#!/usr/bin/env python

"""Create generalized TF-activity summary tables and manuscript-oriented plots."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind


PROJECT_DIR = Path(__file__).resolve().parents[1]
MANIFEST_FILE = PROJECT_DIR / "config" / "sample_manifest.csv"
CONFIG_FILE = PROJECT_DIR / "config" / "analysis_config.json"
KEY_TF_FILE = PROJECT_DIR / "config" / "key_tfs.csv"
CONTRAST_FILE = PROJECT_DIR / "config" / "sample_level_contrasts.csv"

INFERENCE_ROOT = PROJECT_DIR / "results" / "10_tf_activity" / "02_inference"
TABLE_ROOT = PROJECT_DIR / "TABLE" / "TF_activity"
FIGURE_ROOT = PROJECT_DIR / "FIGURE" / "TF_activity"

TABLE_ROOT.mkdir(parents=True, exist_ok=True)
FIGURE_ROOT.mkdir(parents=True, exist_ok=True)


def load_config() -> dict[str, Any]:
    with CONFIG_FILE.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def as_flag(value: Any) -> bool:
    return str(value).strip().lower() in {"true", "t", "1", "yes", "y"}


def sanitize_name(value: str) -> str:
    safe = "".join(
        character if character.isalnum() or character in "._-" else "_"
        for character in str(value)
    )
    while "__" in safe:
        safe = safe.replace("__", "_")
    return safe


def row_z_score(matrix: pd.DataFrame) -> pd.DataFrame:
    values = matrix.astype(float).copy()
    means = values.mean(axis=1)
    standard_deviations = values.std(axis=1, ddof=0).replace(0, np.nan)
    scaled = values.sub(means, axis=0).div(standard_deviations, axis=0)
    return scaled.fillna(0.0)


def bh_adjust(p_values: pd.Series) -> pd.Series:
    values = p_values.to_numpy(dtype=float)
    result = np.full(values.shape, np.nan, dtype=float)

    valid = np.isfinite(values)
    valid_values = values[valid]

    if valid_values.size == 0:
        return pd.Series(result, index=p_values.index)

    order = np.argsort(valid_values)
    ranked = valid_values[order]
    n = len(ranked)

    adjusted_ranked = ranked * n / np.arange(1, n + 1)
    adjusted_ranked = np.minimum.accumulate(adjusted_ranked[::-1])[::-1]
    adjusted_ranked = np.clip(adjusted_ranked, 0, 1)

    adjusted_valid = np.empty(n, dtype=float)
    adjusted_valid[order] = adjusted_ranked
    result[valid] = adjusted_valid

    return pd.Series(result, index=p_values.index)


def holm_adjust(p_values: pd.Series) -> pd.Series:
    values = p_values.to_numpy(dtype=float)
    result = np.full(values.shape, np.nan, dtype=float)

    valid = np.isfinite(values)
    valid_values = values[valid]

    if valid_values.size == 0:
        return pd.Series(result, index=p_values.index)

    order = np.argsort(valid_values)
    ranked = valid_values[order]
    n = len(ranked)

    adjusted_ranked = (n - np.arange(n)) * ranked
    adjusted_ranked = np.maximum.accumulate(adjusted_ranked)
    adjusted_ranked = np.clip(adjusted_ranked, 0, 1)

    adjusted_valid = np.empty(n, dtype=float)
    adjusted_valid[order] = adjusted_ranked
    result[valid] = adjusted_valid

    return pd.Series(result, index=p_values.index)


def adjust_p_values(
    p_values: pd.Series,
    method: str,
) -> pd.Series:
    method = method.lower()

    if method in {"none", "raw"}:
        return p_values
    if method in {"bh", "fdr"}:
        return bh_adjust(p_values)
    if method == "holm":
        return holm_adjust(p_values)

    raise ValueError(f"Unsupported p-value adjustment method: {method}")


def read_key_tfs() -> list[str]:
    if not KEY_TF_FILE.exists():
        return []

    table = pd.read_csv(KEY_TF_FILE)

    if "gene" not in table.columns:
        raise ValueError("key_tfs.csv must contain a gene column.")

    if "include" in table.columns:
        table = table.loc[table["include"].map(as_flag)]

    return table["gene"].dropna().astype(str).drop_duplicates().tolist()


def save_heatmap(
    matrix: pd.DataFrame,
    output_stem: Path,
    title: str,
    colorbar_label: str,
) -> None:
    if matrix.empty:
        return

    n_rows, n_columns = matrix.shape

    figure_width = max(6.0, 0.45 * n_columns + 3.0)
    figure_height = max(4.5, 0.32 * n_rows + 2.2)

    fig, ax = plt.subplots(
        figsize=(figure_width, figure_height)
    )

    image = ax.imshow(
        matrix.to_numpy(dtype=float),
        aspect="auto",
        cmap="RdBu_r",
        vmin=-2.5,
        vmax=2.5,
        interpolation="nearest",
    )

    ax.set_xticks(np.arange(n_columns))
    ax.set_xticklabels(
        matrix.columns.astype(str),
        rotation=45,
        ha="right",
        fontsize=8,
    )

    ax.set_yticks(np.arange(n_rows))
    ax.set_yticklabels(
        matrix.index.astype(str),
        fontsize=9,
    )

    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title(title, fontsize=13, fontweight="bold")

    colorbar = fig.colorbar(
        image,
        ax=ax,
        fraction=0.025,
        pad=0.02,
    )
    colorbar.set_label(colorbar_label, fontsize=9)

    fig.tight_layout()

    fig.savefig(
        output_stem.with_suffix(".pdf"),
        bbox_inches="tight",
    )
    fig.savefig(
        output_stem.with_suffix(".tiff"),
        bbox_inches="tight",
        dpi=600,
    )

    plt.close(fig)


def choose_tfs(
    mean_matrix: pd.DataFrame,
    key_tfs: list[str],
    top_n: int,
) -> list[str]:
    available_key_tfs = [
        tf for tf in key_tfs if tf in mean_matrix.index
    ]

    variable_tfs = (
        mean_matrix.std(axis=1)
        .sort_values(ascending=False)
        .index.astype(str)
        .tolist()
    )

    selected = available_key_tfs + [
        tf for tf in variable_tfs if tf not in available_key_tfs
    ]

    return selected[:top_n]


def load_scores_and_metadata(
    object_name: str,
    grouping_level: str,
    method: str = "ulm",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    result_dir = INFERENCE_ROOT / object_name / grouping_level
    score_file = result_dir / f"{method}_scores.csv"
    metadata_file = result_dir / "pseudobulk_metadata.csv"

    if not score_file.exists() or not metadata_file.exists():
        raise FileNotFoundError(
            f"Missing {method} result files for {object_name}/{grouping_level}"
        )

    scores = pd.read_csv(score_file, index_col=0)
    metadata = pd.read_csv(metadata_file, index_col=0)

    common_index = scores.index.intersection(metadata.index)

    if common_index.empty:
        raise ValueError(
            f"No shared pseudobulk IDs for {object_name}/{grouping_level}"
        )

    return scores.loc[common_index], metadata.loc[common_index]


def infer_group_column(
    manifest_row: pd.Series,
    grouping_level: str,
) -> str:
    if grouping_level == "celltype":
        return str(manifest_row["celltype_col"])
    if grouping_level == "cluster":
        return str(manifest_row["cluster_col"])
    if grouping_level == "overall":
        return "__overall_group"
    return grouping_level


def make_group_heatmaps(
    manifest: pd.DataFrame,
    config: dict[str, Any],
    key_tfs: list[str],
) -> list[dict[str, Any]]:
    summary_rows: list[dict[str, Any]] = []

    for _, row in manifest.iterrows():
        object_name = str(row["object_name"])

        for grouping_level in ("celltype", "cluster"):
            try:
                scores, metadata = load_scores_and_metadata(
                    object_name,
                    grouping_level,
                    method="ulm",
                )

                group_column = infer_group_column(row, grouping_level)

                if group_column not in metadata.columns:
                    raise KeyError(
                        f"Grouping column '{group_column}' was not found in "
                        "pseudobulk metadata."
                    )

                score_long = (
                    scores.assign(
                        pseudobulk_id=scores.index.astype(str)
                    )
                    .melt(
                        id_vars="pseudobulk_id",
                        var_name="tf",
                        value_name="score",
                    )
                    .merge(
                        metadata.reset_index().rename(
                            columns={metadata.index.name or "index": "pseudobulk_id"}
                        ),
                        on="pseudobulk_id",
                        how="left",
                    )
                )

                mean_matrix = (
                    score_long.groupby(
                        ["tf", group_column],
                        observed=False,
                    )["score"]
                    .mean()
                    .unstack(group_column)
                )

                selected_tfs = choose_tfs(
                    mean_matrix,
                    key_tfs=key_tfs,
                    top_n=int(config["top_n_tfs"]),
                )

                plot_matrix = row_z_score(
                    mean_matrix.reindex(selected_tfs)
                )

                output_prefix = (
                    f"{sanitize_name(object_name)}_"
                    f"{sanitize_name(grouping_level)}_ULM_activity"
                )

                plot_matrix.to_csv(
                    TABLE_ROOT / f"{output_prefix}_matrix.csv"
                )

                save_heatmap(
                    matrix=plot_matrix,
                    output_stem=FIGURE_ROOT / output_prefix,
                    title=(
                        f"{row['condition_label']} | {grouping_level} | "
                        "ULM TF activity"
                    ),
                    colorbar_label="Row z-score of mean ULM activity",
                )

                summary_rows.append(
                    {
                        "object_name": object_name,
                        "grouping_level": grouping_level,
                        "status": "completed",
                        "n_pseudobulks": int(scores.shape[0]),
                        "n_tfs": int(scores.shape[1]),
                        "n_groups": int(mean_matrix.shape[1]),
                    }
                )
            except Exception as exc:
                summary_rows.append(
                    {
                        "object_name": object_name,
                        "grouping_level": grouping_level,
                        "status": "failed",
                        "n_pseudobulks": np.nan,
                        "n_tfs": np.nan,
                        "n_groups": np.nan,
                        "error_message": f"{type(exc).__name__}: {exc}",
                    }
                )

    return summary_rows


def build_sample_level_data(
    manifest: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    score_frames: list[pd.DataFrame] = []
    metadata_frames: list[pd.DataFrame] = []

    for _, row in manifest.iterrows():
        object_name = str(row["object_name"])

        scores, metadata = load_scores_and_metadata(
            object_name,
            "overall",
            method="ulm",
        )

        scores = scores.copy()
        metadata = metadata.copy()

        scores.index = [
            f"{object_name}::{index}" for index in scores.index.astype(str)
        ]
        metadata.index = scores.index

        metadata["object_name_manifest"] = object_name
        metadata["condition_label_manifest"] = str(row["condition_label"])

        score_frames.append(scores)
        metadata_frames.append(metadata)

    shared_tfs = set(score_frames[0].columns)

    for frame in score_frames[1:]:
        shared_tfs &= set(frame.columns)

    ordered_shared_tfs = [
        tf for tf in score_frames[0].columns if tf in shared_tfs
    ]

    combined_scores = pd.concat(
        [frame.loc[:, ordered_shared_tfs] for frame in score_frames],
        axis=0,
    )

    combined_metadata = pd.concat(
        metadata_frames,
        axis=0,
        sort=False,
    ).loc[combined_scores.index]

    return combined_scores, combined_metadata


def make_sample_heatmap(
    scores: pd.DataFrame,
    metadata: pd.DataFrame,
    key_tfs: list[str],
    top_n: int,
) -> None:
    selected_tfs = choose_tfs(
        scores.T,
        key_tfs=key_tfs,
        top_n=top_n,
    )

    plot_matrix = row_z_score(
        scores.loc[:, selected_tfs].T
    )

    sample_labels: list[str] = []

    for index in plot_matrix.columns:
        condition = metadata.loc[index].get(
            "condition_label_manifest",
            metadata.loc[index].get("condition_label", ""),
        )

        sample_column_candidates = [
            column for column in metadata.columns
            if column.endswith("sample_id_for_tf")
            or column == "sample_id_for_tf"
            or column == "sample_id"
        ]

        if sample_column_candidates:
            sample_label = str(
                metadata.loc[index, sample_column_candidates[0]]
            )
        else:
            sample_label = str(index).split("::", 1)[-1]

        sample_labels.append(
            f"{sample_label}\n{condition}"
        )

    plot_matrix.columns = sample_labels

    plot_matrix.to_csv(
        TABLE_ROOT / "sample_level_ULM_activity_matrix.csv"
    )

    save_heatmap(
        matrix=plot_matrix,
        output_stem=FIGURE_ROOT / "sample_level_ULM_activity",
        title="Sample-level pseudobulk ULM TF activity",
        colorbar_label="Row z-score of ULM activity",
    )


def run_sample_level_contrasts(
    scores: pd.DataFrame,
    metadata: pd.DataFrame,
) -> pd.DataFrame:
    if not CONTRAST_FILE.exists():
        return pd.DataFrame()

    contrast_plan = pd.read_csv(CONTRAST_FILE)
    contrast_plan = contrast_plan.loc[
        contrast_plan["enabled"].map(as_flag)
    ].copy()

    if contrast_plan.empty:
        return pd.DataFrame()

    result_frames: list[pd.DataFrame] = []

    for _, contrast in contrast_plan.iterrows():
        contrast_id = str(contrast["contrast_id"])
        group_column = str(contrast["group_column"])
        group1 = str(contrast["group1"])
        group2 = str(contrast["group2"])
        stratify_column = str(contrast.get("stratify_column", "")).strip()
        stratify_value = str(contrast.get("stratify_value", "")).strip()
        adjustment_method = str(contrast["p_adjust_method"])
        top_n = int(contrast["top_n"])

        if group_column not in metadata.columns:
            raise KeyError(
                f"Contrast '{contrast_id}' uses missing column: {group_column}"
            )

        keep = pd.Series(True, index=metadata.index)

        if stratify_column and stratify_column.lower() != "nan":
            if stratify_column not in metadata.columns:
                raise KeyError(
                    f"Contrast '{contrast_id}' uses missing stratification "
                    f"column: {stratify_column}"
                )

            keep &= (
                metadata[stratify_column].astype(str)
                == stratify_value
            )

        contrast_metadata = metadata.loc[keep].copy()
        contrast_scores = scores.loc[contrast_metadata.index].copy()

        group1_ids = contrast_metadata.index[
            contrast_metadata[group_column].astype(str) == group1
        ]

        group2_ids = contrast_metadata.index[
            contrast_metadata[group_column].astype(str) == group2
        ]

        if len(group1_ids) < 2 or len(group2_ids) < 2:
            raise ValueError(
                f"Contrast '{contrast_id}' requires at least two independent "
                "biological samples in each group."
            )

        rows: list[dict[str, Any]] = []

        for tf in contrast_scores.columns:
            values1 = contrast_scores.loc[group1_ids, tf].astype(float)
            values2 = contrast_scores.loc[group2_ids, tf].astype(float)

            statistic, p_value = ttest_ind(
                values1,
                values2,
                equal_var=False,
                nan_policy="omit",
            )

            rows.append(
                {
                    "contrast_id": contrast_id,
                    "tf": tf,
                    "group1": group1,
                    "group2": group2,
                    "n_group1": int(values1.notna().sum()),
                    "n_group2": int(values2.notna().sum()),
                    "mean_group1": float(values1.mean()),
                    "mean_group2": float(values2.mean()),
                    "delta_group1_minus_group2": float(
                        values1.mean() - values2.mean()
                    ),
                    "t_statistic": float(statistic),
                    "p_value": float(p_value),
                }
            )

        result = pd.DataFrame(rows)
        result["p_adjusted"] = adjust_p_values(
            result["p_value"],
            adjustment_method,
        )
        result["p_adjust_method"] = adjustment_method

        result.to_csv(
            TABLE_ROOT / f"contrast_{sanitize_name(contrast_id)}.csv",
            index=False,
        )

        top_positive = result.nlargest(
            top_n,
            "delta_group1_minus_group2",
        )
        top_negative = result.nsmallest(
            top_n,
            "delta_group1_minus_group2",
        )

        plot_data = (
            pd.concat([top_positive, top_negative], ignore_index=True)
            .drop_duplicates(subset="tf")
            .sort_values("delta_group1_minus_group2")
        )

        fig, ax = plt.subplots(
            figsize=(7, max(4, 0.28 * len(plot_data) + 1.5))
        )

        bar_colors = [
            "#B2182B" if value > 0 else "#2166AC"
            for value in plot_data["delta_group1_minus_group2"]
        ]

        ax.barh(
            plot_data["tf"],
            plot_data["delta_group1_minus_group2"],
            color=bar_colors,
        )
        ax.axvline(0, color="black", linewidth=0.8)
        ax.set_xlabel(
            f"Mean ULM activity difference ({group1} - {group2})"
        )
        ax.set_ylabel("")
        ax.set_title(contrast_id, fontweight="bold")
        fig.tight_layout()

        output_stem = FIGURE_ROOT / f"contrast_{sanitize_name(contrast_id)}"
        fig.savefig(output_stem.with_suffix(".pdf"), bbox_inches="tight")
        fig.savefig(
            output_stem.with_suffix(".tiff"),
            bbox_inches="tight",
            dpi=600,
        )
        plt.close(fig)

        result_frames.append(result)

    combined = pd.concat(result_frames, ignore_index=True)

    combined.to_csv(
        TABLE_ROOT / "sample_level_contrasts_all.csv",
        index=False,
    )

    return combined


def collect_top_tf_tables() -> pd.DataFrame:
    rows: list[pd.DataFrame] = []

    for path in INFERENCE_ROOT.glob("*/*/*_top_tfs.csv"):
        object_name = path.parents[1].name
        grouping_level = path.parent.name
        method = path.name.replace("_top_tfs.csv", "")

        table = pd.read_csv(path)
        table.insert(0, "method", method)
        table.insert(0, "grouping_level", grouping_level)
        table.insert(0, "object_name", object_name)
        rows.append(table)

    if not rows:
        return pd.DataFrame()

    return pd.concat(rows, ignore_index=True)


def main() -> None:
    config = load_config()
    manifest = pd.read_csv(MANIFEST_FILE, dtype=str)
    key_tfs = read_key_tfs()

    plot_summary = make_group_heatmaps(
        manifest=manifest,
        config=config,
        key_tfs=key_tfs,
    )

    scores, metadata = build_sample_level_data(manifest)

    scores.to_csv(
        TABLE_ROOT / "sample_level_ULM_scores.csv"
    )
    metadata.to_csv(
        TABLE_ROOT / "sample_level_pseudobulk_metadata.csv"
    )

    make_sample_heatmap(
        scores=scores,
        metadata=metadata,
        key_tfs=key_tfs,
        top_n=int(config["top_n_tfs"]),
    )

    contrast_results = run_sample_level_contrasts(
        scores=scores,
        metadata=metadata,
    )

    top_tfs = collect_top_tf_tables()

    if not top_tfs.empty:
        top_tfs.to_csv(
            TABLE_ROOT / "combined_top_TFs.csv",
            index=False,
        )

        recurrent = (
            top_tfs.groupby(["method", "tf"])
            .size()
            .reset_index(name="n_top_hits")
            .sort_values(
                ["method", "n_top_hits", "tf"],
                ascending=[True, False, True],
            )
        )

        recurrent.to_csv(
            TABLE_ROOT / "recurrent_top_TFs.csv",
            index=False,
        )

    pd.DataFrame(plot_summary).to_csv(
        TABLE_ROOT / "TF_report_generation_summary.csv",
        index=False,
    )

    report_lines = [
        "# TF Activity Analysis Report",
        "",
        "## Core principles",
        "",
        "- TF activity was inferred from sample-level pseudobulk expression profiles.",
        "- Independent biological samples, not cells, are the units used for contrasts.",
        "- Celltype and cluster heatmaps show means across sample-level pseudobulks.",
        "- UMAP panels generated by the R script show TF gene expression, not activity.",
        "",
        "## Main outputs",
        "",
        "- `sample_level_ULM_scores.csv`",
        "- `sample_level_pseudobulk_metadata.csv`",
        "- `sample_level_ULM_activity.pdf`",
        "- per-object celltype and cluster TF-activity heatmaps",
    ]

    if not contrast_results.empty:
        report_lines.extend(
            [
                "- explicit sample-level contrast tables and bar plots",
            ]
        )

    (TABLE_ROOT / "TF_activity_report.md").write_text(
        "\n".join(report_lines) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
