#!/usr/bin/env python

"""Run pseudobulk TF-activity inference with decoupler and CollecTRI."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import anndata as ad
import decoupler as dc
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread


PROJECT_DIR = Path(__file__).resolve().parents[1]
MANIFEST_FILE = PROJECT_DIR / "config" / "sample_manifest.csv"
CONFIG_FILE = PROJECT_DIR / "config" / "analysis_config.json"
EXPORT_ROOT = PROJECT_DIR / "results" / "10_tf_activity" / "01_exports"
OUTPUT_ROOT = PROJECT_DIR / "results" / "10_tf_activity" / "02_inference"
TABLE_ROOT = PROJECT_DIR / "TABLE" / "TF_activity"

OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
TABLE_ROOT.mkdir(parents=True, exist_ok=True)


def load_config() -> dict[str, Any]:
    with CONFIG_FILE.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_network(config: dict[str, Any]) -> tuple[pd.DataFrame, Path]:
    network_path = PROJECT_DIR / config["network_file"]

    if network_path.exists():
        network = pd.read_csv(network_path)
    elif bool(config.get("allow_network_download", False)):
        network_path.parent.mkdir(parents=True, exist_ok=True)
        network = dc.op.collectri(
            organism=config["organism"],
            remove_complexes=False,
            license="academic",
            verbose=False,
        )
        network.to_csv(network_path, index=False)
    else:
        raise FileNotFoundError(
            "The cached CollecTRI network was not found at "
            f"{network_path}. Download and archive it before running, or set "
            "allow_network_download=true for the initial retrieval."
        )

    required_columns = {"source", "target", "weight"}
    missing_columns = required_columns.difference(network.columns)

    if missing_columns:
        raise ValueError(
            "The TF network is missing columns: "
            + ", ".join(sorted(missing_columns))
        )

    network = network.loc[:, ["source", "target", "weight"]].copy()
    network["source"] = network["source"].astype(str)
    network["target"] = network["target"].astype(str)
    network["weight"] = pd.to_numeric(network["weight"], errors="raise")

    return network, network_path


def read_exported_object(object_name: str) -> ad.AnnData:
    object_dir = EXPORT_ROOT / object_name

    required_files = [
        object_dir / "counts.mtx",
        object_dir / "genes.tsv",
        object_dir / "barcodes.tsv",
        object_dir / "metadata.csv",
    ]

    missing_files = [str(path) for path in required_files if not path.exists()]

    if missing_files:
        raise FileNotFoundError(
            "Exported input files are missing: " + ", ".join(missing_files)
        )

    counts = mmread(object_dir / "counts.mtx").tocsr().T

    genes = (
        pd.read_csv(
            object_dir / "genes.tsv",
            sep="\t",
            header=None,
            dtype=str,
        )[0]
        .astype(str)
        .tolist()
    )

    barcodes = (
        pd.read_csv(
            object_dir / "barcodes.tsv",
            sep="\t",
            header=None,
            dtype=str,
        )[0]
        .astype(str)
        .tolist()
    )

    metadata = pd.read_csv(
        object_dir / "metadata.csv",
        dtype={"barcode": str},
    )

    if metadata["barcode"].duplicated().any():
        raise ValueError(f"Duplicated metadata barcodes detected: {object_name}")

    metadata = metadata.set_index("barcode").loc[barcodes].copy()

    if counts.shape != (len(barcodes), len(genes)):
        raise ValueError(
            f"Matrix dimensions do not match barcodes and genes for {object_name}"
        )

    adata = ad.AnnData(
        X=counts,
        obs=metadata,
        var=pd.DataFrame(index=pd.Index(genes, name="gene")),
    )

    adata.obs_names = barcodes
    adata.var_names = genes
    adata.var_names_make_unique()

    return adata


def normalize_cpm_log1p(adata: ad.AnnData) -> ad.AnnData:
    adata = adata.copy()
    matrix = adata.X

    if sparse.issparse(matrix):
        library_size = np.asarray(matrix.sum(axis=1)).ravel()
        library_size[library_size == 0] = 1.0
        normalized = matrix.multiply(1e6 / library_size[:, None]).tocsr()
        normalized.data = np.log1p(normalized.data)
        adata.X = normalized
    else:
        library_size = np.asarray(matrix.sum(axis=1)).ravel()
        library_size[library_size == 0] = 1.0
        adata.X = np.log1p((matrix / library_size[:, None]) * 1e6)

    return adata


def build_sample_metadata(
    adata: ad.AnnData,
    sample_col: str,
) -> pd.DataFrame:
    metadata = adata.obs.copy()
    metadata[sample_col] = metadata[sample_col].astype(str)

    rows: list[dict[str, Any]] = []

    for sample_id, frame in metadata.groupby(sample_col, sort=False):
        row: dict[str, Any] = {
            sample_col: str(sample_id),
            "n_cells_total": int(frame.shape[0]),
        }

        for column in metadata.columns:
            values = frame[column].dropna().astype(str).unique().tolist()
            row[column] = values[0] if len(values) == 1 else np.nan

        rows.append(row)

    return pd.DataFrame(rows)


def get_pseudobulk_cell_count_column(pb: ad.AnnData) -> str:
    candidates = (
        "psbulk_n_cells",
        "psbulk_cells",
        "psbulk_n_obs",
    )

    for candidate in candidates:
        if candidate in pb.obs.columns:
            return candidate

    raise KeyError(
        "No pseudobulk cell-count column was found. Available columns: "
        + ", ".join(pb.obs.columns.astype(str))
    )


def run_method(
    pseudobulk: ad.AnnData,
    network: pd.DataFrame,
    method: str,
    min_genes_per_tf: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    if method == "ulm":
        dc.mt.ulm(
            pseudobulk,
            network,
            tmin=min_genes_per_tf,
            verbose=False,
        )
    elif method == "mlm":
        dc.mt.mlm(
            pseudobulk,
            network,
            tmin=min_genes_per_tf,
            verbose=False,
        )
    else:
        raise ValueError(f"Unsupported TF-activity method: {method}")

    score_key = f"score_{method}"
    padj_key = f"padj_{method}"

    if score_key not in pseudobulk.obsm:
        raise KeyError(f"Expected score matrix was not generated: {score_key}")

    if padj_key not in pseudobulk.obsm:
        raise KeyError(f"Expected adjusted-P matrix was not generated: {padj_key}")

    scores = pseudobulk.obsm[score_key].copy()
    adjusted_p = pseudobulk.obsm[padj_key].copy()

    if not isinstance(scores, pd.DataFrame):
        scores = pd.DataFrame(
            scores,
            index=pseudobulk.obs_names,
        )

    if not isinstance(adjusted_p, pd.DataFrame):
        adjusted_p = pd.DataFrame(
            adjusted_p,
            index=pseudobulk.obs_names,
            columns=scores.columns,
        )

    return scores, adjusted_p


def summarize_top_tfs(
    score_matrix: pd.DataFrame,
    top_n: int,
) -> pd.DataFrame:
    return (
        score_matrix.rename_axis("pseudobulk_id")
        .reset_index()
        .melt(
            id_vars="pseudobulk_id",
            var_name="tf",
            value_name="score",
        )
        .sort_values(
            ["pseudobulk_id", "score"],
            ascending=[True, False],
        )
        .groupby("pseudobulk_id", as_index=False)
        .head(top_n)
    )


def infer_grouping(
    adata: ad.AnnData,
    manifest_row: pd.Series,
    grouping_level: str,
) -> tuple[ad.AnnData, str]:
    adata = adata.copy()

    if grouping_level == "celltype":
        grouping_column = str(manifest_row["celltype_col"])
    elif grouping_level == "cluster":
        grouping_column = str(manifest_row["cluster_col"])
    elif grouping_level == "overall":
        grouping_column = "__overall_group"
        adata.obs[grouping_column] = "overall"
    else:
        grouping_column = grouping_level

    if grouping_column not in adata.obs.columns:
        raise KeyError(
            f"Grouping column '{grouping_column}' was not found for "
            f"{manifest_row['object_name']}"
        )

    adata.obs[grouping_column] = adata.obs[grouping_column].astype(str)

    return adata, grouping_column


def run_one_grouping(
    adata: ad.AnnData,
    manifest_row: pd.Series,
    grouping_level: str,
    config: dict[str, Any],
    network: pd.DataFrame,
) -> dict[str, Any]:
    object_name = str(manifest_row["object_name"])
    sample_col = str(manifest_row["sample_col"])

    if sample_col not in adata.obs.columns:
        raise KeyError(
            f"Biological-sample column '{sample_col}' was not found for "
            f"{object_name}"
        )

    working, grouping_column = infer_grouping(
        adata,
        manifest_row,
        grouping_level,
    )

    working.obs[sample_col] = working.obs[sample_col].astype(str)

    pseudobulk = dc.pp.pseudobulk(
        working,
        sample_col=sample_col,
        groups_col=grouping_column,
        mode="sum",
        skip_checks=False,
        verbose=False,
    )

    cell_count_column = get_pseudobulk_cell_count_column(pseudobulk)

    keep = (
        pseudobulk.obs[cell_count_column]
        >= int(config["min_cells_per_pseudobulk"])
    )

    n_before = int(pseudobulk.n_obs)
    pseudobulk = pseudobulk[keep].copy()
    n_after = int(pseudobulk.n_obs)

    output_dir = OUTPUT_ROOT / object_name / grouping_level
    output_dir.mkdir(parents=True, exist_ok=True)

    filtering_table = pd.DataFrame(
        {
            "pseudobulk_id": keep.index.astype(str),
            "n_cells": pseudobulk.obs.reindex(keep.index)[cell_count_column]
            if n_after == n_before
            else np.nan,
            "retained": keep.to_numpy(),
        }
    )

    # Reconstruct a reliable filtering table from the pre-filter object.
    pre_pb = dc.pp.pseudobulk(
        working,
        sample_col=sample_col,
        groups_col=grouping_column,
        mode="sum",
        skip_checks=False,
        verbose=False,
    )
    filtering_table = pd.DataFrame(
        {
            "pseudobulk_id": pre_pb.obs_names.astype(str),
            "n_cells": pre_pb.obs[cell_count_column].to_numpy(),
            "retained": (
                pre_pb.obs[cell_count_column].to_numpy()
                >= int(config["min_cells_per_pseudobulk"])
            ),
        }
    )
    filtering_table.to_csv(
        output_dir / "pseudobulk_filtering.csv",
        index=False,
    )

    if pseudobulk.n_obs == 0:
        return {
            "object_name": object_name,
            "grouping_level": grouping_level,
            "grouping_column": grouping_column,
            "n_pseudobulks_before_filtering": n_before,
            "n_pseudobulks_after_filtering": n_after,
            "status": "no_pseudobulks_retained",
        }

    sample_metadata = build_sample_metadata(
        working,
        sample_col=sample_col,
    )

    if sample_col in pseudobulk.obs.columns:
        pb_sample_ids = pseudobulk.obs[sample_col].astype(str)
    else:
        raise KeyError(
            f"pseudobulk.obs does not contain sample column '{sample_col}'"
        )

    pseudobulk_metadata = (
        pseudobulk.obs.copy()
        .assign(**{sample_col: pb_sample_ids})
        .reset_index(names="pseudobulk_id")
        .merge(
            sample_metadata,
            on=sample_col,
            how="left",
            suffixes=("", "_sample"),
        )
        .set_index("pseudobulk_id")
    )

    pseudobulk.obs = pseudobulk_metadata.copy()

    pseudobulk = normalize_cpm_log1p(pseudobulk)

    pseudobulk.obs.to_csv(
        output_dir / "pseudobulk_metadata.csv",
        index=True,
    )

    method_rows: list[dict[str, Any]] = []

    for method in config["methods"]:
        try:
            scores, adjusted_p = run_method(
                pseudobulk=pseudobulk,
                network=network,
                method=str(method),
                min_genes_per_tf=int(config["min_genes_per_tf"]),
            )

            scores.to_csv(output_dir / f"{method}_scores.csv")
            adjusted_p.to_csv(output_dir / f"{method}_padj.csv")

            summarize_top_tfs(
                scores,
                top_n=int(config["top_n_tfs"]),
            ).to_csv(
                output_dir / f"{method}_top_tfs.csv",
                index=False,
            )

            method_rows.append(
                {
                    "method": method,
                    "status": "completed",
                    "n_tfs": int(scores.shape[1]),
                    "error_message": "",
                }
            )
        except Exception as exc:  # transparent per-method failure
            method_rows.append(
                {
                    "method": method,
                    "status": "failed",
                    "n_tfs": 0,
                    "error_message": f"{type(exc).__name__}: {exc}",
                }
            )

    pd.DataFrame(method_rows).to_csv(
        output_dir / "method_summary.csv",
        index=False,
    )

    return {
        "object_name": object_name,
        "grouping_level": grouping_level,
        "grouping_column": grouping_column,
        "n_pseudobulks_before_filtering": n_before,
        "n_pseudobulks_after_filtering": n_after,
        "status": "completed",
    }


def main() -> None:
    config = load_config()
    manifest = pd.read_csv(MANIFEST_FILE, dtype=str)
    network, network_path = load_network(config)

    network_provenance = pd.DataFrame(
        [
            {
                "network_resource": config.get("network_resource", "CollecTRI"),
                "organism": config["organism"],
                "network_file": str(network_path.relative_to(PROJECT_DIR)),
                "sha256": sha256_file(network_path),
                "n_interactions": int(network.shape[0]),
                "n_tfs": int(network["source"].nunique()),
                "n_targets": int(network["target"].nunique()),
            }
        ]
    )

    network_provenance.to_csv(
        TABLE_ROOT / "TF_network_provenance.csv",
        index=False,
    )

    run_rows: list[dict[str, Any]] = []

    for _, manifest_row in manifest.iterrows():
        object_name = str(manifest_row["object_name"])
        adata = read_exported_object(object_name)

        for grouping_level in config["grouping_levels"]:
            try:
                run_rows.append(
                    run_one_grouping(
                        adata=adata,
                        manifest_row=manifest_row,
                        grouping_level=str(grouping_level),
                        config=config,
                        network=network,
                    )
                )
            except Exception as exc:
                run_rows.append(
                    {
                        "object_name": object_name,
                        "grouping_level": grouping_level,
                        "grouping_column": "",
                        "n_pseudobulks_before_filtering": np.nan,
                        "n_pseudobulks_after_filtering": np.nan,
                        "status": "failed",
                        "error_message": f"{type(exc).__name__}: {exc}",
                    }
                )

    pd.DataFrame(run_rows).to_csv(
        TABLE_ROOT / "TF_inference_run_summary.csv",
        index=False,
    )


if __name__ == "__main__":
    main()
