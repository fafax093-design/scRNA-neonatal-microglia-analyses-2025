# Single-cell RNA-sequencing analyses of neonatal mouse microglia

This repository contains the analysis code supporting a study of early
postnatal microglial states in mice under physiological conditions and following
neonatal lipopolysaccharide challenge. The workflow includes primary
single-cell RNA-sequencing processing, microglial-state annotation, pathway and
trajectory analyses, metabolic-flux post-processing, transcription-factor
activity inference, and integration with an external developmental microglia
dataset.

> **Manuscript title:** Neonatal inflammation disrupts a temporally restricted neurodevelopment-associated microglial program in the mouse brain

## Study overview

The primary single-cell dataset contains 13 neonatal mouse samples collected at
postnatal days 3, 7, and 12 under physiological saline and LPS-exposed
conditions. The repository also includes code for cross-dataset comparison with
the Hammond developmental microglia dataset.

The analysis is organized into 11 numbered stages. Several downstream stages
are configuration driven and can process multiple Seurat objects, annotation
schemes, pathways, or comparisons without editing the main analysis scripts.

## Workflow overview

```text
01  Load and merge primary scRNA-seq data
 ↓
02  Cell-level quality control
 ↓
03  SCTransform, PCA, Harmony integration, UMAP, and clustering
 ↓
04  Marker genes, selected DEG analyses, and reviewed annotations
 ├── 05  Visualization
 ├── 06  Gene Ontology enrichment
 ├── 07  UCell pathway scoring
 ├── 08  Monocle 3 trajectory and pseudotime analysis
 ├── 09  scFEA flux post-processing
 └── 10  Pseudobulk transcription-factor activity analysis

11  Independent cross-dataset label transfer, RPCA integration,
    and Figure I-J generation
```

Stages 01–04 form the principal preprocessing and annotation workflow. Stages
05–10 use annotated objects and may be run independently when their required
inputs and configuration files are available. Stage 11 is a separate
cross-dataset analysis requiring both the physiological reference object and
the Hammond microglia object.

## Repository structure

```text
scRNA-neonatal-microglia-analyses-2025/
├── scripts/
│   ├── 01_load_and_merge_scRNA_data.R
│   ├── 02_quality_control.R
│   ├── 03_integration_and_clustering.R
│   ├── 04_marker_DEG_and_annotation.R
│   ├── 05_visualization.R
│   ├── 06_GO_enrichment.R
│   ├── 07_UCell_pathway_scoring.R
│   ├── 08a_monocle3_trajectory.R
│   ├── 08b_monocle3_graph_test_and_modules.R
│   ├── 08c_monocle3_pseudotime_visualization.R
│   ├── 08d_optional_ClusterGVis_heatmap.R
│   ├── 09_scFEA_flux_postprocessing.R
│   ├── 10a_validate_and_export.R
│   ├── 10b_run_tf_activity.py
│   ├── 10c_make_tf_report.py
│   ├── 10d_make_key_tf_expression_umaps.R
│   ├── 10e_package_tf_results.R
│   └── 11_cross_dataset_integration.R
├── config/
│   ├── 04_marker_annotation_manifest.csv
│   ├── 04_DEG_comparisons.csv
│   ├── 05_visualization_manifest.csv
│   ├── 06_GO_enrichment_manifest.csv
│   ├── 07_UCell_analysis_manifest.csv
│   ├── annotations/
│   ├── visualization/
│   ├── GO/
│   ├── UCell/
│   ├── monocle3/
│   ├── scFEA/
│   ├── sample_manifest.csv
│   ├── analysis_config.json
│   ├── key_tfs.csv
│   └── sample_level_contrasts.csv
├── resources/
│   └── collectri_mouse.csv
├── DATA/
│   ├── INPUT/
│   └── OUTPUT/
├── data/
├── output/
├── TABLE/
├── FIGURE/
├── results/
├── run_tf_pipeline.sh
├── requirements_tf.txt
├── README.md
├── LICENSE
├── CITATION.cff
└── .gitignore
```

Most primary analyses use the uppercase `DATA/`, `TABLE/`, and `FIGURE/`
directories. The cross-dataset script retains the standalone `data/` and
`output/` paths used by that workflow.

Large FASTQ files, Cell Ranger output directories, full processed Seurat
objects, and other large intermediate files are not intended to be committed to
GitHub.

## Analysis stages

| Stage | Script | Description |
|---|---|---|
| 01 | `01_load_and_merge_scRNA_data.R` | Loads preconstructed Seurat objects or 10X Genomics filtered feature-barcode matrices, standardizes metadata, prefixes cell barcodes with sample identifiers, and merges the 13 primary samples. |
| 02 | `02_quality_control.R` | Reconstructs the RNA assay, calculates mitochondrial, ribosomal, and haemoglobin transcript percentages, applies the prespecified QC thresholds, and exports cell-retention and QC summaries. |
| 03 | `03_integration_and_clustering.R` | Performs SCTransform normalization and PCA once, then creates parallel Harmony-integrated and uncorrected branches for UMAP, neighbour-graph construction, and clustering across prespecified resolutions. |
| 04 | `04_marker_DEG_and_annotation.R` | Processes one or more Seurat objects using external manifests, identifies cluster markers, performs configured DEG comparisons, applies reviewed cluster-to-state annotations, and exports annotated objects. |
| 05 | `05_visualization.R` | Generates configurable cluster trees, cell-state UMAPs, group-split UMAPs, feature plots, composition plots, expression boxplots, heatmaps, and marker-gene DotPlots. |
| 06 | `06_GO_enrichment.R` | Performs configurable mouse Gene Ontology enrichment for multiple marker or DEG tables and generates complete enrichment tables and annotated bubble plots. |
| 07 | `07_UCell_pathway_scoring.R` | Calculates UCell pathway scores for multiple objects or pathways, exports cell- and biological-sample-level summaries, and performs sample-level comparisons. |
| 08a | `08a_monocle3_trajectory.R` | Builds Monocle 3 cell-data-set objects, imports the selected Seurat UMAP, learns the principal graph, assigns a reviewed root population, calculates pseudotime, and transfers pseudotime back to Seurat. |
| 08b | `08b_monocle3_graph_test_and_modules.R` | Runs `graph_test()`, filters pseudotime-associated genes, identifies gene modules, optionally applies reviewed module-merging rules, and summarizes module expression by cell state. |
| 08c | `08c_monocle3_pseudotime_visualization.R` | Generates pseudotime expression curves and a module-aware, pseudotime-binned gene-expression heatmap. |
| 08d | `08d_optional_ClusterGVis_heatmap.R` | Optionally generates the ClusterGVis heatmap; its k-means clusters are distinct from Monocle 3 gene modules. |
| 09 | `09_scFEA_flux_postprocessing.R` | Matches per-cell scFEA flux estimates to Seurat metadata, summarizes flux by configured cell states and biological samples, applies reviewed module orders, and generates full and selected-module heatmaps. |
| 10a–10e | TF-activity subworkflow | Exports raw counts and metadata, constructs biological-sample-level pseudobulk profiles, infers TF activity with CollecTRI and decoupler ULM/MLM, generates summary figures, plots selected TF gene expression on UMAPs, and packages verified result bundles. |
| 11 | `11_cross_dataset_integration.R` | Transfers reference microglial-state labels to Hammond cells, integrates the datasets with SCT and reciprocal PCA, generates a shared UMAP, exports composition tables, and produces Figure I-J. |

## Primary sample identifiers

The primary scRNA-seq samples use the following standardized identifiers:

```text
NS_P3_1
NS_P3_2
NS_P3_3
LPS_P3_1
LPS_P3_2
NS_P7_1
NS_P7_2
LPS_P7_1
LPS_P7_2
NS_P12_1
NS_P12_2
LPS_P12_1
LPS_P12_2
```

These identifiers should remain consistent across sequencing metadata, Seurat
metadata, processed objects, tables, and figures.

The technical-batch assignments used in the archived workflow are:

```text
batch1: NS_P3_3

batch2:
  NS_P3_1
  NS_P3_2
  LPS_P3_1
  LPS_P3_2

batch3:
  NS_P12_1
  NS_P12_2
  LPS_P12_1
  LPS_P12_2

batch4:
  NS_P7_1
  NS_P7_2
  LPS_P7_1
  LPS_P7_2
```

A technical batch must not be used as an independent biological-sample
identifier.

## Input data

The primary workflow expects locally available input files under `DATA/INPUT/`.
A representative layout is:

```text
DATA/
├── INPUT/
│   ├── RDS/
│   │   ├── NS_P3_1_seurat.rds
│   │   ├── NS_P3_2_seurat.rds
│   │   ├── NS_P3_3_seurat.rds
│   │   ├── LPS_P3_1_seurat.rds
│   │   └── LPS_P3_2_seurat.rds
│   ├── 10X/
│   │   ├── NS_P7_1/filtered_feature_bc_matrix/
│   │   ├── NS_P7_2/filtered_feature_bc_matrix/
│   │   ├── LPS_P7_1/filtered_feature_bc_matrix/
│   │   ├── LPS_P7_2/filtered_feature_bc_matrix/
│   │   ├── NS_P12_1/filtered_feature_bc_matrix/
│   │   ├── NS_P12_2/filtered_feature_bc_matrix/
│   │   ├── LPS_P12_1/filtered_feature_bc_matrix/
│   │   └── LPS_P12_2/filtered_feature_bc_matrix/
│   ├── external/
│   └── scFEA/
└── OUTPUT/
```

Script-specific input paths and required metadata columns are documented at the
beginning of each script.

## Configuration-driven analyses

Stages 04–10 use external configuration files so that multiple objects,
annotation schemes, pathways, or comparisons can be processed without
hard-coding study-specific values in the main scripts.

### Marker, DEG, and annotation configuration

```text
config/04_marker_annotation_manifest.csv
config/04_DEG_comparisons.csv
config/annotations/
```

Each manifest row identifies one Seurat object, its clustering column,
annotation file, output object, marker assay, and requested analyses.

### Visualization configuration

```text
config/05_visualization_manifest.csv
config/visualization/
```

Object-specific color palettes and feature, boxplot, heatmap, and DotPlot gene
lists are maintained as external tables.

### GO enrichment configuration

```text
config/06_GO_enrichment_manifest.csv
config/GO/
```

Each row defines the input DEG table, grouping and gene columns, significance
thresholds, ontology, optional GO-term annotations, and plot dimensions.

### UCell configuration

```text
config/07_UCell_analysis_manifest.csv
config/UCell/
```

Each row defines the object, pathway, grouping column, independent
biological-sample column, group order, comparison plan, and statistical method.

### Monocle 3 configuration

```text
config/monocle3/
├── 08a_trajectory_manifest.csv
├── 08b_graph_test_manifest.csv
├── 08c_visualization_manifest.csv
└── 08d_ClusterGVis_manifest.csv
```

Root populations, graph-test thresholds, module-merging rules, pseudotime
visualization parameters, and optional ClusterGVis settings are stored outside
the analysis code.

### scFEA configuration

```text
config/scFEA/
├── 09_scFEA_manifest.csv
├── <analysis>_state_annotations.csv
├── <analysis>_module_order.csv
├── <analysis>_selected_modules.csv
└── <analysis>_scFEA_command.txt
```

The R workflow performs scFEA result post-processing only. The exact Python
command, scFEA version or commit, input expression matrix, model files, and
parameters used to generate each flux matrix should be archived separately.

### TF-activity configuration

```text
config/sample_manifest.csv
config/analysis_config.json
config/key_tfs.csv
config/sample_level_contrasts.csv
resources/collectri_mouse.csv
```

The CollecTRI regulatory-network file is cached locally, and its SHA-256 checksum
and dimensions are recorded during the analysis.

## Software requirements

### R

The workflow was developed around Seurat-based single-cell analysis. Major R
packages include:

```text
Seurat
SeuratObject
dplyr
tidyr
tibble
readr
readxl
openxlsx
ggplot2
ggpubr
patchwork
future
harmony
SCP
clustree
ggrastr
ComplexHeatmap
circlize
clusterProfiler
org.Mm.eg.db
UCell
monocle3
SingleCellExperiment
ClusterGVis
pheatmap
RColorBrewer
Matrix
jsonlite
purrr
```

Some scripts contain compatibility handling for Seurat v4 and v5 data access.
The exact environment used for the archived analysis should be reported through
the generated `sessionInfo()` files.

### Python

The TF-activity workflow requires Python and the packages listed in
`requirements_tf.txt`, including:

```text
anndata
decoupler
matplotlib
numpy
pandas
scipy
```

scFEA is an external Python workflow. Its exact software environment and command
must be recorded with the corresponding flux input.

## Running the workflow

Run commands from the repository root.

### Primary preprocessing and annotation

```bash
Rscript scripts/01_load_and_merge_scRNA_data.R
Rscript scripts/02_quality_control.R
Rscript scripts/03_integration_and_clustering.R
Rscript scripts/04_marker_DEG_and_annotation.R
```

Complete all required annotation manifests and reviewed mapping tables before
running stage 04.

### Visualization and downstream analyses

```bash
Rscript scripts/05_visualization.R
Rscript scripts/06_GO_enrichment.R
Rscript scripts/07_UCell_pathway_scoring.R
```

### Monocle 3 trajectory workflow

```bash
Rscript scripts/08a_monocle3_trajectory.R
Rscript scripts/08b_monocle3_graph_test_and_modules.R
Rscript scripts/08c_monocle3_pseudotime_visualization.R
Rscript scripts/08d_optional_ClusterGVis_heatmap.R
```

Stage 08d is optional.

### scFEA post-processing

```bash
Rscript scripts/09_scFEA_flux_postprocessing.R
```

This command requires a completed per-cell flux matrix generated by scFEA.

### TF-activity workflow

```bash
bash run_tf_pipeline.sh
```

The shell workflow runs stages 10a–10e in order. Set the `RSCRIPT` and `PYTHON`
environment variables when non-default interpreters are required.

### Cross-dataset integration

```bash
Rscript scripts/11_cross_dataset_integration.R
```

The script expects the reference and Hammond objects under `data/` and writes
its integrated object, tables, figures, parameters, and `sessionInfo()` under
`output/`.

## Statistical unit and pseudoreplication

Cells from the same animal or biological sample are not independent biological
replicates.

Accordingly:

- UCell inferential comparisons use one summary score per independent sample;
- TF activity is inferred from sample-level pseudobulk profiles;
- scFEA state-level heatmaps are descriptive, while sample-level summaries are
  exported for statistical analysis;
- cell-level UMAPs, feature plots, and score distributions are primarily
  descriptive visualizations.

Metadata columns such as `group`, `treatment`, `age`, or `batch` must not replace
the independent biological-sample identifier unless they genuinely contain
unique sample IDs.

## Important methodological notes

### Harmony and uncorrected branches

Stage 03 generates both:

```text
DATA/OUTPUT/scRNA_HarmonyIntegrated.rds
DATA/OUTPUT/scRNA_NoBatchCorrection.rds
```

The Harmony-integrated object is used for the principal downstream annotation
workflow, while the uncorrected object is retained for diagnostic comparison.

### Reviewed cell-state annotations

Cluster-to-state mappings are stored in external CSV files and validated before
being applied. Annotation labels are not inferred automatically from cluster
numbers alone.

### Monocle 3 input assay

The archived Monocle 3 configuration may use `SCT` counts to reproduce the
original analysis. Raw `RNA` counts are more conventional for a new analysis,
but the assay should not be changed in the archived workflow without confirming
which input generated the reported result.

### Monocle 3 modules and ClusterGVis clusters

Gene modules generated by `monocle3::find_gene_modules()` and k-means clusters
generated by ClusterGVis are separate analytical results and should not be
treated as equivalent.

### TF activity and TF expression

ULM and MLM outputs represent inferred transcription-factor activity.

The UMAP panels generated by `10d_make_key_tf_expression_umaps.R` display RNA
expression of selected TF genes. They are not single-cell TF-activity maps and
should be labeled as expression panels.

### Cross-dataset integration

The integrated assay generated in stage 11 is used for cross-dataset embedding
and visualization. Reference labels are transferred before SCT-RPCA integration,
and prediction-score summaries are exported where available.

## Main outputs

Generated files are organized by analysis stage:

```text
DATA/OUTPUT/       processed Seurat objects and compact result bundles
TABLE/             QC, marker, DEG, enrichment, score, and summary tables
FIGURE/            manuscript-oriented PDF and 600-dpi TIFF figures
results/           intermediate TF-activity exports and inference outputs
output/            cross-dataset integrated object, tables, and Figure I-J
```

Most plotting workflows export both the figure and its underlying numerical
table. Major scripts also record analysis parameters and software versions.

Generated large files should not be committed to GitHub unless required for
reproducibility and compatible with repository file-size limits.

## Reproducibility

The repository follows these principles:

- standardized sample identifiers and sample-prefixed cell barcodes;
- explicit technical-batch assignments;
- input-file and metadata validation;
- fixed random seeds for stochastic analyses;
- configuration files for reviewed annotation and plotting choices;
- preserved underlying data for major figures;
- intermediate objects saved between major stages;
- cached regulatory-network resources with checksums;
- sample-level statistical summaries where biological replication is required;
- `sessionInfo()` or environment records for major analysis stages.

Before creating the archived release, all scripts should be run once from a
clean project directory using the final input files and configuration tables.

## Data availability

The transcriptome sequencing data are available in the BioSample database
under BioProject ID PRJNA1329933.

> **Accession:** BioProject PRJNA1329933

Large sequencing files, Cell Ranger output directories, and processed Seurat
objects are not stored in this GitHub repository.

The external Hammond developmental microglia dataset should be obtained from its
original public source and cited according to the associated publication.

## Code availability

The analysis code is available in this GitHub repository. A permanent versioned
release will be archived in Zenodo and assigned a DOI.

> **Zenodo DOI:** Pending assignment after the first archived release

## Citation

When using this code, cite both the associated manuscript and the archived
software release. Citation metadata will also be provided in `CITATION.cff`.

Recommended software citation format:

```text
Zhu, J. et al. [FINAL MANUSCRIPT TITLE].
scRNA-neonatal-microglia-analyses-2025, version 1.0.0.
Zenodo. DOI to be added after the first archived release (2026).
```

## License

Reuse of the code is governed by the license provided in `LICENSE`.

Before creating the final Zenodo release, confirm the selected license with the
corresponding author and institution.

## Contact

For questions regarding the analysis, contact:

**Jinjin Zhu**

jinjin337@gs.zzu.edu.cn
