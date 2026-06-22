
##############################################

### Description: Single-cell RNA-seq analysis pipeline
### Author:Jinjin zhu

##############################################

### Load required packages
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(ggsci)
library(gplots)

###-----------------------------
### # 1.  Load count matrices and construct Seurat objects #
# Description: # This script loads neonatal mouse brain immune-cell datasets
#from preconstructed # Seurat objects or feature–barcode count matrices,
#adds sample metadata, and
# assigns unique cell barcodes before downstream integration.
###-----------------------------
#   (1) Load packages

# ==============================================================================
# Load raw single-cell RNA-seq data and construct Seurat objects
#
# Description:
# This script loads neonatal mouse brain immune-cell datasets from Seurat RDS
# files or 10X Genomics filtered feature-barcode matrices, adds sample metadata,
# and assigns unique cell barcodes before downstream integration.
# ==============================================================================


# 1 Load packages -------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})



# 2. Define project directories ------------------------------------------------

rds_dir <- "DATA/INPUT/RDS"
matrix_dir <- "DATA/INPUT/10X"
output_dir <- "DATA/OUTPUT"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# 3. Define sample information -------------------------------------------------

sample_information <- data.frame(
  sample_id = c(
    "LPS_P3_1",
    "LPS_P3_2",
    "NS_P3_1",
    "NS_P3_2",
    "NS_P3_3",
    "LPS_P7_1",
    "LPS_P7_2",
    "NS_P7_1",
    "NS_P7_2",
    "LPS_P12_1",
    "LPS_P12_2",
    "NS_P12_1",
    "NS_P12_2"
  ),

  treatment = c(
    "LPS", "LPS", "NS", "NS","NS",
    "LPS", "LPS", "NS", "NS",
    "LPS", "LPS", "NS", "NS"
  ),

  age = c(
    "P3", "P3", "P3", "P3","P3",
    "P7", "P7", "P7", "P7",
    "P12", "P12", "P12", "P12"
  ),

  replicate = c(
    1, 2, 1, 2,
    1, 2, 1, 2,
    1, 2, 1, 2
  ),

  input_type = c(
    "rds", "rds", "rds", "rds",
    "10x", "10x", "10x", "10x",
    "10x", "10x", "10x", "10x"
  ),

  input_path = c(
    file.path(rds_dir, "5L_1_seurat.rds"),
    file.path(rds_dir, "5L_2_seurat.rds"),
    file.path(rds_dir, "NS_1_seurat.rds"),
    file.path(rds_dir, "NS_2_seurat.rds"),

    file.path(
      matrix_dir,
      "2308154_LPS_1_P7",
      "filtered_feature_bc_matrix"
    ),
    file.path(
      matrix_dir,
      "2308154_LPS_2_P7",
      "filtered_feature_bc_matrix"
    ),
    file.path(
      matrix_dir,
      "2308154_NS_1_P7",
      "filtered_feature_bc_matrix"
    ),
    file.path(
      matrix_dir,
      "2308154_NS_2_P7",
      "filtered_feature_bc_matrix"
    ),

    file.path(
      matrix_dir,
      "2308153_5LPS_1_P12",
      "filtered_feature_bc_matrix"
    ),
    file.path(
      matrix_dir,
      "2308153_5LPS_2_P12",
      "filtered_feature_bc_matrix"
    ),
    file.path(
      matrix_dir,
      "2308153_NS_1_P12",
      "filtered_feature_bc_matrix"
    ),
    file.path(
      matrix_dir,
      "2308153_NS_2_P12",
      "filtered_feature_bc_matrix"
    )
  ),

  stringsAsFactors = FALSE
)


# 4. Function for loading one sample -------------------------------------------

load_single_sample <- function(
    sample_id,
    treatment,
    age,
    replicate,
    input_type,
    input_path,
    min_cells = 3,
    min_features = 200
) {

  if (!file.exists(input_path)) {
    stop(
      paste0(
        "Input file or directory not found: ",
        input_path
      )
    )
  }

  if (input_type == "rds") {

    seurat_object <- readRDS(input_path)

    if (!inherits(seurat_object, "Seurat")) {
      stop(
        paste0(
          "The RDS file does not contain a Seurat object: ",
          input_path
        )
      )
    }

  } else if (input_type == "10x") {

    count_matrix <- Read10X(
      data.dir = input_path
    )

    # If Read10X returns several assays, retain the gene-expression matrix.
    if (is.list(count_matrix)) {

      if ("Gene Expression" %in% names(count_matrix)) {
        count_matrix <- count_matrix[["Gene Expression"]]
      } else {
        count_matrix <- count_matrix[[1]]
      }
    }

    seurat_object <- CreateSeuratObject(
      counts = count_matrix,
      project = sample_id,
      min.cells = min_cells,
      min.features = min_features
    )

  } else {

    stop(
      paste0(
        "Unsupported input type for ",
        sample_id,
        ": ",
        input_type
      )
    )
  }

  # Add a unique sample prefix to every cell barcode.
  seurat_object <- RenameCells(
    object = seurat_object,
    add.cell.id = sample_id
  )

  # Add sample-level metadata.
  seurat_object$sample_id <- sample_id
  seurat_object$treatment <- treatment
  seurat_object$age <- age
  seurat_object$replicate <- replicate
  seurat_object$group <- paste(
    treatment,
    age,
    sep = "_"
  )

  return(seurat_object)
}


# 5. Load all samples -----------------------------------------------------------

seurat_list <- lapply(
  seq_len(nrow(sample_information)),
  function(i) {

    load_single_sample(
      sample_id = sample_information$sample_id[i],
      treatment = sample_information$treatment[i],
      age = sample_information$age[i],
      replicate = sample_information$replicate[i],
      input_type = sample_information$input_type[i],
      input_path = sample_information$input_path[i]
    )
  }
)

names(seurat_list) <- sample_information$sample_id


# 6. Check the loaded datasets -------------------------------------------------

sample_summary <- do.call(
  rbind,
  lapply(
    names(seurat_list),
    function(sample_name) {

      current_object <- seurat_list[[sample_name]]

      data.frame(
        sample_id = sample_name,
        cells = ncol(current_object),
        genes = nrow(current_object),
        treatment = unique(current_object$treatment),
        age = unique(current_object$age),
        stringsAsFactors = FALSE
      )
    }
  )
)

print(sample_summary)

write.csv(
  sample_summary,
  file = file.path(
    output_dir,
    "sample_loading_summary.csv"
  ),
  row.names = FALSE
)


# 7. Merge samples --------------------------------------------------------------

scRNA_raw <- merge(
  x = seurat_list[[1]],
  y = seurat_list[-1],
  project = "Neonatal_brain_immune_cells"
)

print(
  table(
    scRNA_raw$age,
    scRNA_raw$treatment
  )
)

print(
  table(
    scRNA_raw$sample_id
  )
)


# 8. Save the merged raw Seurat object -----------------------------------------

saveRDS(
  scRNA_raw,
  file = file.path(
    output_dir,
    "merged_raw_Seurat_object.rds"
  )
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(
    output_dir,
    "Seurat_data_loading_sessionInfo.txt"
  )
)


# ==============================================================================
# (2)  Merge and quality control of single-cell RNA-seq data
#
# Description:
# This script rebuilds a clean Seurat object from the merged raw RNA-count
# matrix, calculates mitochondrial, ribosomal and haemoglobin transcript
# percentages, applies the prespecified cell-level QC thresholds, assigns batch
# labels, exports QC summaries and saves the filtered Seurat object.
# ==============================================================================


# 1. Load packages -------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
})


# 2. Define paths ---------------------------------------------------------------

input_file <- "DATA/OUTPUT/merged_raw_Seurat_object.rds"
output_data_dir <- "DATA/OUTPUT"
output_qc_dir <- "QC"

dir.create(output_data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_qc_dir, recursive = TRUE, showWarnings = FALSE)


# 3. Define QC thresholds -------------------------------------------------------

min_features <- 500
max_features <- 3000
max_counts <- 22000
max_percent_mt <- 10


# 4. Load merged raw Seurat object ---------------------------------------------

merged_seurat <- readRDS(input_file)

if (!inherits(merged_seurat, "Seurat")) {
  stop("The input file must contain a Seurat object.")
}

if (!"RNA" %in% names(merged_seurat@assays)) {
  stop("The input Seurat object does not contain an RNA assay.")
}

if (!"orig.ident" %in% colnames(merged_seurat@meta.data)) {
  stop("The metadata must contain a column named 'orig.ident'.")
}


# 5. Rebuild a clean Seurat object from raw RNA counts -------------------------

# Compatible with Seurat v5 and Seurat v4.
rna_counts <- tryCatch(
  GetAssayData(
    object = merged_seurat,
    assay = "RNA",
    layer = "counts"
  ),
  error = function(e) {
    GetAssayData(
      object = merged_seurat,
      assay = "RNA",
      slot = "counts"
    )
  }
)

# Retain sample-level metadata, but remove assay-derived QC columns that
# CreateSeuratObject recalculates from the raw count matrix.
cell_metadata <- merged_seurat@meta.data[
  colnames(rna_counts),
  ,
  drop = FALSE
]

recalculated_metadata_columns <- intersect(
  c(
    "nCount_RNA",
    "nFeature_RNA",
    "percent.mt",
    "percent.rb",
    "percent.HB"
  ),
  colnames(cell_metadata)
)

cell_metadata <- cell_metadata[
  ,
  setdiff(
    colnames(cell_metadata),
    recalculated_metadata_columns
  ),
  drop = FALSE
]

scRNA <- CreateSeuratObject(
  counts = rna_counts,
  project = "Neonatal_brain_immune_cells",
  meta.data = cell_metadata,
  min.cells = 0,
  min.features = 0
)

rm(merged_seurat, rna_counts, cell_metadata)


# 6. Calculate QC metrics -------------------------------------------------------

# Percentage of transcripts mapped to mouse mitochondrial genes.
scRNA[["percent.mt"]] <- PercentageFeatureSet(
  object = scRNA,
  pattern = "^mt-"
)

# Percentage of transcripts mapped to ribosomal protein genes.
scRNA[["percent.rb"]] <- PercentageFeatureSet(
  object = scRNA,
  pattern = "^Rp[sl]"
)

# Percentage of transcripts mapped to haemoglobin genes.
hb_genes <- c(
  "Hba-a1",
  "Hba-a2",
  "Hbb-b1",
  "Hbb-b2",
  "Hbe1",
  "Hbg1",
  "Hbg2",
  "Hbm",
  "Hbq1",
  "Hbz"
)

hb_genes_present <- CaseMatch(
  search = hb_genes,
  match = rownames(scRNA)
)

if (length(hb_genes_present) == 0) {
  warning("None of the predefined haemoglobin genes were found.")
  scRNA$percent.HB <- 0
} else {
  scRNA[["percent.HB"]] <- PercentageFeatureSet(
    object = scRNA,
    features = hb_genes_present
  )
}


# 7. Assign experimental batches -----------------------------------------------

sample_to_batch <- c(
  "NS" = "batch1",
  "Lps" = "batch1",

  "Lps_5_1" = "batch2",
  "Lps_5_2" = "batch2",
  "NS_1" = "batch2",
  "NS_2" = "batch2",

  "NS_P12_1" = "batch3",
  "NS_P12_2" = "batch3",
  "LPS_5_LPS_P12_1" = "batch3",
  "LPS_5_LPS_P12_2" = "batch3",

  "NS_P7_1" = "batch4",
  "NS_P7_2" = "batch4",
  "LPS_5_P7_1" = "batch4",
  "LPS_5_P7_2" = "batch4"
)

sample_ids <- as.character(scRNA$orig.ident)

unmatched_samples <- setdiff(
  unique(sample_ids),
  names(sample_to_batch)
)

if (length(unmatched_samples) > 0) {
  stop(
    paste0(
      "The following orig.ident values have no batch assignment: ",
      paste(unmatched_samples, collapse = ", ")
    )
  )
}

scRNA$batch <- unname(
  sample_to_batch[sample_ids]
)


# 8. QC-summary function --------------------------------------------------------

summarise_qc <- function(object, stage) {
  object@meta.data |>
    rownames_to_column("cell_barcode") |>
    group_by(orig.ident) |>
    summarise(
      stage = stage,
      n_cells = n(),

      median_nFeature_RNA = median(nFeature_RNA, na.rm = TRUE),
      q1_nFeature_RNA = quantile(nFeature_RNA, 0.25, na.rm = TRUE),
      q3_nFeature_RNA = quantile(nFeature_RNA, 0.75, na.rm = TRUE),

      median_nCount_RNA = median(nCount_RNA, na.rm = TRUE),
      q1_nCount_RNA = quantile(nCount_RNA, 0.25, na.rm = TRUE),
      q3_nCount_RNA = quantile(nCount_RNA, 0.75, na.rm = TRUE),

      median_percent_mt = median(percent.mt, na.rm = TRUE),
      q1_percent_mt = quantile(percent.mt, 0.25, na.rm = TRUE),
      q3_percent_mt = quantile(percent.mt, 0.75, na.rm = TRUE),

      median_percent_rb = median(percent.rb, na.rm = TRUE),
      median_percent_HB = median(percent.HB, na.rm = TRUE),

      .groups = "drop"
    ) |>
    relocate(stage, orig.ident)
}


# 9. Export and visualize pre-QC metrics ----------------------------------------

qc_summary_pre <- summarise_qc(
  object = scRNA,
  stage = "Pre-QC"
)

qc_features <- c(
  "nFeature_RNA",
  "nCount_RNA",
  "percent.mt",
  "percent.rb",
  "percent.HB"
)

qc_plots_pre <- lapply(
  qc_features,
  function(feature_name) {
    VlnPlot(
      object = scRNA,
      group.by = "orig.ident",
      features = feature_name,
      pt.size = 0,
      raster = FALSE
    ) +
      NoLegend() +
      ggtitle(feature_name) +
      theme(
        axis.text.x = element_text(
          angle = 45,
          hjust = 1
        )
      )
  }
)

ggsave(
  filename = file.path(
    output_qc_dir,
    "QC_metrics_before_filtering.pdf"
  ),
  plot = wrap_plots(qc_plots_pre, nrow = 3),
  width = 16,
  height = 8,
  units = "in"
)


# 10. Visualize the prespecified QC thresholds ---------------------------------

p_feature_threshold <- VlnPlot(
  scRNA,
  features = "nFeature_RNA",
  pt.size = 0,
  raster = FALSE
) +
  geom_hline(
    yintercept = c(min_features, max_features),
    linetype = "dashed"
  )

p_count_threshold <- VlnPlot(
  scRNA,
  features = "nCount_RNA",
  pt.size = 0,
  raster = FALSE
) +
  geom_hline(
    yintercept = max_counts,
    linetype = "dashed"
  )

p_mt_threshold <- VlnPlot(
  scRNA,
  features = "percent.mt",
  pt.size = 0,
  raster = FALSE
) +
  geom_hline(
    yintercept = max_percent_mt,
    linetype = "dashed"
  )

ggsave(
  filename = file.path(
    output_qc_dir,
    "QC_filtering_thresholds.pdf"
  ),
  plot = p_feature_threshold +
    p_count_threshold +
    p_mt_threshold,
  width = 12,
  height = 4,
  units = "in"
)


# 11. Apply cell-level QC thresholds -------------------------------------------

cells_before_qc <- table(scRNA$orig.ident)

scRNA <- subset(
  x = scRNA,
  subset =
    nFeature_RNA > min_features &
    nFeature_RNA < max_features &
    nCount_RNA < max_counts &
    percent.mt < max_percent_mt
)

cells_after_qc <- table(scRNA$orig.ident)


# 12. Export cell-retention summary --------------------------------------------

all_samples <- union(
  names(cells_before_qc),
  names(cells_after_qc)
)

cell_retention <- data.frame(
  orig.ident = all_samples,
  cells_before_qc = as.integer(
    cells_before_qc[all_samples]
  ),
  cells_after_qc = as.integer(
    cells_after_qc[all_samples]
  ),
  stringsAsFactors = FALSE
)

cell_retention$cells_before_qc[
  is.na(cell_retention$cells_before_qc)
] <- 0

cell_retention$cells_after_qc[
  is.na(cell_retention$cells_after_qc)
] <- 0

cell_retention$retention_percent <- with(
  cell_retention,
  ifelse(
    cells_before_qc > 0,
    100 * cells_after_qc / cells_before_qc,
    NA_real_
  )
)

write.csv(
  cell_retention,
  file = file.path(
    output_qc_dir,
    "QC_cell_retention_by_sample.csv"
  ),
  row.names = FALSE
)


# 13. Export and visualize post-QC metrics --------------------------------------

qc_summary_post <- summarise_qc(
  object = scRNA,
  stage = "Post-QC"
)

qc_summary <- bind_rows(
  qc_summary_pre,
  qc_summary_post
)

write.csv(
  qc_summary,
  file = file.path(
    output_qc_dir,
    "QC_metric_summary_by_sample.csv"
  ),
  row.names = FALSE
)

qc_plots_post <- lapply(
  qc_features,
  function(feature_name) {
    VlnPlot(
      object = scRNA,
      group.by = "orig.ident",
      features = feature_name,
      pt.size = 0,
      raster = FALSE
    ) +
      NoLegend() +
      ggtitle(
        paste0(
          feature_name,
          " (Post-QC)"
        )
      ) +
      theme(
        axis.text.x = element_text(
          angle = 45,
          hjust = 1
        )
      )
  }
)

ggsave(
  filename = file.path(
    output_qc_dir,
    "QC_metrics_after_filtering.pdf"
  ),
  plot = wrap_plots(qc_plots_post, nrow = 3),
  width = 16,
  height = 8,
  units = "in"
)


# 14. Plot nCount_RNA versus nFeature_RNA --------------------------------------

p_feature_scatter <- FeatureScatter(
  object = scRNA,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA"
)

ggsave(
  filename = file.path(
    output_qc_dir,
    "QC_nCount_vs_nFeature_after_filtering.pdf"
  ),
  plot = p_feature_scatter,
  width = 6,
  height = 5,
  units = "in"
)


# 15. Save the filtered object and analysis information ------------------------

saveRDS(
  scRNA,
  file = file.path(
    output_data_dir,
    "scRNA_after_QC.rds"
  )
)

qc_parameters <- c(
  paste0("Minimum detected genes per cell: > ", min_features),
  paste0("Maximum detected genes per cell: < ", max_features),
  paste0("Maximum UMI count per cell: < ", max_counts),
  paste0("Maximum mitochondrial transcript percentage: < ", max_percent_mt),
  "Ribosomal and haemoglobin transcript percentages were calculated for QC visualization but were not used as filtering criteria."
)

writeLines(
  qc_parameters,
  con = file.path(
    output_qc_dir,
    "QC_filtering_parameters.txt"
  )
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(
    output_qc_dir,
    "QC_sessionInfo.txt"
  )
)

# ==============================================================================
#  (3) Dimensionality reduction, batch integration and clustering
#
# Description:
# This script performs SCTransform normalization, PCA, UMAP and graph-based
# clustering on the QC-filtered Seurat object. Two parallel analysis branches
# are generated:
#   1) Harmony batch integration using the metadata column "batch"
#   2) No batch correction, using PCA directly
#
# Software used in the original analysis:
#   Seurat 4.2.1
#   dplyr 1.1.1
#   harmony 0.1.1
# ==============================================================================


# 1. Load packages -------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(harmony)
  library(ggplot2)
  library(patchwork)
})

set.seed(1234)


# 2. Define input and output paths ---------------------------------------------

input_file <- "DATA/OUTPUT/scRNA_after_QC.rds"

output_data_dir <- "DATA/OUTPUT"
output_figure_dir <- "FIGURE/STEP3_dimension_reduction"

dir.create(
  output_data_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  output_figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# 3. Define analysis parameters ------------------------------------------------

n_pcs <- 50
dims_use <- 1:30

harmony_max_iterations <- 20

harmony_resolutions <- c(
  0.1,
  0.3,
  0.5,
  0.7,
  1.0
)

pca_resolutions <- c(
  0.1,
  0.2,
  0.3,
  0.5,
  0.7,
  1.0
)

random_seed <- 1234


# 4. Load the QC-filtered Seurat object ----------------------------------------

scRNA_qc <- readRDS(input_file)

if (!inherits(scRNA_qc, "Seurat")) {
  stop("The input file must contain a Seurat object.")
}

if (!"RNA" %in% names(scRNA_qc@assays)) {
  stop("The Seurat object does not contain an RNA assay.")
}

required_metadata <- c(
  "orig.ident",
  "group",
  "batch"
)

missing_metadata <- setdiff(
  required_metadata,
  colnames(scRNA_qc@meta.data)
)

if (length(missing_metadata) > 0) {
  stop(
    paste0(
      "The following required metadata columns are missing: ",
      paste(missing_metadata, collapse = ", ")
    )
  )
}

if (anyNA(scRNA_qc$batch)) {
  stop("Missing values were detected in the batch metadata.")
}


# ==============================================================================
# Branch A: Harmony batch integration
# ==============================================================================


# 5. Copy the QC-filtered object -----------------------------------------------

scRNA_harmony <- scRNA_qc


# 6. SCTransform normalization -------------------------------------------------

scRNA_harmony <- SCTransform(
  object = scRNA_harmony,
  assay = "RNA",
  new.assay.name = "SCT",
  verbose = FALSE
)


# 7. Principal-component analysis ---------------------------------------------

scRNA_harmony <- RunPCA(
  object = scRNA_harmony,
  assay = "SCT",
  npcs = n_pcs,
  verbose = FALSE,
  seed.use = random_seed
)

p_elbow_harmony <- ElbowPlot(
  object = scRNA_harmony,
  ndims = n_pcs
)

ggsave(
  filename = file.path(
    output_figure_dir,
    "Harmony_branch_PCA_elbow_plot.pdf"
  ),
  plot = p_elbow_harmony,
  width = 5,
  height = 4,
  units = "in"
)


# 8. Harmony integration -------------------------------------------------------

scRNA_harmony <- RunHarmony(
  object = scRNA_harmony,
  group.by.vars = "batch",
  reduction = "pca",
  assay.use = "SCT",
  max.iter.harmony = harmony_max_iterations,
  reduction.save = "harmony",
  plot_convergence = FALSE,
  verbose = TRUE
)

if (!"harmony" %in% names(scRNA_harmony@reductions)) {
  stop("Harmony did not generate a reduction named 'harmony'.")
}


# 9. UMAP using Harmony embeddings ---------------------------------------------

scRNA_harmony <- RunUMAP(
  object = scRNA_harmony,
  reduction = "harmony",
  dims = dims_use,
  seed.use = random_seed,
  reduction.name = "umap",
  reduction.key = "UMAP_",
  verbose = FALSE
)


# 10. Shared-nearest-neighbour graph -------------------------------------------

scRNA_harmony <- FindNeighbors(
  object = scRNA_harmony,
  reduction = "harmony",
  dims = dims_use,
  verbose = FALSE
)


# 11. Clustering across the prespecified resolution grid -----------------------

for (resolution_value in harmony_resolutions) {

  scRNA_harmony <- FindClusters(
    object = scRNA_harmony,
    resolution = resolution_value,
    random.seed = random_seed,
    verbose = FALSE
  )
}


# 12. Harmony diagnostic UMAPs -------------------------------------------------

p_harmony_batch <- DimPlot(
  object = scRNA_harmony,
  reduction = "umap",
  group.by = "batch",
  pt.size = 0.1
) +
  ggtitle("Harmony-integrated UMAP by batch")

p_harmony_group <- DimPlot(
  object = scRNA_harmony,
  reduction = "umap",
  group.by = "group",
  pt.size = 0.1
) +
  ggtitle("Harmony-integrated UMAP by group")

ggsave(
  filename = file.path(
    output_figure_dir,
    "Harmony_branch_UMAP_diagnostics.pdf"
  ),
  plot = p_harmony_batch + p_harmony_group,
  width = 10,
  height = 4.5,
  units = "in"
)


# 13. Save the Harmony-integrated object ---------------------------------------

saveRDS(
  scRNA_harmony,
  file = file.path(
    output_data_dir,
    "scRNA_HarmonyIntegrated.rds"
  )
)


# ==============================================================================
# Branch B: no batch correction
# ==============================================================================


# 14. Copy the same QC-filtered input object -----------------------------------

scRNA_no_batch <- scRNA_qc


# 15. SCTransform normalization ------------------------------------------------

scRNA_no_batch <- SCTransform(
  object = scRNA_no_batch,
  assay = "RNA",
  new.assay.name = "SCT",
  verbose = FALSE
)


# 16. Principal-component analysis --------------------------------------------

scRNA_no_batch <- RunPCA(
  object = scRNA_no_batch,
  assay = "SCT",
  npcs = n_pcs,
  verbose = FALSE,
  seed.use = random_seed
)

p_elbow_no_batch <- ElbowPlot(
  object = scRNA_no_batch,
  ndims = n_pcs
)

ggsave(
  filename = file.path(
    output_figure_dir,
    "No_batch_correction_PCA_elbow_plot.pdf"
  ),
  plot = p_elbow_no_batch,
  width = 5,
  height = 4,
  units = "in"
)


# 17. UMAP using PCA embeddings ------------------------------------------------

scRNA_no_batch <- RunUMAP(
  object = scRNA_no_batch,
  reduction = "pca",
  dims = dims_use,
  seed.use = random_seed,
  reduction.name = "umap",
  reduction.key = "UMAP_",
  verbose = FALSE
)


# 18. Shared-nearest-neighbour graph -------------------------------------------

scRNA_no_batch <- FindNeighbors(
  object = scRNA_no_batch,
  reduction = "pca",
  dims = dims_use,
  verbose = FALSE
)


# 19. Clustering across the prespecified resolution grid -----------------------

for (resolution_value in pca_resolutions) {

  scRNA_no_batch <- FindClusters(
    object = scRNA_no_batch,
    resolution = resolution_value,
    random.seed = random_seed,
    verbose = FALSE
  )
}


# 20. No-correction diagnostic UMAPs -------------------------------------------

p_no_batch_batch <- DimPlot(
  object = scRNA_no_batch,
  reduction = "umap",
  group.by = "batch",
  pt.size = 0.1
) +
  ggtitle("Uncorrected UMAP by batch")

p_no_batch_group <- DimPlot(
  object = scRNA_no_batch,
  reduction = "umap",
  group.by = "group",
  pt.size = 0.1
) +
  ggtitle("Uncorrected UMAP by group")

ggsave(
  filename = file.path(
    output_figure_dir,
    "No_batch_correction_UMAP_diagnostics.pdf"
  ),
  plot = p_no_batch_batch + p_no_batch_group,
  width = 10,
  height = 4.5,
  units = "in"
)


# 21. Save the uncorrected object ----------------------------------------------

saveRDS(
  scRNA_no_batch,
  file = file.path(
    output_data_dir,
    "scRNA_NoBatchCorrection.rds"
  )
)


# 22. Export clustering metadata -----------------------------------------------

harmony_cluster_columns <- grep(
  pattern = "_snn_res\\.",
  x = colnames(scRNA_harmony@meta.data),
  value = TRUE
)

pca_cluster_columns <- grep(
  pattern = "_snn_res\\.",
  x = colnames(scRNA_no_batch@meta.data),
  value = TRUE
)

harmony_cluster_table <- scRNA_harmony@meta.data |>
  tibble::rownames_to_column("cell_barcode") |>
  dplyr::select(
    cell_barcode,
    orig.ident,
    group,
    batch,
    dplyr::all_of(harmony_cluster_columns)
  )

pca_cluster_table <- scRNA_no_batch@meta.data |>
  tibble::rownames_to_column("cell_barcode") |>
  dplyr::select(
    cell_barcode,
    orig.ident,
    group,
    batch,
    dplyr::all_of(pca_cluster_columns)
  )

write.csv(
  harmony_cluster_table,
  file = file.path(
    output_data_dir,
    "Harmony_cluster_assignments.csv"
  ),
  row.names = FALSE
)

write.csv(
  pca_cluster_table,
  file = file.path(
    output_data_dir,
    "No_batch_correction_cluster_assignments.csv"
  ),
  row.names = FALSE
)


# 23. Save analysis parameters -------------------------------------------------

analysis_parameters <- c(
  paste0("Input file: ", input_file),
  "Normalization: SCTransform using the RNA assay",
  paste0("Number of computed principal components: ", n_pcs),
  paste0(
    "Dimensions used for UMAP and graph construction: ",
    min(dims_use),
    "-",
    max(dims_use)
  ),
  "Harmony integration variable: batch",
  paste0(
    "Harmony maximum iterations: ",
    harmony_max_iterations
  ),
  paste0(
    "Harmony clustering resolutions: ",
    paste(harmony_resolutions, collapse = ", ")
  ),
  paste0(
    "No-correction clustering resolutions: ",
    paste(pca_resolutions, collapse = ", ")
  ),
  paste0("Random seed: ", random_seed)
)

writeLines(
  analysis_parameters,
  con = file.path(
    output_data_dir,
    "Dimension_reduction_clustering_parameters.txt"
  )
)


# 24. Save software information ------------------------------------------------

writeLines(
  capture.output(sessionInfo()),
  con = file.path(
    output_data_dir,
    "Dimension_reduction_clustering_sessionInfo.txt"
  )
)

-------------------------------------------------------------------------------

### ----(4)-Marker Gene Identification and DEG Analysis --------------

#==============================================================================
  # Cluster marker identification, differential expression and cell annotation
  #
  # Description:
  # This script identifies marker genes for graph-based clusters, performs a
  # assigns biological cell-type labels, and identifies marker genes for the
  # annotated cell types.
  # ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(future)
})

set.seed(1234)

input_file <- "DATA/OUTPUT/scRNA_HarmonyIntegrated.rds"
output_data_dir <- "DATA/OUTPUT"
output_table_dir <- "TABLE/DEG"

dir.create(output_data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_table_dir, recursive = TRUE, showWarnings = FALSE)

cluster_column <- "SCT_snn_res.0.5"

marker_assay <- "SCT"
marker_slot <- "data"

deg_assay <- "RNA"

# Keep "counts" only if this was the slot used to generate the reported result.
# For a conventional Wilcoxon-based DEG analysis, the normalized "data" slot is
# generally preferable.
deg_slot <- "counts"

ident_1 <- "0"
ident_2 <- "6"
random_seed <- 1234

workers <- min(
  10L,
  max(
    1L,
    parallel::detectCores(logical = TRUE) - 1L
  )
)

options(future.globals.maxSize = 20 * 1024^3)

future::plan(
  future::multisession,
  workers = workers
)

scRNA <- readRDS(input_file)

if (!inherits(scRNA, "Seurat")) {
  stop("The input file must contain a Seurat object.")
}

if (!cluster_column %in% colnames(scRNA@meta.data)) {
  stop(
    paste0(
      "The clustering column '",
      cluster_column,
      "' was not found in the Seurat metadata."
    )
  )
}

if (!marker_assay %in% names(scRNA@assays)) {
  stop(
    paste0(
      "The marker assay '",
      marker_assay,
      "' was not found."
    )
  )
}

if (!deg_assay %in% names(scRNA@assays)) {
  stop(
    paste0(
      "The DEG assay '",
      deg_assay,
      "' was not found."
    )
  )
}

cluster_ids <- sort(
  unique(
    as.character(
      scRNA@meta.data[[cluster_column]]
    )
  )
)

message(
  "Clusters present in ",
  cluster_column,
  ": ",
  paste(cluster_ids, collapse = ", ")
)

# 1. Marker genes for graph-based clusters -------------------------------------

Idents(scRNA) <- cluster_column

cluster_markers <- FindAllMarkers(
  object = scRNA,
  assay = marker_assay,
  slot = marker_slot,
  test.use = "wilcox",
  logfc.threshold = 0.25,
  min.pct = 0.10,
  only.pos = FALSE,
  random.seed = random_seed,
  verbose = TRUE
)

write.csv(
  cluster_markers,
  file = file.path(
    output_table_dir,
    "Markers_all_clusters_resolution_0.5.csv"
  ),
  row.names = FALSE
)

# 2. DEG between clusters 0 and 6 ----------------------------------------------

missing_deg_clusters <- setdiff(
  c(ident_1, ident_2),
  cluster_ids
)

if (length(missing_deg_clusters) > 0) {
  stop(
    paste0(
      "The following DEG comparison cluster(s) are absent from ",
      cluster_column,
      ": ",
      paste(missing_deg_clusters, collapse = ", ")
    )
  )
}

deg_cluster_0_vs_6 <- FindMarkers(
  object = scRNA,
  ident.1 = ident_1,
  ident.2 = ident_2,
  group.by = cluster_column,
  assay = deg_assay,
  slot = deg_slot,
  test.use = "wilcox",
  logfc.threshold = 0,
  min.pct = 0,
  random.seed = random_seed,
  verbose = TRUE
)

deg_cluster_0_vs_6 <- deg_cluster_0_vs_6 |>
  tibble::rownames_to_column("gene")

write.csv(
  deg_cluster_0_vs_6,
  file = file.path(
    output_table_dir,
    "Cluster_0_vs_6_DEG.csv"
  ),
  row.names = FALSE
)

# 3. Cell-type annotation -------------------------------------------------------

# Replace or extend this mapping so that every cluster has exactly one label.
# The previous label "NDM" has been updated to "Numb+ microglia".

cluster_to_celltype <- c(
  "0" = "Numb+ microglia",
  "1" = "Mg1",
  "2" = "Mg2"
)

unannotated_clusters <- setdiff(
  cluster_ids,
  names(cluster_to_celltype)
)

if (length(unannotated_clusters) > 0) {
  annotation_template <- data.frame(
    cluster = cluster_ids,
    celltype = unname(
      cluster_to_celltype[cluster_ids]
    ),
    stringsAsFactors = FALSE
  )

  write.csv(
    annotation_template,
    file = file.path(
      output_table_dir,
      "Cluster_annotation_template.csv"
    ),
    row.names = FALSE,
    na = ""
  )

  stop(
    paste0(
      "Cell-type labels are missing for cluster(s): ",
      paste(unannotated_clusters, collapse = ", "),
      ". Complete cluster_to_celltype before continuing. ",
      "A template has been written to TABLE/DEG/Cluster_annotation_template.csv."
    )
  )
}

scRNA$celltype <- unname(
  cluster_to_celltype[
    as.character(
      scRNA@meta.data[[cluster_column]]
    )
  ]
)

if (anyNA(scRNA$celltype)) {
  stop("Missing cell-type labels were generated during annotation.")
}

celltype_levels <- unique(
  unname(
    cluster_to_celltype[cluster_ids]
  )
)

scRNA$celltype <- factor(
  scRNA$celltype,
  levels = celltype_levels
)

print(
  table(
    cluster = scRNA@meta.data[[cluster_column]],
    celltype = scRNA$celltype
  )
)

# 4. Marker genes for annotated cell types -------------------------------------

Idents(scRNA) <- "celltype"

celltype_markers <- FindAllMarkers(
  object = scRNA,
  assay = marker_assay,
  slot = marker_slot,
  test.use = "wilcox",
  logfc.threshold = 0.25,
  min.pct = 0.10,
  only.pos = FALSE,
  random.seed = random_seed,
  verbose = TRUE
)

write.csv(
  celltype_markers,
  file = file.path(
    output_table_dir,
    "Markers_all_annotated_celltypes.csv"
  ),
  row.names = FALSE
)

# 5. Export annotation metadata -------------------------------------------------

annotation_metadata <- scRNA@meta.data |>
  tibble::rownames_to_column("cell_barcode") |>
  dplyr::select(
    cell_barcode,
    orig.ident,
    group,
    batch,
    dplyr::all_of(cluster_column),
    celltype
  )

write.csv(
  annotation_metadata,
  file = file.path(
    output_table_dir,
    "Cell_cluster_and_celltype_annotations.csv"
  ),
  row.names = FALSE
)

# 6. Save annotated object ------------------------------------------------------

saveRDS(
  scRNA,
  file = file.path(
    output_data_dir,
    "scRNA_HarmonyIntegrated_annotated.rds"
  )
)

# 7. Save parameters and software information ---------------------------------

analysis_parameters <- c(
  paste0("Input file: ", input_file),
  paste0("Clustering column: ", cluster_column),
  paste0("Marker assay: ", marker_assay),
  paste0("Marker slot: ", marker_slot),
  "Marker test: Wilcoxon rank-sum test",
  "FindAllMarkers logfc.threshold: 0.25",
  "FindAllMarkers min.pct: 0.10",
  "FindAllMarkers only.pos: FALSE",
  paste0("DEG comparison: cluster ", ident_1, " versus cluster ", ident_2),
  paste0("DEG assay: ", deg_assay),
  paste0("DEG slot: ", deg_slot),
  "FindMarkers logfc.threshold: 0",
  "FindMarkers min.pct: 0",
  paste0("Parallel workers: ", workers),
  "future.globals.maxSize: 20 GiB",
  paste0("Random seed: ", random_seed)
)

writeLines(
  analysis_parameters,
  con = file.path(
    output_table_dir,
    "Marker_DEG_annotation_parameters.txt"
  )
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(
    output_table_dir,
    "Marker_DEG_annotation_sessionInfo.txt"
  )
)

future::plan(future::sequential)




###------------- 4.2: DEG Between Specific Clusters-----------------------------


options(future.globals.maxSize = 20000 * 1024^3)
plan(multisession, workers = 20)
# Compare between cluster "0" and "6" under clustering result SCT_snn_res.0.5
deg <- FindMarkers(
  object = scRNA,
  ident.1 = "0",
  ident.2 = "6",
  group.by = "SCT_snn_res.0.5",
  logfc.threshold = 0,
  min.pct = 0,
  assay = "RNA",
  slot = "counts")

# Save DEG result
write.csv(deg, file = "Cluster_0_vs_6_DEG.csv", row.names = TRUE)


--------------------------------------------------------------------------------

###亚群命名#   5. cluster annotation     ###########

metadata <- scRNA@meta.data
metadata$celltype <- recode(metadata$SCT_snn_res.0.5,
                            `0` = "NDM",`1` = "Mg1",`2` = "Mg2")
scRNA@meta.data <- metadata
scRNA$celltype<- factor(scRNA$celltype,levels =
                          c("NDM","Mg1","Mg2"))

--------------------------------------------------------------------------------
###  ------(5).Visualization ----------------------------------------------------
# ==============================================================================
# Visualization of clustering, cell states and marker-gene expression
#
# Description:
# This script generates the cluster tree, UMAPs, cell-composition plot,
# gene-expression plots, per-cell expression heatmap and marker DotPlot used in
# the manuscript.
#
# Software used in the original analysis:
#   Seurat
#   clustree 0.5.0
#   SCP 0.4.2
#   ComplexHeatmap 2.14.0
# ==============================================================================


# 1. Load packages -------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(clustree)
  library(SCP)
  library(dplyr)
  library(readxl)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(RColorBrewer)
})


# 2. Define paths ---------------------------------------------------------------

input_file <- "DATA/OUTPUT/scRNA_HarmonyIntegrated_annotated.rds"

heatmap_gene_file <- "DATA/INPUT/p3_all_cell_genes.xlsx"
dotplot_gene_file <- "DATA/INPUT/dotplot_genes.xlsx"

output_figure_dir <- "FIGURE"
output_umap_dir <- file.path(output_figure_dir, "UMAP")
output_data_dir <- "DATA/OUTPUT"

dir.create(output_figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_umap_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_data_dir, recursive = TRUE, showWarnings = FALSE)


# 3. Load annotated Seurat object ----------------------------------------------

scRNA <- readRDS(input_file)

if (!inherits(scRNA, "Seurat")) {
  stop("The input file must contain a Seurat object.")
}

required_metadata <- c("celltype", "group")

missing_metadata <- setdiff(
  required_metadata,
  colnames(scRNA@meta.data)
)

if (length(missing_metadata) > 0) {
  stop(
    paste0(
      "Missing metadata column(s): ",
      paste(missing_metadata, collapse = ", ")
    )
  )
}

if (!"umap" %in% names(scRNA@reductions)) {
  stop("The Seurat object does not contain a UMAP reduction named 'umap'.")
}

if (!"SCT" %in% names(scRNA@assays)) {
  stop("The Seurat object does not contain an SCT assay.")
}


# 4. Define final cell-state order and colors ----------------------------------

celltype_levels <- levels(
  droplevels(
    factor(scRNA$celltype)
  )
)

if (length(celltype_levels) == 0) {
  celltype_levels <- unique(
    as.character(scRNA$celltype)
  )
}

scRNA$celltype <- factor(
  scRNA$celltype,
  levels = celltype_levels
)

# Replace this vector with the exact color mapping used in the final figures if
# a fixed manuscript-wide palette was applied.
base_celltype_colors <- c(
  "#E5D2DD",
  "#53A85F",
  "#F3B1A0",
  "#FFDD44",
  "#9467BD",
  "#E377C2",
  "#D62728",
  "#8C564B",
  "#2CA02C",
  "#23452F",
  "#1F77B4",
  "#17BECF",
  "#BCBD22",
  "#7F7F7F",
  "#FF7F0E"
)

if (length(celltype_levels) > length(base_celltype_colors)) {
  stop("More cell types were found than colors available in the palette.")
}

celltype_colors <- setNames(
  base_celltype_colors[
    seq_along(celltype_levels)
  ],
  celltype_levels
)


# ==============================================================================
# 5. Cluster-tree visualization
# ==============================================================================

cluster_columns <- grep(
  pattern = "^SCT_snn_res\\.",
  x = colnames(scRNA@meta.data),
  value = TRUE
)

if (length(cluster_columns) < 2) {
  warning(
    "Fewer than two SCT_snn_res.* columns were found; the clustree plot was skipped."
  )
} else {

  p_clustree <- clustree(
    scRNA@meta.data,
    prefix = "SCT_snn_res."
  ) +
    guides(
      edge_colour = "none",
      edge_alpha = "none"
    ) +
    scale_color_brewer(
      palette = "Set1"
    ) +
    scale_edge_color_continuous(
      low = "blue",
      high = "red"
    ) +
    theme(
      legend.position = "bottom"
    )

  ggsave(
    filename = file.path(
      output_figure_dir,
      "ClusterTree.tiff"
    ),
    plot = p_clustree,
    width = 8,
    height = 8,
    units = "in",
    dpi = 600,
    compression = "lzw"
  )

  ggsave(
    filename = file.path(
      output_figure_dir,
      "ClusterTree.pdf"
    ),
    plot = p_clustree,
    width = 8,
    height = 8,
    units = "in"
  )
}


# ==============================================================================
# 6. UMAP by cell type and experimental group
# ==============================================================================

p_umap_celltype <- SCP::CellDimPlot(
  srt = scRNA,
  group.by = "celltype",
  reduction = "umap",
  theme_use = "theme_blank",
  label = TRUE,
  label_insitu = TRUE
)

ggsave(
  filename = file.path(
    output_figure_dir,
    "UMAP_celltype.tiff"
  ),
  plot = p_umap_celltype,
  width = 6,
  height = 4,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  filename = file.path(
    output_figure_dir,
    "UMAP_celltype.pdf"
  ),
  plot = p_umap_celltype,
  width = 4,
  height = 4,
  units = "in"
)

p_umap_group <- SCP::CellDimPlot(
  srt = scRNA,
  group.by = "celltype",
  reduction = "umap",
  theme_use = "theme_blank",
  label = FALSE,
  label_insitu = FALSE,
  show_stat = FALSE,
  split.by = "group"
)

ggsave(
  filename = file.path(
    output_figure_dir,
    "UMAP_split_by_group.tiff"
  ),
  plot = p_umap_group,
  width = 15,
  height = 12,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  filename = file.path(
    output_figure_dir,
    "UMAP_split_by_group.pdf"
  ),
  plot = p_umap_group,
  width = 15,
  height = 12,
  units = "in"
)


# ==============================================================================
# 7. Feature plot
# ==============================================================================

feature_genes <- c("Gpnmb")

missing_feature_genes <- setdiff(
  feature_genes,
  rownames(scRNA)
)

if (length(missing_feature_genes) > 0) {
  warning(
    paste0(
      "Feature gene(s) not found and skipped: ",
      paste(missing_feature_genes, collapse = ", ")
    )
  )
}

feature_genes <- intersect(
  feature_genes,
  rownames(scRNA)
)

for (gene_name in feature_genes) {

  p_feature <- SCP::FeatureDimPlot(
    object = scRNA,
    features = gene_name,
    reduction = "umap",
    cells.highlight = TRUE,
    theme_use = "theme_blank",
    show_stat = FALSE,
    legend.position = "none"
  )

  ggsave(
    filename = file.path(
      output_umap_dir,
      paste0(gene_name, ".tiff")
    ),
    plot = p_feature,
    width = 2,
    height = 2,
    units = "in",
    dpi = 600,
    compression = "lzw"
  )

  ggsave(
    filename = file.path(
      output_umap_dir,
      paste0(gene_name, ".pdf")
    ),
    plot = p_feature,
    width = 2,
    height = 2,
    units = "in"
  )
}


# ==============================================================================
# 8. Cell-state proportion plot
# ==============================================================================

p_cell_ratio <- SCP::CellStatPlot(
  srt = scRNA,
  stat.by = "celltype",
  group.by = "group",
  plot_type = "trend"
)

ggsave(
  filename = file.path(
    output_figure_dir,
    "Cell_state_proportions.pdf"
  ),
  plot = p_cell_ratio,
  width = 4,
  height = 3,
  units = "in"
)

ggsave(
  filename = file.path(
    output_figure_dir,
    "Cell_state_proportions.tiff"
  ),
  plot = p_cell_ratio,
  width = 4,
  height = 3,
  units = "in",
  dpi = 600,
  compression = "lzw"
)


# ==============================================================================
# 9. Gene-expression boxplots
# ==============================================================================

boxplot_genes <- c(
  "Pgk1",
  "Pgam1",
  "Pkm",
  "Ldha",
  "Aif1",
  "P2ry12"
)

boxplot_genes <- intersect(
  boxplot_genes,
  rownames(scRNA)
)

if (length(boxplot_genes) == 0) {
  warning("None of the requested boxplot genes were found.")
} else {

  p_gene_boxplot <- SCP::FeatureStatPlot(
    srt = scRNA,
    stat.by = boxplot_genes,
    fill.by = "group",
    plot_type = "box",
    group.by = "celltype",
    bg.by = "celltype",
    stack = TRUE,
    flip = FALSE
  )

  ggsave(
    filename = file.path(
      output_figure_dir,
      "Gene_expression_boxplots.tiff"
    ),
    plot = p_gene_boxplot,
    width = 8,
    height = 4,
    units = "in",
    dpi = 600,
    compression = "lzw"
  )

  ggsave(
    filename = file.path(
      output_figure_dir,
      "Gene_expression_boxplots.pdf"
    ),
    plot = p_gene_boxplot,
    width = 8,
    height = 4,
    units = "in"
  )
}


# ==============================================================================
# 10. Per-cell gene-expression heatmap
# ==============================================================================

if (!file.exists(heatmap_gene_file)) {
  warning(
    paste0(
      "Heatmap gene file not found; heatmap skipped: ",
      heatmap_gene_file
    )
  )
} else {

  gene_data <- readxl::read_excel(
    heatmap_gene_file,
    sheet = 1
  )

  if (ncol(gene_data) < 2) {
    stop("The heatmap gene file must contain at least two columns.")
  }

  colnames(gene_data)[1:2] <- c(
    "cluster",
    "gene"
  )

  gene_data <- gene_data |>
    dplyr::select(cluster, gene) |>
    dplyr::filter(
      !is.na(gene),
      gene != ""
    ) |>
    dplyr::distinct(gene, .keep_all = TRUE)

  expression_matrix <- GetAssayData(
    object = scRNA,
    assay = "SCT",
    slot = "data"
  )

  valid_gene_data <- gene_data |>
    dplyr::filter(
      gene %in% rownames(expression_matrix)
    )

  if (nrow(valid_gene_data) == 0) {
    stop("None of the requested heatmap genes were found in the SCT assay.")
  }

  heatmap_data <- as.matrix(
    expression_matrix[
      valid_gene_data$gene,
      ,
      drop = FALSE
    ]
  )

  rownames(heatmap_data) <- valid_gene_data$gene

  heatmap_celltypes <- factor(
    as.character(
      scRNA$celltype[
        colnames(heatmap_data)
      ]
    ),
    levels = celltype_levels
  )

  row_clusters <- factor(
    valid_gene_data$cluster,
    levels = unique(
      valid_gene_data$cluster
    )
  )

  min_max_scale <- function(x) {
    x_min <- min(x, na.rm = TRUE)
    x_max <- max(x, na.rm = TRUE)

    if (!is.finite(x_min) || !is.finite(x_max)) {
      return(
        rep(
          NA_real_,
          length(x)
        )
      )
    }

    if (x_max == x_min) {
      return(
        rep(
          0,
          length(x)
        )
      )
    }

    (x - x_min) / (x_max - x_min)
  }

  scaled_heatmap_data <- t(
    apply(
      heatmap_data,
      1,
      min_max_scale
    )
  )

  rownames(scaled_heatmap_data) <- rownames(heatmap_data)
  colnames(scaled_heatmap_data) <- colnames(heatmap_data)

  row_cluster_levels <- levels(row_clusters)

  row_cluster_palette <- setNames(
    grDevices::hcl.colors(
      n = length(row_cluster_levels),
      palette = "Dark 3"
    ),
    row_cluster_levels
  )

  top_annotation <- columnAnnotation(
    celltype = heatmap_celltypes,
    col = list(
      celltype = celltype_colors
    ),
    show_annotation_name = FALSE
  )

  left_annotation <- rowAnnotation(
    cluster = row_clusters,
    col = list(
      cluster = row_cluster_palette
    ),
    show_annotation_name = FALSE
  )

  heatmap_object <- Heatmap(
    scaled_heatmap_data,
    name = "Expression",
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_column_names = FALSE,
    show_row_names = FALSE,
    column_split = heatmap_celltypes,
    row_split = row_clusters,
    top_annotation = top_annotation,
    left_annotation = left_annotation,
    col = circlize::colorRamp2(
      c(0, 0.5, 1),
      c(
        "#4978B3",
        "white",
        "#FF3333"
      )
    ),
    heatmap_legend_param = list(
      at = seq(0, 1, 0.2),
      labels = seq(0, 1, 0.2),
      title = "Expression",
      title_position = "leftcenter-rot"
    ),
    border = TRUE,
    use_raster = TRUE,
    raster_quality = 2,
    column_gap = unit(1, "mm"),
    row_gap = unit(1, "mm")
  )

  pdf(
    file = file.path(
      output_figure_dir,
      "GENE_EXPRESSION.pdf"
    ),
    width = 6.5,
    height = 6,
    useDingbats = FALSE
  )

  draw(
    heatmap_object,
    merge_legends = TRUE
  )

  dev.off()

  write.csv(
    scaled_heatmap_data,
    file = file.path(
      output_data_dir,
      "Heatmap_scaled_expression_matrix.csv"
    ),
    row.names = TRUE
  )

  write.csv(
    valid_gene_data,
    file = file.path(
      output_data_dir,
      "Heatmap_gene_annotations.csv"
    ),
    row.names = FALSE
  )
}


# ==============================================================================
# 11. Marker-gene DotPlot
# ==============================================================================

if (!file.exists(dotplot_gene_file)) {
  warning(
    paste0(
      "DotPlot gene file not found; DotPlot skipped: ",
      dotplot_gene_file
    )
  )
} else {

  marker_table <- readxl::read_excel(
    dotplot_gene_file,
    sheet = 1
  )

  if (!"gene" %in% colnames(marker_table)) {
    stop("The DotPlot gene file must contain a column named 'gene'.")
  }

  marker_genes <- unique(
    as.character(marker_table$gene)
  )

  marker_genes <- marker_genes[
    !is.na(marker_genes) &
      marker_genes != ""
  ]

  missing_marker_genes <- setdiff(
    marker_genes,
    rownames(scRNA)
  )

  if (length(missing_marker_genes) > 0) {
    warning(
      paste0(
        "DotPlot gene(s) not found and skipped: ",
        paste(missing_marker_genes, collapse = ", ")
      )
    )
  }

  marker_genes <- intersect(
    marker_genes,
    rownames(scRNA)
  )

  if (length(marker_genes) == 0) {
    stop("None of the requested DotPlot genes were found.")
  }

  Idents(scRNA) <- "celltype"

  dotplot_base <- DotPlot(
    object = scRNA,
    features = marker_genes,
    assay = "SCT"
  )

  dot_data <- dotplot_base$data

  write.csv(
    dot_data,
    file = file.path(
      output_data_dir,
      "DotPlot_underlying_data.csv"
    ),
    row.names = FALSE
  )

  p_dotplot <- ggplot(
    dot_data,
    aes(
      x = features.plot,
      y = id,
      size = pct.exp,
      fill = avg.exp.scaled
    )
  ) +
    geom_point(
      shape = 21,
      colour = "black",
      stroke = 0.5
    ) +
    guides(
      size = guide_legend(
        override.aes = list(
          shape = 21,
          colour = "black",
          fill = NA
        )
      )
    ) +
    scale_fill_gradientn(
      colours = c(
        "#5749A0",
        "#0F7AB0",
        "#00BBB1",
        "#BEF0B0",
        "#FDF4AF",
        "#F9B64B",
        "#EC840E",
        "#CA443D",
        "#A51A49"
      )
    ) +
    theme(
      panel.background = element_blank(),
      panel.border = element_rect(
        fill = NA
      ),
      panel.grid.major.x = element_line(
        colour = "grey80"
      ),
      panel.grid.major.y = element_line(
        colour = "grey80"
      ),
      axis.title = element_blank(),
      axis.text.y = element_text(
        colour = "black",
        size = 12
      ),
      axis.text.x = element_text(
        colour = "black",
        size = 12,
        angle = 90,
        hjust = 1,
        vjust = 0.5
      )
    )

  ggsave(
    filename = file.path(
      output_figure_dir,
      "DOTPLOT_NOT_SLICED.pdf"
    ),
    plot = p_dotplot,
    width = 6,
    height = 2,
    units = "in"
  )
}


# 12. Save analysis information ------------------------------------------------

analysis_parameters <- c(
  paste0("Input file: ", input_file),
  paste0(
    "Cell-type order: ",
    paste(celltype_levels, collapse = ", ")
  ),
  "Heatmap assay: SCT",
  "Heatmap slot: data",
  "Heatmap scaling: row-wise min-max scaling to 0-1",
  "DotPlot assay: SCT"
)

writeLines(
  analysis_parameters,
  con = file.path(
    output_data_dir,
    "Visualization_parameters.txt"
  )
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(
    output_data_dir,
    "Visualization_sessionInfo.txt"
  )
)

  #==============================================================================
  # (6) Gene Ontology enrichment analysis and bubble-plot visualization
  #
  # Description:
  # This script filters positively regulated genes for each cluster, maps mouse
  # gene symbols to Entrez identifiers, performs GO Biological Process enrichment
  # using compareCluster, exports the complete enrichment results and generates
  # the GO bubble plot used for visualization.
  #
  # Software used in the original analysis:
  #   clusterProfiler 4.2.2
  #   org.Mm.eg.db 3.14.0
  #   GOplot 1.0.2
  # ==============================================================================


# 1. Load packages -------------------------------------------------------------

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(dplyr)
  library(openxlsx)
  library(ggplot2)
  library(cowplot)
  library(aplot)
})


# 2. Define paths ---------------------------------------------------------------

deg_file <- "DATA/INPUT/module_gene.xlsx"

# Optional file containing two columns: Description and Annotation.
# When absent, the bubble plot is saved without the left annotation panel.
term_annotation_file <- "DATA/INPUT/GO_term_annotations.xlsx"

output_dir <- "go_enrichment"
output_figure_dir <- "FIGURE"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_figure_dir, recursive = TRUE, showWarnings = FALSE)


# 3. Define analysis parameters ------------------------------------------------

log2fc_cutoff <- 0
adjusted_p_cutoff <- 0.01

go_ontology <- "BP"
p_adjust_method <- "BH"
enrichment_p_cutoff <- 0.05
enrichment_q_cutoff <- 0.20
min_gene_set_size <- 10
max_gene_set_size <- 500


# 4. Read DEG table -------------------------------------------------------------

deg_table <- openxlsx::read.xlsx(
  deg_file
)

required_columns <- c(
  "cluster",
  "gene",
  "avg_log2FC",
  "p_val_adj"
)

missing_columns <- setdiff(
  required_columns,
  colnames(deg_table)
)

if (length(missing_columns) > 0) {
  stop(
    paste0(
      "The DEG table is missing: ",
      paste(missing_columns, collapse = ", ")
    )
  )
}


# 5. Select significant positively regulated genes -----------------------------

markers <- deg_table |>
  dplyr::filter(
    !is.na(cluster),
    !is.na(gene),
    avg_log2FC > log2fc_cutoff,
    p_val_adj < adjusted_p_cutoff
  ) |>
  dplyr::distinct(
    cluster,
    gene,
    .keep_all = TRUE
  )

if (nrow(markers) == 0) {
  stop("No genes passed the DEG filtering criteria.")
}

write.csv(
  markers,
  file = file.path(
    output_dir,
    "GO_input_genes.csv"
  ),
  row.names = FALSE
)


# 6. Map mouse gene symbols to Entrez identifiers ------------------------------

gene_mapping <- clusterProfiler::bitr(
  unique(markers$gene),
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)

markers_mapped <- markers |>
  dplyr::inner_join(
    gene_mapping,
    by = c(
      "gene" = "SYMBOL"
    )
  )

if (nrow(markers_mapped) == 0) {
  stop("None of the selected genes could be mapped to Entrez identifiers.")
}

write.csv(
  markers_mapped,
  file = file.path(
    output_dir,
    "GO_input_genes_with_Entrez_IDs.csv"
  ),
  row.names = FALSE
)

mapping_summary <- data.frame(
  total_selected_symbols = length(
    unique(markers$gene)
  ),
  mapped_symbols = length(
    unique(markers_mapped$gene)
  ),
  unmapped_symbols = length(
    setdiff(
      unique(markers$gene),
      unique(markers_mapped$gene)
    )
  )
)

write.csv(
  mapping_summary,
  file = file.path(
    output_dir,
    "GO_gene_mapping_summary.csv"
  ),
  row.names = FALSE
)


# 7. GO Biological Process enrichment by cluster ------------------------------

ego <- compareCluster(
  ENTREZID ~ cluster,
  data = markers_mapped,
  fun = "enrichGO",
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = go_ontology,
  pAdjustMethod = p_adjust_method,
  pvalueCutoff = enrichment_p_cutoff,
  qvalueCutoff = enrichment_q_cutoff,
  minGSSize = min_gene_set_size,
  maxGSSize = max_gene_set_size,
  readable = FALSE
)

ego_readable <- setReadable(
  ego,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID"
)

ego_df <- as.data.frame(
  ego_readable
)

if (nrow(ego_df) == 0) {
  stop("GO enrichment returned no terms.")
}

write.csv(
  ego_df,
  file = file.path(
    output_dir,
    "Enrichment_GO_complete.csv"
  ),
  row.names = FALSE
)

saveRDS(
  ego_readable,
  file = file.path(
    output_dir,
    "GO_compareCluster_result.rds"
  )
)


# 8. Prepare bubble-plot data ---------------------------------------------------

plot_data <- ego_df |>
  dplyr::filter(
    !is.na(p.adjust)
  ) |>
  dplyr::mutate(
    adjust_group = cut(
      p.adjust,
      breaks = c(
        -Inf,
        0.0001,
        0.001,
        0.01,
        0.05,
        0.1,
        Inf
      ),
      labels = c(
        "<0.0001",
        "<0.001",
        "<0.01",
        "<0.05",
        "<0.1",
        ">0.1"
      ),
      right = FALSE
    )
  )

adjust_group_levels <- c(
  "<0.0001",
  "<0.001",
  "<0.01",
  "<0.05",
  "<0.1",
  ">0.1"
)

plot_data$adjust_group <- factor(
  plot_data$adjust_group,
  levels = adjust_group_levels
)

# Preserve the cluster order from the input table. Replace this with an explicit
# final order if a different order was used in the published figure.
cluster_order <- unique(
  as.character(markers$cluster)
)

plot_data$Cluster <- factor(
  as.character(plot_data$Cluster),
  levels = cluster_order
)

# Preserve the order of GO terms as they appear in the enrichment output.
description_order <- unique(
  as.character(plot_data$Description)
)

plot_data$Description <- factor(
  plot_data$Description,
  levels = rev(description_order)
)

write.csv(
  plot_data,
  file = file.path(
    output_dir,
    "GO_bubble_plot_data.csv"
  ),
  row.names = FALSE
)


# 9. Generate the main bubble plot ---------------------------------------------

bubble_plot <- ggplot(
  plot_data,
  aes(
    x = Cluster,
    y = Description
  )
) +
  geom_vline(
    xintercept = seq_along(cluster_order),
    colour = "#D3D3D3"
  ) +
  geom_hline(
    yintercept = seq_along(description_order),
    colour = "#E8E8E8"
  ) +
  geom_point(
    aes(
      colour = adjust_group,
      size = Count
    ),
    shape = 19
  ) +
  scale_colour_manual(
    values = c(
      "<0.0001" = "#67000D",
      "<0.001" = "#EF3B2C",
      "<0.01" = "#FB6A4A",
      "<0.05" = "#FC9272",
      "<0.1" = "#FEE0D2",
      ">0.1" = "#FFF5F0"
    ),
    drop = FALSE
  ) +
  cowplot::theme_cowplot() +
  theme(
    panel.grid.major = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 1
    ),
    plot.title = element_text(
      hjust = 0.5
    ),
    legend.direction = "vertical",
    legend.position = "right"
  ) +
  labs(
    x = NULL,
    y = NULL,
    colour = "Adjusted P",
    size = "Gene count"
  ) +
  guides(
    colour = guide_legend(
      override.aes = list(
        size = 4
      ),
      order = 2
    ),
    size = guide_legend(
      order = 1
    )
  ) +
  scale_y_discrete(
    position = "right"
  )


# 10. Optional left annotation panel -------------------------------------------

final_plot <- bubble_plot

if (file.exists(term_annotation_file)) {

  term_annotation <- openxlsx::read.xlsx(
    term_annotation_file
  )

  required_annotation_columns <- c(
    "Description",
    "Annotation"
  )

  missing_annotation_columns <- setdiff(
    required_annotation_columns,
    colnames(term_annotation)
  )

  if (length(missing_annotation_columns) > 0) {
    stop(
      paste0(
        "The GO-term annotation file is missing: ",
        paste(missing_annotation_columns, collapse = ", ")
      )
    )
  }

  term_annotation <- term_annotation |>
    dplyr::filter(
      Description %in% as.character(
        plot_data$Description
      )
    ) |>
    dplyr::distinct(
      Description,
      .keep_all = TRUE
    )

  term_annotation$Description <- factor(
    term_annotation$Description,
    levels = levels(
      plot_data$Description
    )
  )

  annotation_levels <- unique(
    as.character(
      term_annotation$Annotation
    )
  )

  annotation_colors <- setNames(
    grDevices::hcl.colors(
      n = length(annotation_levels),
      palette = "Dark 3"
    ),
    annotation_levels
  )

  left_annotation_plot <- term_annotation |>
    dplyr::mutate(
      panel = ""
    ) |>
    ggplot(
      aes(
        x = panel,
        y = Description,
        fill = Annotation
      )
    ) +
    geom_tile() +
    scale_fill_manual(
      values = annotation_colors
    ) +
    scale_y_discrete(
      position = "right",
      limits = levels(
        plot_data$Description
      )
    ) +
    theme_minimal() +
    labs(
      x = NULL,
      y = NULL,
      fill = "Annotations"
    ) +
    theme(
      axis.text.y = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      legend.position = "top",
      legend.direction = "horizontal"
    )

  final_plot <- bubble_plot |>
    aplot::insert_left(
      left_annotation_plot,
      width = 0.08
    )
}


# 11. Save the final plot -------------------------------------------------------

ggsave(
  filename = file.path(
    output_figure_dir,
    "GO_Bubble.pdf"
  ),
  plot = final_plot,
  width = 18,
  height = 8,
  units = "in"
)

ggsave(
  filename = file.path(
    output_figure_dir,
    "GO_Bubble.tiff"
  ),
  plot = final_plot,
  width = 18,
  height = 8,
  units = "in",
  dpi = 600,
  compression = "lzw"
)


------------------------------------------------------------------------------
  # ==============================================================================
# (7) UCell pathway scoring and visualization
#
# Description:
# This script calculates a UCell score for the gene set "dendrite development"
# and visualizes the score distribution across the specified experimental groups.
#
# Software used in the original analysis:
#   UCell 2.1.1
#   ggrastr 1.0.1
# ==============================================================================


# 1. Load packages -------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(UCell)
  library(ggrastr)
  library(ggplot2)
  library(ggpubr)
  library(dplyr)
})


# 2. Define paths ---------------------------------------------------------------

input_file <- "DATA/OUTPUT/scRNA_HarmonyIntegrated_annotated.rds"
geneset_file <- "DATA/step0-Microglia_Pathways_geneset0.rds"

output_figure_dir <- "FIGURE"
output_data_dir <- "DATA/OUTPUT"

dir.create(output_figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_data_dir, recursive = TRUE, showWarnings = FALSE)


# 3. Define analysis parameters ------------------------------------------------

pathway_name <- "dendrite development"

group_levels <- c(
  "NS_P3",
  "NS_P7",
  "NS_P12"
)

comparisons <- list(
  c("NS_P3", "NS_P7"),
  c("NS_P7", "NS_P12"),
  c("NS_P3", "NS_P12")
)

group_colors <- c(
  "NS_P3" = "#507BA8",
  "NS_P7" = "#F38D37",
  "NS_P12" = "#F1CE60"
)

random_seed <- 1234
set.seed(random_seed)


# 4. Load Seurat object and gene sets ------------------------------------------

scRNA <- readRDS(input_file)
markers <- readRDS(geneset_file)

if (!inherits(scRNA, "Seurat")) {
  stop("The input file must contain a Seurat object.")
}

if (!"group" %in% colnames(scRNA@meta.data)) {
  stop("The Seurat metadata must contain a column named 'group'.")
}

if (!pathway_name %in% names(markers)) {
  stop(
    paste0(
      "The requested pathway was not found in the gene-set object: ",
      pathway_name
    )
  )
}


# 5. Prepare gene set -----------------------------------------------------------

features <- list()
features[[pathway_name]] <- unique(as.character(markers[[pathway_name]]))
features[[pathway_name]] <- features[[pathway_name]][
  !is.na(features[[pathway_name]]) &
    features[[pathway_name]] != ""
]

if (length(features[[pathway_name]]) == 0) {
  stop("The selected pathway gene set is empty.")
}


# 6. Calculate UCell score ------------------------------------------------------

marker_score <- AddModuleScore_UCell(
  object = scRNA,
  features = features
)

ucell_columns <- grep(
  "UCell$",
  colnames(marker_score@meta.data),
  value = TRUE
)

colnames(marker_score@meta.data)[
  match(ucell_columns, colnames(marker_score@meta.data))
] <- gsub("\\.", " ", ucell_columns)

score_column <- paste0(pathway_name, "_UCell")

if (!score_column %in% colnames(marker_score@meta.data)) {
  stop(
    paste0(
      "The expected UCell score column was not found: ",
      score_column
    )
  )
}


# 7. Extract plotting data ------------------------------------------------------

plot_data <- FetchData(
  object = marker_score,
  vars = c("group", score_column)
)

colnames(plot_data) <- c("group", "score")

plot_data <- plot_data |>
  dplyr::filter(
    group %in% group_levels
  )

if (nrow(plot_data) == 0) {
  stop("No cells remained after filtering to the selected groups.")
}

plot_data$group <- factor(
  plot_data$group,
  levels = group_levels
)

cell_number <- plot_data |>
  dplyr::group_by(group) |>
  dplyr::summarise(
    n_cells = n(),
    y = max(score, na.rm = TRUE) + 0.05,
    .groups = "drop"
  ) |>
  dplyr::mutate(
    label = paste0("n=", n_cells)
  )

write.csv(
  plot_data,
  file = file.path(
    output_data_dir,
    "UCell_dendrite_development_scores_by_cell.csv"
  ),
  row.names = FALSE
)

write.csv(
  cell_number,
  file = file.path(
    output_data_dir,
    "UCell_dendrite_development_group_counts.csv"
  ),
  row.names = FALSE
)


# 8. Plot score distributions ---------------------------------------------------

p <- ggplot(
  plot_data,
  aes(
    x = group,
    y = score,
    fill = group,
    color = group
  )
) +
  theme_minimal() +
  theme(
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.text.x = element_text(
      size = 6,
      face = "plain",
      color = "black"
    ),
    axis.text.y = element_text(
      size = 6,
      face = "plain",
      color = "black"
    ),
    plot.title = element_blank(),
    axis.title.y = element_text(
      color = "black",
      size = 8,
      face = "bold",
      vjust = 0.5
    ),
    legend.position = "none"
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = pathway_name
  ) +
  ggrastr::geom_jitter_rast(
    color = "#00000033",
    pch = 19,
    size = 0.8,
    stroke = 0.01,
    position = position_jitter(0.15),
    alpha = 0.3
  ) +
  scale_fill_manual(
    values = group_colors
  ) +
  scale_color_manual(
    values = group_colors
  ) +
  geom_boxplot(
    color = "black",
    outlier.shape = NA,
    alpha = 0.8,
    linewidth = 0.5,
    width = 0.4
  ) +
  ggpubr::stat_compare_means(
    comparisons = comparisons,
    method = "t.test",
    label = "p.signif",
    size = 4,
    vjust = -0.5
  ) +
  geom_text(
    data = cell_number,
    aes(
      x = group,
      y = y,
      label = label
    ),
    inherit.aes = FALSE,
    size = 2.5
  )

ggsave(
  filename = file.path(
    output_figure_dir,
    "dendrite_development_score.pdf"
  ),
  plot = p,
  width = 4.5,
  height = 3.5,
  units = "cm",
  dpi = 600
)

ggsave(
  filename = file.path(
    output_figure_dir,
    "dendrite_development_score.tiff"
  ),
  plot = p,
  width = 4.5,
  height = 3.5,
  units = "cm",
  dpi = 600,
  compression = "lzw"
)


# ==============================================================================
# (8)  Monocle 3 trajectory analysis
#
# Description:
# This script constructs a Monocle 3 cell_data_set from a Seurat object,
# imports the Seurat UMAP coordinates, learns the principal graph, and orders
# cells in pseudotime using PAM microglia as the trajectory root.
# ==============================================================================


# 1. Load packages -------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(monocle3)
  library(dplyr)
  library(ggplot2)
  library(igraph)
})

set.seed(1234)


# 2. Define input and output paths ---------------------------------------------

input_file <- "DATA/INPUT/scRNA2.rds"

output_data_dir <- "DATA/OUTPUT"
output_figure_dir <- "FIGURE/STEP1.5"

dir.create(
  output_data_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  output_figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# 3. Analysis parameters -------------------------------------------------------

num_dimensions <- 50

root_celltypes <- "PAM"

use_partition <- FALSE
close_loop <- TRUE


# 4. Load Seurat object --------------------------------------------------------

scRNA <- readRDS(input_file)

if (!inherits(scRNA, "Seurat")) {
  stop("The input file must contain a Seurat object.")
}

if (!"celltype" %in% colnames(scRNA@meta.data)) {
  stop("The Seurat metadata must contain a column named 'celltype'.")
}

if (!"SCT" %in% names(scRNA@assays)) {
  stop("The Seurat object does not contain an SCT assay.")
}

if (!"umap" %in% names(scRNA@reductions)) {
  stop("The Seurat object does not contain a UMAP reduction.")
}

print(table(scRNA$celltype))

if ("group" %in% colnames(scRNA@meta.data)) {
  print(table(scRNA$group, scRNA$celltype))
}


# 5. Construct the Monocle 3 cell_data_set -------------------------------------

count_matrix <- GetAssayData(
  object = scRNA,
  assay = "SCT",
  slot = "counts"
)

gene_annotation <- data.frame(
  gene_short_name = rownames(count_matrix),
  row.names = rownames(count_matrix)
)

cell_metadata <- scRNA@meta.data

cds <- new_cell_data_set(
  expression_data = count_matrix,
  cell_metadata = cell_metadata,
  gene_metadata = gene_annotation
)


# 6. Preprocess the data --------------------------------------------------------

cds <- preprocess_cds(
  cds,
  num_dim = num_dimensions
)


# 7. Perform UMAP dimensionality reduction -------------------------------------

cds <- reduce_dimension(
  cds,
  preprocess_method = "PCA",
  reduction_method = "UMAP",
  umap.fast_sgd = FALSE,
  cores = 1
)


# 8. Cluster cells --------------------------------------------------------------

# Clustering is required before learning the principal graph.
cds <- cluster_cells(cds)

print(table(cds@clusters$UMAP$partitions))
print(table(cds@clusters$UMAP$clusters))


# 9. Import the Seurat UMAP coordinates ----------------------------------------

seurat_umap <- Embeddings(
  object = scRNA,
  reduction = "umap"
)

if (!all(colnames(cds) %in% rownames(seurat_umap))) {
  stop("Cell names in the Seurat UMAP do not match those in the Monocle object.")
}

seurat_umap <- seurat_umap[
  colnames(cds),
  ,
  drop = FALSE
]

cds@int_colData$reducedDims$UMAP <- seurat_umap


# 10. Learn the principal graph ------------------------------------------------

cds <- learn_graph(
  cds,
  use_partition = use_partition,
  close_loop = close_loop,
  verbose = TRUE
)


# 11. Define the trajectory root -----------------------------------------------

get_root_pr_nodes <- function(
    cds,
    root_celltypes,
    celltype_column = "celltype"
) {

  cell_metadata <- colData(cds)

  root_cell_indices <- which(
    cell_metadata[[celltype_column]] %in% root_celltypes
  )

  if (length(root_cell_indices) == 0) {
    stop(
      paste0(
        "No cells were found for the selected root cell type(s): ",
        paste(root_celltypes, collapse = ", ")
      )
    )
  }

  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$
    pr_graph_cell_proj_closest_vertex

  closest_vertex <- as.matrix(
    closest_vertex[
      colnames(cds),
      ,
      drop = FALSE
    ]
  )

  most_frequent_vertex <- names(
    which.max(
      table(
        closest_vertex[root_cell_indices, 1]
      )
    )
  )

  root_vertex_index <- as.numeric(most_frequent_vertex)

  root_pr_nodes <-
    igraph::V(
      principal_graph(cds)[["UMAP"]]
    )$name[root_vertex_index]

  return(root_pr_nodes)
}

root_pr_nodes <- get_root_pr_nodes(
  cds = cds,
  root_celltypes = root_celltypes
)

message(
  "Selected root principal node: ",
  paste(root_pr_nodes, collapse = ", ")
)


# 12. Order cells in pseudotime ------------------------------------------------

cds <- order_cells(
  cds,
  root_pr_nodes = root_pr_nodes
)


# 13. Define the cell-type colour palette --------------------------------------

celltype_levels <- levels(
  factor(colData(cds)$celltype)
)

celltype_colors <- c(
  "#66C2A5",
  "#FC8D62",
  "#8DA0CB",
  "#E78AC3",
  "#A6D854",
  "#FFD92F",
  "#E5C494",
  "#B3B3B3",
  "#A6CEE3",
  "#1F78B4",
  "#9A3232",
  "#8F6830",
  "#8F8F30",
  "#4F8F30",
  "#308F68"
)

if (length(celltype_levels) > length(celltype_colors)) {
  stop("The colour palette contains fewer colours than cell types.")
}

celltype_colors <- celltype_colors[
  seq_along(celltype_levels)
]

names(celltype_colors) <- celltype_levels


# 14. Plot the trajectory by cell type -----------------------------------------

p_trajectory <- plot_cells(
  cds,
  color_cells_by = "celltype",
  label_groups_by_cluster = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE
) +
  scale_colour_manual(
    values = celltype_colors
  ) +
  theme(
    panel.border = element_rect(
      fill = NA,
      colour = "black",
      linewidth = 0.5,
      linetype = "solid"
    )
  )

ggsave(
  filename = file.path(
    output_figure_dir,
    "UMAP_based_on_Seurat_trajectory.pdf"
  ),
  plot = p_trajectory,
  width = 4,
  height = 4,
  units = "in",
  device = cairo_pdf
)


# 15. Plot pseudotime -----------------------------------------------------------

p_pseudotime <- plot_cells(
  cds,
  color_cells_by = "pseudotime",
  label_groups_by_cluster = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE
) +
  theme(
    panel.border = element_rect(
      fill = NA,
      colour = "black",
      linewidth = 0.5,
      linetype = "solid"
    )
  )

ggsave(
  filename = file.path(
    output_figure_dir,
    "UMAP_pseudotime.pdf"
  ),
  plot = p_pseudotime,
  width = 4,
  height = 4,
  units = "in",
  device = cairo_pdf
)


# 16. Export pseudotime values -------------------------------------------------

pseudotime_values <- monocle3::pseudotime(cds)

pseudotime_table <- data.frame(
  cell_id = names(pseudotime_values),
  celltype = colData(cds)[
    names(pseudotime_values),
    "celltype"
  ],
  pseudotime = as.numeric(pseudotime_values),
  row.names = NULL
)

if ("group" %in% colnames(colData(cds))) {
  pseudotime_table$group <- colData(cds)[
    pseudotime_table$cell_id,
    "group"
  ]
}

write.csv(
  pseudotime_table,
  file = file.path(
    output_data_dir,
    "Monocle3_cell_pseudotime.csv"
  ),
  row.names = FALSE
)


# 17. Save the ordered Monocle object ------------------------------------------

saveRDS(
  cds,
  file = file.path(
    output_data_dir,
    "Monocle3_ordered_cds.rds"
  )
)


# 18. Save analysis parameters -------------------------------------------------

parameter_information <- c(
  paste0("Input file: ", input_file),
  paste0("Expression assay: SCT"),
  paste0("Expression slot: counts"),
  paste0("Number of PCA dimensions: ", num_dimensions),
  paste0(
    "Root cell type(s): ",
    paste(root_celltypes, collapse = ", ")
  ),
  paste0(
    "Root principal node(s): ",
    paste(root_pr_nodes, collapse = ", ")
  ),
  paste0("use_partition: ", use_partition),
  paste0("close_loop: ", close_loop),
  paste0("Random seed: 1234")
)

writeLines(
  parameter_information,
  con = file.path(
    output_data_dir,
    "Monocle3_analysis_parameters.txt"
  )
)


# 19. Save software and package information ------------------------------------

writeLines(
  capture.output(sessionInfo()),
  con = file.path(
    output_data_dir,
    "Monocle3_sessionInfo.txt"
  )
)

# ==============================================================================
# Reorder microglial states and plot the final module-expression heatmap
# ==============================================================================

custom_order <- c(
  "Mg5",
  "Mg6",
  "Mg8",
  "Mg4",
  "Mg3",
  "Mg7",
  "Mg2",
  "Mg1"
)

missing_columns <- setdiff(
  custom_order,
  colnames(aggregate_module_matrix)
)

if (length(missing_columns) > 0) {
  stop(
    paste0(
      "The following cell types are missing from the aggregated matrix: ",
      paste(missing_columns, collapse = ", ")
    )
  )
}

aggregate_module_matrix <- aggregate_module_matrix[
  ,
  custom_order,
  drop = FALSE
]

dir.create(
  "FIGURE",
  recursive = TRUE,
  showWarnings = FALSE
)

pdf(
  file = "FIGURE/CellType_GeneCluster_CombineModule.pdf",
  width = 4.5,
  height = 1.8
)

pheatmap::pheatmap(
  aggregate_module_matrix,
  scale = "row",
  show_rownames = TRUE,
  show_colnames = TRUE,
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  annotation_col = NULL,
  annotation_row = NULL,
  fontsize_col = 11,
  angle_col = 45,
  border_color = "white",
  color = rev(
    RColorBrewer::brewer.pal(
      n = 10,
      name = "RdBu"
    )
  )
)

dev.off()
```
# ============================================================================
# Monocle 3 pseudotime heatmap
#
# Description:
# This script generates the pseudotime expression heatmap used in the
# manuscript. Cells are ordered by Monocle 3 pseudotime, grouped into 0.3-unit
# pseudotime intervals, and trajectory-associated genes are ordered by the first
# interval in which their scaled mean expression reaches at least 99% of the
# gene-specific maximum. Five temporal gene clusters are defined using the
# prespecified peak-interval boundaries used in the original analysis.
# ============================================================================

# 1. Load packages -----------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(monocle3)
  library(dplyr)
  library(Hmisc)
  library(openxlsx)
  library(ComplexHeatmap)
  library(viridisLite)
  library(grid)
})

set.seed(1234)


# 2. Define input and output paths -------------------------------------------

cds_file <- "DATA/OUTPUT/Step1.5_cds.rds"
seurat_file <- "ALL NS_microglia.rds"
gene_module_file <- "DATA/OUTPUT/gene_module_df.csv"
label_gene_file <- "DATA/OUTPUT/gene list.xlsx"

output_data_dir <- "DATA/OUTPUT"
output_figure_dir <- "FIGURE/STEP3-2"

dir.create(output_data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_figure_dir, recursive = TRUE, showWarnings = FALSE)


# 3. Analysis parameters -----------------------------------------------------

pseudotime_bin_width <- 0.3

# Manual boundaries used in the original analysis. These values refer to the
# index of the pseudotime interval at which each gene first reaches >=99% of its
# gene-specific maximum scaled expression.
temporal_cluster_cutoffs <- c(19, 58, 84, 108)

cluster_levels <- paste0("cluster", 1:5)

cluster_colors <- c(
  "cluster1" = "#EDADC5",
  "cluster2" = "#CEAAD0",
  "cluster3" = "#9584C1",
  "cluster4" = "#6CBEC3",
  "cluster5" = "#AAD7C8"
)


# 4. Load input objects ------------------------------------------------------

cds <- readRDS(cds_file)
scRNA <- readRDS(seurat_file)

gene_module_df <- read.csv(
  gene_module_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (!inherits(cds, "cell_data_set")) {
  stop("The cds input must be a Monocle 3 cell_data_set object.")
}

if (!inherits(scRNA, "Seurat")) {
  stop("The scRNA input must be a Seurat object.")
}

if (!"id" %in% colnames(gene_module_df)) {
  stop("The gene-module table must contain a column named 'id'.")
}

if (!"SCT" %in% names(scRNA@assays)) {
  stop("The Seurat object does not contain an SCT assay.")
}

pseudotime_genes <- unique(as.character(gene_module_df$id))
pseudotime_genes <- pseudotime_genes[
  pseudotime_genes %in% rownames(scRNA)
]

if (length(pseudotime_genes) == 0) {
  stop("None of the genes in gene_module_df$id were found in the Seurat object.")
}


# 5. Extract pseudotime and retain finite values -----------------------------

cell_pseudotime <- monocle3::pseudotime(cds)

finite_cells <- names(cell_pseudotime)[
  is.finite(cell_pseudotime)
]

if (length(finite_cells) == 0) {
  stop("No cells have finite pseudotime values.")
}

cds <- cds[, finite_cells]
cell_pseudotime <- cell_pseudotime[finite_cells]

colData(cds)$pseudotime <- as.numeric(cell_pseudotime)


# 6. Divide pseudotime into 0.3-unit intervals -------------------------------

pseudotime_cuts <- seq(
  min(cell_pseudotime),
  max(cell_pseudotime),
  by = pseudotime_bin_width
)

if (tail(pseudotime_cuts, 1) < max(cell_pseudotime)) {
  pseudotime_cuts <- c(
    pseudotime_cuts,
    max(cell_pseudotime)
  )
}

# Hmisc::cut2 is retained because it was used in the original analysis.
colData(cds)$pseudotime_bin <- Hmisc::cut2(
  as.numeric(cell_pseudotime),
  cuts = pseudotime_cuts
)

write.csv(
  as.data.frame(
    table(
      pseudotime_bin = colData(cds)$pseudotime_bin,
      celltype = colData(cds)$celltype
    )
  ),
  file = file.path(
    output_data_dir,
    "Pseudotime_bin_celltype_counts.csv"
  ),
  row.names = FALSE
)


# 7. Transfer pseudotime information to the Seurat object --------------------

missing_cells <- setdiff(
  colnames(cds),
  colnames(scRNA)
)

if (length(missing_cells) > 0) {
  stop(
    paste0(
      "The following Monocle cells are absent from the Seurat object: ",
      paste(head(missing_cells, 10), collapse = ", ")
    )
  )
}

scRNA <- scRNA[, colnames(cds)]

scRNA$pseudotime <- as.numeric(
  colData(cds)$pseudotime
)

scRNA$pseudotime_bin <- colData(cds)$pseudotime_bin


# 8. Calculate mean expression in each pseudotime interval -------------------

Idents(scRNA) <- "pseudotime_bin"

average_expression_list <- AverageExpression(
  object = scRNA,
  assays = "SCT",
  features = pseudotime_genes,
  slot = "data",
  verbose = FALSE
)

average_expression <- as.matrix(
  average_expression_list[["SCT"]]
)

if (nrow(average_expression) == 0 || ncol(average_expression) == 0) {
  stop("AverageExpression returned an empty matrix.")
}


# 9. Scale each gene to the range 0-1 ----------------------------------------

scale_to_unit_interval <- function(x) {
  x_min <- min(x, na.rm = TRUE)
  x_max <- max(x, na.rm = TRUE)

  if (!is.finite(x_min) || !is.finite(x_max)) {
    return(rep(NA_real_, length(x)))
  }

  if (x_max == x_min) {
    return(rep(0, length(x)))
  }

  (x - x_min) / (x_max - x_min)
}

scaled_expression <- t(
  apply(
    average_expression,
    1,
    scale_to_unit_interval
  )
)

rownames(scaled_expression) <- rownames(average_expression)
colnames(scaled_expression) <- colnames(average_expression)

scaled_expression <- scaled_expression[
  rowSums(is.na(scaled_expression)) == 0,
  ,
  drop = FALSE
]

write.csv(
  scaled_expression,
  file = file.path(
    output_data_dir,
    "Pseudotime_bin_scaled_mean_expression.csv"
  ),
  row.names = TRUE
)


# 10. Order genes by the first near-maximum pseudotime interval ---------------

first_near_maximum_bin <- apply(
  scaled_expression,
  1,
  function(x) {
    candidate_bins <- which(x >= 0.99)

    if (length(candidate_bins) == 0) {
      return(which.max(x))
    }

    candidate_bins[1]
  }
)

gene_order <- order(
  first_near_maximum_bin,
  rownames(scaled_expression)
)

ordered_expression <- scaled_expression[
  gene_order,
  ,
  drop = FALSE
]

ordered_peak_bin <- first_near_maximum_bin[
  gene_order
]


# 11. Assign five temporal gene clusters -------------------------------------

temporal_cluster <- cut(
  ordered_peak_bin,
  breaks = c(
    -Inf,
    temporal_cluster_cutoffs,
    Inf
  ),
  labels = cluster_levels,
  ordered_result = TRUE
)

temporal_cluster <- factor(
  temporal_cluster,
  levels = cluster_levels
)

gene_cluster_table <- data.frame(
  gene = rownames(ordered_expression),
  peak_pseudotime_bin = as.integer(ordered_peak_bin),
  cluster = as.character(temporal_cluster),
  stringsAsFactors = FALSE
)

write.csv(
  gene_cluster_table,
  file = file.path(
    output_data_dir,
    "Pseudotime_gene_temporal_clusters.csv"
  ),
  row.names = FALSE
)

openxlsx::write.xlsx(
  gene_cluster_table,
  file = file.path(
    output_data_dir,
    "Pseudotime_gene_temporal_clusters.xlsx"
  ),
  overwrite = TRUE
)


# 12. Read genes selected for labelling --------------------------------------

label_gene_table <- openxlsx::read.xlsx(
  label_gene_file
)

if (!"gene" %in% colnames(label_gene_table)) {
  stop("The label-gene file must contain a column named 'gene'.")
}

label_genes <- unique(
  as.character(label_gene_table$gene)
)

label_positions <- which(
  rownames(ordered_expression) %in% label_genes
)

label_names <- rownames(ordered_expression)[
  label_positions
]


# 13. Define pathway labels shown beside each temporal cluster ---------------

# These labels are display annotations. The enrichment-analysis code and the
# complete enrichment results used to select these terms should be supplied in
# a separate script and output table.

pathway_labels <- c(
  "cluster1" = paste(
    "chromosome segregation",
    "positive regulation of cell cycle",
    "ATP metabolic process",
    sep = "\n"
  ),
  "cluster2" = paste(
    "tumor necrosis factor production",
    "neuron death",
    "lipid localization",
    "synapse pruning",
    sep = "\n"
  ),
  "cluster3" = paste(
    "myeloid leukocyte migration",
    "positive regulation of cytokine production",
    "positive regulation of immune effector process",
    "regulation of vasculature development",
    "regulation of lymphocyte proliferation",
    "regulation of adaptive immune response",
    sep = "\n"
  ),
  "cluster4" = paste(
    "myeloid cell differentiation",
    "lymphocyte differentiation",
    "BMP signaling pathway",
    "response to transforming growth factor beta",
    "actin filament organization",
    sep = "\n"
  ),
  "cluster5" = paste(
    "glial cell migration",
    "gliogenesis",
    "synapse organization",
    "dendrite development",
    "axonogenesis",
    "axon guidance",
    "axon extension",
    "neuron projection guidance",
    sep = "\n"
  )
)

pathway_text_colors <- c(
  "cluster1" = "#009E73",
  "cluster2" = "#E69F00",
  "cluster3" = "#E69F00",
  "cluster4" = "#E69F00",
  "cluster5" = "#E69F00"
)


# 14. Construct the pseudotime heatmap ---------------------------------------

cluster_matrix <- matrix(
  as.character(temporal_cluster),
  ncol = 1,
  dimnames = list(
    rownames(ordered_expression),
    "Temporal cluster"
  )
)

cluster_heatmap <- Heatmap(
  cluster_matrix,
  name = "Temporal cluster",
  col = cluster_colors,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  width = unit(3, "mm")
)

expression_heatmap <- Heatmap(
  ordered_expression,
  name = "%Max",
  col = viridisLite::viridis(256),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  use_raster = FALSE,
  heatmap_legend_param = list(
    title = "%Max"
  )
)

gene_label_annotation <- rowAnnotation(
  gene = anno_mark(
    at = label_positions,
    labels = label_names,
    labels_gp = gpar(fontsize = 8)
  )
)

heatmap_list <- cluster_heatmap +
  expression_heatmap +
  gene_label_annotation


# 15. Save the final heatmap --------------------------------------------------

output_pdf <- file.path(
  output_figure_dir,
  "Pseudotime_heatmap.pdf"
)

pdf(
  output_pdf,
  width = 4,
  height = 6,
  useDingbats = FALSE
)

draw(
  heatmap_list,
  row_split = temporal_cluster,
  column_title = "Heatmap of pseudotime DEGs",
  column_title_gp = gpar(
    fontsize = 12,
    fontface = "bold"
  ),
  merge_legends = TRUE,
  heatmap_legend_side = "right"
)

# Add the pathway terms used as display annotations in the original figure.
# Their positions are figure-specific and should be checked after rendering.
for (i in seq_along(cluster_levels)) {
  current_cluster <- cluster_levels[i]

  decorate_heatmap_body(
    "Temporal cluster",
    {
      grid.text(
        pathway_labels[[current_cluster]],
        x = unit(-70, "npc"),
        y = unit(0, "npc"),
        just = "centre",
        hjust = 0,
        vjust = 0,
        gp = gpar(
          fontsize = 12,
          col = pathway_text_colors[[current_cluster]],
          fontface = "italic"
        )
      )
    },
    row_slice = i
  )
}

dev.off()


# 16. Save analysis parameters and software information ----------------------

analysis_parameters <- c(
  paste0("Pseudotime bin width: ", pseudotime_bin_width),
  paste0(
    "Temporal cluster cutoffs: ",
    paste(temporal_cluster_cutoffs, collapse = ", ")
  ),
  "Gene ordering rule: first pseudotime interval with scaled expression >= 0.99",
  "Expression assay: SCT",
  "Expression slot: data",
  "Gene scaling: row-wise min-max scaling to 0-1",
  "Random seed: 1234"
)

writeLines(
  analysis_parameters,
  con = file.path(
    output_data_dir,
    "Pseudotime_heatmap_analysis_parameters.txt"
  )
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(
    output_data_dir,
    "Pseudotime_heatmap_sessionInfo.txt"
  )
)


  # ==============================================================================
# (9) scFEA flux post-processing and visualization
#
# Description:
# This script imports the per-cell metabolic flux matrix generated by scFEA,
# matches cells to Seurat-defined microglial/macrophage states, calculates mean
# flux for each state, removes modules with no variation across states, orders
# metabolic modules according to the original analysis, and generates the full
# and selected-module heatmaps used for visualization.
#
# Important:
# This script does not run the scFEA model itself. The exact scFEA command used
# to generate adj_flux.csv should be supplied separately.
# ==============================================================================


# 1. Load packages -------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(openxlsx)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})


# 2. Define input and output paths ---------------------------------------------

flux_file <- "DATA/P3 MACROPHAGY/adj_flux.csv"
module_info_file <- "DATA/P3 MACROPHAGY/scFEA.mouse.moduleinfo.csv"
seurat_file <- "DATA/P3 MACROPHAGY/p3 all macrophagy NEW NAME.rds"
selected_module_file <- "DATA/P3 MACROPHAGY/scFEA_Filter.xlsx"

output_data_dir <- "DATA/OUTPUT"
output_figure_dir <- "FIGURE/STEP2"

dir.create(output_data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_figure_dir, recursive = TRUE, showWarnings = FALSE)


# 3. Analysis parameters -------------------------------------------------------

celltype_order <- c(
  "Mg1",
  "Mg2",
  "Cd11c+Mg",
  "PAM",
  "APIM1",
  "APIM2",
  "BAM1",
  "BAM2",
  "BAM3",
  "Inf_BAM"
)

broad_celltype <- c(
  "Mg1" = "Microglia",
  "Mg2" = "Microglia",
  "Cd11c+Mg" = "Microglia",
  "PAM" = "Microglia",
  "APIM1" = "Microglia",
  "APIM2" = "Microglia",
  "BAM1" = "BAM",
  "BAM2" = "BAM",
  "BAM3" = "BAM",
  "Inf_BAM" = "BAM"
)

# Exact row order used in the original analysis after removal of zero-variance
# modules. This assumes that 164 modules remain.
module_order_index <- c(
  1:142,
  163,
  143:159,
  164,
  160:162
)

heatmap_colors <- circlize::colorRamp2(
  c(-2, -1, 0, 1, 2),
  c(
    "#2166AC",
    "#90C0DC",
    "white",
    "#EF8C65",
    "#B2182B"
  )
)


# 4. Read scFEA flux matrix ----------------------------------------------------

flux_raw <- read.csv(
  flux_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (ncol(flux_raw) < 2) {
  stop("The flux file must contain a barcode column and at least one flux column.")
}

colnames(flux_raw)[1] <- "barcode"

if (anyDuplicated(flux_raw$barcode)) {
  stop("Duplicated cell barcodes were detected in the flux file.")
}

rownames(flux_raw) <- flux_raw$barcode

flux <- flux_raw[
  ,
  setdiff(colnames(flux_raw), "barcode"),
  drop = FALSE
]

flux <- as.data.frame(
  lapply(
    flux,
    function(x) as.numeric(as.character(x))
  ),
  row.names = rownames(flux)
)

if (anyNA(flux)) {
  warning("Missing values were detected after converting the flux matrix to numeric.")
}


# 5. Read module annotations ---------------------------------------------------

module_info <- read.csv(
  module_info_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_module_columns <- c(
  "M_id",
  "M_name",
  "SM_anno"
)

missing_module_columns <- setdiff(
  required_module_columns,
  colnames(module_info)
)

if (length(missing_module_columns) > 0) {
  stop(
    paste0(
      "The module-information file is missing: ",
      paste(missing_module_columns, collapse = ", ")
    )
  )
  )


# 6. Read Seurat object and match cell states ----------------------------------

scRNA <- readRDS(seurat_file)

if (!inherits(scRNA, "Seurat")) {
  stop("The input RDS file must contain a Seurat object.")
}

if (!"celltype" %in% colnames(scRNA@meta.data)) {
  stop("The Seurat metadata must contain a column named 'celltype'.")
}

cell_group <- as.character(
  scRNA@meta.data[
    match(
      rownames(flux),
      rownames(scRNA@meta.data)
    ),
    "celltype"
  ]
)

if (anyNA(cell_group)) {
  unmatched_barcodes <- rownames(flux)[is.na(cell_group)]

  stop(
    paste0(
      "Some scFEA barcodes were not found in the Seurat metadata. Examples: ",
      paste(head(unmatched_barcodes, 10), collapse = ", ")
    )
  )
}

cell_metadata <- data.frame(
  barcode = rownames(flux),
  celltype = cell_group,
  stringsAsFactors = FALSE
)

write.csv(
  cell_metadata,
  file = file.path(
    output_data_dir,
    "scFEA_celltype_metadata.csv"
  ),
  row.names = FALSE
)


# 7. Calculate mean flux by cell state -----------------------------------------

flux_with_group <- flux |>
  tibble::rownames_to_column("barcode") |>
  dplyr::left_join(
    cell_metadata,
    by = "barcode"
  )

mean_flux_by_state <- flux_with_group |>
  dplyr::select(-barcode) |>
  dplyr::group_by(celltype) |>
  dplyr::summarise(
    dplyr::across(
      dplyr::everything(),
      ~ mean(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )

missing_celltypes <- setdiff(
  celltype_order,
  mean_flux_by_state$celltype
)

if (length(missing_celltypes) > 0) {
  stop(
    paste0(
      "The following cell states are missing from the averaged flux table: ",
      paste(missing_celltypes, collapse = ", ")
    )
  )
}

mean_flux_by_state <- mean_flux_by_state[
  match(
    celltype_order,
    mean_flux_by_state$celltype
  ),
  ,
  drop = FALSE
]

df_flux <- t(
  as.matrix(
    mean_flux_by_state[
      ,
      setdiff(
        colnames(mean_flux_by_state),
        "celltype"
      ),
      drop = FALSE
    ]
  )
)

colnames(df_flux) <- mean_flux_by_state$celltype
storage.mode(df_flux) <- "numeric"


# 8. Remove modules with zero variance -----------------------------------------

module_sd <- apply(
  df_flux,
  1,
  stats::sd,
  na.rm = TRUE
)

keep_modules <- is.finite(module_sd) & module_sd != 0

df_flux <- df_flux[
  keep_modules,
  ,
  drop = FALSE
]

message(
  "Number of retained metabolic modules: ",
  nrow(df_flux)
)

if (nrow(df_flux) != length(module_order_index)) {
  stop(
    paste0(
      "The original module-order vector expects ",
      length(module_order_index),
      " retained modules, but ",
      nrow(df_flux),
      " were found. Verify the input file and original ordering."
    )
  )
}


# 9. Apply the original metabolic-module order ---------------------------------

df_flux <- df_flux[
  module_order_index,
  ,
  drop = FALSE
]

module_info <- module_info[
  match(
    rownames(df_flux),
    module_info$M_id
  ),
  ,
  drop = FALSE
]

if (anyNA(module_info$M_id)) {
  stop("Some retained flux modules are absent from the module-information file.")
}


# 10. Save processed state-level flux data -------------------------------------

write.csv(
  df_flux,
  file = file.path(
    output_data_dir,
    "scFEA_mean_flux_by_celltype.csv"
  ),
  row.names = TRUE
)

write.csv(
  module_info,
  file = file.path(
    output_data_dir,
    "scFEA_retained_module_information.csv"
  ),
  row.names = FALSE
)

saveRDS(
  list(
    flux_per_cell = flux,
    cell_metadata = cell_metadata,
    mean_flux_by_celltype = df_flux,
    module_info = module_info
  ),
  file = file.path(
    output_data_dir,
    "scFEA_flux_postprocessing_objects.rds"
  )
)


# 11. Read and order selected modules ------------------------------------------

selected_module_table <- openxlsx::read.xlsx(
  selected_module_file
)

if (!"M_id" %in% colnames(selected_module_table)) {
  stop("The selected-module file must contain a column named 'M_id'.")
}

selected_module_ids <- as.character(
  selected_module_table$M_id
)

missing_selected_modules <- setdiff(
  selected_module_ids,
  rownames(df_flux)
)

if (length(missing_selected_modules) > 0) {
  stop(
    paste0(
      "Selected modules missing from the processed flux matrix: ",
      paste(missing_selected_modules, collapse = ", ")
    )
  )
}

df_flux_filter <- df_flux[
  selected_module_ids,
  ,
  drop = FALSE
]


# 12. Prepare column annotations -----------------------------------------------

annotation_col <- data.frame(
  subcelltype = factor(
    celltype_order,
    levels = celltype_order
  ),
  celltype = factor(
    unname(broad_celltype[celltype_order]),
    levels = c("Microglia", "BAM")
  ),
  row.names = celltype_order,
  check.names = FALSE
)

celltype_colors <- c(
  "Microglia" = "#92C5DE",
  "BAM" = "#F4A582"
)

subcelltype_colors <- c(
  "Mg1" = "#D1E5F0",
  "Mg2" = "#ABD9E9",
  "Cd11c+Mg" = "#92C5DE",
  "PAM" = "#6BAED6",
  "APIM1" = "#3182BD",
  "APIM2" = "#08519C",
  "BAM1" = "#FDDBC7",
  "BAM2" = "#FDD49E",
  "BAM3" = "#FDB863",
  "Inf_BAM" = "#F4A582"
)


# 13. Prepare full-matrix row annotations --------------------------------------

annotation_row_full <- data.frame(
  pathway = factor(
    module_info$SM_anno,
    levels = unique(module_info$SM_anno)
  ),
  row.names = module_info$M_name,
  check.names = FALSE
)

df_flux_full_plot <- df_flux
rownames(df_flux_full_plot) <- module_info$M_name

pathway_palette_base <- c(
  "#8DD3C7",
  "#FFFFB3",
  "#BEBADA",
  "#FB8072",
  "#80B1D3",
  "#FDB462",
  "#B3DE69",
  "#FCCDE5",
  "#D9D9D9",
  "#BC80BD",
  "#CCEBC5",
  "#FFED6F"
)

pathway_levels <- levels(annotation_row_full$pathway)

pathway_colors <- setNames(
  rep(
    pathway_palette_base,
    length.out = length(pathway_levels)
  ),
  pathway_levels
)

annotation_colors_full <- list(
  celltype = celltype_colors,
  subcelltype = subcelltype_colors,
  pathway = pathway_colors
)

pathway_run_lengths <- rle(
  as.character(annotation_row_full$pathway)
)$lengths

full_row_gaps <- cumsum(pathway_run_lengths)
full_row_gaps <- full_row_gaps[
  full_row_gaps < nrow(df_flux_full_plot)
]


# 14. Plot the full flux heatmap ------------------------------------------------

full_heatmap <- ComplexHeatmap::pheatmap(
  df_flux_full_plot,
  scale = "row",
  show_rownames = TRUE,
  show_colnames = FALSE,
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  color = heatmap_colors,
  annotation_col = annotation_col,
  annotation_row = annotation_row_full,
  annotation_names_row = FALSE,
  annotation_names_col = FALSE,
  column_title = NULL,
  row_title = NULL,
  fontsize_col = 8,
  fontsize_row = 4,
  annotation_colors = annotation_colors_full,
  legend = FALSE,
  annotation_legend = TRUE,
  border_color = "white",
  gaps_row = full_row_gaps
)

pdf(
  file = file.path(
    output_figure_dir,
    "scFEA_flux_heatmap_full.pdf"
  ),
  width = 12,
  height = 50,
  useDingbats = FALSE
)

ComplexHeatmap::draw(
  full_heatmap,
  merge_legends = TRUE
)

dev.off()


# 15. Prepare selected-module annotations --------------------------------------

selected_module_info <- module_info[
  match(
    selected_module_ids,
    module_info$M_id
  ),
  ,
  drop = FALSE
]

if (anyNA(selected_module_info$M_id)) {
  stop("Some selected modules lack module annotations.")
}

df_flux_filter_plot <- df_flux_filter
rownames(df_flux_filter_plot) <- selected_module_info$M_name

# In the original analysis, the selected pathway groups were supplied in
# column X8 of scFEA_Filter.xlsx.
if (!"X8" %in% colnames(selected_module_table)) {
  stop(
    paste0(
      "The selected-module file must contain column 'X8', ",
      "which stores the display pathway groups used in the original analysis."
    )
  )
}

annotation_row_filter <- data.frame(
  pathway = factor(
    as.character(selected_module_table$X8),
    levels = unique(as.character(selected_module_table$X8))
  ),
  row.names = rownames(df_flux_filter_plot),
  check.names = FALSE
)

filter_pathway_levels <- levels(
  annotation_row_filter$pathway
)

filter_pathway_palette <- c(
  "#709CCC",
  "#8AD2C6",
  "#FB8072"
)

if (length(filter_pathway_levels) > length(filter_pathway_palette)) {
  stop(
    "More selected pathway groups were found than colors defined in the original palette."
  )
}

filter_pathway_colors <- setNames(
  filter_pathway_palette[
    seq_along(filter_pathway_levels)
  ],
  filter_pathway_levels
)

annotation_colors_filter <- list(
  celltype = celltype_colors,
  subcelltype = subcelltype_colors,
  pathway = filter_pathway_colors
)

filter_run_lengths <- rle(
  as.character(annotation_row_filter$pathway)
)$lengths

filter_row_gaps <- cumsum(filter_run_lengths)
filter_row_gaps <- filter_row_gaps[
  filter_row_gaps < nrow(df_flux_filter_plot)
]


# 16. Plot the selected-module flux heatmap ------------------------------------

selected_heatmap <- ComplexHeatmap::pheatmap(
  df_flux_filter_plot,
  scale = "row",
  show_rownames = TRUE,
  show_colnames = FALSE,
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  color = heatmap_colors,
  annotation_col = annotation_col,
  annotation_row = annotation_row_filter,
  annotation_names_row = FALSE,
  annotation_names_col = FALSE,
  column_title = NULL,
  row_title = NULL,
  fontsize_col = 8,
  annotation_colors = annotation_colors_filter,
  legend = FALSE,
  annotation_legend = TRUE,
  border_color = "white",
  gaps_row = filter_row_gaps
)

pdf(
  file = file.path(
    output_figure_dir,
    "scFEA_flux_heatmap_selected_modules.pdf"
  ),
  width = 7,
  height = 5,
  useDingbats = FALSE
)

ComplexHeatmap::draw(
  selected_heatmap,
  merge_legends = TRUE
)

dev.off()


# 17. Export selected flux data and annotations --------------------------------

write.csv(
  df_flux_filter,
  file = file.path(
    output_data_dir,
    "scFEA_selected_module_mean_flux.csv"
  ),
  row.names = TRUE
)

write.csv(
  data.frame(
    M_id = selected_module_info$M_id,
    M_name = selected_module_info$M_name,
    pathway_group = as.character(annotation_row_filter$pathway),
    stringsAsFactors = FALSE
  ),
  file = file.path(
    output_data_dir,
    "scFEA_selected_module_annotations.csv"
  ),
  row.names = FALSE
)


# 18. Save analysis parameters and software information ------------------------

analysis_parameters <- c(
  paste0(
    "Cell-state order: ",
    paste(celltype_order, collapse = ", ")
  ),
  paste0(
    "Modules retained after zero-variance filtering: ",
    nrow(df_flux)
  ),
  "Flux summary: arithmetic mean across cells within each Seurat-defined state",
  "Heatmap scaling: row-wise z-score",
  "Column clustering: disabled",
  "Row clustering: disabled",
  "Module order: original prespecified order retained"
)

writeLines(
  analysis_parameters,
  con = file.path(
    output_data_dir,
    "scFEA_flux_postprocessing_parameters.txt"
  )
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(
    output_data_dir,
    "scFEA_flux_postprocessing_sessionInfo.txt"
  )
)


# ==============================================================================
# (10) Pseudobulk transcription-factor activity analysis
#
# Description:
# This script performs the complete transcription-factor activity analysis for
# neonatal mouse microglia using two preconstructed Seurat objects.
#
# The workflow:
#   1. loads and validates the LPS and physiological Seurat objects;
#   2. constructs sample-level pseudobulk profiles by cell type and cluster;
#   3. infers TF activity with CollecTRI and decoupleR using ULM and MLM;
#   4. generates summary tables and manuscript-oriented figures;
#   5. saves a compact RDS result bundle and optional augmented Seurat objects.
#
# Important:
# The metadata column specified as sample_col must identify independent
# biological samples. Do not use a treatment, age, or pooled group label as the
# pseudobulk sample identifier.
# ==============================================================================


# 1. Load packages -------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(decoupleR)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
})


# 2. Define project directories ------------------------------------------------

input_dir <- "DATA/INPUT/RDS"
resource_dir <- "DATA/RESOURCE"
output_dir <- "DATA/OUTPUT/TF_ACTIVITY"

qc_dir <- file.path(output_dir, "01_QC")
table_dir <- file.path(output_dir, "02_TABLES")
figure_dir <- file.path(output_dir, "03_FIGURES")
rds_output_dir <- file.path(output_dir, "04_RDS")

dir.create(resource_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_output_dir, recursive = TRUE, showWarnings = FALSE)


# 3. Define analysis parameters ------------------------------------------------

set.seed(11)

min_cells_per_pseudobulk <- 30
min_targets_per_tf <- 5
normalization_target_sum <- 1e6

tf_methods <- c("ulm", "mlm")
top_tfs_per_profile <- 15
top_tfs_for_heatmap <- 18
top_tfs_for_sample_heatmap <- 20

condition_levels <- c("LPS", "NS")

key_tfs <- c(
  "Stat1",
  "Rela",
  "Spi1",
  "Trp53",
  "Hif1a",
  "Nfe2l2"
)

analysis_options <- list(
  run_celltype_analysis = TRUE,
  run_cluster_analysis = TRUE,
  run_sample_level_analysis = TRUE,
  make_activity_expression_heatmap = TRUE,
  make_key_tf_umap = TRUE,
  save_augmented_seurat_objects = TRUE
)


# 4. Define input-object information -------------------------------------------

# IMPORTANT:
# sample_col must identify biological replicates. The current script assumes
# that the "batch" column corresponds to independent samples. Replace it with
# the correct animal/sample identifier when this assumption is not valid.

sample_information <- tribble(
  ~object_name, ~condition, ~file_path, ~assay,
  ~sample_col, ~group_col, ~celltype_col, ~cluster_col,

  "lps_mg", "LPS",
  file.path(input_dir, "All_lps_mg0326.rds"),
  "RNA", "batch", "group", "celltype", "seurat_clusters",

  "ns_mg", "NS",
  file.path(input_dir, "All_ns_Mg_0326.rds"),
  "RNA", "batch", "group", "celltype", "seurat_clusters"
)


# 5. Define helper functions ----------------------------------------------------

get_assay_matrix <- function(
    object,
    assay,
    layer_name = c("counts", "data")
) {
  layer_name <- match.arg(layer_name)

  matrix_out <- tryCatch(
    GetAssayData(
      object = object,
      assay = assay,
      layer = layer_name
    ),
    error = function(e) {
      GetAssayData(
        object = object,
        assay = assay,
        slot = layer_name
      )
    }
  )

  return(matrix_out)
}


validate_sample_information <- function(sample_information) {

  required_columns <- c(
    "object_name",
    "condition",
    "file_path",
    "assay",
    "sample_col",
    "group_col",
    "celltype_col",
    "cluster_col"
  )

  missing_columns <- setdiff(
    required_columns,
    colnames(sample_information)
  )

  if (length(missing_columns) > 0) {
    stop(
      "sample_information is missing columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  if (anyDuplicated(sample_information$object_name)) {
    stop("Duplicated object names were detected in sample_information.")
  }

  missing_files <- sample_information$file_path[
    !file.exists(sample_information$file_path)
  ]

  if (length(missing_files) > 0) {
    stop(
      "The following input files do not exist:\n",
      paste(missing_files, collapse = "\n")
    )
  }
}


validate_seurat_object <- function(object, sample_row) {

  object_name <- sample_row$object_name[[1]]
  assay_name <- sample_row$assay[[1]]

  if (!inherits(object, "Seurat")) {
    stop(
      "The input file for ",
      object_name,
      " does not contain a Seurat object."
    )
  }

  if (!assay_name %in% Assays(object)) {
    stop(
      "Assay '",
      assay_name,
      "' is absent from object '",
      object_name,
      "'."
    )
  }

  required_metadata <- c(
    sample_row$sample_col[[1]],
    sample_row$group_col[[1]],
    sample_row$celltype_col[[1]],
    sample_row$cluster_col[[1]]
  )

  missing_metadata <- setdiff(
    required_metadata,
    colnames(object@meta.data)
  )

  if (length(missing_metadata) > 0) {
    stop(
      "Object '",
      object_name,
      "' is missing metadata columns: ",
      paste(missing_metadata, collapse = ", ")
    )
  }

  sample_col <- sample_row$sample_col[[1]]
  group_col <- sample_row$group_col[[1]]

  if (sample_col == group_col) {
    stop(
      "For object '",
      object_name,
      "', sample_col and group_col are identical. ",
      "sample_col must identify biological replicates."
    )
  }

  sample_group_map <- object@meta.data %>%
    transmute(
      sample_id = as.character(.data[[sample_col]]),
      experimental_group = as.character(.data[[group_col]])
    ) %>%
    distinct()

  duplicated_sample_ids <- sample_group_map %>%
    count(sample_id, name = "n_group_values") %>%
    filter(n_group_values > 1)

  if (nrow(duplicated_sample_ids) > 0) {
    stop(
      "At least one sample identifier in object '",
      object_name,
      "' maps to more than one experimental group. ",
      "Check sample_col and group_col."
    )
  }

  invisible(TRUE)
}


build_object_summary <- function(
    object,
    sample_row
) {
  tibble(
    object_name = sample_row$object_name[[1]],
    condition = sample_row$condition[[1]],
    file_path = sample_row$file_path[[1]],
    n_cells = ncol(object),
    n_features = nrow(object),
    assays = paste(Assays(object), collapse = ", "),
    default_assay = DefaultAssay(object),
    sample_col = sample_row$sample_col[[1]],
    group_col = sample_row$group_col[[1]],
    celltype_col = sample_row$celltype_col[[1]],
    cluster_col = sample_row$cluster_col[[1]]
  )
}


count_metadata_levels <- function(
    object,
    object_name,
    condition,
    metadata_column,
    metadata_role
) {
  object@meta.data %>%
    transmute(
      level = as.character(.data[[metadata_column]])
    ) %>%
    count(level, name = "n_cells") %>%
    mutate(
      object_name = object_name,
      condition = condition,
      metadata_role = metadata_role,
      metadata_column = metadata_column,
      .before = 1
    )
}


build_pseudobulk <- function(
    object,
    object_name,
    condition,
    assay,
    sample_col,
    group_col,
    grouping_col,
    min_cells = 30
) {

  counts <- get_assay_matrix(
    object = object,
    assay = assay,
    layer_name = "counts"
  )

  metadata <- object@meta.data %>%
    rownames_to_column("barcode")

  metadata <- metadata[
    match(colnames(counts), metadata$barcode),
    ,
    drop = FALSE
  ]

  if (!identical(metadata$barcode, colnames(counts))) {
    stop(
      "Cell barcodes in metadata and the count matrix are not aligned for ",
      object_name,
      "."
    )
  }

  sample_values <- as.character(metadata[[sample_col]])
  experimental_group_values <- as.character(metadata[[group_col]])
  grouping_values <- as.character(metadata[[grouping_col]])

  invalid_cells <- (
    is.na(sample_values) |
      sample_values == "" |
      is.na(experimental_group_values) |
      experimental_group_values == "" |
      is.na(grouping_values) |
      grouping_values == ""
  )

  if (any(invalid_cells)) {
    stop(
      "Missing sample, group, or grouping annotations were detected in ",
      object_name,
      "."
    )
  }

  pseudobulk_ids <- paste(
    sample_values,
    grouping_values,
    sep = "__"
  )

  pseudobulk_levels <- unique(pseudobulk_ids)
  pseudobulk_index <- match(
    pseudobulk_ids,
    pseudobulk_levels
  )

  design_matrix <- sparseMatrix(
    i = seq_along(pseudobulk_index),
    j = pseudobulk_index,
    x = 1,
    dims = c(
      length(pseudobulk_index),
      length(pseudobulk_levels)
    ),
    dimnames = list(
      colnames(counts),
      pseudobulk_levels
    )
  )

  pseudobulk_counts <- counts %*% design_matrix

  first_cell_index <- match(
    pseudobulk_levels,
    pseudobulk_ids
  )

  pseudobulk_metadata <- tibble(
    pseudobulk_id = pseudobulk_levels,
    object_name = object_name,
    condition_label = condition,
    sample_id = sample_values[first_cell_index],
    experimental_group = experimental_group_values[first_cell_index],
    grouping_level = grouping_col,
    grouping_value = grouping_values[first_cell_index],
    n_cells = tabulate(
      pseudobulk_index,
      nbins = length(pseudobulk_levels)
    )
  )

  keep_profiles <- (
    pseudobulk_metadata$n_cells >= min_cells
  )

  pseudobulk_counts <- pseudobulk_counts[
    ,
    keep_profiles,
    drop = FALSE
  ]

  pseudobulk_metadata <- pseudobulk_metadata[
    keep_profiles,
    ,
    drop = FALSE
  ]

  if (ncol(pseudobulk_counts) == 0) {
    warning(
      "No pseudobulk profiles passed the minimum-cell threshold for ",
      object_name,
      " grouped by ",
      grouping_col,
      "."
    )
  }

  return(
    list(
      counts = pseudobulk_counts,
      metadata = pseudobulk_metadata
    )
  )
}


normalize_pseudobulk <- function(
    count_matrix,
    target_sum = 1e6
) {

  if (ncol(count_matrix) == 0) {
    return(count_matrix)
  }

  library_size <- Matrix::colSums(count_matrix)

  keep_profiles <- library_size > 0

  count_matrix <- count_matrix[
    ,
    keep_profiles,
    drop = FALSE
  ]

  library_size <- library_size[keep_profiles]

  normalized_matrix <- count_matrix %*%
    Diagonal(
      x = target_sum / library_size
    )

  normalized_matrix <- as(
    normalized_matrix,
    "dgCMatrix"
  )

  normalized_matrix@x <- log1p(
    normalized_matrix@x
  )

  return(normalized_matrix)
}


run_tf_method <- function(
    normalized_matrix,
    network,
    method,
    min_targets = 5
) {

  if (ncol(normalized_matrix) == 0) {
    return(tibble())
  }

  input_matrix <- as.matrix(normalized_matrix)

  if (method == "ulm") {

    result <- decoupleR::run_ulm(
      mat = input_matrix,
      network = network,
      .source = "source",
      .target = "target",
      .mor = "mor",
      minsize = min_targets
    )

  } else if (method == "mlm") {

    result <- decoupleR::run_mlm(
      mat = input_matrix,
      network = network,
      .source = "source",
      .target = "target",
      .mor = "mor",
      minsize = min_targets
    )

  } else {

    stop(
      "Unsupported TF-activity method: ",
      method
    )
  }

  result <- result %>%
    group_by(condition) %>%
    mutate(
      p_adjusted = p.adjust(
        p_value,
        method = "BH"
      )
    ) %>%
    ungroup()

  return(result)
}


activity_to_matrix <- function(activity_result) {

  if (nrow(activity_result) == 0) {
    return(matrix(numeric(0), nrow = 0, ncol = 0))
  }

  activity_result %>%
    select(
      source,
      condition,
      score
    ) %>%
    pivot_wider(
      names_from = condition,
      values_from = score
    ) %>%
    column_to_rownames("source") %>%
    as.matrix()
}


summarize_top_tfs <- function(
    activity_result,
    top_n = 15
) {

  if (nrow(activity_result) == 0) {
    return(tibble())
  }

  activity_result %>%
    mutate(
      absolute_score = abs(score),
      activity_direction = if_else(
        score >= 0,
        "activated",
        "repressed"
      )
    ) %>%
    group_by(condition) %>%
    slice_max(
      order_by = absolute_score,
      n = top_n,
      with_ties = FALSE
    ) %>%
    ungroup() %>%
    arrange(
      condition,
      desc(absolute_score)
    )
}


write_matrix_csv <- function(
    matrix_object,
    file_path,
    row_name = "feature"
) {

  matrix_object %>%
    as.data.frame() %>%
    rownames_to_column(row_name) %>%
    write_csv(file_path)
}


write_tf_results <- function(
    result_list,
    output_subdir
) {

  dir.create(
    output_subdir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  write_csv(
    result_list$metadata,
    file.path(
      output_subdir,
      "pseudobulk_metadata.csv"
    )
  )

  for (method in names(result_list$activity)) {

    activity_result <- result_list$activity[[method]]

    write_csv(
      activity_result,
      file.path(
        output_subdir,
        paste0(method, "_activity_long.csv")
      )
    )

    write_matrix_csv(
      activity_to_matrix(activity_result),
      file.path(
        output_subdir,
        paste0(method, "_activity_matrix.csv")
      ),
      row_name = "tf"
    )

    write_csv(
      summarize_top_tfs(
        activity_result,
        top_n = top_tfs_per_profile
      ),
      file.path(
        output_subdir,
        paste0(method, "_top_tfs.csv")
      )
    )
  }
}


collect_activity_results <- function(
    tf_results,
    grouping_name,
    method_name
) {

  map_dfr(
    names(tf_results),
    function(object_name) {

      object_result <- tf_results[[object_name]][[grouping_name]]

      if (
        is.null(object_result) ||
        is.null(object_result$activity[[method_name]])
      ) {
        return(tibble())
      }

      activity_result <- object_result$activity[[method_name]]
      pseudobulk_metadata <- object_result$metadata

      activity_result %>%
        left_join(
          pseudobulk_metadata,
          by = c(
            "condition" = "pseudobulk_id"
          )
        )
    }
  )
}


z_score_vector <- function(x) {

  x_sd <- sd(
    x,
    na.rm = TRUE
  )

  if (
    length(x) <= 1 ||
    is.na(x_sd) ||
    x_sd == 0
  ) {
    return(
      rep(0, length(x))
    )
  }

  (
    x - mean(x, na.rm = TRUE)
  ) / x_sd
}


save_ggplot <- function(
    plot_object,
    file_stem,
    width,
    height
) {

  ggsave(
    filename = file.path(
      figure_dir,
      paste0(file_stem, ".pdf")
    ),
    plot = plot_object,
    width = width,
    height = height,
    units = "in"
  )

  ggsave(
    filename = file.path(
      figure_dir,
      paste0(file_stem, ".tiff")
    ),
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    dpi = 600,
    compression = "lzw"
  )
}


make_group_activity_heatmap <- function(
    tf_results,
    cell_count_table,
    grouping_name,
    method_name = "ulm",
    top_n = 18,
    figure_stem,
    figure_title
) {

  activity_data <- collect_activity_results(
    tf_results = tf_results,
    grouping_name = grouping_name,
    method_name = method_name
  )

  if (nrow(activity_data) == 0) {
    warning(
      "No activity data were available for ",
      grouping_name,
      "."
    )

    return(NULL)
  }

  mean_activity <- activity_data %>%
    group_by(
      object_name,
      condition_label,
      grouping_value,
      source
    ) %>%
    summarise(
      mean_score = mean(
        score,
        na.rm = TRUE
      ),
      .groups = "drop"
    )

  selected_tfs <- mean_activity %>%
    group_by(source) %>%
    summarise(
      activity_sd = sd(
        mean_score,
        na.rm = TRUE
      ),
      .groups = "drop"
    ) %>%
    mutate(
      activity_sd = replace_na(
        activity_sd,
        0
      )
    ) %>%
    arrange(
      desc(activity_sd),
      source
    ) %>%
    slice_head(n = top_n) %>%
    pull(source)

  plot_data <- mean_activity %>%
    filter(
      source %in% selected_tfs
    ) %>%
    left_join(
      cell_count_table %>%
        filter(
          .data$grouping_name == .env$grouping_name
        ),
      by = c(
        "object_name",
        "condition_label",
        "grouping_value"
      )
    ) %>%
    group_by(
      condition_label,
      source
    ) %>%
    mutate(
      z_score = z_score_vector(
        mean_score
      )
    ) %>%
    ungroup() %>%
    mutate(
      condition_label = factor(
        condition_label,
        levels = condition_levels
      ),
      grouping_label = paste0(
        grouping_value,
        "\n(n=",
        n_cells,
        ")"
      ),
      source = factor(
        source,
        levels = rev(selected_tfs)
      )
    )

  grouping_order <- plot_data %>%
    distinct(
      condition_label,
      grouping_label,
      n_cells
    ) %>%
    arrange(
      condition_label,
      desc(n_cells),
      grouping_label
    ) %>%
    pull(grouping_label) %>%
    unique()

  plot_data$grouping_label <- factor(
    plot_data$grouping_label,
    levels = grouping_order
  )

  heatmap_plot <- ggplot(
    plot_data,
    aes(
      x = grouping_label,
      y = source,
      fill = z_score
    )
  ) +
    geom_tile(
      linewidth = 0.2,
      colour = "grey85"
    ) +
    facet_grid(
      . ~ condition_label,
      scales = "free_x",
      space = "free_x"
    ) +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-2.5, 2.5),
      oob = scales::squish,
      name = "Row z-score"
    ) +
    labs(
      title = figure_title,
      x = NULL,
      y = "Transcription factor"
    ) +
    theme_classic(base_size = 11) +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1
      ),
      axis.ticks = element_blank(),
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        linewidth = 0.5
      ),
      strip.background = element_blank(),
      strip.text = element_text(
        face = "bold",
        size = 11
      ),
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      )
    )

  write_csv(
    plot_data,
    file.path(
      table_dir,
      paste0(
        figure_stem,
        "_plot_data.csv"
      )
    )
  )

  save_ggplot(
    plot_object = heatmap_plot,
    file_stem = figure_stem,
    width = 14,
    height = 8
  )

  return(
    list(
      plot = heatmap_plot,
      plot_data = plot_data,
      selected_tfs = selected_tfs
    )
  )
}


build_sample_level_pseudobulk <- function(
    object,
    sample_row
) {

  object_copy <- object
  object_copy$overall_group <- "All_cells"

  build_pseudobulk(
    object = object_copy,
    object_name = sample_row$object_name[[1]],
    condition = sample_row$condition[[1]],
    assay = sample_row$assay[[1]],
    sample_col = sample_row$sample_col[[1]],
    group_col = sample_row$group_col[[1]],
    grouping_col = "overall_group",
    min_cells = min_cells_per_pseudobulk
  )
}


extract_timepoint <- function(x) {

  numeric_value <- gsub(
    "[^0-9]",
    "",
    as.character(x)
  )

  numeric_value[numeric_value == ""] <- NA_character_

  as.numeric(numeric_value)
}


compute_mean_expression_by_group <- function(
    object,
    assay,
    grouping_col,
    features
) {

  object_copy <- object
  DefaultAssay(object_copy) <- assay

  expression_matrix <- tryCatch(
    get_assay_matrix(
      object = object_copy,
      assay = assay,
      layer_name = "data"
    ),
    error = function(e) {
      matrix(
        numeric(0),
        nrow = 0,
        ncol = 0
      )
    }
  )

  if (
    nrow(expression_matrix) == 0 ||
    ncol(expression_matrix) == 0
  ) {
    object_copy <- NormalizeData(
      object = object_copy,
      assay = assay,
      verbose = FALSE
    )

    expression_matrix <- get_assay_matrix(
      object = object_copy,
      assay = assay,
      layer_name = "data"
    )
  }

  available_features <- intersect(
    features,
    rownames(expression_matrix)
  )

  if (length(available_features) == 0) {
    return(tibble())
  }

  expression_matrix <- expression_matrix[
    available_features,
    ,
    drop = FALSE
  ]

  grouping_values <- as.character(
    object_copy@meta.data[
      colnames(expression_matrix),
      grouping_col
    ]
  )

  grouping_levels <- unique(
    grouping_values
  )

  grouping_index <- match(
    grouping_values,
    grouping_levels
  )

  design_matrix <- sparseMatrix(
    i = seq_along(grouping_index),
    j = grouping_index,
    x = 1,
    dims = c(
      length(grouping_index),
      length(grouping_levels)
    ),
    dimnames = list(
      colnames(expression_matrix),
      grouping_levels
    )
  )

  group_sum <- expression_matrix %*%
    design_matrix

  group_size <- tabulate(
    grouping_index,
    nbins = length(grouping_levels)
  )

  group_mean <- group_sum %*%
    Diagonal(
      x = 1 / group_size
    )

  group_mean %>%
    as.matrix() %>%
    as.data.frame() %>%
    rownames_to_column("tf") %>%
    pivot_longer(
      cols = -tf,
      names_to = "grouping_value",
      values_to = "mean_value"
    )
}


# 6. Load and validate Seurat objects ------------------------------------------

validate_sample_information(
  sample_information
)

seurat_objects <- map(
  sample_information$file_path,
  readRDS
)

names(seurat_objects) <- sample_information$object_name

walk2(
  seurat_objects,
  split(
    sample_information,
    seq_len(nrow(sample_information))
  ),
  validate_seurat_object
)

object_summary <- map_dfr(
  seq_len(nrow(sample_information)),
  function(i) {
    build_object_summary(
      object = seurat_objects[[sample_information$object_name[[i]]]],
      sample_row = sample_information[i, ]
    )
  }
)

metadata_level_counts <- map_dfr(
  seq_len(nrow(sample_information)),
  function(i) {

    sample_row <- sample_information[i, ]
    object_name <- sample_row$object_name[[1]]
    condition <- sample_row$condition[[1]]
    object <- seurat_objects[[object_name]]

    bind_rows(
      count_metadata_levels(
        object = object,
        object_name = object_name,
        condition = condition,
        metadata_column = sample_row$sample_col[[1]],
        metadata_role = "sample"
      ),
      count_metadata_levels(
        object = object,
        object_name = object_name,
        condition = condition,
        metadata_column = sample_row$group_col[[1]],
        metadata_role = "experimental_group"
      ),
      count_metadata_levels(
        object = object,
        object_name = object_name,
        condition = condition,
        metadata_column = sample_row$celltype_col[[1]],
        metadata_role = "celltype"
      ),
      count_metadata_levels(
        object = object,
        object_name = object_name,
        condition = condition,
        metadata_column = sample_row$cluster_col[[1]],
        metadata_role = "seurat_cluster"
      )
    )
  }
)

write_csv(
  object_summary,
  file.path(
    qc_dir,
    "object_summary.csv"
  )
)

write_csv(
  metadata_level_counts,
  file.path(
    qc_dir,
    "metadata_level_counts.csv"
  )
)


# 7. Load the mouse CollecTRI network ------------------------------------------

network_file <- file.path(
  resource_dir,
  "collectri_mouse.csv"
)

if (file.exists(network_file)) {

  collectri_network <- read_csv(
    network_file,
    show_col_types = FALSE
  )

} else {

  collectri_network <- decoupleR::get_collectri(
    organism = "mouse",
    split_complexes = FALSE
  )

  write_csv(
    collectri_network,
    network_file
  )
}

if (
  "weight" %in% colnames(collectri_network) &&
  !"mor" %in% colnames(collectri_network)
) {
  collectri_network <- collectri_network %>%
    rename(
      mor = weight
    )
}

required_network_columns <- c(
  "source",
  "target",
  "mor"
)

missing_network_columns <- setdiff(
  required_network_columns,
  colnames(collectri_network)
)

if (length(missing_network_columns) > 0) {
  stop(
    "The CollecTRI network is missing columns: ",
    paste(
      missing_network_columns,
      collapse = ", "
    )
  )
}

collectri_network <- collectri_network %>%
  select(
    source,
    target,
    mor
  ) %>%
  filter(
    !is.na(source),
    !is.na(target),
    !is.na(mor)
  ) %>%
  distinct()


# 8. Construct pseudobulk profiles and infer TF activity -----------------------

tf_results <- list()
method_errors <- tibble()

grouping_plan <- list()

if (analysis_options$run_celltype_analysis) {
  grouping_plan$celltype <- "celltype_col"
}

if (analysis_options$run_cluster_analysis) {
  grouping_plan$seurat_clusters <- "cluster_col"
}

for (i in seq_len(nrow(sample_information))) {

  sample_row <- sample_information[i, ]
  object_name <- sample_row$object_name[[1]]
  object <- seurat_objects[[object_name]]

  message(
    "Processing object: ",
    object_name
  )

  tf_results[[object_name]] <- list()

  for (grouping_name in names(grouping_plan)) {

    grouping_column <- sample_row[[grouping_plan[[grouping_name]]]][[1]]

    message(
      "  Grouping by: ",
      grouping_column
    )

    pseudobulk <- build_pseudobulk(
      object = object,
      object_name = object_name,
      condition = sample_row$condition[[1]],
      assay = sample_row$assay[[1]],
      sample_col = sample_row$sample_col[[1]],
      group_col = sample_row$group_col[[1]],
      grouping_col = grouping_column,
      min_cells = min_cells_per_pseudobulk
    )

    normalized_matrix <- normalize_pseudobulk(
      count_matrix = pseudobulk$counts,
      target_sum = normalization_target_sum
    )

    activity_results <- map(
      tf_methods,
      function(method_name) {

        message(
          "    Running ",
          toupper(method_name)
        )

        tryCatch(
          run_tf_method(
            normalized_matrix = normalized_matrix,
            network = collectri_network,
            method = method_name,
            min_targets = min_targets_per_tf
          ),
          error = function(e) {

            method_errors <<- bind_rows(
              method_errors,
              tibble(
                object_name = object_name,
                grouping_name = grouping_name,
                method = method_name,
                error_message = conditionMessage(e)
              )
            )

            warning(
              "TF activity inference failed for ",
              object_name,
              " / ",
              grouping_name,
              " / ",
              method_name,
              ": ",
              conditionMessage(e)
            )

            tibble()
          }
        )
      }
    )

    names(activity_results) <- tf_methods

    tf_results[[object_name]][[grouping_name]] <- list(
      metadata = pseudobulk$metadata,
      activity = activity_results
    )

    result_subdir <- file.path(
      table_dir,
      object_name,
      grouping_name
    )

    write_tf_results(
      result_list = tf_results[[object_name]][[grouping_name]],
      output_subdir = result_subdir
    )
  }
}


# 9. Build combined TF summary tables ------------------------------------------

combined_top_tfs <- map_dfr(
  names(tf_results),
  function(object_name) {

    map_dfr(
      names(tf_results[[object_name]]),
      function(grouping_name) {

        map_dfr(
          names(
            tf_results[[object_name]][[grouping_name]]$activity
          ),
          function(method_name) {

            summarize_top_tfs(
              tf_results[[object_name]][[grouping_name]]$activity[[method_name]],
              top_n = top_tfs_per_profile
            ) %>%
              mutate(
                object_name = object_name,
                grouping_name = grouping_name,
                method = method_name,
                .before = 1
              )
          }
        )
      }
    )
  }
)

recurrent_top_tfs <- combined_top_tfs %>%
  group_by(
    grouping_name,
    method,
    source
  ) %>%
  summarise(
    n_top_hits = n(),
    mean_absolute_score = mean(
      absolute_score,
      na.rm = TRUE
    ),
    maximum_absolute_score = max(
      absolute_score,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  arrange(
    grouping_name,
    method,
    desc(n_top_hits),
    desc(mean_absolute_score),
    source
  )

write_csv(
  combined_top_tfs,
  file.path(
    table_dir,
    "combined_top_tfs.csv"
  )
)

write_csv(
  recurrent_top_tfs,
  file.path(
    table_dir,
    "recurrent_top_tfs.csv"
  )
)

if (nrow(method_errors) > 0) {
  write_csv(
    method_errors,
    file.path(
      table_dir,
      "method_errors.csv"
    )
  )
}


# 10. Prepare cell-count tables for heatmap labels -----------------------------

cell_count_table <- map_dfr(
  seq_len(nrow(sample_information)),
  function(i) {

    sample_row <- sample_information[i, ]
    object_name <- sample_row$object_name[[1]]
    condition <- sample_row$condition[[1]]
    object <- seurat_objects[[object_name]]

    celltype_counts <- object@meta.data %>%
      transmute(
        grouping_value = as.character(
          .data[[sample_row$celltype_col[[1]]]]
        )
      ) %>%
      count(
        grouping_value,
        name = "n_cells"
      ) %>%
      mutate(
        object_name = object_name,
        condition_label = condition,
        grouping_name = "celltype",
        .before = 1
      )

    cluster_counts <- object@meta.data %>%
      transmute(
        grouping_value = as.character(
          .data[[sample_row$cluster_col[[1]]]]
        )
      ) %>%
      count(
        grouping_value,
        name = "n_cells"
      ) %>%
      mutate(
        object_name = object_name,
        condition_label = condition,
        grouping_name = "seurat_clusters",
        .before = 1
      )

    bind_rows(
      celltype_counts,
      cluster_counts
    )
  }
)

write_csv(
  cell_count_table,
  file.path(
    qc_dir,
    "cell_counts_for_heatmaps.csv"
  )
)


# 11. Generate cell-type TF activity heatmap -----------------------------------

celltype_heatmap_result <- NULL

if (analysis_options$run_celltype_analysis) {

  celltype_heatmap_result <- make_group_activity_heatmap(
    tf_results = tf_results,
    cell_count_table = cell_count_table,
    grouping_name = "celltype",
    method_name = "ulm",
    top_n = top_tfs_for_heatmap,
    figure_stem = "celltype_tf_activity_heatmap",
    figure_title = paste0(
      "TF activity across major microglial states"
    )
  )
}


# 12. Generate cluster-level TF activity heatmap -------------------------------

cluster_heatmap_result <- NULL

if (analysis_options$run_cluster_analysis) {

  cluster_heatmap_result <- make_group_activity_heatmap(
    tf_results = tf_results,
    cell_count_table = cell_count_table,
    grouping_name = "seurat_clusters",
    method_name = "ulm",
    top_n = top_tfs_for_heatmap,
    figure_stem = "cluster_tf_activity_heatmap",
    figure_title = paste0(
      "TF activity across Seurat clusters"
    )
  )
}


# 13. Run sample-level TF activity analysis ------------------------------------

sample_level_results <- list()
sample_level_activity <- tibble()

if (analysis_options$run_sample_level_analysis) {

  for (i in seq_len(nrow(sample_information))) {

    sample_row <- sample_information[i, ]
    object_name <- sample_row$object_name[[1]]
    object <- seurat_objects[[object_name]]

    sample_pseudobulk <- build_sample_level_pseudobulk(
      object = object,
      sample_row = sample_row
    )

    sample_normalized <- normalize_pseudobulk(
      count_matrix = sample_pseudobulk$counts,
      target_sum = normalization_target_sum
    )

    sample_ulm <- run_tf_method(
      normalized_matrix = sample_normalized,
      network = collectri_network,
      method = "ulm",
      min_targets = min_targets_per_tf
    )

    sample_level_results[[object_name]] <- list(
      metadata = sample_pseudobulk$metadata,
      activity = sample_ulm
    )

    sample_level_activity <- bind_rows(
      sample_level_activity,
      sample_ulm %>%
        left_join(
          sample_pseudobulk$metadata,
          by = c(
            "condition" = "pseudobulk_id"
          )
        )
    )
  }

  write_csv(
    sample_level_activity,
    file.path(
      table_dir,
      "sample_level_ulm_activity.csv"
    )
  )


  # 13.1 Sample-level TF activity heatmap --------------------------------------

  selected_sample_tfs <- sample_level_activity %>%
    group_by(source) %>%
    summarise(
      activity_sd = sd(
        score,
        na.rm = TRUE
      ),
      .groups = "drop"
    ) %>%
    mutate(
      activity_sd = replace_na(
        activity_sd,
        0
      )
    ) %>%
    arrange(
      desc(activity_sd),
      source
    ) %>%
    slice_head(
      n = top_tfs_for_sample_heatmap
    ) %>%
    pull(source)

  sample_heatmap_data <- sample_level_activity %>%
    filter(
      source %in% selected_sample_tfs
    ) %>%
    group_by(source) %>%
    mutate(
      z_score = z_score_vector(score)
    ) %>%
    ungroup() %>%
    mutate(
      condition_label = factor(
        condition_label,
        levels = condition_levels
      ),
      timepoint = extract_timepoint(
        experimental_group
      ),
      sample_label = paste(
        condition_label,
        experimental_group,
        sample_id,
        sep = " | "
      ),
      source = factor(
        source,
        levels = rev(selected_sample_tfs)
      )
    )

  sample_order <- sample_heatmap_data %>%
    distinct(
      condition_label,
      timepoint,
      sample_id,
      sample_label
    ) %>%
    arrange(
      condition_label,
      timepoint,
      sample_id
    ) %>%
    pull(sample_label)

  sample_heatmap_data$sample_label <- factor(
    sample_heatmap_data$sample_label,
    levels = sample_order
  )

  sample_heatmap_plot <- ggplot(
    sample_heatmap_data,
    aes(
      x = sample_label,
      y = source,
      fill = z_score
    )
  ) +
    geom_tile(
      linewidth = 0.2,
      colour = "grey85"
    ) +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-2.5, 2.5),
      oob = scales::squish,
      name = "Row z-score"
    ) +
    labs(
      title = "Sample-level pseudobulk TF activity",
      x = NULL,
      y = "Transcription factor"
    ) +
    theme_classic(base_size = 11) +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1
      ),
      axis.ticks = element_blank(),
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        linewidth = 0.5
      ),
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      )
    )

  write_csv(
    sample_heatmap_data,
    file.path(
      table_dir,
      "sample_level_tf_heatmap_plot_data.csv"
    )
  )

  save_ggplot(
    plot_object = sample_heatmap_plot,
    file_stem = "sample_level_tf_activity_heatmap",
    width = 12,
    height = 7
  )


  # 13.2 Descriptive LPS-versus-NS TF activity contrast ------------------------

  tf_contrast <- sample_level_activity %>%
    group_by(
      source,
      condition_label
    ) %>%
    summarise(
      mean_score = mean(
        score,
        na.rm = TRUE
      ),
      n_samples = n_distinct(
        sample_id
      ),
      .groups = "drop"
    ) %>%
    select(
      source,
      condition_label,
      mean_score,
      n_samples
    ) %>%
    pivot_wider(
      names_from = condition_label,
      values_from = c(
        mean_score,
        n_samples
      ),
      values_fill = 0
    ) %>%
    mutate(
      delta_LPS_minus_NS = (
        mean_score_LPS -
          mean_score_NS
      )
    )

  tf_contrast_statistics <- sample_level_activity %>%
    group_by(source) %>%
    summarise(
      p_value = tryCatch(
        {
          condition_count <- n_distinct(
            condition_label
          )

          if (
            condition_count == 2 &&
            all(
              table(condition_label) >= 2
            )
          ) {
            t.test(
              score ~ condition_label
            )$p.value
          } else {
            NA_real_
          }
        },
        error = function(e) {
          NA_real_
        }
      ),
      .groups = "drop"
    ) %>%
    mutate(
      p_adjusted = p.adjust(
        p_value,
        method = "BH"
      )
    )

  tf_contrast <- tf_contrast %>%
    left_join(
      tf_contrast_statistics,
      by = "source"
    )

  top_positive <- tf_contrast %>%
    arrange(
      desc(delta_LPS_minus_NS)
    ) %>%
    slice_head(n = 10)

  top_negative <- tf_contrast %>%
    arrange(
      delta_LPS_minus_NS
    ) %>%
    slice_head(n = 10)

  contrast_plot_data <- bind_rows(
    top_positive,
    top_negative
  ) %>%
    distinct(source, .keep_all = TRUE) %>%
    arrange(
      delta_LPS_minus_NS
    ) %>%
    mutate(
      source = factor(
        source,
        levels = source
      ),
      direction = if_else(
        delta_LPS_minus_NS >= 0,
        "Higher in LPS",
        "Higher in NS"
      )
    )

  contrast_plot <- ggplot(
    contrast_plot_data,
    aes(
      x = delta_LPS_minus_NS,
      y = source,
      fill = direction
    )
  ) +
    geom_col(
      width = 0.75
    ) +
    geom_vline(
      xintercept = 0,
      linewidth = 0.6
    ) +
    scale_fill_manual(
      values = c(
        "Higher in LPS" = "#B2182B",
        "Higher in NS" = "#2166AC"
      )
    ) +
    labs(
      title = "Descriptive LPS versus NS TF activity contrast",
      x = "Mean ULM activity difference (LPS - NS)",
      y = NULL,
      fill = NULL
    ) +
    theme_classic(base_size = 11) +
    theme(
      legend.position = "top",
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      )
    )

  write_csv(
    tf_contrast,
    file.path(
      table_dir,
      "sample_level_tf_activity_contrast.csv"
    )
  )

  save_ggplot(
    plot_object = contrast_plot,
    file_stem = "sample_level_tf_activity_contrast",
    width = 8,
    height = 6
  )
}


# 14. Compare TF activity with TF expression -----------------------------------

activity_expression_result <- NULL

if (
  analysis_options$make_activity_expression_heatmap &&
  analysis_options$run_celltype_analysis
) {

  mean_activity_key_tfs <- collect_activity_results(
    tf_results = tf_results,
    grouping_name = "celltype",
    method_name = "ulm"
  ) %>%
    filter(
      source %in% key_tfs
    ) %>%
    group_by(
      object_name,
      condition_label,
      grouping_value,
      source
    ) %>%
    summarise(
      mean_value = mean(
        score,
        na.rm = TRUE
      ),
      .groups = "drop"
    ) %>%
    transmute(
      object_name,
      condition_label,
      grouping_value,
      tf = source,
      measure = "TF activity",
      mean_value
    )

  mean_expression_key_tfs <- map_dfr(
    seq_len(nrow(sample_information)),
    function(i) {

      sample_row <- sample_information[i, ]
      object_name <- sample_row$object_name[[1]]

      compute_mean_expression_by_group(
        object = seurat_objects[[object_name]],
        assay = sample_row$assay[[1]],
        grouping_col = sample_row$celltype_col[[1]],
        features = key_tfs
      ) %>%
        mutate(
          object_name = object_name,
          condition_label = sample_row$condition[[1]],
          measure = "TF expression",
          .before = 1
        )
    }
  )

  activity_expression_data <- bind_rows(
    mean_activity_key_tfs,
    mean_expression_key_tfs
  ) %>%
    group_by(
      condition_label,
      measure,
      tf
    ) %>%
    mutate(
      z_score = z_score_vector(
        mean_value
      )
    ) %>%
    ungroup() %>%
    mutate(
      condition_label = factor(
        condition_label,
        levels = condition_levels
      ),
      measure = factor(
        measure,
        levels = c(
          "TF activity",
          "TF expression"
        )
      ),
      tf = factor(
        tf,
        levels = rev(key_tfs)
      )
    )

  celltype_order <- activity_expression_data %>%
    distinct(grouping_value) %>%
    arrange(grouping_value) %>%
    pull(grouping_value)

  activity_expression_data$grouping_value <- factor(
    activity_expression_data$grouping_value,
    levels = celltype_order
  )

  activity_expression_plot <- ggplot(
    activity_expression_data,
    aes(
      x = grouping_value,
      y = tf,
      fill = z_score
    )
  ) +
    geom_tile(
      linewidth = 0.2,
      colour = "grey85"
    ) +
    facet_grid(
      measure ~ condition_label,
      scales = "free_x",
      space = "free_x"
    ) +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-2.5, 2.5),
      oob = scales::squish,
      name = "Row z-score"
    ) +
    labs(
      title = "TF activity versus TF expression across microglial states",
      x = NULL,
      y = NULL
    ) +
    theme_classic(base_size = 11) +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1
      ),
      axis.ticks = element_blank(),
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        linewidth = 0.5
      ),
      strip.background = element_blank(),
      strip.text = element_text(
        face = "bold"
      ),
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      )
    )

  write_csv(
    activity_expression_data,
    file.path(
      table_dir,
      "key_tf_activity_expression_plot_data.csv"
    )
  )

  save_ggplot(
    plot_object = activity_expression_plot,
    file_stem = "key_tf_activity_vs_expression",
    width = 14,
    height = 7
  )

  activity_expression_result <- list(
    plot = activity_expression_plot,
    plot_data = activity_expression_data
  )
}


# 15. Generate key TF expression UMAP panels -----------------------------------

key_tf_umap_panels <- list()

if (analysis_options$make_key_tf_umap) {

  for (i in seq_len(nrow(sample_information))) {

    sample_row <- sample_information[i, ]
    object_name <- sample_row$object_name[[1]]
    condition <- sample_row$condition[[1]]
    assay_name <- sample_row$assay[[1]]

    object <- seurat_objects[[object_name]]
    DefaultAssay(object) <- assay_name

    if (!"umap" %in% Reductions(object)) {
      warning(
        "Object '",
        object_name,
        "' has no UMAP reduction. UMAP panels were skipped."
      )

      next
    }

    expression_matrix <- tryCatch(
      get_assay_matrix(
        object = object,
        assay = assay_name,
        layer_name = "data"
      ),
      error = function(e) {
        matrix(
          numeric(0),
          nrow = 0,
          ncol = 0
        )
      }
    )

    if (
      nrow(expression_matrix) == 0 ||
      ncol(expression_matrix) == 0
    ) {
      object <- NormalizeData(
        object = object,
        assay = assay_name,
        verbose = FALSE
      )
    }

    available_tfs <- intersect(
      key_tfs,
      rownames(object)
    )

    if (length(available_tfs) == 0) {
      warning(
        "None of the requested TFs were found in object '",
        object_name,
        "'."
      )

      next
    }

    feature_plots <- map(
      available_tfs,
      function(tf_name) {

        FeaturePlot(
          object = object,
          features = tf_name,
          reduction = "umap",
          order = TRUE,
          min.cutoff = "q05",
          max.cutoff = "q95"
        ) +
          ggtitle(tf_name) +
          theme(
            plot.title = element_text(
              size = 11,
              face = "bold",
              hjust = 0.5
            ),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            legend.title = element_text(
              size = 8
            ),
            legend.text = element_text(
              size = 7
            )
          )
      }
    )

    panel <- wrap_plots(
      feature_plots,
      ncol = 3
    ) +
      plot_annotation(
        title = paste0(
          condition,
          " microglia: key TF expression"
        ),
        theme = theme(
          plot.title = element_text(
            size = 14,
            face = "bold",
            hjust = 0.5
          )
        )
      )

    key_tf_umap_panels[[object_name]] <- panel

    save_ggplot(
      plot_object = panel,
      file_stem = paste0(
        "key_tf_umap_",
        object_name
      ),
      width = 12,
      height = 8
    )
  }

  if (length(key_tf_umap_panels) > 1) {

    combined_umap_panel <- wrap_plots(
      key_tf_umap_panels,
      ncol = 1
    )

    save_ggplot(
      plot_object = combined_umap_panel,
      file_stem = "key_tf_umap_combined",
      width = 12,
      height = 16
    )
  }
}


# 16. Save compact and augmented RDS outputs -----------------------------------

lightweight_results <- list(
  created_at = as.character(Sys.time()),
  analysis_parameters = list(
    min_cells_per_pseudobulk = min_cells_per_pseudobulk,
    min_targets_per_tf = min_targets_per_tf,
    normalization_target_sum = normalization_target_sum,
    tf_methods = tf_methods,
    top_tfs_per_profile = top_tfs_per_profile,
    top_tfs_for_heatmap = top_tfs_for_heatmap,
    key_tfs = key_tfs
  ),
  sample_information = sample_information,
  object_summary = object_summary,
  metadata_level_counts = metadata_level_counts,
  network_summary = tibble(
    n_interactions = nrow(
      collectri_network
    ),
    n_tfs = n_distinct(
      collectri_network$source
    ),
    n_targets = n_distinct(
      collectri_network$target
    )
  ),
  tf_results = tf_results,
  combined_top_tfs = combined_top_tfs,
  recurrent_top_tfs = recurrent_top_tfs,
  method_errors = method_errors,
  sample_level_results = sample_level_results
)

saveRDS(
  lightweight_results,
  file = file.path(
    rds_output_dir,
    "tf_activity_analysis_results.rds"
  ),
  compress = "gzip"
)

if (analysis_options$save_augmented_seurat_objects) {

  for (i in seq_len(nrow(sample_information))) {

    sample_row <- sample_information[i, ]
    object_name <- sample_row$object_name[[1]]
    object <- seurat_objects[[object_name]]

    object@misc$tf_activity_analysis <- list(
      created_at = as.character(
        Sys.time()
      ),
      condition = sample_row$condition[[1]],
      analysis_parameters = lightweight_results$analysis_parameters,
      celltype_results = tf_results[[object_name]][["celltype"]],
      cluster_results = tf_results[[object_name]][["seurat_clusters"]],
      sample_level_results = sample_level_results[[object_name]],
      recurrent_top_tfs = recurrent_top_tfs
    )

    output_file <- file.path(
      rds_output_dir,
      paste0(
        object_name,
        "_with_tf_activity.rds"
      )
    )

    saveRDS(
      object,
      file = output_file,
      compress = "gzip"
    )

    # Verify that the saved RDS file can be read.
    verification_object <- tryCatch(
      readRDS(output_file),
      error = function(e) {
        NULL
      }
    )

    if (is.null(verification_object)) {
      stop(
        "RDS verification failed for: ",
        output_file
      )
    }

    rm(verification_object)
  }
}


# 17. Save session information -------------------------------------------------

capture.output(
  sessionInfo(),
  file = file.path(
    output_dir,
    "sessionInfo.txt"
  )
)

message(
  "TF activity analysis completed successfully."
)

message(
  "Main output directory: ",
  normalizePath(
    output_dir,
    mustWork = FALSE
  )
)


# ==============================================================================
# (11) Cross-dataset integration of developmental microglial states
#
# Description:
# This script transfers microglial-state annotations from a neonatal mouse
# reference dataset to the Hammond developmental microglia dataset, integrates
# the two datasets using SCTransform and reciprocal PCA, and generates the three
# panels used for cross-dataset visualization:
#
#   1. integrated UMAP colored by developmental stage;
#   2. integrated UMAP colored by microglial state;
#   3. microglial-state composition across developmental stages.
#
# The script intentionally excludes marker-gene analysis, RNA velocity,
# trajectory inference, and additional feature plots because they are not
# required for these three panels.
# ==============================================================================


# 1. Load packages -------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
})


# 2. Define project directories ------------------------------------------------

input_dir <- "DATA/INPUT/RDS"
output_dir <- "DATA/OUTPUT/CROSS_DATASET_INTEGRATION"

figure_dir <- file.path(output_dir, "FIGURES")
table_dir <- file.path(output_dir, "TABLES")
rds_dir <- file.path(output_dir, "RDS")

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  table_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  rds_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# 3. Define input files and analysis parameters --------------------------------

reference_file <- file.path(
  input_dir,
  "All ns Mg 0326.rds"
)

hammond_file <- file.path(
  input_dir,
  "ham_microglia_res0.3_annotated.rds"
)

reference_celltype_col <- "celltype"
reference_group_col <- "group"
hammond_group_col <- "group"

rna_assay <- "RNA"
integration_assay <- "SCTint"

n_integration_features <- 3000
n_pcs <- 30
integration_dims <- 1:30
umap_seed <- 11

# Re-running SCTransform creates one new SCT model per input object and avoids
# problems caused by multiple pre-existing SCT models in merged Seurat objects.
rerun_sctransform <- TRUE

# The transferred Hammond labels are used without a confidence cutoff by
# default. Set a numeric value such as 0.40 to label lower-confidence cells as
# "Low confidence".
minimum_prediction_score <- NULL

set.seed(umap_seed)


# 4. Define developmental-stage and cell-state order ---------------------------

developmental_stage_levels <- c(
  "E14",
  "P3",
  "P4/P5",
  "P7",
  "P12",
  "P30"
)

celltype_levels <- c(
  "NDM",
  "Mg1",
  "Mg2",
  "PEM",
  "PAM",
  "Cd74+Mg",
  "Cd11c+Mg",
  "Pf4+Mg"
)


# 5. Define plotting colors -----------------------------------------------------

developmental_stage_colors <- c(
  "E14"   = "#3C78B4",
  "P3"    = "#F4B183",
  "P4/P5" = "#7FB8B5",
  "P7"    = "#69B34C",
  "P12"   = "#F2C66D",
  "P30"   = "#FF7F00"
)

microglial_state_colors <- c(
  "NDM"      = "#E31A1C",
  "Mg1"      = "#A6CEE3",
  "Mg2"      = "#1F78B4",
  "PEM"      = "#B2DF8A",
  "PAM"      = "#33A02C",
  "Cd74+Mg"  = "#FDBF6F",
  "Cd11c+Mg" = "#FF7F00",
  "Pf4+Mg"   = "#FB9A99"
)


# 6. Define helper functions ----------------------------------------------------

standardize_stage_names <- function(x) {

  recode(
    as.character(x),
    "NS_P3" = "P3",
    "NS_P7" = "P7",
    "NS_P12" = "P12",
    "P4_P5" = "P4/P5",
    .default = as.character(x)
  )
}


validate_input_object <- function(
    object,
    object_name,
    group_col,
    celltype_col = NULL
) {

  if (!inherits(object, "Seurat")) {
    stop(
      "The input file for ",
      object_name,
      " does not contain a Seurat object."
    )
  }

  if (!rna_assay %in% Assays(object)) {
    stop(
      "The RNA assay is absent from object: ",
      object_name
    )
  }

  required_metadata <- group_col

  if (!is.null(celltype_col)) {
    required_metadata <- c(
      required_metadata,
      celltype_col
    )
  }

  missing_metadata <- setdiff(
    required_metadata,
    colnames(object@meta.data)
  )

  if (length(missing_metadata) > 0) {
    stop(
      "Object '",
      object_name,
      "' is missing metadata columns: ",
      paste(
        missing_metadata,
        collapse = ", "
      )
    )
  }

  invisible(TRUE)
}


run_integration_sctransform <- function(
    object,
    new_assay_name
) {

  DefaultAssay(object) <- rna_assay

  variables_to_regress <- if (
    "percent.mt" %in% colnames(object@meta.data)
  ) {
    "percent.mt"
  } else {
    NULL
  }

  object <- SCTransform(
    object = object,
    assay = rna_assay,
    new.assay.name = new_assay_name,
    vars.to.regress = variables_to_regress,
    variable.features.n = n_integration_features,
    verbose = FALSE
  )

  DefaultAssay(object) <- new_assay_name

  return(object)
}


save_plot <- function(
    plot_object,
    file_stem,
    width,
    height
) {

  ggsave(
    filename = file.path(
      figure_dir,
      paste0(
        file_stem,
        ".pdf"
      )
    ),
    plot = plot_object,
    width = width,
    height = height,
    units = "in"
  )

  ggsave(
    filename = file.path(
      figure_dir,
      paste0(
        file_stem,
        ".tiff"
      )
    ),
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    dpi = 600,
    compression = "lzw"
  )
}


# 7. Load and validate the input objects ---------------------------------------

if (!file.exists(reference_file)) {
  stop(
    "Reference file does not exist: ",
    reference_file
  )
}

if (!file.exists(hammond_file)) {
  stop(
    "Hammond file does not exist: ",
    hammond_file
  )
}

reference <- readRDS(reference_file)
hammond <- readRDS(hammond_file)

validate_input_object(
  object = reference,
  object_name = "reference",
  group_col = reference_group_col,
  celltype_col = reference_celltype_col
)

validate_input_object(
  object = hammond,
  object_name = "Hammond",
  group_col = hammond_group_col
)

# Prefix cell names to ensure uniqueness after integration.
reference <- RenameCells(
  object = reference,
  add.cell.id = "REF"
)

hammond <- RenameCells(
  object = hammond,
  add.cell.id = "HAM"
)

reference$dataset <- "Reference"
hammond$dataset <- "Hammond"


# 8. Prepare one SCT model per input object ------------------------------------

if (rerun_sctransform) {

  message(
    "Re-running SCTransform for the reference object."
  )

  reference <- run_integration_sctransform(
    object = reference,
    new_assay_name = integration_assay
  )

  message(
    "Re-running SCTransform for the Hammond object."
  )

  hammond <- run_integration_sctransform(
    object = hammond,
    new_assay_name = integration_assay
  )

} else {

  if (!"SCT" %in% Assays(reference) ||
      !"SCT" %in% Assays(hammond)) {
    stop(
      "Both objects must contain an SCT assay when ",
      "rerun_sctransform is FALSE."
    )
  }

  integration_assay <- "SCT"

  DefaultAssay(reference) <- integration_assay
  DefaultAssay(hammond) <- integration_assay
}


# 9. Define shared genes and transfer features ---------------------------------

shared_rna_genes <- intersect(
  rownames(reference[[rna_assay]]),
  rownames(hammond[[rna_assay]])
)

if (length(shared_rna_genes) == 0) {
  stop(
    "No shared RNA features were detected between the two objects."
  )
}

transfer_features <- intersect(
  VariableFeatures(reference),
  rownames(hammond[[integration_assay]])
)

transfer_features <- intersect(
  transfer_features,
  shared_rna_genes
)

if (length(transfer_features) < 1000) {
  warning(
    "Fewer than 1,000 shared variable features are available ",
    "for label transfer: ",
    length(transfer_features)
  )
}

feature_summary <- tibble(
  statistic = c(
    "Reference RNA genes",
    "Hammond RNA genes",
    "Shared RNA genes",
    "Transfer features"
  ),
  value = c(
    nrow(reference[[rna_assay]]),
    nrow(hammond[[rna_assay]]),
    length(shared_rna_genes),
    length(transfer_features)
  )
)

write_csv(
  feature_summary,
  file.path(
    table_dir,
    "feature_summary.csv"
  )
)


# 10. Transfer reference cell-state labels to Hammond cells --------------------

DefaultAssay(reference) <- integration_assay
DefaultAssay(hammond) <- integration_assay

reference <- RunPCA(
  object = reference,
  assay = integration_assay,
  features = transfer_features,
  npcs = n_pcs,
  verbose = FALSE
)

transfer_anchors <- FindTransferAnchors(
  reference = reference,
  query = hammond,
  normalization.method = "SCT",
  reference.assay = integration_assay,
  query.assay = integration_assay,
  reduction = "pcaproject",
  reference.reduction = "pca",
  features = transfer_features,
  dims = integration_dims
)

celltype_predictions <- TransferData(
  anchorset = transfer_anchors,
  refdata = reference@meta.data[
    [reference_celltype_col]
  ],
  dims = integration_dims
)

hammond <- AddMetaData(
  object = hammond,
  metadata = celltype_predictions
)

hammond$predicted_celltype <- as.character(
  hammond$predicted.id
)

hammond$predicted_celltype_score <- hammond$prediction.score.max

if (!is.null(minimum_prediction_score)) {

  hammond$predicted_celltype <- ifelse(
    hammond$predicted_celltype_score >= minimum_prediction_score,
    hammond$predicted_celltype,
    "Low confidence"
  )
}

label_transfer_summary <- hammond@meta.data %>%
  count(
    .data[[hammond_group_col]],
    predicted_celltype,
    name = "n_cells"
  ) %>%
  rename(
    original_group = 1
  ) %>%
  group_by(original_group) %>%
  mutate(
    percentage = 100 * n_cells / sum(n_cells)
  ) %>%
  ungroup()

write_csv(
  label_transfer_summary,
  file.path(
    table_dir,
    "hammond_label_transfer_summary.csv"
  )
)

prediction_score_summary <- hammond@meta.data %>%
  summarise(
    n_cells = n(),
    minimum_score = min(
      predicted_celltype_score,
      na.rm = TRUE
    ),
    first_quartile = quantile(
      predicted_celltype_score,
      0.25,
      na.rm = TRUE
    ),
    median_score = median(
      predicted_celltype_score,
      na.rm = TRUE
    ),
    mean_score = mean(
      predicted_celltype_score,
      na.rm = TRUE
    ),
    third_quartile = quantile(
      predicted_celltype_score,
      0.75,
      na.rm = TRUE
    ),
    maximum_score = max(
      predicted_celltype_score,
      na.rm = TRUE
    )
  )

write_csv(
  prediction_score_summary,
  file.path(
    table_dir,
    "hammond_prediction_score_summary.csv"
  )
)


# 11. Standardize developmental-stage and cell-state metadata ------------------

reference$developmental_stage <- factor(
  standardize_stage_names(
    reference@meta.data[
      [reference_group_col]
    ]
  ),
  levels = developmental_stage_levels
)

hammond$developmental_stage <- factor(
  standardize_stage_names(
    hammond@meta.data[
      [hammond_group_col]
    ]
  ),
  levels = developmental_stage_levels
)

reference$microglial_state <- factor(
  as.character(
    reference@meta.data[
      [reference_celltype_col]
    ]
  ),
  levels = celltype_levels
)

if (is.null(minimum_prediction_score)) {

  hammond$microglial_state <- factor(
    hammond$predicted_celltype,
    levels = celltype_levels
  )

} else {

  hammond$microglial_state <- factor(
    hammond$predicted_celltype,
    levels = c(
      celltype_levels,
      "Low confidence"
    )
  )
}

unexpected_reference_stages <- setdiff(
  unique(
    as.character(
      reference@meta.data[
        [reference_group_col]
      ]
    )
  ),
  c(
    "NS_P3",
    "NS_P7",
    "NS_P12",
    developmental_stage_levels,
    "P4_P5"
  )
)

unexpected_hammond_stages <- setdiff(
  unique(
    as.character(
      hammond@meta.data[
        [hammond_group_col]
      ]
    )
  ),
  c(
    developmental_stage_levels,
    "P4_P5",
    "NS_P3",
    "NS_P7",
    "NS_P12"
  )
)

if (length(unexpected_reference_stages) > 0) {
  warning(
    "Unexpected developmental-stage labels were found in the ",
    "reference object: ",
    paste(
      unexpected_reference_stages,
      collapse = ", "
    )
  )
}

if (length(unexpected_hammond_stages) > 0) {
  warning(
    "Unexpected developmental-stage labels were found in the ",
    "Hammond object: ",
    paste(
      unexpected_hammond_stages,
      collapse = ", "
    )
  )
}

if (anyNA(reference$developmental_stage) ||
    anyNA(hammond$developmental_stage)) {
  stop(
    "At least one developmental-stage label could not be standardized."
  )
}

if (anyNA(reference$microglial_state)) {
  stop(
    "At least one reference cell-state label is absent from ",
    "celltype_levels."
  )
}

if (is.null(minimum_prediction_score) &&
    anyNA(hammond$microglial_state)) {
  stop(
    "At least one transferred Hammond label is absent from ",
    "celltype_levels."
  )
}


# 12. Select SCT integration features ------------------------------------------

object_list <- list(
  Reference = reference,
  Hammond = hammond
)

integration_features_raw <- SelectIntegrationFeatures(
  object.list = object_list,
  nfeatures = n_integration_features,
  assay = rep(
    integration_assay,
    length(object_list)
  )
)

integration_features <- intersect(
  integration_features_raw,
  shared_rna_genes
)

if (length(integration_features) < 1000) {
  warning(
    "Fewer than 1,000 shared integration features were retained: ",
    length(integration_features)
  )
}

write_csv(
  tibble(
    feature = integration_features
  ),
  file.path(
    table_dir,
    "integration_features.csv"
  )
)


# 13. Prepare SCT objects and identify RPCA anchors ----------------------------

object_list <- PrepSCTIntegration(
  object.list = object_list,
  assay = rep(
    integration_assay,
    length(object_list)
  ),
  anchor.features = integration_features
)

object_list <- lapply(
  object_list,
  function(object) {

    DefaultAssay(object) <- integration_assay

    object <- RunPCA(
      object = object,
      assay = integration_assay,
      features = integration_features,
      npcs = n_pcs,
      verbose = FALSE
    )

    return(object)
  }
)

integration_anchors <- FindIntegrationAnchors(
  object.list = object_list,
  assay = rep(
    integration_assay,
    length(object_list)
  ),
  normalization.method = "SCT",
  anchor.features = integration_features,
  reduction = "rpca",
  reference = 1,
  dims = integration_dims
)


# 14. Integrate the datasets and generate a shared UMAP ------------------------

integrated <- IntegrateData(
  anchorset = integration_anchors,
  normalization.method = "SCT",
  dims = integration_dims
)

DefaultAssay(integrated) <- "integrated"

integrated <- RunPCA(
  object = integrated,
  assay = "integrated",
  npcs = n_pcs,
  verbose = FALSE
)

integrated <- RunUMAP(
  object = integrated,
  reduction = "pca",
  dims = integration_dims,
  seed.use = umap_seed,
  verbose = FALSE
)

# Reapply factor levels after integration.
integrated$developmental_stage <- factor(
  as.character(
    integrated$developmental_stage
  ),
  levels = developmental_stage_levels
)

if (is.null(minimum_prediction_score)) {

  integrated$microglial_state <- factor(
    as.character(
      integrated$microglial_state
    ),
    levels = celltype_levels
  )

} else {

  integrated$microglial_state <- factor(
    as.character(
      integrated$microglial_state
    ),
    levels = c(
      celltype_levels,
      "Low confidence"
    )
  )
}


# 15. Export integrated-object summaries ---------------------------------------

integrated_object_summary <- integrated@meta.data %>%
  count(
    dataset,
    developmental_stage,
    microglial_state,
    name = "n_cells",
    .drop = FALSE
  )

write_csv(
  integrated_object_summary,
  file.path(
    table_dir,
    "integrated_object_summary.csv"
  )
)


# 16. Calculate developmental cell-state composition ---------------------------

composition_table <- integrated@meta.data %>%
  count(
    developmental_stage,
    microglial_state,
    name = "n_cells",
    .drop = FALSE
  ) %>%
  group_by(
    developmental_stage
  ) %>%
  mutate(
    total_cells = sum(
      n_cells
    ),
    percentage = 100 * n_cells / total_cells
  ) %>%
  ungroup()

write_csv(
  composition_table,
  file.path(
    table_dir,
    "microglial_state_composition.csv"
  )
)


# 17. Plot the integrated UMAP by developmental stage --------------------------

p_stage <- DimPlot(
  object = integrated,
  reduction = "umap",
  group.by = "developmental_stage",
  cols = developmental_stage_colors,
  label = TRUE,
  repel = TRUE,
  label.size = 4,
  raster = TRUE
) +
  labs(
    title = NULL,
    colour = "Developmental stage"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_text(
      size = 9
    ),
    legend.text = element_text(
      size = 8
    ),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.5
    )
  )

save_plot(
  plot_object = p_stage,
  file_stem = "Figure_I_developmental_stage_UMAP",
  width = 5,
  height = 4
)


# 18. Plot the integrated UMAP by microglial state -----------------------------

plot_celltype_colors <- microglial_state_colors

if (!is.null(minimum_prediction_score)) {
  plot_celltype_colors <- c(
    plot_celltype_colors,
    "Low confidence" = "grey75"
  )
}

p_celltype <- DimPlot(
  object = integrated,
  reduction = "umap",
  group.by = "microglial_state",
  cols = plot_celltype_colors,
  label = TRUE,
  repel = TRUE,
  label.size = 4,
  raster = TRUE
) +
  labs(
    title = NULL,
    colour = "Microglial state"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_text(
      size = 9
    ),
    legend.text = element_text(
      size = 8
    ),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.5
    )
  )

save_plot(
  plot_object = p_celltype,
  file_stem = "Figure_I_microglial_state_UMAP",
  width = 5,
  height = 4
)


# 19. Plot microglial-state composition across developmental stages ------------

p_composition <- ggplot(
  composition_table,
  aes(
    x = developmental_stage,
    y = percentage,
    fill = microglial_state
  )
) +
  geom_col(
    width = 0.75,
    colour = "black",
    linewidth = 0.2
  ) +
  scale_fill_manual(
    values = plot_celltype_colors,
    limits = names(
      plot_celltype_colors
    ),
    breaks = names(
      plot_celltype_colors
    ),
    drop = FALSE
  ) +
  scale_y_continuous(
    limits = c(
      0,
      100
    ),
    breaks = seq(
      0,
      100,
      by = 20
    ),
    expand = expansion(
      mult = c(
        0,
        0.02
      )
    )
  ) +
  labs(
    x = NULL,
    y = "Percentage of cells",
    fill = "Microglial state"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    legend.title = element_text(
      size = 9
    ),
    legend.text = element_text(
      size = 8
    )
  )

save_plot(
  plot_object = p_composition,
  file_stem = "Figure_J_microglial_state_composition",
  width = 4,
  height = 3.5
)


# 20. Assemble and save the combined figure ------------------------------------

p_umap_pair <- (
  p_stage |
    p_celltype
) +
  plot_annotation(
    title = paste0(
      "Integrated UMAP of physiological microglial states ",
      "across datasets"
    ),
    theme = theme(
      plot.title = element_text(
        size = 12,
        face = "plain",
        hjust = 0.5
      )
    )
  )

p_combined <- (
  p_umap_pair |
    p_composition
) +
  plot_layout(
    widths = c(
      2.2,
      1
    )
  )

save_plot(
  plot_object = p_combined,
  file_stem = "Figure_IJ_combined",
  width = 12,
  height = 4
)


# 21. Save the integrated Seurat object ----------------------------------------

saveRDS(
  integrated,
  file = file.path(
    rds_dir,
    "reference_hammond_integrated.rds"
  ),
  compress = "gzip"
)

saveRDS(
  hammond,
  file = file.path(
    rds_dir,
    "hammond_with_transferred_celltypes.rds"
  ),
  compress = "gzip"
)


# 22. Save session information -------------------------------------------------

capture.output(
  sessionInfo(),
  file = file.path(
    output_dir,
    "sessionInfo.txt"
  )
)

message(
  "Cross-dataset integration completed successfully."
)

message(
  "Main output directory: ",
  normalizePath(
    output_dir,
    mustWork = FALSE
  )
)
