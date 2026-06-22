# ==============================================================================
# 11. Cross-dataset label transfer, RPCA integration, and Figure I-J generation
#
# Description:
# This script performs the complete cross-dataset microglia analysis:
#
#   1. loads the physiological reference microglia and Hammond datasets;
#   2. transfers reference microglial-state labels to the Hammond dataset;
#   3. standardizes developmental-stage and microglial-state labels;
#   4. integrates the two datasets using SCT normalization and reciprocal PCA;
#   5. generates a shared PCA and UMAP embedding;
#   6. exports label-transfer and composition tables;
#   7. generates the developmental-stage UMAP, microglial-state UMAP,
#      composition plot, and combined Figure I-J layout;
#   8. saves the integrated Seurat object and software information.
#
# Inputs:
#   data/All ns Mg 0326.rds
#   data/ham_microglia_res0.3_annotated.rds
#
# Main outputs:
#   output/microglia_reference_hammond_integrated.rds
#   output/label_transfer_summary.csv
#   output/celltype_composition.csv
#   output/Figure_I_stage_UMAP.pdf
#   output/Figure_I_celltype_UMAP.pdf
#   output/Figure_J_celltype_composition.pdf
#   output/Figure_IJ_combined.pdf
#
# Software used in the original analysis:
#   Seurat
#   SCP
#   dplyr
#   ggplot2
#   patchwork
# ==============================================================================


# 1. Load packages -------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(SCP)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

random_seed <- 11
set.seed(random_seed)


# 2. Define paths ---------------------------------------------------------------

input_dir <- "data"
output_dir <- "output"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

reference_file <- file.path(
  input_dir,
  "All ns Mg 0326.rds"
)

hammond_file <- file.path(
  input_dir,
  "ham_microglia_res0.3_annotated.rds"
)

integrated_file <- file.path(
  output_dir,
  "microglia_reference_hammond_integrated.rds"
)


# 3. Define analysis parameters ------------------------------------------------

n_transfer_pcs <- 30
n_integration_features <- 3000
integration_dims <- 1:30

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


# 4. Define helper functions ----------------------------------------------------

standardize_group_names <- function(x) {

  dplyr::recode(
    as.character(
      x
    ),
    "NS_P3"  = "P3",
    "NS_P7"  = "P7",
    "NS_P12" = "P12",
    "P4_P5"  = "P4/P5",
    .default = as.character(
      x
    )
  )
}


save_plot <- function(
  plot_object,
  file_stem,
  width,
  height
) {

  ggsave(
    filename = file.path(
      output_dir,
      paste0(
        file_stem,
        ".pdf"
      )
    ),
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    useDingbats = FALSE
  )

  ggsave(
    filename = file.path(
      output_dir,
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


validate_seurat_object <- function(
  object,
  object_name,
  required_metadata
) {

  if (!inherits(
    object,
    "Seurat"
  )) {

    stop(
      object_name,
      " must be a Seurat object."
    )
  }

  missing_metadata <- setdiff(
    required_metadata,
    colnames(
      object@meta.data
    )
  )

  if (length(
    missing_metadata
  ) > 0) {

    stop(
      object_name,
      " is missing metadata columns: ",
      paste(
        missing_metadata,
        collapse = ", "
      )
    )
  }

  if (!"SCT" %in% Assays(
    object
  )) {

    stop(
      object_name,
      " does not contain an SCT assay."
    )
  }

  invisible(
    TRUE
  )
}


# 5. Load and validate the input objects ---------------------------------------

if (!file.exists(
  reference_file
)) {

  stop(
    "Reference file does not exist: ",
    reference_file
  )
}

if (!file.exists(
  hammond_file
)) {

  stop(
    "Hammond file does not exist: ",
    hammond_file
  )
}

reference <- readRDS(
  reference_file
)

hammond <- readRDS(
  hammond_file
)

validate_seurat_object(
  object = reference,
  object_name = "Reference object",
  required_metadata = c(
    "group",
    "celltype"
  )
)

validate_seurat_object(
  object = hammond,
  object_name = "Hammond object",
  required_metadata = "group"
)


# 6. Prepare objects for label transfer and integration ------------------------

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

DefaultAssay(
  reference
) <- "SCT"

DefaultAssay(
  hammond
) <- "SCT"


# 7. Transfer reference microglial-state labels to Hammond cells ---------------

transfer_features <- intersect(
  VariableFeatures(
    reference
  ),
  rownames(
    hammond[
      ["SCT"]
    ]
  )
)

if (length(
  transfer_features
) == 0) {

  stop(
    "No shared variable features were available for label transfer."
  )
}

if (length(
  transfer_features
) < 1000) {

  warning(
    "Fewer than 1,000 shared variable features were available for label ",
    "transfer: ",
    length(
      transfer_features
    )
  )
}

reference <- RunPCA(
  object = reference,
  assay = "SCT",
  features = transfer_features,
  npcs = n_transfer_pcs,
  seed.use = random_seed,
  verbose = FALSE
)

transfer_anchors <- FindTransferAnchors(
  reference = reference,
  query = hammond,
  normalization.method = "SCT",
  reference.assay = "SCT",
  query.assay = "SCT",
  reduction = "pcaproject",
  reference.reduction = "pca",
  features = transfer_features,
  dims = integration_dims
)

predictions <- TransferData(
  anchorset = transfer_anchors,
  refdata = reference$celltype,
  dims = integration_dims
)

hammond <- AddMetaData(
  object = hammond,
  metadata = predictions
)

if (!"predicted.id" %in% colnames(
  hammond@meta.data
)) {

  stop(
    "TransferData did not generate a predicted.id metadata column."
  )
}

hammond$predicted.celltype <- as.character(
  hammond$predicted.id
)

label_transfer_summary <- hammond@meta.data %>%
  count(
    group,
    predicted.celltype,
    name = "cell_number"
  ) %>%
  group_by(
    group
  ) %>%
  mutate(
    percentage = 100 *
      cell_number /
      sum(
        cell_number
      )
  ) %>%
  ungroup() %>%
  arrange(
    group,
    predicted.celltype
  )

write.csv(
  label_transfer_summary,
  file = file.path(
    output_dir,
    "label_transfer_summary.csv"
  ),
  row.names = FALSE
)

if ("prediction.score.max" %in% colnames(
  hammond@meta.data
)) {

  prediction_score_summary <- hammond@meta.data %>%
    transmute(
      group = as.character(
        group
      ),
      predicted.celltype = as.character(
        predicted.celltype
      ),
      prediction.score.max = as.numeric(
        prediction.score.max
      )
    ) %>%
    group_by(
      group,
      predicted.celltype
    ) %>%
    summarise(
      n_cells = n(),
      mean_prediction_score = mean(
        prediction.score.max,
        na.rm = TRUE
      ),
      median_prediction_score = median(
        prediction.score.max,
        na.rm = TRUE
      ),
      .groups = "drop"
    )

  write.csv(
    prediction_score_summary,
    file = file.path(
      output_dir,
      "label_transfer_prediction_score_summary.csv"
    ),
    row.names = FALSE
  )
}


# 8. Standardize developmental-stage and microglial-state labels ---------------

reference$group_standardized <- factor(
  standardize_group_names(
    reference$group
  ),
  levels = developmental_stage_levels
)

hammond$group_standardized <- factor(
  standardize_group_names(
    hammond$group
  ),
  levels = developmental_stage_levels
)

reference$celltype_standardized <- factor(
  as.character(
    reference$celltype
  ),
  levels = celltype_levels
)

hammond$celltype_standardized <- factor(
  as.character(
    hammond$predicted.celltype
  ),
  levels = celltype_levels
)

if (anyNA(
  reference$group_standardized
) ||
    anyNA(
      hammond$group_standardized
    )) {

  unknown_reference_groups <- setdiff(
    unique(
      as.character(
        reference$group
      )
    ),
    names(
      developmental_stage_colors
    )
  )

  unknown_hammond_groups <- setdiff(
    unique(
      standardize_group_names(
        hammond$group
      )
    ),
    developmental_stage_levels
  )

  stop(
    "At least one developmental-stage label was not included in ",
    "developmental_stage_levels. Reference: ",
    paste(
      unknown_reference_groups,
      collapse = ", "
    ),
    "; Hammond: ",
    paste(
      unknown_hammond_groups,
      collapse = ", "
    )
  )
}

if (anyNA(
  reference$celltype_standardized
) ||
    anyNA(
      hammond$celltype_standardized
    )) {

  unknown_reference_celltypes <- setdiff(
    unique(
      as.character(
        reference$celltype
      )
    ),
    celltype_levels
  )

  unknown_hammond_celltypes <- setdiff(
    unique(
      as.character(
        hammond$predicted.celltype
      )
    ),
    celltype_levels
  )

  stop(
    "At least one microglial-state label was not included in celltype_levels. ",
    "Reference: ",
    paste(
      unknown_reference_celltypes,
      collapse = ", "
    ),
    "; Hammond: ",
    paste(
      unknown_hammond_celltypes,
      collapse = ", "
    )
  )
}


# 9. Integrate reference and Hammond datasets using SCT and RPCA ----------------

object_list <- list(
  Reference = reference,
  Hammond = hammond
)

integration_features <- SelectIntegrationFeatures(
  object.list = object_list,
  nfeatures = n_integration_features,
  assay = rep(
    "SCT",
    length(
      object_list
    )
  )
)

object_list <- PrepSCTIntegration(
  object.list = object_list,
  assay = rep(
    "SCT",
    length(
      object_list
    )
  ),
  anchor.features = integration_features
)

object_list <- lapply(
  object_list,
  function(object) {

    DefaultAssay(
      object
    ) <- "SCT"

    RunPCA(
      object = object,
      assay = "SCT",
      features = integration_features,
      npcs = max(
        integration_dims
      ),
      seed.use = random_seed,
      verbose = FALSE
    )
  }
)

integration_anchors <- FindIntegrationAnchors(
  object.list = object_list,
  assay = rep(
    "SCT",
    length(
      object_list
    )
  ),
  normalization.method = "SCT",
  anchor.features = integration_features,
  reduction = "rpca",
  reference = 1,
  dims = integration_dims
)

integrated <- IntegrateData(
  anchorset = integration_anchors,
  normalization.method = "SCT",
  dims = integration_dims
)


# 10. Generate shared PCA and UMAP embeddings ----------------------------------

DefaultAssay(
  integrated
) <- "integrated"

integrated <- RunPCA(
  object = integrated,
  assay = "integrated",
  npcs = max(
    integration_dims
  ),
  seed.use = random_seed,
  verbose = FALSE
)

integrated <- RunUMAP(
  object = integrated,
  reduction = "pca",
  dims = integration_dims,
  seed.use = random_seed,
  reduction.name = "umap",
  reduction.key = "UMAP_",
  verbose = FALSE
)

if (!"umap" %in% Reductions(
  integrated
)) {

  stop(
    "The integrated object does not contain a UMAP reduction."
  )
}

integrated$group_standardized <- factor(
  as.character(
    integrated$group_standardized
  ),
  levels = developmental_stage_levels
)

integrated$celltype_standardized <- factor(
  as.character(
    integrated$celltype_standardized
  ),
  levels = celltype_levels
)


# 11. Export integrated-object summaries ---------------------------------------

integrated_cell_counts <- integrated@meta.data %>%
  transmute(
    dataset = as.character(
      dataset
    ),
    group_standardized = as.character(
      group_standardized
    ),
    celltype_standardized = as.character(
      celltype_standardized
    )
  ) %>%
  count(
    dataset,
    group_standardized,
    celltype_standardized,
    name = "cell_number"
  )

write.csv(
  integrated_cell_counts,
  file = file.path(
    output_dir,
    "integrated_cell_counts_by_dataset_stage_and_state.csv"
  ),
  row.names = FALSE
)

composition_table <- integrated@meta.data %>%
  count(
    group_standardized,
    celltype_standardized,
    name = "cell_number",
    .drop = FALSE
  ) %>%
  group_by(
    group_standardized
  ) %>%
  mutate(
    percentage = 100 *
      cell_number /
      sum(
        cell_number
      )
  ) %>%
  ungroup()

write.csv(
  composition_table,
  file = file.path(
    output_dir,
    "celltype_composition.csv"
  ),
  row.names = FALSE
)

saveRDS(
  integrated,
  file = integrated_file,
  compress = "gzip"
)


# 12. Generate Figure I: UMAP by developmental stage ---------------------------

p_stage <- CellDimPlot(
  srt = integrated,
  group.by = "group_standardized",
  reduction = "umap",
  theme_use = "theme_blank",
  label = TRUE,
  label_insitu = TRUE,
  label.size = 4,
  show_stat = TRUE
) +
  scale_colour_manual(
    values = developmental_stage_colors,
    limits = developmental_stage_levels,
    breaks = developmental_stage_levels,
    drop = FALSE
  ) +
  labs(
    colour = "Developmental stage"
  ) +
  theme(
    legend.title = element_text(
      size = 8
    ),
    legend.text = element_text(
      size = 7
    )
  )


# 13. Generate Figure I: UMAP by microglial state ------------------------------

p_celltype <- CellDimPlot(
  srt = integrated,
  group.by = "celltype_standardized",
  reduction = "umap",
  theme_use = "theme_blank",
  label = TRUE,
  label_insitu = TRUE,
  label.size = 4,
  show_stat = TRUE
) +
  scale_colour_manual(
    values = microglial_state_colors,
    limits = celltype_levels,
    breaks = celltype_levels,
    drop = FALSE
  ) +
  labs(
    colour = "Microglial state"
  ) +
  theme(
    legend.title = element_text(
      size = 8
    ),
    legend.text = element_text(
      size = 7
    )
  )


# 14. Generate Figure J: cell-state composition --------------------------------

p_composition <- CellStatPlot(
  srt = integrated,
  stat.by = "celltype_standardized",
  group.by = "group_standardized",
  plot_type = "trend"
) +
  scale_fill_manual(
    values = microglial_state_colors,
    limits = celltype_levels,
    breaks = celltype_levels,
    drop = FALSE
  ) +
  scale_colour_manual(
    values = microglial_state_colors,
    limits = celltype_levels,
    breaks = celltype_levels,
    drop = FALSE
  ) +
  labs(
    x = NULL,
    y = "Percentage",
    fill = "Microglial state",
    colour = "Microglial state"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    legend.title = element_text(
      size = 8
    ),
    legend.text = element_text(
      size = 7
    )
  )


# 15. Save individual panels ---------------------------------------------------

save_plot(
  plot_object = p_stage,
  file_stem = "Figure_I_stage_UMAP",
  width = 5,
  height = 4
)

save_plot(
  plot_object = p_celltype,
  file_stem = "Figure_I_celltype_UMAP",
  width = 5,
  height = 4
)

save_plot(
  plot_object = p_composition,
  file_stem = "Figure_J_celltype_composition",
  width = 4,
  height = 3.5
)


# 16. Assemble and save the combined Figure I-J --------------------------------

p_umap_pair <- (
  p_stage |
    p_celltype
) +
  plot_annotation(
    title = paste(
      "Integrated UMAP of physiological microglial states",
      "across datasets"
    )
  ) &
  theme(
    plot.title = element_text(
      size = 11,
      face = "plain",
      hjust = 0.5
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


# 17. Save analysis parameters and software information ------------------------

analysis_parameters <- data.frame(
  parameter = c(
    "Reference file",
    "Hammond file",
    "Random seed",
    "Transfer PCA dimensions",
    "Integration features",
    "Integration dimensions",
    "Integration method",
    "Normalization method",
    "Reference dataset index",
    "UMAP reduction"
  ),
  value = c(
    reference_file,
    hammond_file,
    as.character(
      random_seed
    ),
    paste0(
      "1-",
      n_transfer_pcs
    ),
    as.character(
      n_integration_features
    ),
    paste0(
      min(
        integration_dims
      ),
      "-",
      max(
        integration_dims
      )
    ),
    "Reciprocal PCA",
    "SCT",
    "1",
    "PCA-derived UMAP"
  ),
  stringsAsFactors = FALSE
)

write.csv(
  analysis_parameters,
  file = file.path(
    output_dir,
    "cross_dataset_integration_parameters.csv"
  ),
  row.names = FALSE
)

capture.output(
  sessionInfo(),
  file = file.path(
    output_dir,
    "sessionInfo_cross_dataset_integration.txt"
  )
)

message(
  "Cross-dataset label transfer, integration, and Figure I-J generation ",
  "completed successfully."
)

message(
  "Integrated object saved to: ",
  integrated_file
)
