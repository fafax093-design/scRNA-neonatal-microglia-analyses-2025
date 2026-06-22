#!/usr/bin/env Rscript

# ==============================================================================
# 10d. Selected TF expression on Seurat UMAP embeddings
#
# Description:
# This script visualizes expression of selected TF genes on existing Seurat UMAP
# embeddings. These panels display TF gene expression and must not be described
# as inferred TF-activity maps.
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)

script_path <- if (length(file_arg) > 0) {
  sub("^--file=", "", file_arg[[1]])
} else {
  normalizePath(".", mustWork = TRUE)
}

project_dir <- normalizePath(
  file.path(dirname(script_path), ".."),
  mustWork = FALSE
)

manifest_file <- file.path(
  project_dir,
  "config",
  "sample_manifest.csv"
)

key_tf_file <- file.path(
  project_dir,
  "config",
  "key_tfs.csv"
)

figure_dir <- file.path(
  project_dir,
  "FIGURE",
  "TF_activity",
  "TF_expression_UMAP"
)

table_dir <- file.path(
  project_dir,
  "TABLE",
  "TF_activity"
)

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)


as_flag <- function(x) {
  tolower(trimws(as.character(x))) %in% c(
    "true",
    "t",
    "1",
    "yes",
    "y"
  )
}


manifest <- suppressMessages(
  readr::read_csv(
    manifest_file,
    show_col_types = FALSE
  )
)

key_tf_table <- suppressMessages(
  readr::read_csv(
    key_tf_file,
    show_col_types = FALSE
  )
)

if (!"gene" %in% colnames(key_tf_table)) {
  stop("key_tfs.csv must contain a gene column.")
}

if ("include" %in% colnames(key_tf_table)) {
  key_tf_table <- key_tf_table %>%
    filter(as_flag(include))
}

key_tfs <- unique(
  as.character(
    key_tf_table$gene
  )
)

all_panels <- list()
availability_rows <- list()

for (i in seq_len(nrow(manifest))) {
  row <- manifest[i, , drop = FALSE]

  object_name <- row$object_name[[1]]
  object_file <- file.path(
    project_dir,
    row$file_path[[1]]
  )

  object <- readRDS(object_file)

  if (!inherits(object, "Seurat")) {
    stop(
      "Input file does not contain a Seurat object: ",
      object_file
    )
  }

  assay_name <- row$assay[[1]]
  reduction_name <- row$reduction[[1]]

  if (!assay_name %in% Assays(object)) {
    stop(
      "Assay was not found in object '",
      object_name,
      "': ",
      assay_name
    )
  }

  if (!reduction_name %in% Reductions(object)) {
    warning(
      "Reduction '",
      reduction_name,
      "' was not found for object '",
      object_name,
      "'. UMAP plotting was skipped."
    )
    next
  }

  DefaultAssay(object) <- assay_name

  layer_names <- tryCatch(
    Layers(object = object, assay = assay_name),
    error = function(e) character(0)
  )

  if (length(grep("^data", layer_names)) > 1) {
    object <- JoinLayers(
      object = object,
      assay = assay_name
    )
  }

  normalized_data <- tryCatch(
    GetAssayData(
      object = object,
      assay = assay_name,
      layer = "data"
    ),
    error = function(e) {
      GetAssayData(
        object = object,
        assay = assay_name,
        slot = "data"
      )
    }
  )

  if (ncol(normalized_data) == 0) {
    object <- NormalizeData(
      object = object,
      assay = assay_name,
      verbose = FALSE
    )
  }

  available_tfs <- key_tfs[
    key_tfs %in% rownames(object)
  ]

  missing_tfs <- setdiff(
    key_tfs,
    available_tfs
  )

  availability_rows[[object_name]] <- tibble(
    object_name = object_name,
    gene = key_tfs,
    available = key_tfs %in% available_tfs
  )

  if (length(available_tfs) == 0) {
    warning(
      "No requested TF genes were found for object: ",
      object_name
    )
    next
  }

  plot_list <- lapply(
    available_tfs,
    function(tf_gene) {
      FeaturePlot(
        object = object,
        features = tf_gene,
        reduction = reduction_name,
        order = TRUE,
        min.cutoff = "q05",
        max.cutoff = "q95",
        raster = TRUE
      ) +
        ggtitle(tf_gene) +
        theme(
          plot.title = element_text(
            size = 11,
            face = "bold",
            hjust = 0.5
          ),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 7)
        )
    }
  )

  panel <- wrap_plots(
    plot_list,
    ncol = 3
  ) +
    plot_annotation(
      title = paste0(
        row$condition_label[[1]],
        ": selected TF gene expression"
      ),
      subtitle = "RNA expression on the existing Seurat UMAP"
    )

  output_stem <- file.path(
    figure_dir,
    paste0(
      object_name,
      "_key_TF_expression_UMAP"
    )
  )

  ggsave(
    filename = paste0(output_stem, ".pdf"),
    plot = panel,
    width = 12,
    height = max(
      4,
      ceiling(length(plot_list) / 3) * 3.5
    ),
    units = "in"
  )

  ggsave(
    filename = paste0(output_stem, ".tiff"),
    plot = panel,
    width = 12,
    height = max(
      4,
      ceiling(length(plot_list) / 3) * 3.5
    ),
    units = "in",
    dpi = 600,
    compression = "lzw"
  )

  all_panels[[object_name]] <- panel
}

availability_table <- bind_rows(
  availability_rows
)

write_csv(
  availability_table,
  file.path(
    table_dir,
    "key_TF_expression_availability.csv"
  )
)

if (length(all_panels) > 1) {
  combined_panel <- wrap_plots(
    all_panels,
    ncol = 1
  )

  ggsave(
    filename = file.path(
      figure_dir,
      "combined_key_TF_expression_UMAP.pdf"
    ),
    plot = combined_panel,
    width = 12,
    height = 9 * length(all_panels),
    units = "in"
  )

  ggsave(
    filename = file.path(
      figure_dir,
      "combined_key_TF_expression_UMAP.tiff"
    ),
    plot = combined_panel,
    width = 12,
    height = 9 * length(all_panels),
    units = "in",
    dpi = 600,
    compression = "lzw"
  )
}

writeLines(
  capture.output(sessionInfo()),
  con = file.path(
    table_dir,
    "10d_TF_expression_UMAP_sessionInfo.txt"
  )
)
