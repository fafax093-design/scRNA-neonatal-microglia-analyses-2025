#!/usr/bin/env Rscript

# ==============================================================================
# 10a. Validate Seurat objects and export count matrices for TF analysis
#
# Description:
# This script validates one or more Seurat objects listed in
# config/sample_manifest.csv and exports the raw count matrix and metadata needed
# by the Python decoupler workflow.
#
# Critical requirement:
# sample_col must identify independent biological samples. It must not be a
# treatment label, age label, pooled group label, or sequencing-batch label.
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(purrr)
  library(tibble)
  library(Matrix)
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

manifest_path <- file.path(
  project_dir,
  "config",
  "sample_manifest.csv"
)

export_root <- file.path(
  project_dir,
  "results",
  "10_tf_activity",
  "01_exports"
)

table_root <- file.path(
  project_dir,
  "TABLE",
  "TF_activity"
)

dir.create(export_root, recursive = TRUE, showWarnings = FALSE)
dir.create(table_root, recursive = TRUE, showWarnings = FALSE)


as_flag <- function(x) {
  tolower(trimws(as.character(x))) %in% c(
    "true",
    "t",
    "1",
    "yes",
    "y"
  )
}


get_count_matrix <- function(object, assay_name) {
  if (!assay_name %in% Assays(object)) {
    stop("Assay was not found: ", assay_name)
  }

  layer_names <- tryCatch(
    Layers(object = object, assay = assay_name),
    error = function(e) character(0)
  )

  count_layers <- grep(
    pattern = "^counts",
    x = layer_names,
    value = TRUE
  )

  if (length(count_layers) > 1) {
    object <- JoinLayers(
      object = object,
      assay = assay_name
    )
  }

  counts <- tryCatch(
    GetAssayData(
      object = object,
      assay = assay_name,
      layer = "counts"
    ),
    error = function(e) {
      GetAssayData(
        object = object,
        assay = assay_name,
        slot = "counts"
      )
    }
  )

  if (nrow(counts) == 0 || ncol(counts) == 0) {
    stop("The selected count matrix is empty.")
  }

  list(
    object = object,
    counts = counts
  )
}


validate_sample_column <- function(
  object,
  sample_col,
  condition_label,
  object_name
) {
  sample_values <- as.character(
    object@meta.data[[sample_col]]
  )

  if (any(is.na(sample_values) | sample_values == "")) {
    stop(
      "Missing values were detected in sample_col for object: ",
      object_name
    )
  }

  n_samples <- dplyr::n_distinct(sample_values)

  if (n_samples < 2) {
    stop(
      "sample_col contains fewer than two independent samples for object: ",
      object_name
    )
  }

  suspicious_names <- c(
    "group",
    "condition",
    "treatment",
    "age",
    "timepoint",
    "batch"
  )

  if (tolower(sample_col) %in% suspicious_names) {
    warning(
      "The configured sample_col is named '",
      sample_col,
      "' for object '",
      object_name,
      "'. Verify that it identifies independent biological samples rather ",
      "than a pooled experimental label."
    )
  }

  if (all(sample_values == condition_label)) {
    stop(
      "sample_col is identical to the condition label for object: ",
      object_name
    )
  }

  invisible(n_samples)
}


if (!file.exists(manifest_path)) {
  stop("Sample manifest does not exist: ", manifest_path)
}

manifest <- suppressMessages(
  readr::read_csv(
    manifest_path,
    show_col_types = FALSE
  )
)

required_manifest_columns <- c(
  "object_name",
  "file_path",
  "condition_label",
  "organism",
  "assay",
  "sample_col",
  "batch_col",
  "celltype_col",
  "cluster_col",
  "reduction",
  "save_augmented_object"
)

missing_manifest_columns <- setdiff(
  required_manifest_columns,
  colnames(manifest)
)

if (length(missing_manifest_columns) > 0) {
  stop(
    "Sample manifest is missing columns: ",
    paste(missing_manifest_columns, collapse = ", ")
  )
}

if (anyDuplicated(manifest$object_name)) {
  stop("Duplicated object_name values were detected in the manifest.")
}

object_summary <- vector("list", nrow(manifest))
sample_summary <- vector("list", nrow(manifest))
grouping_summary <- vector("list", nrow(manifest))

for (i in seq_len(nrow(manifest))) {
  row <- manifest[i, , drop = FALSE]

  object_name <- row$object_name[[1]]
  object_file <- file.path(project_dir, row$file_path[[1]])

  message("Validating and exporting: ", object_name)

  if (!file.exists(object_file)) {
    stop("Seurat object does not exist: ", object_file)
  }

  object <- readRDS(object_file)

  if (!inherits(object, "Seurat")) {
    stop("Input file does not contain a Seurat object: ", object_file)
  }

  required_metadata <- unique(
    c(
      row$sample_col[[1]],
      row$celltype_col[[1]],
      row$cluster_col[[1]]
    )
  )

  batch_col <- trimws(as.character(row$batch_col[[1]]))

  if (!is.na(batch_col) && batch_col != "") {
    required_metadata <- unique(c(required_metadata, batch_col))
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
      paste(missing_metadata, collapse = ", ")
    )
  }

  validate_sample_column(
    object = object,
    sample_col = row$sample_col[[1]],
    condition_label = row$condition_label[[1]],
    object_name = object_name
  )

  if (!row$reduction[[1]] %in% Reductions(object)) {
    warning(
      "Reduction '",
      row$reduction[[1]],
      "' was not found in object '",
      object_name,
      "'. TF inference can continue, but UMAP-expression plotting will fail."
    )
  }

  count_result <- get_count_matrix(
    object = object,
    assay_name = row$assay[[1]]
  )

  object <- count_result$object
  counts <- count_result$counts

  metadata <- object@meta.data[
    colnames(counts),
    ,
    drop = FALSE
  ] %>%
    rownames_to_column("barcode") %>%
    mutate(
      object_name = object_name,
      condition_label = row$condition_label[[1]]
    )

  if (!identical(metadata$barcode, colnames(counts))) {
    stop(
      "Count-matrix columns and metadata barcodes are not aligned for object: ",
      object_name
    )
  }

  output_dir <- file.path(
    export_root,
    object_name
  )

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  Matrix::writeMM(
    counts,
    file.path(output_dir, "counts.mtx")
  )

  write_tsv(
    tibble(gene = rownames(counts)),
    file.path(output_dir, "genes.tsv"),
    col_names = FALSE
  )

  write_tsv(
    tibble(barcode = colnames(counts)),
    file.path(output_dir, "barcodes.tsv"),
    col_names = FALSE
  )

  write_csv(
    metadata,
    file.path(output_dir, "metadata.csv")
  )

  sample_col <- row$sample_col[[1]]

  sample_metadata <- metadata %>%
    mutate(
      sample_id_for_tf = as.character(.data[[sample_col]])
    ) %>%
    group_by(sample_id_for_tf) %>%
    summarise(
      n_cells = n(),
      condition_label = dplyr::first(condition_label),
      object_name = dplyr::first(object_name),
      across(
        everything(),
        ~ {
          values <- unique(as.character(.x))
          values <- values[!is.na(values) & values != ""]
          if (length(values) == 1) values[[1]] else NA_character_
        }
      ),
      .groups = "drop"
    )

  write_csv(
    sample_metadata,
    file.path(output_dir, "sample_metadata.csv")
  )

  object_summary[[i]] <- tibble(
    object_name = object_name,
    file_path = row$file_path[[1]],
    condition_label = row$condition_label[[1]],
    class = paste(class(object), collapse = ", "),
    n_cells = ncol(object),
    n_features = nrow(object),
    assays = paste(Assays(object), collapse = ", "),
    default_assay = DefaultAssay(object),
    assay_exported = row$assay[[1]],
    sample_col = sample_col,
    n_samples = dplyr::n_distinct(metadata[[sample_col]]),
    celltype_col = row$celltype_col[[1]],
    cluster_col = row$cluster_col[[1]],
    reduction = row$reduction[[1]]
  )

  sample_summary[[i]] <- metadata %>%
    count(
      sample_id = .data[[sample_col]],
      name = "n_cells"
    ) %>%
    mutate(
      object_name = object_name,
      condition_label = row$condition_label[[1]],
      .before = 1
    )

  grouping_summary[[i]] <- bind_rows(
    metadata %>%
      count(
        level = .data[[row$celltype_col[[1]]]],
        name = "n_cells"
      ) %>%
      mutate(
        object_name = object_name,
        grouping_level = "celltype",
        metadata_column = row$celltype_col[[1]],
        .before = 1
      ),
    metadata %>%
      count(
        level = .data[[row$cluster_col[[1]]]],
        name = "n_cells"
      ) %>%
      mutate(
        object_name = object_name,
        grouping_level = "cluster",
        metadata_column = row$cluster_col[[1]],
        .before = 1
      )
  )
}

write_csv(
  bind_rows(object_summary),
  file.path(table_root, "TF_input_object_summary.csv")
)

write_csv(
  bind_rows(sample_summary),
  file.path(table_root, "TF_sample_cell_counts.csv")
)

write_csv(
  bind_rows(grouping_summary),
  file.path(table_root, "TF_grouping_cell_counts.csv")
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(table_root, "10a_export_sessionInfo.txt")
)

message("Seurat export completed: ", export_root)
