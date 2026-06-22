#!/usr/bin/env Rscript

# ==============================================================================
# 10e. Package TF-analysis outputs into verified RDS files
#
# Description:
# This script creates a compact TF-analysis result bundle and optionally stores
# TF-analysis references in the misc slot of each input Seurat object.
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(readr)
  library(dplyr)
  library(tibble)
  library(jsonlite)
  library(purrr)
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

config_file <- file.path(
  project_dir,
  "config",
  "analysis_config.json"
)

inference_root <- file.path(
  project_dir,
  "results",
  "10_tf_activity",
  "02_inference"
)

table_root <- file.path(
  project_dir,
  "TABLE",
  "TF_activity"
)

output_root <- file.path(
  project_dir,
  "DATA",
  "OUTPUT",
  "TF_activity"
)

dir.create(output_root, recursive = TRUE, showWarnings = FALSE)


as_flag <- function(x) {
  tolower(trimws(as.character(x))) %in% c(
    "true",
    "t",
    "1",
    "yes",
    "y"
  )
}


read_optional_csv <- function(path, row_names = FALSE) {
  if (!file.exists(path)) {
    return(NULL)
  }

  if (row_names) {
    table <- read.csv(
      path,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    rownames(table) <- table[[1]]
    table[[1]] <- NULL
    return(table)
  }

  suppressMessages(
    readr::read_csv(
      path,
      show_col_types = FALSE
    )
  )
}


safe_save_rds <- function(object, path, description) {
  temporary_file <- tempfile(
    tmpdir = dirname(path),
    fileext = ".rds"
  )

  on.exit(
    if (file.exists(temporary_file)) {
      unlink(temporary_file)
    },
    add = TRUE
  )

  saveRDS(
    object,
    temporary_file,
    compress = "gzip"
  )

  verified <- tryCatch(
    {
      test_object <- readRDS(temporary_file)
      !is.null(test_object)
    },
    error = function(e) FALSE
  )

  if (!verified) {
    stop(
      "RDS verification failed for ",
      description,
      ": ",
      path
    )
  }

  copied <- file.copy(
    temporary_file,
    path,
    overwrite = TRUE
  )

  if (!copied) {
    stop(
      "Failed to copy verified RDS file to: ",
      path
    )
  }
}


collect_object_results <- function(object_name) {
  object_root <- file.path(
    inference_root,
    object_name
  )

  grouping_levels <- c(
    "celltype",
    "cluster",
    "overall"
  )

  grouping_results <- setNames(
    lapply(
      grouping_levels,
      function(grouping_level) {
        result_dir <- file.path(
          object_root,
          grouping_level
        )

        list(
          pseudobulk_metadata = read_optional_csv(
            file.path(result_dir, "pseudobulk_metadata.csv")
          ),
          pseudobulk_filtering = read_optional_csv(
            file.path(result_dir, "pseudobulk_filtering.csv")
          ),
          ulm_scores = read_optional_csv(
            file.path(result_dir, "ulm_scores.csv")
          ),
          ulm_adjusted_p = read_optional_csv(
            file.path(result_dir, "ulm_padj.csv")
          ),
          ulm_top_tfs = read_optional_csv(
            file.path(result_dir, "ulm_top_tfs.csv")
          ),
          mlm_scores = read_optional_csv(
            file.path(result_dir, "mlm_scores.csv")
          ),
          mlm_adjusted_p = read_optional_csv(
            file.path(result_dir, "mlm_padj.csv")
          ),
          mlm_top_tfs = read_optional_csv(
            file.path(result_dir, "mlm_top_tfs.csv")
          ),
          method_summary = read_optional_csv(
            file.path(result_dir, "method_summary.csv")
          )
        )
      }
    ),
    grouping_levels
  )

  grouping_results
}


manifest <- suppressMessages(
  readr::read_csv(
    manifest_file,
    show_col_types = FALSE
  )
)

analysis_config <- jsonlite::fromJSON(
  config_file,
  simplifyVector = TRUE
)

object_results <- setNames(
  lapply(
    manifest$object_name,
    collect_object_results
  ),
  manifest$object_name
)

global_results <- list(
  network_provenance = read_optional_csv(
    file.path(table_root, "TF_network_provenance.csv")
  ),
  inference_run_summary = read_optional_csv(
    file.path(table_root, "TF_inference_run_summary.csv")
  ),
  report_generation_summary = read_optional_csv(
    file.path(table_root, "TF_report_generation_summary.csv")
  ),
  sample_level_ulm_scores = read_optional_csv(
    file.path(table_root, "sample_level_ULM_scores.csv")
  ),
  sample_level_pseudobulk_metadata = read_optional_csv(
    file.path(table_root, "sample_level_pseudobulk_metadata.csv")
  ),
  sample_level_contrasts = read_optional_csv(
    file.path(table_root, "sample_level_contrasts_all.csv")
  ),
  combined_top_tfs = read_optional_csv(
    file.path(table_root, "combined_top_TFs.csv")
  ),
  recurrent_top_tfs = read_optional_csv(
    file.path(table_root, "recurrent_top_TFs.csv")
  )
)

result_bundle <- list(
  created_at = as.character(Sys.time()),
  analysis_description = paste(
    "Pseudobulk transcription-factor activity inferred with",
    "decoupler and the CollecTRI mouse regulatory network."
  ),
  config = analysis_config,
  manifest = manifest,
  global = global_results,
  objects = object_results
)

bundle_file <- file.path(
  output_root,
  "tf_activity_results_lightweight.rds"
)

safe_save_rds(
  result_bundle,
  bundle_file,
  "lightweight TF-analysis result bundle"
)

verification_rows <- list(
  tibble(
    file_type = "lightweight_results",
    file_name = basename(bundle_file),
    path = bundle_file,
    verified = TRUE
  )
)

for (i in seq_len(nrow(manifest))) {
  if (!as_flag(manifest$save_augmented_object[[i]])) {
    next
  }

  object_name <- manifest$object_name[[i]]
  input_file <- file.path(
    project_dir,
    manifest$file_path[[i]]
  )

  object <- readRDS(input_file)

  object@misc$tf_activity_analysis <- list(
    created_at = as.character(Sys.time()),
    object_name = object_name,
    condition_label = manifest$condition_label[[i]],
    config = analysis_config,
    results = object_results[[object_name]],
    global_context = global_results
  )

  output_file <- file.path(
    output_root,
    paste0(
      object_name,
      "_with_tf_activity_results.rds"
    )
  )

  safe_save_rds(
    object,
    output_file,
    paste0(
      object_name,
      " augmented Seurat object"
    )
  )

  verification_rows[[length(verification_rows) + 1]] <- tibble(
    file_type = "augmented_seurat",
    file_name = basename(output_file),
    path = output_file,
    verified = TRUE
  )
}

verification_table <- bind_rows(
  verification_rows
)

write_csv(
  verification_table,
  file.path(
    output_root,
    "tf_activity_result_manifest.csv"
  )
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(
    table_root,
    "10e_TF_result_packaging_sessionInfo.txt"
  )
)
