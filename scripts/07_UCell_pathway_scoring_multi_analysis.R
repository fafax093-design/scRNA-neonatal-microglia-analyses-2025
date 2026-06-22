# ==============================================================================
# 07. UCell pathway scoring and visualization across multiple analyses
#
# Description:
# This script calculates UCell pathway scores for one or more annotated Seurat
# objects defined in an external manifest. Each analysis may use a different
# Seurat object, pathway gene set, grouping variable, biological-sample column,
# group order, color palette, and statistical comparison plan.
#
# The script exports:
#   1. cell-level UCell scores;
#   2. sample-level score summaries;
#   3. cell-level descriptive plots;
#   4. sample-level inferential plots;
#   5. pairwise statistical-test results;
#   6. optionally, the Seurat object containing the UCell score.
#
# Statistical principle:
# Inferential comparisons are performed using independent biological samples,
# not individual cells. Cell-level plots are descriptive and do not display
# cell-level hypothesis tests.
#
# Author:
# Jinjin Zhu
#
# Configuration file:
#   config/07_UCell_analysis_manifest.csv
#
# Main output directories:
#   TABLE/UCell/<analysis_id>/
#   FIGURE/UCell/<analysis_id>/
#   DATA/OUTPUT/UCell/
# ==============================================================================


# 1. Load packages -------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(UCell)
  library(ggrastr)
  library(ggplot2)
  library(ggpubr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readxl)
})


# 2. Define project paths -------------------------------------------------------

analysis_manifest_file <- "config/07_UCell_analysis_manifest.csv"

output_table_root <- "TABLE/UCell"
output_figure_root <- "FIGURE/UCell"
output_object_root <- "DATA/OUTPUT/UCell"

dir.create(
  output_table_root,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  output_figure_root,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  output_object_root,
  recursive = TRUE,
  showWarnings = FALSE
)


# 3. Define global settings -----------------------------------------------------

random_seed <- 1234

set.seed(
  random_seed
)


# 4. Define helper functions ----------------------------------------------------

as_flag <- function(x) {

  x <- tolower(
    trimws(
      as.character(
        x
      )
    )
  )

  x %in% c(
    "true",
    "t",
    "1",
    "yes",
    "y"
  )
}


sanitize_name <- function(x) {

  x <- gsub(
    pattern = "[^A-Za-z0-9._-]+",
    replacement = "_",
    x = as.character(
      x
    )
  )

  gsub(
    pattern = "_+",
    replacement = "_",
    x = x
  )
}


read_external_table <- function(
  file,
  required = TRUE
) {

  file <- trimws(
    as.character(
      file
    )
  )

  if (is.na(file) ||
      file == "") {

    if (required) {
      stop(
        "A required file path is empty."
      )
    }

    return(
      NULL
    )
  }

  if (!file.exists(
    file
  )) {

    if (required) {
      stop(
        "Input file does not exist: ",
        file
      )
    }

    return(
      NULL
    )
  }

  extension <- tolower(
    tools::file_ext(
      file
    )
  )

  if (extension == "csv") {

    table <- read.csv(
      file,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

  } else if (extension %in% c(
    "tsv",
    "txt"
  )) {

    table <- read.delim(
      file,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

  } else if (extension %in% c(
    "xlsx",
    "xls"
  )) {

    table <- readxl::read_excel(
      file,
      sheet = 1
    ) %>%
      as.data.frame(
        stringsAsFactors = FALSE
      )

  } else {

    stop(
      "Unsupported table format: ",
      file
    )
  }

  table
}


read_pathway_gene_set <- function(
  file,
  pathway_name
) {

  file <- trimws(
    as.character(
      file
    )
  )

  if (!file.exists(
    file
  )) {

    stop(
      "Gene-set file does not exist: ",
      file
    )
  }

  extension <- tolower(
    tools::file_ext(
      file
    )
  )

  if (extension == "rds") {

    gene_set_object <- readRDS(
      file
    )

    if (is.list(
      gene_set_object
    ) &&
        !is.data.frame(
          gene_set_object
        )) {

      if (!pathway_name %in% names(
        gene_set_object
      )) {

        stop(
          "Pathway '",
          pathway_name,
          "' was not found in gene-set file: ",
          file
        )
      }

      genes <- gene_set_object[
        [pathway_name]
      ]

    } else if (is.data.frame(
      gene_set_object
    )) {

      gene_set_table <- gene_set_object

      if (!all(
        c(
          "pathway",
          "gene"
        ) %in% colnames(
          gene_set_table
        )
      )) {

        stop(
          "A data-frame RDS gene-set file must contain columns: pathway, gene."
        )
      }

      genes <- gene_set_table %>%
        filter(
          pathway == !!pathway_name
        ) %>%
        pull(
          gene
        )

    } else {

      stop(
        "Unsupported RDS gene-set structure in file: ",
        file
      )
    }

  } else {

    gene_set_table <- read_external_table(
      file = file,
      required = TRUE
    )

    if (!all(
      c(
        "pathway",
        "gene"
      ) %in% colnames(
        gene_set_table
      )
    )) {

      stop(
        "Tabular gene-set files must contain columns: pathway, gene."
      )
    }

    genes <- gene_set_table %>%
      filter(
        pathway == !!pathway_name
      ) %>%
      pull(
        gene
      )
  }

  genes <- unique(
    trimws(
      as.character(
        genes
      )
    )
  )

  genes <- genes[
    !is.na(
      genes
    ) &
      genes != ""
  ]

  if (length(
    genes
  ) == 0) {

    stop(
      "The selected pathway gene set is empty: ",
      pathway_name
    )
  }

  genes
}


read_group_palette <- function(
  file,
  observed_groups
) {

  palette_table <- read_external_table(
    file = file,
    required = TRUE
  )

  required_columns <- c(
    "group",
    "color"
  )

  missing_columns <- setdiff(
    required_columns,
    colnames(
      palette_table
    )
  )

  if (length(
    missing_columns
  ) > 0) {

    stop(
      "Group palette file is missing columns: ",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  }

  palette_table <- palette_table %>%
    transmute(
      group = trimws(
        as.character(
          group
        )
      ),
      color = trimws(
        as.character(
          color
        )
      )
    ) %>%
    filter(
      !is.na(
        group
      ),
      group != "",
      !is.na(
        color
      ),
      color != ""
    )

  if (anyDuplicated(
    palette_table$group
  )) {

    stop(
      "Duplicated groups were detected in the group palette file."
    )
  }

  missing_groups <- setdiff(
    unique(
      observed_groups
    ),
    palette_table$group
  )

  if (length(
    missing_groups
  ) > 0) {

    warning(
      "The following object groups are not included in the selected plotting ",
      "palette and will be excluded: ",
      paste(
        missing_groups,
        collapse = ", "
      )
    )
  }

  list(
    group_levels = palette_table$group,
    group_colors = setNames(
      palette_table$color,
      palette_table$group
    ),
    table = palette_table
  )
}


validate_manifest <- function(manifest) {

  required_columns <- c(
    "analysis_id",
    "input_file",
    "geneset_file",
    "pathway_name",
    "group_column",
    "sample_column",
    "assay",
    "group_palette_file",
    "comparison_file",
    "sample_summary_method",
    "test_method",
    "p_adjust_method",
    "run_cell_plot",
    "run_sample_plot",
    "save_scored_object",
    "plot_width_cm",
    "plot_height_cm"
  )

  missing_columns <- setdiff(
    required_columns,
    colnames(
      manifest
    )
  )

  if (length(
    missing_columns
  ) > 0) {

    stop(
      "The UCell analysis manifest is missing columns: ",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  }

  manifest <- manifest %>%
    mutate(
      across(
        c(
          analysis_id,
          input_file,
          geneset_file,
          pathway_name,
          group_column,
          sample_column,
          assay,
          group_palette_file,
          comparison_file,
          sample_summary_method,
          test_method,
          p_adjust_method
        ),
        ~ trimws(
          as.character(
            .x
          )
        )
      ),
      run_cell_plot = as_flag(
        run_cell_plot
      ),
      run_sample_plot = as_flag(
        run_sample_plot
      ),
      save_scored_object = as_flag(
        save_scored_object
      ),
      plot_width_cm = as.numeric(
        plot_width_cm
      ),
      plot_height_cm = as.numeric(
        plot_height_cm
      )
    )

  if (any(
    is.na(
      manifest$analysis_id
    ) |
      manifest$analysis_id == ""
  )) {

    stop(
      "Every manifest row must contain a non-empty analysis_id."
    )
  }

  if (anyDuplicated(
    manifest$analysis_id
  )) {

    stop(
      "Duplicated analysis_id values were detected."
    )
  }

  valid_summary_methods <- c(
    "mean",
    "median"
  )

  invalid_summary_methods <- setdiff(
    unique(
      manifest$sample_summary_method
    ),
    valid_summary_methods
  )

  if (length(
    invalid_summary_methods
  ) > 0) {

    stop(
      "Unsupported sample_summary_method values: ",
      paste(
        invalid_summary_methods,
        collapse = ", "
      )
    )
  }

  valid_test_methods <- c(
    "t.test",
    "wilcox.test"
  )

  invalid_test_methods <- setdiff(
    unique(
      manifest$test_method
    ),
    valid_test_methods
  )

  if (length(
    invalid_test_methods
  ) > 0) {

    stop(
      "Unsupported test_method values: ",
      paste(
        invalid_test_methods,
        collapse = ", "
      )
    )
  }

  manifest
}


calculate_pairwise_statistics <- function(
  sample_data,
  comparison_file,
  test_method,
  p_adjust_method
) {

  comparisons <- read_external_table(
    file = comparison_file,
    required = TRUE
  )

  required_columns <- c(
    "group1",
    "group2"
  )

  missing_columns <- setdiff(
    required_columns,
    colnames(
      comparisons
    )
  )

  if (length(
    missing_columns
  ) > 0) {

    stop(
      "Comparison file is missing columns: ",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  }

  comparisons <- comparisons %>%
    transmute(
      group1 = trimws(
        as.character(
          group1
        )
      ),
      group2 = trimws(
        as.character(
          group2
        )
      )
    ) %>%
    filter(
      group1 != "",
      group2 != ""
    ) %>%
    distinct()

  results <- lapply(
    seq_len(
      nrow(
        comparisons
      )
    ),
    function(i) {

      group1 <- comparisons$group1[
        [i]
      ]

      group2 <- comparisons$group2[
        [i]
      ]

      values1 <- sample_data %>%
        filter(
          group == group1
        ) %>%
        pull(
          sample_score
        )

      values2 <- sample_data %>%
        filter(
          group == group2
        ) %>%
        pull(
          sample_score
        )

      if (length(
        values1
      ) < 2 ||
          length(
            values2
          ) < 2) {

        return(
          tibble(
            group1 = group1,
            group2 = group2,
            n_group1 = length(
              values1
            ),
            n_group2 = length(
              values2
            ),
            statistic = NA_real_,
            p_value = NA_real_,
            error_message = "At least two biological samples are required per group."
          )
        )
      }

      test_result <- tryCatch(
        {
          if (test_method == "t.test") {

            stats::t.test(
              values1,
              values2,
              paired = FALSE
            )

          } else {

            stats::wilcox.test(
              values1,
              values2,
              paired = FALSE,
              exact = FALSE
            )
          }
        },
        error = function(e) {
          e
        }
      )

      if (inherits(
        test_result,
        "error"
      )) {

        tibble(
          group1 = group1,
          group2 = group2,
          n_group1 = length(
            values1
          ),
          n_group2 = length(
            values2
          ),
          statistic = NA_real_,
          p_value = NA_real_,
          error_message = conditionMessage(
            test_result
          )
        )

      } else {

        tibble(
          group1 = group1,
          group2 = group2,
          n_group1 = length(
            values1
          ),
          n_group2 = length(
            values2
          ),
          statistic = unname(
            test_result$statistic
          ),
          p_value = test_result$p.value,
          error_message = NA_character_
        )
      }
    }
  )

  statistics <- bind_rows(
    results
  )

  statistics$p_adjusted <- p.adjust(
    statistics$p_value,
    method = p_adjust_method
  )

  statistics$p_significance <- case_when(
    is.na(
      statistics$p_adjusted
    ) ~ "NA",
    statistics$p_adjusted <= 0.0001 ~ "****",
    statistics$p_adjusted <= 0.001 ~ "***",
    statistics$p_adjusted <= 0.01 ~ "**",
    statistics$p_adjusted <= 0.05 ~ "*",
    TRUE ~ "ns"
  )

  statistics
}


save_plot <- function(
  plot_object,
  output_dir,
  file_stem,
  width_cm,
  height_cm
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
    width = width_cm,
    height = height_cm,
    units = "cm",
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
    width = width_cm,
    height = height_cm,
    units = "cm",
    dpi = 600,
    compression = "lzw"
  )
}


# 5. Define the per-analysis UCell function ------------------------------------

run_one_ucell_analysis <- function(
  manifest_row
) {

  analysis_id <- manifest_row$analysis_id[
    [1]
  ]

  analysis_id_safe <- sanitize_name(
    analysis_id
  )

  table_dir <- file.path(
    output_table_root,
    analysis_id_safe
  )

  figure_dir <- file.path(
    output_figure_root,
    analysis_id_safe
  )

  dir.create(
    table_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  dir.create(
    figure_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  input_file <- manifest_row$input_file[
    [1]
  ]

  geneset_file <- manifest_row$geneset_file[
    [1]
  ]

  pathway_name <- manifest_row$pathway_name[
    [1]
  ]

  group_column <- manifest_row$group_column[
    [1]
  ]

  sample_column <- manifest_row$sample_column[
    [1]
  ]

  assay_name <- manifest_row$assay[
    [1]
  ]

  message(
    "\nRunning UCell analysis: ",
    analysis_id,
    " | pathway: ",
    pathway_name
  )

  if (!file.exists(
    input_file
  )) {

    stop(
      "Input Seurat file does not exist: ",
      input_file
    )
  }

  object <- readRDS(
    input_file
  )

  if (!inherits(
    object,
    "Seurat"
  )) {

    stop(
      "The input file does not contain a Seurat object: ",
      input_file
    )
  }

  required_metadata <- c(
    group_column,
    sample_column
  )

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
      "The Seurat object is missing metadata columns: ",
      paste(
        missing_metadata,
        collapse = ", "
      )
    )
  }

  if (!assay_name %in% Assays(
    object
  )) {

    stop(
      "Assay '",
      assay_name,
      "' was not found."
    )
  }

  pathway_genes <- read_pathway_gene_set(
    file = geneset_file,
    pathway_name = pathway_name
  )

  genes_present <- intersect(
    pathway_genes,
    rownames(
      object[
        [assay_name]
      ]
    )
  )

  genes_missing <- setdiff(
    pathway_genes,
    genes_present
  )

  if (length(
    genes_present
  ) == 0) {

    stop(
      "None of the pathway genes were found in assay: ",
      assay_name
    )
  }

  gene_coverage <- tibble(
    pathway = pathway_name,
    n_genes_in_gene_set = length(
      pathway_genes
    ),
    n_genes_present = length(
      genes_present
    ),
    n_genes_missing = length(
      genes_missing
    ),
    coverage_percentage = 100 *
      n_genes_present /
      n_genes_in_gene_set
  )

  write.csv(
    gene_coverage,
    file = file.path(
      table_dir,
      "UCell_gene_set_coverage.csv"
    ),
    row.names = FALSE
  )

  write.csv(
    tibble(
      gene = genes_missing
    ),
    file = file.path(
      table_dir,
      "UCell_missing_genes.csv"
    ),
    row.names = FALSE
  )

  DefaultAssay(
    object
  ) <- assay_name

  features <- setNames(
    list(
      genes_present
    ),
    pathway_name
  )

  metadata_columns_before <- colnames(
    object@meta.data
  )

  object <- UCell::AddModuleScore_UCell(
    obj = object,
    features = features
  )

  metadata_columns_after <- colnames(
    object@meta.data
  )

  new_score_columns <- setdiff(
    metadata_columns_after,
    metadata_columns_before
  )

  ucell_score_columns <- new_score_columns[
    grepl(
      "_UCell$",
      new_score_columns
    )
  ]

  if (length(
    ucell_score_columns
  ) != 1) {

    stop(
      "Expected exactly one new UCell score column, but detected: ",
      paste(
        ucell_score_columns,
        collapse = ", "
      )
    )
  }

  raw_score_column <- ucell_score_columns[
    [1]
  ]

  standardized_score_column <- paste0(
    sanitize_name(
      pathway_name
    ),
    "_UCell"
  )

  object[
    [standardized_score_column]
  ] <- object@meta.data[
    [raw_score_column]
  ]

  cell_scores <- object@meta.data %>%
    rownames_to_column(
      "cell_barcode"
    ) %>%
    transmute(
      cell_barcode = cell_barcode,
      sample_id = as.character(
        .data[
          [sample_column]
        ]
      ),
      group = as.character(
        .data[
          [group_column]
        ]
      ),
      score = as.numeric(
        .data[
          [standardized_score_column]
        ]
      )
    )

  palette_info <- read_group_palette(
    file = manifest_row$group_palette_file[
      [1]
    ],
    observed_groups = unique(
      cell_scores$group
    )
  )

  group_levels <- palette_info$group_levels
  group_colors <- palette_info$group_colors

  cell_scores <- cell_scores %>%
    filter(
      group %in% group_levels
    ) %>%
    mutate(
      group = factor(
        group,
        levels = group_levels
      )
    )

  if (nrow(
    cell_scores
  ) == 0) {

    stop(
      "No cells remained after filtering to the configured groups."
    )
  }

  sample_group_check <- cell_scores %>%
    distinct(
      sample_id,
      group
    ) %>%
    count(
      sample_id,
      name = "n_groups"
    ) %>%
    filter(
      n_groups > 1
    )

  if (nrow(
    sample_group_check
  ) > 0) {

    stop(
      "At least one biological sample maps to more than one group."
    )
  }

  write.csv(
    cell_scores,
    file = file.path(
      table_dir,
      "UCell_scores_by_cell.csv"
    ),
    row.names = FALSE
  )

  write.csv(
    palette_info$table,
    file = file.path(
      table_dir,
      "UCell_group_palette_used.csv"
    ),
    row.names = FALSE
  )

  sample_summary_method <- manifest_row$sample_summary_method[
    [1]
  ]

  sample_scores <- cell_scores %>%
    group_by(
      sample_id,
      group
    ) %>%
    summarise(
      n_cells = n(),
      mean_score = mean(
        score,
        na.rm = TRUE
      ),
      median_score = median(
        score,
        na.rm = TRUE
      ),
      sd_score = sd(
        score,
        na.rm = TRUE
      ),
      .groups = "drop"
    ) %>%
    mutate(
      sample_score = if (
        sample_summary_method == "mean"
      ) {
        mean_score
      } else {
        median_score
      },
      group = factor(
        group,
        levels = group_levels
      )
    )

  write.csv(
    sample_scores,
    file = file.path(
      table_dir,
      "UCell_scores_by_sample.csv"
    ),
    row.names = FALSE
  )

  group_summary <- sample_scores %>%
    group_by(
      group
    ) %>%
    summarise(
      n_samples = n(),
      mean_sample_score = mean(
        sample_score,
        na.rm = TRUE
      ),
      sd_sample_score = sd(
        sample_score,
        na.rm = TRUE
      ),
      median_sample_score = median(
        sample_score,
        na.rm = TRUE
      ),
      q1_sample_score = quantile(
        sample_score,
        0.25,
        na.rm = TRUE
      ),
      q3_sample_score = quantile(
        sample_score,
        0.75,
        na.rm = TRUE
      ),
      .groups = "drop"
    )

  write.csv(
    group_summary,
    file = file.path(
      table_dir,
      "UCell_group_summary.csv"
    ),
    row.names = FALSE
  )

  statistics <- calculate_pairwise_statistics(
    sample_data = sample_scores,
    comparison_file = manifest_row$comparison_file[
      [1]
    ],
    test_method = manifest_row$test_method[
      [1]
    ],
    p_adjust_method = manifest_row$p_adjust_method[
      [1]
    ]
  )

  score_range <- range(
    sample_scores$sample_score,
    na.rm = TRUE
  )

  score_span <- diff(
    score_range
  )

  if (!is.finite(
    score_span
  ) ||
      score_span == 0) {

    score_span <- 0.05
  }

  statistics <- statistics %>%
    mutate(
      y.position =
        max(
          sample_scores$sample_score,
          na.rm = TRUE
        ) +
        seq_len(
          n()
        ) *
        score_span *
        0.12
    )

  write.csv(
    statistics,
    file = file.path(
      table_dir,
      "UCell_pairwise_statistics.csv"
    ),
    row.names = FALSE
  )


  # 5.1 Cell-level descriptive plot -------------------------------------------

  if (manifest_row$run_cell_plot[
    [1]
  ]) {

    cell_counts <- cell_scores %>%
      count(
        group,
        name = "n_cells"
      )

    label_y <- max(
      cell_scores$score,
      na.rm = TRUE
    ) +
      0.03

    cell_counts$y <- label_y
    cell_counts$label <- paste0(
      "cells=",
      cell_counts$n_cells
    )

    p_cell <- ggplot(
      cell_scores,
      aes(
        x = group,
        y = score,
        fill = group,
        color = group
      )
    ) +
      ggrastr::geom_jitter_rast(
        color = "#00000033",
        pch = 19,
        size = 0.8,
        stroke = 0.01,
        position = position_jitter(
          width = 0.15,
          seed = random_seed
        ),
        alpha = 0.3
      ) +
      geom_boxplot(
        color = "black",
        outlier.shape = NA,
        alpha = 0.8,
        linewidth = 0.5,
        width = 0.4
      ) +
      geom_text(
        data = cell_counts,
        aes(
          x = group,
          y = y,
          label = label
        ),
        inherit.aes = FALSE,
        size = 2.5
      ) +
      scale_fill_manual(
        values = group_colors,
        limits = group_levels,
        breaks = group_levels,
        drop = FALSE
      ) +
      scale_color_manual(
        values = group_colors,
        limits = group_levels,
        breaks = group_levels,
        drop = FALSE
      ) +
      labs(
        x = NULL,
        y = paste0(
          pathway_name,
          " UCell score"
        ),
        title = NULL
      ) +
      theme_minimal() +
      theme(
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(
          color = "black"
        ),
        axis.ticks = element_line(
          color = "black",
          linewidth = 0.5
        ),
        axis.text.x = element_text(
          size = 6,
          color = "black"
        ),
        axis.text.y = element_text(
          size = 6,
          color = "black"
        ),
        axis.title.y = element_text(
          color = "black",
          size = 8,
          face = "bold"
        ),
        legend.position = "none"
      )

    save_plot(
      plot_object = p_cell,
      output_dir = figure_dir,
      file_stem = "UCell_cell_level_distribution",
      width_cm = manifest_row$plot_width_cm[
        [1]
      ],
      height_cm = manifest_row$plot_height_cm[
        [1]
      ]
    )
  }


  # 5.2 Sample-level inferential plot -----------------------------------------

  if (manifest_row$run_sample_plot[
    [1]
  ]) {

    sample_counts <- sample_scores %>%
      count(
        group,
        name = "n_samples"
      )

    sample_label_y <- max(
      c(
        sample_scores$sample_score,
        statistics$y.position
      ),
      na.rm = TRUE
    ) +
      score_span *
      0.08

    sample_counts$y <- sample_label_y
    sample_counts$label <- paste0(
      "n=",
      sample_counts$n_samples
    )

    p_sample <- ggplot(
      sample_scores,
      aes(
        x = group,
        y = sample_score,
        fill = group,
        color = group
      )
    ) +
      geom_boxplot(
        color = "black",
        outlier.shape = NA,
        alpha = 0.5,
        linewidth = 0.5,
        width = 0.45
      ) +
      geom_point(
        shape = 21,
        size = 2.2,
        stroke = 0.5,
        position = position_jitter(
          width = 0.08,
          seed = random_seed
        )
      ) +
      geom_text(
        data = sample_counts,
        aes(
          x = group,
          y = y,
          label = label
        ),
        inherit.aes = FALSE,
        size = 2.5
      ) +
      ggpubr::stat_pvalue_manual(
        statistics,
        label = "p_significance",
        xmin = "group1",
        xmax = "group2",
        y.position = "y.position",
        tip.length = 0.01,
        size = 3,
        hide.ns = FALSE
      ) +
      scale_fill_manual(
        values = group_colors,
        limits = group_levels,
        breaks = group_levels,
        drop = FALSE
      ) +
      scale_color_manual(
        values = group_colors,
        limits = group_levels,
        breaks = group_levels,
        drop = FALSE
      ) +
      labs(
        x = NULL,
        y = paste0(
          pathway_name,
          " UCell score\n(",
          sample_summary_method,
          " per biological sample)"
        ),
        title = NULL
      ) +
      theme_minimal() +
      theme(
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(
          color = "black"
        ),
        axis.ticks = element_line(
          color = "black",
          linewidth = 0.5
        ),
        axis.text.x = element_text(
          size = 6,
          color = "black"
        ),
        axis.text.y = element_text(
          size = 6,
          color = "black"
        ),
        axis.title.y = element_text(
          color = "black",
          size = 8,
          face = "bold"
        ),
        legend.position = "none"
      )

    save_plot(
      plot_object = p_sample,
      output_dir = figure_dir,
      file_stem = "UCell_sample_level_comparison",
      width_cm = manifest_row$plot_width_cm[
        [1]
      ],
      height_cm = manifest_row$plot_height_cm[
        [1]
      ]
    )
  }


  # 5.3 Save scored object and parameters -------------------------------------

  if (manifest_row$save_scored_object[
    [1]
  ]) {

    output_object_file <- file.path(
      output_object_root,
      paste0(
        analysis_id_safe,
        "_UCell_scored.rds"
      )
    )

    saveRDS(
      object,
      file = output_object_file,
      compress = "gzip"
    )

  } else {

    output_object_file <- NA_character_
  }

  analysis_parameters <- tibble(
    parameter = c(
      "analysis_id",
      "input_file",
      "geneset_file",
      "pathway_name",
      "group_column",
      "sample_column",
      "assay",
      "score_column",
      "sample_summary_method",
      "test_method",
      "p_adjust_method",
      "random_seed",
      "saved_scored_object"
    ),
    value = c(
      analysis_id,
      input_file,
      geneset_file,
      pathway_name,
      group_column,
      sample_column,
      assay_name,
      standardized_score_column,
      sample_summary_method,
      manifest_row$test_method[
        [1]
      ],
      manifest_row$p_adjust_method[
        [1]
      ],
      as.character(
        random_seed
      ),
      output_object_file
    )
  )

  write.csv(
    analysis_parameters,
    file = file.path(
      table_dir,
      "UCell_analysis_parameters.csv"
    ),
    row.names = FALSE
  )

  tibble(
    analysis_id = analysis_id,
    status = "completed",
    input_file = input_file,
    pathway_name = pathway_name,
    n_cells = nrow(
      cell_scores
    ),
    n_samples = n_distinct(
      sample_scores$sample_id
    ),
    n_groups = n_distinct(
      sample_scores$group
    ),
    n_pathway_genes = length(
      pathway_genes
    ),
    n_pathway_genes_present = length(
      genes_present
    )
  )
}


# 6. Load and validate the UCell analysis manifest -----------------------------

if (!file.exists(
  analysis_manifest_file
)) {

  stop(
    "UCell analysis manifest does not exist: ",
    analysis_manifest_file,
    "\nCreate it from the provided template before running this script."
  )
}

analysis_manifest <- read.csv(
  analysis_manifest_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

analysis_manifest <- validate_manifest(
  analysis_manifest
)


# 7. Run all configured UCell analyses -----------------------------------------

analysis_results <- vector(
  mode = "list",
  length = nrow(
    analysis_manifest
  )
)

for (i in seq_len(
  nrow(
    analysis_manifest
  )
)) {

  manifest_row <- analysis_manifest[
    i,
    ,
    drop = FALSE
  ]

  analysis_id <- manifest_row$analysis_id[
    [1]
  ]

  analysis_results[
    [i]
  ] <- tryCatch(
    run_one_ucell_analysis(
      manifest_row = manifest_row
    ),
    error = function(e) {

      message(
        "UCell analysis failed: ",
        analysis_id,
        "\nReason: ",
        conditionMessage(
          e
        )
      )

      tibble(
        analysis_id = analysis_id,
        status = "failed",
        input_file = manifest_row$input_file[
          [1]
        ],
        pathway_name = manifest_row$pathway_name[
          [1]
        ],
        n_cells = NA_integer_,
        n_samples = NA_integer_,
        n_groups = NA_integer_,
        n_pathway_genes = NA_integer_,
        n_pathway_genes_present = NA_integer_,
        error_message = conditionMessage(
          e
        )
      )
    }
  )
}

analysis_summary <- bind_rows(
  analysis_results
)

write.csv(
  analysis_summary,
  file = file.path(
    output_table_root,
    "07_UCell_analysis_summary.csv"
  ),
  row.names = FALSE
)


# 8. Save software information -------------------------------------------------

writeLines(
  capture.output(
    sessionInfo()
  ),
  con = file.path(
    output_table_root,
    "07_UCell_sessionInfo.txt"
  )
)

message(
  "Configured UCell analyses have finished."
)

print(
  analysis_summary
)
