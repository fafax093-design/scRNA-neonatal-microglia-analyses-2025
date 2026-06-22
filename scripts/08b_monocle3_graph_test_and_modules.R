# ==============================================================================
# 08b. Monocle 3 graph test and pseudotime-associated gene modules
#
# Description:
# This script performs graph_test() on one or more ordered Monocle 3 CDS objects,
# selects genes using q-value and Moran's I thresholds, identifies gene modules,
# optionally merges modules according to a reviewed mapping table, aggregates
# module expression by a metadata-defined cell group, and generates a module
# heatmap.
#
# Configuration:
#   config/monocle3/08b_graph_test_manifest.csv
# ==============================================================================

suppressPackageStartupMessages({
  library(monocle3)
  library(dplyr)
  library(tibble)
  library(pheatmap)
  library(RColorBrewer)
  library(grid)
})

manifest_file <- "config/monocle3/08b_graph_test_manifest.csv"

output_data_root <- "DATA/OUTPUT/Monocle3"
output_table_root <- "TABLE/Monocle3"
output_figure_root <- "FIGURE/Monocle3"

dir.create(output_data_root, recursive = TRUE, showWarnings = FALSE)
dir.create(output_table_root, recursive = TRUE, showWarnings = FALSE)
dir.create(output_figure_root, recursive = TRUE, showWarnings = FALSE)

random_seed <- 1234
set.seed(random_seed)


sanitize_name <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", as.character(x))
  gsub("_+", "_", x)
}


read_optional_table <- function(file) {
  file <- trimws(as.character(file))

  if (is.na(file) || file == "" || !file.exists(file)) {
    return(NULL)
  }

  extension <- tolower(tools::file_ext(file))

  if (extension == "csv") {
    return(
      read.csv(
        file,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    )
  }

  if (extension %in% c("tsv", "txt")) {
    return(
      read.delim(
        file,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    )
  }

  stop("Optional table must be CSV, TSV, or TXT: ", file)
}


apply_module_mapping <- function(module_table, mapping_file) {
  mapping <- read_optional_table(mapping_file)

  if (is.null(mapping)) {
    module_table$module_original <- module_table$module
    return(module_table)
  }

  if (!all(c("old_module", "new_module") %in% colnames(mapping))) {
    stop("Module mapping file must contain old_module and new_module.")
  }

  mapping <- mapping |>
    transmute(
      old_module = as.character(old_module),
      new_module = as.character(new_module)
    )

  module_table <- module_table |>
    mutate(
      module_original = as.character(module),
      module = mapping$new_module[
        match(as.character(module), mapping$old_module)
      ],
      module = ifelse(
        is.na(module),
        module_original,
        module
      )
    )

  module_table
}


draw_pheatmap <- function(
  matrix,
  output_file,
  width,
  height,
  device = c("pdf", "tiff")
) {
  device <- match.arg(device)

  p <- pheatmap::pheatmap(
    matrix,
    scale = "row",
    cluster_cols = FALSE,
    cluster_rows = TRUE,
    show_rownames = TRUE,
    show_colnames = TRUE,
    border_color = "white",
    color = rev(RColorBrewer::brewer.pal(10, "RdBu")),
    silent = TRUE
  )

  if (device == "pdf") {
    pdf(output_file, width = width, height = height, useDingbats = FALSE)
  } else {
    tiff(
      output_file,
      width = width,
      height = height,
      units = "in",
      res = 600,
      compression = "lzw"
    )
  }

  grid::grid.newpage()
  grid::grid.draw(p$gtable)
  dev.off()
}


run_one_module_analysis <- function(row) {
  analysis_id <- row$analysis_id[[1]]
  analysis_id_safe <- sanitize_name(analysis_id)

  table_dir <- file.path(output_table_root, analysis_id_safe)
  figure_dir <- file.path(output_figure_root, analysis_id_safe)

  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

  cds_file <- row$cds_file[[1]]

  if (!file.exists(cds_file)) {
    stop("CDS file does not exist: ", cds_file)
  }

  cds <- readRDS(cds_file)

  if (!inherits(cds, "cell_data_set")) {
    stop("Input file does not contain a Monocle 3 cell_data_set.")
  }

  cell_group_column <- row$cell_group_column[[1]]

  if (!cell_group_column %in% colnames(colData(cds))) {
    stop("Cell-group column was not found: ", cell_group_column)
  }

  message("\nRunning graph_test(): ", analysis_id)

  graph_test_result <- graph_test(
    cds,
    neighbor_graph = "principal_graph",
    cores = row$cores[[1]],
    verbose = TRUE
  )

  graph_test_result <- as.data.frame(graph_test_result)

  write.csv(
    graph_test_result,
    file = file.path(table_dir, "graph_test_complete.csv"),
    row.names = TRUE
  )

  required_columns <- c("q_value", "morans_I")
  missing_columns <- setdiff(required_columns, colnames(graph_test_result))

  if (length(missing_columns) > 0) {
    stop(
      "graph_test result is missing columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  selected_gene_ids <- rownames(graph_test_result)[
    !is.na(graph_test_result$q_value) &
      !is.na(graph_test_result$morans_I) &
      graph_test_result$q_value < row$q_value_cutoff[[1]] &
      graph_test_result$morans_I > row$morans_I_cutoff[[1]]
  ]

  exclude_regex <- trimws(as.character(row$exclude_regex[[1]]))

  if (!is.na(exclude_regex) && exclude_regex != "") {
    gene_names <- rowData(cds)[selected_gene_ids, "gene_short_name"]
    excluded <- grepl(exclude_regex, gene_names)
    excluded[is.na(excluded)] <- FALSE
    selected_gene_ids <- selected_gene_ids[!excluded]
  }

  if (length(selected_gene_ids) == 0) {
    stop("No genes passed the graph_test filtering criteria.")
  }

  selected_result <- graph_test_result[
    selected_gene_ids,
    ,
    drop = FALSE
  ] |>
    rownames_to_column("gene_id") |>
    arrange(desc(morans_I))

  write.csv(
    selected_result,
    file = file.path(table_dir, "graph_test_selected_genes.csv"),
    row.names = FALSE
  )

  set.seed(random_seed)

  gene_module_df <- find_gene_modules(
    cds[selected_gene_ids, ],
    resolution = row$module_resolution[[1]],
    cores = row$module_cores[[1]]
  )

  gene_module_df <- as.data.frame(gene_module_df)
  gene_module_df <- apply_module_mapping(
    gene_module_df,
    row$module_merge_file[[1]]
  )

  write.csv(
    gene_module_df,
    file = file.path(table_dir, "gene_module_assignments.csv"),
    row.names = FALSE
  )

  selected_result$module <- gene_module_df$module[
    match(selected_result$gene_id, gene_module_df$id)
  ]

  if ("module_original" %in% colnames(gene_module_df)) {
    selected_result$module_original <-
      gene_module_df$module_original[
        match(selected_result$gene_id, gene_module_df$id)
      ]
  }

  write.csv(
    selected_result,
    file = file.path(table_dir, "graph_test_selected_genes_with_modules.csv"),
    row.names = FALSE
  )

  module_for_aggregation <- gene_module_df |>
    select(id, module)

  cell_group_df <- tibble(
    cell = rownames(colData(cds)),
    cell_group = as.character(
      colData(cds)[[cell_group_column]]
    )
  )

  aggregate_matrix <- aggregate_gene_expression(
    cds,
    module_for_aggregation,
    cell_group_df,
    scale_agg_values = FALSE
  )

  rownames(aggregate_matrix) <- paste0(
    "Module ",
    rownames(aggregate_matrix)
  )

  group_order_table <- read_optional_table(
    row$group_order_file[[1]]
  )

  if (!is.null(group_order_table)) {
    if (!"group" %in% colnames(group_order_table)) {
      stop("Group-order file must contain a column named group.")
    }

    requested_order <- as.character(group_order_table$group)
    requested_order <- requested_order[
      requested_order %in% colnames(aggregate_matrix)
    ]

    remaining_groups <- setdiff(
      colnames(aggregate_matrix),
      requested_order
    )

    aggregate_matrix <- aggregate_matrix[
      ,
      c(requested_order, remaining_groups),
      drop = FALSE
    ]
  }

  write.csv(
    aggregate_matrix,
    file = file.path(table_dir, "module_expression_by_cell_group.csv"),
    row.names = TRUE
  )

  draw_pheatmap(
    aggregate_matrix,
    file.path(figure_dir, "module_expression_heatmap.pdf"),
    row$heatmap_width[[1]],
    row$heatmap_height[[1]],
    "pdf"
  )

  draw_pheatmap(
    aggregate_matrix,
    file.path(figure_dir, "module_expression_heatmap.tiff"),
    row$heatmap_width[[1]],
    row$heatmap_height[[1]],
    "tiff"
  )

  saveRDS(
    gene_module_df,
    file = file.path(
      output_data_root,
      paste0(analysis_id_safe, "_gene_modules.rds")
    ),
    compress = "gzip"
  )

  parameters <- tibble(
    parameter = c(
      "analysis_id",
      "cds_file",
      "cell_group_column",
      "q_value_cutoff",
      "morans_I_cutoff",
      "exclude_regex",
      "module_resolution",
      "graph_test_cores",
      "module_cores",
      "random_seed",
      "n_selected_genes",
      "n_modules"
    ),
    value = c(
      analysis_id,
      cds_file,
      cell_group_column,
      as.character(row$q_value_cutoff[[1]]),
      as.character(row$morans_I_cutoff[[1]]),
      exclude_regex,
      as.character(row$module_resolution[[1]]),
      as.character(row$cores[[1]]),
      as.character(row$module_cores[[1]]),
      as.character(random_seed),
      as.character(length(selected_gene_ids)),
      as.character(dplyr::n_distinct(gene_module_df$module))
    )
  )

  write.csv(
    parameters,
    file = file.path(table_dir, "graph_test_module_parameters.csv"),
    row.names = FALSE
  )

  tibble(
    analysis_id = analysis_id,
    status = "completed",
    cds_file = cds_file,
    n_selected_genes = length(selected_gene_ids),
    n_modules = dplyr::n_distinct(gene_module_df$module)
  )
}


if (!file.exists(manifest_file)) {
  stop("Graph-test manifest does not exist: ", manifest_file)
}

manifest <- read.csv(
  manifest_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
) |>
  mutate(
    across(
      c(
        analysis_id,
        cds_file,
        cell_group_column,
        exclude_regex,
        module_merge_file,
        group_order_file
      ),
      ~ trimws(as.character(.x))
    ),
    q_value_cutoff = as.numeric(q_value_cutoff),
    morans_I_cutoff = as.numeric(morans_I_cutoff),
    module_resolution = as.numeric(module_resolution),
    cores = as.integer(cores),
    module_cores = as.integer(module_cores),
    heatmap_width = as.numeric(heatmap_width),
    heatmap_height = as.numeric(heatmap_height)
  )

required_columns <- c(
  "analysis_id",
  "cds_file",
  "cell_group_column",
  "q_value_cutoff",
  "morans_I_cutoff",
  "exclude_regex",
  "module_resolution",
  "module_merge_file",
  "group_order_file",
  "cores",
  "module_cores",
  "heatmap_width",
  "heatmap_height"
)

missing_columns <- setdiff(required_columns, colnames(manifest))

if (length(missing_columns) > 0) {
  stop(
    "Graph-test manifest is missing columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

results <- vector("list", nrow(manifest))

for (i in seq_len(nrow(manifest))) {
  current_id <- manifest$analysis_id[[i]]

  results[[i]] <- tryCatch(
    run_one_module_analysis(manifest[i, , drop = FALSE]),
    error = function(e) {
      message(
        "Graph-test analysis failed: ",
        current_id,
        "\nReason: ",
        conditionMessage(e)
      )

      tibble(
        analysis_id = current_id,
        status = "failed",
        cds_file = manifest$cds_file[[i]],
        n_selected_genes = NA_integer_,
        n_modules = NA_integer_,
        error_message = conditionMessage(e)
      )
    }
  )
}

summary_table <- bind_rows(results)

write.csv(
  summary_table,
  file = file.path(output_table_root, "08b_graph_test_summary.csv"),
  row.names = FALSE
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(output_table_root, "08b_graph_test_sessionInfo.txt")
)

print(summary_table)
