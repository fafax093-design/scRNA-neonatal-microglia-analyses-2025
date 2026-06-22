# ==============================================================================
# 08a. Monocle 3 trajectory construction and pseudotime assignment
#
# Description:
# This script builds Monocle 3 cell_data_set objects for one or more Seurat
# objects defined in an external manifest. It imports the selected Seurat UMAP,
# learns a principal graph, assigns a root from a reviewed cell-state population,
# calculates pseudotime, transfers pseudotime back to the Seurat object, exports
# cell-level results, and saves trajectory figures.
#
# Author:
# Jinjin Zhu
#
# Configuration:
#   config/monocle3/08a_trajectory_manifest.csv
#
# Main outputs:
#   DATA/OUTPUT/Monocle3/<analysis_id>_cds.rds
#   DATA/OUTPUT/Monocle3/<analysis_id>_Seurat_with_pseudotime.rds
#   TABLE/Monocle3/<analysis_id>/
#   FIGURE/Monocle3/<analysis_id>/
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(monocle3)
  library(SingleCellExperiment)
  library(dplyr)
  library(tibble)
  library(ggplot2)
})

manifest_file <- "config/monocle3/08a_trajectory_manifest.csv"

output_data_root <- "DATA/OUTPUT/Monocle3"
output_table_root <- "TABLE/Monocle3"
output_figure_root <- "FIGURE/Monocle3"

dir.create(output_data_root, recursive = TRUE, showWarnings = FALSE)
dir.create(output_table_root, recursive = TRUE, showWarnings = FALSE)
dir.create(output_figure_root, recursive = TRUE, showWarnings = FALSE)

random_seed <- 1234
set.seed(random_seed)


as_flag <- function(x) {
  tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes", "y")
}


sanitize_name <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", as.character(x))
  gsub("_+", "_", x)
}


split_values <- function(x) {
  x <- trimws(as.character(x))
  if (is.na(x) || x == "") {
    return(character(0))
  }
  values <- unlist(strsplit(x, ";", fixed = TRUE))
  trimws(values[values != ""])
}


read_palette <- function(file, observed_levels) {
  if (is.na(file) || trimws(file) == "" || !file.exists(file)) {
    colors <- grDevices::hcl.colors(
      length(observed_levels),
      palette = "Dark 3"
    )
    return(setNames(colors, observed_levels))
  }

  palette_table <- read.csv(
    file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (!all(c("label", "color") %in% colnames(palette_table))) {
    stop("Palette file must contain columns: label and color.")
  }

  palette_table <- palette_table |>
    transmute(
      label = trimws(as.character(label)),
      color = trimws(as.character(color))
    ) |>
    filter(label != "", color != "")

  missing_levels <- setdiff(observed_levels, palette_table$label)

  if (length(missing_levels) > 0) {
    stop(
      "Palette file is missing labels: ",
      paste(missing_levels, collapse = ", ")
    )
  }

  setNames(
    palette_table$color[match(observed_levels, palette_table$label)],
    observed_levels
  )
}


save_plot <- function(plot_object, output_dir, file_stem, width, height) {
  ggsave(
    filename = file.path(output_dir, paste0(file_stem, ".pdf")),
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    useDingbats = FALSE
  )

  ggsave(
    filename = file.path(output_dir, paste0(file_stem, ".tiff")),
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    dpi = 600,
    compression = "lzw"
  )
}


get_count_matrix <- function(object, assay_name) {
  if (!assay_name %in% Assays(object)) {
    stop("Count assay was not found: ", assay_name)
  }

  layer_names <- tryCatch(
    Layers(object = object, assay = assay_name),
    error = function(e) character(0)
  )

  count_layers <- grep("^counts", layer_names, value = TRUE)

  if (length(count_layers) > 1) {
    object <- JoinLayers(object = object, assay = assay_name)
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

  list(object = object, counts = counts)
}


get_root_principal_node <- function(cds, root_column, root_values) {
  if (!root_column %in% colnames(colData(cds))) {
    stop("Root metadata column was not found: ", root_column)
  }

  root_cells <- which(
    as.character(colData(cds)[[root_column]]) %in% root_values
  )

  if (length(root_cells) == 0) {
    stop(
      "No cells matched the requested root population: ",
      paste(root_values, collapse = ", ")
    )
  }

  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex

  closest_vertex <- as.matrix(
    closest_vertex[colnames(cds), , drop = FALSE]
  )

  closest_vertex_for_roots <- closest_vertex[root_cells, 1]

  most_common_vertex <- names(
    which.max(
      table(closest_vertex_for_roots)
    )
  )

  graph_vertex_names <- igraph::V(
    principal_graph(cds)[["UMAP"]]
  )$name

  root_pr_node <- graph_vertex_names[
    as.numeric(most_common_vertex)
  ]

  if (length(root_pr_node) != 1 || is.na(root_pr_node)) {
    stop("A unique principal-graph root node could not be determined.")
  }

  root_pr_node
}


validate_manifest <- function(manifest) {
  required_columns <- c(
    "analysis_id",
    "input_file",
    "count_assay",
    "seurat_reduction",
    "celltype_column",
    "group_column",
    "root_column",
    "root_values",
    "num_dim",
    "cores",
    "use_partition",
    "close_loop",
    "palette_file",
    "plot_width",
    "plot_height"
  )

  missing_columns <- setdiff(required_columns, colnames(manifest))

  if (length(missing_columns) > 0) {
    stop(
      "Trajectory manifest is missing columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  manifest |>
    mutate(
      across(
        c(
          analysis_id,
          input_file,
          count_assay,
          seurat_reduction,
          celltype_column,
          group_column,
          root_column,
          root_values,
          palette_file
        ),
        ~ trimws(as.character(.x))
      ),
      num_dim = as.integer(num_dim),
      cores = as.integer(cores),
      use_partition = as_flag(use_partition),
      close_loop = as_flag(close_loop),
      plot_width = as.numeric(plot_width),
      plot_height = as.numeric(plot_height)
    )
}


run_one_trajectory <- function(manifest_row) {
  analysis_id <- manifest_row$analysis_id[[1]]
  analysis_id_safe <- sanitize_name(analysis_id)

  input_file <- manifest_row$input_file[[1]]
  count_assay <- manifest_row$count_assay[[1]]
  seurat_reduction <- manifest_row$seurat_reduction[[1]]
  celltype_column <- manifest_row$celltype_column[[1]]
  group_column <- manifest_row$group_column[[1]]
  root_column <- manifest_row$root_column[[1]]
  root_values <- split_values(manifest_row$root_values[[1]])

  table_dir <- file.path(output_table_root, analysis_id_safe)
  figure_dir <- file.path(output_figure_root, analysis_id_safe)

  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

  message("\nBuilding Monocle 3 trajectory: ", analysis_id)

  if (!file.exists(input_file)) {
    stop("Input Seurat file does not exist: ", input_file)
  }

  object <- readRDS(input_file)

  if (!inherits(object, "Seurat")) {
    stop("Input file does not contain a Seurat object: ", input_file)
  }

  required_metadata <- unique(
    c(celltype_column, group_column, root_column)
  )

  missing_metadata <- setdiff(
    required_metadata,
    colnames(object@meta.data)
  )

  if (length(missing_metadata) > 0) {
    stop(
      "Required metadata columns are missing: ",
      paste(missing_metadata, collapse = ", ")
    )
  }

  if (!seurat_reduction %in% Reductions(object)) {
    stop("Seurat reduction was not found: ", seurat_reduction)
  }

  if (length(root_values) == 0) {
    stop("At least one root value must be specified.")
  }

  count_result <- get_count_matrix(object, count_assay)
  object <- count_result$object
  counts <- count_result$counts

  metadata <- object@meta.data[
    colnames(counts),
    ,
    drop = FALSE
  ]

  if (!identical(rownames(metadata), colnames(counts))) {
    stop("Cell metadata and count matrix are not aligned.")
  }

  gene_metadata <- data.frame(
    gene_short_name = rownames(counts),
    row.names = rownames(counts),
    stringsAsFactors = FALSE
  )

  cds <- new_cell_data_set(
    expression_data = counts,
    cell_metadata = metadata,
    gene_metadata = gene_metadata
  )

  cds <- preprocess_cds(
    cds,
    num_dim = manifest_row$num_dim[[1]]
  )

  cds <- reduce_dimension(
    cds,
    preprocess_method = "PCA",
    reduction_method = "UMAP",
    umap.fast_sgd = FALSE,
    cores = manifest_row$cores[[1]]
  )

  seurat_umap <- Embeddings(
    object = object,
    reduction = seurat_reduction
  )

  missing_cells <- setdiff(colnames(cds), rownames(seurat_umap))

  if (length(missing_cells) > 0) {
    stop(
      "The Seurat UMAP is missing cells present in the CDS: ",
      paste(head(missing_cells, 10), collapse = ", ")
    )
  }

  seurat_umap <- seurat_umap[
    colnames(cds),
    1:2,
    drop = FALSE
  ]

  colnames(seurat_umap) <- c("UMAP_1", "UMAP_2")

  reducedDims(cds)$UMAP <- seurat_umap

  cds <- cluster_cells(
    cds,
    reduction_method = "UMAP"
  )

  cds <- learn_graph(
    cds,
    use_partition = manifest_row$use_partition[[1]],
    close_loop = manifest_row$close_loop[[1]],
    verbose = TRUE
  )

  root_pr_node <- get_root_principal_node(
    cds = cds,
    root_column = root_column,
    root_values = root_values
  )

  cds <- order_cells(
    cds,
    root_pr_nodes = root_pr_node
  )

  pseudotime_values <- monocle3::pseudotime(cds)

  pseudotime_table <- tibble(
    cell_barcode = names(pseudotime_values),
    pseudotime = as.numeric(pseudotime_values),
    finite_pseudotime = is.finite(pseudotime_values),
    celltype = as.character(
      colData(cds)[names(pseudotime_values), celltype_column]
    ),
    group = as.character(
      colData(cds)[names(pseudotime_values), group_column]
    )
  )

  write.csv(
    pseudotime_table,
    file = file.path(table_dir, "pseudotime_by_cell.csv"),
    row.names = FALSE
  )

  object$pseudotime <- pseudotime_values[colnames(object)]
  object$finite_pseudotime <- is.finite(object$pseudotime)

  celltype_levels <- unique(
    as.character(colData(cds)[[celltype_column]])
  )

  celltype_colors <- read_palette(
    file = manifest_row$palette_file[[1]],
    observed_levels = celltype_levels
  )

  p_trajectory <- plot_cells(
    cds,
    reduction_method = "UMAP",
    color_cells_by = celltype_column,
    show_trajectory_graph = TRUE,
    label_cell_groups = FALSE,
    label_groups_by_cluster = FALSE,
    label_branch_points = FALSE,
    label_roots = TRUE,
    label_leaves = FALSE,
    label_principal_points = FALSE,
    cell_size = 0.5
  ) +
    scale_colour_manual(
      values = celltype_colors,
      breaks = celltype_levels,
      limits = celltype_levels,
      drop = FALSE
    ) +
    theme(
      panel.border = element_rect(
        fill = NA,
        color = "black",
        linewidth = 0.5
      )
    )

  p_pseudotime <- plot_cells(
    cds,
    reduction_method = "UMAP",
    color_cells_by = "pseudotime",
    show_trajectory_graph = TRUE,
    trajectory_graph_color = "black",
    trajectory_graph_segment_size = 1.0,
    label_cell_groups = FALSE,
    label_groups_by_cluster = FALSE,
    label_branch_points = FALSE,
    label_roots = TRUE,
    label_leaves = TRUE,
    label_principal_points = FALSE,
    cell_size = 0.5
  ) +
    theme(
      panel.border = element_rect(
        fill = NA,
        color = "black",
        linewidth = 0.5
      )
    )

  p_group <- plot_cells(
    cds,
    reduction_method = "UMAP",
    color_cells_by = group_column,
    show_trajectory_graph = FALSE,
    label_cell_groups = FALSE,
    label_groups_by_cluster = FALSE,
    label_branch_points = FALSE,
    label_roots = FALSE,
    label_leaves = FALSE,
    label_principal_points = FALSE,
    cell_size = 0.5
  ) +
    theme(
      panel.border = element_rect(
        fill = NA,
        color = "black",
        linewidth = 0.5
      )
    )

  save_plot(
    p_trajectory,
    figure_dir,
    "trajectory_by_celltype",
    manifest_row$plot_width[[1]],
    manifest_row$plot_height[[1]]
  )

  save_plot(
    p_pseudotime,
    figure_dir,
    "trajectory_by_pseudotime",
    manifest_row$plot_width[[1]],
    manifest_row$plot_height[[1]]
  )

  save_plot(
    p_group,
    figure_dir,
    "UMAP_by_group",
    manifest_row$plot_width[[1]],
    manifest_row$plot_height[[1]]
  )

  cds_output_file <- file.path(
    output_data_root,
    paste0(analysis_id_safe, "_cds.rds")
  )

  seurat_output_file <- file.path(
    output_data_root,
    paste0(analysis_id_safe, "_Seurat_with_pseudotime.rds")
  )

  saveRDS(cds, cds_output_file, compress = "gzip")
  saveRDS(object, seurat_output_file, compress = "gzip")

  parameters <- tibble(
    parameter = c(
      "analysis_id",
      "input_file",
      "count_assay",
      "seurat_reduction",
      "celltype_column",
      "group_column",
      "root_column",
      "root_values",
      "root_principal_node",
      "num_dim",
      "cores",
      "use_partition",
      "close_loop",
      "random_seed",
      "n_cells",
      "n_finite_pseudotime_cells"
    ),
    value = c(
      analysis_id,
      input_file,
      count_assay,
      seurat_reduction,
      celltype_column,
      group_column,
      root_column,
      paste(root_values, collapse = ";"),
      root_pr_node,
      as.character(manifest_row$num_dim[[1]]),
      as.character(manifest_row$cores[[1]]),
      as.character(manifest_row$use_partition[[1]]),
      as.character(manifest_row$close_loop[[1]]),
      as.character(random_seed),
      as.character(ncol(cds)),
      as.character(sum(is.finite(pseudotime_values)))
    )
  )

  write.csv(
    parameters,
    file = file.path(table_dir, "trajectory_parameters.csv"),
    row.names = FALSE
  )

  tibble(
    analysis_id = analysis_id,
    status = "completed",
    input_file = input_file,
    cds_file = cds_output_file,
    seurat_output_file = seurat_output_file,
    n_cells = ncol(cds),
    n_finite_pseudotime_cells = sum(is.finite(pseudotime_values)),
    root_principal_node = root_pr_node
  )
}


if (!file.exists(manifest_file)) {
  stop(
    "Trajectory manifest does not exist: ",
    manifest_file
  )
}

manifest <- read.csv(
  manifest_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

manifest <- validate_manifest(manifest)

results <- vector("list", nrow(manifest))

for (i in seq_len(nrow(manifest))) {
  current_id <- manifest$analysis_id[[i]]

  results[[i]] <- tryCatch(
    run_one_trajectory(manifest[i, , drop = FALSE]),
    error = function(e) {
      message(
        "Trajectory analysis failed: ",
        current_id,
        "\nReason: ",
        conditionMessage(e)
      )

      tibble(
        analysis_id = current_id,
        status = "failed",
        input_file = manifest$input_file[[i]],
        cds_file = NA_character_,
        seurat_output_file = NA_character_,
        n_cells = NA_integer_,
        n_finite_pseudotime_cells = NA_integer_,
        root_principal_node = NA_character_,
        error_message = conditionMessage(e)
      )
    }
  )
}

summary_table <- bind_rows(results)

write.csv(
  summary_table,
  file = file.path(output_table_root, "08a_trajectory_summary.csv"),
  row.names = FALSE
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(output_table_root, "08a_trajectory_sessionInfo.txt")
)

print(summary_table)
