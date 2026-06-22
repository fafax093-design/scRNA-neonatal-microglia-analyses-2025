# ==============================================================================
# 08c. Monocle 3 pseudotime-gene visualization
#
# Description:
# This script generates pseudotime expression curves and a continuous
# pseudotime heatmap for one or more Monocle 3 analyses. It uses the ordered CDS
# from 08a, graph-test and module assignments from 08b, and the Seurat object
# containing barcode-matched pseudotime values.
#
# Configuration:
#   config/monocle3/08c_visualization_manifest.csv
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(monocle3)
  library(dplyr)
  library(tibble)
  library(Matrix)
  library(ComplexHeatmap)
  library(circlize)
  library(viridisLite)
  library(grid)
})

manifest_file <- "config/monocle3/08c_visualization_manifest.csv"

output_table_root <- "TABLE/Monocle3"
output_figure_root <- "FIGURE/Monocle3"

dir.create(output_table_root, recursive = TRUE, showWarnings = FALSE)
dir.create(output_figure_root, recursive = TRUE, showWarnings = FALSE)


sanitize_name <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", as.character(x))
  gsub("_+", "_", x)
}


read_gene_list <- function(file) {
  file <- trimws(as.character(file))

  if (is.na(file) || file == "" || !file.exists(file)) {
    return(character(0))
  }

  table <- read.csv(
    file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (!"gene" %in% colnames(table)) {
    stop("Selected-gene file must contain a column named gene.")
  }

  genes <- unique(trimws(as.character(table$gene)))
  genes[!is.na(genes) & genes != ""]
}


read_palette <- function(file, levels) {
  file <- trimws(as.character(file))

  if (is.na(file) || file == "" || !file.exists(file)) {
    return(
      setNames(
        grDevices::hcl.colors(length(levels), "Dark 3"),
        levels
      )
    )
  }

  table <- read.csv(
    file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (!all(c("label", "color") %in% colnames(table))) {
    stop("Palette file must contain label and color.")
  }

  missing_levels <- setdiff(levels, table$label)

  if (length(missing_levels) > 0) {
    stop(
      "Palette file is missing labels: ",
      paste(missing_levels, collapse = ", ")
    )
  }

  setNames(
    table$color[match(levels, table$label)],
    levels
  )
}


get_expression_matrix <- function(object, assay, layer) {
  if (!assay %in% Assays(object)) {
    stop("Expression assay was not found: ", assay)
  }

  tryCatch(
    GetAssayData(
      object = object,
      assay = assay,
      layer = layer
    ),
    error = function(e) {
      GetAssayData(
        object = object,
        assay = assay,
        slot = layer
      )
    }
  )
}


scale_01 <- function(x) {
  x_min <- min(x, na.rm = TRUE)
  x_max <- max(x, na.rm = TRUE)

  if (!is.finite(x_min) || !is.finite(x_max) || x_max == x_min) {
    return(rep(0, length(x)))
  }

  (x - x_min) / (x_max - x_min)
}


make_pseudotime_bins <- function(values, method, bin_width, n_bins) {
  if (method == "fixed_width") {
    breaks <- seq(
      min(values, na.rm = TRUE),
      max(values, na.rm = TRUE) + bin_width,
      by = bin_width
    )

    return(
      cut(
        values,
        breaks = breaks,
        include.lowest = TRUE,
        ordered_result = TRUE
      )
    )
  }

  if (method == "equal_cell") {
    probabilities <- seq(0, 1, length.out = n_bins + 1)
    breaks <- unique(
      stats::quantile(
        values,
        probs = probabilities,
        na.rm = TRUE,
        names = FALSE
      )
    )

    if (length(breaks) < 3) {
      stop("Too few unique pseudotime values to construct equal-cell bins.")
    }

    return(
      cut(
        values,
        breaks = breaks,
        include.lowest = TRUE,
        ordered_result = TRUE
      )
    )
  }

  stop("Unsupported bin_method: ", method)
}


save_complex_heatmap <- function(
  heatmap_object,
  output_file,
  width,
  height,
  device = c("pdf", "tiff")
) {
  device <- match.arg(device)

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

  draw(
    heatmap_object,
    merge_legends = TRUE,
    heatmap_legend_side = "right"
  )

  dev.off()
}


run_one_visualization <- function(row) {
  analysis_id <- row$analysis_id[[1]]
  analysis_id_safe <- sanitize_name(analysis_id)

  table_dir <- file.path(output_table_root, analysis_id_safe)
  figure_dir <- file.path(output_figure_root, analysis_id_safe)

  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

  cds <- readRDS(row$cds_file[[1]])
  object <- readRDS(row$seurat_file[[1]])

  if (!inherits(cds, "cell_data_set")) {
    stop("CDS file does not contain a Monocle 3 cell_data_set.")
  }

  if (!inherits(object, "Seurat")) {
    stop("Seurat file does not contain a Seurat object.")
  }

  if (!"pseudotime" %in% colnames(object@meta.data)) {
    stop("The Seurat object does not contain a pseudotime column.")
  }

  graph_test_result <- read.csv(
    row$graph_test_file[[1]],
    row.names = 1,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  module_table <- read.csv(
    row$module_file[[1]],
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (!all(c("id", "module") %in% colnames(module_table))) {
    stop("Module file must contain columns: id and module.")
  }

  selected_genes <- read_gene_list(
    row$selected_gene_file[[1]]
  )

  top_genes <- graph_test_result |>
    rownames_to_column("gene_id") |>
    filter(!is.na(morans_I)) |>
    arrange(desc(morans_I)) |>
    slice_head(n = row$top_n_moran[[1]])

  top_gene_names <- rowData(cds)[top_genes$gene_id, "gene_short_name"]
  top_gene_names <- unique(
    as.character(top_gene_names[!is.na(top_gene_names)])
  )

  celltype_column <- row$celltype_column[[1]]

  if (!celltype_column %in% colnames(colData(cds))) {
    stop("Cell-type column was not found in CDS: ", celltype_column)
  }

  celltype_levels <- unique(
    as.character(colData(cds)[[celltype_column]])
  )

  celltype_colors <- read_palette(
    row$palette_file[[1]],
    celltype_levels
  )

  if (length(top_gene_names) > 0) {
    top_gene_ids <- rownames(cds)[
      rowData(cds)$gene_short_name %in% top_gene_names
    ]

    p_top <- plot_genes_in_pseudotime(
      cds[top_gene_ids, ],
      color_cells_by = celltype_column,
      min_expr = row$min_expr[[1]],
      ncol = row$curve_ncol[[1]],
      cell_size = row$curve_cell_size[[1]]
    ) +
      scale_color_manual(
        values = celltype_colors,
        limits = celltype_levels,
        breaks = celltype_levels,
        drop = FALSE
      )

    ggsave(
      file.path(figure_dir, "pseudotime_curves_top_Moran_genes.pdf"),
      p_top,
      width = row$curve_width[[1]],
      height = row$curve_height[[1]],
      units = "in",
      useDingbats = FALSE
    )

    ggsave(
      file.path(figure_dir, "pseudotime_curves_top_Moran_genes.tiff"),
      p_top,
      width = row$curve_width[[1]],
      height = row$curve_height[[1]],
      units = "in",
      dpi = 600,
      compression = "lzw"
    )
  }

  valid_selected_genes <- intersect(
    selected_genes,
    as.character(rowData(cds)$gene_short_name)
  )

  if (length(valid_selected_genes) > 0) {
    selected_gene_ids <- rownames(cds)[
      rowData(cds)$gene_short_name %in% valid_selected_genes
    ]

    p_selected <- plot_genes_in_pseudotime(
      cds[selected_gene_ids, ],
      color_cells_by = celltype_column,
      min_expr = row$min_expr[[1]],
      ncol = row$curve_ncol[[1]],
      cell_size = row$curve_cell_size[[1]]
    ) +
      scale_color_manual(
        values = celltype_colors,
        limits = celltype_levels,
        breaks = celltype_levels,
        drop = FALSE
      )

    ggsave(
      file.path(figure_dir, "pseudotime_curves_selected_genes.pdf"),
      p_selected,
      width = row$curve_width[[1]],
      height = row$curve_height[[1]],
      units = "in",
      useDingbats = FALSE
    )

    ggsave(
      file.path(figure_dir, "pseudotime_curves_selected_genes.tiff"),
      p_selected,
      width = row$curve_width[[1]],
      height = row$curve_height[[1]],
      units = "in",
      dpi = 600,
      compression = "lzw"
    )
  }

  expression_matrix <- get_expression_matrix(
    object,
    row$expression_assay[[1]],
    row$expression_layer[[1]]
  )

  module_genes <- unique(as.character(module_table$id))

  gene_id_to_symbol <- setNames(
    as.character(rowData(cds)$gene_short_name),
    rownames(cds)
  )

  heatmap_gene_symbols <- gene_id_to_symbol[module_genes]
  names(heatmap_gene_symbols) <- module_genes

  valid_module_rows <- !is.na(heatmap_gene_symbols) &
    heatmap_gene_symbols %in% rownames(expression_matrix)

  module_table <- module_table[valid_module_rows, , drop = FALSE]
  module_table$gene <- unname(
    heatmap_gene_symbols[module_table$id]
  )

  finite_cells <- colnames(object)[
    is.finite(object$pseudotime)
  ]

  finite_cells <- intersect(
    finite_cells,
    colnames(expression_matrix)
  )

  if (length(finite_cells) == 0) {
    stop("No cells with finite pseudotime were available.")
  }

  pseudotime_values <- object$pseudotime[
    match(finite_cells, colnames(object))
  ]

  bins <- make_pseudotime_bins(
    values = pseudotime_values,
    method = row$bin_method[[1]],
    bin_width = row$bin_width[[1]],
    n_bins = row$n_bins[[1]]
  )

  valid_bins <- !is.na(bins)
  finite_cells <- finite_cells[valid_bins]
  bins <- droplevels(bins[valid_bins])

  heatmap_expression <- expression_matrix[
    module_table$gene,
    finite_cells,
    drop = FALSE
  ]

  average_expression <- sapply(
    levels(bins),
    function(bin_label) {
      cells_in_bin <- finite_cells[bins == bin_label]
      Matrix::rowMeans(
        heatmap_expression[
          ,
          cells_in_bin,
          drop = FALSE
        ]
      )
    }
  )

  if (is.null(dim(average_expression))) {
    average_expression <- matrix(
      average_expression,
      ncol = 1,
      dimnames = list(
        rownames(heatmap_expression),
        levels(bins)[1]
      )
    )
  }

  rownames(average_expression) <- module_table$gene

  scaled_expression <- t(
    apply(
      average_expression,
      1,
      scale_01
    )
  )

  peak_bin <- max.col(
    scaled_expression,
    ties.method = "first"
  )

  order_index <- order(
    as.character(module_table$module),
    peak_bin,
    module_table$gene
  )

  scaled_expression <- scaled_expression[
    order_index,
    ,
    drop = FALSE
  ]

  module_table <- module_table[
    order_index,
    ,
    drop = FALSE
  ]

  module_factor <- factor(
    as.character(module_table$module),
    levels = unique(as.character(module_table$module))
  )

  write.csv(
    average_expression,
    file = file.path(table_dir, "pseudotime_bin_mean_expression.csv"),
    row.names = TRUE
  )

  write.csv(
    scaled_expression,
    file = file.path(table_dir, "pseudotime_bin_scaled_expression.csv"),
    row.names = TRUE
  )

  write.csv(
    module_table,
    file = file.path(table_dir, "heatmap_gene_order_and_modules.csv"),
    row.names = FALSE
  )

  marked_genes <- intersect(
    selected_genes,
    rownames(scaled_expression)
  )

  heatmap_object <- Heatmap(
    scaled_expression,
    name = "%Max",
    col = circlize::colorRamp2(
      c(0, 0.5, 1),
      viridisLite::viridis(3)
    ),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    row_split = module_factor,
    use_raster = TRUE,
    border = TRUE
  )

  if (length(marked_genes) > 0) {
    marked_positions <- which(
      rownames(scaled_expression) %in% marked_genes
    )

    heatmap_object <- heatmap_object +
      rowAnnotation(
        genes = anno_mark(
          at = marked_positions,
          labels = rownames(scaled_expression)[marked_positions],
          labels_gp = gpar(fontsize = 8)
        )
      )
  }

  save_complex_heatmap(
    heatmap_object,
    file.path(figure_dir, "pseudotime_gene_heatmap.pdf"),
    row$heatmap_width[[1]],
    row$heatmap_height[[1]],
    "pdf"
  )

  save_complex_heatmap(
    heatmap_object,
    file.path(figure_dir, "pseudotime_gene_heatmap.tiff"),
    row$heatmap_width[[1]],
    row$heatmap_height[[1]],
    "tiff"
  )

  tibble(
    analysis_id = analysis_id,
    status = "completed",
    n_heatmap_genes = nrow(scaled_expression),
    n_pseudotime_bins = ncol(scaled_expression),
    n_selected_curve_genes = length(valid_selected_genes),
    n_top_moran_genes = length(top_gene_names)
  )
}


if (!file.exists(manifest_file)) {
  stop("Visualization manifest does not exist: ", manifest_file)
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
        seurat_file,
        graph_test_file,
        module_file,
        selected_gene_file,
        celltype_column,
        palette_file,
        expression_assay,
        expression_layer,
        bin_method
      ),
      ~ trimws(as.character(.x))
    ),
    top_n_moran = as.integer(top_n_moran),
    min_expr = as.numeric(min_expr),
    curve_ncol = as.integer(curve_ncol),
    curve_cell_size = as.numeric(curve_cell_size),
    curve_width = as.numeric(curve_width),
    curve_height = as.numeric(curve_height),
    bin_width = as.numeric(bin_width),
    n_bins = as.integer(n_bins),
    heatmap_width = as.numeric(heatmap_width),
    heatmap_height = as.numeric(heatmap_height)
  )

required_columns <- c(
  "analysis_id",
  "cds_file",
  "seurat_file",
  "graph_test_file",
  "module_file",
  "selected_gene_file",
  "celltype_column",
  "palette_file",
  "expression_assay",
  "expression_layer",
  "top_n_moran",
  "min_expr",
  "curve_ncol",
  "curve_cell_size",
  "curve_width",
  "curve_height",
  "bin_method",
  "bin_width",
  "n_bins",
  "heatmap_width",
  "heatmap_height"
)

missing_columns <- setdiff(required_columns, colnames(manifest))

if (length(missing_columns) > 0) {
  stop(
    "Visualization manifest is missing columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

results <- vector("list", nrow(manifest))

for (i in seq_len(nrow(manifest))) {
  current_id <- manifest$analysis_id[[i]]

  results[[i]] <- tryCatch(
    run_one_visualization(manifest[i, , drop = FALSE]),
    error = function(e) {
      message(
        "Pseudotime visualization failed: ",
        current_id,
        "\nReason: ",
        conditionMessage(e)
      )

      tibble(
        analysis_id = current_id,
        status = "failed",
        n_heatmap_genes = NA_integer_,
        n_pseudotime_bins = NA_integer_,
        n_selected_curve_genes = NA_integer_,
        n_top_moran_genes = NA_integer_,
        error_message = conditionMessage(e)
      )
    }
  )
}

summary_table <- bind_rows(results)

write.csv(
  summary_table,
  file = file.path(output_table_root, "08c_visualization_summary.csv"),
  row.names = FALSE
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(output_table_root, "08c_visualization_sessionInfo.txt")
)

print(summary_table)
