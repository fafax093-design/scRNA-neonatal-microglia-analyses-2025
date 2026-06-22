# ==============================================================================
# 05. Visualization of clustering, cell states, and marker-gene expression
#     across multiple Seurat objects
#
# Description:
# This script generates manuscript-oriented visualizations for one or more
# annotated Seurat objects defined in an external manifest. Each object may use
# different metadata columns, reductions, assays, clustering prefixes, gene
# lists, cell-state orders, and color palettes.
#
# Available outputs include:
#   1. cluster trees across clustering resolutions;
#   2. UMAPs colored by annotated cell state;
#   3. UMAPs split by experimental group;
#   4. feature-expression UMAPs;
#   5. cell-state composition plots;
#   6. gene-expression boxplots;
#   7. per-cell expression heatmaps;
#   8. marker-gene DotPlots.
#
# Author:
# Jinjin Zhu
#
# Configuration file:
#   config/05_visualization_manifest.csv
#
# Main output directories:
#   FIGURE/STEP5_visualization/<analysis_id>/
#   TABLE/STEP5_visualization/<analysis_id>/
# ==============================================================================


# 1. Load packages -------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(clustree)
  library(SCP)
  library(dplyr)
  library(tibble)
  library(readxl)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})


# 2. Define project paths -------------------------------------------------------

visualization_manifest_file <- "config/05_visualization_manifest.csv"

output_figure_root <- "FIGURE/STEP5_visualization"
output_table_root <- "TABLE/STEP5_visualization"

dir.create(
  output_figure_root,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  output_table_root,
  recursive = TRUE,
  showWarnings = FALSE
)


# 3. Define global plotting parameters -----------------------------------------

random_seed <- 1234

feature_plot_width <- 2
feature_plot_height <- 2

umap_width <- 6
umap_height <- 4

split_umap_width <- 15
split_umap_height <- 12

cell_ratio_width <- 4
cell_ratio_height <- 3

boxplot_width <- 8
boxplot_height <- 4

heatmap_width <- 6.5
heatmap_height <- 6

dotplot_width <- 6
dotplot_height <- 2

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
  required = FALSE
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
        "A required external table path is empty."
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
        "Required external table does not exist: ",
        file
      )
    }

    warning(
      "Optional external table does not exist and was skipped: ",
      file
    )

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


read_gene_list <- function(
  file,
  gene_column = "gene"
) {

  table <- read_external_table(
    file = file,
    required = FALSE
  )

  if (is.null(
    table
  )) {

    return(
      character(0)
    )
  }

  if (!gene_column %in% colnames(
    table
  )) {

    stop(
      "Gene-list file '",
      file,
      "' must contain a column named '",
      gene_column,
      "'."
    )
  }

  genes <- unique(
    trimws(
      as.character(
        table[
          [gene_column]
        ]
      )
    )
  )

  genes[
    !is.na(
      genes
    ) &
      genes != ""
  ]
}


read_palette <- function(
  file,
  observed_celltypes
) {

  palette_table <- read_external_table(
    file = file,
    required = FALSE
  )

  if (is.null(
    palette_table
  )) {

    observed_celltypes <- unique(
      as.character(
        observed_celltypes
      )
    )

    generated_colors <- grDevices::hcl.colors(
      n = length(
        observed_celltypes
      ),
      palette = "Dark 3"
    )

    return(
      list(
        levels = observed_celltypes,
        colors = setNames(
          generated_colors,
          observed_celltypes
        ),
        table = tibble(
          celltype = observed_celltypes,
          color = generated_colors
        )
      )
    )
  }

  required_columns <- c(
    "celltype",
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
      "Palette file '",
      file,
      "' is missing columns: ",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  }

  palette_table <- palette_table %>%
    transmute(
      celltype = trimws(
        as.character(
          celltype
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
        celltype
      ),
      celltype != "",
      !is.na(
        color
      ),
      color != ""
    )

  if (anyDuplicated(
    palette_table$celltype
  )) {

    stop(
      "Duplicated cell-state names were detected in palette file: ",
      file
    )
  }

  missing_celltypes <- setdiff(
    unique(
      as.character(
        observed_celltypes
      )
    ),
    palette_table$celltype
  )

  if (length(
    missing_celltypes
  ) > 0) {

    stop(
      "Palette file '",
      file,
      "' does not contain colors for: ",
      paste(
        missing_celltypes,
        collapse = ", "
      )
    )
  }

  palette_table <- palette_table %>%
    filter(
      celltype %in% unique(
        as.character(
          observed_celltypes
        )
      )
    )

  list(
    levels = palette_table$celltype,
    colors = setNames(
      palette_table$color,
      palette_table$celltype
    ),
    table = palette_table
  )
}


save_ggplot <- function(
  plot_object,
  output_dir,
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
    units = "in"
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


get_assay_matrix <- function(
  object,
  assay,
  data_slot
) {

  if (!assay %in% Assays(
    object
  )) {

    stop(
      "Assay '",
      assay,
      "' was not found."
    )
  }

  matrix <- tryCatch(
    GetAssayData(
      object = object,
      assay = assay,
      layer = data_slot
    ),
    error = function(e) {
      GetAssayData(
        object = object,
        assay = assay,
        slot = data_slot
      )
    }
  )

  matrix
}


min_max_scale <- function(x) {

  x_min <- min(
    x,
    na.rm = TRUE
  )

  x_max <- max(
    x,
    na.rm = TRUE
  )

  if (!is.finite(
    x_min
  ) ||
      !is.finite(
        x_max
      )) {

    return(
      rep(
        NA_real_,
        length(
          x
        )
      )
    )
  }

  if (x_max == x_min) {

    return(
      rep(
        0,
        length(
          x
        )
      )
    )
  }

  (
    x -
      x_min
  ) /
    (
      x_max -
        x_min
    )
}


validate_manifest <- function(manifest) {

  required_columns <- c(
    "analysis_id",
    "input_file",
    "celltype_column",
    "group_column",
    "reduction",
    "expression_assay",
    "expression_slot",
    "cluster_prefix",
    "palette_file",
    "feature_gene_file",
    "boxplot_gene_file",
    "heatmap_gene_file",
    "dotplot_gene_file",
    "run_clustree",
    "run_umap",
    "run_feature_plots",
    "run_cell_proportions",
    "run_boxplots",
    "run_heatmap",
    "run_dotplot"
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
      "The visualization manifest is missing columns: ",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  }

  character_columns <- c(
    "analysis_id",
    "input_file",
    "celltype_column",
    "group_column",
    "reduction",
    "expression_assay",
    "expression_slot",
    "cluster_prefix",
    "palette_file",
    "feature_gene_file",
    "boxplot_gene_file",
    "heatmap_gene_file",
    "dotplot_gene_file"
  )

  flag_columns <- c(
    "run_clustree",
    "run_umap",
    "run_feature_plots",
    "run_cell_proportions",
    "run_boxplots",
    "run_heatmap",
    "run_dotplot"
  )

  manifest <- manifest %>%
    mutate(
      across(
        all_of(
          character_columns
        ),
        ~ trimws(
          as.character(
            .x
          )
        )
      ),
      across(
        all_of(
          flag_columns
        ),
        as_flag
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

    duplicated_ids <- unique(
      manifest$analysis_id[
        duplicated(
          manifest$analysis_id
        )
      ]
    )

    stop(
      "Duplicated analysis_id values were detected: ",
      paste(
        duplicated_ids,
        collapse = ", "
      )
    )
  }

  manifest
}


# 5. Define the per-object visualization function ------------------------------

visualize_one_object <- function(
  manifest_row
) {

  analysis_id <- manifest_row$analysis_id[
    [1]
  ]

  analysis_id_safe <- sanitize_name(
    analysis_id
  )

  input_file <- manifest_row$input_file[
    [1]
  ]

  celltype_column <- manifest_row$celltype_column[
    [1]
  ]

  group_column <- manifest_row$group_column[
    [1]
  ]

  reduction_name <- manifest_row$reduction[
    [1]
  ]

  expression_assay <- manifest_row$expression_assay[
    [1]
  ]

  expression_slot <- manifest_row$expression_slot[
    [1]
  ]

  cluster_prefix <- manifest_row$cluster_prefix[
    [1]
  ]

  figure_dir <- file.path(
    output_figure_root,
    analysis_id_safe
  )

  table_dir <- file.path(
    output_table_root,
    analysis_id_safe
  )

  feature_figure_dir <- file.path(
    figure_dir,
    "feature_plots"
  )

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
    feature_figure_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  message(
    "\nGenerating visualizations for: ",
    analysis_id
  )

  if (!file.exists(
    input_file
  )) {

    stop(
      "Input file does not exist: ",
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
    celltype_column,
    group_column
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
      "Analysis '",
      analysis_id,
      "' is missing metadata columns: ",
      paste(
        missing_metadata,
        collapse = ", "
      )
    )
  }

  if (!reduction_name %in% Reductions(
    object
  )) {

    stop(
      "Reduction '",
      reduction_name,
      "' was not found for analysis '",
      analysis_id,
      "'."
    )
  }

  if (!expression_assay %in% Assays(
    object
  )) {

    stop(
      "Expression assay '",
      expression_assay,
      "' was not found for analysis '",
      analysis_id,
      "'."
    )
  }

  observed_celltypes <- unique(
    as.character(
      object@meta.data[
        [celltype_column]
      ]
    )
  )

  observed_celltypes <- observed_celltypes[
    !is.na(
      observed_celltypes
    ) &
      observed_celltypes != ""
  ]

  if (length(
    observed_celltypes
  ) == 0) {

    stop(
      "No valid cell-state labels were found in column: ",
      celltype_column
    )
  }

  palette_info <- read_palette(
    file = manifest_row$palette_file[
      [1]
    ],
    observed_celltypes = observed_celltypes
  )

  celltype_levels <- palette_info$levels
  celltype_colors <- palette_info$colors

  object[
    [celltype_column]
  ] <- factor(
    as.character(
      object@meta.data[
        [celltype_column]
      ]
    ),
    levels = celltype_levels
  )

  write.csv(
    palette_info$table,
    file = file.path(
      table_dir,
      "celltype_palette_used.csv"
    ),
    row.names = FALSE
  )

  object_summary <- tibble(
    analysis_id = analysis_id,
    input_file = input_file,
    n_cells = ncol(
      object
    ),
    n_features = nrow(
      object
    ),
    n_celltypes = length(
      celltype_levels
    ),
    n_groups = dplyr::n_distinct(
      object@meta.data[
        [group_column]
      ]
    ),
    reduction = reduction_name,
    expression_assay = expression_assay,
    expression_slot = expression_slot
  )

  write.csv(
    object_summary,
    file = file.path(
      table_dir,
      "visualization_input_summary.csv"
    ),
    row.names = FALSE
  )


  # 5.1 Cluster-tree visualization --------------------------------------------

  if (manifest_row$run_clustree[
    [1]
  ]) {

    if (cluster_prefix == "") {

      warning(
        "Cluster-tree plotting was requested for '",
        analysis_id,
        "', but cluster_prefix is empty."
      )

    } else {

      cluster_pattern <- paste0(
        "^",
        gsub(
          pattern = "([.])",
          replacement = "\\\\\\1",
          x = cluster_prefix
        )
      )

      cluster_columns <- grep(
        pattern = cluster_pattern,
        x = colnames(
          object@meta.data
        ),
        value = TRUE
      )

      if (length(
        cluster_columns
      ) < 2) {

        warning(
          "Fewer than two clustering columns matched prefix '",
          cluster_prefix,
          "' for analysis '",
          analysis_id,
          "'. The cluster tree was skipped."
        )

      } else {

        p_clustree <- clustree::clustree(
          object@meta.data,
          prefix = cluster_prefix
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

        save_ggplot(
          plot_object = p_clustree,
          output_dir = figure_dir,
          file_stem = "ClusterTree",
          width = 8,
          height = 8
        )
      }
    }
  }


  # 5.2 UMAPs by cell state and group -----------------------------------------

  if (manifest_row$run_umap[
    [1]
  ]) {

    p_umap_celltype <- SCP::CellDimPlot(
      srt = object,
      group.by = celltype_column,
      reduction = reduction_name,
      theme_use = "theme_blank",
      label = TRUE,
      label_insitu = TRUE
    ) +
      scale_colour_manual(
        values = celltype_colors,
        limits = celltype_levels,
        breaks = celltype_levels,
        drop = FALSE
      )

    save_ggplot(
      plot_object = p_umap_celltype,
      output_dir = figure_dir,
      file_stem = "UMAP_celltype",
      width = umap_width,
      height = umap_height
    )

    p_umap_group <- SCP::CellDimPlot(
      srt = object,
      group.by = celltype_column,
      reduction = reduction_name,
      theme_use = "theme_blank",
      label = FALSE,
      label_insitu = FALSE,
      show_stat = FALSE,
      split.by = group_column
    ) +
      scale_colour_manual(
        values = celltype_colors,
        limits = celltype_levels,
        breaks = celltype_levels,
        drop = FALSE
      )

    save_ggplot(
      plot_object = p_umap_group,
      output_dir = figure_dir,
      file_stem = "UMAP_split_by_group",
      width = split_umap_width,
      height = split_umap_height
    )
  }


  # 5.3 Feature-expression UMAPs ----------------------------------------------

  if (manifest_row$run_feature_plots[
    [1]
  ]) {

    feature_genes <- read_gene_list(
      file = manifest_row$feature_gene_file[
        [1]
      ]
    )

    available_feature_genes <- intersect(
      feature_genes,
      rownames(
        object
      )
    )

    missing_feature_genes <- setdiff(
      feature_genes,
      available_feature_genes
    )

    if (length(
      missing_feature_genes
    ) > 0) {

      warning(
        "Feature genes absent from analysis '",
        analysis_id,
        "': ",
        paste(
          missing_feature_genes,
          collapse = ", "
        )
      )
    }

    for (gene_name in available_feature_genes) {

      p_feature <- SCP::FeatureDimPlot(
        srt = object,
        features = gene_name,
        reduction = reduction_name,
        cells.highlight = TRUE,
        theme_use = "theme_blank",
        show_stat = FALSE,
        legend.position = "none"
      )

      save_ggplot(
        plot_object = p_feature,
        output_dir = feature_figure_dir,
        file_stem = sanitize_name(
          gene_name
        ),
        width = feature_plot_width,
        height = feature_plot_height
      )
    }
  }


  # 5.4 Cell-state composition plot -------------------------------------------

  if (manifest_row$run_cell_proportions[
    [1]
  ]) {

    p_cell_ratio <- SCP::CellStatPlot(
      srt = object,
      stat.by = celltype_column,
      group.by = group_column,
      plot_type = "trend"
    )

    save_ggplot(
      plot_object = p_cell_ratio,
      output_dir = figure_dir,
      file_stem = "Cell_state_proportions",
      width = cell_ratio_width,
      height = cell_ratio_height
    )

    cell_composition <- object@meta.data %>%
      transmute(
        group = as.character(
          .data[
            [group_column]
          ]
        ),
        celltype = as.character(
          .data[
            [celltype_column]
          ]
        )
      ) %>%
      count(
        group,
        celltype,
        name = "n_cells"
      ) %>%
      group_by(
        group
      ) %>%
      mutate(
        total_cells = sum(
          n_cells
        ),
        percentage = 100 *
          n_cells /
          total_cells
      ) %>%
      ungroup()

    write.csv(
      cell_composition,
      file = file.path(
        table_dir,
        "Cell_state_composition.csv"
      ),
      row.names = FALSE
    )
  }


  # 5.5 Gene-expression boxplots ----------------------------------------------

  if (manifest_row$run_boxplots[
    [1]
  ]) {

    boxplot_genes <- read_gene_list(
      file = manifest_row$boxplot_gene_file[
        [1]
      ]
    )

    available_boxplot_genes <- intersect(
      boxplot_genes,
      rownames(
        object
      )
    )

    missing_boxplot_genes <- setdiff(
      boxplot_genes,
      available_boxplot_genes
    )

    if (length(
      missing_boxplot_genes
    ) > 0) {

      warning(
        "Boxplot genes absent from analysis '",
        analysis_id,
        "': ",
        paste(
          missing_boxplot_genes,
          collapse = ", "
        )
      )
    }

    if (length(
      available_boxplot_genes
    ) > 0) {

      p_gene_boxplot <- SCP::FeatureStatPlot(
        srt = object,
        stat.by = available_boxplot_genes,
        fill.by = group_column,
        plot_type = "box",
        group.by = celltype_column,
        bg.by = celltype_column,
        stack = TRUE,
        flip = FALSE
      )

      save_ggplot(
        plot_object = p_gene_boxplot,
        output_dir = figure_dir,
        file_stem = "Gene_expression_boxplots",
        width = boxplot_width,
        height = boxplot_height
      )
    }
  }


  # 5.6 Per-cell expression heatmap -------------------------------------------

  if (manifest_row$run_heatmap[
    [1]
  ]) {

    heatmap_gene_table <- read_external_table(
      file = manifest_row$heatmap_gene_file[
        [1]
      ],
      required = TRUE
    )

    if ("cluster" %in% colnames(
      heatmap_gene_table
    ) &&
        !"module" %in% colnames(
          heatmap_gene_table
        )) {

      heatmap_gene_table <- heatmap_gene_table %>%
        rename(
          module = cluster
        )
    }

    required_heatmap_columns <- c(
      "module",
      "gene"
    )

    missing_heatmap_columns <- setdiff(
      required_heatmap_columns,
      colnames(
        heatmap_gene_table
      )
    )

    if (length(
      missing_heatmap_columns
    ) > 0) {

      stop(
        "Heatmap gene file must contain columns: module and gene."
      )
    }

    heatmap_gene_table <- heatmap_gene_table %>%
      transmute(
        module = trimws(
          as.character(
            module
          )
        ),
        gene = trimws(
          as.character(
            gene
          )
        )
      ) %>%
      filter(
        !is.na(
          gene
        ),
        gene != ""
      ) %>%
      distinct(
        gene,
        .keep_all = TRUE
      )

    expression_matrix <- get_assay_matrix(
      object = object,
      assay = expression_assay,
      data_slot = expression_slot
    )

    valid_heatmap_gene_table <- heatmap_gene_table %>%
      filter(
        gene %in% rownames(
          expression_matrix
        )
      )

    missing_heatmap_genes <- setdiff(
      heatmap_gene_table$gene,
      valid_heatmap_gene_table$gene
    )

    if (length(
      missing_heatmap_genes
    ) > 0) {

      warning(
        "Heatmap genes absent from analysis '",
        analysis_id,
        "': ",
        paste(
          missing_heatmap_genes,
          collapse = ", "
        )
      )
    }

    if (nrow(
      valid_heatmap_gene_table
    ) == 0) {

      stop(
        "None of the requested heatmap genes were found for analysis: ",
        analysis_id
      )
    }

    heatmap_data <- as.matrix(
      expression_matrix[
        valid_heatmap_gene_table$gene,
        ,
        drop = FALSE
      ]
    )

    rownames(
      heatmap_data
    ) <- valid_heatmap_gene_table$gene

    heatmap_celltypes <- factor(
      as.character(
        object@meta.data[
          colnames(
            heatmap_data
          ),
          celltype_column
        ]
      ),
      levels = celltype_levels
    )

    column_order <- order(
      heatmap_celltypes
    )

    heatmap_data <- heatmap_data[
      ,
      column_order,
      drop = FALSE
    ]

    heatmap_celltypes <- heatmap_celltypes[
      column_order
    ]

    row_modules <- factor(
      valid_heatmap_gene_table$module,
      levels = unique(
        valid_heatmap_gene_table$module
      )
    )

    scaled_heatmap_data <- t(
      apply(
        heatmap_data,
        1,
        min_max_scale
      )
    )

    rownames(
      scaled_heatmap_data
    ) <- rownames(
      heatmap_data
    )

    colnames(
      scaled_heatmap_data
    ) <- colnames(
      heatmap_data
    )

    row_module_levels <- levels(
      row_modules
    )

    row_module_colors <- setNames(
      grDevices::hcl.colors(
        n = length(
          row_module_levels
        ),
        palette = "Dark 3"
      ),
      row_module_levels
    )

    top_annotation <- ComplexHeatmap::HeatmapAnnotation(
      celltype = heatmap_celltypes,
      col = list(
        celltype = celltype_colors
      ),
      show_annotation_name = FALSE
    )

    left_annotation <- ComplexHeatmap::rowAnnotation(
      module = row_modules,
      col = list(
        module = row_module_colors
      ),
      show_annotation_name = FALSE
    )

    heatmap_object <- ComplexHeatmap::Heatmap(
      scaled_heatmap_data,
      name = "Expression",
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      show_column_names = FALSE,
      show_row_names = FALSE,
      column_split = heatmap_celltypes,
      row_split = row_modules,
      top_annotation = top_annotation,
      left_annotation = left_annotation,
      col = circlize::colorRamp2(
        c(
          0,
          0.5,
          1
        ),
        c(
          "#4978B3",
          "white",
          "#FF3333"
        )
      ),
      heatmap_legend_param = list(
        at = seq(
          0,
          1,
          0.2
        ),
        labels = seq(
          0,
          1,
          0.2
        ),
        title = "Expression",
        title_position = "leftcenter-rot"
      ),
      border = TRUE,
      use_raster = TRUE,
      raster_quality = 2,
      column_gap = unit(
        1,
        "mm"
      ),
      row_gap = unit(
        1,
        "mm"
      )
    )

    pdf(
      file = file.path(
        figure_dir,
        "GENE_EXPRESSION.pdf"
      ),
      width = heatmap_width,
      height = heatmap_height,
      useDingbats = FALSE
    )

    ComplexHeatmap::draw(
      heatmap_object,
      merge_legends = TRUE
    )

    dev.off()

    tiff(
      filename = file.path(
        figure_dir,
        "GENE_EXPRESSION.tiff"
      ),
      width = heatmap_width,
      height = heatmap_height,
      units = "in",
      res = 600,
      compression = "lzw"
    )

    ComplexHeatmap::draw(
      heatmap_object,
      merge_legends = TRUE
    )

    dev.off()

    write.csv(
      scaled_heatmap_data,
      file = file.path(
        table_dir,
        "Heatmap_scaled_expression_matrix.csv"
      ),
      row.names = TRUE
    )

    write.csv(
      valid_heatmap_gene_table,
      file = file.path(
        table_dir,
        "Heatmap_gene_annotations.csv"
      ),
      row.names = FALSE
    )
  }


  # 5.7 Marker-gene DotPlot ----------------------------------------------------

  if (manifest_row$run_dotplot[
    [1]
  ]) {

    marker_genes <- read_gene_list(
      file = manifest_row$dotplot_gene_file[
        [1]
      ]
    )

    available_marker_genes <- intersect(
      marker_genes,
      rownames(
        object
      )
    )

    missing_marker_genes <- setdiff(
      marker_genes,
      available_marker_genes
    )

    if (length(
      missing_marker_genes
    ) > 0) {

      warning(
        "DotPlot genes absent from analysis '",
        analysis_id,
        "': ",
        paste(
          missing_marker_genes,
          collapse = ", "
        )
      )
    }

    if (length(
      available_marker_genes
    ) > 0) {

      Idents(
        object
      ) <- celltype_column

      dotplot_base <- DotPlot(
        object = object,
        features = available_marker_genes,
        assay = expression_assay
      )

      dot_data <- dotplot_base$data

      write.csv(
        dot_data,
        file = file.path(
          table_dir,
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

      save_ggplot(
        plot_object = p_dotplot,
        output_dir = figure_dir,
        file_stem = "DOTPLOT_NOT_SLICED",
        width = dotplot_width,
        height = dotplot_height
      )
    }
  }


  # 5.8 Save analysis information --------------------------------------------

  analysis_parameters <- tibble(
    parameter = c(
      "analysis_id",
      "input_file",
      "celltype_column",
      "group_column",
      "reduction",
      "expression_assay",
      "expression_slot",
      "cluster_prefix",
      "palette_file",
      "feature_gene_file",
      "boxplot_gene_file",
      "heatmap_gene_file",
      "dotplot_gene_file",
      "random_seed"
    ),
    value = c(
      analysis_id,
      input_file,
      celltype_column,
      group_column,
      reduction_name,
      expression_assay,
      expression_slot,
      cluster_prefix,
      manifest_row$palette_file[
        [1]
      ],
      manifest_row$feature_gene_file[
        [1]
      ],
      manifest_row$boxplot_gene_file[
        [1]
      ],
      manifest_row$heatmap_gene_file[
        [1]
      ],
      manifest_row$dotplot_gene_file[
        [1]
      ],
      as.character(
        random_seed
      )
    )
  )

  write.csv(
    analysis_parameters,
    file = file.path(
      table_dir,
      "Visualization_parameters.csv"
    ),
    row.names = FALSE
  )

  tibble(
    analysis_id = analysis_id,
    status = "completed",
    input_file = input_file,
    n_cells = ncol(
      object
    ),
    n_features = nrow(
      object
    ),
    n_celltypes = length(
      celltype_levels
    ),
    n_groups = dplyr::n_distinct(
      object@meta.data[
        [group_column]
      ]
    )
  )
}


# 6. Load and validate the visualization manifest -----------------------------

if (!file.exists(
  visualization_manifest_file
)) {

  stop(
    "Visualization manifest does not exist: ",
    visualization_manifest_file,
    "\nCreate it from the provided template before running this script."
  )
}

visualization_manifest <- read.csv(
  visualization_manifest_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

visualization_manifest <- validate_manifest(
  visualization_manifest
)


# 7. Process all configured Seurat objects -------------------------------------

visualization_results <- vector(
  mode = "list",
  length = nrow(
    visualization_manifest
  )
)

for (i in seq_len(
  nrow(
    visualization_manifest
  )
)) {

  manifest_row <- visualization_manifest[
    i,
    ,
    drop = FALSE
  ]

  analysis_id <- manifest_row$analysis_id[
    [1]
  ]

  visualization_results[
    [i]
  ] <- tryCatch(
    visualize_one_object(
      manifest_row = manifest_row
    ),
    error = function(e) {

      message(
        "Visualization failed: ",
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
        n_cells = NA_integer_,
        n_features = NA_integer_,
        n_celltypes = NA_integer_,
        n_groups = NA_integer_,
        error_message = conditionMessage(
          e
        )
      )
    }
  )
}

visualization_summary <- bind_rows(
  visualization_results
)

write.csv(
  visualization_summary,
  file = file.path(
    output_table_root,
    "05_visualization_summary.csv"
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
    "05_visualization_sessionInfo.txt"
  )
)

message(
  "Configured visualization analyses have finished."
)

print(
  visualization_summary
)
