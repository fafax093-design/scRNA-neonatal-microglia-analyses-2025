# ==============================================================================
# 06. Gene Ontology enrichment analysis and bubble-plot visualization
#     across multiple DEG tables
#
# Description:
# This script performs Gene Ontology enrichment for one or more differential-
# expression tables defined in an external manifest. For each analysis, it:
#
#   1. filters positively regulated genes using configurable thresholds;
#   2. maps mouse gene symbols to Entrez identifiers;
#   3. performs GO enrichment separately for each configured group;
#   4. exports complete enrichment and gene-mapping results;
#   5. generates a manuscript-oriented GO bubble plot;
#   6. optionally adds a left annotation panel for selected GO terms.
#
# Different analyses may use different input files, column names, thresholds,
# ontologies, annotation files, and figure dimensions.
#
# Author:
# Jinjin Zhu
#
# Configuration file:
#   config/06_GO_enrichment_manifest.csv
#
# Main output directories:
#   TABLE/GO/<analysis_id>/
#   FIGURE/GO/<analysis_id>/
# ==============================================================================


# 1. Load packages -------------------------------------------------------------

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(dplyr)
  library(tibble)
  library(readxl)
  library(openxlsx)
  library(ggplot2)
  library(cowplot)
  library(aplot)
})


# 2. Define project paths -------------------------------------------------------

analysis_manifest_file <- "config/06_GO_enrichment_manifest.csv"

output_table_root <- "TABLE/GO"
output_figure_root <- "FIGURE/GO"

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


# 3. Define global settings -----------------------------------------------------

adjust_group_levels <- c(
  "<0.0001",
  "<0.001",
  "<0.01",
  "<0.05",
  "<0.1",
  ">0.1"
)

adjust_group_colors <- c(
  "<0.0001" = "#67000D",
  "<0.001" = "#EF3B2C",
  "<0.01" = "#FB6A4A",
  "<0.05" = "#FC9272",
  "<0.1" = "#FEE0D2",
  ">0.1" = "#FFF5F0"
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
      "Unsupported input format: ",
      file
    )
  }

  table
}


validate_manifest <- function(manifest) {

  required_columns <- c(
    "analysis_id",
    "deg_file",
    "group_column",
    "gene_column",
    "logfc_column",
    "adjusted_p_column",
    "logfc_cutoff",
    "adjusted_p_cutoff",
    "go_ontology",
    "p_adjust_method",
    "enrichment_p_cutoff",
    "enrichment_q_cutoff",
    "min_gene_set_size",
    "max_gene_set_size",
    "term_annotation_file",
    "restrict_plot_to_annotation_terms",
    "plot_width",
    "plot_height"
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
      "The GO analysis manifest is missing columns: ",
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
          deg_file,
          group_column,
          gene_column,
          logfc_column,
          adjusted_p_column,
          go_ontology,
          p_adjust_method,
          term_annotation_file
        ),
        ~ trimws(
          as.character(
            .x
          )
        )
      ),
      logfc_cutoff = as.numeric(
        logfc_cutoff
      ),
      adjusted_p_cutoff = as.numeric(
        adjusted_p_cutoff
      ),
      enrichment_p_cutoff = as.numeric(
        enrichment_p_cutoff
      ),
      enrichment_q_cutoff = as.numeric(
        enrichment_q_cutoff
      ),
      min_gene_set_size = as.integer(
        min_gene_set_size
      ),
      max_gene_set_size = as.integer(
        max_gene_set_size
      ),
      restrict_plot_to_annotation_terms = as_flag(
        restrict_plot_to_annotation_terms
      ),
      plot_width = as.numeric(
        plot_width
      ),
      plot_height = as.numeric(
        plot_height
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

  valid_ontologies <- c(
    "BP",
    "CC",
    "MF",
    "ALL"
  )

  invalid_ontologies <- setdiff(
    unique(
      manifest$go_ontology
    ),
    valid_ontologies
  )

  if (length(
    invalid_ontologies
  ) > 0) {

    stop(
      "Unsupported GO ontology values were detected: ",
      paste(
        invalid_ontologies,
        collapse = ", "
      )
    )
  }

  manifest
}


prepare_term_annotations <- function(
  file,
  available_descriptions,
  restrict_plot
) {

  term_annotations <- read_external_table(
    file = file,
    required = FALSE
  )

  if (is.null(
    term_annotations
  )) {

    return(
      list(
        table = NULL,
        descriptions = available_descriptions
      )
    )
  }

  required_columns <- c(
    "Description",
    "Annotation"
  )

  missing_columns <- setdiff(
    required_columns,
    colnames(
      term_annotations
    )
  )

  if (length(
    missing_columns
  ) > 0) {

    stop(
      "The GO-term annotation file is missing columns: ",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  }

  term_annotations <- term_annotations %>%
    mutate(
      Description = trimws(
        as.character(
          Description
        )
      ),
      Annotation = trimws(
        as.character(
          Annotation
        )
      )
    ) %>%
    filter(
      !is.na(
        Description
      ),
      Description != "",
      !is.na(
        Annotation
      ),
      Annotation != ""
    ) %>%
    distinct(
      Description,
      .keep_all = TRUE
    )

  if ("order" %in% colnames(
    term_annotations
  )) {

    term_annotations$order <- as.numeric(
      term_annotations$order
    )

    term_annotations <- term_annotations %>%
      arrange(
        order
      )
  }

  matched_descriptions <- intersect(
    term_annotations$Description,
    available_descriptions
  )

  if (length(
    matched_descriptions
  ) == 0) {

    stop(
      "None of the terms in the annotation file were found in the ",
      "enrichment result."
    )
  }

  if (restrict_plot) {

    descriptions <- matched_descriptions

  } else {

    descriptions <- available_descriptions
  }

  list(
    table = term_annotations,
    descriptions = descriptions
  )
}


save_plot <- function(
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


# 5. Define the per-analysis GO function ---------------------------------------

run_one_go_analysis <- function(
  manifest_row
) {

  analysis_id <- manifest_row$analysis_id[
    [1]
  ]

  analysis_id_safe <- sanitize_name(
    analysis_id
  )

  output_table_dir <- file.path(
    output_table_root,
    analysis_id_safe
  )

  output_figure_dir <- file.path(
    output_figure_root,
    analysis_id_safe
  )

  dir.create(
    output_table_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  dir.create(
    output_figure_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  message(
    "\nRunning GO analysis: ",
    analysis_id
  )

  deg_file <- manifest_row$deg_file[
    [1]
  ]

  group_column <- manifest_row$group_column[
    [1]
  ]

  gene_column <- manifest_row$gene_column[
    [1]
  ]

  logfc_column <- manifest_row$logfc_column[
    [1]
  ]

  adjusted_p_column <- manifest_row$adjusted_p_column[
    [1]
  ]

  deg_table <- read_external_table(
    file = deg_file,
    required = TRUE
  )

  required_deg_columns <- c(
    group_column,
    gene_column,
    logfc_column,
    adjusted_p_column
  )

  missing_deg_columns <- setdiff(
    required_deg_columns,
    colnames(
      deg_table
    )
  )

  if (length(
    missing_deg_columns
  ) > 0) {

    stop(
      "DEG table '",
      deg_file,
      "' is missing columns: ",
      paste(
        missing_deg_columns,
        collapse = ", "
      )
    )
  }

  markers <- deg_table %>%
    transmute(
      group = trimws(
        as.character(
          .data[
            [group_column]
          ]
        )
      ),
      gene = trimws(
        as.character(
          .data[
            [gene_column]
          ]
        )
      ),
      avg_log2FC = as.numeric(
        .data[
          [logfc_column]
        ]
      ),
      p_val_adj = as.numeric(
        .data[
          [adjusted_p_column]
        ]
      )
    ) %>%
    filter(
      !is.na(
        group
      ),
      group != "",
      !is.na(
        gene
      ),
      gene != "",
      !is.na(
        avg_log2FC
      ),
      !is.na(
        p_val_adj
      ),
      avg_log2FC >
        manifest_row$logfc_cutoff[
          [1]
        ],
      p_val_adj <
        manifest_row$adjusted_p_cutoff[
          [1]
        ]
    ) %>%
    distinct(
      group,
      gene,
      .keep_all = TRUE
    )

  if (nrow(
    markers
  ) == 0) {

    stop(
      "No genes passed the DEG filtering criteria for analysis: ",
      analysis_id
    )
  }

  group_order <- unique(
    markers$group
  )

  write.csv(
    markers,
    file = file.path(
      output_table_dir,
      "GO_input_genes.csv"
    ),
    row.names = FALSE
  )

  selected_gene_counts <- markers %>%
    count(
      group,
      name = "n_selected_genes"
    )

  write.csv(
    selected_gene_counts,
    file = file.path(
      output_table_dir,
      "GO_input_gene_counts_by_group.csv"
    ),
    row.names = FALSE
  )


  # 5.1 Map mouse gene symbols to Entrez identifiers --------------------------

  gene_mapping <- clusterProfiler::bitr(
    unique(
      markers$gene
    ),
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Mm.eg.db
  )

  markers_mapped <- markers %>%
    inner_join(
      gene_mapping,
      by = c(
        "gene" = "SYMBOL"
      )
    ) %>%
    distinct(
      group,
      gene,
      ENTREZID,
      .keep_all = TRUE
    )

  if (nrow(
    markers_mapped
  ) == 0) {

    stop(
      "None of the selected genes could be mapped to Entrez identifiers."
    )
  }

  write.csv(
    markers_mapped,
    file = file.path(
      output_table_dir,
      "GO_input_genes_with_Entrez_IDs.csv"
    ),
    row.names = FALSE
  )

  unmapped_symbols <- setdiff(
    unique(
      markers$gene
    ),
    unique(
      markers_mapped$gene
    )
  )

  write.csv(
    tibble(
      gene = unmapped_symbols
    ),
    file = file.path(
      output_table_dir,
      "GO_unmapped_gene_symbols.csv"
    ),
    row.names = FALSE
  )

  mapping_summary <- tibble(
    total_selected_symbols = length(
      unique(
        markers$gene
      )
    ),
    mapped_symbols = length(
      unique(
        markers_mapped$gene
      )
    ),
    unmapped_symbols = length(
      unmapped_symbols
    ),
    mapping_percentage = 100 *
      mapped_symbols /
      total_selected_symbols
  )

  write.csv(
    mapping_summary,
    file = file.path(
      output_table_dir,
      "GO_gene_mapping_summary.csv"
    ),
    row.names = FALSE
  )


  # 5.2 Perform GO enrichment by group ----------------------------------------

  gene_clusters <- split(
    markers_mapped$ENTREZID,
    markers_mapped$group
  )

  gene_clusters <- lapply(
    gene_clusters,
    unique
  )

  gene_clusters <- gene_clusters[
    group_order[
      group_order %in% names(
        gene_clusters
      )
    ]
  ]

  ego <- compareCluster(
    geneCluster = gene_clusters,
    fun = "enrichGO",
    OrgDb = org.Mm.eg.db,
    keyType = "ENTREZID",
    ont = manifest_row$go_ontology[
      [1]
    ],
    pAdjustMethod = manifest_row$p_adjust_method[
      [1]
    ],
    pvalueCutoff = manifest_row$enrichment_p_cutoff[
      [1]
    ],
    qvalueCutoff = manifest_row$enrichment_q_cutoff[
      [1]
    ],
    minGSSize = manifest_row$min_gene_set_size[
      [1]
    ],
    maxGSSize = manifest_row$max_gene_set_size[
      [1]
    ],
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

  if (nrow(
    ego_df
  ) == 0) {

    stop(
      "GO enrichment returned no terms for analysis: ",
      analysis_id
    )
  }

  write.csv(
    ego_df,
    file = file.path(
      output_table_dir,
      "Enrichment_GO_complete.csv"
    ),
    row.names = FALSE
  )

  saveRDS(
    ego_readable,
    file = file.path(
      output_table_dir,
      "GO_compareCluster_result.rds"
    ),
    compress = "gzip"
  )


  # 5.3 Prepare bubble-plot data ----------------------------------------------

  available_descriptions <- unique(
    as.character(
      ego_df$Description
    )
  )

  annotation_info <- prepare_term_annotations(
    file = manifest_row$term_annotation_file[
      [1]
    ],
    available_descriptions = available_descriptions,
    restrict_plot =
      manifest_row$restrict_plot_to_annotation_terms[
        [1]
      ]
  )

  plot_data <- ego_df %>%
    filter(
      !is.na(
        p.adjust
      ),
      Description %in%
        annotation_info$descriptions
    ) %>%
    mutate(
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
        labels = adjust_group_levels,
        right = FALSE
      )
    )

  if (nrow(
    plot_data
  ) == 0) {

    stop(
      "No GO terms remained for bubble-plot generation."
    )
  }

  plot_data$adjust_group <- factor(
    plot_data$adjust_group,
    levels = adjust_group_levels
  )

  plot_data$Cluster <- factor(
    as.character(
      plot_data$Cluster
    ),
    levels = group_order
  )

  if (!is.null(
    annotation_info$table
  ) &&
      manifest_row$restrict_plot_to_annotation_terms[
        [1]
      ]) {

    description_order <- annotation_info$table$Description[
      annotation_info$table$Description %in%
        unique(
          as.character(
            plot_data$Description
          )
        )
    ]

  } else {

    description_order <- unique(
      as.character(
        plot_data$Description
      )
    )
  }

  plot_data$Description <- factor(
    plot_data$Description,
    levels = rev(
      description_order
    )
  )

  write.csv(
    plot_data,
    file = file.path(
      output_table_dir,
      "GO_bubble_plot_data.csv"
    ),
    row.names = FALSE
  )


  # 5.4 Generate the GO bubble plot -------------------------------------------

  bubble_plot <- ggplot(
    plot_data,
    aes(
      x = Cluster,
      y = Description
    )
  ) +
    geom_vline(
      xintercept = seq_along(
        group_order
      ),
      colour = "#D3D3D3"
    ) +
    geom_hline(
      yintercept = seq_along(
        description_order
      ),
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
      values = adjust_group_colors,
      limits = adjust_group_levels,
      breaks = adjust_group_levels,
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

  final_plot <- bubble_plot


  # 5.5 Add the optional left annotation panel --------------------------------

  if (!is.null(
    annotation_info$table
  )) {

    term_annotation <- annotation_info$table %>%
      filter(
        Description %in%
          unique(
            as.character(
              plot_data$Description
            )
          )
      ) %>%
      distinct(
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
        n = length(
          annotation_levels
        ),
        palette = "Dark 3"
      ),
      annotation_levels
    )

    left_annotation_plot <- term_annotation %>%
      mutate(
        panel = ""
      ) %>%
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

    final_plot <- aplot::insert_left(
      bubble_plot,
      left_annotation_plot,
      width = 0.08
    )

    write.csv(
      term_annotation,
      file = file.path(
        output_table_dir,
        "GO_term_annotations_used.csv"
      ),
      row.names = FALSE
    )
  }

  save_plot(
    plot_object = final_plot,
    output_dir = output_figure_dir,
    file_stem = "GO_Bubble",
    width = manifest_row$plot_width[
      [1]
    ],
    height = manifest_row$plot_height[
      [1]
    ]
  )


  # 5.6 Save analysis parameters ----------------------------------------------

  analysis_parameters <- tibble(
    parameter = c(
      "analysis_id",
      "deg_file",
      "group_column",
      "gene_column",
      "logfc_column",
      "adjusted_p_column",
      "logfc_cutoff",
      "adjusted_p_cutoff",
      "go_ontology",
      "p_adjust_method",
      "enrichment_p_cutoff",
      "enrichment_q_cutoff",
      "min_gene_set_size",
      "max_gene_set_size",
      "term_annotation_file",
      "restrict_plot_to_annotation_terms",
      "plot_width",
      "plot_height"
    ),
    value = c(
      analysis_id,
      deg_file,
      group_column,
      gene_column,
      logfc_column,
      adjusted_p_column,
      as.character(
        manifest_row$logfc_cutoff[
          [1]
        ]
      ),
      as.character(
        manifest_row$adjusted_p_cutoff[
          [1]
        ]
      ),
      manifest_row$go_ontology[
        [1]
      ],
      manifest_row$p_adjust_method[
        [1]
      ],
      as.character(
        manifest_row$enrichment_p_cutoff[
          [1]
        ]
      ),
      as.character(
        manifest_row$enrichment_q_cutoff[
          [1]
        ]
      ),
      as.character(
        manifest_row$min_gene_set_size[
          [1]
        ]
      ),
      as.character(
        manifest_row$max_gene_set_size[
          [1]
        ]
      ),
      manifest_row$term_annotation_file[
        [1]
      ],
      as.character(
        manifest_row$restrict_plot_to_annotation_terms[
          [1]
        ]
      ),
      as.character(
        manifest_row$plot_width[
          [1]
        ]
      ),
      as.character(
        manifest_row$plot_height[
          [1]
        ]
      )
    )
  )

  write.csv(
    analysis_parameters,
    file = file.path(
      output_table_dir,
      "GO_analysis_parameters.csv"
    ),
    row.names = FALSE
  )

  tibble(
    analysis_id = analysis_id,
    status = "completed",
    deg_file = deg_file,
    n_selected_genes = nrow(
      markers
    ),
    n_mapped_genes = n_distinct(
      markers_mapped$gene
    ),
    n_groups = length(
      gene_clusters
    ),
    n_enriched_terms = nrow(
      ego_df
    ),
    n_plotted_terms = n_distinct(
      plot_data$Description
    )
  )
}


# 6. Load and validate the GO analysis manifest --------------------------------

if (!file.exists(
  analysis_manifest_file
)) {

  stop(
    "GO analysis manifest does not exist: ",
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


# 7. Run all configured GO analyses --------------------------------------------

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
    run_one_go_analysis(
      manifest_row = manifest_row
    ),
    error = function(e) {

      message(
        "GO analysis failed: ",
        analysis_id,
        "\nReason: ",
        conditionMessage(
          e
        )
      )

      tibble(
        analysis_id = analysis_id,
        status = "failed",
        deg_file = manifest_row$deg_file[
          [1]
        ],
        n_selected_genes = NA_integer_,
        n_mapped_genes = NA_integer_,
        n_groups = NA_integer_,
        n_enriched_terms = NA_integer_,
        n_plotted_terms = NA_integer_,
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
    "06_GO_analysis_summary.csv"
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
    "06_GO_enrichment_sessionInfo.txt"
  )
)

message(
  "Configured GO enrichment analyses have finished."
)

print(
  analysis_summary
)
