# ==============================================================================
# 04. Marker-gene identification, differential expression, and annotation
#     across multiple Seurat objects
#
# Description:
# This script processes one or more Seurat objects defined in an external
# analysis manifest. For each object, it can:
#
#   1. identify marker genes for graph-based clusters;
#   2. apply a reviewed cluster-to-cell-state annotation table;
#   3. identify marker genes for the annotated cell states;
#   4. perform one or more prespecified differential-expression comparisons;
#   5. export cell-level annotations and save the annotated Seurat object.
#
# Different objects may use different clustering columns, annotation files,
# marker assays, output files, and DEG comparisons.
#
# Author:
# Jinjin Zhu
#
# Configuration files:
#   config/04_marker_annotation_manifest.csv
#   config/04_DEG_comparisons.csv
#
# Annotation files:
#   config/annotations/*.csv
#
# Main outputs:
#   DATA/OUTPUT/*.rds
#   TABLE/DEG/<analysis_id>/
# ==============================================================================


# 1. Load packages -------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(future)
})


# 2. Define project paths -------------------------------------------------------

analysis_manifest_file <- "config/04_marker_annotation_manifest.csv"
deg_comparison_file <- "config/04_DEG_comparisons.csv"

output_table_root <- "TABLE/DEG"

dir.create(
  output_table_root,
  recursive = TRUE,
  showWarnings = FALSE
)


# 3. Define global analysis parameters -----------------------------------------

random_seed <- 1234

cluster_marker_logfc_threshold <- 0.25
cluster_marker_min_pct <- 0.10
cluster_marker_only_positive <- FALSE

celltype_marker_logfc_threshold <- 0.25
celltype_marker_min_pct <- 0.10
celltype_marker_only_positive <- FALSE

maximum_future_size_gib <- 20

workers <- min(
  10L,
  max(
    1L,
    parallel::detectCores(
      logical = TRUE
    ) - 1L
  )
)

set.seed(
  random_seed
)

options(
  future.globals.maxSize =
    maximum_future_size_gib *
    1024^3
)

future::plan(
  future::multisession,
  workers = workers
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


sort_cluster_ids <- function(cluster_ids) {

  cluster_ids <- unique(
    as.character(
      cluster_ids
    )
  )

  numeric_ids <- suppressWarnings(
    as.numeric(
      cluster_ids
    )
  )

  if (all(
    !is.na(
      numeric_ids
    )
  )) {

    return(
      cluster_ids[
        order(
          numeric_ids
        )
      ]
    )
  }

  sort(
    cluster_ids
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


write_annotation_template <- function(
  cluster_ids,
  output_file
) {

  template <- tibble(
    cluster = cluster_ids,
    celltype = ""
  )

  dir.create(
    dirname(
      output_file
    ),
    recursive = TRUE,
    showWarnings = FALSE
  )

  write.csv(
    template,
    file = output_file,
    row.names = FALSE,
    na = ""
  )

  invisible(
    template
  )
}


validate_annotation_table <- function(
  annotation_table,
  cluster_ids
) {

  required_columns <- c(
    "cluster",
    "celltype"
  )

  missing_columns <- setdiff(
    required_columns,
    colnames(
      annotation_table
    )
  )

  if (length(
    missing_columns
  ) > 0) {

    stop(
      "The annotation table is missing columns: ",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  }

  annotation_table <- annotation_table %>%
    transmute(
      cluster = trimws(
        as.character(
          cluster
        )
      ),
      celltype = trimws(
        as.character(
          celltype
        )
      )
    )

  if (anyDuplicated(
    annotation_table$cluster
  )) {

    duplicated_clusters <- unique(
      annotation_table$cluster[
        duplicated(
          annotation_table$cluster
        )
      ]
    )

    stop(
      "Duplicated cluster IDs were detected: ",
      paste(
        duplicated_clusters,
        collapse = ", "
      )
    )
  }

  empty_labels <- annotation_table %>%
    filter(
      is.na(
        celltype
      ) |
        celltype == ""
    ) %>%
    pull(
      cluster
    )

  if (length(
    empty_labels
  ) > 0) {

    stop(
      "Empty cell-state labels were detected for cluster(s): ",
      paste(
        empty_labels,
        collapse = ", "
      )
    )
  }

  missing_clusters <- setdiff(
    cluster_ids,
    annotation_table$cluster
  )

  if (length(
    missing_clusters
  ) > 0) {

    stop(
      "The annotation table does not contain cluster(s): ",
      paste(
        missing_clusters,
        collapse = ", "
      )
    )
  }

  extra_clusters <- setdiff(
    annotation_table$cluster,
    cluster_ids
  )

  if (length(
    extra_clusters
  ) > 0) {

    warning(
      "The annotation table contains cluster IDs absent from the object: ",
      paste(
        extra_clusters,
        collapse = ", "
      )
    )
  }

  annotation_table %>%
    filter(
      cluster %in% cluster_ids
    )
}


get_required_manifest_columns <- function() {

  c(
    "analysis_id",
    "input_file",
    "output_object_file",
    "cluster_column",
    "annotation_file",
    "celltype_column",
    "marker_assay",
    "marker_slot",
    "prepare_sct_find_markers",
    "run_cluster_markers",
    "run_celltype_markers"
  )
}


validate_manifest <- function(manifest) {

  missing_columns <- setdiff(
    get_required_manifest_columns(),
    colnames(
      manifest
    )
  )

  if (length(
    missing_columns
  ) > 0) {

    stop(
      "The analysis manifest is missing columns: ",
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
          output_object_file,
          cluster_column,
          annotation_file,
          celltype_column,
          marker_assay,
          marker_slot
        ),
        ~ trimws(
          as.character(
            .x
          )
        )
      ),
      prepare_sct_find_markers = as_flag(
        prepare_sct_find_markers
      ),
      run_cluster_markers = as_flag(
        run_cluster_markers
      ),
      run_celltype_markers = as_flag(
        run_celltype_markers
      )
    )

  if (any(
    manifest$analysis_id == "" |
      is.na(
        manifest$analysis_id
      )
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


read_deg_plan <- function(file) {

  if (!file.exists(
    file
  )) {

    return(
      tibble(
        analysis_id = character(),
        comparison_id = character(),
        identity_column = character(),
        ident_1 = character(),
        ident_2 = character(),
        assay = character(),
        slot = character(),
        test_use = character(),
        logfc_threshold = numeric(),
        min_pct = numeric()
      )
    )
  }

  deg_plan <- read.csv(
    file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  required_columns <- c(
    "analysis_id",
    "comparison_id",
    "identity_column",
    "ident_1",
    "ident_2",
    "assay",
    "slot",
    "test_use",
    "logfc_threshold",
    "min_pct"
  )

  missing_columns <- setdiff(
    required_columns,
    colnames(
      deg_plan
    )
  )

  if (length(
    missing_columns
  ) > 0) {

    stop(
      "The DEG comparison file is missing columns: ",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  }

  deg_plan %>%
    mutate(
      across(
        c(
          analysis_id,
          comparison_id,
          identity_column,
          ident_1,
          ident_2,
          assay,
          slot,
          test_use
        ),
        ~ trimws(
          as.character(
            .x
          )
        )
      ),
      logfc_threshold = as.numeric(
        logfc_threshold
      ),
      min_pct = as.numeric(
        min_pct
      )
    ) %>%
    filter(
      analysis_id != "",
      comparison_id != ""
    )
}


run_deg_comparisons <- function(
  object,
  analysis_id,
  comparison_plan,
  output_dir
) {

  current_plan <- comparison_plan %>%
    filter(
      analysis_id == !!analysis_id
    )

  if (nrow(
    current_plan
  ) == 0) {

    message(
      "No DEG comparisons were configured for analysis: ",
      analysis_id
    )

    return(
      invisible(
        NULL
      )
    )
  }

  comparison_summary <- vector(
    mode = "list",
    length = nrow(
      current_plan
    )
  )

  for (i in seq_len(
    nrow(
      current_plan
    )
  )) {

    comparison_row <- current_plan[
      i,
      ,
      drop = FALSE
    ]

    identity_column <- comparison_row$identity_column[
      [1]
    ]

    if (!identity_column %in% colnames(
      object@meta.data
    )) {

      stop(
        "Identity column '",
        identity_column,
        "' was not found for comparison '",
        comparison_row$comparison_id[
          [1]
        ],
        "'."
      )
    }

    assay_name <- comparison_row$assay[
      [1]
    ]

    if (!assay_name %in% Assays(
      object
    )) {

      stop(
        "Assay '",
        assay_name,
        "' was not found for comparison '",
        comparison_row$comparison_id[
          [1]
        ],
        "'."
      )
    }

    Idents(
      object
    ) <- identity_column

    available_identities <- levels(
      Idents(
        object
      )
    )

    requested_identities <- c(
      comparison_row$ident_1[
        [1]
      ],
      comparison_row$ident_2[
        [1]
      ]
    )

    missing_identities <- setdiff(
      requested_identities,
      available_identities
    )

    if (length(
      missing_identities
    ) > 0) {

      stop(
        "Comparison '",
        comparison_row$comparison_id[
          [1]
        ],
        "' refers to absent identity level(s): ",
        paste(
          missing_identities,
          collapse = ", "
        )
      )
    }

    result <- FindMarkers(
      object = object,
      ident.1 = comparison_row$ident_1[
        [1]
      ],
      ident.2 = comparison_row$ident_2[
        [1]
      ],
      assay = assay_name,
      slot = comparison_row$slot[
        [1]
      ],
      test.use = comparison_row$test_use[
        [1]
      ],
      logfc.threshold = comparison_row$logfc_threshold[
        [1]
      ],
      min.pct = comparison_row$min_pct[
        [1]
      ],
      random.seed = random_seed,
      verbose = TRUE
    ) %>%
      rownames_to_column(
        "gene"
      )

    comparison_id_safe <- sanitize_name(
      comparison_row$comparison_id[
        [1]
      ]
    )

    write.csv(
      result,
      file = file.path(
        output_dir,
        paste0(
          "DEG_",
          comparison_id_safe,
          ".csv"
        )
      ),
      row.names = FALSE
    )

    comparison_summary[
      [i]
    ] <- tibble(
      comparison_id = comparison_row$comparison_id[
        [1]
      ],
      identity_column = identity_column,
      ident_1 = comparison_row$ident_1[
        [1]
      ],
      ident_2 = comparison_row$ident_2[
        [1]
      ],
      assay = assay_name,
      slot = comparison_row$slot[
        [1]
      ],
      test_use = comparison_row$test_use[
        [1]
      ],
      logfc_threshold = comparison_row$logfc_threshold[
        [1]
      ],
      min_pct = comparison_row$min_pct[
        [1]
      ],
      n_genes_returned = nrow(
        result
      )
    )
  }

  comparison_summary <- bind_rows(
    comparison_summary
  )

  write.csv(
    comparison_summary,
    file = file.path(
      output_dir,
      "DEG_comparison_summary.csv"
    ),
    row.names = FALSE
  )

  invisible(
    comparison_summary
  )
}


process_one_object <- function(
  manifest_row,
  deg_plan
) {

  analysis_id <- manifest_row$analysis_id[
    [1]
  ]

  analysis_id_safe <- sanitize_name(
    analysis_id
  )

  output_dir <- file.path(
    output_table_root,
    analysis_id_safe
  )

  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  input_file <- manifest_row$input_file[
    [1]
  ]

  output_object_file <- manifest_row$output_object_file[
    [1]
  ]

  cluster_column <- manifest_row$cluster_column[
    [1]
  ]

  annotation_file <- manifest_row$annotation_file[
    [1]
  ]

  celltype_column <- manifest_row$celltype_column[
    [1]
  ]

  marker_assay <- manifest_row$marker_assay[
    [1]
  ]

  marker_slot <- manifest_row$marker_slot[
    [1]
  ]

  prepare_sct <- manifest_row$prepare_sct_find_markers[
    [1]
  ]

  run_cluster_markers <- manifest_row$run_cluster_markers[
    [1]
  ]

  run_celltype_markers <- manifest_row$run_celltype_markers[
    [1]
  ]

  message(
    "\nProcessing analysis: ",
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

  if (!cluster_column %in% colnames(
    object@meta.data
  )) {

    stop(
      "Cluster column '",
      cluster_column,
      "' was not found in analysis '",
      analysis_id,
      "'."
    )
  }

  if (!marker_assay %in% Assays(
    object
  )) {

    stop(
      "Marker assay '",
      marker_assay,
      "' was not found in analysis '",
      analysis_id,
      "'."
    )
  }

  cluster_ids <- sort_cluster_ids(
    object@meta.data[
      [cluster_column]
    ]
  )

  cluster_size_table <- object@meta.data %>%
    transmute(
      cluster = as.character(
        .data[
          [cluster_column]
        ]
      )
    ) %>%
    count(
      cluster,
      name = "n_cells"
    ) %>%
    mutate(
      cluster = factor(
        cluster,
        levels = cluster_ids
      )
    ) %>%
    arrange(
      cluster
    ) %>%
    mutate(
      cluster = as.character(
        cluster
      )
    )

  write.csv(
    cluster_size_table,
    file = file.path(
      output_dir,
      "cluster_sizes.csv"
    ),
    row.names = FALSE
  )

  if (prepare_sct &&
      marker_assay == "SCT") {

    object <- PrepSCTFindMarkers(
      object = object,
      assay = marker_assay
    )
  }

  if (run_cluster_markers) {

    Idents(
      object
    ) <- cluster_column

    cluster_markers <- FindAllMarkers(
      object = object,
      assay = marker_assay,
      slot = marker_slot,
      test.use = "wilcox",
      logfc.threshold = cluster_marker_logfc_threshold,
      min.pct = cluster_marker_min_pct,
      only.pos = cluster_marker_only_positive,
      random.seed = random_seed,
      verbose = TRUE
    )

    write.csv(
      cluster_markers,
      file = file.path(
        output_dir,
        "cluster_markers.csv"
      ),
      row.names = FALSE
    )
  }

  generated_template_file <- file.path(
    output_dir,
    paste0(
      analysis_id_safe,
      "_cluster_annotation_template.csv"
    )
  )

  if (!file.exists(
    annotation_file
  )) {

    write_annotation_template(
      cluster_ids = cluster_ids,
      output_file = generated_template_file
    )

    stop(
      "Annotation file does not exist: ",
      annotation_file,
      "\nA template was generated at: ",
      generated_template_file
    )
  }

  annotation_table <- read.csv(
    annotation_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  annotation_table <- validate_annotation_table(
    annotation_table = annotation_table,
    cluster_ids = cluster_ids
  )

  annotation_table$cluster <- factor(
    annotation_table$cluster,
    levels = cluster_ids
  )

  annotation_table <- annotation_table %>%
    arrange(
      cluster
    ) %>%
    mutate(
      cluster = as.character(
        cluster
      )
    )

  write.csv(
    annotation_table,
    file = file.path(
      output_dir,
      "cluster_annotations_used.csv"
    ),
    row.names = FALSE
  )

  cluster_to_celltype <- setNames(
    annotation_table$celltype,
    annotation_table$cluster
  )

  object[
    [celltype_column]
  ] <- unname(
    cluster_to_celltype[
      as.character(
        object@meta.data[
          [cluster_column]
        ]
      )
    ]
  )

  if (anyNA(
    object@meta.data[
      [celltype_column]
    ]
  )) {

    stop(
      "Missing labels were generated in column: ",
      celltype_column
    )
  }

  celltype_levels <- unique(
    annotation_table$celltype
  )

  object[
    [celltype_column]
  ] <- factor(
    object@meta.data[
      [celltype_column]
    ],
    levels = celltype_levels
  )

  cluster_celltype_counts <- table(
    cluster = object@meta.data[
      [cluster_column]
    ],
    celltype = object@meta.data[
      [celltype_column]
    ]
  )

  write.csv(
    as.data.frame(
      cluster_celltype_counts
    ),
    file = file.path(
      output_dir,
      "cluster_to_celltype_cell_counts.csv"
    ),
    row.names = FALSE
  )

  if (run_celltype_markers) {

    Idents(
      object
    ) <- celltype_column

    celltype_markers <- FindAllMarkers(
      object = object,
      assay = marker_assay,
      slot = marker_slot,
      test.use = "wilcox",
      logfc.threshold = celltype_marker_logfc_threshold,
      min.pct = celltype_marker_min_pct,
      only.pos = celltype_marker_only_positive,
      random.seed = random_seed,
      verbose = TRUE
    )

    write.csv(
      celltype_markers,
      file = file.path(
        output_dir,
        "celltype_markers.csv"
      ),
      row.names = FALSE
    )
  }

  run_deg_comparisons(
    object = object,
    analysis_id = analysis_id,
    comparison_plan = deg_plan,
    output_dir = output_dir
  )

  preferred_metadata <- c(
    "sample_id",
    "orig.ident",
    "treatment",
    "age",
    "group",
    "batch",
    cluster_column,
    celltype_column
  )

  metadata_columns <- intersect(
    preferred_metadata,
    colnames(
      object@meta.data
    )
  )

  cell_annotation_metadata <- object@meta.data %>%
    rownames_to_column(
      "cell_barcode"
    ) %>%
    select(
      cell_barcode,
      all_of(
        metadata_columns
      )
    )

  write.csv(
    cell_annotation_metadata,
    file = file.path(
      output_dir,
      "cell_cluster_and_annotation_metadata.csv"
    ),
    row.names = FALSE
  )

  dir.create(
    dirname(
      output_object_file
    ),
    recursive = TRUE,
    showWarnings = FALSE
  )

  saveRDS(
    object,
    file = output_object_file,
    compress = "gzip"
  )

  parameter_table <- tibble(
    parameter = c(
      "analysis_id",
      "input_file",
      "output_object_file",
      "cluster_column",
      "annotation_file",
      "celltype_column",
      "marker_assay",
      "marker_slot",
      "prepare_sct_find_markers",
      "run_cluster_markers",
      "run_celltype_markers",
      "cluster_marker_logfc_threshold",
      "cluster_marker_min_pct",
      "celltype_marker_logfc_threshold",
      "celltype_marker_min_pct",
      "random_seed"
    ),
    value = c(
      analysis_id,
      input_file,
      output_object_file,
      cluster_column,
      annotation_file,
      celltype_column,
      marker_assay,
      marker_slot,
      as.character(
        prepare_sct
      ),
      as.character(
        run_cluster_markers
      ),
      as.character(
        run_celltype_markers
      ),
      as.character(
        cluster_marker_logfc_threshold
      ),
      as.character(
        cluster_marker_min_pct
      ),
      as.character(
        celltype_marker_logfc_threshold
      ),
      as.character(
        celltype_marker_min_pct
      ),
      as.character(
        random_seed
      )
    )
  )

  write.csv(
    parameter_table,
    file = file.path(
      output_dir,
      "analysis_parameters.csv"
    ),
    row.names = FALSE
  )

  tibble(
    analysis_id = analysis_id,
    status = "completed",
    input_file = input_file,
    output_object_file = output_object_file,
    n_cells = ncol(
      object
    ),
    n_features = nrow(
      object
    ),
    n_clusters = length(
      cluster_ids
    ),
    n_celltypes = length(
      celltype_levels
    )
  )
}


# 5. Load and validate the analysis manifest -----------------------------------

if (!file.exists(
  analysis_manifest_file
)) {

  future::plan(
    future::sequential
  )

  stop(
    "Analysis manifest does not exist: ",
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

deg_plan <- read_deg_plan(
  deg_comparison_file
)


# 6. Process all configured Seurat objects -------------------------------------

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
    process_one_object(
      manifest_row = manifest_row,
      deg_plan = deg_plan
    ),
    error = function(e) {

      message(
        "Analysis failed: ",
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
        output_object_file = manifest_row$output_object_file[
          [1]
        ],
        n_cells = NA_integer_,
        n_features = NA_integer_,
        n_clusters = NA_integer_,
        n_celltypes = NA_integer_,
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
    "04_analysis_summary.csv"
  ),
  row.names = FALSE
)


# 7. Save software information -------------------------------------------------

writeLines(
  capture.output(
    sessionInfo()
  ),
  con = file.path(
    output_table_root,
    "04_marker_DEG_annotation_sessionInfo.txt"
  )
)

future::plan(
  future::sequential
)

message(
  "Configured marker, DEG, and annotation analyses have finished."
)

print(
  analysis_summary
)
