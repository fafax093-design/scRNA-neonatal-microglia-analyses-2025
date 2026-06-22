```r
# ==============================================================================
# 01. Load and merge single-cell RNA-sequencing datasets
#
# Description:
# This script loads neonatal mouse brain immune-cell datasets from
# preconstructed Seurat objects or 10X Genomics filtered feature-barcode
# matrices. It adds standardized sample metadata, assigns unique cell barcodes,
# merges all samples, and saves the merged raw Seurat object for downstream
# quality-control analysis.
#
# Author:
# Jinjin Zhu
#
# Inputs:
#   DATA/INPUT/RDS/
#   DATA/INPUT/10X/
#
# Outputs:
#   DATA/OUTPUT/sample_information.csv
#   DATA/OUTPUT/sample_loading_summary.csv
#   DATA/OUTPUT/merged_raw_Seurat_object.rds
#   DATA/OUTPUT/01_data_loading_sessionInfo.txt
# ==============================================================================


# 1. Load packages -------------------------------------------------------------

suppressPackageStartupMessages({
    library(Seurat)
    library(dplyr)
    library(tibble)
})


# 2. Define project directories ------------------------------------------------

rds_dir <- "DATA/INPUT/RDS"
matrix_dir <- "DATA/INPUT/10X"
output_dir <- "DATA/OUTPUT"

dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
)


# 3. Define sample information -------------------------------------------------

# Sample information is entered row by row to prevent mismatched vector lengths.
#
# IMPORTANT:
# Replace "REPLACE_WITH_NS_P3_3_FILE.rds" with the actual input file for
# sample NS_P3_3 before running or publishing this script.

```r
# 3. Define sample information -------------------------------------------------

sample_information <- tribble(
    ~sample_id,  ~treatment, ~age,  ~replicate, ~input_type, ~input_path,

    "LPS_P3_1",  "LPS",      "P3",  1L,         "rds",
    file.path(
        rds_dir,
        "LPS_P3_1_seurat.rds"
    ),

    "LPS_P3_2",  "LPS",      "P3",  2L,         "rds",
    file.path(
        rds_dir,
        "LPS_P3_2_seurat.rds"
    ),

    "NS_P3_1",   "NS",       "P3",  1L,         "rds",
    file.path(
        rds_dir,
        "NS_P3_1_seurat.rds"
    ),

    "NS_P3_2",   "NS",       "P3",  2L,         "rds",
    file.path(
        rds_dir,
        "NS_P3_2_seurat.rds"
    ),

    "NS_P3_3",   "NS",       "P3",  3L,         "rds",
    file.path(
        rds_dir,
        "NS_P3_3_seurat.rds"
    ),

    "LPS_P7_1",  "LPS",      "P7",  1L,         "10x",
    file.path(
        matrix_dir,
        "LPS_P7_1",
        "filtered_feature_bc_matrix"
    ),

    "LPS_P7_2",  "LPS",      "P7",  2L,         "10x",
    file.path(
        matrix_dir,
        "LPS_P7_2",
        "filtered_feature_bc_matrix"
    ),

    "NS_P7_1",   "NS",       "P7",  1L,         "10x",
    file.path(
        matrix_dir,
        "NS_P7_1",
        "filtered_feature_bc_matrix"
    ),

    "NS_P7_2",   "NS",       "P7",  2L,         "10x",
    file.path(
        matrix_dir,
        "NS_P7_2",
        "filtered_feature_bc_matrix"
    ),

    "LPS_P12_1", "LPS",      "P12", 1L,         "10x",
    file.path(
        matrix_dir,
        "LPS_P12_1",
        "filtered_feature_bc_matrix"
    ),

    "LPS_P12_2", "LPS",      "P12", 2L,         "10x",
    file.path(
        matrix_dir,
        "LPS_P12_2",
        "filtered_feature_bc_matrix"
    ),

    "NS_P12_1",  "NS",       "P12", 1L,         "10x",
    file.path(
        matrix_dir,
        "NS_P12_1",
        "filtered_feature_bc_matrix"
    ),

    "NS_P12_2",  "NS",       "P12", 2L,         "10x",
    file.path(
        matrix_dir,
        "NS_P12_2",
        "filtered_feature_bc_matrix"
    )
)
```



# 4. Validate sample information -----------------------------------------------

required_sample_columns <- c(
    "sample_id",
    "treatment",
    "age",
    "replicate",
    "input_type",
    "input_path"
)

missing_sample_columns <- setdiff(
    required_sample_columns,
    colnames(sample_information)
)

if (length(missing_sample_columns) > 0) {
    stop(
        "sample_information is missing the following columns: ",
        paste(
            missing_sample_columns,
            collapse = ", "
        )
    )
}

if (anyDuplicated(sample_information$sample_id)) {
    duplicated_sample_ids <- unique(
        sample_information$sample_id[
            duplicated(sample_information$sample_id)
        ]
    )

    stop(
        "Duplicated sample IDs were detected: ",
        paste(
            duplicated_sample_ids,
            collapse = ", "
        )
    )
}

valid_input_types <- c(
    "rds",
    "10x"
)

invalid_input_types <- setdiff(
    unique(
        tolower(sample_information$input_type)
    ),
    valid_input_types
)

if (length(invalid_input_types) > 0) {
    stop(
        "Unsupported input types were detected: ",
        paste(
            invalid_input_types,
            collapse = ", "
        )
    )
}

if (anyNA(sample_information$input_path) ||
    any(sample_information$input_path == "")) {
    stop("Missing input paths were detected in sample_information.")
}

missing_input_paths <- sample_information %>%
    filter(
        !file.exists(input_path)
    )

if (nrow(missing_input_paths) > 0) {
    print(
        missing_input_paths %>%
            select(
                sample_id,
                input_type,
                input_path
            )
    )

    stop(
        "One or more input files or directories do not exist. ",
        "Check the paths shown above."
    )
}

write.csv(
    sample_information,
    file = file.path(
        output_dir,
        "sample_information.csv"
    ),
    row.names = FALSE
)


# 5. Define a function for loading one sample ----------------------------------

load_single_sample <- function(
        sample_row,
        min_cells = 3,
        min_features = 200
) {

    sample_id <- sample_row$sample_id[[1]]
    treatment <- sample_row$treatment[[1]]
    age <- sample_row$age[[1]]
    replicate <- sample_row$replicate[[1]]
    input_type <- tolower(
        sample_row$input_type[[1]]
    )
    input_path <- sample_row$input_path[[1]]

    message(
        "Loading sample: ",
        sample_id
    )

    if (input_type == "rds") {

        seurat_object <- readRDS(
            input_path
        )

        if (!inherits(seurat_object, "Seurat")) {
            stop(
                "The RDS file for sample '",
                sample_id,
                "' does not contain a Seurat object."
            )
        }

    } else if (input_type == "10x") {

        count_matrix <- Read10X(
            data.dir = input_path
        )

        # Read10X may return a list when several feature types are present.
        if (is.list(count_matrix)) {

            if ("Gene Expression" %in% names(count_matrix)) {

                count_matrix <- count_matrix[
                    ["Gene Expression"]
                ]

            } else {

                warning(
                    "No feature type named 'Gene Expression' was found for ",
                    sample_id,
                    ". The first matrix returned by Read10X was used."
                )

                count_matrix <- count_matrix[[1]]
            }
        }

        seurat_object <- CreateSeuratObject(
            counts = count_matrix,
            project = sample_id,
            min.cells = min_cells,
            min.features = min_features
        )

    } else {

        stop(
            "Unsupported input type for sample '",
            sample_id,
            "': ",
            input_type
        )
    }

    if (ncol(seurat_object) == 0) {
        stop(
            "No cells were loaded for sample: ",
            sample_id
        )
    }

    if (!"RNA" %in% Assays(seurat_object)) {
        stop(
            "The RNA assay is absent from sample: ",
            sample_id
        )
    }

    # Prefix cell barcodes so that every cell name remains unique after merging.
    seurat_object <- RenameCells(
        object = seurat_object,
        add.cell.id = sample_id
    )

    # Add standardized sample-level metadata.
    seurat_object$orig.ident <- sample_id
    seurat_object$sample_id <- sample_id
    seurat_object$treatment <- treatment
    seurat_object$age <- age
    seurat_object$replicate <- as.integer(
        replicate
    )
    seurat_object$group <- paste(
        treatment,
        age,
        sep = "_"
    )

    return(seurat_object)
}


# 6. Load all samples -----------------------------------------------------------

seurat_list <- lapply(
    seq_len(
        nrow(sample_information)
    ),
    function(i) {

        load_single_sample(
            sample_row = sample_information[
                i,
                ,
                drop = FALSE
            ]
        )
    }
)

names(seurat_list) <- sample_information$sample_id


# 7. Validate the loaded Seurat objects ----------------------------------------

loaded_sample_ids <- names(
    seurat_list
)

if (!identical(
    loaded_sample_ids,
    sample_information$sample_id
)) {
    stop(
        "The loaded-object order does not match sample_information."
    )
}

all_cell_barcodes <- unlist(
    lapply(
        seurat_list,
        colnames
    ),
    use.names = FALSE
)

if (anyDuplicated(all_cell_barcodes)) {
    stop(
        "Duplicated cell barcodes remain after sample-prefix assignment."
    )
}


# 8. Generate the sample-loading summary ---------------------------------------

sample_summary <- bind_rows(
    lapply(
        names(seurat_list),
        function(sample_name) {

            current_object <- seurat_list[
                [sample_name]
            ]

            tibble(
                sample_id = sample_name,
                treatment = unique(
                    as.character(
                        current_object$treatment
                    )
                ),
                age = unique(
                    as.character(
                        current_object$age
                    )
                ),
                replicate = unique(
                    current_object$replicate
                ),
                n_cells = ncol(
                    current_object
                ),
                n_features = nrow(
                    current_object
                )
            )
        }
    )
)

print(
    sample_summary
)

write.csv(
    sample_summary,
    file = file.path(
        output_dir,
        "sample_loading_summary.csv"
    ),
    row.names = FALSE
)


# 9. Merge all samples ----------------------------------------------------------

scRNA_raw <- merge(
    x = seurat_list[[1]],
    y = seurat_list[-1],
    project = "Neonatal_brain_immune_cells"
)

if (ncol(scRNA_raw) != sum(sample_summary$n_cells)) {
    stop(
        "The number of cells in the merged object does not equal ",
        "the sum of cells in the individual objects."
    )
}

message(
    "Merged object: ",
    ncol(scRNA_raw),
    " cells and ",
    nrow(scRNA_raw),
    " features."
)

print(
    table(
        scRNA_raw$age,
        scRNA_raw$treatment
    )
)

print(
    table(
        scRNA_raw$sample_id
    )
)


# 10. Save the merged object and software information --------------------------

saveRDS(
    scRNA_raw,
    file = file.path(
        output_dir,
        "merged_raw_Seurat_object.rds"
    ),
    compress = "gzip"
)

writeLines(
    capture.output(
        sessionInfo()
    ),
    con = file.path(
        output_dir,
        "01_data_loading_sessionInfo.txt"
    )
)

message(
    "Data loading and merging completed successfully."
)

message(
    "Merged Seurat object saved to: ",
    file.path(
        output_dir,
        "merged_raw_Seurat_object.rds"
    )
)
```
