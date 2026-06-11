library(limma)
library(readr)
library(dplyr)
library(tibble)
library(jsonlite)
source(file.path(getwd(), "helpers", "preprocessing.r"))
source(file.path(getwd(), "helpers", "da_analysis.r"))

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
workdir <- args[1]
metadata_file <- file.path(workdir, "metadata.tsv")
results_file <- file.path(workdir, "results.tsv")
data_file <- file.path(workdir, "data.tsv")
options_file <- file.path(workdir, "options.json")

# Read data and metadata
data <- read_tsv(data_file)
metadata <- read_tsv(metadata_file)
options <- fromJSON(options_file)
print(names(options))
# Preprocess data (log, normalize, impute, batch correct)
data_matrix <- as.data.frame(data[,-1]) # remove feature column for processing
rownames(data_matrix) <- data[[1]] # set feature names as rownames

metadata_matrix <- as.data.frame(metadata[, -1]) # assuming first column is sample names
rownames(metadata_matrix) <- metadata[[1]] # set sample names as rownames

# keep reference category for later - this should be the last category listed
ref_cat<-tail(metadata_matrix$condition, 1))

# reorder metadata rows to match data matrix columns
metadata_matrix <- metadata_matrix[match(colnames(data_matrix), rownames(metadata_matrix)), ,drop = FALSE]

# transpose data_matrix as proteomic_data_preprocessing expects samples as rows, features as columns
data_matrix <- t(data_matrix)

# process options - batch correct method must be NULL for differential expression analysis
options <- replace_required(options, "batch_correction_method", NULL)

# preprocess data
processed_data <- preprocess_data(data=data_matrix, 
                              batch_vector=metadata_matrix$batch,
                              class_labels=metadata_matrix$condition, 
                              options = options)

# check if batch column is empty, all NA, or single-level; if so set batch_vector to NULL for da_analysis
if (any(is.na(metadata_matrix$batch)) || any(metadata_matrix$batch == "") || length(unique(metadata_matrix$batch)) < 2) {
  batch_vector <- NULL
} else {
  batch_vector <- as.factor(metadata_matrix$batch)
}

# Get results
# relevel ensures that the last (second) category is always the reference category
results <- da_analysis(t(processed_data), condition=relevel(as.factor(metadata_matrix$condition), ref = as.character(ref_cat), batch=batch_vector)
#rownames(results) <- colnames(processed_data) # ensure feature names are preserved in results
# Write results to file
results <- rownames_to_column(results, var = "feature")
write_tsv(results, results_file)

