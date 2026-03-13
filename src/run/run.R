library(limma)
library(readr)
library(dplyr)
library(tibble)
library(jsonlite)
source(file.path(getwd(),"src", "helpers", "preprocessing.r"))
source(file.path(getwd(),"src", "helpers", "da_analysis.r"))

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
data_file <- args[1]
metadata_file <- args[2]
results_file <- args[3]
options_file <- args[4]

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

# reorder metadata rows to match data matrix columns
metadata_matrix <- metadata_matrix[match(colnames(data_matrix), rownames(metadata_matrix)), ,drop = FALSE]

# transpose data_matrix as proteomic_data_preprocessing expects samples as rows, features as columns
data_matrix <- t(data_matrix)

# process options - batch correct method must be NULL for differential expression analysis
options <- replace_required(options, "batch_correction_method", NULL)

# if column "batch" is present in metadata, use it for batch correction, otherwise set to NULL
if ("batch" %in% colnames(metadata_matrix)) {
  batch_vector <- as.factor(metadata_matrix$batch)
} else {
  batch_vector <- NULL
}

# preprocess data
processed_data <- preprocess_data(data=data_matrix, 
                              batch_vector=batch_vector,
                              class_labels=as.factor(metadata_matrix$condition), 
                              options = options)


# Get results
results <- da_analysis(t(processed_data), condition=as.factor(metadata_matrix$condition), batch=batch_vector)
#rownames(results) <- colnames(processed_data) # ensure feature names are preserved in results
# Write results to file
results <- rownames_to_column(results, var = "feature")
write_tsv(results, results_file)
