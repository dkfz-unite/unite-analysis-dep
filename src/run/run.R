library(limma)
library(readr)
library(dplyr)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
data_file <- args[1]
metadata_file <- args[2]
results_file <- args[3]

# Read data and metadata
data <- read_tsv(data_file)
metadata <- read_tsv(metadata_file)

# Prepare data for analysis
# First column is feature names, rest are samples
features <- data[[1]]
data_matrix <- as.matrix(data[,-1])
rownames(data_matrix) <- features

# Replace zeros with NA (missing values)
data_matrix[data_matrix == 0] <- NA

# Filter proteins: keep only those with values in at least 50% of samples
valid_proteins <- rowSums(!is.na(data_matrix)) >= (ncol(data_matrix) * 0.5)
data_matrix <- data_matrix[valid_proteins, ]
features <- features[valid_proteins]

# Replace remaining NA values with a small value (half of minimum non-zero value)
min_value <- min(data_matrix, na.rm = TRUE)
data_matrix[is.na(data_matrix)] <- min_value / 2

# Log2 transformation
data_matrix <- log2(data_matrix)

# Quantile normalization
data_matrix <- normalizeBetweenArrays(data_matrix, method = "quantile")

# Align metadata rows to match column order of data matrix
metadata <- metadata[match(colnames(data_matrix), metadata[[1]]), ]

# Create design matrix for Limma
# First column is sample names, second column is condition
conditions <- factor(metadata[[2]])
design <- model.matrix(~ 0 + conditions)
colnames(design) <- levels(conditions)

# Fit the linear model
fit <- lmFit(data_matrix, design)

# Create contrast matrix for comparison (exactly 2 conditions)
contrast_string <- paste0(levels(conditions)[1], "-", levels(conditions)[2])
contrast_matrix <- makeContrasts(contrasts = contrast_string, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# Get results
results <- topTable(fit2, number = Inf)

# Add feature names to results
results <- results %>%
  mutate(Feature = rownames(results)) %>%
  select(Feature, everything())

# Write results to file
write_tsv(results, results_file)
