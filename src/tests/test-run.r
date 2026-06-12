library(testthat)
library(readr)
library(jsonlite)
library(dplyr)

run_dir <- getwd() # assume tests are executed from inside "run"
run_script <- file.path(run_dir,"run.R")
fixture_dir <- file.path(run_dir, "..", "tests", "test_data", "test_run")

expect_true(file.exists(run_script))
expect_true(file.exists(file.path(fixture_dir, "data.tsv")))
expect_true(file.exists(file.path(fixture_dir, "metadata.tsv")))
expect_true(file.exists(file.path(fixture_dir, "options.json")))

run_analysis <- function(data_df, metadata_df, options_list) {
	work_dir <- tempfile("run_r_test_")
	dir.create(work_dir, recursive = TRUE)

	data_path <- file.path(work_dir, "data.tsv")
	metadata_path <- file.path(work_dir, "metadata.tsv")
	options_path <- file.path(work_dir, "options.json")
	results_path <- file.path(work_dir, "results.tsv")

	write_tsv(data_df, data_path)
	write_tsv(metadata_df, metadata_path)
	write_json(options_list, options_path, auto_unbox = TRUE, pretty = TRUE, null = "null")

	old_wd <- getwd()
	on.exit(setwd(old_wd), add = TRUE)
	setwd(run_dir)

	output <- system2(
		command = "Rscript",
		args = c(run_script, work_dir),
		stdout = TRUE,
		stderr = TRUE
	)
	status <- attr(output, "status")
	if (is.null(status)) {
		status <- 0
	}

	list(
		status = status,
		output = output,
		results_path = results_path
	)
}

required_result_columns <- c("feature", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "contrast")

test_that("run.R reproduces expected bundled results with batch metadata", {
	data_df <- read_tsv(file.path(fixture_dir, "data.tsv"), show_col_types = FALSE)
	metadata_df <- read_tsv(file.path(fixture_dir, "metadata.tsv"), show_col_types = FALSE)
	options_list <- fromJSON(file.path(fixture_dir, "options.json"), simplifyVector = TRUE)
	expected_df <- read_tsv(file.path(fixture_dir, "results.tsv"), show_col_types = FALSE)

	run <- run_analysis(data_df, metadata_df, options_list)

	expect_equal(run$status, 0, info = paste(run$output, collapse = "\n"))
	expect_true(file.exists(run$results_path))

	result_df <- read_tsv(run$results_path, show_col_types = FALSE)
	expect_equal(colnames(result_df), required_result_columns)
	expect_equal(result_df, expected_df)
})

test_that("run.R works when metadata has empty batch column", {
	data_df <- read_tsv(file.path(fixture_dir, "data.tsv"), show_col_types = FALSE)
	metadata_df <- read_tsv(file.path(fixture_dir, "metadata.tsv"), show_col_types = FALSE)
	# make batch column all empty values
	metadata_no_batch <- metadata_df %>%
		mutate(batch = "")

	options_list <- fromJSON(file.path(fixture_dir, "options.json"), simplifyVector = TRUE)
	options_list$normalization_method <- "quantile"
	options_list$imputation_method <- "MinDet"
	options_list$stratify_imputation_by_batch <- FALSE
	options_list$min_non_missing_fraction <- 0.4
	options_list$require_min_fraction_one_class <- TRUE

	run <- run_analysis(data_df, metadata_no_batch, options_list)

	expect_equal(run$status, 0, info = paste(run$output, collapse = "\n"))
	expect_true(file.exists(run$results_path))

	result_df <- read_tsv(run$results_path, show_col_types = FALSE)
	expect_equal(colnames(result_df), required_result_columns)
	expect_gt(nrow(result_df), 0)
	expect_false(any(is.na(result_df$feature)))
})

test_that("last-occurring condition in metadata is the reference category", {
	# feature1 is higher in A, feature2 higher in B; metadata ends with "B".
	# With B as reference: logFC = A - B, so feature1 logFC > 0, feature2 logFC < 0.
	# this test rquires the reference category be set explicitly in the script as the defualt 
	# of as.factor is to set the first alphabetically as the ref category
	# If the reference were A instead, both signs would flip.
	data_df <- tibble::tibble(
		feature = c("feature1", "feature2", "feature3"),
		s1 = c(100, 1, 10),
		s2 = c(100, 1, 10),
		s3 = c(1, 100, 10),
		s4 = c(1, 100, 10)
	)
	metadata_df <- tibble::tibble(
		sample     = c("s1", "s2", "s3", "s4"),
		condition  = c("A",  "A",  "B",  "B"),
		batch      = c("1",  "1",  "1",  "1")
	)
	options_list <- list(
		normalization_method          = "median",
		normalization_log_offset      = 1,
		imputation_method             = "MinDet",
		stratify_imputation_by_batch  = FALSE,
		batch_correction_method       = NULL,
		min_non_missing_fraction      = 0,
		require_min_fraction_one_class = FALSE
	)

	run <- run_analysis(data_df, metadata_df, options_list)

	expect_equal(run$status, 0, info = paste(run$output, collapse = "\n"))
	result_df <- read_tsv(run$results_path, show_col_types = FALSE)

	f1_logFC <- result_df$logFC[result_df$feature == "feature1"]
	f2_logFC <- result_df$logFC[result_df$feature == "feature2"]

	# B is reference (last in metadata), so logFC = A - B
	expect_gt(f1_logFC, 0)  # feature1 higher in A
	expect_lt(f2_logFC, 0)  # feature2 higher in B
})

test_that("run.R succeeds when all samples share the same batch value", {
	# Regression test: a single-level batch factor must not be passed to da_analysis,
	# as model.matrix cannot compute contrasts for a factor with only one level.
	data_df <- tibble::tibble(
		feature = c("feature1", "feature2", "feature3"),
		s1 = c(100, 1, 10),
		s2 = c(100, 1, 10),
		s3 = c(1, 100, 10),
		s4 = c(1, 100, 10)
	)
	metadata_df <- tibble::tibble(
		sample    = c("s1", "s2", "s3", "s4"),
		condition = c("A",  "A",  "B",  "B"),
		batch     = c("batch1", "batch1", "batch1", "batch1")
	)
	options_list <- list(
		normalization_method          = "median",
		normalization_log_offset      = 1,
		imputation_method             = "MinDet",
		stratify_imputation_by_batch  = FALSE,
		batch_correction_method       = NULL,
		min_non_missing_fraction      = 0,
		require_min_fraction_one_class = FALSE
	)

	run <- run_analysis(data_df, metadata_df, options_list)

	expect_equal(run$status, 0, info = paste(run$output, collapse = "\n"))
	result_df <- read_tsv(run$results_path, show_col_types = FALSE)
	expect_equal(colnames(result_df), required_result_columns)
	expect_gt(nrow(result_df), 0)
})

test_that("run.R works with alternative options configuration", {
	data_df <- read_tsv(file.path(fixture_dir, "data.tsv"), show_col_types = FALSE)
	metadata_df <- read_tsv(file.path(fixture_dir, "metadata.tsv"), show_col_types = FALSE)

	options_list <- fromJSON(file.path(fixture_dir, "options.json"), simplifyVector = TRUE)
	options_list$normalization_method <- "median"
	options_list$imputation_method <- "MinProb"
	options_list$stratify_imputation_by_batch <- FALSE
	options_list$normalization_log_offset <- 1
	options_list$min_non_missing_fraction <- 0.5
	options_list$require_min_fraction_one_class <- FALSE

	run <- run_analysis(data_df, metadata_df, options_list)

	expect_equal(run$status, 0, info = paste(run$output, collapse = "\n"))
	expect_true(file.exists(run$results_path))

	result_df <- read_tsv(run$results_path, show_col_types = FALSE)
	expect_equal(colnames(result_df), required_result_columns)
	expect_gt(nrow(result_df), 0)
	expect_true(all(is.finite(result_df$P.Value)))
})
