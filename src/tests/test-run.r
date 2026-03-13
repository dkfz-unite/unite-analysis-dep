library(testthat)
library(readr)
library(jsonlite)

repo_root <- getwd()
run_script <- file.path(repo_root, "src", "run", "run.R")
fixture_dir <- file.path(repo_root, "src", "tests", "test_data", "test_run")

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
	setwd(repo_root)

	output <- system2(
		command = "Rscript",
		args = c(run_script, data_path, metadata_path, results_path, options_path),
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

required_result_columns <- c("feature", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")

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

test_that("run.R works when metadata has no batch column", {
	data_df <- read_tsv(file.path(fixture_dir, "data.tsv"), show_col_types = FALSE)
	metadata_df <- read_tsv(file.path(fixture_dir, "metadata.tsv"), show_col_types = FALSE)
	metadata_no_batch <- metadata_df[, c("sample", "condition")]

	options_list <- fromJSON(file.path(fixture_dir, "options.json"), simplifyVector = TRUE)
	options_list$normalization_method <- "quantile"
	options_list$imputation_method <- "MinDet"
	options_list$stratify_imputation_by_batch <- FALSE
	options_list$min_non_na_fraction <- 0.4
	options_list$min_frac_in_one_class <- TRUE

	run <- run_analysis(data_df, metadata_no_batch, options_list)

	expect_equal(run$status, 0, info = paste(run$output, collapse = "\n"))
	expect_true(file.exists(run$results_path))

	result_df <- read_tsv(run$results_path, show_col_types = FALSE)
	expect_equal(colnames(result_df), required_result_columns)
	expect_gt(nrow(result_df), 0)
	expect_false(any(is.na(result_df$feature)))
})

test_that("run.R works with alternative options configuration", {
	data_df <- read_tsv(file.path(fixture_dir, "data.tsv"), show_col_types = FALSE)
	metadata_df <- read_tsv(file.path(fixture_dir, "metadata.tsv"), show_col_types = FALSE)

	options_list <- fromJSON(file.path(fixture_dir, "options.json"), simplifyVector = TRUE)
	options_list$normalization_method <- "median"
	options_list$imputation_method <- "MinProb"
	options_list$stratify_imputation_by_batch <- FALSE
	options_list$log_offset <- 1
	options_list$min_non_na_fraction <- 0.5
	options_list$min_frac_in_one_class <- FALSE

	run <- run_analysis(data_df, metadata_df, options_list)

	expect_equal(run$status, 0, info = paste(run$output, collapse = "\n"))
	expect_true(file.exists(run$results_path))

	result_df <- read_tsv(run$results_path, show_col_types = FALSE)
	expect_equal(colnames(result_df), required_result_columns)
	expect_gt(nrow(result_df), 0)
	expect_true(all(is.finite(result_df$P.Value)))
})
