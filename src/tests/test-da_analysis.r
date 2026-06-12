library(testthat)
library(limma)

source(file.path(getwd() ,"helpers", "da_analysis.r")) # run test from "run" directory

# ─── helpers ──────────────────────────────────────────────────────────────────

make_da_data <- function(n_features = 20, n_samples_per_cond = 5,
                         de_features = 5, effect_size = 5, seed = 42) {
  set.seed(seed)
  n <- n_samples_per_cond * 2
  mat <- matrix(rnorm(n_features * n, mean = 10, sd = 1),
                nrow = n_features,
                dimnames = list(paste0("prot", seq_len(n_features)),
                                paste0("s", seq_len(n))))
  condition <- rep(c("A", "B"), each = n_samples_per_cond)
  # add a clear effect to the first de_features proteins in condition B
  mat[seq_len(de_features), condition == "B"] <-
    mat[seq_len(de_features), condition == "B"] + effect_size
  list(data = mat, condition = condition)
}

# ─── return structure ─────────────────────────────────────────────────────────

test_that("da_analysis returns a data.frame", {
  d <- make_da_data()
  result <- da_analysis(d$data, d$condition)
  expect_s3_class(result, "data.frame")
})


test_that("da_analysis rownames are feature names", {
  d <- make_da_data(n_features = 10)
  result <- da_analysis(d$data, d$condition)
  expect_equal(sort(rownames(result)), sort(rownames(d$data)))
})


# ─── condition validation ─────────────────────────────────────────────────────

test_that("da_analysis errors with more than two conditions", {
  d <- make_da_data()
  condition_three <- c(d$condition, rep("C", 2))
  data_three <- cbind(d$data, matrix(rnorm(nrow(d$data) * 2, mean = 10, sd = 1),
                                     nrow = nrow(d$data),
                                     dimnames = list(rownames(d$data), c("s11", "s12"))))
  expect_error(da_analysis(data_three, condition_three))
})

test_that("da_analysis errors with only one condition", {
  d <- make_da_data()
  expect_error(da_analysis(d$data, rep("A", ncol(d$data))))
})

# ─── reference category ───────────────────────────────────────────────────────

test_that("da_analysis ref_category changes the sign of logFC", {
  d <- make_da_data()
  result_ab <- da_analysis(d$data, d$condition, ref_category = "A")
  result_ba <- da_analysis(d$data, d$condition, ref_category = "B")
  expect_equal(result_ab[rownames(result_ba), "logFC"], -result_ba[, "logFC"])
})

test_that("da_analysis errors when ref_category is not one of the conditions", {
  d <- make_da_data()
  expect_error(da_analysis(d$data, d$condition, ref_category = "Z"))
})

# ─── batch handling ───────────────────────────────────────────────────────────

test_that("da_analysis with batch returns same structure as without", {
  d <- make_da_data(n_samples_per_cond = 6)
  batch <- rep(c("1", "2", "3"), times = 4)  # 12 samples, 3 batches
  result_no_batch   <- da_analysis(d$data, d$condition)
  result_with_batch <- da_analysis(d$data, d$condition, batch = batch)
  expect_equal(nrow(result_with_batch), nrow(result_no_batch))
  expect_equal(colnames(result_with_batch), colnames(result_no_batch))
})
