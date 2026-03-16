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


# ─── batch handling ───────────────────────────────────────────────────────────

test_that("da_analysis with batch returns same structure as without", {
  d <- make_da_data(n_samples_per_cond = 6)
  batch <- rep(c("1", "2", "3"), times = 4)  # 12 samples, 3 batches
  result_no_batch   <- da_analysis(d$data, d$condition)
  result_with_batch <- da_analysis(d$data, d$condition, batch = batch)
  expect_equal(nrow(result_with_batch), nrow(result_no_batch))
  expect_equal(colnames(result_with_batch), colnames(result_no_batch))
})
