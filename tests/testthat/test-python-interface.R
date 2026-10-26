test_that("run_genecluster validates matrix dimensions before Python starts", {
  invalid <- matrix(1, nrow = 2, ncol = 3,
                    dimnames = list(c("g1", "g2"), c("c1", "c2", "c3")))
  expect_error(
    run_genecluster(invalid, tempfile(), n_clusters = 2),
    "at least three genes"
  )
})

test_that("run_genecluster validates pseudo-cluster count before Python starts", {
  valid <- matrix(1, nrow = 3, ncol = 3,
                  dimnames = list(paste0("g", 1:3), paste0("c", 1:3)))
  expect_error(
    run_genecluster(valid, tempfile(), n_clusters = 3),
    "between 2 and nrow"
  )
})
