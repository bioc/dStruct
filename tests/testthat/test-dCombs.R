
test_that("There is one d-score per nucleotide", {
  reacs <- data.frame(matrix(runif(30, 0, 10), 10, 3))
  colnames(reacs) <- c("A1", "A2", "B1")

  expect_true(length(dCombs(reacs, getCombs(2, 1)$between_combs)) ==
                nrow(reacs))
})

