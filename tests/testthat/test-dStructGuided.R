data(wan2014)

test_that("Previous result is reproduced when testing a region.", {
  #Previously obtained result for 1__NM_000146__357__367
  nowResult <- dStructGuided(wan2014[[1]],
                              reps_A = 2, reps_B = 1)

  expect_true(all(names(nowResult) == c("pval", "del_d")))

  expect_equal(nowResult[1], c(pval = 0.16015625), tolerance = 10^-1)
  expect_equal(nowResult[2], c(del_d = 0.02291958), tolerance = 10^-1)
})
