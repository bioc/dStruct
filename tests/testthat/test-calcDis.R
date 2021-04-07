
test_that("For numeric vector input, return is a single value", {
  expect_equal(NaN, calcDis(c(0, 0)))
  expect_equal(1, calcDis(c(1, -1)))
  expect_equal(0, calcDis(c(1, 1)))
  expect_true(
    length(
      calcDis(
        runif(10, -10, 10)
        )
      ) == 1
    )
})

test_that("Output is always between 0 and 1", {
  expect_true(
    all(
      replicate(10,
                calcDis(runif(10,
                              runif(1, -20, 0),
                              runif(1, 0, 20)
                )
                )
      ) <= 1
    )
  )

  expect_true(
    all(
      replicate(10,
                calcDis(runif(10,
                              runif(1, -20, 0),
                              runif(1, 0, 20)
                )
                )
      ) >= 0
    )
  )

})

test_that("Vector of length the same as number of rows of matrix input is returned", {
  reacs <- data.frame(matrix(runif(30, 0, 10), 10, 3))
  expect_true(length(calcDis(reacs)) == nrow(reacs))
})
