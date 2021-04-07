
test_that("Return is a list", {
  expect_true(is.list(getCombs(2, 1)) & !is.data.frame(getCombs(2, 1)))
  expect_true(is.list(getCombs(3, 4, batches = TRUE)) &
                !is.data.frame(getCombs(3, 4, batches = TRUE)))
  expect_true(is.list(getCombs(2, 1, batches = TRUE)) &
                !is.data.frame(getCombs(2, 1,
                                        between_combs = data.frame("A1", "B1"),
                                        within_combs = data.frame("A1", "A2"))))
  expect_true(all(names(getCombs(2, 4)) == c("between_combs", "within_combs")))
})

test_that("Elements of the returned list have equal number of rows", {
  tester <- function(test_next) expect_true(nrow(test_next[[1]]) == nrow(test_next[[2]]))

  test_next <- getCombs(2, 1)
  tester(test_next)

  test_next <- getCombs(4, 3, batches = TRUE)
  tester(test_next)

  test_next <- getCombs(7, 7)
  tester(test_next)
})

test_that("Homogeneous sets have samples of the same group", {
  tester <- function(test_next) expect_true(all(
    apply(test_next$within_combs, 2,
          function(x) length(unique(substr(x, 1, 1)))
          == 1)))

  test_next <- getCombs(2, 1)
  tester(test_next)

  test_next <- getCombs(4, 3, batches = TRUE)
  tester(test_next)

  test_next <- getCombs(7, 7)
  tester(test_next)
})

test_that("Heterogeneous sets have samples of both groups", {
  tester <- function(test_next) expect_true(all(
    apply(test_next$between_combs, 2,
          function(x) length(unique(substr(x, 1, 1)))
          == 2)))

  test_next <- getCombs(2, 1)
  tester(test_next)

  test_next <- getCombs(4, 3, batches = TRUE)
  tester(test_next)

  test_next <- getCombs(7, 7)
  tester(test_next)
})
