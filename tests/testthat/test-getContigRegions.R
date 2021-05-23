
test_that("Contiguous region boundaries are properly identified.", {
  tester <- function(x, expected_starts, expected_stops) {
    all(c(IRanges::start(x) == expected_starts,
          IRanges::end(x) == expected_stops))
  }

  test_next <- getContigRegions(c(1:10, 11:20))
  expect_true(tester(test_next, 1, 20))

  test_next <- getContigRegions(c(1:10, 12:20))
  expect_true(tester(test_next, c(1, 12), c(10, 20)))

  test_next <- getContigRegions(c(1:10, 12:20), gap = 1)
  expect_true(tester(test_next, 1, 20))

  test_next <- getContigRegions(c(1:10, 13:20), gap = 1)
  expect_true(tester(test_next, c(1, 13), c(10, 20)))

  test_next <- getContigRegions(c(1:8, 15:20))
  expect_true(tester(test_next, c(1, 15), c(8, 20)))
})


