
test_that("Returns expected values of normalized reactivities", {
  expect_true(all(two.eight.normalize(c(NA, 1:20, NA, -999)) ==
                c(NA, (1:20)/19, NA, -999), na.rm = TRUE))
  expect_true(all(is.na(two.eight.normalize(c(NA, 1:20, NA, -999))[c(1, 22)])))

  expect_true(all(two.eight.normalize(c(NA, 1:10, NA, -999, NaN)) ==
                    c(NA, (1:10)/9.5, NA, -999, NaN), na.rm = TRUE))
  expect_true(all(is.na(two.eight.normalize(c(NA, 1:10, NA, -999, NaN))[c(1, 12, 14)])))

  expect_true(all(is.nan(two.eight.normalize(c(NA, 1:10, NA, -999, NaN))[14])))
})
