
test_that("Returns expected values of normalization constants", {
  expect_true(normalizer(c(NA, 1:20, NA, -999)) == 19)

  expect_true(normalizer(c(NA, 1:10, NA, NaN, -999)) == 9.5)
})
