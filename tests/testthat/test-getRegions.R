data("lai2019")

test_that("Previous result is reproduced.", {
  #Previously obtained results for YAL042W

  combs <- list(between_combs = data.frame(c("A3", "B1", "B2")),
                within_combs = data.frame(c("A1", "A2", "A3")))
  dWithin <- dCombs(rl[["YAL042W"]], combs$within_combs)
  dBetween <- dCombs(rl[["YAL042W"]], combs$between_combs)
  regs <- getContigRegions(getRegions(dWithin,
                              dBetween, rl[["YAL042W"]], min_length = 21))

  expect_true(all(regs$Start == c(19, 530, 583, 1122, 1187)))
  expect_true(all(regs$Stop == c(146, 553, 632, 1177, 1218)))
})
