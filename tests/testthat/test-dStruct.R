data("lai2019")

test_that("Previous result is reproduced when testing individual regions.", {
  #Previously obtained results for YAL042W
  nowResult <- suppressWarnings(dStruct(rdf = lai2019[["YAL042W"]],
                       reps_A = 3, reps_B = 2,
                       batches = TRUE, min_length = 21,
                       between_combs = data.frame(c("A3", "B1", "B2")),
                       within_combs = data.frame(c("A1", "A2", "A3")),
                       ind_regions = TRUE))

  expect_true(all(nowResult$Start == c(19, 530, 583, 1122, 1187)))
  expect_true(all(nowResult$Stop == c(146, 553, 632, 1177, 1218)))

  pvals_to_compare <- data.frame(nowResult = nowResult$pval,
                                 previous = c(0.0001419124,
                                              0.0767881983, 0.0228964197,
                                              0.0353501385, 0.0164846906))
  apply(pvals_to_compare, 1, function(x) {
    names(x) <- NULL
    expect_equal(x[1], x[2],
                 tolerance = 10^-1)
  })

  delD_to_compare <- data.frame(nowResult = nowResult$del_d,
                                previous = c(0.11071808,
                                             0.07939060, 0.01748423,
                                             0.02047984, 0.09202706))
  apply(delD_to_compare, 1, function(x) {
    names(x) <- NULL
    expect_equal(x[1], x[2],
                 tolerance = 10^-1)
  })

  fdr_to_compare <- data.frame(nowResult = nowResult$FDR,
                                previous = c(0.0007095618,
                                             0.0767881983, 0.0381606995,
                                             0.0441876731, 0.0381606995))
  apply(fdr_to_compare, 1, function(x) {
    names(x) <- NULL
    expect_equal(x[1], x[2],
                 tolerance = 10^-1)
  })

})


test_that("Previous result is reproduced when testing regions of a transcript collectively.", {
  #Previously obtained results for YAL042W
  nowResult <- suppressWarnings(dStruct(rdf = lai2019[["YAL042W"]],
                                        reps_A = 3, reps_B = 2,
                                        batches = TRUE, min_length = 21,
                                        between_combs = data.frame(c("A3", "B1", "B2")),
                                        within_combs = data.frame(c("A1", "A2", "A3")),
                                        ind_regions = FALSE))

  expect_true(is.list(nowResult))
  expect_true(all(names(nowResult) == c("regions", "pval", "del_d")))

  expect_true(all(nowResult$regions$Start == c(19, 530, 583, 1122, 1187)))
  expect_true(all(nowResult$regions$Stop == c(146, 553, 632, 1177, 1218)))

  expect_equal(nowResult$pval, 5.997025e-08,
               tolerance = 10^-1)

  expect_equal(nowResult$del_d, 0.09044452,
               tolerance = 10^-1)

})
