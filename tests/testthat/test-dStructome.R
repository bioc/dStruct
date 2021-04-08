data("lai2019")
data("wan2014")

test_that("Previous result is reproduced in de novo mode", {
  transcripts_previously_top <- c("YJR045C", "YOR383C", "YHL033C",
                                  "YMR120C", "YJR009C", "YJL189W")
  lai2019 <- lai2019[c(transcripts_previously_top,
                       sample(names(lai2019[!(names(lai2019) %in% transcripts_previously_top)]), 6))]
  nowResult <- suppressWarnings(dStructome(lai2019, 3, 2, batches= TRUE,
                                           min_length = 21,
                                           between_combs = data.frame(c("A3", "B1", "B2")),
                                           within_combs = data.frame(c("A1", "A2", "A3")),
                                           ind_regions = TRUE, processes = 1))
  nowResult <- nowResult[order(nowResult$FDR), ]
  nowTop <- head(nowResult, n = 6)
  expect_true(all(nowTop$t ==
                    transcripts_previously_top))
  expect_true(all(nowTop$Start == c(15, 76, 19, 1557, 1073, 115)))
  expect_true(all(nowTop$Stop == c(203, 250, 107, 1741, 1237, 203)))

  pvals_to_compare <- data.frame(nowResult = nowTop$pval,
                                 previous = c(2.286121e-10, 1.762228e-09,
                                              1.925648e-08, 1.662185e-08,
                                              1.062762e-07, 1.723124e-07))
  apply(pvals_to_compare, 1, function(x) {
    names(x) <- NULL
    expect_equal(x[1], x[2],
                 tolerance = 10^-6)
  })

  delD_to_compare <- data.frame(nowResult = nowTop$del_d,
                                previous = c(0.09315098,
                                             0.16484672, 0.17185732,
                                             0.12450723, 0.11853069,
                                             0.14850375))
  apply(delD_to_compare, 1, function(x) {
    names(x) <- NULL
    expect_equal(x[1], x[2],
                 tolerance = 10^-6)
  })

})

test_that("Previous result is reproduced in guided mode", {
  previous_result <- data.frame(t = c("1__NM_000146__357__367",
                                      "2__ENST00000527880.1__281__291",
                                      "3__NM_020820__2667__2677",
                                      "1__NM_002482__1624__1634",
                                      "3__NM_032855__1812__1822",
                                      "2__NM_199440__385__395",
                                      "3__NM_004517__241__251",
                                      "3__NM_004068__1337__1347",
                                      "3__NM_002349__4109__4119"),
                                pval = c(0.16015625, 0.20654296875,
                                         0.08740234375, 0.07373046875,
                                         0.04150390625, 0.05078125,
                                         0.00927734375, 0.00341796875,
                                         0.0379634814912788),
                                del_d = c(0.0229195757396641,
                                          0.0906367583431284, 0.0721098887050756,
                                          0.0198376721176467, 0.0770827662548119,
                                          0.0917944930071851, 0.221674090640147,
                                          0.168966119900885, 0.164412994040479),
                                FDR = c(0.18017578125, 0.20654296875,
                                        0.112374441964286, 0.110595703125,
                                        0.09140625, 0.09140625,
                                        0.041748046875, 0.03076171875,
                                        0.09140625))

  nowResult <- dStructome(wan2014,
                          reps_A = 2, reps_B = 1, method = "guided",
                          processes = 1)

  apply(cbind(previous_result$t, nowResult$t),
        1,
        function(x) expect_true(x[1] == x[2]))

  apply(cbind(previous_result$pval, nowResult$pval),
        1,
        function(x) expect_equal(x[1], x[2], tolerance = 10^-6))

  apply(cbind(previous_result$del_d, nowResult$del_d),
        1,
        function(x) expect_equal(x[1], x[2], tolerance = 10^-6))

  apply(cbind(previous_result$FDR, nowResult$FDR),
        1,
        function(x) expect_equal(x[1], x[2], tolerance = 10^-6))
})
