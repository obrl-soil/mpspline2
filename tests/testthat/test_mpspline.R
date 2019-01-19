context('mpspline tests')


test_that("mpspline_conv works for numeric matrices",
          c(
            obj <- matrix(c( 1,  1,  1,  1,   2,   2,   2,   2,
                            NA, 20, 40, 60,  -1,  45,  15,  80,
                            10, 30, 50, 70,   5,  60,  30,  NA,
                             6,  4,  3, 10, 0.1, 0.9, 2.5,   6),
                          ncol = 4, byrow = FALSE),
            fixed <- mpspline_conv(obj),
            expect_is(fixed, 'data.frame'),
            expect_equal(nrow(fixed), 8)
          )
        )

test_that("mpspline_conv works for data frames",
          c(
            obj <- data.frame("SID" = c( 1,  1,  1,  1,   2,   2,   2,   2),
                              "UD"  = c(NA, 20, 40, 60,  -1,  45,  15,  80),
                              "LD"  = c(10, 30, 50, 70,   5,  60,  30,  NA),
                              "VAL" = c( 6,  4,  3, 10, 0.1, 0.9, 2.5,   6),
                              stringsAsFactors = FALSE),
            fixed <- mpspline_conv(obj),
            expect_is(fixed, 'data.frame'),
            expect_equal(obj, fixed)
          )
        )


test_that("mpspline_conv works for SoilProfileCollections",
          c(library(aqp),
            obj <-  data.frame("SID" = c( 1,  1,  1,  1,   2,   2,   2,   2),
                               "UD"  = c(NA, 20, 40, 60,  -1,  45,  15,  80),
                               "LD"  = c(10, 30, 50, 70,   5,  60,  30,  NA),
                               "VAL" = c( 6,  4,  3, 10, 0.1, 0.9, 2.5,   6),
                               stringsAsFactors = FALSE),
            depths(obj) <- SID ~ UD + LD,
            fixed <- mpspline_conv(obj),
            expect_is(fixed, 'data.frame'),
            expect_equal(nrow(fixed), 8)
          )
        )

test_that("mpspline_datchk does what it oughta",
          c(
            obj <-  data.frame("SID" = c( 1,  1,  1,  1,   2,   2,   2,   2),
                               "UD"  = c(NA, 20, 40, 60,  -1,  45,  15,  80),
                               "LD"  = c(10, 30, 50, 70,   5,  60,  30,  NA),
                               "VAL" = c( 6,  4,  3, 10, 0.1, 0.9, 2.5,   6),
                               stringsAsFactors = FALSE),
            obj <- mpspline_conv(obj),
            sites <- split(obj, as.factor(obj[[1]])),
            chkd <- mpspline_datchk(sites),
            expect_is(chkd, 'list'),
            expect_equal(length(chkd), 2),
            expect_is(chkd[[1]], 'data.frame'),
            expect_equal(nrow(chkd[[1]]), 4),
            expect_equal(nrow(chkd[[2]]), 3),  # -ve hor dropped
            expect_equal(chkd[[1]][[2]][1], 0), # surface fixed
            expect_equal(chkd[[2]][[3]][3], 90), # last fixed
            # all na
            allna <- list("A" = data.frame("SID" = "A",
                                           "UD"  = c(0, 30, 60, 90),
                                           "LD"  = c(30, 60, 90 ,120),
                                           "VAL" = NA_real_,
                                           stringsAsFactors = FALSE)),
            expect_error(mpspline_datchk(allna)),
            # overlap
            ols <- list("A" = data.frame("SID" = "A",
                                         "UD"  = c(0, 30, 50, 90),
                                         "LD"  = c(30, 60, 100 ,120),
                                         "VAL" = c(1,2,3,4),
                                         stringsAsFactors = FALSE)),
            expect_error(mpspline_datchk(ols))
          ))
