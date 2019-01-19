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

