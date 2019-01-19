context('mpspline tests')


test_that("mpspline_convert works for numeric matrices",
          c(
            obj <- matrix(c( 1,  1,  1,  1,   2,   2,   2,   2,
                             0, 20, 40, 60,  -1,  45,  15,  80,
                            10, 30, 50, 70,   5,  60,  30, 100,
                             6,  4,  3, 10, 0.1, 0.9, 2.5,   6),
                          ncol = 4, byrow = FALSE),
            fixed <- mpspline_prep(obj),
            expect_is(fixed, 'data.frame'),
            expect_equal(nrow(fixed), 7),
            expect_equal(fixed[[3]][5], 30) # -1 dropped, site sorted
          )
)

test_that("mpspline_convert works for data frames",
          c(
            obj <- data.frame("SID" = c( 1,  1,  1,  1,   2,   2,   2,   2),
                              "UD"  = c( 0, 20, 40, 60,  -1,  45,  15,  80),
                              "LD"  = c(10, 30, 50, 70,   5,  60,  30, 100),
                              "VAL" = c( 6,  4,  3, 10, 0.1, 0.9, 2.5,   6),
                              stringsAsFactors = FALSE),
            fixed <- mpspline_prep(obj),
            expect_is(fixed, 'data.frame'),
            expect_equal(nrow(fixed), 7),
            expect_equal(fixed[[3]][5], 30) # -1 dropped, site sorted
          ))

test_that("mpspline_convert works for SoilProfileCollections",
          c(library(aqp),
            obj <-  data.frame("SID" = c( 1,  1,  1,  1,   2,   2,   2,   2),
                               "UD"  = c( 0, 20, 40, 60,  -1,  45,  15,  80),
                               "LD"  = c(10, 30, 50, 70,   5,  60,  30, 100),
                               "VAL" = c( 6,  4,  3, 10, 0.1, 0.9, 2.5,   6),
                               stringsAsFactors = FALSE),
            depths(obj) <- SID ~ UD + LD,
            fixed <- mpspline_prep(obj),
            expect_is(fixed, 'data.frame'),
            expect_equal(nrow(fixed), 7),
            expect_equal(fixed[[3]][5], 30)
          ))
