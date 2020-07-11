context('wrapper tests')

test_that("mpspline_compact",
          c( s1 <-
               data.frame("SID" = c( 1,  1,  1,  1,    2,   2,   2,  2,   3,    4,  5,  5),
                          "UD"  = c( 0, 20, 30, 50,   -1,  45,  15, 80,   0,   30,  0, 30),
                          "LD"  = c(20, 30, 50, 70,    5,  60,  30, NA,  10,   50, 10, 50),
                          "VAL" = c( 6,  4,  3, 100, 0.1, 0.9, 2.5,  6, 3.5, 10.4, NA, NA),
                          stringsAsFactors = FALSE),
             m1 <-  mpspline_compact(s1, var_name = 'VAL', lam = 0.1,
                                     d = c(0, 5, 15, 30, 60, 100, 200),
                                     vhigh = 14, vlow = 0),
             expect_is(m1, 'list'),
             expect_equal(length(m1), 4),
             expect_is(m1[[1]], 'matrix'),
             expect_is(m1[[2]], 'matrix'),
             expect_is(m1[[3]], 'matrix'),
             expect_is(m1[[4]], 'matrix'),
             # var missing
             m2 <-  mpspline_compact(s1, lam = 0.1,
                                     d = c(0, 5, 15, 30, 60, 100, 200),
                                     vhigh = 14, vlow = 0),
             expect_identical(m1, m2)

          )
)


test_that("mpspline_tidy",
          c( s1 <-
               data.frame("SID" = c( 1,  1,  1,  1,    2,   2,   2,  2,   3,    4,  5,  5),
                          "UD"  = c( 0, 20, 30, 50,   -1,  45,  15, 80,   0,   30,  0, 30),
                          "LD"  = c(20, 30, 50, 70,    5,  60,  30, NA,  10,   50, 10, 50),
                          "VAL" = c( 6,  4,  3, 100, 0.1, 0.9, 2.5,  6, 3.5, 10.4, NA, NA),
                          stringsAsFactors = FALSE),
             m1 <-  mpspline_tidy(s1, var_name = 'VAL', lam = 0.1,
                                     d = c(0, 5, 15, 30, 60, 100, 200),
                                     vhigh = 14, vlow = 0),
             expect_is(m1, 'list'),
             expect_equal(length(m1), 4),
             expect_is(m1[[1]], 'data.frame'),
             expect_is(m1[[2]], 'data.frame'),
             expect_is(m1[[3]], 'data.frame'),
             expect_is(m1[[4]], 'data.frame'),
             expect_is(m1[[4]][[2]], 'character'),
             # var missing,
             m2 <-  mpspline_tidy(s1, lam = 0.1,
                                     d = c(0, 5, 15, 30, 60, 100, 200),
                                     vhigh = 14, vlow = 0),
             expect_identical(m1, m2)

          )
)
