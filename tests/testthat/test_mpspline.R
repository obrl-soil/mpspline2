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
            expect_equal(nrow(fixed), 8),
            # whoops, site IDs aren't numeric!
            obj <- matrix(c("a", "a", "a", 0, 10, 20, 10, 20, 30, 4.5, 6, 7.8),
                          ncol = 4, byrow = FALSE),
            fx <- mpspline_conv(obj),
            expect_is(fx, 'data.frame'),
            expect_equal(nrow(fx), 3),
            expect_is(fx[[2]], 'numeric')
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
            chkd <- mpspline_datchk(sites, 'VAL'),
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
            expect_message(mpspline_datchk(allna, 'VAL')),
            expect_equal(mpspline_datchk(allna, 'VAL')[[1]], NA),
            # overlap
            ols <- list("A" = data.frame("SID" = "A",
                                         "UD"  = c(0, 30, 50, 90),
                                         "LD"  = c(30, 60, 100 ,120),
                                         "VAL" = c(1,2,3,4),
                                         stringsAsFactors = FALSE)),
            expect_message(mpspline_datchk(ols, 'VAL')),
            expect_equal(mpspline_datchk(ols, 'VAL')[[1]], NA),
            # single horizons
            sh <-  list("A" = data.frame("SID" = "A",
                                         "UD"  =   0,
                                         "LD"  =  10,
                                         "VAL" =   4,
                                         stringsAsFactors = FALSE),
                        "B" = data.frame("SID" = "B",
                                         "UD"  =  60,
                                         "LD"  =  80,
                                         "VAL" =   5,
                                         stringsAsFactors = FALSE)),
            chkd <- mpspline_datchk(sh, 'VAL'),
            expect_is(chkd, 'list'),
            expect_equal(length(chkd), 2),
            expect_equal(nrow(chkd[[1]]), 1)
          ))

test_that("mpspline_est1 does the thing",
          c(
            # normal
            s <- data.frame("SID" = "A",
                            "UD"  = c(0,10,20,30),
                            "LD"  = c(10,20,30,40),
                            "VAL" = c(5.4, 5.3, 5.6, 7.0),
                            stringsAsFactors = FALSE),
            spar <- mpspline_est1(s, 'VAL', lam = 0.1),
            # lock this down
            expect_is(spar, 'list'),
            expect_equal(length(spar), 6),
            expect_equal(names(spar), c("s_bar", "b0", "b1", "gamma", "alfa", "Z")),
            expect_is(spar[[1]], 'numeric'),
            expect_equal(length(spar[[1]]), 4),
            expect_equivalent(spar[[1]][1], 5.3938765799840498),
            expect_is(spar[[2]], 'numeric'),
            expect_equal(length(spar[[2]]), 4),
            expect_equivalent(spar[[2]][1], 0),
            expect_is(spar[[3]], 'numeric'),
            expect_equal(length(spar[[3]]), 4),
            expect_equivalent(spar[[3]][1], -0.015308550039877671),
            expect_is(spar[[4]], 'numeric'),
            expect_equal(length(spar[[4]]), 4),
            expect_equivalent(spar[[4]][1], -0.00076542750199388358),
            expect_is(spar[[5]], 'numeric'),
            expect_equal(length(spar[[5]]), 4),
            expect_equivalent(spar[[5]][1], 5.4193908300505127),
            expect_is(spar[[6]], 'matrix'),
            expect_equal(nrow(spar[[6]]), 4),
            expect_equivalent(spar[[6]][1], 1.064285714),
            # one horizon
            s <- data.frame("SID" = "A", "UD" = 0, "LD"  = 10, "VAL" = 5.4,
                            stringsAsFactors = FALSE),
            spar <- mpspline_est1(s, 'VAL', lam = 0.1),
            expect_equal(spar, NA_real_)
            )
          )

test_that("mpspline_fit1 does the thing",
          c(
            # single horizon
            s1 <- data.frame("SID" = "A", "UD" = 0, "LD" = 10, "VAL" = 5.4),
            p1 <- mpspline_est1(s1, 'VAL', lam = 0.1),
            expect_message(mpspline_fit1(s = s1, p = p1, var_name = 'VAL',
                                          d = c(0, 5, 15, 30, 60, 100, 200),
                                          vhigh = 14, vlow = 0)),
            f1 <- mpspline_fit1(s = s1, p = p1, var_name = 'VAL',
                                 d = c(0, 5, 15, 30, 60, 100, 200),
                                 vhigh = 14, vlow = 0),
            expect_is(f1, 'list'),
            expect_equal(length(f1), 2),
            expect_equal(f1[[1]][1], s1[[4]]),
            expect_equal(f1[[1]][11], NA_real_),
            expect_equal(f1[[2]][1], s1[[4]]),
            expect_equal(f1[[2]][3], NA_real_),
            # normal no gaps
            s2 <- data.frame("SID" = c( 1,  1,  1,  1),
                             "UD"  = c( 0, 20, 30, 50),
                             "LD"  = c(20, 30, 50, 70),
                             "VAL" = c( 6,  4,  3, 10),
                              stringsAsFactors = FALSE),
            p2 <- mpspline_est1(s2, 'VAL', lam = 0.1),
            f2 <- mpspline_fit1(s = s2, p = p2, var_name = 'VAL',
                                 d = c(0, 5, 15, 30, 60, 100, 200),
                                 vhigh = 14, vlow = 0),
            expect_is(f2, 'list'),
            expect_is(f2[[1]], 'numeric'),
            expect_is(f2[[2]], 'numeric'),
            expect_equal(length(f2), 2),
            expect_equal(length(f2[[1]]), 200),
            expect_true(sum(is.na(f2[[1]])) == 130),
            expect_true(max(f2[[1]], na.rm = TRUE) <= 14),
            expect_true(min(f2[[1]], na.rm = TRUE) >= 0),
            expect_equal(length(f2[[2]]), 6),
            expect_true(sum(is.na(f2[[2]])) == 1),
            # gaps
            s3 <- data.frame("SID" = c( 1,  1,  1),
                             "UD"  = c( 0, 30, 50),
                             "LD"  = c(20, 50, 70),
                             "VAL" = c( 6,  3, 10),
                             stringsAsFactors = FALSE),
            p3 <- mpspline_est1(s3, 'VAL', lam = 0.1),
            f3 <- mpspline_fit1(s = s3, p = p3, var_name = 'VAL',
                                d = c(0, 5, 15, 30, 60, 100, 200),
                                vhigh = 14, vlow = 0),
            # still getting predictions in the same depth range as s1:
            expect_true(sum(is.na(f3[[1]])) == 130),
            expect_true(max(f3[[1]], na.rm = TRUE) <= 14),
            expect_true(min(f3[[1]], na.rm = TRUE) >= 0),
            expect_equal(length(f3[[2]]), 6),
            expect_true(sum(is.na(f3[[2]])) == 1),
            # deeper than max d - should just truncate
            f4 <- mpspline_fit1(s = s2, p = p2, var_name = 'VAL',
                                d = c(0, 5, 15, 30),
                                vhigh = 14, vlow = 0),
            expect_equal(length(f4[[1]]), 30),
            expect_true(sum(is.na(f4[[1]])) == 0),
            expect_equal(length(f4[[2]]), 3),
            expect_true(sum(is.na(f4[[2]])) == 0),
            expect_equal(f2[[2]][1:3], f4[[2]]),
            # starts below surface
            s4 <- data.frame("SID" = c( 1,  1,  1),
                             "UD"  = c(20, 30, 50),
                             "LD"  = c(30, 50, 70),
                             "VAL" = c( 4,  3, 10),
                             stringsAsFactors = FALSE),
            p4 <- mpspline_est1(s4, 'VAL', lam = 0.1),
            f5 <- mpspline_fit1(s = s4, p = p4, var_name = 'VAL',
                                d = c(0, 5, 15, 30, 60, 100),
                                vhigh = 14, vlow = 0),
            # no extrapolation, only interpolation!!
            expect_true(all(is.na(f5[[1]][1:20]))),
            expect_true(all(is.na(f5[[1]][71:100]))) # x[70] == 69-70cm
          ))

test_that("mpspline_tmse1 does the thing",
          c( s1 <- data.frame("SID" = c( 1,  1,  1,  1),
                              "UD"  = c( 0, 20, 30, 50),
                              "LD"  = c(20, 30, 50, 70),
                              "VAL" = c( 6,  4,  3, 10),
                              stringsAsFactors = FALSE),
             p1 <- mpspline_est1(s1, 'VAL', lam = 0.1),
             f1 <- mpspline_fit1(s = s1, p = p1, var_name = 'VAL',
                                 d = c(0, 5, 15, 30, 60, 100, 200),
                                 vhigh = 14, vlow = 0),
             s_hat_5 <- (0.05 * stats::sd(s1[[4]], na.rm = TRUE))^2,
             var_5 <- s_hat_5^2,
             t1 <- mpspline_tmse1(s1, p1, var_name = 'VAL', s2 = var_5),
             expect_equal(t1, 0.036610703692220463),
             p <- list("s_bar" = NA, "b0" = NA, "b1" = NA, "gamma" = NA,
                       "alfa" = NA, "Z" = NA),
             expect_equal(mpspline_tmse1(s1, p, var_name = 'VAL',
                                         s2 = var_5), NA_real_)
          )
        )

test_that("mpspline_tmse1 does the thing",
          c( s1 <-
               data.frame("SID" = c( 1,  1,  1,  1,    2,   2,   2,  2,   3,    4,  5,  5),
                          "UD"  = c( 0, 20, 30, 50,   -1,  45,  15, 80,   0,   30,  0, 30),
                          "LD"  = c(20, 30, 50, 70,    5,  60,  30, NA,  10,   50, 10, 50),
                          "VAL" = c( 6,  4,  3, 100, 0.1, 0.9, 2.5,  6, 3.5, 10.4, NA, NA),
                              stringsAsFactors = FALSE),
             m1 <- mpspline(s1, var_name = 'VAL',
                            d = c(0, 5, 15, 30, 60, 100, 200),
                            vhigh = 14, vlow = 0),
             # var name skipped
             expect_message(mpspline(s1, d = c(0, 5, 15, 30, 60, 100, 200),
                                     vhigh = 14, vlow = 0)),
             m2 <- mpspline(s1, d = c(0, 5, 15, 30, 60, 100, 200),
                            vhigh = 14, vlow = 0),
             expect_identical(m1, m2)
          )
)
