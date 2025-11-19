context('mpspline multivariable tests')

test_data_multi_het <- data.frame(
  SID = c(1, 1, 1, 1,    2, 2, 2, 2,    3, 3,    4, 4),
  UD  = c(0, 20, 40, 60, 0, 15, 45, 80, 0, 10,   0, 0),
  LD  = c(10, 30, 50, 70, 5, 30, 60, 100, 5, 20, 10, 5),
  VAL1 = c(6, 4, 3, 10,    0.1, 0.9, 2.5, 6,    NA, 100,  NA, NA),
  VAL2 = c(5, 3, 2, 9,     NA, 1.0, 2.0, 5,    10, 15,   NA, NA),
  VAL3 = c(NA, NA, NA, NA, 1, 2, 3, 4,    NA, NA,    50, 60),
  VAL4 = c(1, 2, 3, 4,     5, 6, 7, 8,    9, 10,   11, NA),
  stringsAsFactors = FALSE
)

test_that("Multi-variable wrappers produce correct output structure", c(
  # Test compact wrapper
  compact_multi <- mpspline_compact(test_data_multi_het, var_name = c("VAL1", "VAL2", "VAL3", "VAL4")),
  expect_equal(names(compact_multi), c("VAL1", "VAL2", "VAL3", "VAL4")),
  expect_is(compact_multi$VAL1, 'list'),
  expect_equal(names(compact_multi$VAL1), c("est_icm", "est_1cm", "est_dcm", "tmse")),
  expect_is(compact_multi$VAL1$est_icm, 'matrix'),
  expect_equal(nrow(compact_multi$VAL1$est_icm), 3),  # VAL1: 3 sites with data
  expect_equal(nrow(compact_multi$VAL3$est_icm), 1),  # VAL3: 1 site with data

  # Test tidy wrapper
  tidy_multi <- mpspline_tidy(test_data_multi_het, var_name = c("VAL1", "VAL2", "VAL3", "VAL4")),
  expect_equal(names(tidy_multi), c("est_icm", "est_1cm", "est_dcm", "tmse")),
  expect_true("VARIABLE" %in% colnames(tidy_multi$est_icm)),
  expect_equal(sort(unique(tidy_multi$est_icm$VARIABLE)), c("VAL1", "VAL2", "VAL3", "VAL4")),
  expect_equal(sort(unique(subset(tidy_multi$est_dcm, VARIABLE == "VAL1")$SID)), c(1, 2, 3)),
  expect_equal(sort(unique(subset(tidy_multi$est_dcm, VARIABLE == "VAL3")$SID)), c(2)),

  # Sanity check: multi-variable results match single-variable runs
  compact_val1_single <- mpspline_compact(test_data_multi_het, var_name = "VAL1"),
  expect_equivalent(compact_multi$VAL1$est_1cm, compact_val1_single$est_1cm)
))

test_that("Multi-variable handles variable-specific NA patterns", c(
  result_multi <- mpspline(test_data_multi_het, var_name = c("VAL1", "VAL2", "VAL3", "VAL4")),

  # VAL3 only has data for site 2
  expect_false("3" %in% names(result_multi$VAL3)),

  # VAL4 has partial data for site 4 (single value at depth 0-0)
  expect_true("4" %in% names(result_multi$VAL4)),
  val4_s4 <- result_multi$VAL4[["4"]]$est_1cm,
  expect_equal(val4_s4[1:10], rep(11, 10)),
  expect_true(all(is.na(val4_s4[11:length(val4_s4)])))
))

test_that("Multi-variable respects custom depth intervals and parameters", c(
  # Custom depth intervals
  custom_d <- c(0, 10, 20, 50, 100),
  result_custom <- mpspline_tidy(test_data_multi_het, var_name = c("VAL1", "VAL2"), d = custom_d),
  expect_equal(length(unique(result_custom$est_dcm$UD)), length(custom_d) - 1),
  expect_equal(sort(unique(result_custom$est_dcm$UD)), custom_d[1:(length(custom_d)-1)]),

  # Value constraints applied across variables
  result_constrained <- mpspline_tidy(test_data_multi_het, var_name = c("VAL1", "VAL2"), vlow = 5, vhigh = 8),
  expect_true(all(na.omit(result_constrained$est_1cm$SPLINED_VALUE) >= 5)),
  expect_true(all(na.omit(result_constrained$est_1cm$SPLINED_VALUE) <= 8))
))
