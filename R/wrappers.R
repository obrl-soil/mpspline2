#' Unified wrapper base for spline output formatting
#' @keywords internal
#' @noRd
.mpspline_wrapper <- function(obj, var_name, lam, d, vlow, vhigh, formatter_fn, combiner_fn = NULL) {
  obj <- mpspline_conv(obj)

  if (is.null(var_name)) {
    message("Parameter var_name not supplied, assuming target data is in column 4.")
    var_name <- names(obj)[4]
  } else if (is.numeric(var_name)) {
    var_name <- names(obj)[var_name]
  }

  is_multi_var <- length(var_name) > 1
  spline_results <- mpspline(obj = obj, var_name = var_name, lam = lam,
                             d = d, vlow = vlow, vhigh = vhigh)

  if (is_multi_var) {
    formatted_output <- list()
    for (v_name in var_name) {
      v_results <- spline_results[[v_name]]
      if (is.null(v_results)) {
        v_results <- list()
      }
      formatted_output[[v_name]] <- formatter_fn(v_results, v_name, d, obj, TRUE)
    }
    if (!is.null(combiner_fn)) {
      return(combiner_fn(formatted_output, var_name))
    }
    return(formatted_output)
  }
  formatter_fn(spline_results, var_name, d, obj, FALSE)
}

#' Format spline results as matrices
#' @keywords internal
#' @noRd
.format_compact <- function(splined_results, var_name, d, input_obj, include_variable_col) {
  splined_results <- splined_results[!is.na(splined_results)]
  if (length(splined_results) == 0) {
    return(list('est_icm' = matrix(nrow = 0, ncol = 0),
                'est_1cm' = matrix(nrow = 0, ncol = 0),
                'est_dcm' = matrix(nrow = 0, ncol = 0),
                'tmse' = matrix(nrow = 0, ncol = 0)))
  }

  mh <- max(sapply(splined_results, function(i) length(i[[2]])), na.rm = TRUE)

  sids <- sapply(splined_results, `[[`, 1)

  compact_list <- list(
    'est_icm' = t(sapply(splined_results, function(i) {
      x <- rep(NA, mh)
      x[seq_along(i[[2]])] <- i[[2]]
      x
    }, USE.NAMES = FALSE)),
    'est_1cm' = t(sapply(splined_results, function(i) { i[[3]] })),
    'est_dcm' = {
      x <- t(sapply(splined_results, function(i) { i[[4]] }))
      colnames(x) <- .format_depth_interval_names(d[-length(d)], d[-1])
      depth <- tapply(input_obj[[3]], input_obj[[1]], max)[as.character(sids)]
      cbind(x, depth)
    },
    'tmse' = t(sapply(splined_results, function(i) { i[[5]] }))
  )
  for (name in names(compact_list)) {
    rownames(compact_list[[name]]) <- sids
  }
  return(compact_list)
}

#' Combine multi-variable tidy results
#' @keywords internal
#' @noRd
.combine_tidy_results <- function(formatted_list, var_names) {
  combined <- list()
  for (nm in names(formatted_list[[1]])) {
    combined[[nm]] <- do.call('rbind', lapply(formatted_list, `[[`, nm))
    rownames(combined[[nm]]) <- NULL
  }
  combined
}

#' Format spline results as tidy data frames
#' @keywords internal
#' @noRd
.build_tidy_df <- function(splined_results, col_name, value_extractor, sid_col_name) {
  df <- lapply(splined_results, function(i) {
    d_df <- data.frame("SID" = i[[1]],
                       value_extractor(i),
                       stringsAsFactors = FALSE)
    names(d_df)[1] <- sid_col_name
    d_df
  })
  do.call('rbind', df)
}

.format_tidy <- function(splined_results, var_name, d, input_obj, include_variable_col) {
  splined_results <- splined_results[!is.na(splined_results)]
  if (length(splined_results) == 0) {
    return(list('est_icm' = data.frame(stringsAsFactors = FALSE),
                'est_1cm' = data.frame(stringsAsFactors = FALSE),
                'est_dcm' = data.frame(stringsAsFactors = FALSE),
                'tmse' = data.frame(stringsAsFactors = FALSE)))
  }

  sid_col_name <- names(splined_results[[1]][1])[1]

  # Build ICM (input depth ranges)
  icm <- .build_tidy_df(splined_results, 'DEPTH',
    function(i) {
      list("DEPTH" = names(i$est_icm),
           "SPLINED_VALUE" = as.vector(i$est_icm))
    }, sid_col_name)
  icm$UD <- as.numeric(substr(icm$DEPTH, 1, 3))
  icm$LD <- as.numeric(substr(icm$DEPTH, 5, 7))
  icm <- icm[, c(names(icm[1])[1], 'UD', 'LD', 'SPLINED_VALUE')]

  # Build NCM (1cm predictions)
  max_d <- max(d, na.rm = TRUE)
  ncm <- .build_tidy_df(splined_results, NULL,
    function(i) {
      list("UD" = seq_len(max_d) - 1L,
           "LD" = seq_len(max_d),
           "SPLINED_VALUE" = as.vector(i$est_1cm))
    }, sid_col_name)
  ncm <- ncm[!is.na(ncm$SPLINED_VALUE), ]

  # Build DCM (output depth ranges)
  dcm <- .build_tidy_df(splined_results, 'DEPTH',
    function(i) {
      list("DEPTH" = names(i$est_dcm),
           "SPLINED_VALUE" = as.vector(i$est_dcm))
    }, sid_col_name)
  dcm$UD <- as.numeric(substr(dcm$DEPTH, 1, 3))
  dcm$LD <- as.numeric(substr(dcm$DEPTH, 5, 7))
  dcm <- dcm[, c(names(dcm)[1], 'UD', 'LD', 'SPLINED_VALUE')]
  dcm <- dcm[!is.na(dcm$SPLINED_VALUE), ]

  # Build TMSE (error estimates)
  tmse <- .build_tidy_df(splined_results, NULL,
    function(i) {
      list("ERROR_TYPE" = names(i$est_err),
           "ERROR_VALUE" = as.vector(i$est_err))
    }, sid_col_name)

  output <- list(
    'est_icm' = icm,
    'est_1cm' = ncm,
    'est_dcm' = dcm,
    'tmse' = tmse
  )

  if (include_variable_col) {
    for (name in names(output)) {
      if (nrow(output[[name]]) > 0) {
        output[[name]]$VARIABLE <- var_name
      }
    }
  }

  output
}

#' Spline discrete soils data - multiple sites
#'
#' These functions implement the mass-preserving spline method of Bishop \emph{et
#' al} (1999) (\doi{10.1016/S0016-7061(99)00003-8}) for interpolating between
#' measured soil attributes down a soil profile, across multiple sites' worth of
#' data. \code{mpspline_compact()} returns results as matrices while
#' \code{mpspline_tidy()} returns results as data frames.
#'
#' @param obj data.frame or matrix. Column 1 must contain site identifiers.
#'   Columns 2 and 3 must contain upper and lower sample depths, respectively,
#'   measured in centimeters. Subsequent columns will contain measured values
#'   for those depths.
#' @param var_name character or integer vector denoting the column(s) in
#'   \code{obj} in which target data is stored. If not supplied, the fourth
#'   column of the input object is assumed to contain the target data.
#' @param lam number; smoothing parameter for spline. Defaults to 0.1.
#' @param d sequential integer vector; denotes the output depth ranges in cm.
#'   Defaults to \code{c(0, 5, 15, 30, 60, 100, 200)} after the GlobalSoilMap
#'   specification, giving output predictions over intervals 0-5cm, 5-15cm,
#'   etc.
#' @param vlow numeric; constrains the minimum predicted value to a realistic
#'   number. Defaults to 0.
#' @param vhigh numeric; constrains the maximum predicted value to a realistic
#'   number. Defaults to 1000.
#'
#' @return \code{mpspline_compact()} returns a four-item list containing matrices:
#'   predicted values over input depth ranges, output depth ranges, 1cm predictions,
#'   and RMSE values. Site identifiers are in rownames.
#'
#'   \code{mpspline_tidy()} returns a four-item list of data frames with the
#'   same predictions but in tidy format, with an added VARIABLE column when
#'   processing multiple variables.
#' @examples
#' dat <- data.frame("SID" = c( 1,  1,  1,  1,   2,   2,   2,   2),
#'                    "UD" = c( 0, 20, 40, 60,   0,  15,  45,  80),
#'                    "LD" = c(10, 30, 50, 70,   5,  30,  60, 100),
#'                   "VAL" = c( 6,  4,  3, 10, 0.1, 0.9, 2.5,   6),
#'                    stringsAsFactors = FALSE)
#' # single variable
#' result <- mpspline_compact(obj = dat, var_name = 'VAL')
#'
#' # multiple variables
#' dat_multi <- data.frame( "SID" = c( 1,  1,  1,  1,   2,   2,   2,   2),
#'                           "UD" = c( 0, 20, 40, 60,   0,  15,  45,  80),
#'                           "LD" = c(10, 30, 50, 70,   5,  30,  60, 100),
#'                         "VAL1" = c( 6,  4,  3, 10, 0.1, 0.9, 2.5,   6),
#'                         "VAL2" = c( 5,  3,  2,  9, 0.2, 1.0, 2.0,   5),
#'                         stringsAsFactors = FALSE)
#' result_multi <- mpspline_compact(obj = dat_multi, var_name = c('VAL1', 'VAL2'))
#' @rdname mpspline-wrappers
#' @export
#'
mpspline_compact <- function(obj = NULL, var_name = NULL, lam = 0.1,
                             d = c(0, 5, 15, 30, 60, 100, 200),
                             vlow = 0, vhigh = 1000) {
  .mpspline_wrapper(obj, var_name, lam, d, vlow, vhigh, .format_compact)
}

#' @rdname mpspline-wrappers
#' @examples
#' \dontshow{
#' # Reuse example data from mpspline_compact
#' dat <- data.frame("SID" = c( 1,  1,  1,  1,   2,   2,   2,   2),
#'                    "UD" = c( 0, 20, 40, 60,   0,  15,  45,  80),
#'                    "LD" = c(10, 30, 50, 70,   5,  30,  60, 100),
#'                   "VAL" = c( 6,  4,  3, 10, 0.1, 0.9, 2.5,   6),
#'                    stringsAsFactors = FALSE)
#' }
#' # Single variable with tidy output
#' result <- mpspline_tidy(obj = dat, var_name = 'VAL')
#'
#' # Multiple variables
#' dat_multi <- data.frame( "SID" = c( 1,  1,  1,  1,   2,   2,   2,   2),
#'                           "UD" = c( 0, 20, 40, 60,   0,  15,  45,  80),
#'                           "LD" = c(10, 30, 50, 70,   5,  30,  60, 100),
#'                         "VAL1" = c( 6,  4,  3, 10, 0.1, 0.9, 2.5,   6),
#'                         "VAL2" = c( 5,  3,  2,  9, 0.2, 1.0, 2.0,   5),
#'                         stringsAsFactors = FALSE)
#' result_multi <- mpspline_tidy(obj = dat_multi, var_name = c('VAL1', 'VAL2'))
#' subset(result_multi$est_dcm, VARIABLE == 'VAL1')
#' @export
#'
mpspline_tidy <- function(obj = NULL, var_name = NULL, lam = 0.1,
                          d = c(0, 5, 15, 30, 60, 100, 200),
                          vlow = 0, vhigh = 1000) {
  .mpspline_wrapper(obj, var_name, lam, d, vlow, vhigh, .format_tidy, .combine_tidy_results)
}

