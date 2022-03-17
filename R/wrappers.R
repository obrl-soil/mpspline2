#' Spline discrete soils data - multiple sites, compact output
#'
#' This function implements the mass-preserving spline method of Bishop \emph{et
#' al} (1999) (\doi{10.1016/S0016-7061(99)00003-8}) for interpolating between
#' measured soil attributes down a soil profile, across multiple sites' worth of
#' data. It returns a more compact output object than
#' \code{\link[mpspline2:mpspline]{mpspline()}}.
#' @param obj data.frame or matrix. Column 1 must contain site identifiers.
#'   Columns 2 and 3 must contain upper and lower sample depths, respectively.
#'   Subsequent columns will contain measured values for those depths.
#' @param var_name length-1 character or length-1 integer denoting the column in
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
#' @return A four-item list containing a matrix of predicted values over the
#'   input depth ranges, a matrix of predicted values over the output depth
#'   ranges, a matrix of 1cm predictions, and a matrix of RMSE and IQR-scaled
#'   RMSE values. Site identifiers are in rownames attributes.
#' @examples
#' dat <- data.frame("SID" = c( 1,  1,  1,  1,   2,   2,   2,   2),
#'                    "UD" = c( 0, 20, 40, 60,   0,  15,  45,  80),
#'                    "LD" = c(10, 30, 50, 70,   5,  30,  60, 100),
#'                   "VAL" = c( 6,  4,  3, 10, 0.1, 0.9, 2.5,   6),
#'                    stringsAsFactors = FALSE)
#' mpspline_compact(obj = dat, var_name = 'VAL')
#' @export
#'
mpspline_compact <- function(obj = NULL, var_name = NULL, lam = 0.1,
                             d = c(0, 5, 15, 30, 60, 100, 200),
                             vlow = 0, vhigh = 1000) {

  obj <- mpspline_conv(obj)

  if(is.null(var_name)) {
    message("Parameter var_name not supplied, assuming target data is in column 4.")
    var_name <- names(obj)[4]
  }

  sites <- split(obj, as.factor(obj[[1]]))

  splined <- lapply(sites, mpspline_one, var_name = var_name, lam = lam,
                    d = d, vlow = vlow, vhigh = vhigh)
  keep <- which(!is.na(splined))
  splined <- splined[keep]

  mh <- max(sapply(splined, function(i) length(i[[2]])), na.rm = TRUE)

  list(# don't need an ID - see dimnames(x) on each of the following:
       'est_icm' = t(sapply(splined, function(i) {
         x <- rep(NA, mh)
         x[1:length(i[[2]])] <- i[[2]]
         x
       }, USE.NAMES = FALSE)),
       'est_1cm' = t(sapply(splined, function(i) { i[[3]] })),
       'est_dcm' = {
         x <- t(sapply(splined, function(i) { i[[4]] }))
         names(x) <- mapply(function(u, l) {
           paste0(sprintf('%03d', u), '_', sprintf('%03d', l), '_cm')
         }, u = d[1:(length(d) - 1)], l = d[2:length(d)])
         depth <- sapply(sites[keep], function(i) max(i[[3]]))
         cbind(x, depth)
       },
       'tmse' = t(sapply(splined, function(i) { i[[5]] }))
  )

}

#' Spline discrete soils data - multiple sites, tidy output
#'
#' This function implements the mass-preserving spline method of Bishop \emph{et
#' al} (1999) (\doi{10.1016/S0016-7061(99)00003-8}) for interpolating between
#' measured soil attributes down a soil profile, across multiple sites' worth of
#' data. It returns an output object with tidy data formatting.
#' @param obj data.frame or matrix. Column 1 must contain site identifiers.
#'   Columns 2 and 3 must contain upper and lower sample depths, respectively,
#'   and be measured in centimeters. Subsequent columns will contain measured
#'   values for those depths.
#' @param var_name length-1 character or length-1 integer denoting the column in
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
#' @return A four-item list containing data frames of predicted values over the
#'   input depth ranges, the output depth ranges, 1cm-increment predictions, and
#'   RMSE and IQR-scaled RMSE values.
#' @examples
#' dat <- data.frame("SID" = c( 1,  1,  1,  1,   2,   2,   2,   2),
#'                    "UD" = c( 0, 20, 40, 60,   0,  15,  45,  80),
#'                    "LD" = c(10, 30, 50, 70,   5,  30,  60, 100),
#'                   "VAL" = c( 6,  4,  3, 10, 0.1, 0.9, 2.5,   6),
#'                    stringsAsFactors = FALSE)
#' mpspline_tidy(obj = dat, var_name = 'VAL')
#' @export
#'
mpspline_tidy <- function(obj = NULL, var_name = NULL, lam = 0.1,
                          d = c(0, 5, 15, 30, 60, 100, 200),
                          vlow = 0, vhigh = 1000) {

  obj <- mpspline_conv(obj)

  if(is.null(var_name)) {
    message("Parameter var_name not supplied, assuming target data is in column 4.")
    var_name <- names(obj)[4]
  }

  sites <- split(obj, as.factor(obj[[1]]))

  splined <- lapply(sites, mpspline_one, var_name = var_name, lam = lam,
                    d = d, vlow = vlow, vhigh = vhigh)
  keep <- which(!is.na(splined))
  splined <- splined[keep]

  icm <- lapply(splined, function(i) {
    d <- data.frame("SID" = i[[1]],
                    "DEPTH" = names(i$est_icm),
                    "SPLINED_VALUE" = as.vector(i$est_icm),
                    stringsAsFactors = FALSE)
    names(d)[1] <- names(i[1])[1]
    d
  })
  icm <- do.call('rbind', icm)
  icm$UD <- as.numeric(substr(icm$DEPTH, 1, 3))
  icm$LD <- as.numeric(substr(icm$DEPTH, 5, 7))
  icm <- icm[, c(names(icm[1])[1], 'UD', 'LD', 'SPLINED_VALUE')]

  ncm <- lapply(splined, function(i) {
    d <- data.frame("SID" = i[[1]],
                    "UD" = as.numeric(seq(200)),
                    "LD" = seq(200) + 1,
                    "SPLINED_VALUE" = as.vector(i$est_1cm),
                    stringsAsFactors = FALSE)
    names(d)[1] <- names(i[1])[1]
    d
  })
  ncm <- do.call('rbind', ncm)
  ncm <- ncm[!is.na(ncm$SPLINED_VALUE), ]

  dcm <- lapply(splined, function(i) {
    d <- data.frame("SID" = i[[1]],
                    "DEPTH" = names(i$est_dcm),
                    "SPLINED_VALUE" = as.vector(i$est_dcm),
                    stringsAsFactors = FALSE)
    names(d)[1] <- names(i[1])[1]
    d
  })
  dcm <- do.call('rbind', dcm)
  dcm$UD <- as.numeric(substr(dcm$DEPTH, 1, 3))
  dcm$LD <- as.numeric(substr(dcm$DEPTH, 5, 7))
  dcm <- dcm[, c(names(dcm)[1], 'UD', 'LD', 'SPLINED_VALUE')]
  dcm <- dcm[!is.na(dcm$SPLINED_VALUE), ]

  tmse <- lapply(splined, function(i) {
    d <- data.frame("SID" = i[[1]],
                    "ERROR_TYPE" = names(i$est_err),
                    "ERROR_VALUE" = as.vector(i$est_err),
                    stringsAsFactors = FALSE)
    names(d)[1] <- names(i[1])[1]
    d
  })
  tmse <- do.call('rbind', tmse)

  list(
    'est_icm' = icm,
    'est_1cm' = ncm,
    'est_dcm' = dcm,
    'tmse'    = tmse
  )

}

