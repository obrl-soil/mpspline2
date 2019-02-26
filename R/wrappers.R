#' Spline discrete soils data - multiple sites, compact output
#'
#' This function implements the mass-preserving spline method of (Bishop et al
#' (1999))[http://dx.doi.org/10.1016/S0016-7061(99)00003-8] for interpolating
#' between measured soil attributes down a soil profile, across multiple sites'
#' worth of data. It returns a compact output object similar to
#'  \code{\link[GSIF:mpspline]{GSIF::mpspline()}}
#' @param obj Object of class SoilProfileCollection (see package 'aqp') or data
#'   frame or matrix. For data frames and matrices, column 1 must contain site
#'   identifiers. Columns 2 and 3 must contain upper and lower sample depths,
#'   respectively. Subsequent columns will contain measured values for those
#'   depths. For SoilProfileCollections, the `@horizons` slot must be similarly
#'   arranged, and the `@idcol` and `@depthcol` slots must be correctly defined.
#' @param var_name length-1 character or length-1 integer denoting the column in
#'   `obj` in which target data is stored. If not supplied, the fourth column of
#'   the input object is assumed to contain the target data.
#' @param lam number; smoothing parameter for spline. Defaults to 0.1.
#' @param d sequential integer vector; denotes the output depth ranges in cm.
#'   Defaults to `c(0, 5, 15, 30, 60, 100, 200)` after the globalsoilmap.net
#'   specification, giving output predictions over intervals 0-5cm, 5-15cm,
#'   etc.
#' @param vlow numeric; constrains the minimum predicted value to a realistic
#'   number. Defaults to 0.
#' @param vhigh numeric; constrains the maximum predicted value to a realistic
#'   number. Defaults to 1000.
#' @return A four-item list containing a matrix of
#'   predicted values over the input depth ranges, a matrix of predicted
#'   values over the output depth ranges, a matrix of 1cm predictions, and a
#'   matrix of RMSE and IQR-scaled RMSE values. Site identifiers are in rownames
#'   attributes.
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

#' Spline discrete soils data - multiple sites, SoilProfileCollection output
#'
#' This function implements the mass-preserving spline method of (Bishop et al
#' (1999))[http://dx.doi.org/10.1016/S0016-7061(99)00003-8] for interpolating
#' between measured soil attributes down a soil profile, across multiple sites'
#' worth of data. It returns a SoilProfileCollection object containing all
#' outputs.
#' @param obj Object of class SoilProfileCollection (see package 'aqp') or data
#'   frame or matrix. For data frames and matrices, column 1 must contain site
#'   identifiers. Columns 2 and 3 must contain upper and lower sample depths,
#'   respectively. Subsequent columns will contain measured values for those
#'   depths. For SoilProfileCollections, the `@horizons` slot must be similarly
#'   arranged, and the `@idcol` and `@depthcol` slots must be correctly defined.
#' @param var_name length-1 character or length-1 integer denoting the column in
#'   `obj` in which target data is stored. If not supplied, the fourth column of
#'   the input object is assumed to contain the target data.
#' @param lam number; smoothing parameter for spline. Defaults to 0.1.
#' @param d sequential integer vector; denotes the output depth ranges in cm.
#'   Defaults to `c(0, 5, 15, 30, 60, 100, 200)` after the globalsoilmap.net
#'   specification, giving output predictions over intervals 0-5cm, 5-15cm,
#'   etc.
#' @param vlow numeric; constrains the minimum predicted value to a realistic
#'   number. Defaults to 0.
#' @param vhigh numeric; constrains the maximum predicted value to a realistic
#'   number. Defaults to 1000.
#' @return A SoilProfileCollection object containing spline predictions over
#' the input depth ranges, the output depth ranges, and at 1cm intervals, along
#' with root mean squared error (RMSE) and IQR-scaled RMSE.
#' @examples
#' dat <- data.frame("SID" = c( 1,  1,  1,  1,   2,   2,   2,   2),
#'                    "UD" = c( 0, 20, 40, 60,   0,  15,  45,  80),
#'                    "LD" = c(10, 30, 50, 70,   5,  30,  60, 100),
#'                   "VAL" = c( 6,  4,  3, 10, 0.1, 0.9, 2.5,   6),
#'                    stringsAsFactors = FALSE)
#' mpspline_spc(obj = dat, var_name = 'VAL')
#' @importFrom aqp depths<- horizons<- site
#' @importFrom stats as.formula sd
#' @export
#'
mpspline_spc <- function(obj = NULL, var_name = NULL, lam = 0.1,
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

  # hating myself a little here b/c redundant but I need the input
  # horizons post-cleaning
  orig <- suppressMessages(
    lapply(sites, mpspline_datchk, var_name = var_name)
  )
  orig <- orig[keep]

  sidnm <- names(obj)[1]
  all_df <- mapply(function(orig, spln) {
    df_1cm <- data.frame(spln[[1]], '1cm',
                         seq(spln[[3]]) - 1, seq(spln[[3]]),
                         spln[[3]], spln[[5]][[1]], spln[[5]][[2]])
    names(df_1cm) <- c(sidnm, 'est', 'UD_cm', 'LD_cm', var_name,
                       'RMSE', 'RMSE_IQR')
    df_1cm <- df_1cm[!is.na(df_1cm[[var_name]]), ]

    df_dcm <- data.frame(spln[[1]], 'dcm',
                         d[1:(length(d) - 1)], d[2:length(d)],
                         spln[[4]], spln[[5]][[1]], spln[[5]][[2]],
                         row.names = NULL)
    names(df_dcm) <- names(df_1cm)
    df_dcm <- df_dcm[!is.na(df_dcm[[var_name]]), ]

    df_icm <-  data.frame(spln[[1]], 'icm',
                          orig[[2]], orig[[3]], # see whining above
                          spln[[2]], spln[[5]][[1]], spln[[5]][[2]],
                          row.names = NULL)
    names(df_icm) <- names(df_dcm)
    rbind(df_icm, df_dcm, df_1cm)
  },
  orig = orig, spln = splined, SIMPLIFY = FALSE)
  all_df <- do.call('rbind', all_df)
  rownames(all_df) <- seq(dim(all_df)[1])
  all_df <- cbind(paste0(all_df[[1]], '_', all_df[[2]]),
                  all_df)
  names(all_df)[1] <- paste0(sidnm, "_est")
  nm <- names(all_df)
  fm <- as.formula(sprintf("%s ~ %s + %s", nm[1], nm[4], nm[5]))
  suppressWarnings(depths(all_df) <- fm)
  nso <- horizons(all_df)[, seq(3)]
  nso <- base::unique(nso)
  nso$est <- base::ordered(nso$est,
                           levels = c('icm', '1cm', 'dcm')) # overkill? nahhhh
  rownames(nso) <- NULL
  all_df@site <- nso

  all_df
}
