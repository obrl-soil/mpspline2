#' Convert data for splining
#'
#' Generate a consistent input object for splining
#' @param obj data.frame or matrix. Column 1 must contain site identifiers.
#'   Columns 2 and 3 must contain upper and lower sample depths, respectively.
#'   Subsequent columns will contain measured values for those depths.
#' @return data frame, sorted by site ID, upper and lower depth.
#' @keywords internal
#' @rdname mpspline_conv
#'
mpspline_conv <- function(obj = NULL) {
  UseMethod('mpspline_conv')
}

# not that I think this is common but jic
#' @rdname mpspline_conv
#' @inherit mpspline_conv return
#' @method mpspline_conv matrix
#'
mpspline_conv.matrix <- function(obj = NULL) {
  # return df with some names set
  out <- data.frame(obj, stringsAsFactors = FALSE)
  if(typeof(obj) == 'character') {
    out[, c(2:ncol(out))] <- as.numeric(unlist(out[, c(2:ncol(out))]))
  }
  names(out)[1:3] <- c('SID', 'UD', 'LD')
  out
}

#' @rdname mpspline_conv
#' @inherit mpspline_conv return
#' @method mpspline_conv data.frame
#'
mpspline_conv.data.frame <- function(obj = NULL) {
  obj # >.>
}

#' pre-spline data checks
#'
#' Runs a few data quality checks and makes some repairs where possible.
#' @param s data frame, input data for a single soil profile.
#' @param var_name length-1 character or length-1 integer denoting the column in
#'   \code{site} in which target data is stored. If not supplied, the fourth
#'   column of the input object is assumed to contain the target data.
#' @return If data passes checks it is returned unchanged. Sites with no data to
#'   spline and sites with overlapping input depth ranges return NA.
#' @keywords internal
#'
mpspline_datchk <- function(s = NULL, var_name = NULL) {
  sid <- s[[1]][1]
  nrs <- nrow(s)

  # Drop sites where target param has no data at all
  if(all(is.na(s[[var_name]]))) {
    message('No data values are present for site ', sid, '.')
    return(NA)
  }

  # remove any horizons where target value missing
  n_mv <- length(which(is.na(s[[var_name]])))
  if(n_mv > 0) {
    message(n_mv, " depth range(s) with missing analytical data removed from ", sid, '.')
    s <- s[!is.na(s[[var_name]]), ]
    nrs <- nrs - n_mv
  }

  # replace any missing surface value for upper depth
  if(is.na(s[[2]][1])) {
    message("Missing surface upper depth replaced with 0 in ", sid, '.')
    s[[2]][1] <- 0
  }

  # replace any missing max input depth value
  if(is.na(s[[3]][nrs])) {
    # NB more conservative than existing approach (stretches from last ud to
    # either 150 or 200cm), but ud + 10cm is more realistic
    message("Missing deepest lower depth replaced with (deepest upper depth + 10) in ", sid, ".")
    s[[3]][nrs] <- s[[2]][nrs] + 10
  }

  # remove any horizons with -ve depths (e.g. organic horizons)
  n_negd <- length(which((s[[2]] < 0 | s[[3]] < 0)))
  if(n_negd > 0) {
    message(n_negd, " depth range(s) with negative depths removed from ", sid, '.')
    s <- s[!(s[[2]] < 0 | s[[3]] < 0), ]
  }

  # remove any horizons with 0-thickness depths (same UD and LD)
  n_0th <- length(which(s[[2]] == s[[3]]))
  if(n_0th > 0) {
    message(n_0th, " depth range(s) with thickness of 0 removed from ", sid, '.')
    s <- s[!(s[[2]] == s[[3]]), ]
  }

  # remove any horizons where upper depth > lower depth
  bkwrds <- length(which(s[[2]] > s[[3]]))
  if(bkwrds > 0) {
    message(bkwrds, " depth ranges removed from ", sid, " as upper depth exceeded lower depth.")
    s <- s[!(s[[2]] > s[[3]]), ]
  }

  # sort by ud, ld
  s <- s[order(s[[2]], s[[3]]), ]
  rownames(s) <- NULL

  # Drop sites with overlapping data depth ranges (user should fix)
  if(nrs > 1) {
    if(any(s[[2]][2:nrs] < s[[3]][1:(nrs - 1)], na.rm = TRUE)) {
      message("Overlapping depth ranges detected in site ", sid, '.')
      return(NA)
    }
  }
  s
}

#' Estimate spline parameters
#'
#' Estimate key parameters for building a mass-preserving spline across a single
#' soil profile
#' @param s data.frame containing a single profile's worth of soil info
#' @param var_name length-1 character or length-1 integer denoting the column in
#'   \code{site} in which target data is stored. If not supplied, the fourth
#'   column of the input object is assumed to contain the target data.
#' @param lam number; smoothing parameter for spline. Defaults to 0.1.
#' @return A list of parameters used for spline fitting.
#' @keywords internal
#'
mpspline_est1 <- function(s = NULL, var_name = NULL, lam = NULL) {

  n <- dim(s)[1]
  # don't bother when only one layer's worth of data
  if(n == 1) {
    return(NA_real_) # NB message later
  }

  # find boundaries, layer thicknesses & gap thicknesses (if present)
  nb <- n - 1
  th <- s[[3]] - s[[2]]
  # (ud[2] - ld[1]...ud[nrow] - ld[nrow - 1]):
  gp <- s[[2]][2:n] - s[[3]][1:nb]

  # (n-1) x (n-1) matrix R with diag elements R[i, i] = 2(x[i + 1] - x[i - 1])
  # and off-diag elements R[i + 1, i] = R[i, i + 1] = x[i + 1] - x[i]
  # x refers to the boundaries of the n horizons x[0] - x[n] where x[0] is
  # usually the surface.
  # (this way is not usually much faster than previous, just more concise)
  n1_mat <- diag(th[2:n] * 2, nrow = nb, ncol = nb)
  ud <- which(n1_mat > 0) - 1
  ud <- ud[which(ud > 0)]
  ld <- which(n1_mat > 0) + 1
  ld <- ld[which(ld <= nb^2)]
  n1_mat[ud] <- n1_mat[ld] <- th[2:(n - 1)]
  th_mat <- diag(th, ncol = nb, nrow = nb) # samp ranges on diag
  gp_mat <- diag(gp, ncol = nb, nrow = nb) # gap ranges on diag
  R <- n1_mat + 2 * th_mat + 6 * gp_mat

  # (n-1) x n matrix Q with Q[i, i] = -1, Q[i, i + 1] = 1 and
  # Q[i, j] = 0 otherwise
  Q <- diag(-1, ncol = n, nrow = nb)
  ud <- which(Q == -1) - 1
  ud <- ud[which(ud > 0)]
  Q[c(ud, n * nb)] <- 1

  # also need inverse of R
  R_inv <- try(solve(R), TRUE)
  stopifnot(is.matrix(R_inv)) # failboat sometimes

  # create the matrix coefficent Z (part of eq 7 in paper) -
  # [I + 6 * n * lam * Q^t * R^-1 * Q] (I = diag(n))
  pr_mat <- matrix(6 * n * lam, ncol = nb, nrow = n)
  f_dub <- pr_mat * t(Q) %*% R_inv
  Z <- diag(n) + f_dub %*% Q

  # solve Z for the input data values
  s_bar <- solve(Z, s[[var_name]])

  # name s_bar with input depth ranges
  names(s_bar) <- mapply(function(u, l) {
    paste0(sprintf('%03d', u), '_', sprintf('%03d', l), '_cm')
  }, u = as.integer(s[[2]]), l = as.integer(s[[3]]))

  # calculate the fitted value at the knots (middle of each input range)
  b     <- as.vector(6 * R_inv %*% Q %*% s_bar)
  b0    <- c(0, b)
  b1    <- c(b, 0)
  gamma <- (b1 - b0) / (th * 2)
  alfa  <- s_bar - b0 * th / 2 - gamma * th^2/3 # nb s_bar names inherit

  # just return the stuff needed for subsequent steps
  list("s_bar" = s_bar, "b0" = b0, "b1" = b1,
       "gamma" = gamma, "alfa" = alfa, "Z" = Z)
}

#' Fit spline parameters
#'
#' Fit spline parameters to data for a single site.
#' @param s data.frame; data for one site
#' @param p list; estimated spline parameters for one site from
#'   \code{\link[mpspline2:mpspline_est1]{mpspline_est1}}
#' @inheritParams mpspline
#' @return list of two vectors: fitted values at 1cm intervals and the average
#'   of same over the requested depth ranges.
#' @keywords internal
#'
mpspline_fit1 <- function(s = NULL, p = NULL, var_name = NULL,
                          d = NULL, vhigh = NULL, vlow = NULL) {
  md    <- max(d)
  # single horizon - no splining for you!
  if(all(is.na(p))) {
    ud <- s[[2]]
    ld <- s[[3]]
    not_1cm <- rep(NA_real_, times = md)
    not_1cm[ud:ld] <- s[[var_name]]
    not_dcm <- rep(NA_real_, times = length(d) - 1)
    not_dcm[which(d >= ud & d <= ld)] <- s[[var_name]]
    # constrain input data to supplied limits
    not_1cm[which(not_1cm > vhigh)] <- vhigh
    not_1cm[which(not_1cm < vlow)]  <- vlow
    not_dcm[which(not_dcm > vhigh)] <- vhigh
    not_dcm[which(not_dcm < vlow)]  <- vlow
    message(paste0("Site ", s[[1]],
                   " only has data for one depth range and was not splined."))
    return(list("est_1cm" = not_1cm, "est_dcm" = not_dcm))
  }

  b0    <- p[['b0']]
  b1    <- p[['b1']]
  gamma <- p[['gamma']]
  alfa  <- p[['alfa']]
  names(alfa) <- NULL

  nj <- max(s[[3]])
  if (nj > md) { nj <- md } # if profile > max d, ignore the deeper part

  cm_ds <- (1:nj) - 1 # need 0-index to match splinetool.exe *sigh*

  # tag each depth with its membership layer number, NA for gaps
  cm_h <- sapply(cm_ds, function(i) {
    x <- which(s[[2]] <= i & s[[3]] > i)
    if(length(x) == 0) { NA_integer_ } else { x }
    })

  est_1cm <- sapply(cm_ds, function(k) {  # for every cm to max(d); 0 = 0 - 1 cm
    # if input data from profile starts below surface, return NA until s[[2]][1]
    # is reached.
    if (k < s[[2]][1]) {
      return(NA_real_)
    }
    # Re: above, splinetool.exe continues to extrapolate outside the input depth
    # ranges which can lead to substantial over-/underpredictions in the
    # surface. GSIF::mpspline appears to attempt to compensate by repeating the
    # prediction for the shallowest 'known' 1cm slice back up to the surface.
    # This is more conservative but still carries out *extrapolation* with an
    # *interpolation* algorithm. I argue that this is not appropriate. If a
    # surface value is missing, either go get it measured or explicitly insert
    # an expert-knowledge-based or PTF'd surface value into the data before
    # splining.

    h <- cm_h[k + 1] # uuugghhhhhhh

    # if not in a gap, predict using the params from that layer:
    if(!is.na(h)) {
      alfa[h] + b0[h] * (k - s[[2]][h]) + gamma[h] * (k - s[[2]][h])^2
    } else {
      # infill with reference to p data from layers either side of the gap:
      bh  <- which(s[[3]] <= k)[length(which(s[[3]] <= k))] # prev sample range
      nh  <- which(s[[2]] >= k)[1] # next sample range
      phi <- alfa[nh] - b1[bh] * (s[[2]][nh] - s[[3]][bh])
      phi + b1[bh] * (k - s[[3]][bh])
    }
  }, USE.NAMES = FALSE)
  #names(est_1cm) <- seq(nj)

  # constrain 1cm estimates to supplied limits
  est_1cm[which(est_1cm > vhigh)] <- vhigh
  est_1cm[which(est_1cm < vlow)]  <- vlow

  # pad vec to max(d)
  if(nj < md) {
    est_1cm <- rep(est_1cm, length.out = md)
    est_1cm[(nj + 1):md] <- NA_real_
  }

  # average est_1cm over the output depths
  # NB doesn't match current mpspline but does match splinetool.exe
  est_dcm <- mapply(function(u, l) { mean(est_1cm[u:l], na.rm = TRUE) },
                    u = d[1:length(d) - 1] + 1, l = d[2:length(d)])
  est_dcm[is.nan(est_dcm)] <- NA_real_

  # return list with 1cm and dcm as vectors
  list("est_1cm" = est_1cm, "est_dcm" = est_dcm)
}

#' calculate RMSE
#'
#' Calculates Root Mean Squared Error (RMSE) for estimates on a single site
#' @param s data.frame; data for one site
#' @param p list; estimated spline parameters for one site from
#'   \code{\link[mpspline2:mpspline_est1]{mpspline_est1}}
#' @param var_name length-1 character or length-1 integer denoting the column in
#'   \code{site} in which target data is stored. If not supplied, the fourth
#'   column of the input object is assumed to contain the target data.
#' @return length-2 named numeric - RMSE and RMSE scaled against input data's
#'   interquartile range.
#' @keywords internal
#' @note Useful for comparing the results of varying parameter \code{lam}.
#' @importFrom stats IQR
#'
mpspline_rmse1 <- function(s = NULL, p = NULL, var_name = NULL) {
  # single layer in input
  if(all(is.na(p))) {
    return(c('RMSE' = NA_real_, 'RMSE_IQR' = NA_real_))
  }

  rmse <- sqrt(sum((s[[var_name]] - p[['s_bar']])^2) / dim(s)[1] )
  rmse_iqr <- rmse / stats::IQR(s[[var_name]])
  c('RMSE' = rmse, 'RMSE_IQR' = rmse_iqr)
}

#' Spline discrete soils data - single site
#'
#' This function implements the mass-preserving spline method of Bishop \emph{et
#' al} (1999) (\doi{10.1016/S0016-7061(99)00003-8}) for interpolating between
#' measured soil attributes down a single soil profile.
#' @param site data frame containing data for a single soil profile.
#'   Column 1 must contain site identifiers. Columns 2 and 3 must contain upper
#'   and lower sample depths, respectively, measured in centimeters. Subsequent
#'   columns will contain measured values for those depths.
#' @param var_name length-1 character or length-1 integer denoting the column in
#'   \code{site} in which target data is stored. If not supplied, the fourth
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
#' @return A list with the following elements: Site ID, vector of predicted
#'   values over input intervals, vector of predicted values for each cm down
#'   the profile to \code{max(d)}, vector of predicted values over \code{d}
#'   (output) intervals, and root mean squared error.
#' @examples
#' dat <- data.frame("SID" = c( 1,  1,  1,  1),
#'                    "UD" = c( 0, 20, 40, 60),
#'                    "LD" = c(10, 30, 50, 70),
#'                   "VAL" = c( 6,  4,  3, 10),
#'                    stringsAsFactors = FALSE)
#' mpspline_one(site = dat, var_name = 'VAL')
#' @importFrom stats sd
#' @export
#'
mpspline_one <- function(site = NULL, var_name = NULL, lam = 0.1,
                         d = c(0, 5, 15, 30, 60, 100, 200),
                         vlow = 0, vhigh = 1000) {

  if(!inherits(site, 'data.frame')) {
    site <- mpspline_conv(site)
  }

  # assume if varname missing (or conv from mat etc)
  if(is.null(var_name)) {
    message("Parameter var_name not supplied, assuming target data is in column 4.")
    var_name <- names(site)[4]
  }

  # do some checks and tidying up
  s <- mpspline_datchk(site, var_name = var_name)
  if(!inherits(s, 'data.frame')) { return(NA) }

  # estimate spline parameters for each site and fit
  p <- mpspline_est1(s, var_name = var_name, lam = lam)
  e <- mpspline_fit1(s, p, var_name = var_name,
                     d = d, vhigh = vhigh, vlow = vlow)

  # estimate error
  t <- mpspline_rmse1(s, p, var_name = var_name)

  if(is.factor(s[[1]])) {
    s[[1]] <- as.character(s[[1]])
  }

  if(all(is.na(p))) { icm <- s[[var_name]][1]
    names(icm) <- paste0(sprintf('%03d', as.integer(s[[2]][1])), '_',
                         sprintf('%03d', as.integer(s[[3]][1])), '_cm')
    }

  # arrange outputs
  splined <-
    list("ID" = s[[1]][1],
         # below matches splinetool.exe - should this return alfa tho??
         "est_icm" = if(all(is.na(p))) { icm } else { p[['s_bar']] },
         "est_1cm" = e[[1]],
         "est_dcm" = e[[2]],
         "est_err" = t
         )
  # preserve input SID column name
  names(splined)[1] <- names(site)[1]
  # name est_dcm like est_ins:
  names(splined[['est_dcm']]) <- mapply(function(u, l) {
    paste0(sprintf('%03d', u), '_', sprintf('%03d', l), '_cm')
    }, u = d[1:(length(d) - 1)], l = d[2:length(d)])
  splined
}

#' Spline discrete soils data - multiple sites
#'
#' This function implements the mass-preserving spline method of Bishop \emph{et
#' al} (1999) (\doi{10.1016/S0016-7061(99)00003-8}) for interpolating between
#' measured soil attributes down a soil profile, across multiple sites' worth of
#' data.
#' @param obj data.frame or matrix. Column 1 must contain site
#'   identifiers. Columns 2 and 3 must contain upper and lower sample depths,
#'   respectively. Subsequent columns will contain measured values for those
#'   depths.
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
#' @return A nested list of data for each input site. List elements are: Site
#'   ID, vector of predicted values over input intervals, vector of predicted
#'   values for each cm down the profile to \code{max(d)}, vector of predicted
#'   values over \code{d} (output) intervals, and root mean squared error.
#' @examples
#' dat <- data.frame("SID" = c( 1,  1,  1,  1,   2,   2,   2,   2),
#'                    "UD" = c( 0, 20, 40, 60,   0,  15,  45,  80),
#'                    "LD" = c(10, 30, 50, 70,   5,  30,  60, 100),
#'                   "VAL" = c( 6,  4,  3, 10, 0.1, 0.9, 2.5,   6),
#'                    stringsAsFactors = FALSE)
#' m1 <- mpspline(obj = dat, var_name = 'VAL')
#' @export
#'
mpspline <- function(obj = NULL, var_name = NULL, lam = 0.1,
                     d = c(0, 5, 15, 30, 60, 100, 200),
                     vlow = 0, vhigh = 1000) {

  obj <- mpspline_conv(obj)

  # assume if varname missing (or conv from mat etc)
  if(is.null(var_name)) {
    message("Parameter var_name not supplied, assuming target data is in column 4.")
    var_name <- names(obj)[4]
  }

  sites <- split(obj, as.factor(obj[[1]]))

  splined <- lapply(sites, mpspline_one, var_name = var_name, lam = lam,
                    d = d, vlow = vlow, vhigh = vhigh)

  splined[which(!is.na(splined))]
}

