#' Convert data for splining
#'
#' Generate a consistent input object for splining
#' @param obj Object of class SoilProfileCollection (aqp) or data frame or
#'   matrix. For data frames and matrices, column 1 must contain site
#'   identifiers. Columns 2 and 3 must contain upper and lower sample depths,
#'   respectively. Subsequent columns will contain measured values for those
#'   depths. For SoilProfileCollections, the `@horizons` slot must be similarly
#'   arranged, and the `@idcol` and `@depthcol` slots must be correctly defined.
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
  # check numeric
  if(!typeof(obj) %in% c('double', 'integer')) {
    stop('This is not a numeric matrix')
    }

  # return df with some names set
  obj <- as.data.frame(obj)
  names(obj)[1:3] <- c('SID', 'UD', 'LD')
  obj
}

#' @rdname mpspline_conv
#' @inherit mpspline_conv return
#' @method mpspline_conv data.frame
#'
mpspline_conv.data.frame <- function(obj = NULL) {
  obj # >.>
}

#' @rdname mpspline_conv
#' @inherit mpspline_conv return
#' @method mpspline_conv SoilProfileCollection
#'
mpspline_conv.SoilProfileCollection <- function(obj = NULL) {
  dc <- obj@depthcols
  ic <- obj@idcol
  ac <- names(obj@horizons)[-which(names(obj@horizons) %in% c(dc, ic))]
  data.frame(c(obj@horizons[ic],
               obj@horizons[dc],
               obj@horizons[ac]), stringsAsFactors = FALSE)
}

#' pre-spline data checks
#'
#' Runs a few data quality checks and makes some repairs where possible.
#' @param sites list of appropriately formatted site data. Sites with no data to
#'   spline and sites with overlapping depth ranges return NA.
#' @param var_name target variable
#' @keywords internal
#'
mpspline_datchk <- function(sites = NULL, var_name = NULL) {
  # nb seq_along() allows access to list names (SIDs)
  lapply(seq_along(sites), function(i) {
    s <- sites[[i]] # just for readability

    # Drop sites where target param has no data at all
    if(all(is.na(s[[var_name]]))) {
      message('No data values are present for site ', names(sites)[i], '.')
      return(NA)
    }

    # remove any horizons where val == NA
    s <- s[!is.na(s[[var_name]]), ]

    # replace any missing surface value
    if(is.na(s[[2]][1])) {
      s[[2]][1] <- 0
    }

    # replace any missing max input depth value
    if(is.na(s[[3]][nrow(s)])) {
      # NB more conservative than existing approach (stretches from last ud to
      # either 150 or 200cm), but ud + 10cm is more realistic
      s[[3]][nrow(s)] <- s[[2]][nrow(s)] + 10
    }

    # remove any horizons with -ve depths
    s <- s[!(s[[2]] < 0 | s[[3]] < 0), ]

    # sort by cols 1, 2, 3 asc
    s <- s[order(s[[1]], s[[2]], s[[3]]), ]
    rownames(s) <- NULL

    # Drop sites with overlapping data depth ranges
    if(any(diff(as.vector(rbind(s[[2]], s[[3]]))) < 0)) {
      message("Overlapping depth ranges detected in site ", names(sites)[i], '.')
      return(NA)
    }

    s
    })
}

#' Estimate spline parameters
#'
#' estimate key parameters for building a mass-preserving spline across a single
#' profile
#' @param s data.frame containing a single profile's worth of soil info
#' @param var_name target variable.
#' @param lam number; smoothing parameter for spline. Defaults to 0.1.
#' @return A list of parameters used for spline fitting.
#' @keywords internal
#'
mpspline_est1 <- function(s = NULL, var_name = NULL, lam = NULL) {

  # don't bother when only one layer's worth of data
  if(dim(s)[1] == 1) {
    return(NA) # NB message later
  }

  # find n inputs, boundaries, layer thicknesses & gap thicknesses (if present)
  n <- dim(s)[1]
  nb <- n - 1
  th <- s[[3]] - s[[2]]
  # (ud[2] - ld[1]...ud[nrow] - ld[nrow - 1]):
  gp <- as.vector(rbind(s[[2]][2:length(s[[2]])], s[[3]][1:length(s[[3]])-1]))
  # http://r.789695.n4.nabble.com/odd-even-indices-of-a-vector-td4694337.html:
  gp <- gp[c(TRUE, FALSE)] - gp[c(FALSE, TRUE)]

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
  s_bar <- solve(Z, as.matrix(s[[var_name]]))

  # calculate the fitted value at the knots
  b  <- 6 * R_inv %*% Q %*% s_bar
  b0 <- rbind(0, b) # add a row to top = 0
  b1 <- rbind(b, 0) # add a row to bottom = 0
  gamma <- (b1 - b0) / as.matrix(th * 2)
  alfa <- s_bar - b0 * as.matrix(th) / 2 - gamma * as.matrix(th)^2/3

  # just return the stuff needed for subsequent steps (decompose to vec or nah?)
  list("s_bar" = s_bar, "b0" = b0, "b1" = b1,
       "gamma" = gamma, "alfa" = alfa, "Z" = Z)
}

#' Fit spline parameters
#'
#' Fit spline parameters to data for a single site.
#' @param s data for one site
#' @param p estimated spline parameters for one site
#' @inheritParams mpspline
#' @return list of two vectors: fitted values at 1cm intervals and the average
#'   of same over the requested depth ranges.
#' @keywords internal
#' @importFrom utils tail
#'
mpspline_fit1 <- function(s = NULL, p = NULL, var_name = NULL,
                          d = NULL, vhigh = NULL, vlow = NULL) {

  # single horizon - no splining for you!
  if(all(is.na(p))) {
    ud <- s[[2]]
    ld <- s[[3]]
    not_1cm <- rep(NA_real_, times = max(d))
    not_1cm[seq(ud, ld)] <- s[[var_name]]
    not_dcm <- rep(NA_real_, times = length(d))
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

  nj <- max(s[[3]])
  if (nj > max(d)) { nj <- max(d) } # if profile > max d, ignore the deeper part

  cm_ds <- seq(nj) - 1 # need 0-index to match splinetool.exe *sigh*

  # name each depth with its membership layer, NA for gaps
  names(cm_ds) <- sapply(cm_ds, function(i) {
    x <- which(s[[2]] <= i & s[[3]] > i)
    if(length(x) == 0) { NA_integer_ } else { x }
  })

  est_1cm <- sapply(cm_ds, function(k) {  # for every cm to max(d); 0 = 0 - 1 cm
    # if input data from profile starts below surface, return alfa of first
    # available horizon for all depths above s[[2]][1]:
    if (k < s[[2]][1]) {
      return(alfa[1])
    }
    # NB splinetool.exe does something different, but not better IMO - looks
    # like it can over or underestimate where a 0-x sample is not present.
    # Would be more conservative to refuse to predict above the first sample.

    h <- as.integer(names(cm_ds)[k + 1]) # uuugghhhhhhh

    # if not in a gap, predict using the params from that layer:
    if(!is.na(h)) {
      alfa[h] + b0[h] * (k - s[[2]][h]) + gamma[h] * (k - s[[2]][h])^2
    } else {
      # infill with reference to params from layers either side of the gap:
      bh  <- utils::tail(which(s[[3]] <= k), 1) # prev sample range
      nh  <- which(s[[2]] >= k)[1] # next sample range
      phi <- alfa[nh] - b1[bh] * (s[[2]][nh] - s[[3]][bh])
      phi + b1[bh] * (k - s[[3]][bh])
    }
  }, USE.NAMES = FALSE)
 # est_1cm <- unlist(est_1cm, use.names = FALSE) # sometimes sapply isn't?!
  names(est_1cm) <- seq(nj)

  # constrain 1cm estimates to supplied limits
  est_1cm[which(est_1cm > vhigh)] <- vhigh
  est_1cm[which(est_1cm < vlow)]  <- vlow

  # pad vec to max(d)
  if(nj < max(d)) {
    est_1cm <- rep(est_1cm, length.out = max(d))
    est_1cm[(nj + 1):max(d)] <- NA_real_
  }

  # average est_1cm over the output depths
  # NB doesn't match current mpspline but does match splinetool.exe
  od_mat <- matrix(c(d[1:length(d) - 1],
                     d[2:length(d)]), ncol = 2, byrow = FALSE)
  # will need a tweak if anyone tries to go past 999cm:
  nms <- apply(od_mat, 1, function(i) paste0(sprintf('%03d', i[1]), '_',
                                             sprintf('%03d', i[2]), '_cm'))
  od_mat[, 1] <- od_mat[, 1] + 1 # e.g. 1 - 5, = >= 0cm and < 5cm etc
  est_dcm <- apply(od_mat, 1, FUN = function(i) {
    mean(est_1cm[i[1]:i[2]], na.rm = TRUE)
  })
  names(est_dcm) <- nms
  est_dcm[is.nan(est_dcm)] <- NA_real_

  # return list with 1cm and dcm as vectors
  list("est_1cm" = est_1cm, "est_dcm" = est_dcm)
}

#' calculate TSME
#'
#' Calculates Total Mean Squared Error (TMSE) for a single site
#' @param s site data frame
#' @param p estimated spline params for site
#' @param var_name target variable
#' @param s2 numeric, 5\% of the variance for the parent dataset
#' @return numeric, tmse
#' @keywords internal
#'
mpspline_tmse1 <- function(s = NULL, p = NULL, var_name = NULL, s2 = NULL) {
  # single layer in input
  if(all(is.na(p))) {
    return(NA_real_)
  }

  s_bar <- p[['s_bar']]
  Z     <- p[['Z']]
  n     <- dim(s)[1]

  ssq <- sum((as.matrix(s[[var_name]]) - s_bar)^2)
  g <- solve(Z)
  ei <- eigen(g)$values
  df <- n - sum(ei)
  sig2w <- ssq / df
  ## calculate the Carter and Eagleson estimate of residual variance
  dfc <- n - 2 * sum(ei) + sum(ei^2)
  sig2c <- ssq / dfc
  ## calculate the estimate of the true mean squared error
  ssq / n - 2 * s2 * df / n + s2
}

#' Spline discrete soils data
#'
#' This function implements the mass-preserving spline method of Bishop et al
#' (1999) - http://dx.doi.org/10.1016/S0016-7061(99)00003-8] for interpolating
#' between measured soils data parameters down a profile.
#' @param obj Object of class SoilProfileCollection (aqp) or data frame or
#'   matrix. For data frames and matrices, column 1 must contain site
#'   identifiers. Columns 2 and 3 must contain upper and lower sample depths,
#'   respectively. Subsequent columns will contain measured values for those
#'   depths. For SoilProfileCollections, the `@horizons` slot must be similarly
#'   arranged, and the `@idcol` and `@depthcol` slots must be correctly defined.
#' @param var_name length-1 character denoting the column in `obj` in which
#'   target data is stored.
#' @param lam number; smoothing parameter for spline. Defaults to 0.1.
#' @param d sequential integer vector; denotes the output depth ranges in cm.
#'   Defaults to `c(0, 5, 15, 30, 60, 100, 200)` after the globalsoilmap.net
#'   specification, giving output predictions over intervals 0-5cm, 5-15cm,
#'   etc.
#' @param vlow numeric; constrains the minimum predicted value to a realistic
#'   number. Defaults to 0.
#' @param vhigh numeric; constrains the maximum predicted value to a realistic
#'   number. Defaults to 1000.
#' @return list of five data elements for each site - input data, predicted
#'   values at each cm down the profile, predicted values over `d` intervals,
#'   and TMSE.
#' @examples
#' dat <- data.frame("SID" = c( 1,  1,  1,  1,   2,   2,   2,   2),
#'                    "UD" = c( 0, 20, 40, 60,   0,  15,  45,  80),
#'                    "LD" = c(10, 30, 50, 70,   5,  30,  60, 100),
#'                   "VAL" = c( 6,  4,  3, 10, 0.1, 0.9, 2.5,   6),
#'                    stringsAsFactors = FALSE)
#' mpspline(obj = dat, var_name = 'VAL')
#' @importFrom stats sd
#' @export
#'
mpspline <- function(obj = NULL, var_name = NULL, lam = 0.1,
                     d = c(0,5,15,30,60,100,200),
                     vlow = 0, vhigh = 1000) {

  nice_obj <- mpspline_conv(obj)

  # split input data into a list by site
  sites <- split(nice_obj, as.factor(nice_obj[[1]]))

  # do some checks and tidying up
  sites <- mpspline_datchk(sites, var_name = var_name)
  sites <- sites[!is.na(sites)]

  # find the max number of samples/horizons across all input sites post-clean
  # tbh not sure this is really necessary
  max_dat <- max(table(do.call('rbind', sites)[[1]]))

  # *shrug* it could happen I guess??
  if(max_dat == 1) {
    stop('Supplied sites all have one depth range, please check inputs.')
  }

  # estimate spline parameters for each site # prog bar it later maybe
  message("Estimating spline parameters at ", Sys.time())
  params <- lapply(sites, function(i) {
    mpspline_est1(i, var_name = var_name, lam = lam)
    })

  # fit spline to each site
  message("Fitting spline parameters at ", Sys.time())
  splined <- mapply(function(s, p) {
    mpspline_fit1(s, p, var_name = var_name,
                  d = d, vhigh = vhigh, vlow = vlow) },
    s = sites, p = params, SIMPLIFY = FALSE)

  ### TO DO: TESTS FOR MPSPLINE_FIT1

  # cl_dat is all the data values in obj[[var_name]] post-clean
  cl_dat <- unlist(sapply(sites, function(i) i[[var_name]]), use.names = FALSE)
  s_hat_5 <- (0.05 * stats::sd(cl_dat, na.rm = TRUE))^2
  var_5 <- s_hat_5^2

  tmses <- mapply(function(s, p) {
    mpspline_tmse1(s, p, var_name = var_name, s2 = var_5)
    },
    s = sites, p = params)

  ## TO DO: tests for mpspline_tsme1

  # OUTPUTS - make a condensed option? List of dataframes

  # list-based
  out <- mapply(function(s, e, tmse, lam) {
    list("inputs"  = s,
         "est_1cm" = e[[1]],
         "est_dcm" = e[[2]],
         "tmse"    = tmse,
         "lam"     = lam)
  },
  s = sites, e = splined, tmse = tmses, lam = lam, SIMPLIFY = FALSE)

  # mat_based next
  #out_mat <- list("inputs" = do.call('rbind', sites),
  #                "est_1cm" = )
#
  out

}

