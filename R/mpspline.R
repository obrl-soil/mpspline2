#' convert data for splining
#'
#' generate a consistent input object for splining
#'
#' @param obj Object of class SoilProfileCollection (aqp) or data frame or
#'   matrix. For data frames and matrices, column 1 must contain site
#'   identifiers. Columns 2 and 3 must contain upper and lower sample depths,
#'   respectively. Subsequent columns will contain measured values for those
#'   depths. For SoilProfileCollections, the `@horizons` slot must be similarly
#'   arranged, and the `@idcol` and `@depthcol` slots must be correctly defined.
#' @return data frame, sorted by site ID, upper and lower depth.
#' @rdname mpspline
#'
mpspline_conv <- function(obj = NULL) {
  UseMethod('mpsline_conv')
}

# not that I think this is common but jic
#' @rdname mpsline_conv
#' @inherit mpsline_conv return
#' @method mpsline_conv matrix
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

#' @rdname mpsline_conv
#' @inherit mpsline_conv return
#' @method mpsline_conv data.frame
#'
mpspline_conv.data.frame <- function(obj = NULL) {
  obj # >.>
}

#' @rdname mpsline_conv
#' @inherit mpsline_conv return
#' @method mpsline_conv SoilProfileCollection
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
#' @param sites list of appropriately formatted site data
#' @keywords internal
#'
mpspline_datchk <- function(sites = NULL) {
  # nb seq_along() allows access to list names (SIDs)
  lapply(seq_along(sites), function(i) {
    s <- sites[[i]] # just for readability

    # make sure at least some data is present
    if(all(is.na(s[[4]]))) {
      stop('No analytical values are present for site ', names(sites)[i], '.')
    }

    # replace any missing surface values
    if(is.na(s[[2]][1])) {
      s[[2]][1] <- 0
    }

    # replace any missing max input depth values
    if(is.na(s[[3]][nrow(s)])) {
      # NB more conservative than existing approach (stretches from last ud to
      # either 150 or 200cm), but I would argue ud + 10 is more realistic
      s[[3]][nrow(s)] <- s[[2]][nrow(s)] + 10
    }

    # remove horizons with -ve depths
    s <- s[!(s[[2]] < 0 | s[[3]] < 0), ]

    # sort by cols 1, 2, 3 asc
    s <- s[order(s[[1]], s[[2]], s[[3]]), ]
    rownames(s) <- NULL

    # Fail if overlapping data depth ranges
    if(any(diff(as.vector(rbind(s[[2]], s[[3]]))) < 0)) {
      stop("Overlapping horizons detected in site ", names(sites)[i], '.')
    }

    s
    })
}

#' Estimate spline parameters
#'
#' estimate key parameters for building a mass-preserving spline across a site's
#' worth of data
#' @param site data.frame containing a site's worth of soil info
#' @param lam smoothing thingo
#' @return A list of parameters used for spline fitting
#' @keywords internal
#'
mpspline_est1 <- function(s = NULL, lam = NULL) {

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
  s_bar <- solve(Z, as.matrix(s[[4]])) # tsme comes from comparing these with input

  # calculate the fitted value at the knots
  b  <- 6 * R_inv %*% Q %*% s_bar
  b0 <- rbind(0, b) # add a row to top = 0
  b1 <- rbind(b, 0) # add a row to bottom = 0
  gamma <- (b1 - b0) / as.matrix(th * 2)
  alfa <- s_bar - b0 * as.matrix(th) / 2 - gamma * as.matrix(th)^2/3

  # just return the stuff needed for subsequent steps (decompose to vec or nah?)
  list("s_bar" = s_bar, "b0" = b0, "b1" = b1, "gamma" = gamma, "alfa" = alfa)
}

# Note: Mass-preserving spline explained in detail in [http://dx.doi.org/10.1016/S0016-7061(99)00003-8];

# Spline fitting for horizon data (created by Brendan Malone; adjusted by T. Hengl)
mpspline <- function(obj = NULL, var.name = NULL, lam = 0.1,
                     d = t(c(0,5,15,30,60,100,200)),
                     vlow = 0, vhigh = 1000, show_progress = TRUE) {

  nice_obj <- mpspline_conv(obj)

  # find the max number of samples/horizons across all input sites
  max_dat <- max(table(nice_obj[[1]]))

  # *shrug* it could happen I guess??
  if(max_dat == 1) {
    stop('Supplied sites all have one horizon, please check inputs.')
  }

  ## organize the data:
  ndata <- nrow(objd)
  mxd   <- max(d)

  # 5% of the standard deviation of the target attribute
  s <- 0.05 * sd(unlist(unclass(objd[, svar.lst])), na.rm = TRUE)
  s2 <- s*s   # overall variance of soil attribute

  # split input data into a list by site
  sites <- split(nice_obj, as.factor(nice_obj[[1]]))

  # do some checks and tidying up
  sites <- mpspline_datchk(sites)


  # estimate spline parameters for each site  # prog bar it later maybe
  params <- lapply(sites[1:10], function(i) mpspline_est1(i, lam = lam))



                  ## fit the spline
                  ## spline will be interpolated onto these depths (1cm res)
                  xfit <- as.matrix(t(c(1:mxd)))
                  nj   <- max(v)
                  if (nj > mxd) { nj <- mxd }
                  yfit<- xfit
                  for (k in 1:nj) {
                    xd <- xfit[k]
                    if (xd < u[1]) {
                      p <- alfa[1]
                    } else {
                      for(its in 1:n) {
                        if(its < n) {
                          tf2 = as.numeric(xd > v[its] & xd < u[its + 1])
                        } else {
                          tf2 <- 0
                        }

                        if (xd >= u[its] & xd <= v[its]) {
                          p = alfa[its] + b0[its] * (xd - u[its]) + gamma[its] * (xd - u[its]) ^ 2
                        } else if(tf2) {
                          phi = alfa[its + 1] - b1[its] * (u[its + 1] - v[its])
                          p = phi + b1[its] * (xd - v[its])
                        }
                      }
                    }
                    yfit[k] = p
                  }
                  if (nj < mxd) { yfit[, (nj+1):mxd] = NA }

                  yfit[which(yfit > vhigh)] <- vhigh
                  yfit[which(yfit < vlow)]  <- vlow
                  m_fyfit[st, ] <- yfit

                  ## Averages of the spline at specified depths
                  nd <- length(d) - 1  # number of depth intervals
                  dl <- d + 1     #  increase d by 1
                  for (cj in 1:nd) {
                    xd1 <- dl[cj]
                    xd2 <- dl[cj + 1] - 1
                    if (nj >= xd1 & nj <= xd2) {
                      xd2 <- nj - 1
                      yave[st, cj] <- mean(yfit[, xd1:xd2])
                      } else {
                        yave[st, cj] <- mean(yfit[, xd1:xd2])
                        } # average of the spline at the specified depth intervals
                    yave[st, cj + 1] <- max(v)
                    } #maximum soil depth

                  ## Spline estimates at observed depths
                  dave[st, 1:n] <- sbar

                  ## CALCULATION OF THE ERROR BETWEEN OBSERVED AND FITTED VALUES
                  ## calculate Wahba's estimate of the residual variance sigma^2
                  ssq <- sum((t(y) - sbar)^2)
                  g <- solve(z)
                  ei <- eigen(g)
                  ei <- ei$values
                  df <- n - sum(ei)
                  sig2w <- ssq / df
                  ## calculate the Carter and Eagleson estimate of residual variance
                  dfc <- n - 2 * sum(ei) + sum(ei^2)
                  sig2c <- ssq / dfc
                  ## calculate the estimate of the true mean squared error
                  tmse <- ssq / n - 2 * s2 * df/ n + s2
                  sset[st] <- tmse

                }
              }

              if (show.progress) { setTxtProgressBar(pb, st)  }
            }
            if (show.progress) {
              close(pb)
              #cat(st, "\r")  ## TH: Suggested by D. Rossiter but not required
              #flush.console()
            }

            ## asthetics for output
            ## yave
            yave<- as.data.frame(yave)
            jmat<- matrix(NA, ncol = 1, nrow = length(d))
            for (i in 1:length(d) - 1) {
              a1 <- paste(d[i], d[ i + 1], sep= "-")
              a1 <- paste(a1, "cm", sep=" ")
              jmat[i] <- a1
              }
            jmat[length(d)] <- "soil depth"
            for (jj in 1:length(jmat)) {
              names(yave)[jj] <- jmat[jj]
            }

            retval <- list(idcol = objd_m[, 1],
                           var.fitted = dave,
                           var.std = yave,
                           var.1cm = t(m_fyfit))

            return(retval)
}

# end of script;

