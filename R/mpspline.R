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
  # data quality checks. seq_along() allows access to list names (SIDs)
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
    s <- s[-which(s[[2]] < 0 | s[[3]] < 0), ]

    # sort by cols 1, 2, 3 asc
    s <- s[order(s[[1]], s[[2]], s[[3]]), ]
    rownames(s) <- NULL

    # warn for overlapping data depth ranges
    if(any(diff(as.vector(rbind(s[[2]], s[[3]]))) < 0)) {
      stop("Overlapping horizons detected in site ", names(sites)[i], '.')
    }

    s
    })
}

# Note: Mass-preserving spline explained in detail in [http://dx.doi.org/10.1016/S0016-7061(99)00003-8];

# Spline fitting for horizon data (created by Brendan Malone; adjusted by T. Hengl)
mpspline <- function(obj = NULL, var.name = NULL, lam = 0.1,
                     d = t(c(0,5,15,30,60,100,200)),
                     vlow = 0, vhigh = 1000, show.progress=TRUE) {

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




            ## Fit splines profile by profile:
            message("Fitting mass preserving splines per profile...")
            if (show.progress)
              pb <- txtProgressBar(min = 0,
                                   max = length(sel),
                                   style = 3)
            for(st in as.vector(which(sel))) {
              subs <- matrix(unlist(c(1:np,
                                      as.vector(objd_m[st, upperb.lst]),
                                      as.vector(objd_m[st, lowerb.lst]),
                                      as.vector(objd_m[st, svar.lst]))), ncol = 4)
              d.ho <-
                rowMeans(data.frame(x = subs[, 2],
                                    y = c(NA, subs[1:(nrow(subs) - 1), 3])),
                         na.rm = TRUE)
              ## mask out missing values
              if (ncol(as.matrix(subs[!is.na(subs[, 2]) &
                                      !is.na(subs[, 3]) & !is.na(subs[, 4]), ]))==1)
              {
                subs = t(as.matrix(subs[!is.na(subs[, 2]) &
                                          !is.na(subs[, 3]) & !is.na(subs[, 4]), ]))
              }
              else {
                subs <- subs[!is.na(subs[, 2]) & !is.na(subs[, 3]) &
                               !is.na(subs[, 4]), ]
              }

              ## manipulate the profile data to the required form
              ir <- c(1:length(subs[,1]))
              ir <- as.matrix(t(ir))
              u <- subs[ir,2]
              u <- as.matrix(t(u))   # upper
              v <- subs[ir,3]
              v <- as.matrix(t(v))   # lower
              y <- subs[ir,4]
              y <- as.matrix(t(y))   # concentration
              n <- length(y)       # number of observations in the profile

              ############################################################################################################################################################
              ## routine for handling profiles with one observation
              if (n == 1){
                message(paste("Spline not fitted to profile:",
                              objd_m[st, 1], sep = " "))
                ## spline will be interpolated onto these depths (1cm res)
                xfit<- as.matrix(t(c(1:mxd)))
                nj<- max(v)
                if (nj > mxd) { nj <- mxd }
                yfitv <- xfit
                yfit[, 1:nj] <- y   ## values extrapolated onto yfit
                if (nj < mxd) { yfit[, (nj + 1):mxd] = NA }
                m_fyfit[st, ] <- yfit

                ## Averages of the spline at specified depths
                nd <- length(d) - 1  ## number of depth intervals
                dl <- d + 1     ##  increase d by 1
                for (cj in 1:nd) {
                  xd1 <- dl[cj]
                  xd2 <- dl[cj + 1] - 1
                  if (nj >= xd1 & nj <= xd2) {
                    xd2 <- nj - 1
                    yave[st, cj] <- mean(yfit[, xd1:xd2])
                    } else {
                      # average of the spline at the specified depth intervals
                      yave[st, cj] <- mean(yfit[,xd1:xd2])
                    }
                  yave[st, cj + 1] <- max(v)
                  } # maximum soil depth
              }

              ## End of single observation profile routine
              ###############################################################################################################################################################

              ## Start of routine for fitting spline to profiles with multiple observations

              else  {
                ###############################################################################################################################################################
                ## ESTIMATION OF SPLINE PARAMETERS
                np1 <- n + 1  # number of interval boundaries
                nm1 <- n - 1
                delta <- v - u  # depths of each layer
                del <- c(u[2:n], u[n]) - v   # del is (u1-v0,u2-v1, ...)

                ## create the (n-1)x(n-1) matrix r; first create r with 1's on
                ## the diagonal and upper diagonal, and 0's elsewhere
                r <- matrix(0, ncol = nm1, nrow = nm1)
                for(dig in 1:nm1){
                  r[dig, dig] <- 1
                }
                for(udig in 1:nm1 - 1){
                  r[udig, udig + 1] <- 1
                }

                ## then create a diagonal matrix d2 of differences to
                ## premultiply the current r
                d2 <- matrix(0, ncol = nm1, nrow = nm1)
                diag(d2) <- delta[2:n]  # delta = depth of each layer

                ## then premultiply and add the transpose; this gives half of r
                r <- d2 %*% r
                r <- r + t(r)

                ## then create a new diagonal matrix for differences to add to
                ## the diagonal
                d1 <- matrix(0, ncol = nm1, nrow = nm1)
                diag(d1) <- delta[1:nm1]  # delta = depth of each layer

                d3 <- matrix(0, ncol = nm1, nrow = nm1)
                diag(d3) <- del[1:nm1]  # del =  differences

                r <- r+2*d1 + 6*d3

                ## create the (n-1)xn matrix q
                q <- matrix(0, ncol = n, nrow = n)
                for (dig in 1:n){
                  q[dig, dig] <- -1
                }
                for (udig in 1:n - 1){
                  q[udig, udig + 1] <- 1
                }
                q <- q[1:nm1, 1:n]
                dim.mat <- matrix(q[], ncol = length(1:n),nrow = length(1:nm1))

                ## inverse of r
                rinv <- try(solve(r), TRUE)

                ## Note: in same cases this will fail due to singular matrix
                ## problems, hence you need to check if the object is
                ## meaningfull:
                if(is.matrix(rinv)){
                  ## identity matrix i
                  ind <- diag(n)

                  ## create the matrix coefficent z

                  pr.mat <- matrix(0, ncol = length(1:nm1), nrow = length(1:n))
                  pr.mat[] <- 6 * n * lam
                  fdub <- pr.mat * t(dim.mat) %*% rinv
                  z <- fdub %*% dim.mat + ind

                  ## solve for the fitted layer means
                  sbar <- solve(z, t(y))

                  ## calculate the fitted value at the knots
                  b <- 6 * rinv %*% dim.mat %*% sbar
                  b0 <- rbind(0, b) # add a row to top = 0
                  b1 <- rbind(b, 0) # add a row to bottom = 0
                  gamma <- (b1 - b0) / t(2 * delta)
                  alfa <- sbar - b0 * t(delta) / 2 - gamma * t(delta)^2/3

                  ## END ESTIMATION OF SPLINE PARAMETERS
                  ##############################################################


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

