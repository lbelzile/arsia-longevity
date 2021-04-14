# Additional functions used in 4.3-EnglandWales_analysis.R

#' @inheritParams nll_elife
#' @param mle an object of class \code{elife_par}
#' @return a smooth spline.spline object
prof_exp_scale <- function(mle = NULL,
                           time,
                           time2 = NULL,
                           event = NULL,
                           thresh = 0,
                           ltrunc = NULL,
                           rtrunc = NULL,
                           type = c("right","left","interval","interval2"),
                           weights = NULL){
  type <- match.arg(type)
  stopifnot("Provide a single threshold" = length(thresh) == 1L)
  if(is.null(weights)){
    weights <- rep(1, length(time))
  }
  if(!is.null(mle)){
    stopifnot("`mle` should be an object of class `elife_par` as returned by `fit_elife`" =  inherits(mle, "elife_par"))
  } else{
    mle <- fit_elife(time = time,
                     time2 = time2,
                     event = event,
                     thresh = thresh,
                     ltrunc = ltrunc,
                     rtrunc = rtrunc,
                     type = type,
                     family = "exp",
                     weights = weights)
  }
  psi <- mle$par + seq(pmax(-mle$par + 1e-4, -4*mle$std.error), 4*mle$std.error, length.out = 201)
  pll <- sapply(psi, function(scale){
    nll_elife(par = scale,
              time = time,
              time2 = time2,
              event = event,
              thresh = thresh,
              type = type,
              ltrunc = ltrunc,
              rtrunc = rtrunc,
              family = "exp",
              weights = weights)
  })
  conf_interv_fn <- function(object,
                             ...){
    args <- list(...)
    if ("warn" %in% names(args) && is.logical(args$warn)) {
      warn <- args$warn
    }  else {
      warn <- TRUE
    }
    if (is.null(object$pll) && is.null(object$r)) {
      stop("Object should contain arguments `pll` or `r` in order to compute confidence intervals.")
    }
    if (is.null(object$r)) {
      object$r <- sign(object$psi.max - object$psi) *
        sqrt(2 * (object$maxpll - object$pll))
    } else {
      object$r[is.infinite(object$r)] <- NA
    }
    if (is.null(object$normal)) {
      object$normal <- c(object$psi.max, object$std.error)
    }
    fit.r <- stats::smooth.spline(x = na.omit(cbind(object$r,
                                                    object$psi)), cv = FALSE)
    return(fit.r)
  }
  conf_interv_fn <- conf_interv_fn(list(psi = psi,
                                        pll = -2*(mle$loglik + pll),
                                        maxpll = 0,
                                        psi.max = mle$par,
                                        std.error = mle$std.error,
                                        mle = mle$par))
  return(conf_interv_fn)
}

pred_fun <- function(x, fit.spline = profile_exp){
  predict(fit.spline, qnorm(x))$y
}
#' Uncertainty quantification for quantile-quantile plots
#' @export
#' @param B number of bootstrap samples
#' @param dat vector of data
#' @param par parameter estimates of the model
#' @param lower lower bounds (truncation or lowest possible value)
#' @param upper upper bound for right-censoring or right-truncation
#' @param level level of the confidence intervals
#' @inheritParams nll_elife
#' @keywords internal
uq_qqplot_elife_exp <-
  function(B = 9999L,
           n,
           par,
           lower,
           upper
  ){
    family <- "exp"
    type2 <- "ltrt"
    if(!missing(lower) && !missing(upper)){
      if(length(upper) != 1 && length(upper) != 1){
        stopifnot( "`upper` and `lower` should be vectors of the same length." = length(lower) == length(upper),
                   "`Length of data `dat` does not match vector of lower and upper bounds." = n == length(upper))
      }
    }
    if(missing(lower)){
      ltrunc <- 0
    } else{
      ltrunc <- lower
    }
    if(missing(upper)){
      rtrunc <- Inf
    } else{
      rtrunc <- upper
    }
    stopifnot("`lower` should be smaller than `upper`." = isTRUE(all(lower <= upper)),
              "The matrix of parameters should have B rows." = length(par) == B)
    # parametric bootstrap samples
    # - simulate new data with the same sampling scheme
    # - estimate parameters of the distribution
    # - compute quantiles corresponding to plotting positions
    xppos <- matrix(NA, nrow = B, ncol = n)
    yppos <- matrix(NA, nrow = B, ncol = n)
    for(b in seq_len(B)){
      boot_dat <- samp_elife(n = n,
                             scale = par[b],
                             lower = lower,
                             upper = upper,
                             family = family,
                             type2 = type2)
      yppos[b,] <- sort(boot_dat)
      fit_boot <- fit_elife(time = boot_dat,
                            ltrunc = lower,
                            rtrunc = upper,
                            family = family)
      np_boot <- npsurv(time = boot_dat,
                        type = "interval",
                        event = rep(1L, n),
                        ltrunc = lower,
                        rtrunc = upper)$cdf
      xpos <- n / (n + 1) * (np_boot(boot_dat) - np_boot(ltrunc))/(np_boot(rtrunc) - np_boot(ltrunc))
      xppos[b,] <- qelife(p = pmax(0, pelife(q = ltrunc,
                                             scale = fit_boot$par,
                                             family = family)*(1-xpos) +
                                     xpos * pelife(q = rtrunc,
                                                   scale = fit_boot$par,
                                                   family = family)),
                          scale = fit_boot$par,
                          family = family)[order(boot_dat)]
    }
    return(list(x = xppos, y = yppos))
    # return(boot::envelope(mat = ppos, level = level))
  }
