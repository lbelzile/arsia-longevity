#####################################################
##### Section 4.3 - Is there a Cap on Longevity? ####
#####################################################
#### Data for England & Wales from IDL + ONS     ####
#### the database "englandwales" is described    ####
#### in the longevity package; type              ####
#### ?longevity::englandwales for a description  ####
#####################################################

# Load R packages
library(longevity) # analysis of excess lifetime
library(dplyr) #data manipulation
library(lubridate) #dates
library(ggdist) #Bayesian visualization
library(rust) # ratio-of-uniform sampling algorithm
library(ggplot2) # grammar of graphics
library(patchwork) # combine ggplot objects 
library(cobs) # robust spline
# Logical switch for uncertainty quantification
uncertainty <- FALSE
theme_set(theme_bw() +
            theme(axis.line = element_line(colour = "black"),
                  # panel.grid.minor = element_blank(),
                  panel.border = element_blank(),
                  panel.background = element_blank()))

# Combine databases (lazy-loaded)
ew_dat <- englandwales %>%
    mutate(nyears = ndays / 365.25,
           ltrunc = ltrunc / 365.25,
           rtrunc = rtrunc / 365.25)
# Maximum records is (to rounding) 115 years
max(ew_dat$nyears)
# Threshold stability plots
thresh_ew_dat <- with(ew_dat,
                  tstab(time = nyears,
                        thresh = 105:111,
                        ltrunc = ltrunc,
                        rtrunc = rtrunc,
                        family = "gp",
                        plot.type = "ggplot",
                        method = "profile")
)
# Non-exponential behaviour for 105-108 and 110
fit_109_exp <- with(ew_dat,
               fit_elife(time = nyears,
                         thresh = 109,
                         ltrunc = ltrunc,
                         rtrunc = rtrunc,
                         family = "exp",
                         export = TRUE))

source("EnglandWales_analysis_supp.R")
# Compute profile likelihood for exponential scale parameter
profile_exp <-
  with(fit_109_exp,
       prof_exp_scale(
         mle = fit_109_exp,
         time = time,
         thresh = 0,
         ltrunc = ltrunc,
         rtrunc = rtrunc)
  )

# How do we get uncertainty quantification?
# There are two sources: parameter estimates
# and sampling uncertainty.
# For the latter, we sample parameter values
# from the normal approximation to the profile log-likelihood
# For the former, we re-estimate the model on bootstrap samples
# whose starting values are drawn from the approximation
# 
# Bootstrap for simulating data from the fitted model
# For each replicate dataset, we re-estimate the 
# model parameters and the plotting positions 
# (nonparametric) MLE.
  set.seed(202104)
  B <- 9999L
  uq_qqplot <- with(fit_109_exp,
    uq_qqplot_elife_exp(B = B,
                        n = fit_109_exp$nexc,
                        par = pred_fun(runif(B)),
                        lower = fit_109_exp$ltrunc,
                        upper = fit_109_exp$rtrunc)
  )
  # Compute pointwise 90% confidence intervals via cubic spline
  lb <- cobs::cobs(x = c(uq_qqplot$x), y = c(uq_qqplot$y), constraint = "concave", tau = 0.05)
  ub <- cobs::cobs(x = c(uq_qqplot$x), y = c(uq_qqplot$y), constraint = "increase", tau = 0.95)
  save(uq_qqplot, fit_109_exp, lb, ub, file = "bootstrap_qqplot_exp_109_EW.RData")


# Check the empirical coverage of these approximate intervals
qqplot <- plot(fit_109_exp, plot.type = "ggplot", plot = FALSE)
mean(predict(ub, qqplot$qq$data[,1])[,2] > qqplot$qq$data[,2] & predict(lb, qqplot$qq$data[,1])[,2] < qqplot$qq$data[,2])
# Q-Q plot with 90% confidence intervals
g1 <- ggplot(data = qqplot$qq$data,
             aes(x = x, y = y)) +
  geom_abline(slope = 1, intercept = 0, color = "grey", alpha = 0.5) +
  geom_point(pch = 20) +
  geom_line(
  data = as.data.frame(predict(lb, seq(0, 7, by= 0.1))),
  aes(x = z, y = fit), alpha = 0.5) +
  geom_line(
    data = as.data.frame(predict(ub, seq(0, 7, by= 0.1))),
    aes(x = z, y = fit), alpha = 0.5) +
  scale_x_continuous("theoretical quantiles",
                     breaks = c(0,2,4,6),
                     labels = c(109,111,113,115),
                     limits = c(0,6.5),
                     expand = c(0,0)) +
  scale_y_continuous("empirical quantiles",
                     breaks = c(0,2,4,6),
                     labels = c(109,111,113,115),
                     limits = c(0,6.5),
                     expand = c(0,0))

# Fit various parametric models to the data 
fit_exp <- list()
fit_gp <- list()
fit_gomp <- list()
thresh <- 105:111
for(i in seq_along(thresh)){
  mod1 <- with(ew_dat,
                 fit_elife(time = nyears,
                           thresh = thresh[i],
                           ltrunc = ltrunc,
                           rtrunc = rtrunc,
                           family = "exp",
                           export = TRUE))
  fit_exp[[i]] <- mod1
  mod2 <- with(ew_dat,
                 fit_elife(time = nyears,
                           thresh = thresh[i],
                           ltrunc = ltrunc,
                           rtrunc = rtrunc,
                           family = "gomp",
                           export = TRUE))
  fit_gomp[[i]] <- mod2
  mod3 <- with(ew_dat,
                   fit_elife(time = nyears,
                             thresh = thresh[i],
                             ltrunc = ltrunc,
                             rtrunc = rtrunc,
                             family = "gp",
                             export = TRUE))
  try(print(anova(mod1, mod2))) # Gompertz vs exp
  try(print(anova(mod1, mod3))) # GP vs exp
}
# Bundle all fit together
loglik_df <- data.frame(
  threshold = thresh,
  gp = as.numeric(lapply(fit_gp, logLik)),
  gomp = as.numeric(lapply(fit_gomp, logLik)),
  exp = as.numeric(lapply(fit_exp, logLik)))

#########################################################
# Hazard plots (right panel of Figure 5)
# Create a wrapper function for log posterior
logpost_fn <- function(par, thresh, family){
    with(ew_dat, lpost_elife(par = par,
                                 time = nyears,
                                 ltrunc = ltrunc,
                                 rtrunc = rtrunc,
                                 family = family,
                                 thresh = thresh))
  }
# Set number of simulations and the thresholds
nB <- 1e4L
th_haz <- 109
set.seed(202103)
# The ratio-of-uniform sampling algorithm requires 
# bounded log-posterior densities
# (potentially after Box-Cox transform)
# Gompertz, exponential and generalized Pareto
postsamp_gomp <- rust::ru(logf = logpost_fn,
                          thresh = th_haz,
                          d = 2,
                          family = "gomp",
                          n = nB,
                          init = boxcox_transfo(fit_gomp[[4]]$par, lambda = c(1, 0.3)),
                          lower = c(1e-5, 1e-5),
                          a_method = "N",
                          b_method = "N",
                          trans = "BC",
                          lambda = c(1,0.3))
postsamp_gp <- rust::ru(logf = logpost_fn,
                        family = "gp",
                        thresh = th_haz,
                        d = 2,
                        n = nB,
                        init = as.numeric(fit_gp[[4]]$par),
                        lower = c(0, -1),
                        trans = "none")
postsamp_exp <- rust::ru(logf = logpost_fn,
                         family = "exp",
                         thresh = th_haz,
                         d = 1,
                         n = nB,
                         init = as.numeric(fit_exp[[4]]$par),
                         lower =  0,
                         trans = "none")

# Create a data frame combining the hazard estimates pointwise
# for each posterior sample
xseq <- seq(0.01, 15, length.out = 1001L)
exp_haz_df <- tibble(
  draw = rep(1:nB, each = length(xseq)),
  x = rep(xseq+th_haz, length.out = length(xseq)*nB),
  y = c(sapply(1:nB, function(ind){hazard_fn_elife(x = xseq, par = postsamp_exp$sim_vals[ind,1], family = "exp")})))
gomp_haz_df <- tibble(
  draw = rep(1:nB, each = length(xseq)),
  x = rep(xseq+th_haz, length.out = length(xseq)*nB),
  y = c(sapply(1:nB, function(ind){hazard_fn_elife(x = xseq, par = postsamp_gomp$sim_vals[ind,], family = "gomp")})))
gp_haz_df <- tibble(
  draw = rep(1:nB, each = length(xseq)),
  x = rep(xseq+th_haz, length.out = length(xseq)*nB),
  y = c(sapply(1:nB, function(ind){hazard_fn_elife(x = xseq, par = postsamp_gp$sim_vals[ind,], family = "gp")})))

# Plot the hazards with curvewise bands 
exp_curvwise_conf <- exp_haz_df %>%
  group_by(x) %>%
  curve_interval(y, .width = .5) %>%
  mutate(distribution = factor("exponential"))
gomp_curvwise_conf <- gomp_haz_df %>%
  group_by(x) %>%
  curve_interval(y, .width = c(.5)) %>%
  mutate(distribution = factor("Gompertz"))
gp_curvwise_conf <- gp_haz_df %>%
  group_by(x) %>%
  curve_interval(y, .width = c(.5), na.rm = TRUE) %>%
  mutate(distribution = factor("generalized Pareto"))
curvwise <- bind_rows(exp_curvwise_conf, gomp_curvwise_conf, gp_curvwise_conf)

g2 <-  ggplot(data = curvwise, aes(x = x, y = y, 
                                   color = distribution, 
                                   fill = distribution)) +
  geom_lineribbon(data = curvwise,
                  aes(ymin = .lower, ymax = .upper),
                  alpha = 0.1, show.legend = FALSE) +
  geom_line(aes(lty = distribution), lwd = 1.5) +
  coord_cartesian(ylim = c(0, 4),
                  expand = FALSE,
                  xlim = th_haz + c(0, 8.5)) +
  labs(x = "age",
       y = "hazard") +
  scale_color_manual(values=c("#060606", "#E69F00", "#56B4E9")) +
  scale_fill_manual(values=c("#060606", "#E69F00", "#56B4E9")) +
  theme(legend.position = c(0.25,0.8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
# Add the two figures together
g1 + g2

# Likelihood ratio tests for differences between gender
# Fit different parametric models, with same parameters for
# both gender (H0) and same distributions, but different 
# parameters (H1)
th <- 105:110
gender_ew_dat <- list()
for(i in seq_along(th)){
  gender_ew_dat[[i]] <-
    data.frame(thresh = th[i],
               gp = with(ew_dat,
                         test_elife(covariate = gender,
                                    time = nyears,
                                    ltrunc = ltrunc,
                                    rtrunc = rtrunc,
                                    family =  "gp",
                                    thresh = th[i])$pval),
               gomp = with(ew_dat,
                       test_elife(covariate = gender,
                                  time = nyears,
                                 ltrunc = ltrunc,
                                 rtrunc = rtrunc,
                                 family = "gomp",
                                  thresh = th[i])$pval),
               exp = with(ew_dat,
                          test_elife(covariate = gender,
                                     time = nyears,
                                     ltrunc = ltrunc,
                                     rtrunc = rtrunc,
                                     family = "exp",
                                     thresh = th[i])$pval)
    )
}
pvals_gender_df <- dplyr::bind_rows(gender_ew_dat)


# England and Wales data
# Illustration of likelihood ratio tests for
# nested models with non-regular asymptotic
# null distributions
B <- 9999L
set.seed(202103)
boot_test_list <- list()
for(i in 1:4){
  j <- i + 3
  # Extract truncation bounds
  boot_ltrunc <- with(fit_exp[[j]], ltrunc)
  boot_rtrunc <- with(fit_exp[[j]], rtrunc)
  mod_0 <- fit_exp[[j]]
  mod_1 <- fit_gomp[[j]]
  lrt_stat <- vector(mode = "numeric", length = B + 1)
  lrt_stat[1] <- 2*(mod_1$loglik - mod_0$loglik)
  for(b in 1:B){
    # Simulate new bootstrap sample from null model (exponential)
    boot_samp <- with(fit_exp[[j]],
      samp_elife(n = nexc,
                 lower = ltrunc,
                 upper = rtrunc,
                 scale = par,
                 family = "exp",
                 type2 = "ltrt")
    )
    # Fit both null (exponential) and alternative (Gompertz) models
    boot_fit_gomp <- fit_elife(time = boot_samp,
                          ltrunc = boot_ltrunc,
                          rtrunc = boot_rtrunc,
                          thresh = 0,
                          family = "gomp")
    boot_fit_exp <- fit_elife(time = boot_samp,
                              ltrunc = boot_ltrunc,
                              rtrunc = boot_rtrunc,
                              thresh = 0,
                              family = "exp")
    # Compute the likelihood ratio statistic - can be "zero"
    lrt_stat[b+1] <- 2*pmax(0,boot_fit_gomp$loglik - boot_fit_exp$loglik)
    }
    # Bootstrap p-value is the proportion as extreme as the value calculated
    # which is stored in the first row of the tibble
    boot_test_list[[i]] <- tibble(
      threshold = mod_0$thresh,
      stat = lrt_stat[1],
      nexc = mod_0$nexc,
      boot = 1-rank(lrt_stat)[1]/(B+1),
      asympt = anova(mod_0, mod_1)$`Pr(>Chisq)`[2])
}
boot_test_df <- dplyr::bind_rows(boot_test_list)
