#####################################################
##### Section 4.1 - Is there a Cap on Longevity? ####
#####################################################
#### Data for Netherland from Einmahl et al.     ####
#### the database "dutch" is available from the  ####
#### longevity package; see ?longevity::dutch    ####
#####################################################

## Packages and graphical options
library(longevity) # analysis of excess lifetime
library(tidyverse) # tidy data
library(lubridate) # dates
library(ggdist) # Bayesian visualization
library(rust) # ratio-of-uniform sampling algorithm
library(ggplot2) # grammar of graphics
library(patchwork) # combine ggplot objects
# Set theme for figures
theme_set(theme_bw() +
  theme(
    axis.line = element_line(colour = "black"),
    panel.border = element_blank(),
    panel.background = element_blank()
  ))
# Get data and rescale observations to years
# to avoid having parameters whose scale is
# order of magnitude different
data(dutch, package = "longevity")
yr_samp <- year(attr(x = dutch, which = "sampling_frame"))
dutch1 <- dutch %>%
  subset(!is.na(ndays)) %>%
  # Remove interval censored data
  mutate(
    time = ndays / 365.25,
    time2 = time,
    ltrunc = ltrunc / 365.25,
    rtrunc = rtrunc / 365.25,
    event = 1
  ) %>%
  subset(time > 98) %>%
  select(time, time2, ltrunc, rtrunc, event, gender, byear)
# Subset all interval-censored interval-truncated records
dutch2 <- dutch %>%
  subset(is.na(ndays)) %>%
  mutate(
    time2 = ceiling_date(dmy(paste("01-", dmonth, "-", dyear)), unit = "month") - 1 -
      dmy(paste("01-01-", byear)),
    time = dmy(paste("01-", dmonth, "-", dyear)) - dmy(paste("31-12-", byear)),
    ltrunc = dmy(paste("01-01-1986")) - dmy(paste("31-12-", byear)),
    rtrunc = dmy(paste("31-12-2015")) - dmy(paste("01-01-", byear))
  ) %>%
  select(time, time2, ltrunc, rtrunc, gender, byear) %>%
  mutate(
    time = as.numeric(time) / 365.25,
    time2 = as.numeric(time2) / 365.25,
    ltrunc = as.numeric(ltrunc) / 365.25,
    rtrunc = as.numeric(rtrunc) / 365.25,
    event = 3
  ) %>%
  subset(time > 98)
# Combine databases
dutch_data <- rbind(dutch1, dutch2)


# Threshold stability plot for the Dutch
# lifetime data with 95% profile-based
# pointwise confidence intervals.
# WARNING: computationally intensive for profile lik
tstab_c <- with(
  dutch_data,
  tstab(
    time = time,
    time2 = time2,
    event = event,
    type = "interval",
    thresh = 98:108,
    ltrunc = ltrunc,
    rtrunc = rtrunc,
    family = "gp",
    method = "profile",
    which.plot = "shape",
    plot = FALSE
  )
)
# Create plot
tstab_cp <- plot(tstab_c,
  plot.type = "ggplot",
  which.plot = "shape",
  plot = FALSE
)
save(tstab_cp, tstab_c, file = "dutch_qqplot.RData")


# Fit the Northrop-Coleman model
# the difference between the thresholds
# for the latter must be large enough for
# us to reliably estimate the shapes
# Exploratory analysis (not shown) suggested
# 97.5 with 2.5 years increments
# Specifically, halving the number of thresholds lead to a
# reduction of ~50% of the standard errors for mid-range shape

# Likelihood ratio tests to compare H0: GP versus
# H1: Northrop & Coleman model via anova()
thr <- seq(98, 107, by = 3)
lr_test <- list()
fit_gppiece <- list()
for (i in seq_len(length(thr) - 1L)) {
  dutch_nc_fit <- with(
    dutch_data,
    fit_elife(
      time = time,
      time2 = time2,
      event = event,
      type = "interval",
      thresh = thr[i:length(thr)],
      ltrunc = ltrunc,
      rtrunc = rtrunc,
      family = "gppiece",
      restart = TRUE
    )
  )
  fit_gppiece[[i]] <- dutch_nc_fit
  dutch_gp_fit <- with(
    dutch_data,
    fit_elife(
      time = time,
      time2 = time2,
      event = event,
      type = "interval",
      thresh = thr[i],
      ltrunc = ltrunc,
      rtrunc = rtrunc,
      family = "gp"
    )
  )
  # Likelihood ratio test
  anov <- try(anova(dutch_gp_fit, dutch_nc_fit))
  lr_test[[i]] <- tibble(
    thresh = thr[i],
    stat = anov$Chisq[2],
    pvalue = anov$`Pr(>Chisq)`[2]
  )
}
bind_rows(lr_test)

# Fit generalized Pareto and Gompertz models
# over a range of threshold
fit_gp <- list()
fit_gomp <- list()
th <- 98:108
for (i in seq_along(th)) {
  fit_gomp[[i]] <- with(
    dutch_data,
    fit_elife(
      time = time,
      time2 = time2,
      event = event,
      type = "interval",
      thresh = th[i],
      ltrunc = ltrunc,
      rtrunc = rtrunc,
      family = "gomp"
    )
  )
  fit_gp[[i]] <- with(
    dutch_data,
    fit_elife(
      time = time,
      time2 = time2,
      event = event,
      type = "interval",
      thresh = th[i],
      ltrunc = ltrunc,
      rtrunc = rtrunc,
      family = "gp"
    )
  )
}
# Extract some information from the fitted object
# and create a tibble
fit_df <- bind_rows(
  lapply(c(fit_gomp, fit_gp), function(x) {
    tibble(
      family = x$family,
      thresh = x$thresh,
      scale = x$par[1],
      shape = x$par[2],
      loglik = x$loglik
    )
  })
)

# Apply to selected thresholds to avoid visual
# clutter- leave space so there are enough data
# to reliably estimate shapes in intervals
dutch_nc_fit <- with(
  dutch_data,
  fit_elife(
    time = time,
    time2 = time2,
    event = event,
    type = "interval",
    thresh = c(98, 101, 104, 107),
    ltrunc = ltrunc,
    rtrunc = rtrunc,
    family = "gppiece",
    restart = TRUE
  )
)


# Penultimate shapes from the Northrop and Coleman model
# with Wald 95% confidence intervals (~ +/- 2 standard errors).
nc_shape_df <- with(
  dutch_nc_fit,
  tibble(
    xlow = thresh,
    xhigh = c(thresh[-1], Inf),
    shape = par[-1],
    ylow = par[-1] - qnorm(0.975) * std.error[-1],
    yhigh = par[-1] + qnorm(0.975) * std.error[-1]
  )
)

# Color-blind friendly palette
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

g1 <- tstab_cp$g2 +
  geom_hline(yintercept = 0, alpha = 0.1, lwd = 0.1) +
  geom_segment(
    data = nc_shape_df,
    aes(
      x = xlow,
      xend = xhigh,
      y = shape,
      yend = shape
    ), lty = 2
  ) +
  geom_segment(
    data = nc_shape_df,
    aes(
      x = xlow,
      xend = xhigh,
      y = ylow,
      yend = ylow
    ), alpha = 0.4, lty = 2
  ) +
  geom_segment(
    data = nc_shape_df,
    aes(
      x = xlow,
      xend = xhigh,
      y = yhigh,
      yend = yhigh
    ), alpha = 0.4, lty = 2
  ) +
  scale_y_continuous(
    breaks = seq(-0.2, 0.5, length.out = 8),
    minor_breaks = seq(-0.2, 0.5, length.out = 15),
    limits = c(-0.25, 0.5),
    labels = paste0("$", sprintf(-2:5 / 10, fmt = "%.1f"), "$")
  ) +
  scale_x_continuous(
    breaks = seq(98, 108, by = 2L),
    minor_breaks = seq(99L, 107, by = 2L),
    limits = c(97, 108.3), expand = c(0, 0)
  ) +
  theme(
    panel.grid = element_blank(),
    legend.position = c(0.1, 0.6)
  ) +
  coord_cartesian(
    xlim = c(97.5, 109),
    ylim = c(-0.25, 0.55),
    expand = FALSE,
    default = FALSE,
    clip = "on"
  ) +
  labs(x = "threshold", y = "shape")


# Bayesian analysis using ratio-of-uniform sampling
# and MDI priors
logpost_fn <- function(par, thresh, family) {
  with(dutch_data, lpost_elife(par,
    time = time,
    time2 = time2,
    event = event,
    type = "interval",
    thresh = thresh,
    ltrunc = ltrunc,
    rtrunc = rtrunc,
    family = family
  ))
}

# Bayesian analysis - more or less the same as frequentist
# given the prior is noninformative, but easier to sample from
set.seed(202104)
nB <- 1e3
postsamp_gomp <- list()
for (i in 1:length(th)) {
  postsamp_gomp[[i]] <-
    rust::ru(
      logf = logpost_fn,
      thresh = th[i],
      family = "gomp", d = 2,
      n = nB,
      init = boxcox_transfo(par = fit_gomp[[i]]$par, lambda = c(1, 0.33)),
      trans = "BC",
      rotate = FALSE,
      lambda = c(1, 0.33)
    )
}
# Create a wrapper to compute the penultimate
# approximation fn for the Gompertz model
gomp_penult_fn <- function(x, par, thresh) {
  -par[2] * exp(-par[2] * (x - thresh) / par[1])
}
# Combine lists into vectors
gomp_penult_df <- list()
for (i in 1:4) {
  j <- c(1L, 4L, 7L, 10L)[i]
  gomp_penult_df[[i]] <- tibble(
    draw = rep(1:nB, each = 101),
    x = rep(seq(from = th[j], to = 109, length.out = 101), length.out = 101L * nB),
    y = c(sapply(1:nB, function(ind) {
      gomp_penult_fn(
        x = seq(from = th[j], to = 109, length.out = 101),
        par = postsamp_gomp[[j]]$sim_vals[ind, ],
        thresh = th[j]
      )
    })),
    threshold = factor(th[j])
  )
}
gomp_penult_df <- bind_rows(gomp_penult_df)
gomp_penult_summary_df <-
  group_by(gomp_penult_df, x, threshold) %>%
  point_interval(
    .interval = qi,
    .point = median,
    .width = 0.95,
    .exclude = c("draw")
  )

g2 <- gomp_penult_summary_df %>%
  ggplot() +
  geom_lineribbon(
    mapping = aes(
      x = x,
      y = y,
      ymin = .lower,
      ymax = .upper,
      col = threshold,
      fill = threshold
    ),
    alpha = 0.3
  ) +
  scale_colour_manual(values = cbPalette) +
  scale_fill_manual(values = cbPalette) +
  geom_hline(yintercept = 0, col = "gray", alpha = 0.5) +
  geom_text(
    data = tibble(
      n = unlist(lapply(fit_gomp, nobs)),
      x = th,
      y = rep(c(0.45, 0.4), length.out = length(th))
    ),
    aes(x = x, y = y, label = n), size = 3
  ) +
  scale_y_continuous(
    breaks = seq(-0.2, 0.5, length.out = 8),
    minor_breaks = seq(-0.2, 0.5, length.out = 15),
    labels = paste0("$", sprintf(-2:5 / 10, fmt = "%.1f"), "$")
  ) +
  scale_x_continuous(
    breaks = seq(98, 108, by = 2L),
    minor_breaks = seq(99L, 107, by = 2L)
  ) +
  theme(
    panel.grid = element_blank(),
    legend.position = c(0.2, 0.6)
  ) +
  coord_cartesian(
    xlim = c(97.5, 109),
    ylim = c(-0.25, 0.5),
    expand = FALSE,
    default = FALSE,
    clip = "on"
  ) +
  labs(x = "threshold", y = "shape")

# Output plots for the paper
g1 + g2

# Differences between gender
th <- 100:108
gender_dutch <- list()
for (i in seq_along(th)) {
  gender_dutch[[i]] <-
    data.frame(
      thresh = th[i],
      gp = with(
        dutch_data,
        test_elife(
          covariate = gender,
          time = time,
          time2 = time2,
          event = event,
          type = "interval",
          ltrunc = ltrunc,
          rtrunc = rtrunc,
          family = "gp",
          thresh = th[i]
        )$pval
      ),
      gomp = with(
        dutch_data,
        test_elife(
          covariate = gender,
          time = time,
          time2 = time2,
          event = event,
          type = "interval",
          ltrunc = ltrunc,
          rtrunc = rtrunc,
          family = "gomp",
          thresh = th[i]
        )$pval
      ),
      exp = with(
        dutch_data,
        test_elife(
          covariate = gender,
          time = time,
          time2 = time2,
          event = event,
          type = "interval",
          ltrunc = ltrunc,
          rtrunc = rtrunc,
          family = "exp",
          thresh = th[i]
        )$pval
      )
    )
}
pvals_gender_df <- dplyr::bind_rows(gender_dutch)

th <- 105:108
barplot(table(with(dutch_data, byear[time > 105])))
cohorts_cut <- c(1860, 1885, 1892, 1900, 1911)
bcohort_dutch <- list()
for (i in seq_along(th)) {
  # Cohort effect
  bcohort_dutch[[i]] <- with(
    dutch_data %>% filter(time > th[i]),
    test_elife(
      covariate = factor(cut(byear, cohorts_cut)),
      time = time,
      time2 = time2,
      event = event,
      type = "interval",
      ltrunc = ltrunc,
      rtrunc = rtrunc,
      family = "gp",
      thresh = th[i]
    )$pval
  )
}
# P-values for four cohorts based on LRT above 105,...
# 0.10 0.12 0.41 0.33
unlist(bcohort_dutch)
