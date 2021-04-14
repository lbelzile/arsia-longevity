#####################################################
##### Section 4.2 - Is there a Cap on Longevity? ####
#####################################################
#### Data for Japan from Hanamaya & Sibuya; the  ####
#### database "japanese" is available from the   ####
#### longevity package; see ?longevity::japanese ####
#####################################################
library(longevity)
library(lubridate) # data manipulation
library(ggplot2)   # grammar of graphics
library(patchwork) # combine plots
library(tidybayes) # plots
library(tidyr)     # modify data frames (tidyverse)
library(dplyr)     # database manipulation
library(rust)      # ratio-of-uniform

# Specification for the graphics
theme_set(theme_bw() +
            theme(axis.line = element_line(colour = "black"),
                  # panel.grid.minor = element_blank(),
                  panel.border = element_blank(),
                  panel.background = element_blank()))

# Load data
data(japanese, package = "longevity")
# Birth cohorts of centenarian are included from
# 1947 onward (thus 1947) up until 2012 (male) and 2014 (female)
japan_dat <- with(japanese,
                data.frame(
                  time = age,
                  time2 = age + 364/365,
                  rtrunc = ifelse(gender == "female", 2015, 2013) - byear,
                  gender =  gender,
                  weight = count,
                  event = 3
                )
            )
# Two individuals originally in the table are excluded
# because they are not truncated and thus do not comply
# with the sampling mechanism
# Since they are named individuals, we use the full info
tdeath <- as.numeric(c(ymd("2015-04-01") - ymd("1898-03-05"),
                       ymd("2013-06-12") - ymd("1897-04-19")))/365.25
# Merge databases
japan_dat <- rbind(japan_dat,
      data.frame(time = tdeath,
           time2 = tdeath,
           rtrunc = Inf,
           gender = c("female","male"),
           weight = 1L,
           event = 1))


#####################################################
# Threshold stability plots
tstab <- with(japan_dat,
       tstab(time = time,
               time2 = time2,
               event = event,
               weights = weight,
               rtrunc = rtrunc,
               type = "interval2",
               family = "gp",
               method = "wald",
               plot.type = "ggplot",
               thresh = 100:112))
# We can use this output to calculate the various endpoints
# The scale parameters are transformed for the purpose of
# threshold stability plot, so we do not add the thresholds
-tstab$scale[,"estimate"]/tstab$shape[,"estimate"]

##############################################################
# Fit different parametric models for comparison purposes
fit_gp <- list()
fit_exp <- list()
fit_gomp <- list()
fit_extgp <- list()
th <- 100:112
for(i in seq_along(th)){
  fit_gp_i <- with(japan_dat,
     fit_elife(time = time,
          time2 = time2,
          status = event,
          weights = weight,
          rtrunc = rtrunc,
          type = "interval2",
          family = "gp",
          thresh = th[i]))
  fit_exp_i <- with(japan_dat,
                  fit_elife(time = time,
                            time2 = time2,
                            status = event,
                            weights = weight,
                            rtrunc = rtrunc,
                            type = "interval2",
                            family = "exp",
                            thresh = th[i]))
  fit_gomp_i <- with(japan_dat,
                      fit_elife(time = time,
                                time2 = time2,
                                status = event,
                                weights = weight,
                                rtrunc = rtrunc,
                                type = "interval2",
                                family = "gomp",
                                thresh = th[i]))
  fit_extgp_i <- with(japan_dat,
                    fit_elife(time = time,
                              time2 = time2,
                              status = event,
                              weights = weight,
                              rtrunc = rtrunc,
                              type = "interval2",
                              family = "exp",
                              thresh = th[i]))
  fit_gp[[i]] <- fit_gp_i
  fit_exp[[i]] <- fit_exp_i
  fit_gomp[[i]] <- fit_gomp_i
  fit_extgp[[i]] <- fit_extgp_i
}
# No evidence of exponentiality until 111 onwards
for(i in 1:length(fit_gp)){
  mod1 <- fit_gp[[i]]
  mod0 <-  fit_exp[[i]]
  print(anova(mod1, mod0)$`Pr(>Chisq)`[2])
}
# Same holds for Gompertz model (not shown)
# This is due to loss of power or to lack of difference
# Here none of the estimates of the Gompertz shape are null

########################################################################
# Profile log likelihood for the endpoint
prof_endpt_japan <- list()
levels <- c(0.5,0.75,0.9)
for(i in 1:3){
  prof_endpt_japan[[i]] <- data.frame(t(sapply(th, function(th){
    with(japan_dat,
         prof_gp_endpt(time = time,
                       time2 = time2,
                       event = event,
                       weights = weight,
                       rtrunc = rtrunc,
                       type = "interval2",
                       psi = seq(118L, 180L, by = 1L),
                       thresh = th,
                       level = levels[i],
                       confint = TRUE)
    )
  })), threshold = th, level = factor(levels[i])
  )
}
prof_endpt_japan_df <- dplyr::bind_rows(prof_endpt_japan)

g1 <- ggplot(data = prof_endpt_japan_df) +
  geom_interval(aes(x = threshold,
                    y = Estimate,
                    col = level,
                    ymin = Lower.CI,
                    ymax = Upper.CI),
                show_point = FALSE) +
  geom_point(data = prof_endpt_japan_df %>%
               filter(level == "0.5"),
             aes(x = threshold,
                 y = Estimate),
             pch = 3, size = 1.5) +
  scale_color_brewer(direction = -1) +
  scale_y_continuous(
    breaks = seq(120L, 150L, by = 10L),
    minor_breaks = NULL,
    limits = c(115,150),
    oob = scales::oob_keep,
    labels = seq(120L, 150L, by = 10L)) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = "none") +
  labs(x = "threshold",
       y = "human lifespan") +
  scale_x_continuous(
    breaks = seq(min(th), max(th), by = 2L),
    minor_breaks = NULL,
    limits = range(th) + c(-0.5,0.5)) +
  coord_cartesian(
    ylim = c(115L,150L),
    clip = "on"
  )

# Simulation study to assess loss of information due to
# interval censoring (privacy preserving, but wasteful)
# At the same time, we can compute bootstrap p-values
# for the different models
th_sim <- 110
# Keep the right truncation limits, but
# remove the two untruncated records
jap_sim <- japanese[japanese$age >= th_sim,]
byear_sim <- with(jap_sim, byear)
n_sim <- sum(jap_sim$count)
rtrunc_r_sim <- rep.int(x = ifelse(jap_sim$gender == "male",
                                 2013, 2015) - byear_sim,
          times = jap_sim$count)
param_sim <- fit_gp[[which(th_sim == 100:112)]]$par

set.seed(202103)
B <- 10000L # start small for now
boot <- matrix(nrow = B, ncol = 8)
for(i in seq_len(B)){
 bday <- sample.int(n = 364, size = n_sim, replace = TRUE)
 rtrunc_sim <- with(jap_sim,
                    ifelse(rep(gender == "male", times = count),
                       2013, 2015) - (bday / 365 + rep(byear, times = count)))
 samp_sim <- samp_elife(n = n_sim,
             scale = param_sim["scale"],
             shape = param_sim["shape"],
             upper = rtrunc_sim,
             lower = rep(0, n_sim),
             family = "gp",
             type2 = "ltrt")
 fit_f <- fit_elife(time = samp_sim,
                    rtrunc = rtrunc_sim,
                    family = "gp")
 fit_c <- fit_elife(time = floor(samp_sim),
                    time2 = ceiling(samp_sim),
                    rtrunc = rtrunc_r_sim,
                    event = rep(3L, n_sim),
                    type = "interval2",
                    family = "gp")
 boot[i,] <- c(fit_f$par, fit_f$std.error, fit_c$par, fit_c$std.error)
}
boot <- as.data.frame(boot)
colnames(boot) <- c("f_par_scale","f_par_shape",
                    "f_std_scale","f_std_shape",
                    "c_par_scale","c_par_shape",
                    "c_std_scale","c_std_shape")
save(boot, file = "japanese_sim.RData")
# What is the sampling variability of full relative to censored
# estimation for the endpoint?
# Super careful because there is positive shape estimates
# and potentially thus -scale/shape is negative -> inflates variance
with(boot, sqrt(var((f_par_scale/f_par_shape)[f_par_shape < 0]) / var((c_par_scale/c_par_shape)[c_par_shape < 0])))
endpt_sim <- th_sim - fit_gp[[which(th == th_sim)]]$par[1]/fit_gp[[which(th == th_sim)]]$par[2]
# Very small median negative bias
with(boot, median(th_sim - (f_par_scale/f_par_shape)[f_par_shape < 0])) - endpt
with(boot, median(th_sim - (c_par_scale/c_par_shape)[c_par_shape < 0])) - endpt
# Proportion of shape estimates that are positive
sum(boot$f_par_shape > 0)/B
sum(boot$c_par_shape > 0)/B
#Mean bias - trimmed mean,
with(boot, th_sim - mean((f_par_scale/f_par_shape)[f_par_shape < 0]) - endpt)
with(boot, th_sim - mean((c_par_scale/c_par_shape)[c_par_shape < 0]) - endpt)


# In contrast, there is barely any difference for the shape parameters
with(boot, sqrt(var(f_par_shape[f_par_shape < 0])/var(c_par_shape[c_par_shape < 0])))
with(boot, mean(f_std_shape/c_std_shape))

with(boot, sqrt(var(f_par_scale)/var(c_par_scale)))
with(boot, mean(f_std_scale/c_std_scale))

# Plot of the sampling distribution with full info and interval-censored records
g2 <- ggplot(data = data.frame(endpoint = with(boot, c(ifelse(f_par_shape > 0, Inf, th_sim - f_par_scale/f_par_shape),
                                                 ifelse(c_par_shape > 0, Inf, th_sim - c_par_scale/c_par_shape))),
                         type = factor(rep(c("full","interval censored"), each = B)))) +
    stat_sample_slabinterval(aes(y = endpoint, x = type),
                             slab_type = "histogram",
                             breaks = c(seq(115,160, by= 0.75), Inf),
                             .width = c(0.75,0.95)) +
  geom_hline(yintercept = endpt_sim) +
  scale_y_continuous(
    name = 'human lifespan',
    breaks = seq(120L, 150L, by = 10L),
    minor_breaks = NULL,
    limits = c(115,160),
    expand = c(0,0),
    oob = scales::squish_infinite,
    labels = seq(120L, 150L, by = 10L)) +
  coord_cartesian(
    ylim = c(115L,150L),
    clip = "on"
  )


# Output plots for paper
g1 + g2


# Differences between gender
gender_japan <- list()
for(i in seq_along(th)){
gender_japan[[i]] <-
  data.frame(thresh = th[i],
  gp = with(japan_dat,
  test_elife(covariate = gender,
             time = time,
             time2 = time2,
             event = event,
             weights = weight,
             rtrunc = rtrunc,
             type = "interval2",
             family = "gp",
             thresh = th[i])$pval
  ),
  gomp = with(japan_dat,
       test_elife(covariate = gender,
                  time = time,
                  time2 = time2,
                  event = event,
                  weights = weight,
                  rtrunc = rtrunc,
                  type = "interval2",
                  family = "gomp",
                  thresh = th[i])$pval
  ),
  exp = with(japan_dat,
       test_elife(covariate = gender,
                  time = time,
                  time2 = time2,
                  event = event,
                  weights = weight,
                  rtrunc = rtrunc,
                  type = "interval2",
                  family = "exp",
                  thresh = th[i])$pval
  )
  )
}
pvals_gender_df <- dplyr::bind_rows(gender_japan)
