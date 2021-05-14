# The tables in the Supplementary material were created
# using the Matlab toolbox LATool using only records of
# supercentenarian; as such, the truncation bounds were
# calculated using the date at which data collection starts
# and ends for supercentenarians (110+)
library(dplyr)
library(lubridate)
library(longevity)

# Database of March 27th, including all first three waves
data(idl2021, package = "longevity")
# Threshold 110, but exact number of days can vary
# from one individual to the next because of leap years
u110d <- 40176L
u110 <- 110
idl <- idl2021 %>%
  filter(ageyear >= 110) %>%
  mutate(
    ltrunc = pmax(u110d, as.integer(c1 - bdate)) / 365.25,
    rtrunc = as.integer(c2 - bdate) / 365.25,
    nyears = ndays / 365.25
  )
# Define regions for various tables
north_europe <- idl %>% filter(country %in% c("OS", "BE", "DN", "FI", "DE", "NO", "SV", "EW"))
south_europe <- idl %>% filter(country %in% c("FR", "SP"))
europe <- bind_rows(north_europe, south_europe)
north_america <- idl %>% filter(country %in% c("US", "QC"))
world <- bind_rows(europe, north_america)

# Counts by gender and region
rbind(
  with(north_europe, table(gender)),
  with(south_europe, table(gender)),
  with(europe, table(gender)),
  with(north_america, table(gender)),
  with(world, table(gender))
)
#---------------------------------------------------------
# Likelihood ratio tests comparing different
# generalized Pareto and different exponential models
# for women and men against the null hypothesis
# of the same distribution for women and men
u110y <- u110d / 365.25
NE_gender_gp <- with(
  north_europe,
  test_elife(
    time = nyears,
    thresh = u110y,
    ltrunc = ltrunc,
    rtrunc = rtrunc,
    family = "gp",
    covariate = gender
  )
)
# Summary of likelihood ratio test
print(NE_gender_gp)
# Same thing, but with exponential distribution
NE_gender_exp <-
  with(
    north_europe,
    test_elife(
      time = nyears,
      thresh = u110,
      ltrunc = ltrunc,
      rtrunc = rtrunc,
      family = "exp",
      covariate = gender
    )
  )
NE_gender_exp
#-------------------------------------------------------
# Likelihood ratio tests comparing different generalized
# Pareto and exponential distributions between groups of
# countries for IDL supercentenarian against the null of
# common distribution for both groups of countries.
test_country_group_gp <-
  with(
    bind_rows(
      europe %>% mutate(group = "europe"),
      north_america %>% mutate(group = "america")
    ),
    test_elife(
      time = nyears,
      thresh = u110,
      ltrunc = ltrunc,
      rtrunc = rtrunc,
      family = "gp",
      covariate = group
    )
  )
test_country_group_gp

#-------------------------------------------------------
# Comparison of generalized Pareto
# versus exponential distribution
fit_NE_gp <- with(
  north_europe,
  fit_elife(
    time = nyears,
    thresh = u110,
    ltrunc = ltrunc,
    rtrunc = rtrunc,
    family = "gp"
  )
)
fit_NE_gp
fit_NE_exp <- with(
  north_europe,
  fit_elife(
    time = nyears,
    thresh = u110,
    ltrunc = ltrunc,
    rtrunc = rtrunc,
    family = "exp"
  )
)
# Model comparison (likelihood ratio test)
anova(fit_NE_gp, fit_NE_exp)

#-------------------------------------------------------
# Fit different models for different death intervals:
# this requires changing the observation period and
# thus the truncation.
north_europe_split <-
  north_europe %>%
  mutate(
    deathgr = ddate > ymd("2006-12-31"),
    c1 = if_else(ddate > ymd("2006-12-31"), pmax(c1, ymd("2007-01-01")), c1),
    c2 = if_else(ddate <= ymd("2006-12-31"), pmin(c2, ymd("2006-12-31")), c2),
    ltrunc = as.integer(c1 - bdate) / 365.25,
    rtrunc = as.integer(c2 - bdate) / 365.25,
    nyears = ndays / 365.25
  )

world_split <-
  world %>%
  mutate(
    deathgr = ddate > ymd("2001-12-31"),
    c1 = if_else(ddate > ymd("2001-12-31"), pmax(c1, ymd("2002-01-01")), c1),
    c2 = if_else(ddate <= ymd("2001-12-31"), pmin(c2, ymd("2001-12-31")), c2),
    ltrunc = as.integer(c1 - bdate) / 365.25,
    rtrunc = as.integer(c2 - bdate) / 365.25,
    nyears = ndays / 365.25
  )
# In the paper, this comparison is done using a Wald test statistic
# We use instead the likelihood ratio test below
death_split_NE <-
  with(
    north_europe_split,
    test_elife(
      time = nyears,
      thresh = u110,
      ltrunc = ltrunc,
      rtrunc = rtrunc,
      family = "exp",
      covariate = deathgr
    )
  )

death_split_NE_1 <-
  with(
    north_europe_split %>% dplyr::filter(deathgr == TRUE),
    fit_elife(
      time = nyears,
      thresh = u110,
      ltrunc = ltrunc,
      rtrunc = rtrunc,
      family = "exp"
    )
  )

death_split_NE_2 <-
  with(
    north_europe_split %>% dplyr::filter(deathgr == FALSE),
    fit_elife(
      time = nyears,
      thresh = u110,
      ltrunc = ltrunc,
      rtrunc = rtrunc,
      family = "exp"
    )
  )
death_split_NE
death_split_NE_1
death_split_NE_2
wald_NE <- (death_split_NE_2$par - death_split_NE_1$par) / sqrt(death_split_NE_1$std.error^2 + death_split_NE_2$std.error^2)
pval_wald_NE <- 2*(1-pnorm(wald_NE))
death_split_WR <-
  with(
    world_split,
    test_elife(
      time = nyears,
      thresh = u110,
      ltrunc = ltrunc,
      rtrunc = rtrunc,
      family = "exp",
      covariate = deathgr
    )
  )

death_split_WR_1 <-
  with(
    world_split %>% dplyr::filter(deathgr == TRUE),
    fit_elife(
      time = nyears,
      thresh = u110,
      ltrunc = ltrunc,
      rtrunc = rtrunc,
      family = "exp"
    )
  )

death_split_WR_2 <-
  with(
    world_split %>% dplyr::filter(deathgr == FALSE),
    fit_elife(
      time = nyears,
      thresh = u110,
      ltrunc = ltrunc,
      rtrunc = rtrunc,
      family = "exp"
    )
  )
death_split_WR
death_split_WR_1
death_split_WR_2
wald_WR <- (death_split_WR_2$par - death_split_WR_1$par) / sqrt(death_split_WR_1$std.error^2 + death_split_WR_2$std.error^2)
pval_wald_WR <- 2*(1-pnorm(wald_WR))


#-------------------------------------------------------
# Analysis of data for countries other than
# France, Italy and England & Wales in the IDL
# that also include both semisupercentenarian
# and supercentenarians
#
# In the Supplementary material and for simplicity,
# the sampling frame is defined to be the smallest
# interval for which people can be seen. This is
# wasteful, but simplifies the analysis.
semi <- idl2021 %>%
  filter(country %in% c("OS", "BE", "DN", "DE", "QC")) %>%
  mutate(
    lb = pmax(c1, d1),
    ub = pmin(c2, d2),
    ltrunc = as.integer(lb - bdate) / 365.25,
    rtrunc = as.integer(ub - bdate) / 365.25,
    nyears = ndays / 365.25
  ) %>%
  filter(ddate >= lb & ddate <= ub) %>%
  select(nyears, ltrunc, rtrunc, country)
thyr <- 105:110
vfit <- tibble(
  par = numeric(6),
  std = numeric(6),
  nexc = integer(6),
  pval = numeric(6)
)
for (i in seq_along(thyr)) {
  mod0 <- with(
    semi,
    fit_elife(
      time = nyears,
      ltrunc = ltrunc,
      rtrunc = rtrunc,
      family = "exp",
      thresh = thyr[i]
    )
  )
  mod1 <- with(
    semi,
    fit_elife(
      time = nyears,
      ltrunc = ltrunc,
      rtrunc = rtrunc,
      family = "gp",
      thresh = thyr[i]
    )
  )
  vfit$par[i] <- coef(mod0)
  vfit$std[i] <- mod0$std.error
  vfit$nexc[i] <- nobs(mod1)
  vfit$pval[i] <- anova(mod0, mod1)$`Pr(>Chisq)`[2]
}
# @longevity can handle the double interval-truncation
semi_full <- idl2021 %>%
  filter(country %in% c("OS", "BE", "DN", "DE", "QC")) %>%
  mutate(
    ltrunc1 = ltrunc1 / 365.25,
    rtrunc1 = rtrunc1 / 365.25,
    ltrunc2 = ltrunc2 / 365.25,
    rtrunc2 = rtrunc2 / 365.25,
    nyears = ndays / 365.25
  )
for (i in seq_along(thyr)) {
  mod0 <- with(
    semi_full,
    fit_elife(
      time = nyears,
      ltrunc = cbind(ltrunc1, ltrunc2),
      rtrunc = cbind(rtrunc1, rtrunc2),
      family = "exp",
      thresh = thyr[i]
    )
  )
  mod1 <- with(
    semi_full,
    fit_elife(
      time = nyears,
      ltrunc = cbind(ltrunc1, ltrunc2),
      rtrunc = cbind(rtrunc1, rtrunc2),
      family = "gp",
      thresh = thyr[i]
    )
  )
  vfit$par[i] <- coef(mod0)
  vfit$std[i] <- mod0$std.error
  vfit$nexc[i] <- nobs(mod1)
  vfit$pval[i] <- anova(mod0, mod1)$`Pr(>Chisq)`[2]
}
