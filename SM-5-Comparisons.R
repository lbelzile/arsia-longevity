#####################################################
#####     Supplementary material section 5       ####
#####       Is there a Cap on Longevity?         ####
#####################################################
#### Data for IDL can be downloaded from         ####
#### supercentenarian.org, but the database      ####
#### is subject to change. Most data cannot      ####
#### be shared due to privacy requirements       ####
#####################################################

library(lubridate)
library(tidyverse)
library(longevity)
library(patchwork)
library(ggplot2)

library(tikzDevice)
options(tikzLatexPackages =
          c("\\usepackage{tikz}\n",
            "\\usepackage[active,tightpage,psfixbb]{preview}\n",
            "\\usepackage{amsmath}",
            "\\PreviewEnvironment{pgfpicture}\n",
            "\\setlength\\PreviewBorder{0pt}\n",
            "\\usepackage{fourier}\n"
          ))
# Threshold 110, but exact number of days can vary
# from one individual to the next because of leap years
# Pick the lowest value for which declared 110 for all 
# countries in the IDL as threshold 
u110d <- 40175L
u110y <- u110d / 365.25
# but filter by those categorized as 110
idl <- idl2021 %>%
  filter(ageyear >= 110) %>%
  mutate(
    ltrunc = pmax(u110d, as.integer(c1 - bdate)) / 365.25,
    rtrunc = as.integer(c2 - bdate) / 365.25,
    nyears = ndays / 365.25
  )
# Define regions for various tables
north_europe <- idl %>% 
  filter(country %in% c("OS", "BE", "DN", "FI", "DE", "NO", "SV", "EW"))
south_europe <- idl %>% 
  filter(country %in% c("FR", "ES"))
europe <- bind_rows(north_europe, south_europe)
north_america <- idl %>% 
  filter(country %in% c("US", "QC"))
world <- bind_rows(europe, north_america)

# Counts by gender and region
counts_gender <- rbind(
  with(north_europe, table(gender)),
  with(south_europe, table(gender)),
  with(europe, table(gender)),
  with(north_america, table(gender)),
  with(world, table(gender))
)
region <- c("North Europe", 
           "South Europe", 
           "Europe", 
           "North America",
           "world")
table_counts <- tibble(
  region = region,
  "list of countries" = 
    c("Austria, Belgium, Denmark, Finland, 
      Germany, Norway, Sweden, England and Wales",
                "France, Spain",
                "North and South Europe",
                "USA, Quebec",
                "Europe and North America"),
  women = counts_gender[,1],
  men = counts_gender[,2]
)
# Create list of databases for looping
dbs <- list(north_europe, 
            south_europe, 
            europe,
            north_america,
            world)

#---------------------------------------------------------
# Likelihood ratio tests comparing different
# generalized Pareto and different exponential models
# for women and men against the null hypothesis
# of the same distribution for women and men

lrt_gender <-
  sapply(dbs, function(db){
# Likelihood ratio test assuming GP
# with different parameters for women/men
gp_lrt <- with(
  db,
  test_elife(
    time = nyears,
    thresh = u110y,
    ltrunc = ltrunc,
    rtrunc = rtrunc,
    family = "gp",
    covariate = gender
  )
)
# Same thing, but with exponential distribution
exp_lrt <- with(
  db,
  test_elife(
    time = nyears,
    thresh = u110y,
    ltrunc = ltrunc,
    rtrunc = rtrunc,
    family = "exp",
    covariate = gender
  )
)
# Return only p-values
  return(c(gp_lrt$pval, exp_lrt$pval))
})
# Summary of likelihood ratio test
lrt_gender_tab <- 
  tibble(region = region,
       "generalized Pareto" = lrt_gender[1,],
       "exponential" = lrt_gender[2,])

#-------------------------------------------------------
# Likelihood ratio tests comparing different generalized
# Pareto and exponential distributions between groups of
# countries for IDL supercentenarian against the null of
# common distribution for both groups of countries.

test_country_group1_gp <-
  with(
    bind_rows(
        north_europe %>% mutate(group = "north"),
        south_europe %>% mutate(group = "south")
      ),
    test_elife(
      time = nyears,
      thresh = u110y,
      ltrunc = ltrunc,
      rtrunc = rtrunc,
      family = "gp",
      covariate = group
    )
  )

test_country_group1_exp <-
  with(
    bind_rows(
      north_europe %>% mutate(group = "north"),
      south_europe %>% mutate(group = "south")
    ),
    test_elife(
      time = nyears,
      thresh = u110y,
      ltrunc = ltrunc,
      rtrunc = rtrunc,
      family = "exp",
      covariate = group
    )
  )
test_country_group2_gp <-
  with(
    bind_rows(
      europe %>% mutate(group = "europe"),
      north_america %>% mutate(group = "america")
    ),
    test_elife(
      time = nyears,
      thresh = u110y,
      ltrunc = ltrunc,
      rtrunc = rtrunc,
      family = "gp",
      covariate = group
    )
  )

test_country_group2_exp <-
  with(
    bind_rows(
      europe %>% mutate(group = "europe"),
      north_america %>% mutate(group = "america")
    ),
    test_elife(
      time = nyears,
      thresh = u110y,
      ltrunc = ltrunc,
      rtrunc = rtrunc,
      family = "exp",
      covariate = group
    )
  )

table_diff_regions <- 
  tibble(region = c("North Europe vs South Europe",
                  "Europe vs North America"),
       "generalized Pareto" = c(
         test_country_group1_gp$pval,
         test_country_group2_gp$pval),
       exponential = c(
         test_country_group1_exp$pval,
         test_country_group2_exp$pval)
)

#-------------------------------------------------------
# Comparison of generalized Pareto
# versus exponential distribution
shape_lrt_exp_vs_gp <- 
  t(sapply(dbs, 
         function(db){
  fit_gp <- with(
    db,
    fit_elife(
      time = nyears,
      thresh = u110y,
      ltrunc = ltrunc,
      rtrunc = rtrunc,
      family = "gp"
    )
  )
  fit_exp <- with(
    db,
    fit_elife(
      time = nyears,
      thresh = u110y,
      ltrunc = ltrunc,
      rtrunc = rtrunc,
      family = "exp"
    )
  )
  # Return information for table
  # Point estimates for shape
  c(shape = as.numeric(fit_gp$par["shape"]),
    # Wald-based confidence intervals
    lci = as.numeric(fit_gp$par["shape"] + qnorm(0.025)*fit_gp$std.error["shape"]),
    uci = as.numeric(fit_gp$par["shape"] + qnorm(0.975)*fit_gp$std.error["shape"]),
    # Model comparison (likelihood ratio test)
    pval = anova(fit_gp, fit_exp)[2,"Pr(>Chisq)"])
}))

table_shape_expvsgp <- cbind(region = region, 
                             as_tibble(shape_lrt_exp_vs_gp))

#-------------------------------------------------------
# Analysis of data for countries other than
# France, Italy and England & Wales in the IDL
# that also include both semisupercentenarian
# and supercentenarians
countries_leftover <- c("OS","BE","DN","DE","QC")
country_leftover_comparisons <-
  t(sapply(105:110, function(th){
    gp_model <-
      with(idl2021 %>% 
           filter(country %in% countries_leftover),
         fit_elife(time = ndays / 365.25,
                   thresh = th,
                   ltrunc = cbind(ltrunc1, ltrunc2) / 365.25,
                   rtrunc = cbind(rtrunc1, rtrunc2) / 365.25,
                   family = "gp"))
    exp_model <- 
      with(idl2021 %>% 
             filter(country %in% countries_leftover),
           fit_elife(time = ndays / 365.25,
                     thresh = th,
                     ltrunc = cbind(ltrunc1, ltrunc2) / 365.25,
                     rtrunc = cbind(rtrunc1, rtrunc2) / 365.25,
                     family = "exp"))
    c(par = as.numeric(exp_model$par),
      se = as.numeric(exp_model$std.error),
      pval = anova(exp_model, gp_model)[2,"Pr(>Chisq)"],
      nexc = exp_model$nexc)
}))
tab_leftover <- as_tibble(country_leftover_comparisons) %>%
  mutate(threshold = 105:110,
         nexc = as.integer(nexc))

#-------------------------------------------------------
# Compare early/late death cohort through changepoint
# with time-varying hazard (time being calendar and 
# assumed common to every country/subnational entry)
# Modify the likelihood

breaks <- lubridate::ymd(paste0(1995:2010, "-12-31"))

Wald_hazard_ratio <- function(data, m = ymd("2001-12-31")){
  likelihood_split <- function(pars, 
                               data, 
                               m = ymd("2001-12-31"), 
                               u110d = 40175L,
                               transfo = TRUE){
    if(pars[1] < 0){
      return(1e10)
    }
    scale1 <- pars[1]
    if(transfo){
      scale2 <- exp(pars[2])*scale1
    } else{
      scale2 <- pars[2]
    }
    loglik <- with(data,
                   ifelse(ddate > m,
                          dexp(ndays - u110d, rate = 1/scale2, log = TRUE),
                          dexp(ndays - u110d, rate = 1/scale1, log = TRUE))) -
      log(pexp(q = with(data, pmax(0, as.numeric(m - bdate) - u110d)),
               rate = 1/scale1) -
            pexp(q = with(data, pmax(0, as.numeric(c1 - bdate) - u110d)),
                 rate = 1/scale1) +
            pexp(q = with(data, pmax(0, as.numeric(c2 - bdate) - u110d)),
                 rate = 1/scale2) -
            pexp(q = with(data, pmax(0, as.numeric(m - bdate) - u110d)),
                 rate = 1/scale2))
    -sum(loglik)
  }
  sep_fit <- optim(par = c(500, 0),
                   fn = likelihood_split,
                   method = "N",
                   m = m,
                   data = data,
                   hessian = TRUE)
  loghazard = sep_fit$par[2]
  n1 = sum(with(data, ddate <= m))
  n2 = sum(with(data, ddate > m))
  se = sqrt(solve(sep_fit$hessian)[2,2])
  
  scales_mle <- c(sep_fit$par[1],
                  exp(sep_fit$par[2])*sep_fit$par[1])
  # Standard errors of scale
  obs_info <- numDeriv::hessian(func = function(x){
    likelihood_split(pars = x, 
                     m = m,
                     data = data,
                     transfo = FALSE)
    },
    x = scales_mle)
  var_scales <- diag(solve(obs_info))
  wald_stat <- diff(scales_mle)^2/sum(var_scales)
  wald_pval <- pchisq(wald_stat, df = 1, lower.tail = FALSE)
  psi <- loghazard + seq(-3, 3, length.out = 101L)*se
  pll <- sapply(psi, function(psi_i){
    opt_prof <- optim(par = sep_fit$par[1], 
                      fn = function(lambda){
                        likelihood_split(pars = c(lambda, psi_i), 
                                         data = data, m = m)}, 
                      method = "Brent", 
                      lower = 200, 
                      upper = 800)
    opt_prof$value
    })
  opt_prof_zero <- optim(par = sep_fit$par[1], 
                    fn = function(lambda){
                      likelihood_split(pars = c(lambda, 0), 
                                       data = data, m = m)}, 
                    method = "Brent", 
                    lower = 200, 
                    upper = 800)
  object <- list(pll = -pll, 
                 mle = sep_fit$par[2:1],
                 maxpll = -pll[which.min(pll)],
                 psi.max = psi[which.min(pll)],
                 std.error = se,
                 psi = psi)
  pval_lrt <- pchisq(2*(object$maxpll + opt_prof_zero$value), 
                     df = 1, 
                     lower.tail = FALSE)
  object$r <- with(object, sign(psi.max-psi)*sqrt(2*(maxpll - pll)))
  class(object) <- "eprof"
  profile_confint_95 <- mev:::confint.eprof(object, 
                                         parm = "profile")
  profile_confint_90 <- mev:::confint.eprof(object, 
                                            parm = "profile",
                                            level = 0.9)
  
  wald = loghazard / se
  pval = 2*pnorm(abs(wald), lower.tail = FALSE)
  list(loghazard = loghazard,
       scale1 = sep_fit$par[1],
       scale2 = exp(loghazard)*sep_fit$par[1],
       se = se, 
       wald_scale = wald_stat,
       wald_pval = wald_pval,
       wald = wald, 
       pval = pval, 
       n1 = n1, 
       n2 = n2,
       pval_lrt = pval_lrt,
       lci90 = profile_confint_90[2],
       uci90 = profile_confint_90[3],
       lci95 = profile_confint_95[2],
       uci95 = profile_confint_95[3])
}

# Log hazards
wald_table <- cbind(
  north_europe = sapply(breaks, function(br){
    Wald_hazard_ratio(north_europe, br)
  }),
  south_europe = sapply(breaks, function(br){
    Wald_hazard_ratio(south_europe, br)
  }),
  europe = sapply(breaks, function(br){
    Wald_hazard_ratio(europe, br)
  }),
  north_america = sapply(breaks, function(br){
    Wald_hazard_ratio(north_america, br)
  }),
  world = sapply(breaks, function(br){
    Wald_hazard_ratio(world, br)
  })
)
# Without Jeanne Calment
wald_table_mJC <- cbind(
  north_europe = sapply(breaks, function(br){
    Wald_hazard_ratio(north_europe, br)
  }),
  south_europe = sapply(breaks, function(br){
    Wald_hazard_ratio(south_europe %>% filter(ndays != max(ndays)), br)
  }),
  europe = sapply(breaks, function(br){
    Wald_hazard_ratio(europe %>% filter(ndays != max(ndays)), br)
  }),
  north_america = sapply(breaks, function(br){
    Wald_hazard_ratio(north_america, br)
  }),
  world = sapply(breaks, function(br){
    Wald_hazard_ratio(world %>% filter(ndays != max(ndays)), br)
  })
)
wald_table2 <- round(matrix(unlist(wald_table),
                           ncol = nrow(wald_table), byrow = TRUE),
                    digits = 3)
colnames(wald_table2) <- rownames(wald_table)

wald_tab <- data.frame(
  group = rep(region,
              each = nrow(wald_table2)/length(region)),
  year = rep(as.integer(year(breaks)), length.out = nrow(wald_table2)),
  wald_table2)

wald_table2_mJC <- round(matrix(unlist(wald_table_mJC),
                            ncol = nrow(wald_table_mJC), byrow = TRUE),
                     digits = 3)
colnames(wald_table2_mJC) <- rownames(wald_table_mJC)

wald_tab_mJC <- data.frame(
  group = rep(region,
              each = nrow(wald_table2_mJC)/length(region)),
  year = rep(as.integer(year(breaks)), length.out = nrow(wald_table2_mJC)),
  wald_table2_mJC)
# Remove intermediate tables
rm(wald_table_mJC, wald_table2_mJC, wald_table2, wald_table)




# Create plots for PDF format (with LaTeX)
g1 <- ggplot(data = wald_tab %>% filter(group != "Europe",
                                  year <= 2007 & year > 1996),
       mapping = aes(x = year,
                     y = loghazard,
                     color = group, 
                     group = group,
                     ymin = lci95,
                     ymax = uci95)) +
  geom_hline(yintercept = 0, alpha = 0.2) +
  geom_line(position = position_dodge(width = 0.5)) + 
  geom_linerange(position = position_dodge(width = 0.5),
                size = 1.5) +
  geom_point(
    col = "white", 
    pch = 15,
    size = 1.25,
    position = position_dodge(width = 0.5)) +
  viridis::scale_color_viridis(discrete = TRUE, 
                               option = "G",
                               end = 0.8,
                               direction = -1) + 
  theme_classic() + 
  scale_x_continuous(breaks = seq(1996L, 2006L, by = 2L)) +
  scale_y_continuous(limits = c(-0.78,0.5), expand = c(0,0),
                     breaks = seq(-0.75, 0.5, by = 0.25),
                     labels = paste0("$",seq(-0.75, 0.5, by = 0.25),"$"))  +
  labs(y = "",
       subtitle = "all supercentenarians",
       color = "region") +
  theme(legend.position = "bottom",
        plot.subtitle = element_text(hjust = 1)) + 
  guides(size = "none")
g2 <- ggplot(data = wald_tab_mJC %>% filter(group != "Europe",
                                        year <= 2007 & year > 1996),
             mapping = aes(x = year,
                           y = loghazard,
                           color = group, 
                           group = group,
                           ymin = lci95,
                           ymax = uci95)) +
  geom_hline(yintercept = 0, alpha = 0.2) +
  geom_line(position = position_dodge(width = 0.5)) + 
  geom_linerange(position = position_dodge(width = 0.5),
                 size = 1.5) +
  geom_point(
    col = "white", 
    pch = 15,
    size = 1.25,
    position = position_dodge(width = 0.5)) +
  viridis::scale_color_viridis(discrete = TRUE, 
                               option = "G", 
                               end = 0.8,
                               direction = -1) + 
  theme_classic() + 
  scale_x_continuous(breaks = seq(1996L, 2006L, by = 2L)) +
  scale_y_continuous(limits = c(-0.78,0.5), expand = c(0,0),
                     breaks = seq(-0.75, 0.5, by = 0.25),
                     labels = paste0("$",seq(-0.75, 0.5, by = 0.25),"$"), 
                     minor_breaks = seq(-0.45, 0.45, by = 0.05)) +
  labs(y = "",
       subtitle = "excluding Jeanne Calment",
       color = "region") +
  theme(legend.position = "bottom",
        plot.subtitle = element_text(hjust = 1)) + 
  guides(size = "none")


# tikz(file = "SM5-loghazard.tex",
#      width = 8, 
#      height = 4.5, 
#      standAlone = TRUE)
patchwork <- g1 + g2 + plot_layout(guides = 'collect') & 
  theme(legend.position = "bottom")
patchwork + 
  plot_annotation(title = "$\\log(\\sigma_2/\\sigma_1)$")
# dev.off()

# P-values quoted for transition 2002-2003
wald_tab %>% filter(year == 2002) %>% 
  select(group, pval_lrt)

# Export tables

system("mkdir tables")
# Table SM 1 - counts
table_counts %>%
  knitr::kable(format = 'latex', 
        booktabs = TRUE,
        align = c("l","p{6.8cm}","c","c"),
        caption = "\\textsf{IDL} dataset regions used in the analysis, with sample size broken down by gender.",
        label = "groups",
        position = "htbp!") %>% 
  cat(file = "tables/Table_SM1.tex", 
      sep = "\n",
      append = FALSE)
  

# Table SM 2 - gender
lrt_gender_tab  %>%
  knitr::kable(format = 'latex', 
               booktabs = TRUE,
               align = "lcc",
               digits = 2,
               caption = "$p$-values for likelihood ratio tests comparing different generalized Pareto models (middle) and different exponential models (right) for women and men against the null hypothesis of the same distribution for women and men, for \\textsf{IDL} supercentenarian from different regions.",
               label = "women-men",
               position = "htbp!") %>% 
  cat(file = "tables/Table_SM2.tex", 
      sep = "\n",
      append = FALSE)

# Table SM 3 - death cohort (changepoint)

table_changepoint <- 
  wald_tab %>%
  mutate(ndiff = abs(n2 - n1)) %>%
  group_by(group) %>%
  filter(ndiff == min(ndiff)) %>%
  bind_cols(
    data.frame(t(sapply(dbs, function(db){
      db %>% 
        summarize(minyear = year(min(c1)), 
                  maxyear = year(max(c2)))
    })))) %>% ungroup()
  
table_changepoint %>%
  transmute(region = group,
            ddate1 = paste0(minyear, "--", 
                            year - 1),
            n1 = n1,
            scale1 = scale1/365.25,
            ddate2 = paste0(year, "--", 
                            maxyear),
            n2 = n2,
            scale2 = scale2/365.25,
            pval = wald_pval) %>%
   knitr::kable(format = 'latex', 
               booktabs = TRUE,
               align = "lccccccc",
               digits = c(0,0,0,2,0,0,2,3),
               col.names = c("region",
                             "early",
                             "$n_1$",
                             "$\\widehat{\\sigma}_1$",
                             "late",
                             "$n_2$",
                             "$\\widehat{\\sigma}_2$",
                             "$p$-value"),
               caption = "$p$-values for Wald tests of the null hypothesis $\\mathscr{H}_0: \\sigma_1 = \\sigma_2$ for different exponential hazards for early $(\\sigma_1)$ and late $(\\sigma_2)$ periods, with sample sizes and years of death for the observation scheme.",
               label = "early-late",
               position = "htbp!",
               escape = FALSE) %>% 
  cat(file = "tables/Table_SM3.tex", 
      sep = "\n",
      append = FALSE)


# Table SM 4 - regions

table_diff_regions  %>%
  knitr::kable(format = 'latex', 
               booktabs = TRUE,
               align = "lcc",
               digits = c(0,2,2),
               caption = "$p$-values for the likelihood ratio tests comparing different generalized Pareto (middle) and exponential (right) distributions between groups of countries for \\textsf{IDL} supercentenarians against the null of common distribution for both groups of countries.",
               label = "countries",
               position = "htbp!") %>% 
  cat(file = "tables/Table_SM4.tex", 
      sep = "\n",
      append = FALSE)

# Table SM 5 - exponential vs GP by threshold

table_shape_expvsgp %>%
  transmute(region = region, 
            shape = paste0("$", 
                           ifelse(shape < 0, "\\hphantom{-}",""),
                           formatC(shape, digits = 2, format = "f"),
                           " (",
                           formatC(lci, digits = 2, format = "f"),
                           ", ",
                           formatC(uci, digits = 2, format = "f"),
                           ")$"),
            pval = pval) %>%
  knitr::kable(format = 'latex', 
               booktabs = TRUE,
               align = "lrc",
               digits = c(0,2,2),
               escape = FALSE,
               col.names = c("region", "shape (95\\% Wald CI)","$p$-value"),
               caption = "Estimates (95\\% Wald-based confidence intervals) of the shape parameter of the generalized Pareto distribution, with $p$-values for likelihood ratio test of the null hypothesis of exponential distribution ($\\xi = 0$).",
               label = "test-exp",
               position = "htbp!") %>% 
  cat(file = "tables/Table_SM5.tex", 
      sep = "\n",
      append = FALSE)
  
  

# Table SM 6 - other countries

cbind(rownames = c("$\\sigma$","$p$-value","$n_u$"),
      tab_leftover %>%
  transmute(sigma = paste0("$", 
                           formatC(par, 
                                   digits = 2, 
                                   format = "f"),
                           "\\ (", 
                           formatC(se, 
                                   digits = 2,
                                   format = "f"),")$"),
            pval = paste0("$", formatC(pval,digits = 2), "$"),
            nexc = paste0("$", formatC(nexc, format = "d"),"$")) %>%
  t()) %>%
  knitr::kable(format = 'latex', 
               booktabs = TRUE,
               align = "lcccccc",
               row.names = FALSE,
               col.names = c("threshold",105:110),
               escape = FALSE,
               caption = "Estimates (standard errors) of the scale parameter $\\sigma$ of fitted exponential distributions and $p$-value for likelihood ratio tests of the exponential model as a submodel of the generalized Pareto model, for the \\textsf{IDL} 2021 semi- and supercentenarian datasets for Austria, Belgium, Denmark, Germany and Quebec, as a function of threshold, with number of threshold exceedances ($n_u$).",
               label = "thresholds",
               position = "htbp!") %>% 
  cat(file = "tables/Table_SM6.tex", 
      sep = "\n",
      append = FALSE)
  