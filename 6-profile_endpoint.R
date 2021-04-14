#####################################################
#####  Section 6 - Is there a Cap on Longevity?  ####
#####################################################
#### Data for IDL can be downloaded from         ####
#### supercentenarian.org, but the database      ####
#### is subject to change. Most data cannot      ####
#### be shared due to privacy requirements       ####                                ####
#####################################################


# Compute the profile log likelihood of the GP endpoint
# for various countries and pool them, akin to Fig 1
# in Davison's discussion of Rootzén and Zholud.

library(tidyverse) # data manipulation
library(lubridate) # dates
library(longevity) # model excess lifetime
library(ggplot2)   # grammar of graphics
theme_set(theme_bw() +
            theme(axis.line = element_line(colour = "black"),
                  panel.grid.minor = element_blank(),
                  panel.border = element_blank(),
                  panel.background = element_blank()))
# For each data, we transform them to yearly scale
# Dutch data
data(dutch, package = "longevity")
yr_samp <- year(attr(x = dutch, which = "sampling_frame"))
dutch1 <- dutch %>%
  subset(!is.na(ndays)) %>%
  mutate(time = ndays / 365.25,
         time2 = time,
         ltrunc = ltrunc / 365.25,
         rtrunc = rtrunc / 365.25,
         event = 1) %>%
  subset(time > 98) %>%
  select(time, time2, ltrunc, rtrunc, event, gender)
dutch2 <- dutch %>%
  subset(is.na(ndays)) %>%
  mutate(time2 = ceiling_date(dmy(paste("01-", dmonth, "-", dyear)), unit = "month") - 1 -
           dmy(paste("01-01-",byear)),
         time = dmy(paste("01-", dmonth, "-", dyear)) - dmy(paste("31-12-", byear)),
         ltrunc = dmy(paste("01-01-1986")) - dmy(paste("31-12-", byear)),
         rtrunc = dmy(paste("31-12-2015")) - dmy(paste("01-01-", byear))
  ) %>%
  select(time, time2, ltrunc, rtrunc, gender) %>%
  mutate(time = as.numeric(time) / 365.25,
         time2 = as.numeric(time2) / 365.25,
         ltrunc = as.numeric(ltrunc) / 365.25,
         rtrunc = as.numeric(rtrunc) / 365.25,
         event = 3) %>%
  subset(time > 98)
ne_dat <- rbind(dutch1, dutch2)

# Japanese data
data(japanese, package = "longevity")
jp_dat <- with(japanese,
                  data.frame(
                    time = age,
                    time2 = age + 364/365,
                    rtrunc = ifelse(gender == "female", 2015, 2013) - byear,
                    gender =  gender,
                    weight = count,
                    event = 3
                  )
)
tdeath <- as.numeric(c(ymd("2015-04-01") - ymd("1898-03-05"),
                       ymd("2013-06-12") - ymd("1897-04-19")))/365.25
jp_dat <- rbind(jp_dat,
                   data.frame(time = tdeath,
                              time2 = tdeath,
                              rtrunc = Inf,
                              gender = c("female","male"),
                              weight = 1L,
                              event = 1))
# England and Wales data
ew_dat <- idl %>%
  filter(country == "EW") %>%
  mutate(nyears = ndays / 365.25,
         ltrunc1 = ltrunc1 / 365.25,
         rtrunc1 = rtrunc1 / 365.25,
         ltrunc2 = ltrunc2 / 365.25,
         rtrunc2 = rtrunc2 / 365.25)

# French data
fr_dat <- idl %>%
  filter(country == "FR") %>%
  mutate(nyears = ndays / 365.25,
         ltrunc1 = ltrunc1 / 365.25,
         rtrunc1 = rtrunc1 / 365.25,
         ltrunc2 = ltrunc2 / 365.25,
         rtrunc2 = rtrunc2 / 365.25)
# European data
europe_dat <- idl %>%
  filter(country %in% c("DE","NO","FI","CH","SV","OS","BE","DN","ES")) %>%
  mutate(nyears = ndays / 365.25,
         ltrunc1 = ltrunc1 / 365.25,
         rtrunc1 = rtrunc1 / 365.25,
         ltrunc2 = ltrunc2 / 365.25,
         rtrunc2 = rtrunc2 / 365.25)

idl_dat <- idl %>%
  mutate(nyears = ndays / 365.25,
         ltrunc1 = ltrunc1 / 365.25,
         rtrunc1 = rtrunc1 / 365.25,
         ltrunc2 = ltrunc2 / 365.25,
         rtrunc2 = rtrunc2 / 365.25)

# Italian data
it_dat <- italian %>%
  mutate(nyears = ndays / 365.25,
         ltrunc = ltrunc / 365.25)
rm(dutch, dutch1, dutch2, japanese, tdeath, yr_samp)


# Profile log likelihood for the endpoint
psi <- seq(115, 180, by = 0.5)
thresh <- c(105, 108L, 110L)
prof_df <- list()
nexc_df <- list()
for(i in seq_along(thresh)){
prof_jp_endpt <- with(jp_dat,
       prof_gp_endpt(time = time,
                     time2 = time2,
                     event = event,
                     weights = weight,
                     rtrunc = rtrunc,
                     type = "interval2",
                     psi = psi[psi > max(time2)],
                     thresh = thresh[i],
                     confint = FALSE)
  )
prof_ne_endpt <- with(ne_dat,
                      prof_gp_endpt(time = time,
                                    time2 = time2,
                                    event = event,
                                    ltrunc = ltrunc,
                                    rtrunc = rtrunc,
                                    type = "interval2",
                                    psi = psi[psi > max(time)],
                                    thresh = thresh[i],
                                    confint = FALSE)
)
prof_ew_endpt <- with(ew_dat,
                      prof_gp_endpt(time = nyears,
                                            ltrunc = cbind(ltrunc1, ltrunc2),
                                            rtrunc = cbind(rtrunc1, rtrunc2),
                                    psi = psi[psi > max(nyears)],
                                    thresh = (38350+ c(0, 365*3, 365*5+1))[i]/365.25,#thresh[i],
                                    confint = FALSE)
)
prof_europe_endpt <- with(europe_dat,
                      prof_gp_endpt(time = nyears,
                                            ltrunc = cbind(ltrunc1, ltrunc2),
                                            rtrunc = cbind(rtrunc1, rtrunc2),
                                            psi = psi[psi > max(nyears)],
                                            thresh = (38350+ c(0, 365*3, 365*5+1))[i]/365.25,#thresh[i],
                                            confint = FALSE)
)
prof_idl_endpt <- with(idl_dat,
                          prof_gp_endpt(time = nyears,
                                                ltrunc = cbind(ltrunc1, ltrunc2),
                                                rtrunc = cbind(rtrunc1, rtrunc2),
                                                psi = psi[psi > max(nyears)],
                                                thresh = (38350+ c(0, 365*3, 365*5+1))[i]/365.25,#thresh[i],
                                                confint = FALSE)
)
prof_it_endpt <- with(it_dat,
                      prof_gp_endpt(time = nyears,
                                    event = event,
                                    ltrunc = ltrunc,
                                    rtrunc = NULL,
                                    type = "right",
                                    psi = psi[psi > max(nyears)],
                                    thresh = thresh[i],
                                    confint = FALSE)
)
prof_fr_endpt <- with(fr_dat,
                      prof_gp_endpt(time = nyears,
                                            ltrunc = cbind(ltrunc1, ltrunc2),
                                            rtrunc = cbind(rtrunc1, rtrunc2),
                                    psi = psi[psi > max(nyears)],
                                    thresh = (38350+ c(0, 365*3, 365*5+1))[i]/365.25,#thresh[i],
                                    confint = FALSE)
)
# The combined log likelihood include non-overlapping samples
# We find the maximum by simply looking on the grid of psi values
pll_full_endpt <-
  with(prof_it_endpt, pll[psi > max(fr_dat$nyears)]) +
  with(prof_idl_endpt, pll[psi > max(fr_dat$nyears)]) +
  with(prof_jp_endpt, pll[psi > max(fr_dat$nyears)]) +
  with(prof_ne_endpt, pll[psi > max(fr_dat$nyears)])

# (Messy) code to create a database with the output
# for every curve, every country.
prof_df[[i]] <- tibble(
  threshold = factor(thresh[[i]]),
  country = factor(
    c(rep("Italy", length(prof_it_endpt$psi)),
      rep("England and Wales",length(prof_ew_endpt$psi)),
      rep("Netherlands",length(prof_ne_endpt$psi)),
      rep("Japan",length(prof_jp_endpt$psi)),
      rep("France", each = length(prof_fr_endpt$psi)),
      rep("rest of Europe", each = length(prof_europe_endpt$psi)),
      rep(c("IDL","combined"), each = length(prof_idl_endpt$psi)))),
    endpt = c(prof_it_endpt$psi,
              prof_ew_endpt$psi,
              prof_ne_endpt$psi,
              prof_jp_endpt$psi,
              prof_fr_endpt$psi,
              prof_europe_endpt$psi,
              prof_idl_endpt$psi,
              prof_idl_endpt$psi),
  pll = - c(with(prof_it_endpt, maxpll - pll),
          with(prof_ew_endpt, maxpll - pll),
          with(prof_ne_endpt, maxpll - pll),
          with(prof_jp_endpt, maxpll - pll),
          with(prof_fr_endpt, maxpll - pll),
          with(prof_europe_endpt, maxpll - pll),
          with(prof_idl_endpt, maxpll - pll),
          max(pll_full_endpt) - pll_full_endpt)
  )
# Extract the number of exceedances for each threshold
nexc_v <- c(prof_it_endpt$nexc,
            prof_ew_endpt$nexc,
            prof_ne_endpt$nexc,
            prof_jp_endpt$nexc,
            prof_fr_endpt$nexc,
            prof_europe_endpt$nexc,
            prof_idl_endpt$nexc)
nexc_v <- c(nexc_v, prof_jp_endpt$nexc + prof_idl_endpt$nexc + prof_it_endpt$nexc + prof_ne_endpt$nexc)
nexc_df[[i]] <- tibble(threshold = factor(thresh[[i]]),
                       country = factor(c("Italy",
                                          "England and Wales",
                                          "Netherlands",
                                          "Japan",
                                          "France",
                                          "rest of Europe",
                                          "IDL",
                                          "combined")),
                      nexc = nexc_v)

}
# Unlist the data and merge
prof_df <- bind_rows(prof_df)
nexc_df <- bind_rows(nexc_df)
labels_df <- tibble(country = factor(c("Italy",
                "England and Wales",
                "Netherlands",
                "Japan",
                "France",
                "rest of Europe",
                "IDL",
                "combined")),
  nexc = with(nexc_df, paste(nexc[threshold == "105"],
                             nexc[threshold == "108"],
                             nexc[threshold == "110"],
                             sep = "/"))
  )


g1 <- ggplot(data = prof_df,
             mapping = aes(x = endpt,
                           y = pll)) +
  geom_hline(yintercept = -qchisq(c(0.95,0.99), df = 1)/2, alpha = 0.1) +
  geom_line(aes(linetype = threshold)) +
  scale_linetype_manual(values = 3:1)  +
  labs(x = "lifespan (years)",
       y = "profile log-likelihood") +
  scale_y_continuous(
    breaks = round(seq(-3, 0, by = 1)),
    labels = paste0("$", -3:0,"$"),
    minor_breaks = NULL,
    limits = c(-3.5, 0),
    oob = scales::squish_infinite) +
  scale_x_continuous(
    breaks = seq(120L, 180L, by = 10L),
    minor_breaks = NULL,
    limits = c(115,180),
    expand = c(0,0)) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(color = "#F2F2F2"),
        legend.position = c(0.75,0.2)) +
  facet_wrap(~country, scales = 'free_x') +
  geom_text(
    data    = labels_df,
    hjust   = 1,
    vjust   = -0.1,
    size = 3,
    mapping = aes(x = Inf, y = -Inf, label = nexc)
  )

g1
