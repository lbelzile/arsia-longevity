#################################################
####   Preprocessing of the IDL database   ######
#################################################
# This script shows how to transform the IDL
# comma separated values file "idlcomplete.csv"
# and turn it into a data frame that is suitable
# for other analyses using the R package longevity.

library(lubridate)
library(tidyverse)

# Sampling frame (extracted from the metadata) for the IDL
# at the time of writing the paper
sampling_frame <- matrix(
  ncol = 4L,
  byrow = TRUE,
  data = c(
    "AUT", "110+", "01-01-2005", "31-12-2012",
    "AUT", "105-109", "01-01-2003", "31-12-2014",
    "BEL", "110+", "01-01-1990", "31-12-2015",
    "BEL", "105-109", "01-01-1977", "31-12-2015",
    "CAN", "110+", "01-01-1983", "31-12-2009",
    "CAN", "105-109", "01-01-1985", "31-12-2009",
    "DNK", "110+", "01-01-1996", "31-12-2014",
    "DNK", "105-109", "01-01-1970", "31-12-2014",
    "FIN", "110+", "01-01-1989", "31-12-2006",
    "DEU", "110+", "01-01-1994", "31-12-2005",
    "DEU", "105-109", "01-01-1989", "31-12-2005",
    "NOR", "110+", "01-01-1989", "31-12-2004",
    "NOR", "105-109", "01-01-1986", "31-12-2006",
    "SWE", "110+", "01-01-1986", "31-12-2003",
    "ESP", "110+", "01-01-1989", "31-12-2016",
    "CHE", "110+", "01-01-1993", "31-12-2000",
    "CHE", "105-109", "01-01-1987", "31-12-2005", # 1971
    "ITA", "110+", "01-01-1973", "31-12-2003",
    "FRA", "105-109", "01-01-1978", "31-12-2017",
    "FRA", "110+", "01-01-1987", "31-12-2017",
    "USA", "110+", "01-01-1980", "31-12-2010",
    "UK", "110+", "01-01-1968", "31-12-2017",
    "UK", "105-109", "01-01-2000", "31-12-2014"
  )
) %>%
  as_tibble(.name_repair = "universal") %>%
  rename(
    country = ...1,
    group = ...2,
    ldate = ...3,
    rdate = ...4
  ) %>%
  mutate(
    country = factor(country),
    group = factor(group),
    ldate = dmy(ldate),
    rdate = dmy(rdate)
  )
# Group the countries to obtain the bounds
# [c1, c2] for 110+ and [d1, d2] for [105-110)
sf_wide <- tidyr::pivot_wider(sampling_frame,
  names_from = "group",
  values_from = c("ldate", "rdate")
) %>%
  rename(
    c1 = "ldate_110+",
    d1 = "ldate_105-109",
    c2 = "rdate_110+",
    d2 = "rdate_105-109"
  )

# 1) Read the IDL data
# 2) Select relevant columns
# 3) Rename columns
# 4) Extract country ID from the IDNumber (first three characters)
# 5) Remove Japan (Jdanov, 2021) and USA which only includes birth/death year
# 6) Cast columns to factor or Date, whenever relevant
idl <- read.csv(file = "idl_complete.csv") %>%
  as_tibble() %>%
  select(IDNUMBER, AGEDAYS, AGEYEARS, SEX, BDATE, DDATE, UPDATE, VALIDATION) %>%
  rename(
    ndays = AGEDAYS,
    gender = SEX,
    bdate = BDATE,
    ddate = DDATE,
    wave = UPDATE,
    ageyear = AGEYEARS,
    valid = VALIDATION
  ) %>%
  mutate(country = factor(substr(IDNUMBER, 1, ifelse(substr(IDNUMBER, 1, 2) == "UK", 2, 3)))) %>%
  # Create country manually b/c one record is . (AUSTRIA)
  select(!IDNUMBER) %>%
  filter(!country %in% c("JPN", "USA")) %>%
  # USA doesn't record dates, only years
  mutate(
    bdate = dmy(bdate),
    ddate = dmy(ddate),
    valid = factor(valid),
    gender = factor(gender)
  )

# Select only USA and repeat the treatment
idl_USA <- read.csv(file = "idl_complete.csv") %>%
  as_tibble() %>%
  select(IDNUMBER, AGEDAYS, AGEYEARS, SEX, BDATE, DDATE, UPDATE, VALIDATION) %>%
  rename(
    ndays = AGEDAYS,
    gender = SEX,
    byear = BDATE,
    dyear = DDATE,
    wave = UPDATE,
    ageyear = AGEYEARS,
    valid = VALIDATION
  ) %>%
  mutate(country = factor(substr(IDNUMBER, 1, ifelse(substr(IDNUMBER, 1, 2) == "UK", 2, 3)))) %>%
  # Create country manually b/c one record is . (AUSTRIA)
  select(!IDNUMBER) %>%
  filter(country %in% "USA") %>%
  # USA doesn't record dates, only years
  mutate(
    bdate = pmax(
      dmy(paste0("01-01-", byear)),
      dmy(paste("01-01-", dyear)) - ndays
    ),
    ddate = pmin(
      dmy(paste0("31-12-", dyear)),
      dmy(paste("31-12-", byear)) + ndays
    ),
    valid = factor(valid),
    gender = factor(gender)
  ) %>%
  filter(!ageyear < 110) %>%
  select(!byear, dyear)

# 1) Combine the database
# 2) Exclude observations that are not validated
# 3) Create an indicator for semisupercentenarian and supercentenarian
# 4) Join with the metadata (sampling frame) per country
# 5) Calculate the date at which the individual reaches 105 or 110 (handling leap years)
# 6) Group in by status (semi- or super)
# 7) Change country code to two characters
idl <- bind_rows(idl, idl_USA) %>%
  select(!valid) %>%
  # all validated (flag)
  mutate(
    country = factor(country, exclude = NULL),
    group = factor(ifelse(ageyear < 110, "105-109", "110+"))
  ) %>%
  left_join(y = sf_wide, by = c("country")) %>%
  mutate(
    x110 = if_else(mday(bdate) == 29 & month(bdate) == 2,
      bdate + days(1) + years(110), # leap year and people born on February 29th
      bdate + years(110)
    ),
    x105 = if_else(mday(bdate) == 29 & month(bdate) == 2,
      bdate + days(1) + years(105),
      bdate + years(105)
    )
  ) %>%
  select(!dyear, group) %>%
  filter(if_else(ageyear >= 110,
    ddate >= c1 & ddate <= c2,
    ddate >= d1 & ddate <= d2
  )) %>%
  arrange(country, ndays) %>%
  mutate(gender = factor(ifelse(gender == "F", "female", "male"))) %>%
  mutate(country = fct_recode(country,
    "CH" = "CHE",
    "OS" = "AUT",
    "BE" = "BEL",
    "QC" = "CAN",
    "DE" = "DEU",
    "DN" = "DNK",
    "ES" = "ESP",
    "FI" = "FIN",
    "FR" = "FRA",
    "NO" = "NOR",
    "SV" = "SWE",
    "EW" = "UK",
    "IT" = "ITA",
    "US" = "USA"
  ))
# Pivot the table and create columns with left and right truncation
# potentially double because of interval truncation
idl$ltrunc1 <- NA
idl$rtrunc1 <- NA
idl$ltrunc2 <- NA
idl$rtrunc2 <- NA
# This for loop is inefficient, but the code is more digestable
for (i in seq_along(idl$ndays)) {
  if (is.na(idl$d1[i])) { # No semisupercentenarian
    idl$ltrunc1[i] <- max(idl$x110[i], idl$c1[i]) - as.integer(idl$bdate[i])
  } else if (idl$c1[i] <= idl$d1[i] & idl$d1[i] < idl$x110[i]) {
    idl$ltrunc1[i] <- max(idl$x105[i], idl$d1[i]) - as.integer(idl$bdate[i])
  } else if (idl$c1[i] <= idl$d1[i] & idl$x110[i] <= idl$d1[i]) {
    idl$ltrunc1[i] <- max(idl$x110[i], idl$c1[i]) - as.integer(idl$bdate[i])
  } else if (idl$d1[i] < idl$c1[i] & idl$x110[i] > idl$d1[i]) {
    idl$ltrunc1[i] <- max(idl$x105[i], idl$d1[i]) - as.integer(idl$bdate[i])
  } else if (idl$d1[i] < idl$c1[i] & idl$x110[i] <= idl$d1[i]) {
    idl$ltrunc1[i] <- as.integer(idl$c1[i] - idl$bdate[i])
  }
  # Right truncation
  if (is.na(idl$d2[i])) {
    idl$rtrunc1[i] <- as.integer(idl$c2[i] - idl$bdate[i])
  } else if (idl$x110[i] <= idl$c2[i]) {
    idl$rtrunc1[i] <- as.integer(idl$c2[i] - idl$bdate[i])
  } else if (idl$c2[i] < idl$x110[i] & idl$x110[i] <= idl$d2[i]) {
    idl$rtrunc1[i] <- as.integer(idl$x110[i] - idl$bdate[i])
  } else if (idl$x110[i] > idl$d2[i]) {
    idl$rtrunc1[i] <- as.integer(idl$d2[i] - idl$bdate[i])
  }
  # Double truncation
  if (!is.na(idl$d1[i]) && !is.na(idl$d2[i])) {
    if (idl$d1[i] < idl$c1[i] & idl$x110[i] > idl$d1[i]) {
      if (idl$x110[i] < idl$c1[i]) { # double interval truncation
        idl$ltrunc2[i] <- as.integer(idl$c1[i] - idl$bdate[i])
        idl$rtrunc1[i] <- as.integer(idl$x110[i] - idl$bdate[i])
        idl$rtrunc2[i] <- as.integer(idl$c2[i] - idl$bdate[i])
      }
    } else if (idl$d2[i] < idl$x110[i] & idl$x110[i] < idl$c2[i]) {
      idl$rtrunc1[i] <- as.integer(idl$d2[i] - idl$bdate[i])
      idl$ltrunc2[i] <- as.integer(idl$x110[i] - idl$bdate[i])
      idl$rtrunc2[i] <- as.integer(idl$c2[i] - idl$bdate[i])
    }
  }
}
idl <- idl %>% mutate(
  ltrunc1 = as.integer(ltrunc1),
  rtrunc1 = as.integer(rtrunc1),
  ltrunc2 = as.integer(ltrunc2),
  rtrunc2 = as.integer(rtrunc2)
)
# Fix this operation for the US, which doesn't have known
# birth date and death date, only year
idl <- idl %>%
  mutate(
    ltrunc1 = if_else(country == "US", pmax(40117L, as.integer(c1 - bdate)), ltrunc1),
    rtrunc1 = if_else(country == "US", as.integer(c2 - bdate), rtrunc1)
  )

usethis::use_data(idl, overwrite = TRUE)
