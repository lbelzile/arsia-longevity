setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
code_dir <- getwd()
fig_dir <- paste0(dirname(getwd()), "/figures")
library(tikzDevice)
options(tikzLatexPackages =
          c("\\usepackage{tikz}\n",
            "\\usepackage[active,tightpage,psfixbb]{preview}\n",
            "\\usepackage{amsmath}",
            "\\PreviewEnvironment{pgfpicture}\n",
            "\\setlength\\PreviewBorder{0pt}\n",
            "\\usepackage{fourier}\n",
            "\\DeclareMathAlphabet{\\mathdis}{OT1}{pag}{m}{n}\n"
          )
)
setTikzDefaults(overwrite = FALSE)

set.seed(12345)
c2 <- B <- 20 # finishing date and number of birth cohorts
rate <- 150 # mean annual number of exceedances

trunc.fit <- function(y, U)
{
  nlogL <-
    function(mu, y, U)
      - sum(dexp(y, rate = 1 / mu, log = T) - pexp(U, rate = 1 / mu, log.p =
                                                     TRUE))
  fit <- optimize(nlogL,
                  interval = c(0.01, 2),
                  y = y,
                  U = U)
  fit$minimum
}

trunccens.fit <- function(y, U)
{
  nlogL <-
    function(mu, y, U)
      - sum(ifelse(
        y < U,
        dexp(y, rate = 1 / mu, log = T) - pexp(U, rate = 1 / mu, log.p = TRUE),
        -U / mu
      ))
  fit <- optimize(nlogL,
                  interval = c(0.01, 2),
                  y = y,
                  U = U)
  fit$minimum
}


R <- 1000
colnames <- c("Est 1", "Est 2", "Est 3", "Est 4")
res <- matrix(
  NA,
  ncol = 4,
  nrow = R,
  dimnames = list(NULL, colnames)
)
for (r in 1:R)
{
  N <- rpois(n = 1, lambda = B * rate)  # number of individuals
  x <- sort(runif(N) * B)           # entry times
  b <- 1 + floor(x)               # cohorts
  t <- rexp(n = N) / log(2)          # excess lifetimes

  b.max <- tapply(t + x, b, max)  #
  b.in <- c(1:B)
  b.in <- min(b.in[(b.max > c2)]) - 1
  keep <- (b <= b.in)

  res[r, 1] <- mean(t[keep]) # incorrect fit
  res[r, 2] <- trunc.fit(y = t[keep], U = c2 - x[keep])

  keep <- (x + t < c2)
  res[r, 3] <- trunc.fit(y = t[keep], U = c2 - x[keep])
  res[r, 4] <- trunccens.fit(y = t, U = c2 - x)

}

par(mfrow = c(1, 2),
    pty = "s",
    bty = "l")
boxplot(res[,1:3],
        panel.first = {
          abline(h = 1 / log(2), col = "grey")
        })
# abline(h=1/log(2), col="grey")
# Relative bias
100 * (apply(res, 2, mean) - 1 / log(2)) / (1 / log(2))
# Standard deviation
apply(res, 2, sd)

(apply(res, 2, mean) - 1 / log(2)) * sqrt(R) / apply(res, 2, sd)



library(SMPracticals)
qqexp(
  t[keep],
  panel.first = abline(0, 1 / log(2), col = "grey"),
  pch = 16,
  cex = 0.5,
  ylim = c(0, 10),
  xlim = c(0, 10),
  xaxs = "i",
  yaxs = "i"
)



library(tidyverse)
library(patchwork)

sim_df <-
  res[,1:3] %>%
  as_tibble %>%
  pivot_longer(cols = 1:3,
               names_to = "estimator",
               values_to = "value"
               )
g1 <- ggplot(data = sim_df, aes(x = factor(estimator, labels = c("naive","extinct only","full")),
                          y = value)) +
  geom_hline(yintercept = 1/log(2)) +
  geom_boxplot(alpha = 0.9) +
  theme_classic() +
  labs(x = "estimator",
       y = "scale $\\sigma$")

g2 <- ggplot(tibble(x = qexp(ppoints(sum(keep)),rate = log(2)),
                    y = sort(t[keep])),
             aes(x = x, y = y)) +
  geom_abline(intercept = 0, slope = 1, col = "grey") +
  geom_point(pch = 16, cex = 0.8) +
  scale_y_continuous(limits = c(0,10), expand = c(0,0)) +
  scale_x_continuous(limits = c(0,10), expand = c(0,0)) +
  theme_classic() +
  labs(x = "exponential plotting positions",
       y = "ordered sample")

setwd(fig_dir)
tikz("Extinct_bias.tex", width = 8, height = 4, standAlone = TRUE)
g1 + g2
dev.off()
