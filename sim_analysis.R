## Practical
rm(list = ls())
setwd("")
library(tidyverse)
library(haven)
library(nlme)
library(xtable)

## Original data analysis
d <- read_sas("exam9_1.sas7bdat")
d$SUBJECT <- as.factor(d$SUBJECT)
d$REGIMEN <- relevel(as.factor(d$REGIMEN), ref = "F")
d$PERIOD <- as.factor(d$PERIOD)
d <- mutate(d, baseline = rep(QTCF[TIME == 0], each = 6)) %>%
  arrange(SUBJECT, PERIOD, TIME) %>%
  filter(TIME != 0, !(SUBJECT %in% c(206, 229))) %>%
  mutate(PERIOD = as.factor(rep(rep(1:4, 5), 39)))
d <- group_by(d, SUBJECT) %>% mutate(meanbaseline = rep(mean(baseline), 20))
d$timeF <- factor(d$TIME, labels = c("0.5", "1.0", "1.5", "2.5", "4.0"))
m <- lme(QTCF ~ baseline : timeF + meanbaseline : timeF + PERIOD * timeF + timeF : REGIMEN - 1,
         random = ~ 1 | SUBJECT / PERIOD,
         correlation = corSymm(form = ~ as.integer(timeF) | SUBJECT / PERIOD), weights = varIdent(form=~ 1 | TIME),
         data = d)
beta <- m$coefficients$fixed


## Analysis
simresults <- readRDS("simresults.rds")  # Read results
## Average estimates
ests <- group_by(simresults, effect) %>% # Average estimates
    summarise(m1_UN = mean(est_m1_UN), m1_AR1 = mean(est_m1_AR1), m1_iid = mean(est_m1_iid),
              m2_UN = mean(est_m2_UN), m2_AR1 = mean(est_m2_AR1), m2_iid = mean(est_m2_iid),
              m3_UN = mean(est_m3_UN), m3_AR1 = mean(est_m3_AR1), m3_iid = mean(est_m3_iid),
              m4_UN = mean(est_m4_UN), m4_AR1 = mean(est_m4_AR1), m4_iid = mean(est_m4_iid),
              nonparametric = mean(est_nonparametric))
ests <- arrange(ests, substr(ests$effect,17,17), substr(ests$effect,6,8)) # Sort with treatment first, and then time point
for(i in 1:nrow(ests)) ests[i,2:ncol(ests)] <- ests[i,2:ncol(ests)] - beta[30+i] # Turn average estimates into bias
round(as.matrix(ests[,-1]),2) # Results for bias
print(xtable(as.matrix(ests[,-1])), include.rownames = FALSE)

## Standard deviation of estimates
sds <- group_by(simresults, effect) %>% # Standard deviation of estimates
    summarise(m1_UN = sd(est_m1_UN), m1_AR1 = sd(est_m1_AR1), m1_iid = sd(est_m1_iid),
              m2_UN = sd(est_m2_UN), m2_AR1 = sd(est_m2_AR1), m2_iid = sd(est_m2_iid),
              m3_UN = sd(est_m3_UN), m3_AR1 = sd(est_m3_AR1), m3_iid = sd(est_m3_iid),
              m4_UN = sd(est_m4_UN), m4_AR1 = sd(est_m4_AR1), m4_iid = sd(est_m4_iid),
              nonparametric = sd(est_nonparametric))
sds <- arrange(sds, substr(sds$effect,17,17), substr(sds$effect,6,8)) # Sort with treatment first, and then time point
print(xtable(as.matrix(sds[,-1])), include.rownames = FALSE)

## Average standard errors
ses <- group_by(simresults, effect) %>% # Average estimates
    summarise(m1_UN = mean(se_m1_UN), m1_AR1 = mean(se_m1_AR1), m1_iid = mean(se_m1_iid),
              m2_UN = mean(se_m2_UN), m2_AR1 = mean(se_m2_AR1), m2_iid = mean(se_m2_iid),
              m3_UN = mean(se_m3_UN), m3_AR1 = mean(se_m3_AR1), m3_iid = mean(se_m3_iid),
              m4_UN = mean(se_m4_UN), m4_AR1 = mean(se_m4_AR1), m4_iid = mean(se_m4_iid),
              nonparametric = mean(se_nonparametric))
ses <- arrange(ses, substr(ses$effect,17,17), substr(ses$effect,6,8)) # Sort with treatment first, and then time point
print(xtable(as.matrix(ses[,-1])), include.rownames = FALSE)

## Coverages
coverages <- group_by(simresults, effect) %>% 
    summarise(m1_UN = mean((est_m1_UN - 1.96 * se_m1_UN) < beta & (est_m1_UN + 1.96 * se_m1_UN) > beta),
              m1_AR1 = mean((est_m1_AR1 - 1.96 * se_m1_AR1) < beta & (est_m1_AR1 + 1.96 * se_m1_AR1) > beta),
              m1_iid = mean((est_m1_iid - 1.96 * se_m1_iid) < beta & (est_m1_iid + 1.96 * se_m1_iid) > beta),
              m2_UN = mean((est_m2_UN - 1.96 * se_m2_UN) < beta & (est_m2_UN + 1.96 * se_m2_UN) > beta),
              m2_AR1 = mean((est_m2_AR1 - 1.96 * se_m2_AR1) < beta & (est_m2_AR1 + 1.96 * se_m2_AR1) > beta),
              m2_iid = mean((est_m2_iid - 1.96 * se_m2_iid) < beta & (est_m2_iid + 1.96 * se_m2_iid) > beta),
              m3_UN = mean((est_m3_UN - 1.96 * se_m3_UN) < beta & (est_m3_UN + 1.96 * se_m3_UN) > beta),
              m3_AR1 = mean((est_m3_AR1 - 1.96 * se_m3_AR1) < beta & (est_m3_AR1 + 1.96 * se_m3_AR1) > beta),
              m3_iid = mean((est_m3_iid - 1.96 * se_m3_iid) < beta & (est_m3_iid + 1.96 * se_m3_iid) > beta))
covarages <- arrange(coverages, substr(coverages$effect,17,17), substr(coverages$effect,6,8)) # Sort with treatment first, and then time point
print(xtable(as.matrix(coverages[,-1]), digits = 3), include.rownames = FALSE)
print(xtable(as.matrix(coverages[,-1]), digits = 2), include.rownames = FALSE)

## Adjusted coverages
k <- qt(.975, df = 38)
coverages2 <- group_by(simresults, effect) %>% 
    summarise(m1_UN = mean((est_m1_UN - k * se_m1_UN * sqrt(39/38)) < beta & (est_m1_UN + k * se_m1_UN * sqrt(39/38)) > beta),
              m1_AR1 = mean((est_m1_AR1 - k * se_m1_AR1 * sqrt(39/38)) < beta & (est_m1_AR1 + k * se_m1_AR1 * sqrt(39/38)) > beta),
              m1_iid = mean((est_m1_iid - k * se_m1_iid * sqrt(39/38)) < beta & (est_m1_iid + k * se_m1_iid * sqrt(39/38)) > beta),
              m2_UN = mean((est_m2_UN - k * se_m2_UN * sqrt(39/38)) < beta & (est_m2_UN + k * se_m2_UN * sqrt(39/38)) > beta),
              m2_AR1 = mean((est_m2_AR1 - k * se_m2_AR1 * sqrt(39/38)) < beta & (est_m2_AR1 + k * se_m2_AR1 * sqrt(39/38)) > beta),
              m2_iid = mean((est_m2_iid - k * se_m2_iid * sqrt(39/38)) < beta & (est_m2_iid + k * se_m2_iid * sqrt(39/38)) > beta),
              m3_UN = mean((est_m3_UN - k * se_m3_UN * sqrt(39/38)) < beta & (est_m3_UN + k * se_m3_UN * sqrt(39/38)) > beta),
              m3_AR1 = mean((est_m3_AR1 - k * se_m3_AR1 * sqrt(39/38)) < beta & (est_m3_AR1 + k * se_m3_AR1 * sqrt(39/38)) > beta),
              m3_iid = mean((est_m3_iid - k * se_m3_iid * sqrt(39/38)) < beta & (est_m3_iid + k * se_m3_iid * sqrt(39/38)) > beta))
covarages2 <- arrange(coverages2, substr(coverages2$effect,17,17), substr(coverages2$effect,6,8)) # Sort with treatment first, and then time point
print(xtable(as.matrix(coverages2[,-1]), digits = 3), include.rownames = FALSE)
print(xtable(as.matrix(coverages2[,-1]), digits = 2), include.rownames = FALSE)

