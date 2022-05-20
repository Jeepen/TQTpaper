## Praktisk
rm(list = ls())
set.seed(07022022)
setwd("")
library(haven)
library(LMMstar)
library(tidyverse)
library(MASS)
library(mgcv)
library(xtable)
library(estimatr)
library(clubSandwich)

## Data
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
k <- qt(.975, df = 38)

## Models
## Mean 1, unstructured
m1_UN <- lme(QTCF ~ baseline : timeF + meanbaseline : timeF + PERIOD * timeF + timeF : REGIMEN - 1,
             random = ~ 1 | SUBJECT / PERIOD,
             correlation = corSymm(form = ~ as.integer(timeF) | SUBJECT / PERIOD), weights = varIdent(form=~ 1 | TIME),
             data = d)
beta <- m1_UN$coefficients$fixed[31:45]
round(beta, 2)
tibble(sqrt(diag(vcovCR(m1_UN, cluster = d$SUBJECT, type = "CR1")))[31:45])
cbind(beta - k * sqrt(diag(vcovCR(m1_UN, cluster = d$SUBJECT, type = "CR1")))[31:45],
      beta + k * sqrt(diag(vcovCR(m1_UN, cluster = d$SUBJECT, type = "CR1")))[31:45])

## Mean 1, AR(1)
m1_AR1 <- lme(QTCF ~ baseline : timeF + meanbaseline : timeF + PERIOD * timeF + timeF : REGIMEN - 1,
              random = ~ 1 | SUBJECT / PERIOD,
              correlation = corAR1(form = ~ as.integer(timeF) | SUBJECT / PERIOD), 
              data = d)
beta <- m1_AR1$coefficients$fixed[31:45]
round(beta, 2)
tibble(sqrt(diag(vcovCR(m1_AR1, cluster = d$SUBJECT, type = "CR1")))[31:45])
cbind(beta - k * sqrt(diag(vcovCR(m1_AR1, cluster = d$SUBJECT, type = "CR1")))[31:45],
      beta + k * sqrt(diag(vcovCR(m1_AR1, cluster = d$SUBJECT, type = "CR1")))[31:45])

## Mean 1, independence
m1_iid <- lm_robust(QTCF ~ baseline : timeF + meanbaseline : timeF + PERIOD * timeF + timeF : REGIMEN - 1, data = d, clusters = SUBJECT)
beta <- m1_iid$coefficients[31:45]
round(beta, 2)
round(sqrt(diag(vcov(m1_iid)) * 39/38), 2)[31:45]
cbind(beta - k * sqrt(diag(vcov(m1_iid)) * 39/38)[31:45],
      beta + k * sqrt(diag(vcov(m1_iid)) * 39/38)[31:45])

## Mean 2, unstructured
m2_UN <- lme(QTCF ~ baseline : timeF + PERIOD * timeF + timeF : REGIMEN - 1,
             random = ~ 1 | SUBJECT / PERIOD,
             correlation = corSymm(form = ~ as.integer(timeF) | SUBJECT / PERIOD), weights = varIdent(form=~ 1 | TIME),
             data = d)
beta <- m2_UN$coefficients$fixed[26:40]
round(beta, 2)
tibble(sqrt(diag(vcovCR(m2_UN, cluster = d$SUBJECT, type = "CR1")))[26:40])
cbind(beta - k * sqrt(diag(vcovCR(m2_UN, cluster = d$SUBJECT, type = "CR1")))[26:40],
      beta + k * sqrt(diag(vcovCR(m2_UN, cluster = d$SUBJECT, type = "CR1")))[26:40])

## Mean 2, AR(1)
m2_AR1 <- lme(QTCF ~ baseline : timeF + PERIOD * timeF + timeF : REGIMEN - 1,
              random = ~ 1 | SUBJECT / PERIOD,
              correlation = corAR1(form = ~ as.integer(timeF) | SUBJECT / PERIOD), 
              data = d)
beta <- m2_AR1$coefficients$fixed[26:40]
round(beta, 2)
tibble(sqrt(diag(vcovCR(m2_AR1, cluster = d$SUBJECT, type = "CR1")))[26:40])
cbind(beta - k * sqrt(diag(vcovCR(m2_AR1, cluster = d$SUBJECT, type = "CR1")))[26:40],
      beta + k * sqrt(diag(vcovCR(m2_AR1, cluster = d$SUBJECT, type = "CR1")))[26:40])

## Mean 2, independence
m2_iid <- lm_robust(QTCF ~ baseline : timeF + PERIOD * timeF + timeF : REGIMEN - 1, data = d, clusters = SUBJECT)
beta <- m2_iid$coefficients[26:40]
round(beta, 2)
round(sqrt(diag(vcov(m2_iid)) * 39/38), 2)[26:40]
cbind(beta - k * sqrt(diag(vcov(m2_iid)) * 39/38)[26:40],
      beta + k * sqrt(diag(vcov(m2_iid)) * 39/38)[26:40])

## Mean 3, unstructured
m3_UN <- lme(QTCF ~ baseline + timeF + timeF : REGIMEN - 1,
             random = ~ 1 | SUBJECT / PERIOD,
             correlation = corSymm(form = ~ as.integer(timeF) | SUBJECT / PERIOD), weights = varIdent(form=~ 1 | TIME),
             data = d)
beta <- m3_UN$coefficients$fixed[7:21]
round(beta, 2)
tibble(sqrt(diag(vcovCR(m3_UN, cluster = d$SUBJECT, type = "CR1")))[7:21])
cbind(beta - k * sqrt(diag(vcovCR(m3_UN, cluster = d$SUBJECT, type = "CR1")))[7:21],
      beta + k * sqrt(diag(vcovCR(m3_UN, cluster = d$SUBJECT, type = "CR1")))[7:21])

## Mean 3, AR(1)
m3_AR1 <- lme(QTCF ~ baseline + timeF + timeF : REGIMEN - 1,
              random = ~ 1 | SUBJECT / PERIOD,
              correlation = corAR1(form = ~ as.integer(timeF) | SUBJECT / PERIOD), 
              data = d)
beta <- m3_AR1$coefficients$fixed[7:21]
round(beta, 2)
round(sqrt(diag(vcov(m3_AR1)))[7:21], 2)
tibble(sqrt(diag(vcovCR(m3_AR1, cluster = d$SUBJECT, type = "CR1")))[7:21])
cbind(beta - k * sqrt(diag(vcovCR(m3_AR1, cluster = d$SUBJECT, type = "CR1")))[7:21],
      beta + k * sqrt(diag(vcovCR(m3_AR1, cluster = d$SUBJECT, type = "CR1")))[7:21])

## Mean 3, independence
m3_iid <- lm_robust(QTCF ~ baseline + timeF + timeF : REGIMEN - 1, data = d, clusters = SUBJECT)
beta <- m3_iid$coefficients[7:21]
round(beta, 2)
round(sqrt(diag(vcov(m3_iid)) * 39/38), 2)[7:21]
cbind(beta - k * sqrt(diag(vcov(m3_iid)) * 39/38)[7:21],
      beta + k * sqrt(diag(vcov(m3_iid)) * 39/38)[7:21])

## Mean 4, unstructured
m4_UN <- lme(QTCF ~ baseline : REGIMEN + timeF + timeF : REGIMEN - 1,
             random = ~ 1 | SUBJECT / PERIOD,
             correlation = corSymm(form = ~ as.integer(timeF) | SUBJECT / PERIOD), weights = varIdent(form=~ 1 | TIME),
             data = d)
d <- transform(d, YF = predict(m4_UN, rbind(d, mutate(d, REGIMEN = as.factor("F"))))[(nrow(d)+1):(2*nrow(d))],
               YC = predict(m4_UN, rbind(d, mutate(d, REGIMEN = as.factor("C"))))[(nrow(d)+1):(2*nrow(d))],
               YD = predict(m4_UN, rbind(d, mutate(d, REGIMEN = as.factor("D"))))[(nrow(d)+1):(2*nrow(d))],
               YE = predict(m4_UN, rbind(d, mutate(d, REGIMEN = as.factor("E"))))[(nrow(d)+1):(2*nrow(d))])
est_m4_UN <- as.numeric(as.matrix(group_by(d, TIME) %>% summarise(estC = mean(YC - YF), estD = mean(YD - YF), estE = mean(YE - YF)))[,-1])
round(est_m4_UN, 2)
firstterms <- group_by(d, SUBJECT, TIME) %>% summarise(C = mean(YC - YF),
                                                       D = mean(YD - YF),
                                                       E = mean(YE - YF))
ests <- group_by(d, TIME) %>% summarise(C = mean(YC - YF),
                                        D = mean(YD - YF),
                                        E = mean(YE - YF))
mu3SE <- group_by(d, SUBJECT, TIME) %>% 
    summarise(phi3C = QTCF[REGIMEN == "C"] - QTCF[REGIMEN == "F"] - sum(((REGIMEN == "C") - 1/4) * YC - ((REGIMEN == "F") - 1/4) * YF),
              phi3D = QTCF[REGIMEN == "D"] - QTCF[REGIMEN == "F"] - sum(((REGIMEN == "D") - 1/4) * YD - ((REGIMEN == "F") - 1/4) * YF),
              phi3E = QTCF[REGIMEN == "E"] - QTCF[REGIMEN == "F"] - sum(((REGIMEN == "E") - 1/4) * YE - ((REGIMEN == "F") - 1/4) * YF)) #%>% 
for(t in c("0.5", "1", "1.5", "2.5", "4")){
  for(j in 2:4){
    mu3SE[mu3SE$TIME == t, j+1] <- mu3SE[mu3SE$TIME == t, j+1] - as.numeric(ests[ests$TIME == t, j])
  }
}
mu3SE <- group_by(mu3SE, TIME) %>% summarise(seC = sqrt(mean(phi3C^2) / 38),
                                             seD = sqrt(mean(phi3D^2) / 38),
                                             seE = sqrt(mean(phi3E^2) / 38))
cbind(ests[5,4] - k * mu3SE[5,4], ests[5,4] + k * mu3SE[5,4])

## Mean 4, AR(1)
m4_AR1 <- lme(QTCF ~ baseline : REGIMEN + timeF + timeF : REGIMEN - 1,
              random = ~ 1 | SUBJECT / PERIOD,
              correlation = corAR1(form = ~ as.integer(timeF) | SUBJECT / PERIOD), 
              data = d)
d <- transform(d, YF = predict(m4_AR1, rbind(d, mutate(d, REGIMEN = as.factor("F"))))[(nrow(d)+1):(2*nrow(d))],
               YC = predict(m4_AR1, rbind(d, mutate(d, REGIMEN = as.factor("C"))))[(nrow(d)+1):(2*nrow(d))],
               YD = predict(m4_AR1, rbind(d, mutate(d, REGIMEN = as.factor("D"))))[(nrow(d)+1):(2*nrow(d))],
               YE = predict(m4_AR1, rbind(d, mutate(d, REGIMEN = as.factor("E"))))[(nrow(d)+1):(2*nrow(d))])
est_m4_AR1 <- as.numeric(as.matrix(group_by(d, TIME) %>% summarise(estC = mean(YC - YF), estD = mean(YD - YF), estE = mean(YE - YF)))[,-1])
round(est_m4_AR1, 2)
firstterms <- group_by(d, SUBJECT, TIME) %>% summarise(C = mean(YC - YF),
                                                       D = mean(YD - YF),
                                                       E = mean(YE - YF))
ests <- group_by(d, TIME) %>% summarise(C = mean(YC - YF),
                                        D = mean(YD - YF),
                                        E = mean(YE - YF))
mu3SE <- group_by(d, SUBJECT, TIME) %>% 
    summarise(phi3C = QTCF[REGIMEN == "C"] - QTCF[REGIMEN == "F"] - sum(((REGIMEN == "C") - 1/4) * YC - ((REGIMEN == "F") - 1/4) * YF),
              phi3D = QTCF[REGIMEN == "D"] - QTCF[REGIMEN == "F"] - sum(((REGIMEN == "D") - 1/4) * YD - ((REGIMEN == "F") - 1/4) * YF),
              phi3E = QTCF[REGIMEN == "E"] - QTCF[REGIMEN == "F"] - sum(((REGIMEN == "E") - 1/4) * YE - ((REGIMEN == "F") - 1/4) * YF)) #%>% 
for(t in c("0.5", "1", "1.5", "2.5", "4")){
  for(j in 2:4){
    mu3SE[mu3SE$TIME == t, j+1] <- mu3SE[mu3SE$TIME == t, j+1] - as.numeric(ests[ests$TIME == t, j])
  }
}
mu3SE <- group_by(mu3SE, TIME) %>% summarise(seC = sqrt(mean(phi3C^2) / 38),
                                             seD = sqrt(mean(phi3D^2) / 38),
                                             seE = sqrt(mean(phi3E^2) / 38))
cbind(ests[5,4] - 1.96 * mu3SE[5,4], ests[5,4] + 1.96 * mu3SE[5,4])

## Mean 4, independence
m4_iid <- lm(QTCF ~ baseline : REGIMEN + timeF + timeF : REGIMEN - 1, data = d)
d <- transform(d, YF = predict(m4_iid, mutate(d, REGIMEN = as.factor("F"))),
               YC = predict(m4_iid, mutate(d, REGIMEN = as.factor("C"))),
               YD = predict(m4_iid, mutate(d, REGIMEN = as.factor("D"))),
               YE = predict(m4_iid, mutate(d, REGIMEN = as.factor("E"))))
est_m4_iid <- as.numeric(as.matrix(group_by(d, TIME) %>% summarise(estC = mean(YC - YF), estD = mean(YD - YF), estE = mean(YE - YF)))[,-1])
round(est_m4_iid, 2)
firstterms <- group_by(d, SUBJECT, TIME) %>% summarise(C = mean(YC - YF),
                                                       D = mean(YD - YF),
                                                       E = mean(YE - YF))
ests <- group_by(d, TIME) %>% summarise(C = mean(YC - YF),
                                        D = mean(YD - YF),
                                        E = mean(YE - YF))
mu3SE <- group_by(d, SUBJECT, TIME) %>% 
    summarise(phi3C = QTCF[REGIMEN == "C"] - QTCF[REGIMEN == "F"] - sum(((REGIMEN == "C") - 1/4) * YC - ((REGIMEN == "F") - 1/4) * YF),
              phi3D = QTCF[REGIMEN == "D"] - QTCF[REGIMEN == "F"] - sum(((REGIMEN == "D") - 1/4) * YD - ((REGIMEN == "F") - 1/4) * YF),
              phi3E = QTCF[REGIMEN == "E"] - QTCF[REGIMEN == "F"] - sum(((REGIMEN == "E") - 1/4) * YE - ((REGIMEN == "F") - 1/4) * YF)) #%>% 
for(t in c("0.5", "1", "1.5", "2.5", "4")){
  for(j in 2:4){
    mu3SE[mu3SE$TIME == t, j+1] <- mu3SE[mu3SE$TIME == t, j+1] - as.numeric(ests[ests$TIME == t, j])
  }
}
mu3SE <- group_by(mu3SE, TIME) %>% summarise(seC = sqrt(mean(phi3C^2) / 38),
                                             seD = sqrt(mean(phi3D^2) / 38),
                                             seE = sqrt(mean(phi3E^2) / 38))
cbind(ests[5,4] - k * mu3SE[5,4], ests[5,4] + k * mu3SE[5,4])

## Non-parametric
estNonParametric <- as.numeric(as.matrix(group_by(d, TIME) %>% summarise(meanC = mean(QTCF[REGIMEN == "C"] - QTCF[REGIMEN == "F"]),
                                                                         meanD = mean(QTCF[REGIMEN == "D"] - QTCF[REGIMEN == "F"]),
                                                                         meanE = mean(QTCF[REGIMEN == "E"] - QTCF[REGIMEN == "F"])))[,-1])
round(estNonParametric, 2)
SE_nonparametric <- round(as.numeric(as.matrix(group_by(d, SUBJECT, TIME) %>% summarise(CF = QTCF[REGIMEN == "C"] - QTCF[REGIMEN == "F"],
                                                                                        DF = QTCF[REGIMEN == "D"] - QTCF[REGIMEN == "F"],
                                                                                        EF = QTCF[REGIMEN == "E"] - QTCF[REGIMEN == "F"]) %>% 
                                               group_by(TIME) %>% summarise(seC = sd(CF) / sqrt(with(d, length(unique(SUBJECT)))), 
                                                                            seD = sd(DF) / sqrt(with(d, length(unique(SUBJECT)))),
                                                                            seE = sd(EF) / sqrt(with(d, length(unique(SUBJECT))))))[,-1]), 2)
estNonParametric[15] + c(-k, k) * SE_nonparametric[15]
