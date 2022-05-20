## Praktisk
rm(list = ls())
set.seed(07022022)
setwd("")
library(haven)
library(LMMstar)
library(tidyverse)
library(MASS)
library(mgcv)
library(doSNOW)
library(doRNG)
library(doMC)
library(parallel)
library(xtable)
library(doParallel)
library(clubSandwich)
registerDoParallel(8)
registerDoRNG(seed = 09022022)

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
basemodel <- lmm(baseline ~ 1, structure = "CS", repetition = ~ PERIOD | SUBJECT, data = filter(d, TIME == 0.5))
meanX <- rep(coef(basemodel), 4)
varX <- basemodel$Omega[[1]]

## Model
m <- lme(QTCF ~ baseline : timeF + meanbaseline : timeF + PERIOD * timeF + timeF : REGIMEN - 1,
         random = ~ 1 | SUBJECT / PERIOD,
         correlation = corSymm(form = ~ as.integer(timeF) | SUBJECT / PERIOD), weights = varIdent(form=~ 1 | TIME),
         data = d)
beta <- m$coefficients$fixed
noiseCov <- extract.lme.cov(m)[1:20,1:20]

## Simulation
nsim <- 10000
pb <- txtProgressBar(min = 1, max = nsim, style = 3)
simresults <- foreach(i = 1:nsim, .combine = "rbind", .packages = c("mgcv", "MASS", "dplyr", "LMMstar", "stats", "clubSandwich", "estimatr")) %dorng% {
  cat("\n", "Simulation ", i, "\n")
  ## Simulate treatments
  Z <- unlist(sample(list(c(rep("C", 5), rep("D", 5), rep("F", 5), rep("E", 5)),
                          c(rep("D", 5), rep("E", 5), rep("C", 5), rep("F", 5)),
                          c(rep("E", 5), rep("F", 5), rep("D", 5), rep("C", 5)),
                          c(rep("F", 5), rep("C", 5), rep("E", 5), rep("D", 5))), 39, replace = TRUE))
  X <- rep(c(t(mvrnorm(n = 39, mu = meanX, Sigma = varX))), each = 5)
  data_sim <- data.frame(id = as.factor(rep(1:39, each = 20)),
                         treat = relevel(as.factor(Z), ref = "F"),
                         period = as.factor(rep(1:4, each = 5)),
                         time = c(0.5, 1, 1.5, 2.5, 4),
                         baseline = X)
  data_sim <- group_by(data_sim, id) %>% mutate(meanbaseline = rep(mean(baseline), 20))
  data_sim <- mutate(data_sim, timeF = factor(time, labels = c("0.5", "1.0", "1.5", "2.5", "4.0")))
  
  ## Get mean QTc measurements for each obs
  qtcMean <- numeric(780)
  for(j in 1:780){
    qtcMean[j] <- beta[substr(names(beta), 7, 7) == data_sim$period[j]] + # period effect
      beta[substr(names(beta),1,3) == "bas" & substr(names(beta), 15, 17) == data_sim$timeF[j]] * data_sim$baseline[j] + # baseline
      beta[substr(names(beta), 10, 13) == "mean" & substr(names(beta), 6, 8) == data_sim$timeF[j]] * data_sim$meanbaseline[j] # mean baseline
    if(data_sim$timeF[j] != "0.5"){
      qtcMean[j] <- qtcMean[j] + beta[substr(names(beta), 6, 8) == data_sim$timeF[j] & nchar(substr(names(beta), 1, 100)) == 8] # time point
    }
    if(data_sim$period[j] != 1 & data_sim$timeF[j] != "0.5"){
      qtcMean[j] <- qtcMean[j] + beta[which(substr(names(beta), 16, 16) == data_sim$period[j] & substr(names(beta), 6, 8) == data_sim$timeF[j])] # period : time point
    }
    if(data_sim$treat[j] != "F"){
      qtcMean[j] <- qtcMean[j] +
        beta[which(substr(names(beta),6,8) == data_sim$timeF[j] & substr(names(beta), nchar(names(beta)), nchar(names(beta))) == data_sim$treat[j])] # treatment effect                                 
    }   
  }
  
  ## Get the covariance matrix from the model object and simulate residuals
  noise <- c(t(mvrnorm(n = 39, mu = rep(0,20), Sigma = noiseCov)))
  
  ## Observed measurement equals mean + residual
  qtc <- qtcMean + noise
  data_sim$qtc <- qtc
  data_sim$timeF <- factor(data_sim$time, labels = c("0.5", "1.0", "1.5", "2.5", "4.0"))
  data_sim$treat <- as.factor(data_sim$treat)
  data_sim$treat <- relevel(data_sim$treat, ref = "F")
  
  ## True model
  ## Mean 1, unstructured
  m1_UN <- lme(qtc ~ baseline : timeF + meanbaseline : timeF + period * timeF + timeF : treat - 1,
               random = ~ 1 | id / period,
               correlation = corSymm(form = ~ as.integer(timeF) | id / period), weights = varIdent(form=~ 1 | time),
               control = lmeControl(msMaxIter = 1000, msMaxEval = 1000),
               data = data_sim)
  est_m1_UN <- m1_UN$coefficients$fixed[31:45]
  se_m1_UN <- sqrt(diag(vcovCR(m1_UN, cluster = data_sim$id, type = "CR0")))[31:45]
  
  ## Mean 1, AR(1)
  m1_AR1 <- lme(qtc ~ baseline : timeF + meanbaseline : timeF + period * timeF + timeF : treat - 1,
                random = ~ 1 | id / period,
                correlation = corAR1(form = ~ as.integer(timeF) | id / period),
                control = lmeControl(msMaxIter = 1000, msMaxEval = 1000),
                data = data_sim)
  est_m1_AR1 <- m1_AR1$coefficients$fixed[31:45]
  se_m1_AR1 <- sqrt(diag(vcovCR(m1_AR1, cluster = data_sim$id, type = "CR0")))[31:45]
  
  ## Mean 1, independence
  m1_iid <- lm_robust(qtc ~ baseline : timeF + meanbaseline : timeF + period * timeF + timeF : treat - 1,
                      data = data_sim, clusters = id)
  est_m1_iid <- m1_iid$coefficients[31:45]
  se_m1_iid <- sqrt(diag(vcov(m1_iid)))[31:45]

  ## Mean 2, unstructured
  m2_UN <- lme(qtc ~ baseline : timeF + period * timeF + timeF : treat - 1,
               random = ~ 1 | id / period,
               correlation = corSymm(form = ~ as.integer(timeF) | id / period), weights = varIdent(form=~ 1 | time),
               control = lmeControl(msMaxIter = 1000, msMaxEval = 1000),
               data = data_sim)
  est_m2_UN <- m2_UN$coefficients$fixed[26:40]
  se_m2_UN <- sqrt(diag(vcovCR(m2_UN, cluster = data_sim$id, type = "CR0")))[26:40]
  
  ## Mean 2, AR(1)
  m2_AR1 <- lme(qtc ~ baseline : timeF + period * timeF + timeF : treat - 1,
                random = ~ 1 | id / period,
                correlation = corAR1(form = ~ as.integer(timeF) | id / period),
                control = lmeControl(msMaxIter = 1000, msMaxEval = 1000),
                data = data_sim)
  est_m2_AR1 <- m2_AR1$coefficients$fixed[26:40]
  se_m2_AR1 <- sqrt(diag(vcovCR(m2_AR1, cluster = data_sim$id, type = "CR0")))[26:40]
  
  ## Mean 2, independence
  m2_iid <- lm_robust(qtc ~ baseline : timeF + period * timeF + timeF : treat - 1,
                      data = data_sim, clusters = id)
  est_m2_iid <- m2_iid$coefficients[26:40]
  se_m2_iid <- sqrt(diag(vcov(m2_iid)))[26:40]

  ## Mean 3, unstructured
  m3_UN <- lme(qtc ~ baseline + timeF + timeF : treat - 1,
               random = ~ 1 | id / period,
               correlation = corSymm(form = ~ as.integer(timeF) | id / period), weights = varIdent(form=~ 1 | time),
               control = lmeControl(msMaxIter = 1000, msMaxEval = 1000),
               data = data_sim)
  est_m3_UN <- m3_UN$coefficients$fixed[7:21]
  se_m3_UN <- sqrt(diag(vcovCR(m3_UN, cluster = data_sim$id, type = "CR0")))[7:21]
  
  ## Mean 3, AR(1)
  m3_AR1 <- lme(qtc ~ baseline + timeF + timeF : treat - 1,
                random = ~ 1 | id / period,
                correlation = corAR1(form = ~ as.integer(timeF) | id / period),
                control = lmeControl(msMaxIter = 1000, msMaxEval = 1000),
                data = data_sim)
  est_m3_AR1 <- m3_AR1$coefficients$fixed[7:21]
  se_m3_AR1 <- sqrt(diag(vcovCR(m3_AR1, cluster = data_sim$id, type = "CR0")))[7:21]
  
  ## Mean 3, independence
  m3_iid <- lm_robust(qtc ~ baseline + timeF + timeF : treat - 1,
                      data = data_sim, clusters = id)
  est_m3_iid <- m3_iid$coefficients[7:21]
  se_m3_iid <- sqrt(diag(vcov(m3_iid)))[7:21]

  ## Mean 4, unstructured
  m4_UN <- lme(qtc ~ baseline : treat + timeF + timeF : treat - 1,
               random = ~ 1 | id / period,
               correlation = corSymm(form = ~ as.integer(timeF) | id / period), weights = varIdent(form=~ 1 | time),
               control = lmeControl(msMaxIter = 1000, msMaxEval = 1000),
               data = data_sim)
  data_sim <- transform(data_sim, YF = predict(m4_UN, rbind(data_sim, mutate(data_sim, treat = as.factor("F"))))[(nrow(data_sim)+1):(2*nrow(data_sim))],
                        YC = predict(m4_UN, rbind(data_sim, mutate(data_sim, treat = as.factor("C"))))[(nrow(data_sim)+1):(2*nrow(data_sim))],
                        YD = predict(m4_UN, rbind(data_sim, mutate(data_sim, treat = as.factor("D"))))[(nrow(data_sim)+1):(2*nrow(data_sim))],
                        YE = predict(m4_UN, rbind(data_sim, mutate(data_sim, treat = as.factor("E"))))[(nrow(data_sim)+1):(2*nrow(data_sim))])
  ests <- group_by(data_sim, time) %>% summarise(C = mean(YC - YF),
                                                 D = mean(YD - YF),
                                                 E = mean(YE - YF))
  est_m4_UN <- as.numeric(as.matrix(ests)[,-1])
  mu3SE <- group_by(data_sim, id, time) %>% 
      summarise(phi3C = qtc[treat == "C"] - qtc[treat == "F"] - sum(((treat == "C") - 1/4) * YC - ((treat == "F") - 1/4) * YF),
                phi3D = qtc[treat == "D"] - qtc[treat == "F"] - sum(((treat == "D") - 1/4) * YD - ((treat == "F") - 1/4) * YF),
                phi3E = qtc[treat == "E"] - qtc[treat == "F"] - sum(((treat == "E") - 1/4) * YE - ((treat == "F") - 1/4) * YF)) #%>% 
  for(t in c("0.5", "1", "1.5", "2.5", "4")){
      for(j in 2:4){
          mu3SE[mu3SE$time == t, j+1] <- mu3SE[mu3SE$time == t, j+1] - as.numeric(ests[ests$time == t, j])
      }
  }
  mu3SE <- group_by(mu3SE, time) %>% summarise(seC = sqrt(mean(phi3C^2) / 39),
                                               seD = sqrt(mean(phi3D^2) / 39),
                                               seE = sqrt(mean(phi3E^2) / 39))
  se_m4_UN <- as.numeric(as.matrix(mu3SE)[,-1])  

  ## Mean 4, AR(1)
  m4_AR1 <- lme(qtc ~ baseline : treat + timeF + timeF : treat - 1,
               random = ~ 1 | id / period,
               correlation = corAR1(form = ~ as.integer(timeF) | id / period),
               control = lmeControl(msMaxIter = 1000, msMaxEval = 1000),
               data = data_sim)
  data_sim <- transform(data_sim, YF = predict(m4_AR1, rbind(data_sim, mutate(data_sim, treat = as.factor("F"))))[(nrow(data_sim)+1):(2*nrow(data_sim))],
                        YC = predict(m4_AR1, rbind(data_sim, mutate(data_sim, treat = as.factor("C"))))[(nrow(data_sim)+1):(2*nrow(data_sim))],
                        YD = predict(m4_AR1, rbind(data_sim, mutate(data_sim, treat = as.factor("D"))))[(nrow(data_sim)+1):(2*nrow(data_sim))],
                        YE = predict(m4_AR1, rbind(data_sim, mutate(data_sim, treat = as.factor("E"))))[(nrow(data_sim)+1):(2*nrow(data_sim))])
  ests <- group_by(data_sim, time) %>% summarise(C = mean(YC - YF),
                                                 D = mean(YD - YF),
                                                 E = mean(YE - YF))
  est_m4_AR1 <- as.numeric(as.matrix(ests)[,-1])
  mu3SE <- group_by(data_sim, id, time) %>% 
      summarise(phi3C = qtc[treat == "C"] - qtc[treat == "F"] - sum(((treat == "C") - 1/4) * YC - ((treat == "F") - 1/4) * YF),
                phi3D = qtc[treat == "D"] - qtc[treat == "F"] - sum(((treat == "D") - 1/4) * YD - ((treat == "F") - 1/4) * YF),
                phi3E = qtc[treat == "E"] - qtc[treat == "F"] - sum(((treat == "E") - 1/4) * YE - ((treat == "F") - 1/4) * YF)) #%>% 
  for(t in c("0.5", "1", "1.5", "2.5", "4")){
      for(j in 2:4){
          mu3SE[mu3SE$time == t, j+1] <- mu3SE[mu3SE$time == t, j+1] - as.numeric(ests[ests$time == t, j])
      }
  }
  mu3SE <- group_by(mu3SE, time) %>% summarise(seC = sqrt(mean(phi3C^2) / 39),
                                               seD = sqrt(mean(phi3D^2) / 39),
                                               seE = sqrt(mean(phi3E^2) / 39))
  se_m4_AR1 <- as.numeric(as.matrix(mu3SE)[,-1])

  ## Mean 4, independence
  m4_iid <- lm(qtc ~ baseline : treat + timeF + timeF : treat - 1, data = data_sim)
  data_sim <- transform(data_sim, YF = predict(m4_iid, rbind(data_sim, mutate(data_sim, treat = as.factor("F"))))[(nrow(data_sim)+1):(2*nrow(data_sim))],
                        YC = predict(m4_iid, mutate(data_sim, treat = as.factor("C"))),
                        YD = predict(m4_iid, mutate(data_sim, treat = as.factor("D"))),
                        YE = predict(m4_iid, mutate(data_sim, treat = as.factor("E"))))
  ests <- group_by(data_sim, time) %>% summarise(C = mean(YC - YF),
                                                 D = mean(YD - YF),
                                                 E = mean(YE - YF))
  est_m4_iid <- as.numeric(as.matrix(ests)[,-1])
  mu3SE <- group_by(data_sim, id, time) %>% 
      summarise(phi3C = qtc[treat == "C"] - qtc[treat == "F"] - sum(((treat == "C") - 1/4) * YC - ((treat == "F") - 1/4) * YF),
                phi3D = qtc[treat == "D"] - qtc[treat == "F"] - sum(((treat == "D") - 1/4) * YD - ((treat == "F") - 1/4) * YF),
                phi3E = qtc[treat == "E"] - qtc[treat == "F"] - sum(((treat == "E") - 1/4) * YE - ((treat == "F") - 1/4) * YF)) #%>% 
  for(t in c("0.5", "1", "1.5", "2.5", "4")){
      for(j in 2:4){
          mu3SE[mu3SE$time == t, j+1] <- mu3SE[mu3SE$time == t, j+1] - as.numeric(ests[ests$time == t, j])
      }
  }
  mu3SE <- group_by(mu3SE, time) %>% summarise(seC = sqrt(mean(phi3C^2) / 39),
                                               seD = sqrt(mean(phi3D^2) / 39),
                                               seE = sqrt(mean(phi3E^2) / 39))
  se_m4_iid <- as.numeric(as.matrix(mu3SE)[,-1])
  
  ## Non-parametric estimator  
  est_nonparametric <- as.numeric(as.matrix(group_by(data_sim, time) %>% summarise(meanC = mean(qtc[treat == "C"] - qtc[treat == "F"]),
                                                                                   meanD = mean(qtc[treat == "D"] - qtc[treat == "F"]),
                                                                                   meanE = mean(qtc[treat == "E"] - qtc[treat == "F"])))[,-1])
  se_nonparametric <- round(as.numeric(as.matrix(group_by(data_sim, id, time) %>% summarise(CF = qtc[treat == "C"] - qtc[treat == "F"],
                                                                                            DF = qtc[treat == "D"] - qtc[treat == "F"],
                                                                                            EF = qtc[treat == "E"] - qtc[treat == "F"]) %>% 
                                                 group_by(time) %>% summarise(seC = sd(CF) / sqrt(with(d, length(unique(id)))), 
                                                                              seD = sd(DF) / sqrt(with(d, length(unique(id)))),
                                                                              seE = sd(EF) / sqrt(with(d, length(unique(id))))))[,-1]), 2)

  ## Save results
  setTxtProgressBar(pb, i)
  data.frame(effect = names(beta)[31:45], est_m1_UN, est_m1_AR1, est_m1_iid,
             est_m2_UN, est_m2_AR1, est_m2_iid,
             est_m3_UN, est_m3_AR1, est_m3_iid,
             est_m4_UN, est_m4_AR1, est_m4_iid,
             est_nonparametric,
             se_m1_UN, se_m1_AR1, se_m1_iid,
             se_m2_UN, se_m2_AR1, se_m2_iid,
             se_m3_UN, se_m3_AR1, se_m3_iid,
             se_m4_UN, se_m4_AR1, se_m4_iid,
             se_nonparametric
             )
}
simresults <- transform(simresults, beta = rep(beta[31:45], 1000))
saveRDS(simresults, "")

