## Praktisk
rm(list = ls())
set.seed(07022022)
setwd("~/Dropbox/phd/LSMeans/Code/")
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
registerDoParallel()
registerDoRNG(seed = 09022022)

## Data
d <- read_sas("exam9_1.sas7bdat")
d$SUBJECT <- as.factor(d$SUBJECT)
d$REGIMEN[d$TIME == 0] <- "bas"
d$REGIMEN <- relevel(as.factor(d$REGIMEN), ref = "F")
d$PERIOD <- as.factor(d$PERIOD)
d$timeF <- as.factor(d$TIME)
d <- mutate(d, SubPer = SUBJECT : PERIOD, regimeTime = REGIMEN : timeF, periodTime = PERIOD : timeF) %>% arrange(SUBJECT, PERIOD, TIME)
d$regimeTime <- droplevels(d$regimeTime)
d$periodTime <- droplevels(d$periodTime)

## Model 
m <- lmm(QTCF ~ regimeTime + periodTime, structure = "ID", repetition = ~TIME | SubPer, data = d, df = FALSE, type.information="expected")
mmat <- model.matrix(m)
model <- lme(QTCF ~ mmat[,-1], random=~ 1 | SUBJECT / PERIOD, data = d, na.action = na.omit,
             correlation = corSymm(form = ~ as.integer(timeF) | SUBJECT / PERIOD))
beta <- model$coefficients$fixed
names(beta)[-1] <- substr(names(beta)[-1], 11, 50)

## Simulation
nsim <- 1000
pb <- txtProgressBar(min = 1, max = nsim, style = 3)
simresults <- foreach(i = 1:nsim, .combine = "rbind", .packages = c("mgcv", "MASS", "dplyr", "LMMstar", "stats")) %dopar% {
  cat("\n", "Iteration: ", i, "\n")
  ## Simulate treatments
  Z <- unlist(sample(list(c("bas", rep("C", 5), "bas", rep("D", 5), "bas", rep("F", 5), "bas", rep("E", 5)),
                          c("bas", rep("D", 5), "bas", rep("E", 5), "bas", rep("C", 5), "bas", rep("F", 5)),
                          c("bas", rep("E", 5), "bas", rep("F", 5), "bas", rep("D", 5), "bas", rep("C", 5)),
                          c("bas", rep("F", 5), "bas", rep("C", 5), "bas", rep("E", 5), "bas", rep("D", 5))), 41, replace = TRUE))
  data_sim <- data.frame(id = as.factor(rep(1:41, each = 24)), treat = Z, period = as.factor(rep(1:4, each = 6)), time = c(0, 0.5, 1, 1.5, 2.5, 4))

  ## Get mean QTc measurements for each obs
  qtcMean <- numeric(984)
  for(j in 1:984){
    if(data_sim$treat[j] == "bas"){
      qtcMean[j] <- beta[1]
    }
    else if(data_sim$treat[j] == "F"){
      qtcMean[j] <- beta[1] + beta[which(substr(names(beta), 13, 15) == data_sim$time[j] & substr(names(beta), 11, 11) == data_sim$period[j])]
    }   
    else{
      qtcMean[j] <- beta[1] + beta[which(substr(names(beta), 13, 15) == data_sim$time[j] & substr(names(beta), 11, 11) == data_sim$period[j])] +
        beta[which(substr(names(beta), 13, 15) == data_sim$time[j] & substr(names(beta), 11, 11) == data_sim$treat[j])]
    }
  }

  ## Get the covariance matrix from the model object and simulate residuals
  noiseCov <- extract.lme.cov(model)[1:24,1:24]
  noise <- c(t(mvrnorm(n = 41, mu = rep(0,24), Sigma = noiseCov)))

  ## Observed measurement equals mean + residual
  qtc <- qtcMean + noise
  data_sim$qtc <- qtc
  data_sim$timeF <- as.factor(data_sim$time)
  data_sim$treat <- as.factor(data_sim$treat)
  data_sim$treat <- relevel(data_sim$treat, ref = "F")
  data_sim <- mutate(data_sim, SubPer = id : period, regimeTime = treat : timeF, periodTime = period : timeF) %>% arrange(id, period, time)
  data_sim$regimeTime <- droplevels(data_sim$regimeTime)
  data_sim$periodTime <- droplevels(data_sim$periodTime)

  ## We need the model matrix from the following object to get a design matrix with full rank
  msim <- lmm(qtc ~ regimeTime + periodTime, structure = "ID", repetition = ~time | SubPer, data = data_sim, df = FALSE, type.information="expected",
              method.fit = "ML")
  mmatsim <<- model.matrix(msim)
  
  ## Fit Kenward Roger model
  modelsim <- lme(qtc ~ mmatsim[,-1], random=~ 1 | id / period, data = data_sim,
                  correlation = corSymm(form = ~ as.integer(timeF) | id / period))
  ## Get estimate
  estKR <- modelsim$coefficients$fixed[2:16]

  ## Change data so baseline can be used as covariate
  data_sim <- mutate(data_sim, baseline = rep(qtc[time == 0], each = 6)) %>% filter(time != 0)

  ## Fit mixed model
  modelsim2 <- lme(qtc ~ baseline + period + timeF + period : timeF + treat : timeF,
                   random=~ 1 | id / period, data = data_sim, correlation = corAR1(form = ~ as.integer(timeF) | id / period))
  ## Get estimate
  estEasy <- modelsim2$coefficients$fixed[22:36]
  cat("KR est", estKR, "\n")
  cat("Simple est", estEasy, "\n")
  setTxtProgressBar(pb, i)
  data.frame(effect = substr(names(estKR), 24, 50), estKR, estEasy)
}

ests <- group_by(simresults, effect) %>% summarise(estKR = mean(estKR), estEasy = mean(estEasy))
print(xtable(cbind(ests, model$coefficients$fixed[2:16])), include.rownames=FALSE)
print(xtable(group_by(simresults, effect) %>% summarise(estKR = sd(estKR), estEasy = sd(estEasy))), include.rownames=FALSE)

round(apply(simresults, 2, mean), 2) # estKR = 8.33, estEasy = 8.34 - Average estimate
round(model$coefficients$fixed[16], 2)  # 8.35 - the truth in the simulations
round(apply(simresults, 2, sd), 2) # estKR = 1.5, estEasy = 1.56 - standard error of estimates in the simulations





 

