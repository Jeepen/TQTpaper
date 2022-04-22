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
d$REGIMEN <- relevel(as.factor(d$REGIMEN), ref = "F")
d$PERIOD <- as.factor(d$PERIOD)
d <- mutate(d, baseline = rep(QTCF[TIME == 0], each = 6)) %>%
    arrange(SUBJECT, PERIOD, TIME) %>%
    filter(TIME != 0, !(SUBJECT %in% c(206, 229))) %>%
    mutate(PERIOD = as.factor(rep(rep(1:4, 5), 39)))
d <- group_by(d, SUBJECT) %>% mutate(meanbaseline = rep(mean(baseline), 20))
d$timeF <- factor(d$TIME, labels = c("0.5", "1.0", "1.5", "2.5", "4.0"))
basemodel <- lmm(baseline ~ 1, structure = "UN", repetition = ~ PERIOD | SUBJECT, data = filter(d, TIME == 0.5))
meanX <- rep(coef(basemodel), 4)
varX <- basemodel$Omega[[1]]

## Model
m <- lme(QTCF ~ baseline : timeF + meanbaseline : timeF + PERIOD * timeF + timeF : REGIMEN - 1,
         random = ~ 1 | SUBJECT / PERIOD,
         correlation = corSymm(form = ~ as.integer(timeF) | SUBJECT / PERIOD), weights = varIdent(form=~ 1 | TIME),
         control = list(msMaxIter = 1000, msMaxEval = 1000),
         data = d)
beta <- m$coefficients$fixed
noiseCov <- extract.lme.cov(m)[1:20,1:20]

## Simulation
nsim <- 1000
pb <- txtProgressBar(min = 1, max = nsim, style = 3)
simresults <- foreach(i = 1:nsim, .combine = "rbind", .packages = c("mgcv", "MASS", "dplyr", "LMMstar", "stats")) %dopar% {
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

  ## We need the model matrix from the following object to get a design matrix with full rank
  msim <- lme(qtc ~ baseline : timeF + meanbaseline : timeF + period * timeF + timeF : treat - 1,
              random = ~ 1 | id / period,
              correlation = corSymm(form = ~ as.integer(timeF) | id / period), weights = varIdent(form=~ 1 | time),
              control = list(msMaxIter = 1000, msMaxEval = 1000),
              data = data_sim)
  
  ## Get estimate
  estKR <- msim$coefficients$fixed[31:45]
  estNonParametric <- as.numeric(as.matrix(group_by(data_sim, time) %>% summarise(meanC = mean(qtc[treat == "C"] - qtc[treat == "F"]),
                                                                                  meanD = mean(qtc[treat == "D"] - qtc[treat == "F"]),
                                                                                  meanE = mean(qtc[treat == "E"] - qtc[treat == "F"])))[,-1])

  ## Fit mixed model
  msim2 <- lme(qtc ~ baseline : timeF + timeF * period + timeF : treat - 1,
              random = ~ 1 | id / period,
              correlation = corAR1(form = ~ as.integer(timeF) | id / period), #weights = varIdent(form=~ 1 | time),
              control = list(msMaxIter = 1000, msMaxEval = 1000),
              data = data_sim)
  ## Get estimate
  estEasy <- msim2$coefficients$fixed[26:40]


  ## data_sim <- transform(data_sim, YF = predict(msim2, rbind(data_sim, mutate(data_sim, treat = as.factor("F"))))[(nrow(data_sim)+1):(2*nrow(data_sim))],
  ##             YC = predict(msim2, rbind(data_sim, mutate(data_sim, treat = as.factor("C"))))[(nrow(data_sim)+1):(2*nrow(data_sim))],
  ##             YD = predict(msim2, rbind(data_sim, mutate(data_sim, treat = as.factor("D"))))[(nrow(data_sim)+1):(2*nrow(data_sim))],
  ##             YE = predict(msim2, rbind(data_sim, mutate(data_sim, treat = as.factor("E"))))[(nrow(data_sim)+1):(2*nrow(data_sim))])

  ## proj <- group_by(data_sim, id, time) %>% summarise(phi3C = sum(((treat == "C") - 1/4) * YC - ((treat == "F") - 1/4) * YF),
  ##                                                    phi3D = sum(((treat == "D") - 1/4) * YD - ((treat == "F") - 1/4) * YF),
  ##                                                    phi3E = sum(((treat == "E") - 1/4) * YE - ((treat == "F") - 1/4) * YF)) %>% 
  ## group_by(time) %>% summarise(meanC = mean(phi3C), meanD = mean(phi3D), meanE = mean(phi3E))
  ## semi <- estNonParametric - as.numeric(as.matrix(proj)[,-1])
  
  cat("KR est", estKR, "\n")
  cat("Simple est", estEasy, "\n")
  setTxtProgressBar(pb, i)
  data.frame(effect = names(estKR), estKR, estEasy, estNonParametric)
}

ests <- group_by(simresults, effect) %>% summarise(estKR = mean(estKR), estEasy = mean(estEasy), estNonParametric = mean(estNonParametric))
ests <- arrange(ests, substr(ests$effect,15,15), substr(ests$effect,6,8))
print(xtable(cbind(ests, beta[31:45])), include.rownames=FALSE)
print(xtable(group_by(simresults, effect) %>% summarise(estKR = sd(estKR), estEasy = sd(estEasy), estNonParametric = sd(estNonParametric))), include.rownames=FALSE)

print(xtable(group_by(simresults, effect) %>% summarise(correlation = cor(estKR, estEasy))), include.rownames = FALSE)




 

