#Load in relevant packages
library(magrittr)
library(stats4)
library(bbmle)
library(rAmCharts)
library(LearnClust)
library(dplyr)
library(haven)
library(boot)
library(ggplot2)

######################Data preperation#########################################

#Load in data - called uam

#Restrict to survey years 2011-2020
uam = uam %>%
  subset(year %in% c("2011","2012","2013","2014","2015",
                     "2016","2017","2018","2019","2020"))

#Assume injecting duration of one
uam$injdur[which(uam$injdur == 0)] = 1

#Only keep the complete cases for the variables of interest in this analysis
keep = c("injdur", "hcvdbs", "year")
uam = uam[keep]
uam = uam[complete.cases(uam), ]

#######################Piecewise constant#######################################
#Aggregate data for the PC model
uam1 = aggregate(hcvdbs ~ year + injdur, data = uam, sum, na.rm = T)
colnames(uam1) = c("year", "injdur", "hcv")
uam2 = uam
uam2$hcvdbs = 1
uam2 = aggregate(hcvdbs ~ year + injdur, data = uam2, sum, na.rm = T)
colnames(uam2) = c("year", "injdur", "n")
data = merge(uam1, uam2, by=c("year", "injdur"))

#Load in bespoke functions
setwd("C:/Users/Conor.Egan/Documents/Ross FOI code")
source("fn_mkcut.R")
source("fn_univarATnoint.R")
source("fn_getparamATnoint.R")
source("fn_pred.R")

#Calculate the injecting duration and time contributions
data.hcv = mkcut(data, data$year, data$injdur,
                 c(1900,1980,1990,1995,200,2005,2010,2015,2020),
                 c(0,1,3,5,10,15,25,60))

#Include sensitivity
data.hcv$senshcv = 1

#Initial values for the multiplicative no interaction model
iv.t = rep(-4,8)
iv.a = c(1,0,0,0,0,0)
inits2 = c(iv.t,iv.a)

#Optimise unknown parameters for model
modc2 = optim(inits2, univarATnoint, method = "BFGS", y = data.hcv$hcv,
              n = data.hcv$n, s = data.hcv$senshcv, data = data.hcv,
              control = list(trace = T, maxit = 500), hessian = T)

#Obtain the suvival probability predictions
modc2.params = getparamATnoint(modc2)
modc2.pred = pred(modc2.params$h, data.hcv, data.hcv$hcv, data.hcv$n,
                  data.hcv$senshcv)
hcv.pred = cbind(data.hcv$year, data.hcv$injdur, data.hcv$n, data.hcv$senshcv)
colnames(hcv.pred) = c("year", "injdur", "n", "s")
hcv.pred = cbind(hcv.pred, modc2.pred)
hcv.pred = as.data.frame(hcv.pred)

#Use bootstrapping to get confidence intervals
set.seed(100)

#Make a function that returns the desired values (prob ci or foi ci)
h_CI = function(data, indices) {
  data = data[indices, ]
  #Same process as above
  data1 = aggregate(hcvdbs ~ year + injdur, data = data, sum, na.rm = T)
  colnames(data1) = c("year", "injdur", "hcv")
  data2 = data
  data2$hcvdbs = 1
  data2 = aggregate(hcvdbs ~ year + injdur, data = data2, sum, na.rm = T)
  colnames(data2) = c("year", "injdur", "n")
  data = merge(data1, data2, by=c("year", "injdur"))
  data.hcv = mkcut(data, data$year, data$injdur,
                   c(1900,1980,1990,1995,200,2005,2010,2015,2020),
                   c(0,1,3,5,10,15,25,60))
  data.hcv$senshcv = 1
  iv.t = rep(-4,8)
  iv.a = c(1,0,0,0,0,0)
  inits2 = c(iv.t,iv.a)
  modc2 = optim(inits2, univarATnoint, method = "BFGS", y = data.hcv$hcv,
                n = data.hcv$n, s = data.hcv$senshcv, data = data.hcv,
                control = list(trace = T, maxit = 500), hessian = T)
  modc2.params = getparamATnoint(modc2)
  modc2.pred = pred(modc2.params$h, data.hcv, data.hcv$hcv, data.hcv$n,
                    data.hcv$senshcv)
  hcv.pred = cbind(data.hcv$year, data.hcv$injdur, data.hcv$n, data.hcv$senshcv)
  colnames(hcv.pred) = c("year", "injdur", "n", "s")
  hcv.pred = cbind(hcv.pred, modc2.pred)
  hcv.pred = as.data.frame(hcv.pred)
  
  #This part is to ensure the data frame sizes match the full data
  full.modc2.pred = data.frame(year = c(rep(2011,44), rep(2012,43),
                                        rep(2013,44), rep(2014,45),
                                        rep(2015,48), rep(2016,45),
                                        rep(2017,48), rep(2018,47),
                                        rep(2019,47), rep(2011,45)),
                               injdur = c(1,10,11,12,13,14,15,16,17,18,19,2,20,
                                          21,22,23,24,25,26,27,28,29,3,30,31,32,
                                          33,34,35,36,38,4,40,41,44,45,48,5,51,
                                          52,6,7,8,9,1,10,11,12,13,14,15,16,17,
                                          18,19,2,20,21,22,23,24,25,26,27,28,29,
                                          3,30,31,32,33,34,35,36,37,38,39,4,40,
                                          41,42,43,5,6,7,8,9,1,10,11,12,13,14,
                                          15,16,17,18,19,2,20,21,22,23,24,25,26,
                                          27,28,29,3,30,31,32,33,34,35,36,37,38,
                                          39,4,40,42,43,45,47,5,6,7,8,9,1,10,11,
                                          12,13,14,15,16,17,18,19,2,20,21,22,23,
                                          24,25,26,27,28,29,3,30,31,32,33,34,35,
                                          36,37,38,4,40,41,42,43,47,5,51,52,6,7,
                                          8,9,1,10,11,12,13,14,15,16,17,18,19,2,
                                          20,21,22,23,24,25,26,27,28,29,3,30,31,
                                          32,33,34,35,36,37,38,39,4,40,41,42,44,
                                          45,47,49,5,50,51,6,7,8,9,1,10,11,12,
                                          13,14,15,16,17,18,19,2,20,21,22,23,24,
                                          25,26,27,28,29,3,30,31,32,33,34,35,36,
                                          37,38,39,4,40,42,43,44,45,46,5,6,7,8,
                                          9,1,10,11,12,13,14,15,16,17,18,19,2,
                                          20,21,22,23,24,25,26,27,28,29,3,30,31,
                                          32,33,34,35,36,37,38,39,4,40,41,42,43,
                                          44,45,46,48,5,52,6,7,8,9,1,10,11,12,
                                          13,14,15,16,17,18,19,2,20,21,22,23,24,
                                          25,26,27,28,29,3,30,31,32,33,34,35,36,
                                          37,38,39,4,40,41,42,43,45,46,47,5,50,
                                          6,7,8,9,1,10,11,12,13,14,15,16,17,18,
                                          19,2,20,21,22,23,24,25,26,27,28,29,3,
                                          30,31,32,33,34,35,36,37,38,39,4,40,41,
                                          42,43,44,46,47,5,53,6,7,8,9,1,10,11,
                                          12,13,14,15,16,17,18,19,2,20,21,22,23,
                                          24,25,26,27,28,29,3,30,31,32,33,34,35,
                                          36,37,38,39,4,40,41,42,44,47,5,50,6,7,
                                          8,9))
  #Empty frame of model predictions
  full.modc2.pred$modc2.pred = NA
  
  #Merge with the predictions from the bootstrapped data set
  merged_df = merge(full.modc2.pred, hcv.pred, by = c("year", "injdur"),
                    all.x = T)
  merged_df$modc2.pred.x[is.na(merged_df$modc2.pred.x)] =
    merged_df$modc2.pred.y[is.na(merged_df$modc2.pred.x)]
  results_df = merged_df[, c("year", "injdur", "modc2.pred.x")]
  colnames(results_df) = c("year", "injdur", "modc2.pred")
  results_df = results_df[order(results_df$injdur, results_df$year), ]
  return(modc2.params$h)
}

#Bootstrapping with 1000 replications
h_ci = boot(data = uam, statistic = h_CI, R = 1000)

#Order data frame to match values
hcv.pred = hcv.pred[order(hcv.pred$injdur, hcv.pred$year), ]

#Append the correct hazard values
hcv.pred = hcv.pred %>%
  mutate(hazard = case_when(year <= 2015 & injdur == 1 ~
                              boot.ci(h_ci, type = "perc", index = 7)$t0,
                            year <= 2015 & injdur > 1 & injdur <= 3 ~
                              boot.ci(h_ci, type = "perc", index = 15)$t0,
                            year <= 2015 & injdur > 3 & injdur <= 5 ~
                              boot.ci(h_ci, type = "perc", index = 23)$t0,
                            year <= 2015 & injdur > 5 & injdur <= 10 ~
                              boot.ci(h_ci, type = "perc", index = 31)$t0,
                            year <= 2015 & injdur > 10 & injdur <= 15 ~
                              boot.ci(h_ci, type = "perc", index = 39)$t0,
                            year <= 2015 & injdur > 15 & injdur <= 25 ~
                              boot.ci(h_ci, type = "perc", index = 47)$t0,
                            year <= 2015 & injdur > 25 & injdur <= 60 ~
                              boot.ci(h_ci, type = "perc", index = 55)$t0,
                            year > 2015 & injdur == 1 ~
                              boot.ci(h_ci, type = "perc", index = 8)$t0,
                            year > 2015 & injdur > 1 & injdur <= 3 ~
                              boot.ci(h_ci, type = "perc", index = 16)$t0,
                            year > 2015 & injdur > 3 & injdur <= 5 ~
                              boot.ci(h_ci, type = "perc", index = 24)$t0,
                            year > 2015 & injdur > 5 & injdur <= 10 ~
                              boot.ci(h_ci, type = "perc", index = 32)$t0,
                            year > 2015 & injdur > 10 & injdur <= 15 ~
                              boot.ci(h_ci, type = "perc", index = 40)$t0,
                            year > 2015 & injdur > 15 & injdur <= 25 ~
                              boot.ci(h_ci, type = "perc", index = 48)$t0,
                            year > 2015 & injdur > 25 & injdur <= 60 ~
                              boot.ci(h_ci, type = "perc", index = 56)$t0),
         lower = case_when(year <= 2015 & injdur == 1 ~
                             boot.ci(h_ci, type = "perc", index = 7)$percent[4],
                           year <= 2015 & injdur > 1 & injdur <= 3 ~
                             boot.ci(h_ci, type = "perc", index = 15)$percent[4],
                           year <= 2015 & injdur > 3 & injdur <= 5 ~
                             boot.ci(h_ci, type = "perc", index = 23)$percent[4],
                           year <= 2015 & injdur > 5 & injdur <= 10 ~
                             boot.ci(h_ci, type = "perc", index = 31)$percent[4],
                           year <= 2015 & injdur > 10 & injdur <= 15 ~
                             boot.ci(h_ci, type = "perc", index = 39)$percent[4],
                           year <= 2015 & injdur > 15 & injdur <= 25 ~
                             boot.ci(h_ci, type = "perc", index = 47)$percent[4],
                           year <= 2015 & injdur > 25 & injdur <= 60 ~
                             boot.ci(h_ci, type = "perc", index = 55)$percent[4],
                           year > 2015 & injdur == 1 ~
                             boot.ci(h_ci, type = "perc", index = 8)$percent[4],
                           year > 2015 & injdur > 1 & injdur <= 3 ~
                             boot.ci(h_ci, type = "perc", index = 16)$percent[4],
                           year > 2015 & injdur > 3 & injdur <= 5 ~
                             boot.ci(h_ci, type = "perc", index = 24)$percent[4],
                           year > 2015 & injdur > 5 & injdur <= 10 ~
                             boot.ci(h_ci, type = "perc", index = 32)$percent[4],
                           year > 2015 & injdur > 10 & injdur <= 15 ~
                             boot.ci(h_ci, type = "perc", index = 40)$percent[4],
                           year > 2015 & injdur > 15 & injdur <= 25 ~
                             boot.ci(h_ci, type = "perc", index = 48)$percent[4],
                           year > 2015 & injdur > 25 & injdur <= 60 ~
                             boot.ci(h_ci, type = "perc", index = 56)$percent[4]),
         upper = case_when(year <= 2015 & injdur == 1 ~
                             boot.ci(h_ci, type = "perc", index = 7)$percent[5],
                           year <= 2015 & injdur > 1 & injdur <= 3 ~
                             boot.ci(h_ci, type = "perc", index = 15)$percent[5],
                           year <= 2015 & injdur > 3 & injdur <= 5 ~
                             boot.ci(h_ci, type = "perc", index = 23)$percent[5],
                           year <= 2015 & injdur > 5 & injdur <= 10 ~
                             boot.ci(h_ci, type = "perc", index = 31)$percent[5],
                           year <= 2015 & injdur > 10 & injdur <= 15 ~
                             boot.ci(h_ci, type = "perc", index = 39)$percent[5],
                           year <= 2015 & injdur > 15 & injdur <= 25 ~
                             boot.ci(h_ci, type = "perc", index = 47)$percent[5],
                           year <= 2015 & injdur > 25 & injdur <= 60 ~
                             boot.ci(h_ci, type = "perc", index = 55)$percent[5],
                           year > 2015 & injdur == 1 ~
                             boot.ci(h_ci, type = "perc", index = 8)$percent[5],
                           year > 2015 & injdur > 1 & injdur <= 3 ~
                             boot.ci(h_ci, type = "perc", index = 16)$percent[5],
                           year > 2015 & injdur > 3 & injdur <= 5 ~
                             boot.ci(h_ci, type = "perc", index = 24)$percent[5],
                           year > 2015 & injdur > 5 & injdur <= 10 ~
                             boot.ci(h_ci, type = "perc", index = 32)$percent[5],
                           year > 2015 & injdur > 10 & injdur <= 15 ~
                             boot.ci(h_ci, type = "perc", index = 40)$percent[5],
                           year > 2015 & injdur > 15 & injdur <= 25 ~
                             boot.ci(h_ci, type = "perc", index = 48)$percent[5],
                           year > 2015 & injdur > 25 & injdur <= 60 ~
                             boot.ci(h_ci, type = "perc", index = 56)$percent[5]))

#Plot the FOI
tiff("pc_FOI.tiff", units="in", width=9, height=7, res=300)
ggplot(hcv.pred, aes(x = injdur)) +
  geom_line(aes(y = hazard), size = 1, color = "green") +
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2) +
  xlab("Injecting Duration (years)") +
  ylab("Force of Infection") +
  facet_wrap(.~year, ncol=4, scales = "free_y")
dev.off()

#Obtain years of interest for manuscript
x3 = hcv.pred %>%
  subset(hcv.pred$year == 2011 |hcv.pred$year == 2014 |
           hcv.pred$year == 2017 |hcv.pred$year == 2019)

#Plot the FOI
#tiff("fp_sub_FOI.tiff", units="in", width=9, height=7, res=300)
ggplot(x3, aes(x = injdur)) +
  geom_line(aes(y = hazard), size = 1, color = "green") +
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2) +
  xlab("Injecting Duration") +
  ylab("Force of Infection") +
  facet_wrap(.~year, ncol=2, scales = "free_y")
#dev.off()
