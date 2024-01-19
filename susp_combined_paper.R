#Load in the relevant packages
library(magrittr)
library(dplyr)
library(tidyverse)
library(lmtest)
library(ciTools)
library(ResourceSelection)
library(scatterplot3d)
library(berryFunctions)
library(haven)
library(Hmisc)
library(rms)
library(ggplot2)
library(stats4)
library(bbmle)
library(rAmCharts)
library(LearnClust)
library(boot)

######################Data preparation#########################################

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

#Create a variable that indicates survival i.e., opposite of HCV dbs status
uam$surv = ifelse(uam$hcvdbs == 1, 0, 1)

######################Fractional polynomial#####################################
#Create data frame of variables of interest - order by injecting duration
y = uam$surv[order(uam$injdur)]
a = uam$injdur[order(uam$injdur)]
t = uam$year[order(uam$injdur)]
neg = table(y, round(a), t)[1, , ]
neg = as.vector(t(neg))
pos = table(y, round(a), t)[2, , ]
pos = as.vector(t(pos))
tot = neg + pos
df1 = data.frame(neg,pos,tot)
df1 = subset(df1, tot!=0)
df = unique(data.frame(a,t))

#Best fitting fractional polynomial
model.fpIJ_T = glm(formula = y ~ I(a^(0.4)) + I(a^(2.8)) +
                     I(a^(0.4)):I(t^(0.3)) + I(a^(2.8)):I(t^(0.4)) - 1,
                   family = binomial(link = log))

#Survival probability predictions
p.at = predict(model.fpIJ_T, type = "response", se.fit=T)

#Store in data frame - will combine with other model estimates below
x1 = data.frame(age = df$a, time = df$t, prob.fp = unique(p.at$fit),
                lower.fp = unique(p.at$fit) - 1.96*unique(p.at$se.fit),
                upper.fp = unique(p.at$fit) + 1.96*unique(p.at$se.fit),
                surv = df1$pos/df1$tot, n = df1$tot)

#######################Natural cubic spline#####################################
#Centered this on 2015
t = t-median(t)

#Number of knots for best fitting model in injecting duration and time dims
knots_q = as.vector(quantile(a, c(0.05, 0.23, 0.41, 0.59, 0.77, 0.95)))
knots_q = c(0, knots_q)
knots_qt = as.vector(quantile(t, c(0.1, 0.5, 0.9)))

#Fit the best fitting NCS
fit.rcs.at = glm(y ~ (rcs(a, knots_q) + rcs(a, knots_q)%ia%rcs(t, knots_qt)) - 1,
                 family = binomial(link = "log"))

#Predicted survival probability
p.at = predict(fit.rcs.at, type = "response", se.fit=T)

#Append the results to the table
x1$prob.ncs = unique(p.at$fit)
x1$lower.ncs = unique(p.at$fit) - 1.96*unique(p.at$se.fit)
x1$upper.ncs = unique(p.at$fit) + 1.96*unique(p.at$se.fit)

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
  return(results_df$modc2.pred)
}

#Bootstrapping with 1000 replications
h_ci = boot(data = uam, statistic = h_CI, R = 1000)

#Store the estimated values and CIs in the data frame
for (i in 1:length(modc2.pred)){
  x1$prob.pc[i] = h_ci$t0[i]
  x1$lower.pc[i] = h_ci$t0[i] - 1.96*apply(h_ci$t, 2, sd)[i]
  x1$upper.pc[i] = h_ci$t0[i] + 1.96*apply(h_ci$t, 2, sd)[i]
}

##############################Joint plots#######################################
#Full plot
tiff("susp_3models.tiff", units="in", width=9, height=7, res=300)
ggplot(data = x1, aes(x=age)) +
  geom_point(aes(y=surv, size=n), color = "black")+
  geom_line(aes(y=prob.fp), color="midnightblue", linetype="solid")+
  geom_ribbon(aes(ymin=lower.fp, ymax=upper.fp), alpha=0.2)+
  geom_line(aes(y=prob.ncs), color="red", linetype="solid")+
  geom_ribbon(aes(ymin=lower.ncs, ymax=upper.ncs), alpha=0.2)+
  geom_line(aes(y=1-prob.pc), color="green", linetype="solid")+
  geom_ribbon(aes(ymin=1-lower.pc, ymax=1-upper.pc), alpha=0.2)+
  xlab("Injecting Duration (years)")+
  ylab("Susceptible")+
  facet_wrap(~ time, nrow = 3)
dev.off()

#Obtain years of interest for manuscript
x3 = x1 %>%
  subset(x1$time == 2011 |x1$time == 2014 |x1$time == 2017 |x1$time == 2019)

tiff("sub_susp_3models.tiff", units="in", width=9, height=7, res=300)
ggplot(data = x3, aes(x=age)) +
  geom_point(aes(y=surv, size=n), color = "black")+
  geom_line(aes(y=prob.fp), color="midnightblue", linetype="solid")+
  geom_ribbon(aes(ymin=lower.fp, ymax=upper.fp), alpha=0.2)+
  geom_line(aes(y=prob.ncs), color="red", linetype="solid")+
  geom_ribbon(aes(ymin=lower.ncs, ymax=upper.ncs), alpha=0.2)+
  geom_line(aes(y=1-prob.pc), color="green", linetype="solid")+
  geom_ribbon(aes(ymin=1-lower.pc, ymax=1-upper.pc), alpha=0.2)+
  xlab("Injecting Duration (years)")+
  ylab("Susceptible")+
  facet_wrap(~ time, nrow = 2)
dev.off()

