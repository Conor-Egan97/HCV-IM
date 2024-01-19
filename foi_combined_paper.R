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


#Set the working directory to the folder with all the FOI functions
setwd("C:/Users/Conor.Egan/Documents")
uam = read_dta('Conor Egan dataset_13092022.dta')

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

set.seed(147)

#Make a data frame for boot CI
df_1 = data.frame(y, a, t)

#Make a function that returns the table of plausible values
FOI_CI = function(data, indices) {
  d = data[indices,] # allows boot to select sample
  df  = as.data.frame(expand.grid(1:53, 2011:2020))
  names(df) = c("a", "t")
  fit <- glm(formula = d$y ~ I(d$a^(0.4)) + I(d$a^(2.8)) +
               I(d$a^(0.4)):I(d$t^0.3) + I(d$a^(2.8)):I(d$t^0.4) - 1,
             family = binomial(link = log))
  foi = as.numeric(coef(fit)[1])*0.4*((df$a)^(-0.6)) +
    as.numeric(coef(fit)[2])*2.8*((df$a)^(1.8)) +
    as.numeric(coef(fit)[3])*0.4*((df$a)^(-0.6))*((df$t)^(0.3)) +
    as.numeric(coef(fit)[4])*2.8*((df$a)^(1.8))*((df$t)^(0.4)) 
  foi = -foi
  return(foi)
}

# bootstrapping with 1000 replications
results = boot(data=df_1, statistic=FOI_CI, R=1000)

#Make a table with the plausible FOI values
FOI_table = as.data.frame(expand.grid(unique(a), unique(t)))

for (i in 1:530){
  FOI_table$foi[i] = boot.ci(results, type="perc", index = i)$t0
  FOI_table$foi_LB[i] = boot.ci(results, type="perc", index = i)$percent[4]
  FOI_table$foi_UB[i] = boot.ci(results, type="perc", index = i)$percent[5]
}

names(FOI_table) = c("injdur", "Year", "FOI", "LB", "UB")

#Obtain years of interest for manuscript
x3 = FOI_table %>%
  subset(FOI_table$Year == 2011 |FOI_table$Year == 2014 |
           FOI_table$Year == 2017 |FOI_table$Year == 2019)

#Plot the FOI
#tiff("fp_sub_FOI.tiff", units="in", width=9, height=7, res=300)
fp_foi = ggplot(data=x3, aes(x=injdur, y=FOI)) +
  geom_line(size = 1, color = "midnightblue") +
  geom_ribbon(aes(ymin=LB, ymax=UB), alpha = 0.2) +
  xlab("Injecting Duration (years)") +
  ylab("Force of infection") +
  scale_x_continuous(limits = c(0, 38))+
  facet_wrap(~ Year, nrow = 3)
#dev.off()

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

#Get the force of infection from this model
set.seed(477)

#Make a data frame for boot CI
df_3 = data.frame(y, a, t)

#Make a function that returns the table of plausible values
IR_CI = function(data, indices) {
  d <- data[indices,]
  df  = as.data.frame(expand.grid(1:53, -4:5))
  names(df) = c("a", "t")
  fit = glm(d$y ~ (rcs(d$a, knots_q) + rcs(d$a,knots_q)%ia%rcs(d$t,knots_qt))-1,
            family=binomial(link = "log"))
  
  #Calculate the other two coefficients
  #First for the linear IJ term
  knot_6_Aspline = (coef(fit)[2]*(0-30)/((30)^2)+
                      coef(fit)[3]*(1-30)/((30)^2)+
                      coef(fit)[4]*(5-30)/((30)^2)+
                      coef(fit)[5]*(11-30)/((30)^2)+
                      coef(fit)[6]*(15-30)/((30)^2))/(30-20) 
  
  knot_7_Aspline = -(coef(fit)[2]/((30)^2)+coef(fit)[3]/((30)^2)+
                       coef(fit)[4]/((30)^2)+coef(fit)[5]/((30)^2)+
                       coef(fit)[6]/((30)^2)+knot_6_Aspline)
  
  #Now for the interaction betwen the time spline and the linear IJ term
  knot_2_A_Tspline = (coef(fit)[8]*(-4-4)/((8)^2))/(4) 
  knot_3_A_Tspline = -(coef(fit)[8]/((8)^2)+knot_2_A_Tspline)
  
  #Now for the interaction betwen the IJ spline and the linear time term
  knot_6_T_Aspline = (coef(fit)[9]*(0-30)/((30)^2)+
                        coef(fit)[10]*(1-30)/((30)^2)+
                        coef(fit)[11]*(5-30)/((30)^2) + 
                        coef(fit)[12]*(11-30)/((30)^2) + 
                        coef(fit)[13]*(15-30)/((30)^2))/(30-20)
  
  knot_7_T_Aspline = -(coef(fit)[9]/((30)^2)+
                         coef(fit)[10]/((30)^2)+
                         coef(fit)[11]/((30)^2)+
                         coef(fit)[12]/((30)^2)+
                         coef(fit)[13]/((30)^2)+
                         knot_6_T_Aspline)
  
  foi_A = (coef(fit)[1]) + (3/((30)^2))*(coef(fit)[2])*(pmax(df$a-0,0)^2) +
    (3/((30)^2))*(coef(fit)[3])*(pmax(df$a-1,0)^2) +
    (3/((30)^2))*(coef(fit)[4])*(pmax(df$a-5,0)^2) +
    (3/((30)^2))*(coef(fit)[5])*(pmax(df$a-11,0)^2) +
    (3/((30)^2))*(coef(fit)[6])*(pmax(df$a-15,0)^2) +
    3*(knot_6_Aspline)*(pmax(df$a-20,0)^2) +
    3*(knot_7_Aspline)*(pmax(df$a-30,0)^2)
  
  foi_A_Tspline = (coef(fit)[7])*df$t +
    (1/((8)^2))*(coef(fit)[8])*(pmax(df$t+4,0)^3) +
    (knot_2_A_Tspline)*(pmax(df$t-0,0)^3) +
    (knot_3_A_Tspline)*(pmax(df$t-4,0)^3)
  
  foi_T_Aspline = (df$t)*((3/((30)^2))*(coef(fit)[9])*(pmax(df$a-0,0)^2)  +
                            (3/((30)^2))*(coef(fit)[10])*(pmax(df$a-1,0)^2) +
                            (3/((30)^2))*(coef(fit)[11])*(pmax(df$a-5,0)^2) +
                            (3/((30)^2))*(coef(fit)[12])*(pmax(df$a-11,0)^2) +
                            (3/((30)^2))*(coef(fit)[13])*(pmax(df$a-15,0)^2) +
                            3*(knot_6_T_Aspline)*(pmax(df$a-20,0)^2) +
                            3*(knot_7_T_Aspline)*(pmax(df$a-30,0)^2))
  
  foi = foi_A + foi_A_Tspline + foi_T_Aspline
  foi = -foi
  return(foi)
}

# bootstrapping with 1000 replications
results = boot(data=df_3, statistic=IR_CI, R=1000)

#Make a table with the plausible FOI values
IR_table = as.data.frame(expand.grid(unique(a), unique(t)))

for (i in 1:530){
  IR_table$IR[i] = boot.ci(results, type="perc", index = i)$t0
  IR_table$IR_LB[i] = boot.ci(results, type="perc", index = i)$percent[4]
  IR_table$IR_UB[i] = boot.ci(results, type="perc", index = i)$percent[5]
}

IR_table$Var2 = IR_table$Var2 + 2015

names(IR_table) = c("injdur", "Year", "IR", "LB", "UB")

#Obtain years of interest for manuscript
x3 = IR_table %>%
  subset(IR_table$Year == 2011 |IR_table$Year == 2014 |
           IR_table$Year == 2017 |IR_table$Year == 2019)

#Plot the FOI
ncs_foi = ggplot(data=x3, aes(x=injdur, y=IR)) +
  geom_line(size = 1, color = "red") +
  geom_ribbon(aes(ymin=LB, ymax=UB), alpha = 0.2) +
  xlab("Injecting Duration (years)") +
  ylab("Incidence rate") +
  facet_wrap(~ Year, nrow = 2)
#dev.off()

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

#Obtain years of interest for manuscript
x3 = hcv.pred %>%
  subset(hcv.pred$year == 2011 |hcv.pred$year == 2014 |
           hcv.pred$year == 2017 |hcv.pred$year == 2019)

#Plot the FOI
#tiff("fp_sub_FOI.tiff", units="in", width=9, height=7, res=300)
pc_foi = ggplot(x3, aes(x = injdur)) +
  geom_line(aes(y = hazard), size = 1, color = "green") +
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2) +
  xlab("Injecting Duration (years)") +
  ylab("Force of Infection") +
  facet_wrap(.~year, ncol=2, scales = "free_y")
#dev.off()

##############################Joint plots#######################################
library(gridExtra)

#Full plot
#tiff("foi_3models.tiff", units="in", width=9, height=7, res=300)
grid.arrange(pc_foi, fp_foi, ncs_foi, nrow = 2)
#dev.off()

