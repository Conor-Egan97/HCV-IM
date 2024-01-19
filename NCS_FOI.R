#Load in libraries
library(magrittr)
library(tidyverse)
library(rms)
library(Hmisc)
library(ggplot2)
library(ResourceSelection)
library(boot)
library(haven)

######################Data preparation#########################################
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

tiff("NCS_injdur_time_FOI.tiff", units="in", width=9, height=7, res=300)
ggplot(data=IR_table, aes(x=injdur, y=IR)) +
  geom_line(size = 1, color = "red") +
  geom_ribbon(aes(ymin=LB, ymax=UB), alpha = 0.2) +
  xlab("Injecting Duration (years)") +
  ylab("Incidence rate") +
  facet_wrap(~ Year, nrow = 3)
dev.off()

#Obtain years of interest for manuscript
x3 = IR_table %>%
  subset(IR_table$Year == 2011 |IR_table$Year == 2014 |
           IR_table$Year == 2017 |IR_table$Year == 2019)

#Plot the FOI
ggplot(data=x3, aes(x=injdur, y=IR)) +
  geom_line(size = 1, color = "red") +
  geom_ribbon(aes(ymin=LB, ymax=UB), alpha = 0.2) +
  xlab("Injecting Duration (years)") +
  ylab("Incidence rate") +
  facet_wrap(~ Year, nrow = 2)

#Make a subset of the table that only contains specific age injecting groups
FOI_table_IJ_subset = IR_table %>% 
  subset(IR_table$injdur == 1 | IR_table$injdur == 5 |
           IR_table$injdur == 10 | IR_table$injdur == 15 |
           IR_table$injdur == 20)


#Make injecting duratoin a factor
FOI_table_IJ_subset$injdur = as.factor(FOI_table_IJ_subset$injdur)

#tiff("FP_time_FOI_ALL.tiff", units="in", width=9, height=7, res=300)
ggplot(data=FOI_table_IJ_subset, aes(x=Year, y=IR, colour = injdur)) +
  geom_line() +
  geom_ribbon(aes(ymin=LB, ymax=UB), alpha = 0.2) +
  xlab("Survey year") +
  ylab("Force of infection") +
  scale_x_continuous(breaks=seq(2011,2020,by=1))+
  scale_y_continuous(breaks=seq(0,0.2,by=0.05))+ 
  labs(color='Injecting duration') 
#dev.off()
