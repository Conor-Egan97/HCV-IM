library(magrittr)
library(dplyr)
library(tidyverse)
library(mfp)
library(lmtest)
library(ciTools)
library(ResourceSelection)
library(scatterplot3d)
library(berryFunctions)
library(haven)
library(boot)

######################Data preparation#########################################

#Set the working directory to the folder with all the FOI functions
#setwd("C:/Users/Conor.Egan/Documents")
load("/Users/conoregan/Documents/Cambridge/University/Data/uam_censored.Rdata")

uam = uam_censored

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

#tiff("FP_injdur_time_FOI.tiff", units="in", width=9, height=7, res=300)
ggplot(data=FOI_table, aes(x=injdur, y=FOI)) +
  geom_line(size = 1, color = "midnightblue") +
  geom_ribbon(aes(ymin=LB, ymax=UB), alpha = 0.2) +
  xlab("Injecting Duration (years)") +
  ylab("Force of infection") +
  scale_x_continuous(limits = c(0, 38))+
  facet_wrap(~ Year, nrow = 3)
#dev.off()

#Obtain years of interest for manuscript
x3 = FOI_table %>%
  subset(FOI_table$Year == 2011 |FOI_table$Year == 2014 |
           FOI_table$Year == 2017 |FOI_table$Year == 2019)

#Plot the FOI
#tiff("fp_sub_FOI.tiff", units="in", width=9, height=7, res=300)
ggplot(data=x3, aes(x=injdur, y=FOI)) +
  geom_line(size = 1, color = "midnightblue") +
  geom_ribbon(aes(ymin=LB, ymax=UB), alpha = 0.2) +
  xlab("Injecting Duration (years)") +
  ylab("Force of infection") +
  scale_x_continuous(limits = c(0, 38))+
  facet_wrap(~ Year, nrow = 3)
#dev.off()

#Make a subset of the table that only contains specific age injecting groups
FOI_table_IJ_subset = FOI_table %>% 
  subset(FOI_table$injdur == 1 | FOI_table$injdur == 5 |
           FOI_table$injdur == 10 | FOI_table$injdur == 15 |
           FOI_table$injdur == 20)

#Make injecting duratoin a factor
FOI_table_IJ_subset$injdur = as.factor(FOI_table_IJ_subset$injdur)

#tiff("FP_time_FOI_ALL.tiff", units="in", width=9, height=7, res=300)
ggplot(data=FOI_table_IJ_subset, aes(x=Year, y=FOI, colour = injdur)) +
  geom_line() +
  geom_ribbon(aes(ymin=LB, ymax=UB), alpha = 0.2) +
  xlab("Survey year") +
  ylab("Force of infection") +
  scale_x_continuous(breaks=seq(2011,2020,by=1))+
  scale_y_continuous(breaks=seq(0,0.2,by=0.05))+ 
  labs(color='Injecting duration') 
#dev.off()