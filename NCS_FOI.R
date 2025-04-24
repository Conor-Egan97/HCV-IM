#Load required libraries
library(magrittr)
library(tidyverse)
library(rms)
library(boot)

#####################################
#Data Preparation
#####################################

#Load simulated data
uam = read.csv("~/Desktop/simulated_uam_data.csv")

#Filter data for survey years 2011-2020
uam = uam %>%
  filter(year %in% as.character(2011:2020))

#Assume injecting duration of 1 for individuals reporting 0
uam$injdur[uam$injdur == 0] = 1

#Keep complete cases for relevant variables
uam = uam %>%
  select(injdur, hcvdbs, year) %>%
  drop_na()

#Create binary survival variable: 1 = survived (HCV negative), 0 = not
uam$surv = ifelse(uam$hcvdbs == 1, 0, 1)

#Order variables for modeling
y = uam$surv[order(uam$injdur)]
a = uam$injdur[order(uam$injdur)]
t = uam$year[order(uam$injdur)]
t = as.numeric(t) - 2015  #Center time on 2015

#Create table of outcomes by injecting duration and time
neg = as.vector(t(table(y, round(a), t)[1, , ]))
pos = as.vector(t(table(y, round(a), t)[2, , ]))
tot = neg + pos

df1 = data.frame(neg, pos, tot) %>%
  filter(tot != 0)
df = unique(data.frame(a, t))

#####################################
#Natural Cubic Spline Model
#####################################

#Define knots for spline fitting
knots_q = c(0, 1, 5, 11, 15, 20, 30)
knots_qt = quantile(t, probs = c(0.1, 0.5, 0.9))

#Fit the model
fit.rcs.at = glm(y ~ (rcs(a, knots_q) + rcs(a, knots_q) %ia% rcs(t, knots_qt)) - 1,
                  family = binomial(link = "log"))

#Prepare for bootstrapping
set.seed(477)
df_3 = data.frame(y, a, t)

#Function to compute force of infection (FOI) from bootstrapped model
IR_CI = function(data, indices) {
  d = data[indices, ]
  df = expand.grid(a = 1:53, t = -4:5)
  
  fit = glm(d$y ~ (rcs(d$a, knots_q) + rcs(d$a, knots_q) %ia% rcs(d$t, knots_qt)) - 1,
             family = binomial(link = "log"))
  
  #Compute spline-derived coefficients manually
  coef_vals = coef(fit)
  
  #Age spline (individual components and tail terms)
  a_spline = function(a_val) {
    terms = c(0, 1, 5, 11, 15, 20, 30)
    weights = c(coef_vals[2:6], NA, NA)
    weights[6] = (sum(coef_vals[2:6] * (terms[1:5] - 30) / 30^2)) / (30 - 20)
    weights[7] = -sum(c(coef_vals[2:6], weights[6]) / 30^2)
    sum(3 * weights * (pmax(a_val - terms, 0)^2))
  }
  
  #Time spline
  t_spline = function(t_val) {
    s1 = coef_vals[7] * t_val
    s2 = coef_vals[8] * (pmax(t_val + 4, 0)^3) / (8^2)
    s3 = ((coef_vals[8] * (-8) / 64) / 4) * (pmax(t_val, 0)^3)
    s4 = -(coef_vals[8] / 64 + s3) * (pmax(t_val - 4, 0)^3)
    s1 + s2 + s3 + s4
  }
  
  #Interaction spline
  t_spline_a = function(a_val, t_val) {
    terms = c(0, 1, 5, 11, 15, 20, 30)
    weights = c(coef_vals[9:13], NA, NA)
    weights[6] = (sum(coef_vals[9:13] * (terms[1:5] - 30) / 30^2)) / (30 - 20)
    weights[7] = -sum(c(coef_vals[9:13], weights[6]) / 30^2)
    sum(t_val * 3 * weights * (pmax(a_val - terms, 0)^2))
  }
  
  #Final FOI computation
  foi = mapply(function(a_val, t_val) {
    - (coef_vals[1] + a_spline(a_val) + t_spline(t_val) + t_spline_a(a_val, t_val))
  }, df$a, df$t)
  
  return(foi)
}

#Run bootstrapping
results = boot(data = df_3, statistic = IR_CI, R = 1000)

#Construct FOI table
IR_table = expand.grid(unique(a), unique(t)) %>%
  setNames(c("injdur", "Year"))

for (i in 1:nrow(IR_table)) {
  ci = boot.ci(results, type = "perc", index = i)
  IR_table$IR[i] = ci$t0
  IR_table$LB[i] = ci$percent[4]
  IR_table$UB[i] = ci$percent[5]
}

#Convert centered year back to actual year
IR_table$Year = IR_table$Year + 2015

#####################################
#Plotting
#####################################

#Full FOI over time
ggplot(IR_table, aes(x = injdur, y = IR)) +
  geom_line(size = 1, color = "red") +
  geom_ribbon(aes(ymin = LB, ymax = UB), alpha = 0.2) +
  xlab("Injecting Duration (years)") +
  ylab("Incidence rate") +
  facet_wrap(~ Year, nrow = 3)

#Subset for selected years
x3 = IR_table %>%
  filter(Year %in% c(2011, 2014, 2017, 2019))

ggplot(x3, aes(x = injdur, y = IR)) +
  geom_line(size = 1, color = "red") +
  geom_ribbon(aes(ymin = LB, ymax = UB), alpha = 0.2) +
  xlab("Injecting Duration (years)") +
  ylab("Incidence rate") +
  facet_wrap(~ Year, nrow = 2)

#Subset for specific injecting durations
FOI_table_IJ_subset = IR_table %>%
  filter(injdur %in% c(1, 5, 10, 15, 20)) %>%
  mutate(injdur = as.factor(injdur))

ggplot(FOI_table_IJ_subset, aes(x = Year, y = IR, colour = injdur)) +
  geom_line() +
  geom_ribbon(aes(ymin = LB, ymax = UB), alpha = 0.2) +
  xlab("Survey year") +
  ylab("Force of infection") +
  scale_x_continuous(breaks = 2011:2020) +
  scale_y_continuous(breaks = seq(0, 0.2, by = 0.05)) +
  labs(color = "Injecting duration")
