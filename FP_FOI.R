#Load necessary packages
library(dplyr)
library(boot)
library(ggplot2)

#------------------ Data Preparation ------------------

#Load simulated data
uam = read.csv("~/Desktop/simulated_uam_data.csv")

#Filter for survey years 2011–2020
uam = uam %>% filter(year %in% as.character(2011:2020))

#Replace 0 injecting duration with 1
uam$injdur[uam$injdur == 0] = 1

#Keep only complete cases for relevant variables
uam = uam %>% select(injdur, hcvdbs, year) %>% na.omit()

#Create survival variable: 1 = HCV negative, 0 = HCV positive
uam$surv = ifelse(uam$hcvdbs == 1, 0, 1)

#------------------ Fractional Polynomial Modeling ------------------

#Sort data by injecting duration
y = uam$surv[order(uam$injdur)]
a = uam$injdur[order(uam$injdur)]
t = uam$year[order(uam$injdur)]

#Construct outcome tables
neg = as.vector(t(table(y, round(a), t)[1, , ]))
pos = as.vector(t(table(y, round(a), t)[2, , ]))
tot = neg + pos
df1 = data.frame(neg, pos, tot) %>% filter(tot != 0)
df = unique(data.frame(a, t))

#Fit fractional polynomial model
model.fpIJ_T = glm(
  y ~ I(a^0.4) + I(a^2.8) + I(a^0.4):I(t^0.3) + I(a^2.8):I(t^0.4) - 1,
  family = binomial(link = "log")
)

#Predict survival probabilities
p.at = predict(model.fpIJ_T, type = "response", se.fit = TRUE)

#Create predictions data frame
x1 = data.frame(
  age = df$a,
  time = df$t,
  prob.fp = unique(p.at$fit),
  lower.fp = unique(p.at$fit) - 1.96 * unique(p.at$se.fit),
  upper.fp = unique(p.at$fit) + 1.96 * unique(p.at$se.fit),
  surv = df1$pos / df1$tot,
  n = df1$tot
)

#------------------ Bootstrap for Force of Infection ------------------

set.seed(147)
df_1 = data.frame(y, a, t)

#Function to compute FOI from resampled data
FOI_CI = function(data, indices) {
  d = data[indices, ]
  df = expand.grid(a = 1:30, t = 2011:2020)
  
  fit = tryCatch({
    glm(
      formula = d$y ~ I(d$a^0.4) + I(d$a^2.8) +
        I(d$a^0.4):I(d$t^0.3) + I(d$a^2.8):I(d$t^0.4) - 1,
      family = binomial(link = "log")
    )
  }, error = function(e) return(NULL))
  
  if (is.null(fit)) return(rep(NA, nrow(df)))
  
  coefs = coef(fit)
  if (any(is.na(coefs))) return(rep(NA, nrow(df)))
  
  foi = coefs[1] * 0.4 * (df$a^-0.6) +
    coefs[2] * 2.8 * (df$a^1.8) +
    coefs[3] * 0.4 * (df$a^-0.6) * (df$t^0.3) +
    coefs[4] * 2.8 * (df$a^1.8) * (df$t^0.4)
  
  return(-foi)
}

#Bootstrap FOI estimates
results = boot(data = df_1, statistic = FOI_CI, R = 1000)

#Create table of plausible FOI values
FOI_table = expand.grid(injdur = unique(a), Year = unique(t))

n_rows = nrow(FOI_table)
FOI_table$FOI = NA
FOI_table$LB = NA
FOI_table$UB = NA

for (i in 1:n_rows) {
  ci = tryCatch(boot.ci(results, type = "perc", index = i), error = function(e) return(NULL))
  if (!is.null(ci) && !is.null(ci$percent)) {
    FOI_table$FOI[i] = ci$t0
    FOI_table$LB[i] = ci$percent[4]
    FOI_table$UB[i] = ci$percent[5]
  }
}

#------------------ Plotting ------------------

#Full FOI vs injecting duration plot by year
ggplot(FOI_table, aes(x = injdur, y = FOI)) +
  geom_line(size = 1, color = "midnightblue") +
  geom_ribbon(aes(ymin = LB, ymax = UB), alpha = 0.2) +
  xlab("Injecting Duration (years)") +
  ylab("Force of Infection") +
  scale_x_continuous(limits = c(0, 38)) +
  facet_wrap(~Year, nrow = 3)

#Subset of years of interest
x3 = FOI_table %>% filter(Year %in% c(2011, 2014, 2017, 2019))

#FOI plot for selected years
ggplot(x3, aes(x = injdur, y = FOI)) +
  geom_line(size = 1, color = "midnightblue") +
  geom_ribbon(aes(ymin = LB, ymax = UB), alpha = 0.2) +
  xlab("Injecting Duration (years)") +
  ylab("Force of Infection") +
  scale_x_continuous(limits = c(0, 38)) +
  facet_wrap(~Year, nrow = 3)

#Subset for specific injecting durations
FOI_table_IJ_subset = FOI_table %>%
  filter(injdur %in% c(1, 5, 10, 15, 20)) %>%
  mutate(injdur = as.factor(injdur))

#FOI over time by injecting duration
ggplot(FOI_table_IJ_subset, aes(x = Year, y = FOI, colour = injdur)) +
  geom_line() +
  geom_ribbon(aes(ymin = LB, ymax = UB), alpha = 0.2) +
  xlab("Survey Year") +
  ylab("Force of Infection") +
  scale_x_continuous(breaks = 2011:2020) +
  scale_y_continuous(breaks = seq(0, 0.2, by = 0.05)) +
  labs(color = "Injecting Duration")
