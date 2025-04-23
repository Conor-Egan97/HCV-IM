set.seed(123)

#Total years and samples per year
years = 2011:2020
samples_per_year = 1000

#Create base data
uam = do.call(rbind, lapply(years, function(y) {
  injdur = sample(0:30, samples_per_year, replace = TRUE)
  
  #Simulate probability of HCV positivity based on injecting duration
  prob_hcv = plogis(-2 + 0.1 * injdur)  # logistic function, ~10% at injdur=0 up to ~90% at injdur=30
  hcvdbs = rbinom(samples_per_year, 1, prob = prob_hcv)
  
  data.frame(year = y, injdur = injdur, hcvdbs = hcvdbs)
}))

#View a snapshot
head(uam)
table(uam$year)
summary(uam$injdur)
prop.table(table(uam$hcvdbs))

#Optional: save it if needed
#write.csv(uam, "simulated_uam_data.csv", row.names = FALSE)
