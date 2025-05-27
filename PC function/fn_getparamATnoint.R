getparamATnoint <- function(model){
  ## code is pasted from model - hence still called "inits"
  inits <- model$par
  npar <- length(inits)
  Bt <- array(0,8)
  Ba <- array(0,7)
  h <- array(0,c(8,7))
  logh <- h
  Bt[] <- inits[1:8]
  Ba[] <- c(inits[9:10],0,inits[11:14])
  
  for(t in 1:8){
    for(a in 1:7){
      h[t,a] <- exp(Bt[t]+Ba[a])
      logh[t,a] <- Bt[t]+Ba[a]
    }
  }
  
  results <- list(NA)
  results$time <- round(exp(Bt[]),digits=3)
  results$age <- round(exp(Ba[]),digits=3)
  results$h <- h[]
  results$Bt <- Bt
  results$Ba <- Ba
  results$logh <- logh
  ## deviance? Maybe a different function...
  results$params <- npar
  return(results)
}
