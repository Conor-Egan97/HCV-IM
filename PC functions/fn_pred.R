#params <- mod3.params
#data <- data.hcv

pred <- function(params,data,y,n,s){
  sumht <- rep(0,length(y))
  pi_tru0 <- sumht
  ht<-params
  sumht <- sumht +ht[1,1]*data$x1_1+  ht[1,2]*data$x1_2+  ht[1,3]*data$x1_3+  ht[1,4]*data$x1_4+  ht[1,5]*data$x1_5+  ht[1,6]*data$x1_6+  ht[1,7]*data$x1_7
  sumht <- sumht +ht[2,1]*data$x2_1+  ht[2,2]*data$x2_2+  ht[2,3]*data$x2_3+  ht[2,4]*data$x2_4+  ht[2,5]*data$x2_5+  ht[2,6]*data$x2_6+  ht[2,7]*data$x2_7
  sumht <- sumht +ht[3,1]*data$x3_1+  ht[3,2]*data$x3_2+  ht[3,3]*data$x3_3+  ht[3,4]*data$x3_4+  ht[3,5]*data$x3_5+  ht[3,6]*data$x3_6+  ht[3,7]*data$x3_7
  sumht <- sumht +ht[4,1]*data$x4_1+  ht[4,2]*data$x4_2+  ht[4,3]*data$x4_3+  ht[4,4]*data$x4_4+  ht[4,5]*data$x4_5+  ht[4,6]*data$x4_6+  ht[4,7]*data$x4_7
  sumht <- sumht +ht[5,1]*data$x5_1+  ht[5,2]*data$x5_2+  ht[5,3]*data$x5_3+  ht[5,4]*data$x5_4+  ht[5,5]*data$x5_5+  ht[5,6]*data$x5_6+  ht[5,7]*data$x5_7
  sumht <- sumht +ht[6,1]*data$x6_1+  ht[6,2]*data$x6_2+  ht[6,3]*data$x6_3+  ht[6,4]*data$x6_4+  ht[6,5]*data$x6_5+  ht[6,6]*data$x6_6+  ht[6,7]*data$x6_7
  sumht <- sumht +ht[7,1]*data$x7_1+  ht[7,2]*data$x7_2+  ht[7,3]*data$x7_3+  ht[7,4]*data$x7_4+  ht[7,5]*data$x7_5+  ht[7,6]*data$x7_6+  ht[7,7]*data$x7_7
  sumht <- sumht +ht[8,1]*data$x8_1+  ht[8,2]*data$x8_2+  ht[8,3]*data$x8_3+  ht[8,4]*data$x8_4+  ht[8,5]*data$x8_5+  ht[8,6]*data$x8_6+  ht[8,7]*data$x8_7
  
  pi_tru0 <- exp(-sumht[])
  pi_tru1 <- 1 - pi_tru0
  pi1 <- s*pi_tru1
  pi0 <- 1-pi1
  
  return(pi1)
}

