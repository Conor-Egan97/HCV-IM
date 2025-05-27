## This function calculates for each individual how long they spent in each injecting duration bracket
## revised mkcut function (based on "piecewise constant 01.R" from Analysis 01)
## assign the groupings directly to the data
mkcut = function(data, year, inj, c.year, c.inj){
  dataplus <- data
  #Calculate the number of year intervals
  n.year = length(c.year) - 1
  #Calculate the number of injecting intervals
  n.inj = length(c.inj) - 1
  #This will be the number of groupings as each year is repeated per injecting duration
  n = length(year)
  #This gives the year a person started injecting
  fyear = year-inj
  z = year*0
  #groups <- array(NA, c(n,n.year,n.inj))
  minx = z
  maxx = z
  for(y in 1:n.year){
    for(i in 1:n.inj){
      #c <- (y-1)*n.inj + i
      maxx <- pmin(c.year[y+1]+z, c.inj[i+1]+fyear, year)
      minx <- pmax(c.year[y]+z, c.inj[i]+fyear, fyear)
      #groups[,y,i] <- pmax(0, maxx-minx )
      nam<-paste("x",y,"_",i,sep="")
      dataplus[nam] <- pmax(0, maxx-minx )
    }
  }
  return(dataplus)
}
