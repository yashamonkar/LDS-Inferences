#______________________________________________________________________________#
###Code to implement the K-S Test.
#The K-S test is fit to the aggregated values. 
#The idea is to fit to both wind and solar individually. 
#Note:- This takes either solar or wind.
#If combined fields are to be used must be aggregated initially.


###INPUT
#1. Data
#2. Simulations
#3. Type


###OUTPUT
#1. D-Statistic
#2. P-Values


#______________________________________________________________________________#
get_KS_Test <- function(Dat, Sims, Type){
  
  #Hyper-Parameters
  nsim <- length(Sims)
  
  #Data-Wranging
  Agg_Dat <- rowSums(Dat)

  
  #Kolmogorov Smirnov Tests
  D_Stat <- P_Val <- list()
  for(i in 1:nsim){
    Agg_Sim <- rowSums(Sims[[i]])
    KS <- ks.test(Agg_Dat, Agg_Sim)
    D_Stat[[i]] <- KS$statistic
    P_Val[[i]] <- KS$p.value
  }
  
  #Output in Data-Form
  KS_Results <- data.frame(D_Stat = unlist(D_Stat),
                           P_Value = unlist(P_Val),
                           Type = rep(Type, nsim))
  
  return(KS_Results)
  
  
}