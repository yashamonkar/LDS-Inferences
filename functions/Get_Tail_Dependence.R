#______________________________________________________________________________#
#Function to compute and visualize the tail dependences between variables at
#each Site.

#Input
#1. Data Fields
#2. Simulations.
#Note:- Simulations are a list of dataframes.
#Significane of Tail Dependences. 


#Output
#1. Plot of True and Simulated Correlations for Tail Dependences

#Packages
library(dplyr)

Get_Tail_Dependence <- function(Fld1,Fld2, # Data
                                 Fld1_Sims,Fld2_Sims, #Simulations
                                 Level){
  #Hyper-Parameters
  sims <- length(Fld1_Sims)
  threshold <- 1-Level
  
  #Computing True Data Probabilites
  wgs <- sgw <- list() #Wind Given Solar
  for(i in 1:ncol(Fld1)){
    f1 <- Fld1[,i]
    f2 <- Fld2[,i]
    
    f1_thres <- quantile(f1,threshold)
    f2_thres <- quantile(f2,threshold)
    
    f1b <- ifelse(f1 > f1_thres, 1, 0)
    f2b <- ifelse(f2 > f2_thres, 1, 0)
    
    fb <- as.data.frame(cbind(f1b,f2b))
    
    
    #P(Fld1 > Fld1*|Fld2 > Fld2*) - P(Wind|Solar)
    ft <- fb %>% filter(f2b == 1)
    wgs[[i]] <- sum(ft$f1b)/nrow(ft)
    
    #P(Fld2 > Fld2*|Fld1 > Fld1*) - P(Solar|Wind)
    ft <- fb %>% filter(f1b == 1)
    sgw[[i]] <- sum(ft$f2b)/nrow(ft)
  }

  wgs <- unlist(wgs)
  sgw <- unlist(sgw)
  
  ###Computing Simulation - P(Wind|Solar)
  plot(wgs,wgs, 
     main = paste0("P(Wind |Solar ) ",threshold*100,"th"),
     xlab = "True Data", ylab = "Mean Simulated", type='n',
     pch=19, col='red', ylim = c(0,0.5))
  abline(coef = c(0,1), lwd = 2)

  for(i in 1:ncol(Fld1)){
    wgs_sims <- list()
    
    for(j in 1:sims) {
      
      Fld1_temp <- Fld1_Sims[[j]]
      Fld2_temp <- Fld2_Sims[[j]]
      
      f1 <- Fld1_temp[,i]
      f2 <- Fld2_temp[,i]
    
      f1_thres <- quantile(f1,threshold)
      f2_thres <- quantile(f2,threshold)
    
      f1b <- ifelse(f1 > f1_thres, 1, 0)
      f2b <- ifelse(f2 > f2_thres, 1, 0)
    
      fb <- as.data.frame(cbind(f1b,f2b))

      #P(Fld1 > Fld1*|Fld2 > Fld2*) - P(Wind|Solar)
      ft <- fb %>% filter(f2b == 1)
      wgs_sims[[j]] <- sum(ft$f1b)/nrow(ft)
      }
    
  #Plotting the Results
  tq <- unlist(wgs_sims)
  med <- median(tq)
  quants <- c(0.05,0.25,0.75,0.90)
  qt <- quantile(tq, prob = quants)
  
  segments(wgs[i], qt[1], wgs[i], qt[4], col = "red")
  segments(wgs[i], qt[2], wgs[i], qt[3], col = "blue", lwd = 2)
  points(wgs[i],med, col ="black", pch = 19)
  }
  
  
  
  ###Computing Simulation - P(Solar|Wind)
  plot(sgw,sgw, 
       main = paste0("P(Solar |Wind ) ",threshold*100,"th"),
       xlab = "True Data", ylab = "Mean Simulated", type='n',
       pch=19, col='red', ylim = c(0,0.5))
  abline(coef = c(0,1), lwd = 2)
  
  for(i in 1:ncol(Fld1)){
    sgw_sims <- list()
    
    for(j in 1:sims) {
      
      Fld1_temp <- Fld1_Sims[[j]]
      Fld2_temp <- Fld2_Sims[[j]]
      
      f1 <- Fld1_temp[,i]
      f2 <- Fld2_temp[,i]
      
      f1_thres <- quantile(f1,threshold)
      f2_thres <- quantile(f2,threshold)
      
      f1b <- ifelse(f1 > f1_thres, 1, 0)
      f2b <- ifelse(f2 > f2_thres, 1, 0)
      
      fb <- as.data.frame(cbind(f1b,f2b))
      
      #P(Fld2 > Fld2*|Fld1 > Fld1*) - P(Solar|Wind)
      ft <- fb %>% filter(f1b == 1)
      sgw_sims[[j]] <- sum(ft$f2b)/nrow(ft)
    }
    
    #Plotting the Results
    tq <- unlist(sgw_sims)
    med <- median(tq)
    quants <- c(0.05,0.25,0.75,0.90)
    qt <- quantile(tq, prob = quants)
    
    segments(sgw[i], qt[1], sgw[i], qt[4], col = "red")
    segments(sgw[i], qt[2], sgw[i], qt[3], col = "blue", lwd = 2)
    points(sgw[i],med, col ="black", pch = 19)
  }
 
} 
