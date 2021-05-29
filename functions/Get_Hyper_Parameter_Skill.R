#______________________________________________________________________________#
###Function to compare skill with different hyper-parameters.


###Input
#1. True Data
#2. Simulations - For all hyper-paramaters
#3. Hyper-Parameter Values



###Output
#1. Distribution of D-Scores and p-values - Ggplot2

#Hyper_Par <- c("30 Days", "20 Days", "Base KNN")


#______________________________________________________________________________#

get_hyperparameter_skill <- function(Dat, Com_Sims, Hyper_Par){
  
  #Load Dependencies
  library(ggplot2)
  library(lubridate)
  library(cowplot)
  library(gridExtra)
  
  #Parameters
  N <- length(Com_Sims)
  
  #Data
  Dat <- scale(Dat)
  Dat <- rowMeans(Dat)
  
  #Simulations
  All_D_Stat <- All_P_Val <- list()
  for(i in 1:N){
    
    #Get the current Index  
    temp <- Com_Sims[[i]]  
  
    #Store KS D-Statistic and P-Value
    D_Stat <- P_Val <- list()
  
    for(j in 1:length(temp)){
      tt <- cbind(temp[[j]]$WPnew, temp[[j]]$SSnew)
      tt <- scale(tt)
      tt <- rowMeans(tt)
      ks_test <- ks.test(Dat, tt)
      D_Stat[[j]] <- ks_test$statistic
      P_Val[[j]] <- ks_test$p.value
    }
    
  All_D_Stat[[i]] <- D_Stat
  All_P_Val[[i]] <- P_Val
  }
  
  #Getting the Index
  ind <- list()
  for(i in 1:N){
    ind[[i]] <- rep(Hyper_Par[i], length(Com_Sims[[i]]))
  }
  
  #Convert to Data-Frame Format
  KS_Results <- data.frame(D_Stat = unlist(All_D_Stat),
                           P_Value = unlist(All_P_Val))
  KS_Results$Type <- unlist(ind)
  
  
  #Plotting the Results
  p_main <- ggplot(KS_Results) + geom_point(aes(D_Stat, P_Value, color = Type)) +
    ggtitle("Hyper-Parameter Selection - Two-sample Kolmogorov-Smirnov test") + 
    xlab("D_Stat") + ylab("P_Value") + 
    geom_hline(yintercept = 0.05) +
    theme_bw()
  
  xbox <- axis_canvas(p_main, axis = "x", coord_flip = TRUE) + 
    geom_boxplot(data = KS_Results, aes(y = D_Stat, x = factor(Type), color = factor(Type))) + 
    scale_x_discrete() + coord_flip()
  
  ybox <- axis_canvas(p_main, axis = "y") + 
    geom_boxplot(data = KS_Results, aes(y = P_Value, x = factor(Type), color = factor(Type))) +
    scale_x_discrete()
  
  pnull <- ggdraw() # generate empty plot
  
  p1 <- insert_xaxis_grob(
    insert_xaxis_grob(p_main, xbox, grid::unit(0.6, "in"), position = "top"),
    pnull, grid::unit(0.2, "in"), position = "top")
  
  p2 <- insert_yaxis_grob(
    insert_yaxis_grob(p1, ybox, grid::unit(0.6, "in"), position = "right"),
    pnull, grid::unit(0.2, "in"), position = "right")
  
  p3 <- ggdraw(p2)
  
  print(p3)
  
}
