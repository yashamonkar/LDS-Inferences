#______________________________________________________________________________#
#Function to plot Site Correlation between Data and Simulations. 

#Input
#1. Data Fields
#2. Simulations.
#3. Grid - Grid Box
#Note:- Simulations are a list of dataframes.


#Output
#1. Plot of True and Simulated Correlations at each Site. 
#2. Stick Plot for 1. 
#3. Visual Check for Simulations for both fields - True Data. 


Get_Site_Correlation <- function(Fld1,Fld2, # Data
                                 Fld1_Sims,Fld2_Sims, #Simulations
                                 Grid){
  

  library(ggplot2)
  library(usmap)
  
  #Hyper-Parameters
  sims <- length(Fld1_Sims)


  #______________________________________________________________________________# 
  #Correlation for True Data
  M <- rep(NA, ncol(Fld1))
  for(i in 1:ncol(Fld1)){ 
    M[i] <- cor(Fld1[,i],Fld2[,i])
  }
  
  #Correlations for Simulations
  M_sim <- list()
  for(j in 1:sims){
    temp <- rep(NA, ncol(Fld1))
    for(i in 1:ncol(Fld1)){ 
      temp[i] <- cor(Fld1_Sims[[j]][,i] , 
                     Fld2_Sims[[j]][,i])
    }
    M_sim[[j]] <- temp}
  
  
  par(mfrow = c(1,1), mar = c(6,6,6,3))
  plot(M,M_sim[[j]], main = "Simulation vs Data Correlation",
       xlab = "Reanalysis Data Correlation", ylab = "Simulation Correlation",
       xlim = c(-0.4,0.4),ylim=c(-0.4,0.4),type='n',
       cex.lab = 1.7,
       cex.main = 2)
  abline(coef = c(0,1), lwd = 2)
  abline(v = 0, lty = 2)
  abline(h = 0, lty = 2)
  
  quant <- c(0.05, 0.25,0.75, 0.9) #Region of Interest
  
  #Computing the Quantiles
  for(i in 1:ncol(Fld1)){
    temp <- list()
    for(j in 1:nsim){
      temp[j] <- M_sim[[j]][i]
    }
    temp <- unlist(temp)
    med <- median(temp)
    qt <- quantile(temp, prob = quant)
    segments(M[i], qt[1], M[i], qt[4], col = "red", lwd = 2)
    #segments(M[i], qt[2], M[i], qt[3], col = "blue")
    points(M[i],med, col ="blue", pch = 19) 
  }
  legend('topleft', legend = c("Median Simulation Value","5 - 95 Percentile Range"),
         col = c('blue','red'),
         lty = c(0,1), pch = c(19,NA), cex = 1.5, lwd = c(NA,3))
  
  
  
  #______________________________________________________________________________#
  #True Data Grid Correlation
  world <- map_data("world")
  us <- map_data("state")
  
  Grid$Corr <- M
 
  p <- ggplot() +
    geom_map(data=world, map=world,
             aes(x=long, y=lat, map_id=region),
             fill="#D3D3D3", color="#000000", size=0.15) +
    geom_map(data=us, map=us,
             aes(x=long, y=lat, map_id=region),
             fill="#D3D3D3", color="#000000", size=0.15) +
    geom_tile(data = Grid, aes(x= lon-360, y = lat, fill = Corr)) +
    scale_fill_gradient2(midpoint=0, low="blue", mid="white",high="red", 
                          limits = c(-0.4, 0.4)) +
    labs(title = "Reanalysis Data Correlation") +
    scale_x_continuous(name = "lon", limits = c(-107, -92)) +
    scale_y_continuous(name = "lat", limits = c(25.5, 36.5)) +
    theme_bw() +
    theme(legend.text=element_text(size=15),
          legend.title=element_text(size=15),
          axis.text=element_text(size=0),
          axis.title=element_text(size=0),
          plot.title = element_text(size=20),
          legend.key.height  = unit(1.25, "cm"),
          axis.ticks = element_blank())
    
  print(p)
  
  
  #______________________________________________________________________________#
  #Median Simulations Grid Correlation
  
  #Computing Median Simulation Values
  med_corr <- as.numeric()
  for(i in 1:ncol(Fld1)){
    temp <- list()
    for(j in 1:nsim){
      temp[j] <- M_sim[[j]][i]
    }
    temp <- unlist(temp)
    med_corr[i] <- median(temp)
  }
  
  
  Grid$Corr <- med_corr

  p <- ggplot() +
    geom_map(data=world, map=world,
             aes(x=long, y=lat, map_id=region),
             fill="#D3D3D3", color="#000000", size=0.15) +
    geom_map(data=us, map=us,
             aes(x=long, y=lat, map_id=region),
             fill="#D3D3D3", color="#000000", size=0.15) +
    geom_tile(data = Grid, aes(x= lon-360, y = lat, fill = Corr)) +
    scale_fill_gradient2(midpoint=0, low="blue", mid="white",high="red", 
                         limits = c(-0.4, 0.4)) +
    labs(title = "Median Simulation Correlation") +
    scale_x_continuous(name = "lon", limits = c(-107, -92)) +
    scale_y_continuous(name = "lat", limits = c(25.5, 36.5)) +
    theme_bw() +
    theme(legend.text=element_text(size=15),
          legend.title=element_text(size=15),
          axis.text=element_text(size=0),
          axis.title=element_text(size=0),
          plot.title = element_text(size=20),
          legend.key.height  = unit(1.25, "cm"),
          axis.ticks = element_blank())
  print(p)
  
  
  
  #______________________________________________________________________________#
  #Difference between Data and Median Values
  
  Grid$Corr <- med_corr - M
  p <- ggplot() +
    geom_map(data=world, map=world,
             aes(x=long, y=lat, map_id=region),
             fill="#D3D3D3", color="#000000", size=0.15) +
    geom_map(data=us, map=us,
             aes(x=long, y=lat, map_id=region),
             fill="#D3D3D3", color="#000000", size=0.15) +
    geom_tile(data = Grid, aes(x= lon-360, y = lat, fill = Corr)) +
    scale_fill_gradient2(midpoint=0, low="blue", mid="white",high="red", 
                         limits = c(-0.4, 0.4)) +
    labs(title = "Difference between Renalysis Data and Simulations") +
    scale_x_continuous(name = "lon", limits = c(-107, -92)) +
    scale_y_continuous(name = "lat", limits = c(25.5, 36.5)) +
    theme_bw() +
    theme(legend.text=element_text(size=15),
          legend.title=element_text(size=15),
          axis.text=element_text(size=0),
          axis.title=element_text(size=0),
          plot.title = element_text(size=20),
          legend.key.height  = unit(1.25, "cm"),
          axis.ticks = element_blank())
  print(p)
  
  
  
  
  #______________________________________________________________________________#
  #Simulation Visualization
  n_site <- sample(ncol(Fld1), 1) #Site
  sim <- sample(6,1)
  sel_sim <- sample(sims, 1) #Simulation Run
  
  par(mfrow=c(3,2), mar = c(4,4,2,2))
  for(i in 1:6){
    #Sample Starting Point
    str <- sample(nrow(Fld1)-100,1)
    
    if(i == sim){
      temp <- scale(as.data.frame(Fld1_Sims[[sel_sim]]))
      plot(temp[str:(str+90),n_site],type='l',
           ylab ="Scaled Units", xlab = "Simulated", lwd = 2, ylim = c(-2.5,2.5))
      temp <- scale(as.data.frame(Fld2_Sims[[sel_sim]]))
      lines(1:91,temp[str:(str+90),n_site], lwd = 2, col="red")
    } else{
      temp <- scale(as.data.frame(Fld1[,n_site]))
      plot(temp[str:(str+90),], type='l',
           xlab =" Observed", ylab = " Scaled Units", lwd = 2, ylim = c(-2.5,2.5))
      temp <- scale(as.data.frame(Fld2[,n_site]))
      lines(1:91, temp[str:(str+90),], lwd = 2, col = "red")
    }
  legend("topright", legend = c("Wind","Solar"), lwd = 2, col = c("black","red"))
  }
  
  
  #______________________________________________________________________________#
  #Simulation Visualization  - Scatter Plots
  n_site <- sample(ncol(Fld1), 1) #Site
  sim <- sample(6,1)
  sel_sim <- sample(sims, 1) #Simulation Run
  
  par(mfrow=c(3,2), mar = c(4,4,2,2))
  for(i in 1:6){
    #Sample Starting Point
    str <- sample(nrow(Fld1)-200,1)
    
    if(i == sim){
      temp1 <- as.data.frame(Fld1_Sims[[sel_sim]])
      temp2 <- as.data.frame(Fld2_Sims[[sel_sim]])
      
      plot(temp1[str:(str+190),n_site],temp2[str:(str+190),n_site], pch = 19,
           ylab ="Solar - Simulated", xlab = "Wind - Simulated")
    } else{
      
      temp1 <- as.data.frame(Fld1[,n_site])
      temp2 <- as.data.frame(Fld2[,n_site])
      
      plot(temp1[str:(str+190),], temp2[str:(str+190),], pch = 19,
           ylab ="Solar - Observed", xlab = "Wind - Observed")

    }
  }
  
  
  
  

}