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
  #Compute the True Data Correlation
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
  
  
  #Set-up the storage
  Plt_dataset <- data.frame(Reanalysis = M,
                            Median_Corr = rep(NA, length(M)),
                            upper = rep(NA, length(M)),
                            lower = rep(NA, length(M)),
                            Legend1 = "Median Simulation Value",
                            Legend2 = "5-95 Percentile Range")
  
  for(i in 1:ncol(Fld1)){
    temp <- list()
    for(j in 1:sims){
      temp[j] <- M_sim[[j]][i]
    }
    temp <- unlist(temp)
    Plt_dataset$Median_Corr[i] <- median(temp)
    Plt_dataset$lower[i] <- quantile(temp, prob = 0.05)
    Plt_dataset$upper[i] <- quantile(temp, prob = 0.95)
  }
  
  #Plotting the results
  p1 <- ggplot(Plt_dataset) +
    geom_abline(intercept = 0, slope = 1, size = 1) +
    geom_vline(xintercept = 0, linetype =  "dashed") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_linerange(mapping = aes(x=Reanalysis, ymin = lower, 
                                 ymax = upper), size = 0.5, color = 'red') +
    geom_point(mapping = aes(x = Reanalysis, y = Median_Corr, 
                             size = Legend1), color ='blue', size = 0.75) +
    ggtitle("Simulation vs Data Correlation") + 
    scale_x_continuous(name = "Reanalysis Data Correlation",
                       limits = c(-0.4, 0.4)) +
    scale_y_continuous(name = "Simulation Correlation",
                       limits = c(-0.4, 0.4)) +
    theme_bw() +
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=10),
          plot.title = element_text(size=12),
          legend.position = c(0.24,0.81),
          legend.title=element_text(size=0),
          legend.text =element_text(size=5),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank())


  
  #______________________________________________________________________________#
  #True Data Grid Correlation
  world <- map_data("world")
  us <- map_data("state")
  
  Grid$Corr <- Plt_dataset$Reanalysis
 
  p2 <- ggplot() +
    geom_map(data=world, map=world,
             aes(x=long, y=lat, map_id=region),
             fill="#D3D3D3", color="#000000", size=0.15) +
    geom_map(data=us, map=us,
             aes(x=long, y=lat, map_id=region),
             fill="#D3D3D3", color="#000000", size=0.15) +
    geom_tile(data = Grid, aes(x= lon-360, y = lat, fill = Corr)) +
    scale_fill_gradient2(midpoint=0, low="blue", mid="white",high="red", 
                          limits = c(-0.4, 0.4)) +
    labs(title = "  Reanalysis Data Correlation") +
    scale_x_continuous(name = "lon", limits = c(-107, -92)) +
    scale_y_continuous(name = "lat", limits = c(25.5, 36.5)) +
    theme_bw() +
    theme(legend.text=element_text(size=7),
          legend.title=element_text(size=7),
          axis.text=element_text(size=0),
          axis.title=element_text(size=0),
          plot.title = element_text(size=12),
          legend.key.height  = unit(0.75, "cm"),
          axis.ticks = element_blank())
    
  
  
  
  #______________________________________________________________________________#
  #Median Simulations Grid Correlation
  
  Grid$Corr <- Plt_dataset$Median_Corr

  p3 <- ggplot() +
    geom_map(data=world, map=world,
             aes(x=long, y=lat, map_id=region),
             fill="#D3D3D3", color="#000000", size=0.15) +
    geom_map(data=us, map=us,
             aes(x=long, y=lat, map_id=region),
             fill="#D3D3D3", color="#000000", size=0.15) +
    geom_tile(data = Grid, aes(x= lon-360, y = lat, fill = Corr)) +
    scale_fill_gradient2(midpoint=0, low="blue", mid="white",high="red", 
                         limits = c(-0.4, 0.4)) +
    labs(title = "  Median Simulation Correlation") +
    scale_x_continuous(name = "lon", limits = c(-107, -92)) +
    scale_y_continuous(name = "lat", limits = c(25.5, 36.5)) +
    theme_bw() +
    theme(legend.text=element_text(size=7),
          legend.title=element_text(size=7),
          axis.text=element_text(size=0),
          axis.title=element_text(size=0),
          plot.title = element_text(size=12),
          legend.key.height  = unit(0.75, "cm"),
          axis.ticks = element_blank())
  
  
  
  #______________________________________________________________________________#
  #Difference between Data and Median Values
  
  Grid$Corr <- Plt_dataset$Reanalysis-Plt_dataset$Median_Corr
  
  p4 <- ggplot() +
    geom_map(data=world, map=world,
             aes(x=long, y=lat, map_id=region),
             fill="#D3D3D3", color="#000000", size=0.15) +
    geom_map(data=us, map=us,
             aes(x=long, y=lat, map_id=region),
             fill="#D3D3D3", color="#000000", size=0.15) +
    geom_tile(data = Grid, aes(x= lon-360, y = lat, fill = Corr)) +
    scale_fill_gradient2(midpoint=0, low="blue", mid="white",high="red", 
                         limits = c(-0.4, 0.4)) +
    labs(title = "  Difference between (B) and (C)") +
    scale_x_continuous(name = "lon", limits = c(-107, -92)) +
    scale_y_continuous(name = "lat", limits = c(25.5, 36.5)) +
    theme_bw() +
    theme(legend.text=element_text(size=7),
          legend.title=element_text(size=7),
          axis.text=element_text(size=0),
          axis.title=element_text(size=0),
          plot.title = element_text(size=12),
          legend.key.height  = unit(0.75, "cm"),
          axis.ticks = element_blank())
  
  
  #Plotting the combined plots
  plot_grid(p1,p2,p3,p4,
              ncol = 2,
              labels = c('A', 'B', 'C', 'D'),
              label_size = 12)

  

}