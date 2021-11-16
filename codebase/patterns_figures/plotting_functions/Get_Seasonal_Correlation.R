#______________________________________________________________________________#
#Code to plot the changing correlation structure betweem two feilds for seasons.
#Seasons - JJA, SON, DJF, MAM. 


#Input
#1. Grid Locations. 
#2. True Fields.
#3. Simulations. 
#4. Start Date.
#5. Color Hex Code

#Output
#1. Simulated vs Data plots for each season. 
#2. Simulated vs Data, Map of Mean Simulated, Data and Difference for each season.
#

#______________________________________________________________________________#
Get_Seasonal_Correlation <- function(Fld1,Fld2, # Data
                                     Fld1_Sims,Fld2_Sims, #Simulations
                                     Grid, #Location
                                     start_date,
                                     col_hx){
  #Load Dependencies
  library(ggplot2)
  library(gridExtra)
  library(usmap)
  
  #Hyper-Parameters
  sims <- length(Fld1_Sims)
  
  st_date <- as.Date(start_date, format="%m-%d-%Y")
  time_stamps <- seq(st_date, by = "day", length.out =nrow(Fld1))
  month_stamps <- as.numeric(format(time_stamps,"%m"))
  
  
  #______________________________________________________________________________#
  ###DJF
  season <- c(12,1,2)
  
  Fld_s1 <- as.data.frame(Fld1); Fld_s1$Month <- month_stamps
  Fld_s1 <- filter(Fld_s1, Month %in% season); Fld_s1$Month <- NULL
  
  Fld_s2 <- as.data.frame(Fld2); Fld_s2$Month <- month_stamps
  Fld_s2 <- filter(Fld_s2, Month %in% season); Fld_s2$Month <- NULL
  
  #Data Correlation
  M <- rep(NA, ncol(Fld_s1))
  for(i in 1:ncol(Fld_s1)){ 
    M[i] <- cor(Fld_s1[,i],Fld_s2[,i])
  }
  
  #Correlations for Simulations
  M_sim <- list()
  for(j in 1:sims){
    temp <- rep(NA, ncol(Fld1))
    for(i in 1:ncol(Fld1)){ 
      
      
      Fld_s1 <- as.data.frame(Fld1_Sims[[j]][,i]); Fld_s1$Month <- month_stamps
      Fld_s1 <- filter(Fld_s1, Month %in% season); Fld_s1$Month <- NULL
      
      Fld_s2 <- as.data.frame(Fld2_Sims[[j]][,i]); Fld_s2$Month <- month_stamps
      Fld_s2 <- filter(Fld_s2, Month %in% season); Fld_s2$Month <- NULL
      
      temp[i] <- cor(Fld_s1 , 
                     Fld_s2)
    }
    M_sim[[j]] <- temp
  }
  
  M_sim  <-  as.data.frame(matrix(unlist(M_sim), nrow=length(M_sim[[1]])))
  median <- apply(M_sim, 1, median)
  
  
  
  #Stick Plot
  st_plot <- data.frame(Dat = M, Sims = median)
  st_plot$Diff <- st_plot$Dat-st_plot$Sims
  
  p11 <- ggplot(st_plot, aes(Dat,Sims)) +
    geom_point(col = col_hx, size = 0.75) +
    geom_abline(intercept = 0, slope = 1) +
    geom_vline(xintercept = 0, col ='black', linetype = "dashed", size = 0.5) +
    geom_hline(yintercept = 0, col='black', linetype = "dashed", size = 0.5) +
    ylab("Simulation Correlation") +
    xlab("Data Correlation") +
    scale_x_continuous(limits = c(-0.4, 0.4)) + 
    scale_y_continuous(limits = c(-0.4, 0.4)) +
    labs(title = "Simulation vs Data Correlation",
         subtitle = "Dec-Jan-Feb") +
    theme_bw() +
    theme(axis.text=element_text(size=5),
          axis.title=element_text(size=5),
          plot.title = element_text(size=6),
          plot.subtitle = element_text(size=8),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank())
  
  #Spatial Correlation for Median Simulation
  Grid_Plot <- data.frame(lon = Grid$lon-360,lat = Grid$lat, Corr = st_plot$Sims)
  #True Data Grid Correlation
  world <- map_data("world")
  us <- map_data("state")
  
  p12 <- ggplot() + 
    geom_map(data=world, map=world,
             aes(x=long, y=lat, map_id=region),
             fill="#D3D3D3", color="#000000", size=0.15) +
    geom_map(data=us, map=us,
             aes(x=long, y=lat, map_id=region),
             fill="#D3D3D3", color="#000000", size=0.15) +
    scale_x_continuous(name = "lon", limits = c(-107, -92)) +
    scale_y_continuous(name = "lat", limits = c(25.5, 36.5)) +
    geom_tile(data = Grid_Plot, aes(x= lon, y = lat, fill = Corr)) +
    scale_fill_gradient2(midpoint=0, low="blue", mid="white",high="red", 
                          limits = c(-0.5, 0.5)) +
    labs(title = "Median Simulation Correlation") +
    theme_bw() +
    theme(axis.text=element_text(size=0),
          axis.title=element_text(size=0),
          axis.ticks = element_blank(),
          plot.title = element_text(size=7),
          legend.text=element_text(size=10),
          legend.title=element_text(size=10),
          legend.key.height  = unit(2.0, "cm"),
          legend.key.width  = unit(0.5, "cm"))
  
  #Spatial Correlation for True Data
  Grid_Plot <- data.frame(lon = Grid$lon-360,lat = Grid$lat, Corr = st_plot$Dat)
  us <- map_data("state", region = "Texas")
  p13 <- ggplot() + 
    geom_map(data=world, map=world,
             aes(x=long, y=lat, map_id=region),
             fill="#D3D3D3", color="#000000", size=0.15) +
    geom_map(data=us, map=us,
             aes(x=long, y=lat, map_id=region),
             fill="#D3D3D3", color="#000000", size=0.15) +
    scale_x_continuous(name = "lon", limits = c(-107, -92)) +
    scale_y_continuous(name = "lat", limits = c(25.5, 36.5)) +
    geom_tile(data = Grid_Plot, aes(x= lon, y = lat, fill = Corr)) +
    scale_fill_gradient2(midpoint=0, low="blue", mid="white",high="red", 
                         limits = c(-0.5, 0.5)) +
    labs(title = "Data Correlation") +
    theme_bw() +
    theme(axis.text=element_text(size=0),
          axis.title=element_text(size=0),
          axis.ticks = element_blank(),
          plot.title = element_text(size=8),
          legend.text=element_text(size=5),
          legend.title=element_text(size=5),
          legend.key.height  = unit(0.5, "cm"),
          legend.key.width  = unit(0.25, "cm"))
  
  #Spatial Correlation for difference
  Grid_Plot <- data.frame(lon = Grid$lon-360,lat = Grid$lat, Corr = st_plot$Diff)
  us <- map_data("state", region = "Texas")
  p14 <- ggplot() + 
    geom_map(data=world, map=world,
             aes(x=long, y=lat, map_id=region),
             fill="#D3D3D3", color="#000000", size=0.15) +
    geom_map(data=us, map=us,
             aes(x=long, y=lat, map_id=region),
             fill="#D3D3D3", color="#000000", size=0.15) +
    scale_x_continuous(name = "lon", limits = c(-107, -92)) +
    scale_y_continuous(name = "lat", limits = c(25.5, 36.5)) +
    geom_tile(data = Grid_Plot, aes(x= lon, y = lat, fill = Corr)) +
    scale_fill_gradient2(midpoint=0, low="blue", mid="white",high="red", 
                         limits = c(-0.5, 0.5)) +
    labs(title = "Difference") +
    theme_bw() +
    theme(axis.text=element_text(size=0),
          axis.title=element_text(size=0),
          axis.ticks = element_blank(),
          plot.title = element_text(size=8),
          legend.text=element_text(size=5),
          legend.title=element_text(size=5),
          legend.key.height  = unit(0.5, "cm"),
          legend.key.width  = unit(0.25, "cm"))
  
  #______________________________________________________________________________#
  ###MAM
  season <- c(3,4,5)
  
  Fld_s1 <- as.data.frame(Fld1); Fld_s1$Month <- month_stamps
  Fld_s1 <- filter(Fld_s1, Month %in% season); Fld_s1$Month <- NULL
  
  Fld_s2 <- as.data.frame(Fld2); Fld_s2$Month <- month_stamps
  Fld_s2 <- filter(Fld_s2, Month %in% season); Fld_s2$Month <- NULL
  
  #Data Correlation
  M <- rep(NA, ncol(Fld_s1))
  for(i in 1:ncol(Fld_s1)){ 
    M[i] <- cor(Fld_s1[,i],Fld_s2[,i])
  }
  
  #Correlations for Simulations
  M_sim <- list()
  for(j in 1:sims){
    temp <- rep(NA, ncol(Fld1))
    for(i in 1:ncol(Fld1)){ 
      
      
      Fld_s1 <- as.data.frame(Fld1_Sims[[j]][,i]); Fld_s1$Month <- month_stamps
      Fld_s1 <- filter(Fld_s1, Month %in% season); Fld_s1$Month <- NULL
      
      Fld_s2 <- as.data.frame(Fld2_Sims[[j]][,i]); Fld_s2$Month <- month_stamps
      Fld_s2 <- filter(Fld_s2, Month %in% season); Fld_s2$Month <- NULL
      
      temp[i] <- cor(Fld_s1 , 
                     Fld_s2)
    }
    M_sim[[j]] <- temp
  }
  
  M_sim  <-  as.data.frame(matrix(unlist(M_sim), nrow=length(M_sim[[1]])))
  median <- apply(M_sim, 1, median)
  
  
  
  #Stick Plot
  st_plot <- data.frame(Dat = M, Sims = median)
  st_plot$Diff <- st_plot$Dat-st_plot$Sims
  
  p21 <- ggplot(st_plot, aes(Dat,Sims)) +
    geom_point(col = col_hx, size = 0.75) +
    geom_abline(intercept = 0, slope = 1) +
    geom_vline(xintercept = 0, col ='black', linetype = "dashed") +
    geom_hline(yintercept = 0, col='black', linetype = "dashed") +
    ylab("Simulation Correlation") +
    xlab("Data Correlation") +
    scale_x_continuous(limits = c(-0.4, 0.4)) + 
    scale_y_continuous(limits = c(-0.4, 0.4)) +
    labs(title = "Simulation vs Data Correlation",
         subtitle = "Mar-April-May") +
    theme_bw() +
    theme(axis.text=element_text(size=5),
          axis.title=element_text(size=5),
          plot.title = element_text(size=6),
          plot.subtitle = element_text(size=8),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank())
  
  #Spatial Correlation for Median Simulation
  Grid_Plot <- data.frame(lon = Grid$lon-360,lat = Grid$lat, Corr = st_plot$Sims)
  
  p22 <- ggplot() + 
    geom_map(data=world, map=world,
             aes(x=long, y=lat, map_id=region),
             fill="#D3D3D3", color="#000000", size=0.15) +
    geom_map(data=us, map=us,
             aes(x=long, y=lat, map_id=region),
             fill="#D3D3D3", color="#000000", size=0.15) +
    scale_x_continuous(name = "lon", limits = c(-107, -92)) +
    scale_y_continuous(name = "lat", limits = c(25.5, 36.5)) +
    geom_tile(data = Grid_Plot, aes(x= lon, y = lat, fill = Corr)) +
    scale_fill_gradient2(midpoint=0, low="blue", mid="white",high="red", 
                         limits = c(-0.5, 0.5)) +
    labs(title = "Median Simulation Correlation") +
    theme_bw() +
    theme(axis.text=element_text(size=0),
          axis.title=element_text(size=0),
          axis.ticks = element_blank(),
          plot.title = element_text(size=7),
          legend.text=element_text(size=5),
          legend.title=element_text(size=5),
          legend.key.height  = unit(0.5, "cm"),
          legend.key.width  = unit(0.25, "cm"))
  
  #Spatial Correlation for True Data
  Grid_Plot <- data.frame(lon = Grid$lon-360,lat = Grid$lat, Corr = st_plot$Dat)
  us <- map_data("state", region = "Texas")
  p23 <- ggplot() + 
    geom_map(data=world, map=world,
             aes(x=long, y=lat, map_id=region),
             fill="#D3D3D3", color="#000000", size=0.15) +
    geom_map(data=us, map=us,
             aes(x=long, y=lat, map_id=region),
             fill="#D3D3D3", color="#000000", size=0.15) +
    scale_x_continuous(name = "lon", limits = c(-107, -92)) +
    scale_y_continuous(name = "lat", limits = c(25.5, 36.5)) +
    geom_tile(data = Grid_Plot, aes(x= lon, y = lat, fill = Corr)) +
    scale_fill_gradient2(midpoint=0, low="blue", mid="white",high="red", 
                         limits = c(-0.5, 0.5)) +
    labs(title = "Data Correlation") +
    theme_bw() +
    theme(axis.text=element_text(size=0),
          axis.title=element_text(size=0),
          axis.ticks = element_blank(),
          plot.title = element_text(size=8),
          legend.text=element_text(size=5),
          legend.title=element_text(size=5),
          legend.key.height  = unit(0.5, "cm"),
          legend.key.width  = unit(0.25, "cm"))
  #Spatial Correlation for difference
  Grid_Plot <- data.frame(lon = Grid$lon-360,lat = Grid$lat, Corr = st_plot$Diff)
  us <- map_data("state", region = "Texas")
  p24 <- ggplot() + 
    geom_map(data=world, map=world,
             aes(x=long, y=lat, map_id=region),
             fill="#D3D3D3", color="#000000", size=0.15) +
    geom_map(data=us, map=us,
             aes(x=long, y=lat, map_id=region),
             fill="#D3D3D3", color="#000000", size=0.15) +
    scale_x_continuous(name = "lon", limits = c(-107, -92)) +
    scale_y_continuous(name = "lat", limits = c(25.5, 36.5)) +
    geom_tile(data = Grid_Plot, aes(x= lon, y = lat, fill = Corr)) +
    scale_fill_gradient2(midpoint=0, low="blue", mid="white",high="red", 
                         limits = c(-0.5, 0.5)) +
    labs(title = "Difference") +
    theme_bw() +
    theme(axis.text=element_text(size=0),
          axis.title=element_text(size=0),
          axis.ticks = element_blank(),
          plot.title = element_text(size=8),
          legend.text=element_text(size=5),
          legend.title=element_text(size=5),
          legend.key.height  = unit(0.5, "cm"),
          legend.key.width  = unit(0.25, "cm"))
  
  #______________________________________________________________________________#
  ###JJA
  season <- c(6,7,8)
  
  Fld_s1 <- as.data.frame(Fld1); Fld_s1$Month <- month_stamps
  Fld_s1 <- filter(Fld_s1, Month %in% season); Fld_s1$Month <- NULL
  
  Fld_s2 <- as.data.frame(Fld2); Fld_s2$Month <- month_stamps
  Fld_s2 <- filter(Fld_s2, Month %in% season); Fld_s2$Month <- NULL
  
  #Data Correlation
  M <- rep(NA, ncol(Fld_s1))
  for(i in 1:ncol(Fld_s1)){ 
    M[i] <- cor(Fld_s1[,i],Fld_s2[,i])
  }
  
  #Correlations for Simulations
  M_sim <- list()
  for(j in 1:sims){
    temp <- rep(NA, ncol(Fld1))
    for(i in 1:ncol(Fld1)){ 
      
      
      Fld_s1 <- as.data.frame(Fld1_Sims[[j]][,i]); Fld_s1$Month <- month_stamps
      Fld_s1 <- filter(Fld_s1, Month %in% season); Fld_s1$Month <- NULL
      
      Fld_s2 <- as.data.frame(Fld2_Sims[[j]][,i]); Fld_s2$Month <- month_stamps
      Fld_s2 <- filter(Fld_s2, Month %in% season); Fld_s2$Month <- NULL
      
      temp[i] <- cor(Fld_s1 , 
                     Fld_s2)
    }
    M_sim[[j]] <- temp
  }
  
  M_sim  <-  as.data.frame(matrix(unlist(M_sim), nrow=length(M_sim[[1]])))
  median <- apply(M_sim, 1, median)
  
  #Stick Plot
  st_plot <- data.frame(Dat = M, Sims = median)
  st_plot$Diff <- st_plot$Dat-st_plot$Sims
  
  p31 <- ggplot(st_plot, aes(Dat,Sims)) +
    geom_point(col = col_hx, size = 0.75) +
    geom_abline(intercept = 0, slope = 1) +
    geom_vline(xintercept = 0, col ='black', linetype = "dashed") +
    geom_hline(yintercept = 0, col='black', linetype = "dashed") +
    ylab("Simulation Correlation") +
    xlab("Data Correlation") +
    scale_x_continuous(limits = c(-0.4, 0.4)) + 
    scale_y_continuous(limits = c(-0.4, 0.4)) +
    labs(title = "Simulation vs Data Correlation",
         subtitle = "Jun-July-Aug") +
    theme_bw() +
    theme(axis.text=element_text(size=5),
          axis.title=element_text(size=5),
          plot.title = element_text(size=6),
          plot.subtitle = element_text(size=8),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank())
  
  #Spatial Correlation for Median Simulation
  Grid_Plot <- data.frame(lon = Grid$lon-360,lat = Grid$lat, Corr = st_plot$Sims)
  
  p32 <- ggplot() + 
    geom_map(data=world, map=world,
             aes(x=long, y=lat, map_id=region),
             fill="#D3D3D3", color="#000000", size=0.15) +
    geom_map(data=us, map=us,
             aes(x=long, y=lat, map_id=region),
             fill="#D3D3D3", color="#000000", size=0.15) +
    scale_x_continuous(name = "lon", limits = c(-107, -92)) +
    scale_y_continuous(name = "lat", limits = c(25.5, 36.5)) +
    geom_tile(data = Grid_Plot, aes(x= lon, y = lat, fill = Corr)) +
    scale_fill_gradient2(midpoint=0, low="blue", mid="white",high="red", 
                         limits = c(-0.5, 0.5)) +
    labs(title = "Median Simulation Correlation") +
    theme_bw() +
    theme(axis.text=element_text(size=0),
          axis.title=element_text(size=0),
          axis.ticks = element_blank(),
          plot.title = element_text(size=7),
          legend.text=element_text(size=5),
          legend.title=element_text(size=5),
          legend.key.height  = unit(0.5, "cm"),
          legend.key.width  = unit(0.25, "cm"))
  
  #Spatial Correlation for True Data
  Grid_Plot <- data.frame(lon = Grid$lon-360,lat = Grid$lat, Corr = st_plot$Dat)
  us <- map_data("state", region = "Texas")
  p33 <- ggplot() + 
    geom_map(data=world, map=world,
             aes(x=long, y=lat, map_id=region),
             fill="#D3D3D3", color="#000000", size=0.15) +
    geom_map(data=us, map=us,
             aes(x=long, y=lat, map_id=region),
             fill="#D3D3D3", color="#000000", size=0.15) +
    scale_x_continuous(name = "lon", limits = c(-107, -92)) +
    scale_y_continuous(name = "lat", limits = c(25.5, 36.5)) +
    geom_tile(data = Grid_Plot, aes(x= lon, y = lat, fill = Corr)) +
    scale_fill_gradient2(midpoint=0, low="blue", mid="white",high="red", 
                         limits = c(-0.5, 0.5)) +
    labs(title = "Data Correlation") +
    theme_bw() +
    theme(axis.text=element_text(size=0),
          axis.title=element_text(size=0),
          axis.ticks = element_blank(),
          plot.title = element_text(size=8),
          legend.text=element_text(size=5),
          legend.title=element_text(size=5),
          legend.key.height  = unit(0.5, "cm"),
          legend.key.width  = unit(0.25, "cm"))
  
  #Spatial Correlation for difference
  Grid_Plot <- data.frame(lon = Grid$lon-360,lat = Grid$lat, Corr = st_plot$Diff)
  us <- map_data("state", region = "Texas")
  p34 <- ggplot() + 
    geom_map(data=world, map=world,
             aes(x=long, y=lat, map_id=region),
             fill="#D3D3D3", color="#000000", size=0.15) +
    geom_map(data=us, map=us,
             aes(x=long, y=lat, map_id=region),
             fill="#D3D3D3", color="#000000", size=0.15) +
    scale_x_continuous(name = "lon", limits = c(-107, -92)) +
    scale_y_continuous(name = "lat", limits = c(25.5, 36.5)) +
    geom_tile(data = Grid_Plot, aes(x= lon, y = lat, fill = Corr)) +
    scale_fill_gradient2(midpoint=0, low="blue", mid="white",high="red", 
                         limits = c(-0.5, 0.5)) +
    labs(title = "Difference") +
    theme_bw() +
    theme(axis.text=element_text(size=0),
          axis.title=element_text(size=0),
          axis.ticks = element_blank(),
          plot.title = element_text(size=8),
          legend.text=element_text(size=5),
          legend.title=element_text(size=5),
          legend.key.height  = unit(0.5, "cm"),
          legend.key.width  = unit(0.25, "cm"))
  
  
  #______________________________________________________________________________#
  ###SON
  season <- c(9,10,11)
  
  Fld_s1 <- as.data.frame(Fld1); Fld_s1$Month <- month_stamps
  Fld_s1 <- filter(Fld_s1, Month %in% season); Fld_s1$Month <- NULL
  
  Fld_s2 <- as.data.frame(Fld2); Fld_s2$Month <- month_stamps
  Fld_s2 <- filter(Fld_s2, Month %in% season); Fld_s2$Month <- NULL
  
  #Data Correlation
  M <- rep(NA, ncol(Fld_s1))
  for(i in 1:ncol(Fld_s1)){ 
    M[i] <- cor(Fld_s1[,i],Fld_s2[,i])
  }
  
  #Correlations for Simulations
  M_sim <- list()
  for(j in 1:sims){
    temp <- rep(NA, ncol(Fld1))
    for(i in 1:ncol(Fld1)){ 
      
      
      Fld_s1 <- as.data.frame(Fld1_Sims[[j]][,i]); Fld_s1$Month <- month_stamps
      Fld_s1 <- filter(Fld_s1, Month %in% season); Fld_s1$Month <- NULL
      
      Fld_s2 <- as.data.frame(Fld2_Sims[[j]][,i]); Fld_s2$Month <- month_stamps
      Fld_s2 <- filter(Fld_s2, Month %in% season); Fld_s2$Month <- NULL
      
      temp[i] <- cor(Fld_s1 , 
                     Fld_s2)
    }
    M_sim[[j]] <- temp
  }
  
  M_sim  <-  as.data.frame(matrix(unlist(M_sim), nrow=length(M_sim[[1]])))
  median <- apply(M_sim, 1, median)
  
  #Stick Plot
  st_plot <- data.frame(Dat = M, Sims = median)
  st_plot$Diff <- st_plot$Dat-st_plot$Sims
  
  p41 <- ggplot(st_plot, aes(Dat,Sims)) +
    geom_point(col = col_hx, size = 0.75) +
    geom_abline(intercept = 0, slope = 1) +
    geom_vline(xintercept = 0, col ='black', linetype = "dashed") +
    geom_hline(yintercept = 0, col='black', linetype = "dashed") +
    ylab("Simulation Correlation") +
    xlab("Data Correlation") +
    scale_x_continuous(limits = c(-0.4, 0.4)) + 
    scale_y_continuous(limits = c(-0.4, 0.4)) +
    labs(title = "Simulation vs Data Correlation",
         subtitle = "Sept-Oct-Nov") +
    theme_bw() +
    theme(axis.text=element_text(size=5),
          axis.title=element_text(size=5),
          plot.title = element_text(size=6),
          plot.subtitle = element_text(size=8),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank())
  
  #Spatial Correlation for Median Simulation
  Grid_Plot <- data.frame(lon = Grid$lon-360,lat = Grid$lat, Corr = st_plot$Sims)
  
  p42 <- ggplot() + 
    geom_map(data=world, map=world,
             aes(x=long, y=lat, map_id=region),
             fill="#D3D3D3", color="#000000", size=0.15) +
    geom_map(data=us, map=us,
             aes(x=long, y=lat, map_id=region),
             fill="#D3D3D3", color="#000000", size=0.15) +
    scale_x_continuous(name = "lon", limits = c(-107, -92)) +
    scale_y_continuous(name = "lat", limits = c(25.5, 36.5)) +
    geom_tile(data = Grid_Plot, aes(x= lon, y = lat, fill = Corr)) +
    scale_fill_gradient2(midpoint=0, low="blue", mid="white",high="red", 
                         limits = c(-0.5, 0.5)) +
    labs(title = "Median Simulation Correlation") +
    theme_bw() +
    theme(axis.text=element_text(size=0),
          axis.title=element_text(size=0),
          axis.ticks = element_blank(),
          plot.title = element_text(size=7),
          legend.text=element_text(size=5),
          legend.title=element_text(size=5),
          legend.key.height  = unit(0.5, "cm"),
          legend.key.width  = unit(0.25, "cm"))
  
  #Spatial Correlation for True Data
  Grid_Plot <- data.frame(lon = Grid$lon-360,lat = Grid$lat, Corr = st_plot$Dat)
  us <- map_data("state", region = "Texas")
  p43 <- ggplot() + 
    geom_map(data=world, map=world,
             aes(x=long, y=lat, map_id=region),
             fill="#D3D3D3", color="#000000", size=0.15) +
    geom_map(data=us, map=us,
             aes(x=long, y=lat, map_id=region),
             fill="#D3D3D3", color="#000000", size=0.15) +
    scale_x_continuous(name = "lon", limits = c(-107, -92)) +
    scale_y_continuous(name = "lat", limits = c(25.5, 36.5)) +
    geom_tile(data = Grid_Plot, aes(x= lon, y = lat, fill = Corr)) +
    scale_fill_gradient2(midpoint=0, low="blue", mid="white",high="red", 
                         limits = c(-0.5, 0.5)) +
    labs(title = "Data Correlation") +
    theme_bw() +
    theme(axis.text=element_text(size=0),
          axis.title=element_text(size=0),
          axis.ticks = element_blank(),
          plot.title = element_text(size=8),
          legend.text=element_text(size=5),
          legend.title=element_text(size=5),
          legend.key.height  = unit(0.5, "cm"),
          legend.key.width  = unit(0.25, "cm"))
  
  #Spatial Correlation for difference
  Grid_Plot <- data.frame(lon = Grid$lon-360,lat = Grid$lat, Corr = st_plot$Diff)
  us <- map_data("state", region = "Texas")
  p44 <- ggplot() + 
    geom_map(data=world, map=world,
             aes(x=long, y=lat, map_id=region),
             fill="#D3D3D3", color="#000000", size=0.15) +
    geom_map(data=us, map=us,
             aes(x=long, y=lat, map_id=region),
             fill="#D3D3D3", color="#000000", size=0.15) +
    scale_x_continuous(name = "lon", limits = c(-107, -92)) +
    scale_y_continuous(name = "lat", limits = c(25.5, 36.5)) +
    geom_tile(data = Grid_Plot, aes(x= lon, y = lat, fill = Corr)) +
    scale_fill_gradient2(midpoint=0, low="blue", mid="white",high="red", 
                         limits = c(-0.5, 0.5)) +
    labs(title = "Difference") +
    theme_bw() +
    theme(axis.text=element_text(size=0),
          axis.title=element_text(size=0),
          axis.ticks = element_blank(),
          plot.title = element_text(size=8),
          legend.text=element_text(size=5),
          legend.title=element_text(size=5),
          legend.key.height  = unit(0.5, "cm"),
          legend.key.width  = unit(0.15, "cm"))

  
  #______________________________________________________________________________#
  

  
  p1 <- plot_grid(p11, 
                  p13 + theme(legend.position="none"), 
                  p12 + theme(legend.position="none"), 
                  p14 + theme(legend.position="none"),
                              nrow = 2)
  
  p2 <- plot_grid(p21, 
                  p23 + theme(legend.position="none"), 
                  p22 + theme(legend.position="none"), 
                  p24 + theme(legend.position="none"),
                  nrow = 2)
  
  p3 <- plot_grid(p31, 
                  p33 + theme(legend.position="none"), 
                  p32 + theme(legend.position="none"), 
                  p34 + theme(legend.position="none"),
                  nrow = 2)
  
  p4 <- plot_grid(p41, 
                  p43 + theme(legend.position="none"), 
                  p42 + theme(legend.position="none"), 
                  p44 + theme(legend.position="none"),
                  nrow = 2)
  
  #Combine all of them
  p <- plot_grid(p1,p2,p3,p4,
                 ncol = 2,
                 labels = c('A', 'B', 'C', 'D'), 
                 label_size = 15)
  
  legend <- get_legend(
    p12 + theme(legend.box.margin = margin(0, 0, 0, 12))
  )

  
  p_seasonal <- plot_grid(p, legend, rel_widths = c(3, .4))
  
  print(p_seasonal)
  


  
}
