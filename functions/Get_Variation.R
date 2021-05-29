#______________________________________________________________________________#
#Code to get seasonal variation at each Site across the seasons.
#Seasons - JJA, SON, DJF, MAM. 
#Code is for fields individually.


#Input
#1. Grid Locations. 
#2. True Fields.
#3. Simulations. 
#4. Start Date.
#5. Field_Name

#Output
#1. Mean power at each Site.
#2. Mean power at each Site by season
#3. Variation at each Site.
#4, Variation at each Site by season.




setwd("~/ERCOT") #Code for personal device


#______________________________________________________________________________#
###Load Packages
library(maps)       
library(dplyr)
library(ggplot2)


#______________________________________________________________________________#
###Reading the Data
grid_locs <- read.csv("data/ERCOT_0_5_deg_lat_lon_index_key.csv", 
                      header = TRUE, sep=",")


ssrd <- read.table("data/ERCOT_Solar_Rad_Daily.txt",  
                   sep =" ", 
                   header = TRUE) #Surface Solar Radiation Downwards
WP <- read.table("data/ERCOT_Wind_Power_Daily.txt",  
                 sep =" ", 
                 header = TRUE) #Surface Solar Radiation Downwards


#______________________________________________________________________________#
Get_Variation <- function(Dat, Grid, start_date, frac, Field_Title, leg_title){
  
  #Load Dependencies
  library(ggplot2)
  library(gridExtra)
  library(usmap)
  library(ggpubr)
  
  #Dates
  st_date <- as.Date(start_date, format="%m-%d-%Y")
  time_stamps <- seq(st_date, by = "day", length.out =nrow(Dat))
  month_stamps <- as.numeric(format(time_stamps,"%m"))
  
  #Transform to needed Number.
  Dat <- Dat*frac
  
  #------------------------------------------------------------------------------#
  ###Code to get the color-scale
  #Annual
  mean_cl <- var_cl <-  list()
  for(i in 1:ncol(Dat)){
    mean_cl[[i]] <- mean(Dat[,i])
    var_cl[[i]] <- (quantile(Dat[,i], probs = c(0.9)) -  quantile(Dat[,i], probs = c(0.1)))/mean_cl[[i]]
  }
  
  #DJF
  mean_djf <- var_djf <- list()
  season <- c(12,1,2)
  Fld_s1 <- as.data.frame(Dat); Fld_s1$Month <- month_stamps
  Fld_s1 <- filter(Fld_s1, Month %in% season); Fld_s1$Month <- NULL
  for(i in 1:ncol(Fld_s1)){
    mean_djf[[i]] <- mean(Fld_s1[,i])
    var_djf[[i]] <- (quantile(Fld_s1[,i], probs = c(0.9)) -  quantile(Fld_s1[,i], probs = c(0.1)))/mean_djf[[i]]
  }
  
  #MAM
  mean_mam <- var_mam <- list()
  season <- c(3,4,5)
  Fld_s1 <- as.data.frame(Dat); Fld_s1$Month <- month_stamps
  Fld_s1 <- filter(Fld_s1, Month %in% season); Fld_s1$Month <- NULL
  for(i in 1:ncol(Fld_s1)){
    mean_mam[[i]] <- mean(Fld_s1[,i])
    var_mam[[i]] <- (quantile(Fld_s1[,i], probs = c(0.9)) -  quantile(Fld_s1[,i], probs = c(0.1)))/mean_mam[[i]]
  }
  
  #JJA
  mean_jja <- var_jja <- list()
  season <- c(6,7,8)
  Fld_s1 <- as.data.frame(Dat); Fld_s1$Month <- month_stamps
  Fld_s1 <- filter(Fld_s1, Month %in% season); Fld_s1$Month <- NULL
  for(i in 1:ncol(Fld_s1)){
    mean_jja[[i]] <- mean(Fld_s1[,i])
    var_jja[[i]] <- (quantile(Fld_s1[,i], probs = c(0.9)) -  quantile(Fld_s1[,i], probs = c(0.1)))/mean_jja[[i]]
  }
  
  #SON
  mean_son <- var_son <- list()
  season <- c(9,10,11)
  Fld_s1 <- as.data.frame(Dat); Fld_s1$Month <- month_stamps
  Fld_s1 <- filter(Fld_s1, Month %in% season); Fld_s1$Month <- NULL
  for(i in 1:ncol(Fld_s1)){
    mean_son[[i]] <- mean(Fld_s1[,i])
    var_son[[i]] <- (quantile(Fld_s1[,i], probs = c(0.9)) -  quantile(Fld_s1[,i], probs = c(0.1)))/mean_son[[i]]
  }
  
  #Mean
  mean_site <- c(unlist(mean_cl), unlist(mean_djf), unlist(mean_mam),
                     unlist(mean_jja), unlist(mean_son))
  min_mean <- min(mean_site)
  max_mean <- max(mean_site)
  mid_mean <- (max_mean + min_mean)/2
  
  #Variation
  var_site <- c(unlist(var_cl), unlist(var_djf), unlist(var_mam),
                 unlist(var_jja), unlist(var_son))
  min_var <- min(var_site)
  max_var <- max(var_site)
  mid_var <- (max_var + min_var)/2
  
  

  #------------------------------------------------------------------------------#
  ###ANNUAL
  mean_site <- low_site <- high_site <- list()
  for(i in 1:ncol(Dat)){
    mean_site[[i]] <- mean(Dat[,i])
    low_site[[i]] <- quantile(Dat[,i], probs = c(0.1))
    high_site[[i]] <- quantile(Dat[,i], probs = c(0.9))
  }
  mean_site <- unlist(mean_site)
  low_site <- unlist(low_site)
  high_site <- unlist(high_site)
  
  #Spatial Correlation for Median Simulation
  Grid_Plot <- data.frame(lon = Grid$lon-360,lat = Grid$lat,
                          Mean = mean_site,
                          Variation = (high_site - low_site)/mean_site)
  
  world <- map_data("world")
  us <- map_data("state")

  p1 <- ggplot() +
    geom_map(data=world, map=world,
             aes(x=long, y=lat, map_id=region),
             fill="#D3D3D3", color="#000000", size=0.15) +
    geom_map(data=us, map=us,
             aes(x=long, y=lat, map_id=region),
             fill="#D3D3D3", color="#000000", size=0.15) +
    geom_tile(data = Grid_Plot, aes(x= lon, y = lat, fill = Mean)) +
    scale_fill_gradient2(midpoint=mid_mean, low="blue", mid="white",high="red", 
                         limits = c(min_mean, max_mean)) +
    labs(title = paste0("Mean Production - ", Field_Title)) +
    scale_x_continuous(name = "lon", limits = c(-107, -92)) +
    scale_y_continuous(name = "lat", limits = c(25.5, 36.5)) +
    labs(color=leg_title)  +
    theme_bw() +
    theme(legend.text=element_text(size=20),
          legend.title=element_text(size=20),
          axis.text=element_text(size=0),
          axis.title=element_text(size=0),
          axis.ticks = element_blank(),
          plot.title = element_text(size=20),
          legend.key.height  = unit(2.5, "cm"))
  
  p2 <- ggplot() + 
    geom_map(data=world, map=world,
             aes(x=long, y=lat, map_id=region),
             fill="#D3D3D3", color="#000000", size=0.15) +
    geom_map(data=us, map=us,
             aes(x=long, y=lat, map_id=region),
             fill="#D3D3D3", color="#000000", size=0.15) +
    geom_tile(data = Grid_Plot, aes(x= lon, y = lat, fill = Variation)) +
    scale_fill_gradient2(midpoint=mid_var, low="blue", mid="white",high="red", 
                         limits = c(min_var, max_var)) +
    labs(title = paste0("Variation in Production - ", Field_Title)) +
    scale_x_continuous(name = "lon", limits = c(-107, -92)) +
    scale_y_continuous(name = "lat", limits = c(25.5, 36.6)) +
    labs(color='Variation')  +
    theme_bw() +
    theme(legend.text=element_text(size=20),
          legend.title=element_text(size=20),
          axis.text=element_text(size=0),
          axis.title=element_text(size=0),
          axis.ticks = element_blank(),
          plot.title = element_text(size=20),
          legend.key.height  = unit(2.5, "cm"))

  
  grid.arrange(p1, nrow = 1)
  grid.arrange(p2, nrow = 1)
  
  #______________________________________________________________________________#
  #Seasonal
  
  #------------------------------------------------------------------------------#
  #DJF
  season <- c(12,1,2)
  name <- "DJF"
  
  Fld_s1 <- as.data.frame(Dat); Fld_s1$Month <- month_stamps
  Fld_s1 <- filter(Fld_s1, Month %in% season); Fld_s1$Month <- NULL
  
  #Compute Mean and Quantiles
  mean_site <- low_site <- high_site <- list()
  for(i in 1:ncol(Fld_s1)){
    mean_site[[i]] <- mean(Fld_s1[,i])
    low_site[[i]] <- quantile(Fld_s1[,i], probs = c(0.1))
    high_site[[i]] <- quantile(Fld_s1[,i], probs = c(0.9))
  }
  mean_site <- unlist(mean_site)
  low_site <- unlist(low_site)
  high_site <- unlist(high_site)
  
  #Spatial Correlation for Median Simulation
  Grid_Plot <- data.frame(lon = Grid$lon-360,lat = Grid$lat,
                          Mean = mean_site,
                          Variation = (high_site - low_site)/mean_site)
  
  p11 <- ggplot() +
    geom_map(data=world, map=world,
             aes(x=long, y=lat, map_id=region),
             fill="#D3D3D3", color="#000000", size=0.15) +
    geom_map(data=us, map=us,
             aes(x=long, y=lat, map_id=region),
             fill="#D3D3D3", color="#000000", size=0.15) +
    geom_tile(data = Grid_Plot, aes(x= lon, y = lat, fill = Mean)) +
    scale_fill_gradient2(midpoint=mid_mean, low="blue", mid="white",high="red", 
                         limits = c(min_mean, max_mean)) +
    labs(title = paste0(name)) +
    scale_x_continuous(name = "lon", limits = c(-107, -92)) +
    scale_y_continuous(name = "lat", limits = c(25.5, 36.5)) +
    labs(color=leg_title)  +
    theme_bw() +
    theme(axis.text=element_text(size=0),
          axis.title=element_text(size=0),
          axis.ticks = element_blank(),
          plot.title = element_text(size=20),
          legend.key.height  = unit(1.75, "cm"),
          legend.text=element_text(size=20),
          legend.title=element_text(size=20))
  
  p12 <- ggplot() + 
    geom_map(data=world, map=world,
             aes(x=long, y=lat, map_id=region),
             fill="#D3D3D3", color="#000000", size=0.15) +
    geom_map(data=us, map=us,
             aes(x=long, y=lat, map_id=region),
             fill="#D3D3D3", color="#000000", size=0.15) +
    geom_tile(data = Grid_Plot, aes(x= lon, y = lat, fill = Variation)) +
    scale_fill_gradient2(midpoint=mid_var, low="blue", mid="white",high="red", 
                         limits = c(min_var, max_var)) +
    labs(title = paste0(name)) +
    scale_x_continuous(name = "lon", limits = c(-107, -92)) +
    scale_y_continuous(name = "lat", limits = c(25.5, 36.6)) +
    labs(color='Variation')  +
    theme_bw() +
    theme(axis.text=element_text(size=0),
          axis.title=element_text(size=0),
          axis.ticks = element_blank(),
          plot.title = element_text(size=20),
          legend.key.height  = unit(1.75, "cm"),
          legend.text=element_text(size=20),
          legend.title=element_text(size=20))  #,
  #legend.position = "none")
  
  #------------------------------------------------------------------------------#
  #MAM
  season <- c(3,4,5)
  name <- "MAM"
  
  Fld_s1 <- as.data.frame(Dat); Fld_s1$Month <- month_stamps
  Fld_s1 <- filter(Fld_s1, Month %in% season); Fld_s1$Month <- NULL
  
  #Compute Mean and Quantiles
  mean_site <- low_site <- high_site <- list()
  for(i in 1:ncol(Fld_s1)){
    mean_site[[i]] <- mean(Fld_s1[,i])
    low_site[[i]] <- quantile(Fld_s1[,i], probs = c(0.1))
    high_site[[i]] <- quantile(Fld_s1[,i], probs = c(0.9))
  }
  mean_site <- unlist(mean_site)
  low_site <- unlist(low_site)
  high_site <- unlist(high_site)
  
  #Spatial Correlation for Median Simulation
  Grid_Plot <- data.frame(lon = Grid$lon-360,lat = Grid$lat,
                          Mean = mean_site,
                          Variation = (high_site - low_site)/mean_site)
  
  p21 <- ggplot() +
    geom_map(data=world, map=world,
             aes(x=long, y=lat, map_id=region),
             fill="#D3D3D3", color="#000000", size=0.15) +
    geom_map(data=us, map=us,
             aes(x=long, y=lat, map_id=region),
             fill="#D3D3D3", color="#000000", size=0.15) +
    geom_tile(data = Grid_Plot, aes(x= lon, y = lat, fill = Mean)) +
    scale_fill_gradient2(midpoint=mid_mean, low="blue", mid="white",high="red", 
                         limits = c(min_mean, max_mean)) +
    labs(title = paste0(name)) +
    scale_x_continuous(name = "lon", limits = c(-107, -92)) +
    scale_y_continuous(name = "lat", limits = c(25.5, 36.5)) +
    labs(color=leg_title)  +
    theme_bw() +
    theme(axis.text=element_text(size=0),
          axis.title=element_text(size=0),
          axis.ticks = element_blank(),
          plot.title = element_text(size=20),
          legend.key.height  = unit(1.75, "cm"),
          legend.text=element_text(size=20),
          legend.title=element_text(size=20))  #,
  #legend.position = "none")
  
  p22 <- ggplot() + 
    geom_map(data=world, map=world,
             aes(x=long, y=lat, map_id=region),
             fill="#D3D3D3", color="#000000", size=0.15) +
    geom_map(data=us, map=us,
             aes(x=long, y=lat, map_id=region),
             fill="#D3D3D3", color="#000000", size=0.15) +
    geom_tile(data = Grid_Plot, aes(x= lon, y = lat, fill = Variation)) +
    scale_fill_gradient2(midpoint=mid_var, low="blue", mid="white",high="red", 
                         limits = c(min_var, max_var)) +
    labs(title = paste0(name)) +
    scale_x_continuous(name = "lon", limits = c(-107, -92)) +
    scale_y_continuous(name = "lat", limits = c(25.5, 36.6)) +
    labs(color='Variation')  +
    theme_bw() +
    theme(axis.text=element_text(size=0),
          axis.title=element_text(size=0),
          axis.ticks = element_blank(),
          plot.title = element_text(size=20),
          legend.key.height  = unit(1.75, "cm"),
          legend.text=element_text(size=20),
          legend.title=element_text(size=20))  #,
  #legend.position = "none")
  
  #------------------------------------------------------------------------------#
  #JJA
  season <- c(6,7,8)
  name <- "JJA"
  
  Fld_s1 <- as.data.frame(Dat); Fld_s1$Month <- month_stamps
  Fld_s1 <- filter(Fld_s1, Month %in% season); Fld_s1$Month <- NULL
  
  #Compute Mean and Quantiles
  mean_site <- low_site <- high_site <- list()
  for(i in 1:ncol(Fld_s1)){
    mean_site[[i]] <- mean(Fld_s1[,i])
    low_site[[i]] <- quantile(Fld_s1[,i], probs = c(0.1))
    high_site[[i]] <- quantile(Fld_s1[,i], probs = c(0.9))
  }
  mean_site <- unlist(mean_site)
  low_site <- unlist(low_site)
  high_site <- unlist(high_site)
  
  #Spatial Correlation for Median Simulation
  Grid_Plot <- data.frame(lon = Grid$lon-360,lat = Grid$lat,
                          Mean = mean_site,
                          Variation = (high_site - low_site)/mean_site)
  
  p31 <- ggplot() +
    geom_map(data=world, map=world,
             aes(x=long, y=lat, map_id=region),
             fill="#D3D3D3", color="#000000", size=0.15) +
    geom_map(data=us, map=us,
             aes(x=long, y=lat, map_id=region),
             fill="#D3D3D3", color="#000000", size=0.15) +
    geom_tile(data = Grid_Plot, aes(x= lon, y = lat, fill = Mean)) +
    scale_fill_gradient2(midpoint=mid_mean, low="blue", mid="white",high="red", 
                         limits = c(min_mean, max_mean)) +
    labs(title = paste0(name)) +
    scale_x_continuous(name = "lon", limits = c(-107, -92)) +
    scale_y_continuous(name = "lat", limits = c(25.5, 36.5)) +
    labs(color=leg_title)  +
    theme_bw() +
    theme(axis.text=element_text(size=0),
          axis.title=element_text(size=0),
          axis.ticks = element_blank(),
          plot.title = element_text(size=20),
          legend.key.height  = unit(1.75, "cm"),
          legend.text=element_text(size=20),
          legend.title=element_text(size=20))  #,
  #legend.position = "none")
  
  p32 <- ggplot() + 
    geom_map(data=world, map=world,
             aes(x=long, y=lat, map_id=region),
             fill="#D3D3D3", color="#000000", size=0.15) +
    geom_map(data=us, map=us,
             aes(x=long, y=lat, map_id=region),
             fill="#D3D3D3", color="#000000", size=0.15) +
    geom_tile(data = Grid_Plot, aes(x= lon, y = lat, fill = Variation)) +
    scale_fill_gradient2(midpoint=mid_var, low="blue", mid="white",high="red", 
                         limits = c(min_var, max_var)) +
    labs(title = paste0(name)) +
    scale_x_continuous(name = "lon", limits = c(-107, -92)) +
    scale_y_continuous(name = "lat", limits = c(25.5, 36.6)) +
    labs(color='Variation')  +
    theme_bw() +
    theme(axis.text=element_text(size=0),
          axis.title=element_text(size=0),
          axis.ticks = element_blank(),
          plot.title = element_text(size=20),
          legend.key.height  = unit(1.75, "cm"),
          legend.text=element_text(size=20),
          legend.title=element_text(size=20))  #,
  #legend.position = "none")
  
  #------------------------------------------------------------------------------#
  #DJF
  season <- c(9,10,11)
  name <- "SON"
  
  Fld_s1 <- as.data.frame(Dat); Fld_s1$Month <- month_stamps
  Fld_s1 <- filter(Fld_s1, Month %in% season); Fld_s1$Month <- NULL
  
  #Compute Mean and Quantiles
  mean_site <- low_site <- high_site <- list()
  for(i in 1:ncol(Fld_s1)){
    mean_site[[i]] <- mean(Fld_s1[,i])
    low_site[[i]] <- quantile(Fld_s1[,i], probs = c(0.1))
    high_site[[i]] <- quantile(Fld_s1[,i], probs = c(0.9))
  }
  mean_site <- unlist(mean_site)
  low_site <- unlist(low_site)
  high_site <- unlist(high_site)
  
  #Spatial Correlation for Median Simulation
  Grid_Plot <- data.frame(lon = Grid$lon-360,lat = Grid$lat,
                          Mean = mean_site,
                          Variation = (high_site - low_site)/mean_site)
  
  p41 <- ggplot() +
    geom_map(data=world, map=world,
             aes(x=long, y=lat, map_id=region),
             fill="#D3D3D3", color="#000000", size=0.15) +
    geom_map(data=us, map=us,
             aes(x=long, y=lat, map_id=region),
             fill="#D3D3D3", color="#000000", size=0.15) +
    geom_tile(data = Grid_Plot, aes(x= lon, y = lat, fill = Mean)) +
    scale_fill_gradient2(midpoint=mid_mean, low="blue", mid="white",high="red", 
                         limits = c(min_mean, max_mean)) +
    labs(title = paste0(name)) +
    scale_x_continuous(name = "lon", limits = c(-107, -92)) +
    scale_y_continuous(name = "lat", limits = c(25.5, 36.5)) +
    labs(color=leg_title)  +
    theme_bw() +
    theme(axis.text=element_text(size=0),
          axis.title=element_text(size=0),
          axis.ticks = element_blank(),
          plot.title = element_text(size=20),
          legend.key.height  = unit(1.75, "cm"),
          legend.text=element_text(size=20),
          legend.title=element_text(size=20))  #,
  #legend.position = "none")
  
  p42 <- ggplot() + 
    geom_map(data=world, map=world,
             aes(x=long, y=lat, map_id=region),
             fill="#D3D3D3", color="#000000", size=0.15) +
    geom_map(data=us, map=us,
             aes(x=long, y=lat, map_id=region),
             fill="#D3D3D3", color="#000000", size=0.15) +
    geom_tile(data = Grid_Plot, aes(x= lon, y = lat, fill = Variation)) +
    scale_fill_gradient2(midpoint=mid_var, low="blue", mid="white",high="red", 
                         limits = c(min_var, max_var)) +
    labs(title = paste0(name)) +
    scale_x_continuous(name = "lon", limits = c(-107, -92)) +
    scale_y_continuous(name = "lat", limits = c(25.5, 36.6)) +
    labs(color='Variation')  +
    theme_bw() +
    theme(axis.text=element_text(size=0),
          axis.title=element_text(size=0),
          axis.ticks = element_blank(),
          plot.title = element_text(size=20),
          legend.key.height  = unit(1.75, "cm"),
          legend.text=element_text(size=20),
          legend.title=element_text(size=20))  #,
          #legend.position = "none")
  
  
  #______________________________________________________________________________#
  
  #grid.arrange(p11, p21, p31, p41, nrow = 2)
  print(ggarrange(p11, p21, p31, p41, nrow = 2, ncol = 2, common.legend = TRUE, legend = "right"))
  #grid.arrange(p12, p22, p32, p42, nrow = 2)
  print(ggarrange(p12, p22, p32, p42, nrow = 2, ncol = 2, common.legend = TRUE, legend = "right"))
  

}




#______________________________________________________________________________#
pdf("ERCOT_Variation.pdf")
Get_Variation(Dat = WP, 
              Grid = grid_locs, 
              start_date = "01-01-1979",
              frac = 4*7*90*90/(2*10^6),
              Field_Title = "Wind",
              leg_title = "CF")


Get_Variation(Dat = ssrd,
              Grid = grid_locs, 
              start_date = "01-01-1979",
              frac = 1,
              Field_Title = "Solar",
              leg_title = "W/sq-m")
dev.off()
