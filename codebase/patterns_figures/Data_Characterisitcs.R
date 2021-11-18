#______________________________________________________________________________#
###Updated Code for Patterns Figures

###Wind and Solar Characteristics across the Texas Interconnect###


###Script to generate plots for
#1. PC-1 Eigenvectors for Wind and Solar Fields
#2. ACF/PACF Plots for PC-1 for Wind and Solar Fields
#3. Annual and Seasonal Mean and Variation for Wind-Solar 

###Data Inputs
#1. Wind Data
#2. Solar Data
#3. Grid Locations

###Outputs
#1. PCA Eigenvectors and Auto-Correlations
#2. Gridwise Mean and Variation


#______________________________________________________________________________
###Set-up working directory###
setwd("~/GitHub/LDS-Inferences") #Code for personal device

###Load Packages and Dependencies
library(maps)       
library(dplyr)
library(ggplot2)



#______________________________________________________________________________#
###Reading the Data###

#Load the Grid Locations
grid_locs <- read.csv("data/ERCOT_0_5_deg_lat_lon_index_key.csv", 
                      header = TRUE, sep=",")

#Load the Downward Surface Solar Radiation
ssrd <- read.table("data/ERCOT_Solar_CF_Daily.txt",  
                   sep =" ", 
                   header = TRUE) 

#Load the Wind Power Data
WP <- read.table("data/ERCOT_Wind_CF_Daily.txt",  
                 sep =" ", 
                 header = TRUE) 


#______________________________________________________________________________#
###Get Generation Variation###





###Function 1 
#Objective - Compute Mean and Variation annual and seasonally
###Inputs
#   1. Data Field (Dat)
#   2. Grid Locations (Grid)
#   3. Start Date (start_date)
#   4. Fraction - Transformational Scaling Parameter (frac)
#   5. Field Title (Field_Title)
#   6. Legend Title (leg_title)
###Output
#   1. Daily Mean and Variation in Generation
#   2. Seasonal Mean and Variation in Generation


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
    labs(title = paste0("  Mean - ", Field_Title)) +
    scale_x_continuous(name = "lon", limits = c(-107, -92)) +
    scale_y_continuous(name = "lat", limits = c(25.5, 36.5)) +
    labs(color=leg_title)  +
    theme_bw() +
    theme(legend.text=element_text(size=10),
          legend.title=element_text(size=10),
          axis.text=element_text(size=0),
          axis.title=element_text(size=0),
          axis.ticks = element_blank(),
          plot.title = element_text(size=15),
          legend.key.height  = unit(1, "cm"))
  
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
    labs(title = paste0("  Variation - ", Field_Title)) +
    scale_x_continuous(name = "lon", limits = c(-107, -92)) +
    scale_y_continuous(name = "lat", limits = c(25.5, 36.6)) +
    labs(color='Variation')  +
    theme_bw() +
    theme(legend.text=element_text(size=10),
          legend.title=element_text(size=10),
          axis.text=element_text(size=0),
          axis.title=element_text(size=0),
          axis.ticks = element_blank(),
          plot.title = element_text(size=15),
          legend.key.height  = unit(1, "cm"))

  
  #grid.arrange(p1, nrow = 1)
  #grid.arrange(p2, nrow = 1)
  
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
    labs(title = paste0(" ",name)) +
    scale_x_continuous(name = "lon", limits = c(-107, -92)) +
    scale_y_continuous(name = "lat", limits = c(25.5, 36.5)) +
    labs(color=leg_title)  +
    theme_bw() +
    theme(axis.text=element_text(size=0),
          axis.title=element_text(size=0),
          axis.ticks = element_blank(),
          plot.title = element_text(size=10),
          legend.key.height  = unit(1, "cm"),
          legend.text=element_text(size=10),
          legend.title=element_text(size=10))  
  
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
    labs(title = paste0(" ",name)) +
    scale_x_continuous(name = "lon", limits = c(-107, -92)) +
    scale_y_continuous(name = "lat", limits = c(25.5, 36.6)) +
    labs(color='Variation')  +
    theme_bw() +
    theme(axis.text=element_text(size=0),
          axis.title=element_text(size=0),
          axis.ticks = element_blank(),
          plot.title = element_text(size=10),
          legend.key.height  = unit(1, "cm"),
          legend.text=element_text(size=10),
          legend.title=element_text(size=10))  
  
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
    labs(title = paste0(" ",name)) +
    scale_x_continuous(name = "lon", limits = c(-107, -92)) +
    scale_y_continuous(name = "lat", limits = c(25.5, 36.5)) +
    labs(color=leg_title)  +
    theme_bw() +
    theme(axis.text=element_text(size=0),
          axis.title=element_text(size=0),
          axis.ticks = element_blank(),
          plot.title = element_text(size=10),
          legend.key.height  = unit(1, "cm"),
          legend.text=element_text(size=10),
          legend.title=element_text(size=10))  
  
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
    labs(title = paste0(" ",name)) +
    scale_x_continuous(name = "lon", limits = c(-107, -92)) +
    scale_y_continuous(name = "lat", limits = c(25.5, 36.6)) +
    labs(color='Variation')  +
    theme_bw() +
    theme(axis.text=element_text(size=0),
          axis.title=element_text(size=0),
          axis.ticks = element_blank(),
          plot.title = element_text(size=10),
          legend.key.height  = unit(1, "cm"),
          legend.text=element_text(size=10),
          legend.title=element_text(size=10))  
  
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
    labs(title = paste0(" ",name)) +
    scale_x_continuous(name = "lon", limits = c(-107, -92)) +
    scale_y_continuous(name = "lat", limits = c(25.5, 36.5)) +
    labs(color=leg_title)  +
    theme_bw() +
    theme(axis.text=element_text(size=0),
          axis.title=element_text(size=0),
          axis.ticks = element_blank(),
          plot.title = element_text(size=10),
          legend.key.height  = unit(1, "cm"),
          legend.text=element_text(size=10),
          legend.title=element_text(size=10))  
  
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
    labs(title = paste0(" ",name)) +
    scale_x_continuous(name = "lon", limits = c(-107, -92)) +
    scale_y_continuous(name = "lat", limits = c(25.5, 36.6)) +
    labs(color='Variation')  +
    theme_bw() +
    theme(axis.text=element_text(size=0),
          axis.title=element_text(size=0),
          axis.ticks = element_blank(),
          plot.title = element_text(size=10),
          legend.key.height  = unit(1, "cm"),
          legend.text=element_text(size=10),
          legend.title=element_text(size=10))  
  
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
    labs(title = paste0(" ",name)) +
    scale_x_continuous(name = "lon", limits = c(-107, -92)) +
    scale_y_continuous(name = "lat", limits = c(25.5, 36.5)) +
    labs(color=leg_title)  +
    theme_bw() +
    theme(axis.text=element_text(size=0),
          axis.title=element_text(size=0),
          axis.ticks = element_blank(),
          plot.title = element_text(size=10),
          legend.key.height  = unit(1, "cm"),
          legend.text=element_text(size=10),
          legend.title=element_text(size=10))  
  
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
    labs(title = paste0(" ",name)) +
    scale_x_continuous(name = "lon", limits = c(-107, -92)) +
    scale_y_continuous(name = "lat", limits = c(25.5, 36.6)) +
    labs(color='Variation')  +
    theme_bw() +
    theme(axis.text=element_text(size=0),
          axis.title=element_text(size=0),
          axis.ticks = element_blank(),
          plot.title = element_text(size=10),
          legend.key.height  = unit(1, "cm"),
          legend.text=element_text(size=10),
          legend.title=element_text(size=10))  
  
  
  #______________________________________________________________________________#
  
  #print(ggarrange(p11, p21, p31, p41, nrow = 2, ncol = 2, common.legend = TRUE, legend = "right",
  #                labels = c("a", "b", "c", "d")))
  #print(ggarrange(p12, p22, p32, p42, nrow = 2, ncol = 2, common.legend = TRUE, legend = "right",
  #                labels = c("a", "b", "c", "d")))
  p3 <- ggarrange(p11, p21, p31, p41, nrow = 2, ncol = 2, common.legend = TRUE, 
                  legend = "right")
  p4 <- ggarrange(p12, p22, p32, p42, nrow = 2, ncol = 2, common.legend = TRUE, 
                  legend = "right")
  
  out = list(p1=p1,p2=p2,p3=p3,p4=p4)
  
}


#______________________________________________________________________________#
###Get Generation Variation###

###Function to get the 1st PCA and its PACF
###Input
#   1. Field
#   2. Field Name
###Output
#   1. PACF Plot

get_pca_pacf <- function(Fld, Field_Name){
  
  #PCA on the Data
  pcw <- prcomp(Fld, scale = TRUE, center = TRUE)
  var <- (pcw$sdev^2)
  var <- (var)/sum(var)
  
  #Plotting Partial correlation
  PCA_Dataset <- data.frame(Lag = acf(pcw$x[,1], plot = FALSE)$lag, 
                            ACF = acf(pcw$x[,1], plot = FALSE)$acf,
                            PACF = c(0,pacf(pcw$x[,1], plot = FALSE)$acf))
  Signif_Threshold <- 2/sqrt(length(pcw$x[,1]))
  
  #Plotting ACF
  p1 <- ggplot() +
    geom_segment(PCA_Dataset, mapping = aes(x=Lag, xend=Lag, y=0, yend=ACF)) +
    geom_hline(yintercept = Signif_Threshold, color = 'blue', linetype = 'dashed') +
    geom_hline(yintercept = -Signif_Threshold, color = 'blue', linetype = 'dashed') +
    scale_x_continuous(name = "Lag (Days)", limits = c(0, 30)) +
    scale_y_continuous(name = "ACF", limits = c(-0.1, 1)) +
    labs(title = paste0("  ACF Reanalysis Data PC-1 ",Field_Name)) +
    theme_bw() +
    theme(plot.title = element_text(size=10))  
  
  
  #Plotting ACF
  p2 <- ggplot() +
    geom_segment(PCA_Dataset, mapping = aes(x=Lag, xend=Lag, y=0, yend=PACF)) +
    geom_hline(yintercept = Signif_Threshold, color = 'blue', linetype = 'dashed') +
    geom_hline(yintercept = -Signif_Threshold, color = 'blue', linetype = 'dashed') +
    scale_x_continuous(name = "Lag (Days)", limits = c(0, 30)) +
    scale_y_continuous(name = "PACF", limits = c(-0.1, 1)) +
    labs(title = paste0("  PACF Reanalysis Data PC-1 ",Field_Name)) +
    theme_bw() +
    theme(plot.title = element_text(size=10))  
  
  
  out = list(p1=p1,p2=p2)
  
                            

  
}




#______________________________________________________________________________#
###Function to get the 1st PCA and its PACF
###Input
#   1. Wind Field (Fld1)
#   2. Solar Field (Fld2)
#   3. Grid Locations (Grid)
#   4. Start Date (start_date) 

###Output
#   1. 4x4 map of seasonal correlation plots
#   2. 


get_seasonal_corr <- function(Fld1, Fld2, Grid, start_date){
  
  #Load Packages/Dependencies
  library(ggplot2)
  library(gridExtra)
  library(usmap)
  library(ggpubr)
  
  #Hyper-Parameters
  st_date <- as.Date(start_date, format="%m-%d-%Y")
  time_stamps <- seq(st_date, by = "day", length.out =nrow(Fld1))
  month_stamps <- as.numeric(format(time_stamps,"%m"))
  
  #Spatial-Hyperparameters
  world <- map_data("world")
  us <- map_data("state")
  
  #-----------------------------------------------------------------------#
  ###DJF
  season <- c(12,1,2)
  sn <- "DJF"
  
  Fld_s1 <- as.data.frame(Fld1); Fld_s1$Month <- month_stamps
  Fld_s1 <- filter(Fld_s1, Month %in% season); Fld_s1$Month <- NULL
  
  Fld_s2 <- as.data.frame(Fld2); Fld_s2$Month <- month_stamps
  Fld_s2 <- filter(Fld_s2, Month %in% season); Fld_s2$Month <- NULL
  
  #Data Correlation
  M <- rep(NA, ncol(Fld_s1))
  for(i in 1:ncol(Fld_s1)){ 
    M[i] <- cor(Fld_s1[,i],Fld_s2[,i])
  }
  
  #Spatial Correlation for Median Simulation
  Grid_Plot <- data.frame(lon = Grid$lon-360,lat = Grid$lat, Corr = unlist(M))
  
  p1 <- ggplot() + 
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
    labs(title = paste0("  Data Correlation - ", sn)) +
    theme_bw() +
    theme(axis.text=element_text(size=0),
          axis.title=element_text(size=0),
          axis.ticks = element_blank(),
          legend.key.height  = unit(1.75, "cm"),
          legend.text=element_text(size=15),
          legend.title=element_text(size=15))
  
  #-----------------------------------------------------------------------#
  ###MAM
  season <- c(3,4,5)
  sn <- "MAM"
  
  Fld_s1 <- as.data.frame(Fld1); Fld_s1$Month <- month_stamps
  Fld_s1 <- filter(Fld_s1, Month %in% season); Fld_s1$Month <- NULL
  
  Fld_s2 <- as.data.frame(Fld2); Fld_s2$Month <- month_stamps
  Fld_s2 <- filter(Fld_s2, Month %in% season); Fld_s2$Month <- NULL
  
  #Data Correlation
  M <- rep(NA, ncol(Fld_s1))
  for(i in 1:ncol(Fld_s1)){ 
    M[i] <- cor(Fld_s1[,i],Fld_s2[,i])
  }
  
  #Spatial Correlation for Median Simulation
  Grid_Plot <- data.frame(lon = Grid$lon-360,lat = Grid$lat, Corr = unlist(M))
  
  p2 <- ggplot() + 
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
    labs(title = paste0("  Data Correlation - ", sn)) +
    theme_bw() +
    theme(axis.text=element_text(size=0),
          axis.title=element_text(size=0),
          axis.ticks = element_blank(),
          legend.key.height  = unit(1.75, "cm"),
          legend.text=element_text(size=15),
          legend.title=element_text(size=15))
  
  
  #-----------------------------------------------------------------------#
  ###JJA
  season <- c(6,7,8)
  sn <- "JJA"
  
  Fld_s1 <- as.data.frame(Fld1); Fld_s1$Month <- month_stamps
  Fld_s1 <- filter(Fld_s1, Month %in% season); Fld_s1$Month <- NULL
  
  Fld_s2 <- as.data.frame(Fld2); Fld_s2$Month <- month_stamps
  Fld_s2 <- filter(Fld_s2, Month %in% season); Fld_s2$Month <- NULL
  
  #Data Correlation
  M <- rep(NA, ncol(Fld_s1))
  for(i in 1:ncol(Fld_s1)){ 
    M[i] <- cor(Fld_s1[,i],Fld_s2[,i])
  }
  
  #Spatial Correlation for Median Simulation
  Grid_Plot <- data.frame(lon = Grid$lon-360,lat = Grid$lat, Corr = unlist(M))
  
  p3 <- ggplot() + 
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
    labs(title = paste0("  Data Correlation - ", sn)) +
    theme_bw() +
    theme(axis.text=element_text(size=0),
          axis.title=element_text(size=0),
          axis.ticks = element_blank(),
          legend.key.height  = unit(1.75, "cm"),
          legend.text=element_text(size=15),
          legend.title=element_text(size=15))
  
  #-----------------------------------------------------------------------#
  ###SON
  season <- c(9,10,11)
  sn <- "SON"
  
  Fld_s1 <- as.data.frame(Fld1); Fld_s1$Month <- month_stamps
  Fld_s1 <- filter(Fld_s1, Month %in% season); Fld_s1$Month <- NULL
  
  Fld_s2 <- as.data.frame(Fld2); Fld_s2$Month <- month_stamps
  Fld_s2 <- filter(Fld_s2, Month %in% season); Fld_s2$Month <- NULL
  
  #Data Correlation
  M <- rep(NA, ncol(Fld_s1))
  for(i in 1:ncol(Fld_s1)){ 
    M[i] <- cor(Fld_s1[,i],Fld_s2[,i])
  }
  
  #Spatial Correlation for Median Simulation
  Grid_Plot <- data.frame(lon = Grid$lon-360,lat = Grid$lat, Corr = unlist(M))
  
  p4 <- ggplot() + 
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
    labs(title = paste0("  Data Correlation - ", sn)) +
    theme_bw() +
    theme(axis.text=element_text(size=0),
          axis.title=element_text(size=0),
          axis.ticks = element_blank(),
          legend.key.height  = unit(1.75, "cm"),
          legend.text=element_text(size=15),
          legend.title=element_text(size=15))
  
  
  
  ###Print the plots
  print(ggarrange(p1, p2, p3, p4, 
                  nrow = 2, ncol = 2, 
                  common.legend = TRUE, legend = "right",
                  labels = c("A", "B", "C", "D")))
  
}

#______________________________________________________________________________#
###Running the Functions###



pdf("figures/patterns_figures/Data_Characteristics.pdf", compress = TRUE)

#------------------------------------------------------------------------------#
###Get the Wind Grid-Wise Variation Plots
p_wind <- Get_Variation(Dat = WP, 
                        Grid = grid_locs, 
                        start_date = "01-01-1950",
                        frac = 1,
                        Field_Title = "Wind",
                        leg_title = "CF")


#Get the Solar Grid-Wise Variation
p_solar <- Get_Variation(Dat = ssrd,
                         Grid = grid_locs, 
                         start_date = "01-01-1950",
                         frac = 1,
                         Field_Title = "Solar",
                         leg_title = "CF")

ggarrange(p_wind$p1, p_wind$p2, p_solar$p1, p_solar$p2,
          nrow = 2, ncol = 2,
          labels = c("A", "B", "C", "D"))

ggarrange(p_wind$p3, p_wind$p4, p_solar$p3, p_solar$p4,
          nrow = 2, ncol = 2,
          labels = c("A", "B", "C", "D"))


#------------------------------------------------------------------------------#
###Getting the ACF Correlation Plots
p_wind <- get_pca_pacf(Fld = WP, Field_Name = "Wind")
p_solar <- get_pca_pacf(Fld = ssrd, Field_Name = "Solar")

ggarrange(p_wind$p1, p_wind$p2, p_solar$p1, p_solar$p2,
          nrow = 2, ncol = 2,
          labels = c("A", "B", "C", "D"))




#Get the Seasonal Wind-Solar Correlation
get_seasonal_corr(Fld1 = WP, Fld2 = ssrd,
                  Grid = grid_locs, start_date = "01-01-1950")





dev.off()
