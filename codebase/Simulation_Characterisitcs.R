#______________________________________________________________________________#
###Simulation Characteristics across the Texas Interconnect###

###Script to generate plots for
#1. KDE of Aggregated Wind and Solar across TI.
#2. Energy Droughts at various thresholds.
#3. Annual Exceedances at various severity-duration-thresholds 

###Data Inputs
#1. Wind Data
#2. Solar Data
#3. Grid Locations
#4. KSTS Simulations
#5. KNN Simulations

###Outputs
#1. KDE Plots of Aggregated Production. 
#2. Severity-Duration plots for Energy Droughts.
#3. Annual Exceedance Computations.



#______________________________________________________________________________
###Set-up working directory###
setwd("~/GitHub/LDS-Inferences") #Code for personal device

###Load Packages and Dependencies
library(maps)       
library(dplyr)
library(ggplot2)
library(foreach)    #Parallel Execution
library(doParallel) #Backend to foreach



#______________________________________________________________________________#
###Reading the Data###

#Load the Grid Locations
grid_locs <- read.csv("data/ERCOT_0_5_deg_lat_lon_index_key.csv", 
                      header = TRUE, sep=",")

#Load the Downward Surface Solar Radiation
ssrd <- read.table("data/ERCOT_Solar_Rad_Daily.txt",  
                   sep =" ", 
                   header = TRUE) 

#Load the Wind Power Data
WP <- read.table("data/ERCOT_Wind_Power_Daily.txt",  
                 sep =" ", 
                 header = TRUE) 

###Subset
n <- 2*365
WP <- WP[1:n,]
ssrd <- ssrd[1:n,]


#______________________________________________________________________________#
###Read the simulations###

#KSTS
wind_ksts <- solar_ksts <- list()
load("simulations/KSTS_Joint_Simulations.RData")
nl <- length(ynew_results)
for(i in 1:nl){
  wind_ksts[[i]] <- ynew_results[[i]]$WPnew
  solar_ksts[[i]] <- ynew_results[[i]]$SSnew
  ynew_results[[i]] <- 1
}
ynew_results <- NULL

#Knn
wind_knn <- solar_knn <- list()
load("simulations/KNN_Joint_Simulations.RData")
nl <- length(ynew_results)
for(i in 1:nl){
  wind_knn[[i]] <- ynew_results[[i]]$WPnew
  solar_knn[[i]] <- ynew_results[[i]]$SSnew
  ynew_results[[i]] <- 1
}
ynew_results <- NULL



#______________________________________________________________________________#
###Function 1### 
#Objective - Get the Total Aggregate Energy Production for Wind and Solar.
###Inputs
#   1. Original Field (Dat)
#   2. Field Name (Field)
#   3. KSTS Simulations (ksts)
#   4. KNN Simulations (knn)
#   5. Land Allocation (land_alc)
#   6. Max Capacity (max_cap)
###Output
#   1. All DOYs in the moving window

Get_Total_Energy <- function(Dat, Field, ksts, knn, land_alc, max_cap){
  
  #Hyper-Parameters
  nsim_knn <- length(knn)
  nsim_ksts <- length(ksts)
  
  #COnvert W/sq-m to MWhr
  t_fac <- 24*land_alc/(10^6)   #hrs*sq-m/M
  tx <- rowSums(Dat)*t_fac      #MWhr
  tx <- tx/max_cap              #Division by Max-Capacity
  
  #Compute the KDE on the Reanalysis Data
  og_pdf <- density(tx, from = 0, to = 1)
  
  #KSTS
  sim_ksts <- matrix(NA, ncol = nsim_ksts, nrow = length(og_pdf$x)) #Storing the Simulated CDF's
  for(j in 1:nsim_ksts){
    #Computing each CDF
    sim <- as.data.frame(ksts[[j]])
    sim <- rowSums(sim)*t_fac/max_cap
    pdf_sim <- density(sim, from = 0, to = 1)
    sim_ksts[,j] <- pdf_sim$y
  }
  lower_ksts <- apply(sim_ksts, 1, function(x) quantile(x, probs=.05))
  upper_ksts <- apply(sim_ksts, 1, function(x) quantile(x, probs=.95))
  median_ksts <- apply(sim_ksts, 1, median)
  
  
  #KNN
  sim_knn <- matrix(NA, ncol = nsim_knn, nrow = length(og_pdf$x)) #Storing the Simulated CDF's
  for(j in 1:nsim_knn){
    #Computing each CDF
    sim <- as.data.frame(knn[[j]])
    sim <- rowSums(sim)*t_fac/max_cap
    pdf_sim <- density(sim, from = 0, to = 1)
    sim_knn[,j] <- pdf_sim$y
  }
  lower_knn <- apply(sim_knn, 1, function(x) quantile(x, probs=.05))
  upper_knn <- apply(sim_knn, 1, function(x) quantile(x, probs=.95))
  median_knn <- apply(sim_knn, 1, median)
  
  
  #Plotting the results
  par(mfrow=c(1,1), mar = c(6,6,4,1))
  plot(og_pdf$x, og_pdf$y, type='l',col='red',
       lwd = 3, main = paste0("Daily Generation Potential - ", Field), 
       xlab = "Fraction of Maximum Production", ylab = "Density - f(x) ",
       ylim = c(0, max(upper_knn)),
       xlim = c(0,1),
       cex.lab = 1.5,
       cex.main = 1.75)
  #KSTS
  polygon(c(og_pdf$x,rev(og_pdf$x)),c(lower_ksts,rev(upper_ksts)),col="#af8dc3")
  
  #KNN
  polygon(c(og_pdf$x,rev(og_pdf$x)),c(lower_knn,rev(upper_knn)),col="#7fbf7b")
  
  #Data-Line
  lines(og_pdf$x, og_pdf$y, col = "black", lwd = 3)
  
  legend('topleft', legend = c("Data", "KSTS", "KNN"),
         col = c('black',"#af8dc3", "#7fbf7b"), 
         lty = c(1, NA, NA),  lwd = c(3, NA, NA),
         pch = c(NA, 15, 15), pt.cex = c(NA,2,2),
         cex = 1.25,
         bty = "n")
  
}

#______________________________________________________________________________#
###Running the Total Energy Plots###

Get_Total_Energy(Dat = WP, 
                 Field = "Wind",
                 ksts = wind_ksts,
                 knn = wind_knn,
                 land_alc = 4*7*90*90, 
                 max_cap = 2*24*216) #Installed Turbine Capacity

Get_Total_Energy(Dat = ssrd, 
                 Field = "Solar",
                 ksts = solar_ksts,
                 knn = solar_knn,
                 land_alc = 4*7*90*90, 
                 max_cap = 24*4*7*90*90*max(rowSums(ssrd))/10^6) #Max Solar

#______________________________________________________________________________#
###ENERGY DROUGHT AND EXCEEDANCE PROBABILITY###

#This section requires all the sites to be normalized first.
#Therefore the earlier loaded simulations are discarded and re-loaded again.
solar_knn <- wind_knn <- NULL
solar_ksts <- wind_ksts <- NULL

###Read the simulations
#KSTS
comb_ksts <- list()
load("simulations/KSTS_Joint_Simulations.RData")
nl <-  length(ynew_results)
for(i in 1:nl){
  comb_ksts[[i]] <- cbind(ynew_results[[i]]$WPnew, ynew_results[[i]]$SSnew)
  comb_ksts[[i]] <- apply(comb_ksts[[i]], 2, function(x){x/max(x)})
  ynew_results[[i]] <- 1
  
}
ynew_results <- NULL

#KNN
comb_knn <- list()
load("simulations/KSTS_Joint_Simulations.RData")
nl <- length(ynew_results)
for(i in 1:nl){
  comb_knn[[i]] <- cbind(ynew_results[[i]]$WPnew, ynew_results[[i]]$SSnew)
  comb_knn[[i]] <- apply(comb_knn[[i]], 2, function(x){x/max(x)})
  ynew_results[[i]] <- 1
  
}
ynew_results <- NULL


#Data Normalization
Fld <- as.matrix(cbind(WP,ssrd))
Fld <- apply(Fld, 2, function(x){x/max(x)})



#______________________________________________________________________________#
###Function 2###
#Objective - To get Severity vs Duration plots for a given thresholds
###Inputs
#   1. Original Data (True_Data)
#   2. Field Name (Field_Name)
#   3. KSTS Simulations (Sims_KSTS)
#   4. KNN Simulations (Sims_KNN)
#   5. Threshold (thresh)
#   6. Start_Date

get_energy_droughts <- function(True_Data, 
                                Field_Name,
                                Sims_KSTS,
                                Sims_KNN,
                                thresh,
                                Start_Date)
{
  #Source the function.
  source("functions/Get_Severity_Duration.R")
  
  #Load Dependencies
  library(ggplot2)
  library(lubridate)
  library(cowplot)
  library(gridExtra)
  
  #Hyper-Parameters
  st_date = Start_Date
  nsim_ksts <- length(Sims_KSTS)
  nsim_knn <- length(Sims_KNN)
  mean_daily <- mean(rowSums(True_Data))
  
  #Compute Severity and Duration for Data
  Data_SD <- get_Severity_Duration(Data = True_Data, Thresh = thresh,
                                   start_date = st_date, Type = "Data")
  
  #KSTS
  cores=detectCores()
  registerDoParallel(cores)
  Sim_SD_KSTS <- foreach(m = 1:nsim_ksts, .verbose = TRUE) %dopar% {
    source("functions/Get_Severity_Duration.R")
    get_Severity_Duration(Data = Sims_KSTS[[m]], Thresh = thresh,
                          start_date = st_date, Type = "KSTS")
  }
  stopImplicitCluster()
  Sims_KSTS <- NULL
  
  #KNN
  cores=detectCores()
  registerDoParallel(cores)
  Sim_SD_KNN <- foreach(m = 1:nsim_knn, .verbose = TRUE) %dopar% {
    source("functions/Get_Severity_Duration.R")
    get_Severity_Duration(Data = Sims_KNN[[m]], Thresh = thresh,
                          start_date = st_date, Type = "KNN")
  }
  Sims_KNN <- NULL
  stopImplicitCluster()
  
  #Convert to a Data-Frame Structure
  Sims_KNN <- bind_rows(lapply(Sim_SD_KNN,data.frame))
  Sims_KSTS <- bind_rows(lapply(Sim_SD_KSTS,data.frame))
  Sev_Dur <- rbind(Sims_KSTS,Data_SD,Sims_KNN)
  Sev_Dur$Severity <- Sev_Dur$Severity/mean_daily
  
  #Remove the non-events.
  Sev_Dur<- Sev_Dur %>% filter(Duration > 0)
  
  group.colors <- c(KSTS ="#af8dc3", Data = "#000000", KNN = "#7fbf7b")
  
  p_main <- ggplot(Sev_Dur) + 
    geom_point(aes(log10(Duration), log10(Severity), color = Type), alpha = 0.75) +
    ggtitle(paste0("Duration vs Severity \n Threshold - ", thresh*100,"th Percentile")) + 
    scale_x_continuous(name = "Duration (Days)", labels = scales::math_format(10^.x)) +
    scale_y_continuous(name = "Severity (Multiples of Mean Production)", labels = scales::math_format(10^.x)) +
    theme_bw() +
    theme(legend.text=element_text(size=20),
          legend.title=element_text(size=15),
          axis.text=element_text(size=15),
          axis.title=element_text(size=20),
          plot.title = element_text(size=20)) +
    guides(colour = guide_legend(override.aes = list(size=5))) +
    scale_color_manual(values=group.colors)
  
  xbox <- axis_canvas(p_main, axis = "x", coord_flip = TRUE) + 
    geom_boxplot(data = Sev_Dur, aes(y = log10(Duration), x = factor(Type), color = factor(Type))) + 
    scale_x_discrete() + coord_flip() +
    scale_color_manual(values=group.colors)
  
  ybox <- axis_canvas(p_main, axis = "y") + 
    geom_boxplot(data = Sev_Dur, aes(y = log10(Severity), x = factor(Type), color = factor(Type))) +
    scale_x_discrete() +
    scale_color_manual(values=group.colors)
  
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


#______________________________________________________________________________#
###Running the Energy Droughts Plotting Function at multiple thresholds###

get_energy_droughts(True_Data = Fld,
                    Field_Name = " ",
                    Sims_KSTS = comb_ksts,
                    Sims_KNN = comb_knn,
                    thresh = 0.15,
                    Start_Date = "01-01-1979")


get_energy_droughts(True_Data = Fld,
                    Field_Name = " ",
                    Sims_KSTS = comb_ksts,
                    Sims_KNN = comb_knn,
                    thresh = 0.20,
                    Start_Date = "01-01-1979")

get_energy_droughts(True_Data = Fld,
                    Field_Name = " ",
                    Sims_KSTS = comb_ksts,
                    Sims_KNN = comb_knn,
                    thresh = 0.25,
                    Start_Date = "01-01-1979")

get_energy_droughts(True_Data = Fld,
                    Field_Name = " ",
                    Sims_KSTS = comb_ksts,
                    Sims_KNN = comb_knn,
                    thresh = 0.30,
                    Start_Date = "01-01-1979")

get_energy_droughts(True_Data = Fld,
                    Field_Name = " ",
                    Sims_KSTS = comb_ksts,
                    Sims_KNN = comb_knn,
                    thresh = 0.35,
                    Start_Date = "01-01-1979")

get_energy_droughts(True_Data = Fld,
                    Field_Name = " ",
                    Sims_KSTS = comb_ksts,
                    Sims_KNN = comb_knn,
                    thresh = 0.40,
                    Start_Date = "01-01-1979")

get_energy_droughts(True_Data = Fld,
                    Field_Name = " ",
                    Sims_KSTS = comb_ksts,
                    Sims_KNN = comb_knn,
                    thresh = 0.45,
                    Start_Date = "01-01-1979")

