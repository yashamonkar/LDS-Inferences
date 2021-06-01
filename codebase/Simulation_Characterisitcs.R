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
###Function 1 
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
