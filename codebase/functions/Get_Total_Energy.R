#______________________________________________________________________________#
#Code to get total energy plot for each wind and solar seperately. 
#Include simulations for KNN and KSTS


#Input
#1. Data.
#2. Field Name
#3. KSTS Simulations
#4. KNN Simulations
#5. Land Allocation


#Output
#1. Mean power at each Site.
#2. Mean power at each Site by season
#3. Variation at each Site.
#4, Variation at each Site by season.


setwd("~/ERCOT") #Code for personal device

#______________________________________________________________________________#
#.libPaths("/rigel/cwc/users/yva2000/rpackages/")
###Load Packages
library(maps)       
library(dplyr)
library(ggplot2)
library(logspline)
library(foreach)    #Parallel Execution
library(doParallel) #Backend to foreach


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
###Read the simulations
#KSTS
wind_ksts <- solar_ksts <- list()
load("Simulations/Joint_Raw_Simulations.RData")
nl <-  4 #length(ynew_results)
for(i in 1:nl){
  #wind_ksts[[i]] <- ynew_results[[i]]$WPnew
  solar_ksts[[i]] <- ynew_results[[i]]$SSnew
  ynew_results[[i]] <- 1
}
ynew_results <- NULL

#Knn
wind_knn <- solar_knn <- list()
load("Simulations/KNN_Joint_Raw_Simulations.RData")
nl <- 4 #length(ynew_results)
for(i in 1:nl){
  #wind_knn[[i]] <- ynew_results[[i]]$WPnew
  solar_knn[[i]] <- ynew_results[[i]]$SSnew
  ynew_results[[i]] <- 1
}
ynew_results <- NULL

#______________________________________________________________________________#
Get_Total_Energy <- function(Dat, Field, ksts, knn, land_alc, max_cap){
  
  #Hyper-Parameters
  nsim_knn <- length(knn)
  nsim_ksts <- length(ksts)
  
  #COnvert W/sq-m to MWhr
  t_fac <- 24*land_alc/(10^6)   #hrs*kms/Mega
  
  tx <- rowSums(Dat)*t_fac
  tx <- tx/max_cap
  
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
  #lines(og_pdf$x, median_ksts, lwd = 2)
  #KNN
  polygon(c(og_pdf$x,rev(og_pdf$x)),c(lower_knn,rev(upper_knn)),col="#7fbf7b")
  #lines(og_pdf$x, og_pdf$y, col='red', lwd = 2)
  #Data-Line
  lines(og_pdf$x, og_pdf$y, col = "black", lwd = 3)
  
  legend('topleft', legend = c("Data", "KSTS", "KNN"),
         col = c('black',"#af8dc3", "#7fbf7b"), 
         lty = c(1, NA, NA),  lwd = c(3, NA, NA),
         pch = c(NA, 15, 15), pt.cex = c(NA,2,2),
         cex = 1.25,
         bty = "n")
  
  #legend('topleft', legend = c("Data", "KSTS", "KNN"),
  #      lty = 1, col = c('red','black','blue'), lwd = 3, cex = 1.25)
  
  
  
}





#______________________________________________________________________________#
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
                 land_alc = 10117.1, #2.5 acres in m^2
                 max_cap = 24*10117.1*max(rowSums(ssrd))/10^6) #Max Solar
