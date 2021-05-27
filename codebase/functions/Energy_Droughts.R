#Code to get the severity vs Duration plots for Data, KNN and KSTS.


#setwd("~/ERCOT") #Code for personal device


#______________________________________________________________________________#
.libPaths("/rigel/cwc/users/yva2000/rpackages/")
###Load Packages       
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
comb_ksts <- list()
load("Simulations/Joint_Raw_Simulations.RData")
nl <-  length(ynew_results)
for(i in 1:nl){
  comb_ksts[[i]] <- cbind(ynew_results[[i]]$WPnew, ynew_results[[i]]$SSnew)
  comb_ksts[[i]] <- apply(comb_ksts[[i]], 2, function(x){x/max(x)})
  ynew_results[[i]] <- 1
  
}
ynew_results <- NULL

#Knn
comb_knn <- list()
load("Simulations/KNN_Joint_Raw_Simulations.RData")
nl <- length(ynew_results)
for(i in 1:nl){
  comb_knn[[i]] <- cbind(ynew_results[[i]]$WPnew, ynew_results[[i]]$SSnew)
  comb_knn[[i]] <- apply(comb_knn[[i]], 2, function(x){x/max(x)})
  ynew_results[[i]] <- 1
  
}
ynew_results <- NULL

#______________________________________________________________________________#
###Code to visualize the Energy Droughts###
#Note:- Time Step is Day.


###Input
#1.  Data Field - Fld
#2.  Simulations 
#3.  Start Date
#4.  Threshold
#5.  Field Name


Fld <- as.matrix(cbind(WP,ssrd))
Fld <- apply(Fld, 2, function(x){x/max(x)})

###Output
#1.  ggplot of log(Severity) vs Duration with Marginals


#______________________________________________________________________________#
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
  
  #Sev_Dur$Severity <- log(Sev_Dur$Severity)
  #Sev_Dur$Duration <- log(Sev_Dur$Duration)
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

