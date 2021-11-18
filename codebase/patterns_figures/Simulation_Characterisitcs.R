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
#setwd("~/GitHub/LDS-Inferences") #Code for personal device

###Set-up Library Paths
.libPaths("/rigel/cwc/users/yva2000/rpackages/")

###Load Packages and Dependencies
library(maps)       
library(dplyr)
library(ggplot2)
library(foreach)    #Parallel Execution
library(doParallel) #Backend to foreach
library(gridExtra)
library(cowplot)

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

###Subset
#n <- 5*365
#WP <- WP[1:n,]
#ssrd <- ssrd[1:n,]



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


#Create the Figure PDF
pdf("figures/patterns_figures/Simulation_Characteristics.pdf")

#______________________________________________________________________________#
###Function 1### 
#Objective - Get the Total Aggregate Energy Production for Wind and Solar.
###Inputs
#   1. Original Field (Dat)
#   2. Field Name (Field)
#   3. KSTS Simulations (ksts)
#   4. KNN Simulations (knn)
#   5. Max Capacity (max_cap) (Visual Parameter Limiting the x-axis)
###Output
#   1. All DOYs in the moving window

Get_Total_Energy <- function(Dat, Field, ksts, knn, max_cap){
  
  #Hyper-Parameters
  nsim_knn <- length(knn)
  nsim_ksts <- length(ksts)
  
  #COnvert W/sq-m to MWhr
  tx <- rowSums(Dat)/ncol(Dat)      #Add-Up the capacity factors
  
  #Compute the KDE on the Reanalysis Data
  og_pdf <- density(tx, from = 0, to = 1)
  
  #KSTS
  sim_ksts <- matrix(NA, ncol = nsim_ksts, nrow = length(og_pdf$x)) #Storing the Simulated CDF's
  for(j in 1:nsim_ksts){
    #Computing each CDF
    sim <- as.data.frame(ksts[[j]])
    sim <- rowSums(sim)/ncol(Dat)
    pdf_sim <- density(sim, from = 0, to = 1)
    sim_ksts[,j] <- pdf_sim$y
  }
  lower_ksts <- apply(sim_ksts, 1, function(x) quantile(x, probs=.01))
  upper_ksts <- apply(sim_ksts, 1, function(x) quantile(x, probs=.99))
  median_ksts <- apply(sim_ksts, 1, median)
  
  
  #KNN
  sim_knn <- matrix(NA, ncol = nsim_knn, nrow = length(og_pdf$x)) #Storing the Simulated CDF's
  for(j in 1:nsim_knn){
    #Computing each CDF
    sim <- as.data.frame(knn[[j]])
    sim <- rowSums(sim)/ncol(Dat)
    pdf_sim <- density(sim, from = 0, to = 1)
    sim_knn[,j] <- pdf_sim$y
  }
  lower_knn <- apply(sim_knn, 1, function(x) quantile(x, probs=.01))
  upper_knn <- apply(sim_knn, 1, function(x) quantile(x, probs=.99))
  median_knn <- apply(sim_knn, 1, median)
  
  
  #Plotting Dataset
  Plt_Data <- data.frame(X = rep(og_pdf$x,2),
                         lower = c(lower_knn, lower_ksts),
                         upper = c(upper_knn, upper_ksts),
                         Type = rep(c("KNN","KSTS"), each = length(og_pdf$x)))
  
  OG_PDF <- data.frame(X=og_pdf$x,
                       Y=og_pdf$y)
  
  group.colors <- c(KSTS ="#af8dc3", KNN = "#7fbf7b")
  
  
  
  #Plotting the results
  p <- ggplot(Plt_Data, aes(X)) +
    geom_ribbon(aes(ymin = lower, ymax =  upper, fill = Type), color ='black', alpha = 0.85) +
    geom_line(OG_PDF, mapping = aes(x = X, y = Y), color = 'red', size = 0.9) +
    scale_x_continuous(name = "Fraction of Max Production", limits = c(0, max_cap)) +
    scale_y_continuous(name = "Density - f(x)") +
    labs(title = paste0("  Daily Generation Potential - ", Field)) + 
    theme_bw() +
    theme(legend.text=element_text(size=15),
          legend.title=element_text(size=0),
          plot.title = element_text(size=15),
          legend.key.height  = unit(1.5, "cm"),
          legend.key.width  = unit(1.5, "cm")) +
    guides(fill = guide_legend(override.aes = list(size=2))) +
    scale_fill_manual(values=group.colors)
  
  return(p)
  
  
  
  
}

#______________________________________________________________________________#
###Running the Total Energy Plots###

p1 <- Get_Total_Energy(Dat = WP, 
                 Field = "Wind",
                 ksts = wind_ksts,
                 knn = wind_knn,
                 max_cap = 1)

p2 <- Get_Total_Energy(Dat = ssrd, 
                 Field = "Solar",
                 ksts = solar_ksts,
                 knn = solar_knn,
                 max_cap = 0.4) 

p_total <- plot_grid(p1 + theme(legend.position="none"),
                     p2 + theme(legend.position="none"),
                     nrow =2,
                     labels = c('A', 'B'), 
                     label_size = 12)


legend_b <- get_legend(
  p1 + 
    guides(color = guide_legend(nrow = 1, override.aes = list(size=2))) +
    theme(legend.position = "bottom")
)

# add the legend to the row we made earlier. Give it one-third of 
# the width of one plot (via rel_widths).
#plot_grid(p_total, legend, rel_widths = c(1, .2))
plot_grid(p_total, legend_b, ncol = 1, rel_heights = c(1, .2))

solar_knn <- wind_knn <- NULL
solar_ksts <- wind_ksts <- NULL





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
#   7. Field_type (Joint, Wind or Solar)

###Output
#   1.  Single Severity vs Duration plots with boxplots for marginals.

get_energy_droughts <- function(True_Data, 
                                Field_Name,
                                Sims_KSTS,
                                Sims_KNN,
                                thresh,
                                Start_Date,
                                Field_type)
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
    geom_point(aes(log10(Duration), log10(Severity), color = Type), 
               size = 0.75, alpha = 0.75) +
    ggtitle(paste0("Duration vs Severity \n Threshold - ", thresh*100,"th Percentile")) + 
    scale_x_continuous(name = "Duration (Days)", labels = scales::math_format(10^.x)) +
    scale_y_continuous(name = "Severity", labels = scales::math_format(10^.x)) +
    theme_bw() +
    theme(legend.text=element_text(size=15),
          legend.title=element_text(size=0),
          axis.text=element_text(size=10),
          axis.title=element_text(size=10),
          plot.title = element_text(size=12)) +
    guides(colour = guide_legend(override.aes = list(size=2))) +
    scale_color_manual(values=group.colors)
  
  #Get the legend - Side
  #legend <- get_legend(
  #  p_main + theme(legend.box.margin = margin(0, 0, 0, 12))
  #)
  
  #Get the legend - Bottom
  legend_b <- get_legend(
    p_main + 
      guides(color = guide_legend(nrow = 1, override.aes = list(size=4))) +
      theme(legend.position = "bottom")
  )
  
  
  p_main <- p_main + theme(legend.position="none")
  
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
  
  #print(p3)
  #return(p3)
  out = list(p=p3, legend = legend_b)
  
  
}

#______________________________________________________________________________#
###Function 3
###Function to fit a log-fit to the output of severity duration.
###Inputs
#1. Data (Dat)
#2. KSTS Simulations (Sims)
#3. Threshold Percentile (Thresh)
#4. Desired Severity [2] e.g. c(1.25,1.5)
#5. Desired Duration [3] e.g. c(20,25,30)
#6. Field_type (Joint, Wind or Solar)

###Outputs
#1. Panel Plot


get_aep_plot <- function(Dat, Sims, Thresh, Severities, Durations, 
                         Field_type){
  
  #Hyper-Parameters
  thresh = Thresh
  st_date = "01-01-1950"
  mean_daily <- mean(rowSums(Dat))
  nl <- length(Sims)
  
  #Set-up Grid of Values
  Dur = Durations
  Sev = Severities
  ngrids <- length(Sev)*length(Dur)
  Sev_Dur <- data.frame(Duration = rep(Dur, each = length(Sev)),
                        Severity = rep(Sev, length(Dur)))
  
  
  #---------------------------------------------------------------------------#
  ###Function to fit a log-fit to the output of severity duration.
  ###Inputs of Function
  #1. Severity-Duration for all Data. (ED)
  #2. Desired Duration (Dur)
  #3. Desired Severity (Sev)
  #4. Mean Daily Value to facilitate comparision (MD)
  #5. Number of Years (N_yr)
  
  ###Output 
  #1. Annual Probability
  
  
  get_annual_prob <- function(ED, Dur, Sev, MD, N_yr){
    
    #Load Packages
    library(locfit)
    
    #Cleaning
    ED$Severity <- ED$Severity/MD
    ED$Severity <- log10(ED$Severity)
    ED$Duration <- log10(ED$Duration)
    
    #Compute Exceedances
    n_exceed <- list()
    for(i in 1:nrow(ED)){
      t_count <-  which(ED$Severity > ED$Severity[i] &
                          ED$Duration > ED$Duration[i])
      n_exceed[[i]] <- length(t_count)
    }
    ED$exc <- unlist(n_exceed)
    
    ###Local Regression for the Simulations
    fit <- locfit(exc ~ Duration + Severity,family = "poisson",
                  data=ED, alpha = 0.33)
    
    ###Predict the exceedances
    exceed <- predict(fit, matrix(c(log10(Dur),log10(Sev)), nrow = 1))
    
    ##Compute the probability of annual exceedances
    p_annual <- exceed/N_yr
    
    return(p_annual)
    
  }
  
  
  #---------------------------------------------------------------------------#
  ###Compute the severity and drought for the data.
  #Source the function.
  source("functions/Get_Severity_Duration.R")
  
  #Compute Severity and Duration for Data
  Data_SD <- get_Severity_Duration(Data = Dat, Thresh = thresh,
                                   start_date = st_date, Type = "Data")
  data_pexc <- list()
  for(i in 1:nrow(Sev_Dur)){
    data_pexc[[i]] <- 100*get_annual_prob(ED = Data_SD, 
                                          Dur = Sev_Dur$Duration[i],
                                          Sev = Sev_Dur$Severity[i], 
                                          MD = mean_daily, N_yr = 71)
  }
  data_pexc <- unlist(data_pexc)
  
  
  #---------------------------------------------------------------------------#
  ###Compute the severity and drought for the the simulations.
  cores=detectCores()
  registerDoParallel(cores)
  sim_pexc <- list()
  sim_pexc <- foreach(m = 1:nl, .verbose = TRUE) %dopar% {
    
    ###Load Functions
    source("functions/Get_Severity_Duration.R")
    
    #Compute Severity and Duration for Simulations
    KSTS_SD <- get_Severity_Duration(Data = Sims[[m]], Thresh = thresh,
                                     start_date = st_date, Type = "Sims")
    
    get_annual_prob <- function(ED, Dur, Sev, MD, N_yr){
      
      #Load Packages
      library(locfit)
      
      #Cleaning
      ED$Severity <- ED$Severity/MD
      ED$Severity <- log10(ED$Severity)
      ED$Duration <- log10(ED$Duration)
      
      #Compute Exceedances
      n_exceed <- list()
      for(i in 1:nrow(ED)){
        t_count <-  which(ED$Severity > ED$Severity[i] &
                            ED$Duration > ED$Duration[i])
        n_exceed[[i]] <- length(t_count)
      }
      ED$exc <- unlist(n_exceed)
      
      ###Local Regression for the Simulations
      fit <- locfit(exc ~ Duration + Severity,family = "poisson",
                    data=ED, alpha = 0.33)
      
      ###Predict the exceedances
      exceed <- predict(fit, matrix(c(log10(Dur),log10(Sev)), nrow = 1))
      
      ##Compute the probability of annual exceedances
      p_annual <- exceed/N_yr
      
      return(p_annual)
      
    }
    
    temp_pexc <- list()
    for(i in 1:nrow(Sev_Dur)){
      temp_pexc[[i]] <- 100*get_annual_prob(ED = KSTS_SD, 
                                            Dur = Sev_Dur$Duration[i],
                                            Sev = Sev_Dur$Severity[i], 
                                            MD = mean_daily, N_yr = 71)
    }
    
    sim_pexc[[m]] <- unlist(temp_pexc)
  }
  stopImplicitCluster()
  
  
  #---------------------------------------------------------------------------#
  ###Plotting the results
  
  #Generate Labels
  class_labels <- list()
  t = 1
  for(i in 1:length(Dur)){
    for(j in 1:length(Sev)){
      class_labels[[t]] <- paste0("Dur = ", Dur[i], " days")
      t = t+1
    }
  }
  class_labels <- unlist(class_labels)
  
  #Simulations
  sims_plot <- bind_rows(lapply(sim_pexc,data.frame))
  sims_plot$Class <- as.factor(rep(1:ngrids, nl))
  colnames(sims_plot) <- c("Exceedance", "Class")
  sims_plot$Severity <- 100*rep(rep(Sev, length(Dur)), nl)
  sims_plot$Severity <- paste0("Severity = ", sims_plot$Severity," %")
  sims_plot$Duration <- rep(rep(Dur, each = length(Sev)), nl)
  sims_plot$Duration <- as.factor(sims_plot$Duration)
  
  
  #Data
  data_plot <- data.frame(Exceedance = data_pexc,
                          Class = as.factor(1:ngrids),
                          Severity = 100*rep(Sev, length(Dur)))
  data_plot$Severity <- paste0("Severity = ", data_plot$Severity," %")
  data_plot$Duration <- rep(Dur, each = length(Sev))
  data_plot$Duration <- as.factor(data_plot$Duration)
  
  
  p <- ggplot(sims_plot, aes(x=Duration, y=Exceedance)) + 
    geom_boxplot() +
    ggtitle(paste0("Annual Exceedances - ",Field_type," \n Threshold - ", thresh*100,
                   "th Percentile ")) +
    ylab("Annual Exceedance Probabilities (%)") +
    xlab("Duration (Days)") +
    theme_bw() +
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=10),
          plot.title = element_text(size=12),
          strip.text = element_text(size=10))
  
  p <- p + facet_grid(Severity ~ ., scale = "free")
  p <- p + geom_point(data_plot, mapping = aes(x=Duration, y=Exceedance), color = 'red', size = 3)
  
  
  return(p)
  
  
}


#______________________________________________________________________________#
#______________________________________________________________________________#
#______________________________________________________________________________#
###ENERGY DROUGHT AND EXCEEDANCE PROBABILITY###
#Joint Wind and Solar

#This section requires all the sites to be normalized first.

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
load("simulations/KNN_Joint_Simulations.RData")
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

#------------------------------------------------------------------------------#
###Running the Energy Droughts Plotting Function at multiple thresholds###



p1 <- get_energy_droughts(True_Data = Fld,
                    Field_Name = " ",
                    Sims_KSTS = comb_ksts,
                    Sims_KNN = comb_knn,
                    thresh = 0.25,
                    Start_Date = "01-01-1950",
                    Field_type = "Joint Fields")

p2 <- get_energy_droughts(True_Data = Fld,
                    Field_Name = " ",
                    Sims_KSTS = comb_ksts,
                    Sims_KNN = comb_knn,
                    thresh = 0.30,
                    Start_Date = "01-01-1950",
                    Field_type = "Joint Fields")

p3 <- get_energy_droughts(True_Data = Fld,
                    Field_Name = " ",
                    Sims_KSTS = comb_ksts,
                    Sims_KNN = comb_knn,
                    thresh = 0.35,
                    Start_Date = "01-01-1950",
                    Field_type = "Joint Fields")

p4 <- get_energy_droughts(True_Data = Fld,
                    Field_Name = " ",
                    Sims_KSTS = comb_ksts,
                    Sims_KNN = comb_knn,
                    thresh = 0.40,
                    Start_Date = "01-01-1950",
                    Field_type = "Joint Fields")

p_droughts <- plot_grid(p1$p + theme(legend.position="none"),
                        p2$p+ theme(legend.position="none"),
                        p3$p+ theme(legend.position="none"),
                        p4$p+ theme(legend.position="none"),
                        nrow =2,
                        labels = c('A', 'B', 'C', 'D'), 
                        label_size = 12)


###Adding the legend to the plot
#plot_grid(p_droughts, legend, rel_widths = c(3, .4))
plot_grid(p_droughts, p1$legend, ncol = 1, rel_heights = c(1, .1))


#------------------------------------------------------------------------------#
###Run the function to compute exceedance probabilities###

p1 <- get_aep_plot(Dat = Fld, Sims = comb_ksts, Thresh = 0.25, 
             Severities = c(1.25,1.5), Durations = c(20,25,30),
             Field_type = "Joint Fields")


p2 <- get_aep_plot(Dat = Fld, Sims = comb_ksts, Thresh = 0.30, 
             Severities = c(2,5), Durations = c(30,45,60),
             Field_type = "Joint Fields")


p3 <- get_aep_plot(Dat = Fld, Sims = comb_ksts, Thresh = 0.15, 
             Severities = c(0.75,1), Durations = c(10,15,20),
             Field_type = "Joint Fields")

p4 <- get_aep_plot(Dat = Fld, Sims = comb_ksts, Thresh = 0.20, 
             Severities = c(1,1.25), Durations = c(10,20,30),
             Field_type = "Joint Fields")


plot_grid(p1,NULL,p2,
          nrow =1,
          labels = c('A', '' ,'B'), 
          label_size = 15,
          rel_widths = c(1, 0.1, 1))

plot_grid(p3,NULL,p4,
          nrow =1,
          labels = c('A','' ,'B'), 
          label_size = 15,
          rel_widths = c(1, 0.1, 1))


comb_ksts <- comb_knn <- Fld <- NULL




#______________________________________________________________________________#
#______________________________________________________________________________#
#______________________________________________________________________________#
###Just Wind

###Read the Wind simulations
#KSTS
comb_ksts <- list()
load("simulations/KSTS_Joint_Simulations.RData")
nl <-  length(ynew_results)
for(i in 1:nl){
  comb_ksts[[i]] <- ynew_results[[i]]$WPnew
  ynew_results[[i]] <- 1
  
}
ynew_results <- NULL

#KNN
comb_knn <- list()
load("simulations/KNN_Joint_Simulations.RData")
nl <- length(ynew_results)
for(i in 1:nl){
  comb_knn[[i]] <- ynew_results[[i]]$WPnew
  ynew_results[[i]] <- 1
  
}
ynew_results <- NULL

#Data Normalization
Fld <- as.matrix(WP)

#------------------------------------------------------------------------------#
###Running the Energy Droughts Plotting Function at multiple thresholds###


p1 <- get_energy_droughts(True_Data = Fld,
                    Field_Name = " ",
                    Sims_KSTS = comb_ksts,
                    Sims_KNN = comb_knn,
                    thresh = 0.20,
                    Start_Date = "01-01-1950",
                    Field_type = "Wind")


p2 <- get_aep_plot(Dat = Fld, Sims = comb_ksts, Thresh = 0.20, 
             Severities = c(1,1.25), Durations = c(10,20,30),
             Field_type = "Wind")



comb_ksts <- comb_knn <- Fld <- NULL



#______________________________________________________________________________#
#______________________________________________________________________________#
#______________________________________________________________________________#
###Just Solar

###Read the Wind simulations
#KSTS
comb_ksts <- list()
load("simulations/KSTS_Joint_Simulations.RData")
nl <-  length(ynew_results)
for(i in 1:nl){
  comb_ksts[[i]] <- ynew_results[[i]]$SSnew
  ynew_results[[i]] <- 1
  
}
ynew_results <- NULL

#KNN
comb_knn <- list()
load("simulations/KNN_Joint_Simulations.RData")
nl <- length(ynew_results)
for(i in 1:nl){
  comb_knn[[i]] <- ynew_results[[i]]$SSnew
  ynew_results[[i]] <- 1
  
}
ynew_results <- NULL

#Data Normalization
Fld <- as.matrix(ssrd)

#------------------------------------------------------------------------------#
###Running the Energy Droughts Plotting Function at multiple thresholds###


p3 <- get_energy_droughts(True_Data = Fld,
                    Field_Name = " ",
                    Sims_KSTS = comb_ksts,
                    Sims_KNN = comb_knn,
                    thresh = 0.20,
                    Start_Date = "01-01-1950",
                    Field_type = "Solar")


#------------------------------------------------------------------------------#
###Run the function to compute exceedance probabilities###

p4 <- get_aep_plot(Dat = Fld, Sims = comb_ksts, Thresh = 0.20, 
             Severities = c(2,5), Durations = c(30,45,60),
             Field_type = "Solar")


p_droughts <- plot_grid(p1$p,p3$p,
                    nrow =2,
                    labels = c('A', 'C'),
                    label_size = 12)

p_droughts <- plot_grid(p_droughts, p1$legend, ncol = 1, rel_heights = c(1, .1))

p_excd <- plot_grid(p2,p4,
                    nrow =2,
                    labels = c('B', 'D'),
                    label_size = 12)

plot_grid(p_droughts,p_excd,
          nrow = 1)


comb_ksts <- comb_knn <- Fld <- NULL





#Close the PDF
dev.off()