#Code to plot the probability of occurence
#For each threshold, decide a duration and severity.
#Compute the output probability for data and simulations.




setwd("~/ERCOT") #Code for personal device


#______________________________________________________________________________#
#.libPaths("/rigel/cwc/users/yva2000/rpackages/")
###Load Packages       
library(dplyr)
library(ggplot2)
library(logspline)
library(foreach)    #Parallel Execution
library(doParallel) #Backend to foreach
library(locfit)
library(lubridate)
library(gridExtra)

###Load Functions
source("functions/Get_Severity_Duration.R")

#______________________________________________________________________________#
###Reading the Data and Cleaning it
grid_locs <- read.csv("data/ERCOT_0_5_deg_lat_lon_index_key.csv", 
                      header = TRUE, sep=",")

ssrd <- read.table("data/ERCOT_Solar_Rad_Daily.txt",  
                   sep =" ", 
                   header = TRUE) #Surface Solar Radiation Downwards

WP <- read.table("data/ERCOT_Wind_Power_Daily.txt",  
                 sep =" ", 
                 header = TRUE) #Surface Solar Radiation Downwards

Fld <- as.matrix(cbind(WP,ssrd))
Fld <- apply(Fld, 2, function(x){x/max(x)})

ssrd <- WP <- NULL

#______________________________________________________________________________#
###Read the simulations
load("Simulations/Joint_Raw_Simulations.RData")
nl <- 6 #length(ynew_results)
comb_ksts <- list()
for(i in 1:nl){
  comb_ksts[[i]] <- cbind(ynew_results[[i]]$WPnew, ynew_results[[i]]$SSnew)
  comb_ksts[[i]] <- apply(comb_ksts[[i]], 2, function(x){x/max(x)})
  ynew_results[[i]] <- 1
}
ynew_results <- NULL


#______________________________________________________________________________#
###Function to fit a log-fit to the output of severity duration.
###Inputs
#1. Data
#2. Simulations
#3. Threshold Percentile
#4. Severity (2)
#5. Duration (3)

###Outputs
#1. Panel Plot

#Dat <- Fld
#Sims <- comb_ksts
#Thresh <- 0.25
#Severities <- c(1.25,1.5)
#Durations <- c(20,25,30)

get_aep_plot <- function(Dat, Sims, Thresh, Severities, Durations){
  
  #Hyper-Parameters
  thresh = Thresh
  st_date = "01-01-1979"
  mean_daily <- mean(rowSums(Dat))
  nl <- length(comb_ksts)
  
  #Set-up Grid of Values
  Dur = Durations
  Sev = Severities
  ngrids <- length(Sev)*length(Dur)
  Sev_Dur <- data.frame(Duration = rep(Dur, each = length(Sev)),
                        Severity = rep(Sev, length(Dur)))
  
  
  #---------------------------------------------------------------------------#
  ###Function to fit a log-fit to the output of severity duration.
  ###Inputs of Function
  #1. Severity-Duration for all Data.
  #2. Desired Duration
  #3. Desired Severity
  #4. Mean Daily Value to facilitate comparision
  #5. Number of Years
  
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
                                          MD = mean_daily, N_yr = 40)
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
                                            MD = mean_daily, N_yr = 40)
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
  
  #pdf(paste0("Energy_Droughts_",100*thresh,"_pt_SP.pdf"))
  
  p <- ggplot(sims_plot, aes(x=Duration, y=Exceedance)) + 
    geom_boxplot() +
    ggtitle(paste0("Energy Droughts & Annual Exceedances \n Threshold - ", thresh*100,"th Percentile")) +
    ylab("Annual Exceedance Probabilities (%)") +
    xlab("Duration (Days)") +
    theme_bw() +
    theme(axis.text=element_text(size=20),
          axis.title=element_text(size=20),
          plot.title = element_text(size=20),
          strip.text = element_text(size=20))
  
  p <- p + facet_grid(Severity ~ ., scale = "free")
  p <- p + geom_point(data_plot, mapping = aes(x=Duration, y=Exceedance), color = 'red', size = 2)
  print(p)
  
  
}



#______________________________________________________________________________#
###Function to fit a log-fit to the output of severity duration.
get_aep_plot(Dat = Fld, Sims = comb_ksts, Thresh = 0.21, 
             Severities = c(1.25,1.5), Durations = c(20,25,30))



get_aep_plot(Dat = Fld, Sims = comb_ksts, Thresh = 0.22, 
             Severities = c(1.25,1.5), Durations = c(20,25,30))






















learn_git <- function(){




#______________________________________________________________________________#
###Function to fit a log-fit to the output of severity duration.
###Inputs of Function
#1. Severity-Duration for all Data.
#2. Desired Duration
#3. Desired Severity
#4. Mean Daily Value to facilitate comparision
#5. Number of Years

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


#______________________________________________________________________________#
###Hyper-Parameters
True_Data = Fld
thresh = 0.25
st_date = "01-01-1979"

#Set-up Grid of Values
Dur <- c(20,25,30)
Sev = c(1.25,1.5)
ngrids <- length(Sev)*length(Dur)

Sev_Dur <- data.frame(Duration = rep(Dur, each = length(Sev)),
                         Severity = rep(Sev, length(Dur)))

#______________________________________________________________________________#
###Compute the severity and drought for the data.

#Source the function.
source("functions/Get_Severity_Duration.R")

#Hyper-Parameters
mean_daily <- mean(rowSums(True_Data))

#Compute Severity and Duration for Data
Data_SD <- get_Severity_Duration(Data = True_Data, Thresh = thresh,
                                 start_date = st_date, Type = "Data")

data_pexc <- list()
for(i in 1:nrow(Sev_Dur)){
  data_pexc[[i]] <- 100*get_annual_prob(ED = Data_SD, 
                                        Dur = Sev_Dur$Duration[i],
                                        Sev = Sev_Dur$Severity[i], 
                                        MD = mean_daily, N_yr = 40)
}
data_pexc <- unlist(data_pexc)


#______________________________________________________________________________#
###Compute the severity and drought for the the simulations.
#sim_pexc <- list()
#for(j in 1:nl){
#  
#  #Compute Severity and Duration for Simulations
#  KSTS_SD <- get_Severity_Duration(Data = comb_ksts[[j]], Thresh = thresh,
#                                   start_date = st_date, Type = "Sims")
#  
#  
#  
#  temp_pexc <- list()
#  for(i in 1:nrow(Sev_Dur)){
#    temp_pexc[[i]] <- 100*get_annual_prob(ED = KSTS_SD, 
#                                          Dur = Sev_Dur$Duration[i],
#                                          Sev = Sev_Dur$Severity[i], 
#                                          MD = mean_daily, N_yr = 40)
#  }
#  
#  sim_pexc[[j]] <- unlist(temp_pexc)
#}


#______________________________________________________________________________#
###Compute the severity and drought for the the simulations.
cores=detectCores()
registerDoParallel(cores)
sim_pexc <- list()
sim_pexc <- foreach(m = 1:nl, .verbose = TRUE) %dopar% {
  
  #Compute Severity and Duration for Simulations
  KSTS_SD <- get_Severity_Duration(Data = comb_ksts[[m]], Thresh = thresh,
                                   start_date = st_date, Type = "Sims")
  
  #Function to get the annual exceedances
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
                                          MD = mean_daily, N_yr = 40)
  }
  
  sim_pexc[[m]] <- unlist(temp_pexc)
}
stopImplicitCluster()

#______________________________________________________________________________#
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

#pdf(paste0("Energy_Droughts_",100*thresh,"_pt_SP.pdf"))

p <- ggplot(sims_plot, aes(x=Duration, y=Exceedance)) + 
  geom_boxplot() +
  ggtitle(paste0("Energy Droughts & Annual Exceedances \n Threshold - ", thresh*100,"th Percentile")) +
  ylab("Annual Exceedance Probabilities (%)") +
  xlab("Duration (Days)") +
  theme_bw() +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20),
        plot.title = element_text(size=20),
        strip.text = element_text(size=20))

p <- p + facet_grid(Severity ~ ., scale = "free")
p <- p + geom_point(data_plot, mapping = aes(x=Duration, y=Exceedance), color = 'red', size = 2)
print(p)

#dev.off()

}
