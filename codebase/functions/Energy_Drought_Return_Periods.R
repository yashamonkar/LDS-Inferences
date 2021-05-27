#______________________________________________________________________________#
#Code to get the return periods for Severity and Duartion. 
# P(d > D, s >S)


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


#Set-up Data
Fld <- as.matrix(cbind(WP,ssrd))
Fld <- apply(Fld, 2, function(x){x/max(x)})
mean_daily <- mean(rowSums(Fld))

#______________________________________________________________________________#
###Read the simulations
#KSTS
comb_ksts <- list()
load("Simulations/Joint_Raw_Simulations.RData")
nl <-  length(ynew_results)
for(i in 1:nl){
  comb_ksts[[i]] <- cbind(ynew_results[[i]]$WPnew, ynew_results[[i]]$SSnew)
  ynew_results[[i]] <- 1
  comb_ksts[[i]] <- apply(comb_ksts[[i]], 2, function(x){x/max(x)})
}
ynew_results <- NULL


#______________________________________________________________________________#
###Compute the Severity and Duration for particular drought events
source("functions/Get_Severity_Duration.R")
thresh <- 0.3

#Compute Sev-Dur
cores=detectCores()
registerDoParallel(cores)
Sev_Dur <- foreach(m = 1:nl, .verbose = TRUE) %dopar% {
  source("functions/Get_Severity_Duration.R")
  get_Severity_Duration(Data = comb_ksts[[m]], Thresh = thresh,
                        start_date = "01-01-1979", Type = "Sims")
}
stopImplicitCluster()


#______________________________________________________________________________#
#Stack up all results
Sev_Dur <- bind_rows(lapply(Sev_Dur,data.frame))
Sev_Dur <- Sev_Dur %>% filter(Duration > 3)
Sev_Dur$Severity <- Sev_Dur$Severity/mean_daily

#Events per Duration 
dmin <- min(Sev_Dur$Duration)
dmax <- max(Sev_Dur$Duration)
dur_events <- data.frame(Duration = dmin:dmax,
                            Tally = dmin:dmax)
for(i in 1:nrow(dur_events)){
  t <- Sev_Dur %>% filter(Duration == dur_events$Duration[i])
  dur_events$Tally[i] <- nrow(t)
}
dur_events$Exceedance <- rev(cumsum(rev(dur_events$Tally)))

#Compute the probability of exceedances per-year
dur_events$Probs <-  dur_events$Exceedance/(nl*40)

#Plot for the marginal distribution
temp_events <- dur_events %>% filter(Tally > 0) 
ggplot(temp_events, aes(x = 1/Probs, y = Duration/7)) +
  ggtitle(paste0("Return Period for Energy Drought Duration - ", 
                 thresh*100,"th Percentile")) +
  geom_line(size = 1.5) +
  xlab("Return Period (Years)") +
  ylab("Duration (Weeks)") +
  theme_bw() +
  theme(legend.text=element_text(size=20),
        legend.title=element_text(size=15),
        axis.text=element_text(size=15),
        axis.title=element_text(size=20),
        plot.title = element_text(size=20))

#______________________________________________________________________________#
###Compute the joint return periods. 
#P(Sev,Dur) = P(Sev|Dur)P(Dur)
#P(Dur) has been computed in the steps above.
#Target contours of 0.05, 0.001, 0.0005


###Contour of 0.05/Return Period of 20 years
cont_level <- 0.05
events <- sort(unique(Sev_Dur$Duration))

rperiod <- list()
for(i in 1:length(events)){
  temp_dur <- Sev_Dur %>% filter(Duration == events[i])
  prob <- cont_level/(nrow(temp_dur)/(nl*40))
  if(prob < 0.75) {
    rperiod[[i]] <- quantile(temp_dur$Severity, probs = (1-prob))
    
  } else {
    rperiod[[i]] <- NA
  }
}
p <- cbind(Duration = events,
             Severity = unlist(rperiod))
p_20 <- as.data.frame(p[complete.cases(p), ])


###Contour of 0.02/Return Period of 50 years
cont_level <- 0.02
events <- sort(unique(Sev_Dur$Duration))

rperiod <- list()
for(i in 1:length(events)){
  temp_dur <- Sev_Dur %>% filter(Duration == events[i])
  prob <- cont_level/(nrow(temp_dur)/(nl*40))
  if(prob < 0.75) {
    rperiod[[i]] <- quantile(temp_dur$Severity, probs = (1-prob))
    
  } else {
    rperiod[[i]] <- NA
  }
}
p <- cbind(Duration = events,
           Severity = unlist(rperiod))
p_50 <- as.data.frame(p[complete.cases(p), ])


###Contour of 0.01/Return Period of 100 years
cont_level <- 0.01
events <- sort(unique(Sev_Dur$Duration))

rperiod <- list()
for(i in 1:length(events)){
  temp_dur <- Sev_Dur %>% filter(Duration == events[i])
  prob <- cont_level/(nrow(temp_dur)/(nl*40))
  if(prob < 0.75) {
    rperiod[[i]] <- quantile(temp_dur$Severity, probs = (1-prob))
    
  } else {
    rperiod[[i]] <- NA
  }
}
p <- cbind(Duration = events,
           Severity = unlist(rperiod))
p_100 <- as.data.frame(p[complete.cases(p), ])

###Contour of 0.005/Return Period of 200 years
cont_level <- 0.005
events <- sort(unique(Sev_Dur$Duration))

rperiod <- list()
for(i in 1:length(events)){
  temp_dur <- Sev_Dur %>% filter(Duration == events[i])
  prob <- cont_level/(nrow(temp_dur)/(nl*40))
  if(prob < 0.75) {
    rperiod[[i]] <- quantile(temp_dur$Severity, probs = (1-prob))
    
  } else {
    rperiod[[i]] <- NA
  }
}
p <- cbind(Duration = events,
           Severity = unlist(rperiod))
p_200 <- as.data.frame(p[complete.cases(p), ])


###Contour of 0.002/Return Period of 500 years
cont_level <- 0.002
events <- sort(unique(Sev_Dur$Duration))

rperiod <- list()
for(i in 1:length(events)){
  temp_dur <- Sev_Dur %>% filter(Duration == events[i])
  prob <- cont_level/(nrow(temp_dur)/(nl*40))
  if(prob < 0.75) {
    rperiod[[i]] <- quantile(temp_dur$Severity, probs = (1-prob))
    
  } else {
    rperiod[[i]] <- NA
  }
}
p <- cbind(Duration = events,
           Severity = unlist(rperiod))
p_500 <- as.data.frame(p[complete.cases(p), ])


###Contour of 0.001/Return Period of 1000 years
cont_level <- 0.001
events <- sort(unique(Sev_Dur$Duration))

rperiod <- list()
for(i in 1:length(events)){
  temp_dur <- Sev_Dur %>% filter(Duration == events[i])
  prob <- cont_level/(nrow(temp_dur)/(nl*40))
  if(prob < 0.75) {
    rperiod[[i]] <- quantile(temp_dur$Severity, probs = (1-prob))
    
  } else {
    rperiod[[i]] <- NA
  }
}
p <- cbind(Duration = events,
           Severity = unlist(rperiod))
p_1000 <- as.data.frame(p[complete.cases(p), ])






ggplot() +
  geom_point(data = Sev_Dur, mapping = aes(x = log(Duration), y = log(Severity)),
             alpha = 0.15) +
  geom_line(data = p_20, mapping = aes(x = log(Duration), y = log(Severity)), 
            col ='red', size = 1.5) +
  geom_line(data = p_50, mapping = aes(x = log(Duration), y = log(Severity)), 
            col ='blue', size = 1.5) +
  geom_line(data = p_100, mapping = aes(x = log(Duration), y = log(Severity)), 
            col ='green', size = 1.5) +
  geom_line(data = p_200, mapping = aes(x = log(Duration), y = log(Severity)), 
            col ='orange', size = 1.5) +
  geom_line(data = p_500, mapping = aes(x = log(Duration), y = log(Severity)), 
            col ='violet', size = 1.5) +
  geom_line(data = p_1000, mapping = aes(x = log(Duration), y = log(Severity)), 
            col ='yellow', size = 1.5) +
  ggtitle("Joint Severity and Duration Return Periods")



