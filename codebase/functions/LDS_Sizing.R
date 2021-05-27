#Relationship between sizing of Long Duration Systems(LDS) and Data Record Length
#The LDS size is the largest severity for a particular thershold. 

#Step 1:- Transfrom the data to have max of 1. 
#Step 2:- Compute the maximum severity given a threshold. 
#Step 3:- Divide the Severity by mean daily production.


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
#Hyper-Parametres and source functions
source("functions/Get_Severity_Duration.R")
thresh <- 0.3
data_sizing <-  list()  


#______________________________________________________________________________#
#10-yr break points
start_dates <- c("01-01-1979","01-01-1989","01-01-1999","01-01-2009")
st_date <- as.Date("01-01-1979", format="%m-%d-%Y")
time_stamps <- seq(st_date, by = "day", length.out = nrow(WP))

#First Decade
st <- which(time_stamps == as.Date(start_dates[[1]],format="%m-%d-%Y"))
ed <- which(time_stamps == as.Date(start_dates[[2]],format="%m-%d-%Y"))
Fld <- as.matrix(cbind(WP[(st:(ed-1)),],ssrd[(st:(ed-1)),]))
Fld <- apply(Fld, 2, function(x){x/max(x)})
sev <- get_Severity_Duration(Data = Fld, Thresh = thresh, 
                             start_date = "01-01-1979", Type = "10-yrs")
sev <- sev[which(sev$Severity == max(sev$Severity)),]
sev$Severity <- sev$Severity/mean(rowSums(Fld))
data_sizing[[1]] <- sev

#Second Decade
st <- which(time_stamps == as.Date(start_dates[[2]],format="%m-%d-%Y"))
ed <- which(time_stamps == as.Date(start_dates[[3]],format="%m-%d-%Y"))
Fld <- as.matrix(cbind(WP[(st:(ed-1)),],ssrd[(st:(ed-1)),]))
Fld <- apply(Fld, 2, function(x){x/max(x)})
sev <- get_Severity_Duration(Data = Fld, Thresh = thresh, 
                             start_date = "01-01-1989", Type = "10-yrs")
sev <- sev[which(sev$Severity == max(sev$Severity)),]
sev$Severity <- sev$Severity/mean(rowSums(Fld))
data_sizing[[2]] <- sev


#Third Decade
st <- which(time_stamps == as.Date(start_dates[[3]],format="%m-%d-%Y"))
ed <- which(time_stamps == as.Date(start_dates[[4]],format="%m-%d-%Y"))
Fld <- as.matrix(cbind(WP[(st:(ed-1)),],ssrd[(st:(ed-1)),]))
Fld <- apply(Fld, 2, function(x){x/max(x)})
sev <- get_Severity_Duration(Data = Fld, Thresh = thresh, 
                             start_date = "01-01-1999", Type = "10-yrs")
sev <- sev[which(sev$Severity == max(sev$Severity)),]
sev$Severity <- sev$Severity/mean(rowSums(Fld))
data_sizing[[3]] <- sev


#Last Decade
st <- which(time_stamps == as.Date(start_dates[[4]],format="%m-%d-%Y"))
Fld <- as.matrix(cbind(WP[(st:nrow(WP)),],ssrd[(st:nrow(ssrd)),]))
Fld <- apply(Fld, 2, function(x){x/max(x)})
sev <- get_Severity_Duration(Data = Fld, Thresh = thresh, 
                             start_date = "01-01-2009", Type = "10-yrs")
sev <- sev[which(sev$Severity == max(sev$Severity)),]
sev$Severity <- sev$Severity/mean(rowSums(Fld))
data_sizing[[4]] <- sev

#Bootstrap 10-yr data to avoid start date bias.
n_length <- 365*10
ns <- 48
nsample <- sample(1:(nrow(WP)-n_length-1),ns, replace = TRUE)
data_10_boot <- list()

cores=detectCores()
registerDoParallel(cores)
data_10_boot <- foreach(m = 1:ns, .verbose = TRUE) %dopar% {
  source("functions/Get_Severity_Duration.R")
  st <- nsample[m]
  ed <- n_length + st
  Fld <- as.matrix(cbind(WP[st:ed,],ssrd[st:ed,]))
  Fld <- apply(Fld, 2, function(x){x/max(x)})
  sev <- get_Severity_Duration(Data = Fld, Thresh = thresh, 
                               start_date = time_stamps[nsample[m]],
                               Type = "10-yrs")
  sev <- sev[which(sev$Severity == max(sev$Severity)),]
  sev$Severity <- sev$Severity/mean(rowSums(Fld))
  data_10_boot[[m]] <- sev
  }
stopImplicitCluster()


#______________________________________________________________________________#
#20-yr break points
start_dates <- c("01-01-1979","01-01-1999")
st_date <- as.Date("01-01-1979", format="%m-%d-%Y")
time_stamps <- seq(st_date, by = "day", length.out = nrow(WP))

#First Double Decade
st <- which(time_stamps == as.Date(start_dates[[1]],format="%m-%d-%Y"))
ed <- which(time_stamps == as.Date(start_dates[[2]],format="%m-%d-%Y"))
Fld <- as.matrix(cbind(WP[(st:(ed-1)),],ssrd[(st:(ed-1)),]))
Fld <- apply(Fld, 2, function(x){x/max(x)})
sev <- get_Severity_Duration(Data = Fld, Thresh = thresh, 
                             start_date = "01-01-1979", Type = "20-yrs")
sev <- sev[which(sev$Severity == max(sev$Severity)),]
sev$Severity <- sev$Severity/mean(rowSums(Fld))
data_sizing[[5]] <- sev

#Last Double Decade
st <- which(time_stamps == as.Date(start_dates[[2]],format="%m-%d-%Y"))
Fld <- as.matrix(cbind(WP[(st:nrow(WP)),],ssrd[(st:nrow(ssrd)),]))
Fld <- apply(Fld, 2, function(x){x/max(x)})
sev <- get_Severity_Duration(Data = Fld, Thresh = thresh, 
                             start_date = "01-01-1999", Type = "20-yrs")
sev <- sev[which(sev$Severity == max(sev$Severity)),]
sev$Severity <- sev$Severity/mean(rowSums(Fld))
data_sizing[[6]] <- sev


#Bootstrap 20-yr data
n_length <- 365*20
ns <- 48
nsample <- sample(1:(nrow(WP)-n_length-1),ns, replace = TRUE)
data_20_boot <- list()

cores=detectCores()
registerDoParallel(cores)
data_20_boot <- foreach(m = 1:ns, .verbose = TRUE) %dopar% {
  source("functions/Get_Severity_Duration.R")
  st <- nsample[m]
  ed <- n_length + st
  Fld <- as.matrix(cbind(WP[st:ed,],ssrd[st:ed,]))
  Fld <- apply(Fld, 2, function(x){x/max(x)})
  sev <- get_Severity_Duration(Data = Fld, Thresh = thresh, 
                               start_date = time_stamps[nsample[m]],
                               Type = "20-yrs")
  sev <- sev[which(sev$Severity == max(sev$Severity)),]
  sev$Severity <- sev$Severity/mean(rowSums(Fld))
  data_20_boot[[m]] <- sev
}
stopImplicitCluster()


#______________________________________________________________________________#
#The entire dataset
Fld <- as.matrix(cbind(WP,ssrd))
Fld <- apply(Fld, 2, function(x){x/max(x)})
sev <- get_Severity_Duration(Data = Fld, Thresh = thresh, 
                             start_date = "01-01-1979", Type = "40-yrs")
sev <- sev[which(sev$Severity == max(sev$Severity)),]
sev$Severity <- sev$Severity/mean(rowSums(Fld))
data_sizing[[7]] <- sev


#______________________________________________________________________________#
#______________________________________________________________________________#
#______________________________________________________________________________#
###Simulations###





#______________________________________________________________________________#
###10-yr Simulations###

#Load the simulations
load("Simulations/Joint_Raw_Simulations_10yrs.RData")
ns <- length(ynew_results)
sim <- list()
for(i in 1:ns){
  sim[[i]] <- as.matrix(cbind(ynew_results[[i]]$WPnew,
                         ynew_results[[i]]$SSnew))
  sim[[i]] <- apply(sim[[i]], 2, function(x){x/max(x)})
  ynew_results[[i]] <- 1 
}
ynew_results <- NULL


#Compute the drought severity
sim_10_yrs <- list()
cores=detectCores()
registerDoParallel(cores)
start.time <- Sys.time()
sim_10_yrs <- foreach(m = 1:ns, .verbose = TRUE) %dopar% {
  source("functions/Get_Severity_Duration.R")
  
  Sim_Fld <- sim[[m]]
  Sim_Fld <- apply(Sim_Fld, 2, function(x){x/max(x)})
  sev <- get_Severity_Duration(Data = Sim_Fld, Thresh = thresh, 
                               start_date = "01-01-1979",
                               Type = "10-yrs")
  sev <- sev[which(sev$Severity == max(sev$Severity)),]
  sev$Severity <- sev$Severity/mean(rowSums(Sim_Fld))
  sim_10_yrs[[m]] <- sev
}
end.time <- Sys.time()
time.taken <- end.time - start.time
print(time.taken)
stopImplicitCluster()
sim <- NULL


#______________________________________________________________________________#
###20-yr Simulations###

sim <- list()

#Load the simulations - Part 1
load("Simulations/Joint_Raw_Simulations_20yrs_1.RData")
ns <- length(ynew_results)
for(i in 1:ns){
  sim[[i]] <- as.matrix(cbind(ynew_results[[i]]$WPnew,
                              ynew_results[[i]]$SSnew))
  sim[[i]] <- apply(sim[[i]], 2, function(x){x/max(x)})
  ynew_results[[i]] <- 1 
}
ynew_results <- NULL

#Load the simulations - Part 2
load("Simulations/Joint_Raw_Simulations_20yrs_2.RData")
ns2 <- length(ynew_results)
for(i in 1:ns){
  sim[[ns+i]] <- as.matrix(cbind(ynew_results[[i]]$WPnew,
                              ynew_results[[i]]$SSnew))
  sim[[ns+i]] <- apply(sim[[ns+i]], 2, function(x){x/max(x)})
  ynew_results[[i]] <- 1 
}
ynew_results <- NULL


#Compute the drought severity
sim_20_yrs <- list()
cores=detectCores()
registerDoParallel(cores)
start.time <- Sys.time()
ns <- length(sim)
sim_20_yrs <- foreach(m = 1:ns, .verbose = TRUE) %dopar% {
  source("functions/Get_Severity_Duration.R")
  
  Sim_Fld <- sim[[m]]
  Sim_Fld <- apply(Sim_Fld, 2, function(x){x/max(x)})
  sev <- get_Severity_Duration(Data = Sim_Fld, Thresh = thresh, 
                               start_date = "01-01-1979",
                               Type = "20-yrs")
  sev <- sev[which(sev$Severity == max(sev$Severity)),]
  sev$Severity <- sev$Severity/mean(rowSums(Sim_Fld))
  sim_20_yrs[[m]] <- sev
}
end.time <- Sys.time()
time.taken <- end.time - start.time
print(time.taken)
stopImplicitCluster()
sim <- NULL


#______________________________________________________________________________#
###40-yr Simulations###

#Load the simulations
load("Simulations/Joint_Raw_Simulations.RData")
ns <- length(ynew_results)
sim <- list()
for(i in 1:ns){
  sim[[i]] <- as.matrix(cbind(ynew_results[[i]]$WPnew,
                              ynew_results[[i]]$SSnew))
  ynew_results[[i]] <- 1 
  sim[[i]] <- apply(sim[[i]], 2, function(x){x/max(x)})
  
}
ynew_results <- NULL


#Compute the drought severity
sim_40_yrs <- list()
cores=detectCores()
registerDoParallel(cores)
start.time <- Sys.time()
sim_40_yrs <- foreach(m = 1:ns, .verbose = TRUE) %dopar% {
  source("functions/Get_Severity_Duration.R")
  
  Sim_Fld <- sim[[m]]
  Sim_Fld <- apply(Sim_Fld, 2, function(x){x/max(x)})
  sev <- get_Severity_Duration(Data = Sim_Fld, Thresh = thresh, 
                               start_date = "01-01-1979",
                               Type = "40-yrs")
  sev <- sev[which(sev$Severity == max(sev$Severity)),]
  sev$Severity <- sev$Severity/mean(rowSums(Sim_Fld))
  sim_40_yrs[[m]] <- sev
}
end.time <- Sys.time()
time.taken <- end.time - start.time
print(time.taken)
stopImplicitCluster()
sim <- NULL



#______________________________________________________________________________#
#Convert to a Data-Frame Structure
data_sizing <- bind_rows(lapply(data_sizing,data.frame))
data_10_boot <- bind_rows(lapply(data_10_boot,data.frame))
data_20_boot <- bind_rows(lapply(data_20_boot,data.frame))
sim_10_yrs <- bind_rows(lapply(sim_10_yrs,data.frame))
sim_20_yrs <- bind_rows(lapply(sim_20_yrs,data.frame))
sim_40_yrs <- bind_rows(lapply(sim_40_yrs,data.frame))



###Plot comparing data and simulations.
LDS_Sizing_Sims <- rbind(sim_10_yrs, sim_20_yrs, sim_40_yrs, data_sizing[7,])
ggplot(LDS_Sizing_Sims, aes(x=Type, y=Severity)) + 
  geom_violin() +
  scale_x_discrete(name ="Data/Record Length") +
  scale_y_continuous(name ="LDS Size (fraction of mean production)") +
  geom_point(data_sizing, mapping = aes(x=Type, y=Severity), 
             size = 5, fill = 'red', shape = 23) +
  ggtitle(paste0("Sizing - LDS Systems \n ", thresh*100,"th Percentile")) +
  theme_bw() +
  theme(legend.text=element_text(size=20),
        legend.title=element_text(size=15),
        axis.text=element_text(size=15),
        axis.title=element_text(size=20),
        plot.title = element_text(size=20))


#Plot comparing data bootstrap and simulations
LDS_Sims <- rbind(sim_10_yrs, sim_20_yrs, sim_40_yrs)
LDS_Sims$Class <- "Simulations"
LDS_Data <- rbind(data_10_boot,data_20_boot, data_sizing[7,])
LDS_Data$Class <- "Data"
LDS_Sizing <- rbind(LDS_Sims, LDS_Data)

ggplot(LDS_Sizing, aes(x=Type, y=Severity, fill = Class)) + 
  geom_violin() +
  scale_x_discrete(name ="Data/Record Length") +
  scale_y_continuous(name ="Severity (fraction of mean production)") +
  ggtitle(paste0("Sizing - LDS Systems \n ", thresh*100,"th Percentile")) +
  theme_bw() +
  theme(legend.text=element_text(size=20),
        legend.title=element_text(size=15),
        axis.text=element_text(size=15),
        axis.title=element_text(size=20),
        plot.title = element_text(size=20))

