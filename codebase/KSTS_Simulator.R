#______________________________________________________________________________#
###K-Nearest Neighbor Space-Time Simulator (KSTS)###


###Data Inputs
#1. Wind Data
#2. Solar Data
#3. Grid Location

###Code Inputs
#1. Functions for Simulation Checks

###Outputs
#1. KSTS Simulations
#2. Simulation Checks

 
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
###Set the KSTS Hyper-Parameters###

###Subset
#n <- 5*365
#WP <- WP[1:n,]
#ssrd <- ssrd[1:n,]

###Concatenate the fields
Fld <- as.matrix(cbind(WP,ssrd))
colnames(Fld) <- NULL

###Data Parameters
n_site <- ncol(ssrd)        #Number of Sites 
ngrids <- ncol(Fld)         #Number of grids (Grids = Sites x Fields)
N_valid <- nrow(Fld)        #Number of Time Steps
 
###Embeddings/State Space Parameters
max_embd <- 2  #Max Value of Lagged Embedding
sel_lags <- c(1,2) #Individual Lags Selected
n_lags <- length(sel_lags) #Number of Lags
w <- c(1,0) #Scaling Weights


#______________________________________________________________________________#
###KSTS Algorithm Functions###

###Function 1 
#Objective - Moving Window Indices - Returns the day/hour index of interest.
###Inputs
#   1. Current Day-of-Year (curr)
#   2. Max value of Day-of-Year (max)
#   3. Moving Window Size (window)
###Output
#   1. All DOYs in the moving window

close_ind <- function(curr,max,window){
  if((curr-window) < 1){
    
    indx <- c(tail(1:max,(1+abs(curr-window))),
              1:(curr+window))
    
  } else if((curr+window) > max) {
    
    indx <- c(((curr-window):max),
              (1:((curr+window)-max)))
    
  } else {
    
    indx <- (curr-window):(curr+window)
    
  }
  return(indx)
}



###Function 2 
#Objective - Compute the K-Nearest Neighbors for each site
###Inputs
#   1. Current Feature Vector (x)
#   2. All historic Feature Vectors (xtest)
#   3. Number of Neighbors (nneib)
#   4. Scaling Weights (weights)
###Outputs
#   1. Indices corresponding to the k-nearest neighbors 

knn.sim.index <- function(x,xtest,nneib,weights){
  
  d <- matrix(NA, ncol=ncol(xtest), nrow=nrow(x))
  #Compute Distances
  for(i in 1:ncol(x))
    d[,i] <- weights[i]*(x[,i]-xtest[,i])^2
  d <- rowSums(d)
  sorted.data = sort.int(d,method="quick",index.return=TRUE)
  
  sorted.neib = matrix(sorted.data$ix[1:nneib])
  
  out = list(yknn=sorted.neib
  )
}

###Function 3 
#Objective - KSTS Simulator

###Input
#   1. Concatenated Data (Fld)
#   2. Number of Grid Points (ngrids)
#   3. Record Length (N_valid)
#   4. Number of Nearest Neighbors (nneib)
#   5. Scaling Weights (weights)
#   6. Record Start Date (start_date)
#   7. Moving Window Size (day_mv)
#   8. Maximum Embedding (max_embd)
#   9. Selected Embedding Lags (sel_lags)
#   10. Number of selected lags  (n_lags)

###Output
#   1. A single KSTS Simulation Realization


#Knn with embeddings without climate 
ksts <- function(Fld,ngrids, N_valid, #Data Parameters
                            nneib,weights, #KNN Parameter
                            start_date, day_mv,  #Date Seasonality Parameters
                            max_embd, sel_lags, n_lags) #Embedding Parameters
{ 
  #Load Dependcies
  library(dplyr)
  library(lubridate)
  
  #Day Indices
  st_date <- as.Date(start_date, format="%m-%d-%Y")
  time_stamps <- seq(st_date, by = "day", length.out = N_valid)
  day_index <- yday(time_stamps)

  #Setting up Storage for Simulations
  Xnew = matrix(NA,nr=N_valid,nc=ngrids)
  for(i in 1:max_embd){
    Xnew[i,] = jitter(Fld[i,]) 
  }
  
  #Creating the feature Vector/state space
  X = array(NA, c(N_valid-max_embd,n_lags,ngrids))
  Y = array(NA, c(N_valid-max_embd,1,ngrids))
  for (j in 1:ngrids) 
  {
    #Get Lagged Structure Upto Max Embedding
    x_fld <- embed(Fld[,j],max_embd+1) 
    
    X[,,j] <- x_fld[,sel_lags+1]
    Y[,,j] <- x_fld[,1]
  }
  
  
  #Starting the Simulator
  pb = txtProgressBar(min = (max_embd+1), max = N_valid, initial = 1) 
  for (i in (max_embd+1):N_valid){
    
    setTxtProgressBar(pb,i)
    
    nn_index <- list() #Store all the nearest neighbours.
    day <- day_index[i]
    sel_days <- close_ind(day,max = 366, window = day_mv)
    
    #Subset to the moving window
    indx <- day_index
    indx[i] <- 999  #Remove the current value
    days <- tail(indx, -max_embd)
    X_t <- X[days %in% sel_days,,]
    Y_t <- Y[days %in% sel_days,,]
    
    for (j in 1:ngrids){
      
      #Setting the Test Parameters
      sel_pars <- i-sel_lags
      xtest <- matrix(Xnew[sel_pars,j], nrow =1)
      
      #Running the KNN Algorithm
      nn_index[[j]] <- knn.sim.index(X_t[,,j],xtest,nneib,weights)
    }
    
    #Computing the Resampling Probability
    un_index <- unique(unlist(nn_index))
    un_prob <- rep(NA, length(un_index))
    
    nn_index  <-  matrix(unlist(nn_index), nrow=nneib)
    
    for(k in 1:length(un_index)){
      temp <- which(nn_index == un_index[k]) %% nneib
      for(l in 1:length(temp)){
        if(temp[l]==0) {
          temp[l] = 1/nneib
        } else {
        temp[l]=1/temp[l]
        }
      }
      un_prob[k] <- sum(temp)
    }
  #Computing the Resampling Probablites 
  pj <- data.frame(cbind(un_prob,un_index)) 
  thresh <- apply(pj, 2, FUN = function(x) tail(sort(x), nneib+1)[1])[1]
  pj <- pj %>% filter(un_prob > thresh)  
  pj$un_prob <- pj$un_prob/sum(pj$un_prob)
  ns <- sample(pj$un_index,1,prob = pj$un_prob)
  
  Xnew[i,] <- Y_t[ns,]
  
  }
  n_site <- ngrids/2
  WPnew = Xnew[,1:n_site]
  SSnew = Xnew[,(n_site+1):ncol(Xnew)]
  
  out = list(WPnew=WPnew,SSnew=SSnew)
}


#______________________________________________________________________________#
###Load the Simulation Skill Functions###

#General Simulation Skill Checks
source("functions/Get_Simulation_Skill.R")

#PC-Wavelet Analysis
source("functions/Get_PCWavelet_Sim_Skill.R")

#Annual-Monthly Cycle
source("functions/Get_Annual_Cycle.R")

#Wind-Solar Cross Correlation
source("functions/Get_Site_Correlation.R")

#PC Eigenvector Analysis
source("functions/Get_PCA_ggplot.R")

#Wind-Solar Seasonal Correlation
source('functions/Get_Seasonal_Correlation.R')



###Simulation Hyper-Parameters###
nneib <- 65
nsim <- 48


#______________________________________________________________________________#
###Run the Simulation###

#Multi-Core Processing
cores=detectCores()
registerDoParallel(cores)
start.time <- Sys.time()
ynew_results <- foreach(m = 1:nsim, .verbose = TRUE) %dopar% {
    ksts(Fld,ngrids, N_valid,
         nneib,w,
         start_date = "01-01-1950",day_mv = 30,
         max_embd, sel_lags, n_lags)
  }
end.time <- Sys.time()
time.taken <- end.time - start.time
print(time.taken)
stopImplicitCluster()

#Saving Raw Simulations
save(ynew_results, file = paste0("KSTS_Joint_Simulations.RData") )


#______________________________________________________________________________#
###Plotting the Simulation Skill Metrics###

#Saving Simulation Skill
pdf(paste0("KSTS_Daily_Joint_",nneib,".pdf"))

###Seperating the Simulations
wind_sims <- solar_sims <- list()
for(i in 1:nsim){
  wind_sims[[i]] <- ynew_results[[i]]$WPnew
  solar_sims[[i]] <- ynew_results[[i]]$SSnew
}
ynew_results <- NULL


###Wind
Field <- "Wind"
fld <- Fld[,1:n_site]
Get_Simulation_Skill(True_Data = fld,
                     Simulations = wind_sims,
                     Field_Name = Field)

Get_Annual_Cycle(True_Data = fld, Field_Name = Field,
                 Simulations = wind_sims,
                 Resolution = "daily",
                 Start_Date = "01-01-1950")
  
get_pca_plot(X = fld, 
               Grid = grid_locs, 
               Field = Field,
               Sims = wind_sims)
  
Get_PCWavelet_Sim_Skill(Data_Field = fld,
                          Field_Name = Field,
                          Sims = wind_sims,
                          PCs=2)
  
  
###Solar
Field <- "Solar"
fld <- Fld[,(n_site+1):ncol(Fld)]
Get_Simulation_Skill(True_Data = fld,
                     Simulations = solar_sims,
                     Field_Name = Field)

Get_Annual_Cycle(True_Data = fld, Field_Name = Field,
                 Simulations = solar_sims,
                 Resolution = "daily",
                 Start_Date = "01-01-1950")
  
get_pca_plot(X = fld, 
               Grid = grid_locs, 
               Field = Field,
               Sims = solar_sims)
  
Get_PCWavelet_Sim_Skill(Data_Field = fld,
                          Field_Name = Field,
                          Sims = solar_sims,
                          PCs=2)
  

###Cross - Correlation across Sites. 
Get_Site_Correlation(Fld1 = Fld[,1:n_site], #Wind
                     Fld2 = Fld[,(n_site+1):ncol(Fld)], #Solar
                     Fld1_Sims = wind_sims,
                     Fld2_Sims = solar_sims,
                     Grid = grid_locs)
  
Get_Seasonal_Correlation(Fld1 = Fld[,1:n_site], #Wind
                       Fld2 = Fld[,(n_site+1):ncol(Fld)], #Solar
                       Fld1_Sims = wind_sims,
                       Fld2_Sims = solar_sims,
                       Grid = grid_locs,
                       start_date = "01-01-1950",
                       col_hx = "#af8dc3")

  
#Close the plot
dev.off()


#################END OF SCRIPT#################
#______________________________________________________________________________#
#______________________________________________________________________________#
#______________________________________________________________________________#