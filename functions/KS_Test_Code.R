#______________________________________________________________________________#
###Code to apply the get_KS_Test function to both KSTS and KNN Algorithms. 


###Set the Directory
setwd("~/ERCOT") #Code for personal device
source("functions/Get_KS_Test.R")


###Load the Data
ssrd <- read.table("data/ERCOT_Solar_Rad_Daily.txt",  
                   sep =" ", 
                   header = TRUE) #Surface Solar Radiation Downwards
WP <- read.table("data/ERCOT_Wind_Power_Daily.txt",  
                 sep =" ", 
                 header = TRUE) #Wind Capacity Factors

#n <- 2*365 #Subset to 10-yrs
#WP <- WP[1:n,]
#ssrd <- ssrd[1:n,]


#______________________________________________________________________________#
#####KSTS#####

###Load the Simulations
load("Joint_Raw_Simulations.RData") 
KSTS <- ynew_results
ynew_results <- NULL

#Solar
Sims <- list()
for(i in 1:length(KSTS)){
  Sims[[i]] <- KSTS[[i]]$SSnew
}
KSTS_Solar <- get_KS_Test(Dat = ssrd, Sims = Sims, Type = "KSTS-Solar")
Sims <- NA

#Wind
Sims <- list()
for(i in 1:length(KSTS)){
  Sims[[i]] <- KSTS[[i]]$WPnew
}
KSTS_Wind <- get_KS_Test(Dat = WP, Sims = Sims, Type = "KSTS-Wind")
Sims <- NULL
KSTS <- NULL


#______________________________________________________________________________#
#####KSTS#####

###Load the Simulations
load("KNN_Joint_Raw_Simulations.RData") 
KNN <- ynew_results
ynew_results <- NULL

###Solar
Sims <- list()
for(i in 1:length(KNN)){
  Sims[[i]] <- KNN[[i]]$SSnew
}
KNN_Solar <- get_KS_Test(Dat = ssrd, Sims = Sims, Type = "KNN-Solar")
Sims <- NULL

#KNN
Sims <- list()
for(i in 1:length(KNN)){
  Sims[[i]] <- KNN[[i]]$WPnew
}
KNN_Wind <- get_KS_Test(Dat = WP, Sims = Sims, Type = "KNN-Wind")
Sims <- NULL
KNN <- NULL


#______________________________________________________________________________#
Solar <- cbind(KSTS_Solar, KNN_Solar)
Wind <- cbind(KSTS_Wind, KNN_Wind)

Solar$P_Value <- NULL;Solar$P_Value <- NULL  #Has to be done twice
Wind$P_Value <- NULL;Wind$P_Value <- NULL    #For KSTS and KNN.

print(Wind)
print(Solar)

#Converting to LATEX Output
xtable(Solar, caption = "Kolmogorov-Smirnov (KS) test results on aggregated solar simulations generated using KSTS and KNN.",
       digits = 2)

xtable(Wind, caption = "Kolmogorov-Smirnov (KS) test results on aggregated wind simulations generated using KSTS and KNN.",
       digits = 2)
