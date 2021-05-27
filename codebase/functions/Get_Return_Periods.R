#Code to compute and plot the return periods.
#Take a Single 2000-year Simulation. 
#Compute the exceedances
#Fit a poisson regression using locfit. 
#Predict on a grid. 
#Contour it.


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

Fld <- as.matrix(cbind(WP,ssrd))
Fld <- apply(Fld, 2, function(x){x/max(x)})

ssrd <- WP <- NULL

#______________________________________________________________________________#
###Read the simulations
#KSTS
#load("Simulations/Joint_Raw_Simulations_200yrs.RData")
#comb_ksts <- cbind(ynew_results$WPnew, ynew_results$SSnew)
#ynew_results <- NULL
#comb_ksts <- apply(comb_ksts, 2, function(x){x/max(x)})


###Join the Simulations
#load("Simulations/Joint_Raw_Simulations.RData")
#nl <- length(ynew_results)
#comb_ksts <- list()
#for(i in 1:nl){
#  comb_ksts[[i]] <- cbind(ynew_results[[i]]$WPnew, ynew_results[[i]]$SSnew)
#  ynew_results[[i]] <- 1
#}
#ynew_results <- NULL
#comb_ksts <- bind_rows(lapply(comb_ksts,data.frame))
#comb_ksts <- apply(comb_ksts, 2, function(x){x/max(x)})


#______________________________________________________________________________#
###Hyper-Parameters
True_Data = Fld
thresh = 0.30
Start_Date = "01-01-1979"
#Sims_KSTS = comb_ksts; comb_ksts <- NULL




#______________________________________________________________________________#
###Compute the count exceedances

#Source the function.
source("functions/Get_Severity_Duration.R")

#Hyper-Parameters
st_date = Start_Date
mean_daily <- mean(rowSums(True_Data))

#Compute Severity and Duration for Data
Data_SD <- get_Severity_Duration(Data = True_Data, Thresh = thresh,
                                 start_date = st_date, Type = "Data")

#KSTS
#KSTS_SD <- get_Severity_Duration(Data = Sims_KSTS, Thresh = thresh,
#                                 start_date = st_date, Type = "KSTS")
KSTS_SD <- read.table("Sims_SD.txt", header = TRUE, sep = " ")

#COnvert to log-scale
#KSTS_SD <- KSTS_SD %>% filter(Duration > 1)
#KSTS_SD$Severity <- KSTS_SD$Severity/mean_daily
KSTS_SD$Severity <- log10(KSTS_SD$Severity)
KSTS_SD$Duration <- log10(KSTS_SD$Duration)

#Data_SD <- Data_SD %>% filter(Duration > 1)
#Data_SD$Severity <- Data_SD$Severity/mean_daily
Data_SD$Severity <- log10(Data_SD$Severity)
Data_SD$Duration <- log10(Data_SD$Duration)

Comb_SD <- rbind(KSTS_SD, Data_SD)


#Compute the Exceedances for the Simulations
N_years <- 1920 #round(dim(Sims_KSTS)[1]/365)
n_exceed <- list()
for(i in 1:nrow(KSTS_SD)){
  t_count <-  which(KSTS_SD$Severity > KSTS_SD$Severity[i] &
                      KSTS_SD$Duration > KSTS_SD$Duration[i])
  n_exceed[[i]] <- length(t_count)
}
n_exceed <- unlist(n_exceed)  
KSTS_SD$pexc <- n_exceed   #(n_exceed+1)/N_years
#KSTS_SD <- KSTS_SD %>% filter(pexc < 1)

#Compute the Exceedances for the Data
n_exceed <- list()
for(i in 1:nrow(Data_SD)){
  t_count <-  which(Data_SD$Severity > Data_SD$Severity[i] &
                      Data_SD$Duration > Data_SD$Duration[i])
  n_exceed[[i]] <- length(t_count)
}
n_exceed <- unlist(n_exceed)  
Data_SD$pexc <- n_exceed   #(n_exceed+1)/N_years
#Data_SD <- Data_SD %>% filter(pexc < 1)



#______________________________________________________________________________#
###Local Regression for the Simulations
fit <- locfit(pexc ~ Duration + Severity,family = "poisson",
              data=KSTS_SD, alpha = 0.33)

#plot(KSTS_SD$Duration, KSTS_SD$Severity, pch = 19, cex = 0,
#     xlab = "Duration", ylab = "Severity")
#plot(fit, add = TRUE)
#points(KSTS_SD$Duration, KSTS_SD$Severity, pch = 19, cex = 0.5, col ='blue')


#Local Regression for the Data
fit_data <- locfit(pexc ~ Duration + Severity,family = "poisson",
              data=Data_SD, alpha = 0.33)



#______________________________________________________________________________#
###Make Predictions at points across the grid.

#Set-up a grid
griddf <- expand.grid(Duration = seq(from = min(KSTS_SD$Duration), to = 1.5*max(KSTS_SD$Duration), l = 100),
                      Severity = seq(from = min(KSTS_SD$Severity), to = 1.5*max(KSTS_SD$Severity), l = 100))

#Predict for each point on the grid.
griddf$Pred_KSTS <- griddf$Pred_Data <- NA
for(i in 1:nrow(griddf)){
  griddf$Pred_KSTS[i] <- predict(fit, matrix(c(griddf$Duration[i], griddf$Severity[i]), ncol = 2))
  griddf$Pred_Data[i] <- predict(fit_data, matrix(c(griddf$Duration[i], griddf$Severity[i]), ncol = 2))
  
}

#Compute the exceedances per year
griddf$Ex_year_KSTS <- griddf$Pred_KSTS/(N_years) 
griddf$Ex_year_Data <- griddf$Pred_Data/(40) 



#______________________________________________________________________________#
###Make the contour plot
group.colors <- c(KSTS ="#af8dc3", Data = "#000000")

#pdf("Return_Periods_10.pdf")

ggplot(griddf, aes(Duration, Severity)) +
  geom_contour(aes(z = Ex_year_KSTS), breaks = c(0.1, 0.05, 0.02, 0.01), size = 1, col = group.colors[1]) +
  geom_contour(aes(z = Ex_year_Data), breaks = c(0.1, 0.05), size = 1, linetype = "dashed", col = group.colors[2]) +
  geom_point(data = Comb_SD, mapping = aes(x = Duration, y = Severity, color = Type), alpha = 0.5) +
  #annotate("text", x=1.18, -2.2, label= "10yr", col = 'red', size = 3) +
  #annotate("text", x=1.3, -2.2, label= "20yr", col = 'red', size = 3) +
  #annotate("text", x=1.4, -2.2, label= "50yr", col = 'red', size = 3) +
  #annotate("text", x=1.5, -2.2, label= "100yr", col = 'red', size = 3) +
  ggtitle(paste0("Return Periods for Energy Droughts \n Threshold - ", thresh*100,"th Percentile")) + 
  scale_x_continuous(name = "Duration (Days)", labels = scales::math_format(10^.x)) +
  scale_y_continuous(name = "Severity", labels = scales::math_format(10^.x)) +
  theme_bw() +
  theme(legend.text=element_text(size=20),
        legend.title=element_text(size=15),
        axis.text=element_text(size=15),
        axis.title=element_text(size=20),
        plot.title = element_text(size=20)) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  scale_color_manual(values=group.colors)


ggplot(griddf, aes(Duration, Severity)) +
  geom_contour(aes(z = Ex_year_KSTS), breaks = c(0.1, 0.05, 0.02, 0.01), size = 1, col = "red") +
  geom_point(data = KSTS_SD, mapping = aes(x = Duration, y = Severity, color = Type), alpha = 0.5) +
  #annotate("text", x=1.18, -2.2, label= "10yr", col = 'red', size = 3) +
  #annotate("text", x=1.3, -2.2, label= "20yr", col = 'red', size = 3) +
  #annotate("text", x=1.4, -2.2, label= "50yr", col = 'red', size = 3) +
  #annotate("text", x=1.5, -2.2, label= "100yr", col = 'red', size = 3) +
  ggtitle(paste0("Return Periods for Energy Droughts \n Threshold - ", thresh*100,"th Percentile")) + 
  scale_x_continuous(name = "Duration (Days)", labels = scales::math_format(10^.x)) +
  scale_y_continuous(name = "Severity ", labels = scales::math_format(10^.x)) +
  theme_bw() +
  theme(legend.text=element_text(size=20),
        legend.title=element_text(size=15),
        axis.text=element_text(size=15),
        axis.title=element_text(size=20),
        plot.title = element_text(size=20)) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  scale_color_manual(values=group.colors)


ggplot(griddf, aes(Duration, Severity)) +
  geom_contour(aes(z = Ex_year_Data), breaks = c(0.1, 0.05), size = 1, col = "red") +
  geom_point(data = Data_SD, mapping = aes(x = Duration, y = Severity, color = Type), alpha = 0.5) +
  #annotate("text", x=1.18, -2.2, label= "10yr", col = 'red', size = 3) +
  #annotate("text", x=1.3, -2.2, label= "20yr", col = 'red', size = 3) +
  ggtitle(paste0("Return Periods for Energy Droughts \n Threshold - ", thresh*100,"th Percentile")) + 
  scale_x_continuous(name = "Duration (Days)", labels = scales::math_format(10^.x)) +
  scale_y_continuous(name = "Severity", labels = scales::math_format(10^.x)) +
  theme_bw() +
  theme(legend.text=element_text(size=20),
        legend.title=element_text(size=15),
        axis.text=element_text(size=15),
        axis.title=element_text(size=20),
        plot.title = element_text(size=20)) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  scale_color_manual(values=group.colors)


dev.off()


