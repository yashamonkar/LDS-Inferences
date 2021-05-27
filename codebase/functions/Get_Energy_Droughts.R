#______________________________________________________________________________#
###Code to visualize the Energy Droughts###
#Note:- Time Step is Day.


###Input
#1.  Data Field - Fld
#2.  Simulations 
#3.  Start Date
#4.  Threshold
#5.  Field Name



###Output
#1.  ggplot of log(Severity) vs Duration with Marginals
source("functions/Get_Severity_Duration.R")

#______________________________________________________________________________#
get_energy_droughts <- function(True_Data, Field_Name,  #True Data Parameters
                                Sims,
                                thresh,
                                Start_Date)
  {
  
  #Load Dependencies
  library(ggplot2)
  library(lubridate)
  library(cowplot)
  library(gridExtra)
  
  #Hyper-Parameters
  st_date = Start_Date
  nsim <- length(Sims)
  
  #Compute Severity and Duration for Data
  Data_SD <- get_Severity_Duration(Data = True_Data, Thresh = thresh,
                                   start_date = st_date, Type = "Data")
  
  Sim_SD <- list()
  for(i in 1:nsim){
    Sim_Data <- Sims[[i]]
    Sim_SD[[i]] <- get_Severity_Duration(Data = Sim_Data, Thresh = thresh,
                                         start_date = st_date, Type = "Simulation")
    
  }
  Sims <- NULL
  
  #Convert to a Data-Frame Structure
  Sims <- bind_rows(lapply(Sim_SD,data.frame))
  Sev_Dur <- rbind(Sims,Data_SD)
  
  #Remove the 1 and 2 Day Duration Events
  Sev_Dur<- Sev_Dur %>% filter(Duration > 1)
  
  Sev_Dur$Severity <- log(Sev_Dur$Severity)
  Sev_Dur$Duration <- log(Sev_Dur$Duration)
  
  
  p_main <- ggplot(Sev_Dur) + geom_point(aes(Duration, Severity, color = Type)) +
    ggtitle(paste0("Duration vs Severity - ", Field_Name, "\n Threshold - ", thresh*100,"th Percentile")) + 
    xlab("log(Duration)") + ylab("log(Severity)") + 
    theme_bw()
  
  xbox <- axis_canvas(p_main, axis = "x", coord_flip = TRUE) + 
    geom_boxplot(data = Sev_Dur, aes(y = Duration, x = factor(Type), color = factor(Type))) + 
    scale_x_discrete() + coord_flip()
  
  ybox <- axis_canvas(p_main, axis = "y") + 
    geom_boxplot(data = Sev_Dur, aes(y = Severity, x = factor(Type), color = factor(Type))) +
    scale_x_discrete()
  
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