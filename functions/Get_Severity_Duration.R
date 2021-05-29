#______________________________________________________________________________#
###Function to compute Severity and Duration of Energy Droughts for Simulation and Data. 
###Note:- Coded only for deficit.

#Input:-  1. Data - A single Matrix or Dataframe.
#         2. Threshold e.g. 0.3 [Think of it in terms of daily climatological percentile which has to be supplied.]
#         3. Type  e.g. "Data" or "Simulations"
#         4. Start-Date.

#Output:- Severity and Duration in a data frame form. 

#Takes in either Data or a single Simulation Run.
#Output is Data Frame containing Type, Severity and Threshold.


get_Severity_Duration <- function(Data, Thresh, start_date, Type){
  
  #Load Packages
  library(lubridate)
  library(dplyr)
  
  #Clean up the Data
  colnames(Data) <- NULL
  Data <- as.data.frame(Data)
  
  
  #Get the Day Indices
  st_date <- as.Date(start_date, format="%m-%d-%Y")
  time_stamps <- seq(st_date, by = "day", length.out = nrow(Data))
  Data$Dates <- format(time_stamps,"%m-%d")
  
  #Daily Average
  Data_Thresh <- Data %>% group_by(Dates) %>% summarise_all(funs(Q3 = quantile), 
                                                         probs = Thresh)
  Data_Thresh <- as.data.frame(Data_Thresh)
  
  #Keep the Dates
  Dates_Index <- data.frame(Time = time_stamps, Dates = Data$Dates)
  Data$Dates <- NULL
  Data_Thresh <- merge(Dates_Index, Data_Thresh, by = "Dates")
  Data_Thresh <- Data_Thresh[order(Data_Thresh$Time),]
  Data_Thresh$Dates <- Data_Thresh$Time <- NULL

  
  
  #Deviances from Threshold Percentile
  Data_dev <- matrix(NA, ncol = ncol(Data), nrow = nrow(Data))
  for(i in 1:ncol(Data)){
    Data_dev[,i] <- (Data[,i] - Data_Thresh[,i])
  }
  Data_dev <- rowSums(Data_dev)
  
  
  #Cumulative Deficit
  CD <- rep(NA, length(Data_dev))
  CD[1] <- max(0,-Data_dev[1])
  for(i in 2:length(Data_dev)){
    CD[i] <- max(0, CD[i-1] - Data_dev[i])
  }
  
  if(max(CD) > 0){
    #Creating Group-ID runs for non-zero values.
    x <- match(CD != 0, TRUE)
    grp_run <- with(rle(!is.na(x)), {
      lv <- lengths[values]
      replace(x, !is.na(x), rep(seq_along(lv), lv))
    })
  
    #Compute Severity and Duration
    SnD <- data.frame(Deficit = CD, Type = grp_run)
    SnD <- SnD[complete.cases(SnD), ]
    Sev_Dur <- SnD %>% 
      group_by(Type) %>% 
      summarise(Severity = max(Deficit), Duration = n())
    Sev_Dur$Type <- Type 
    } else {
      Sev_Dur <- data.frame(Type = Type,
                            Severity = 1,
                            Duration = 0)
    
  }
  
  return(Sev_Dur)
  
  
}