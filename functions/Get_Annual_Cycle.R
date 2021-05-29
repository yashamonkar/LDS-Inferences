#______________________________________________________________________________#
###Code to visualize the Simulation Skill for Seasonality###
#Note:- Time-Step is a Day.
#The Code Randomly Selects 9 sites to plot. 


###Input
#1. True_Data:- Original Data Matrix
#2. Field_Name:- Name of the Field. e.g Wind
#2. Start_Date:- Day the data starts. e.g 01-01-1970 (or) "01-01-1970 00:00"
#3. Simulations:- List of List for Simulations. 
#4. Resolution - Time Resolution of the Data. e.g. daily or hourly


###Output
#1. ggplot for the annual cycle.


#______________________________________________________________________________#
Get_Annual_Cycle <- function(True_Data, Simulations, Field_Name, 
                             Resolution, Start_Date){
  
  library(ggplot2)
  library(gridExtra)
  
  ylab.text = expression(paste("W/m"^"2"))
  
  ###Develop a time index
  if(Resolution == "daily"){
    st_date <- as.Date(Start_Date, format="%m-%d-%Y")
    time_stamps <- seq(st_date, by = "day", length.out =nrow(True_Data))
    month_stamps <- as.numeric(format(time_stamps,"%m"))
    
    #Sample Sites
    sites <- sample(1:ncol(True_Data), 8, replace = FALSE)
    

    for(i in 1:8){
      temp <- sites[i]
      
      #Dataframe for True Site Data
      DF_Data <- data.frame(Value = True_Data[,temp],
                            Month = month_stamps,
                            Type = "Data")
    
      #Get the Simulations for the Site
      site_sims <- list()
      for(j in 1:length(Simulations)){
        tt <- as.data.frame(Simulations[[j]])
        site_sims[[j]] <- tt[,temp]
      }
      site_sims <- unlist(site_sims)
    
      #Simulation Time Stamps
      sim_time_stamps <- seq(st_date, by = "day", length.out =nrow(tt))
      sim_month_stamps <- as.numeric(format(sim_time_stamps,"%m"))
      
      
      #Dataframe for Simulations
      DF_Sims <- data.frame(Value = site_sims,
                            Month = rep(sim_month_stamps,j),
                            Type = "Simulations")
      
      #Merge two datasets
      DF_Plot <- rbind(DF_Data, DF_Sims)
                     
      #Generating the plot
      paste0("p",i) 
      p <-  ggplot(data=DF_Plot) + 
      geom_boxplot( aes(x=factor(Month), y=Value, fill=Type), 
                    position=position_dodge(1), outlier.shape = NA) +
      ylab(ylab.text) +
      xlab("Month") + 
      ggtitle(paste0(" ", Field_Name, " Distribution for Site/Grid - ", temp)) +
        theme_bw()
      assign(paste0("p",i),p)
    
    }

    grid.arrange(p1, p2, nrow = 2)
    grid.arrange(p3, p4, nrow = 2)
    grid.arrange(p5, p6, nrow = 2)
    grid.arrange(p7, p8, nrow = 2)
    
  }
  
  
  #------------------------------------------------------------------------------#
  if(Resolution == "hourly"){
    st_date <- as.POSIXct(Start_Date, format="%m-%d-%Y %H:%M")
    time_stamps <- seq(st_date, by = "hour", length.out =nrow(True_Data))
    month_stamps <- as.numeric(format(time_stamps,"%m"))
    hour_stamps <- as.numeric(format(time_stamps,"%H"))
    
    #Sample Sites
    sites <- sample(1:ncol(True_Data), 8, replace = FALSE)
    
    
    for(i in 1:8){
      temp <- sites[i]
      
      #Dataframe for True Site Data
      DF_Data <- data.frame(Value = True_Data[,temp],
                            Month = month_stamps,
                            Hour = hour_stamps,
                            Type = "Data")
      
      #Get the Simulations for the Site
      site_sims <- list()
      for(j in 1:length(Simulations)){
        tt <- as.data.frame(Simulations[[j]])
        site_sims[[j]] <- tt[,temp]
      }
      site_sims <- unlist(site_sims)
      
      #Simulation Time Stamps
      sim_time_stamps <- seq(st_date, by = "hour", length.out =nrow(tt))
      sim_month_stamps <- as.numeric(format(sim_time_stamps,"%m"))
      sim_hour_stamps <- as.numeric(format(sim_time_stamps,"%H"))
      
      
      #Dataframe for Simulations
      DF_Sims <- data.frame(Value = site_sims,
                            Month = rep(sim_month_stamps,j),
                            Hour = rep(sim_hour_stamps, j),
                            Type = "Simulations")
      
      #Merge two datasets
      DF_Plot <- rbind(DF_Data, DF_Sims)
      
      #Generating the plot for Month
      paste0("p",i) 
      p <-  ggplot(data=DF_Plot) + 
        geom_boxplot( aes(x=factor(Month), y=Value, fill=Type), 
                      position=position_dodge(1), outlier.shape = NA) +
        ylab(paste0(Field_Name, " Site Value")) +
        xlab("Month") + 
        ggtitle(paste0("Monthly ", Field_Name, " for Site ", temp)) +
        theme_bw()
      assign(paste0("p",i),p)
      
      
      #Generating the plot for Day
      paste0("d",i) 
      d <-  ggplot(data=DF_Plot) + 
        geom_boxplot( aes(x=factor(Hour), y=Value, fill=Type), 
                      position=position_dodge(1), outlier.shape = NA) +
        ylab(paste0(Field_Name, " Site Value")) +
        xlab("Hour") + 
        ggtitle(paste0("Hourly ", Field_Name, " for Site ", temp)) +
        theme_bw()
      assign(paste0("d",i),d)
      
    }
    
    grid.arrange(p1, p2, p3, p4,nrow = 2)
    grid.arrange(p5, p6, p7, p8,nrow = 2)
    grid.arrange(d1, d2, nrow = 2)
    grid.arrange(d3, d4, nrow = 2)
    grid.arrange(d5, d6, nrow = 2)
    grid.arrange(d7, d8, nrow = 2)
    
  }
}