#______________________________________________________________________________#
###Code to visualize the Simulation Skill for Seasonality###
#Note:- Time-Step is a Day.
#The Code Randomly Selects 9 sites to plot. 


###Input
#1. True_Data:- Original Data Matrix
#2. Field_Name:- Name of the Field. e.g Wind
#2. Start_Date:- Day the data starts. e.g 01-01-1970 (or) "01-01-1970 00:00"
#3. Simulations:- List of List for Simulations. 


###Output
#1. ggplot for the annual cycle.


#______________________________________________________________________________#
Get_Annual_Cycle <- function(True_Data, Simulations, Field_Name, Start_Date){
  
  library(ggplot2)
  library(gridExtra)
  
  ylab.text = expression(paste("CF"))
  
 
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
      ggtitle(Field_Name) +
        theme_bw() +
        theme(axis.text=element_text(size=5),
              axis.title=element_text(size=10),
              plot.title = element_text(size=12),
              legend.title = element_text(size=0),
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank())
      
      assign(paste0("p",i),p)
    
    }

    p <- plot_grid(p1+ theme(legend.position="none"),
                   p2+ theme(legend.position="none"),ncol=1)
    
    #Get Legend
    legend <- get_legend(
      p1 + 
        guides(color = guide_legend(nrow = 1, override.aes = list(size=2))) +
        theme(legend.position = "bottom")
    )
    
    out = list(p=p, legend=legend)
}