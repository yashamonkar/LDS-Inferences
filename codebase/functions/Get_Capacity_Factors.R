#______________________________________________________________________________#
#Function to plot the Capacity Factors for Wind.

#Input
#1. Wind Data. 
#2. Wind Simulations.
#3. Start Date

#Output
#1. 2 x 2 Grids of PDF of Daily Capacity Factors
#2. 2 x 2 Grids of Annual Capacity Factors


#______________________________________________________________________________#
Get_Capacity_Factors <- function(Dat, Sims, start_date){
  
  #Load Dependencies
  library(ggplot2)
  library(gridExtra)
  library(dplyr)
  
  #Land Area Constans
  area <- 226800 #Land Area in m-sq
  inst_capacity <- 2*10^6
  
  #Select the Sites
  n_sites <- sample(1:ncol(Dat), 8, replace = FALSE)
  n_sims <- length(wind_sims)
  
  #----------------------------------------------------------------------------#
  ###DAILY
  for(ns in 1:length(n_sites)){
    
    j <- n_sites[ns]
    
    #Data
    X <- data.frame(Val = Dat[,j], Type = "Data")
    X$Val <- X$Val*area/inst_capacity
    
    #Simulations
    y <- list()
    for(i in 1:n_sims){
      y[[i]] <- wind_sims[[i]][,j]
    }
    Y <- data.frame(Val = unlist(y), Type = "Sims")
    Y$Val <- Y$Val*area/inst_capacity
    
    #Concatenate
    ddf <- rbind(X,Y)
    
    #Plot
    p <- ggplot() + 
      geom_density(data=ddf, aes(x=Val, group=Type, fill=Type),alpha=0.5, adjust=2) + 
      xlab("Daily Capacity Factor") +
      ylab("Density") + 
      ggtitle(label = paste0("CF Site - ",j)) +
      xlim(0, 1) +
      theme_bw() + 
      theme(axis.text.y=element_text(size=5),
            axis.text.x=element_text(size=5),
            axis.title=element_text(size=10),
            plot.title = element_text(size = 15))
    
    assign(paste0("p", ns), p)
      
  }
  
  #Plot the results
  grid.arrange(p1, p2, p3,
               p4, nrow = 2)
  
  grid.arrange(p5, p6, p7,
               p8, nrow = 2)
  
  #----------------------------------------------------------------------------#
  ###ANNUAL
  
  #Dates
  st_date <- as.Date(start_date, format="%m-%d-%Y")
  time_stamps <- seq(st_date, by = "day", length.out =nrow(Dat))
  year_stamps <- as.numeric(format(time_stamps,"%Y"))
  
  for(ns in 1:length(n_sites)){
    j <- n_sites[ns]
    
    #Data
    X <- data.frame(Val = Dat[,j], Year = year_stamps)
    X <- X %>% group_by(Year) %>% summarise(Val = sum(Val))
    X$Val <- X$Val*area/(365*inst_capacity)
    X$Type <- "Data";X$Year <- NULL
    
    #Simulations
    y <- list()
    for(i in 1:n_sims){
      y_temp <- wind_sims[[i]][,j]
      y_temp <- data.frame(Val = unlist(y_temp), Year = year_stamps)
      y_temp <- y_temp %>% group_by(Year) %>% summarise(Val = sum(Val))
      y_temp$Val <- y_temp$Val*area/(365*inst_capacity)
      y[[i]] <- y_temp$Val
    }
    Y <- data.frame(Val = unlist(y), Type = "Sims")
    
    #Concatenate
    ddf <- rbind(X,Y)
    
    #Plot
    q <- ggplot() + 
      geom_density(data=ddf, aes(x=Val, group=Type, fill=Type),alpha=0.5, adjust=2) + 
      xlab("Annual Capacity Factor") +
      ylab("Density") + 
      ggtitle(label = paste0("Annual CF Site - ",j)) +
      xlim(0, 0.5) +
      theme_bw() + 
      theme(axis.text.y=element_text(size=5),
            axis.text.x=element_text(size=10),
            axis.title=element_text(size=10),
            plot.title = element_text(size = 15))
    
    assign(paste0("q", ns), q)
    
  }
  
  #Plot the results
  grid.arrange(q1, q2, q3,
               q4, nrow = 2)
  
  grid.arrange(q5, q6, q7,
               q8, nrow = 2)
  
}