#______________________________________________________________________________#
###Plots for Simulation Checks - Skill Assessments###


###Data Inputs
#1. Wind Data
#2. Solar Data
#3. Grid Locations
#4. KSTS Simulations
#5. KNN Simulations

###Outputs
# All the simulation Check assessments in the patterns paper



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
library(gridExtra)
library(cowplot)

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

###Subset
#n <- 5*365
#WP <- WP[1:n,]
#ssrd <- ssrd[1:n,]



#______________________________________________________________________________#
###Read the simulations###
#setwd("~/GitHub/LDS-Inferences/codebase/patterns_figures")
#KSTS
wind_ksts <- solar_ksts <- list()
load("simulations/KSTS_Joint_Simulations.RData")
nl <- length(ynew_results)
for(i in 1:nl){
  wind_ksts[[i]] <- ynew_results[[i]]$WPnew
  solar_ksts[[i]] <- ynew_results[[i]]$SSnew
  ynew_results[[i]] <- 1
}
ynew_results <- NULL

#Knn
wind_knn <- solar_knn <- list()
load("simulations/KNN_Joint_Simulations.RData")
nl <- length(ynew_results)
for(i in 1:nl){
  wind_knn[[i]] <- ynew_results[[i]]$WPnew
  solar_knn[[i]] <- ynew_results[[i]]$SSnew
  ynew_results[[i]] <- 1
}
ynew_results <- NULL


pdf("Simulation_Checks.pdf")


#______________________________________________________________________________#
### - Plot on boxplots of the base moments
###Input
#1. True Data Field(Fld)
#2. Simulations(Fld_Sims)

###Output
#1. 2x2 Boxplot
get_base_moments <- function(Fld,
                             Fld_Sims){
  
  #Hyper-parameters
  n_sim <- length(Fld_Sims)
  
  #Select the ramdom grid cells
  n_comp <- 20
  n_sel <- sample(1:ncol(Fld), n_comp, replace = FALSE) #Select n random locatins
  
  #Compute the moments
  #Moments
  mean_sim <- sd_sim <- max_sim <- min_sim <- list()
  for(i in 1:n_sim){
    temp <- as.data.frame(Fld_Sims[[i]]) 
    temp <- temp %>% select(n_sel)
    
    mean_sim[[i]] <- apply(temp,2,mean)
    sd_sim[[i]] <- apply(temp,2,sd)
    max_sim[[i]] <- apply(temp,2,max)
    min_sim[[i]] <- apply(temp, 2, min)
  }
  
  #True Values
  true_values <- as.data.frame(Fld) %>% select(n_sel)
  True_Values <- data.frame(ID = 1:n_comp,
                            Mean = colMeans(true_values),
                            SD = apply(true_values,2,sd),
                            Min = apply(true_values,2,min),
                            Max = apply(true_values,2,max))
  
  #Mean
  mean_sim <- matrix(unlist(mean_sim), ncol = n_comp, byrow = TRUE)
  Plt_dataset <- data.frame(Val = unlist(list(mean_sim)),
                            ID = rep(1:n_comp, each = n_sim))
  
  
  p1 <- ggplot(Plt_dataset)+
    geom_boxplot(mapping = aes(x = factor(ID), y=Val)) +
    geom_point(True_Values, 
               mapping = aes(x = factor(ID),y=Mean), 
               color = 'red', size = 1)+
    ggtitle("Mean") + 
    scale_y_continuous(name = "CF") +
    theme_bw() +
    theme(axis.text=element_text(size=5),
          axis.title=element_text(size=5),
          plot.title = element_text(size=10),
          axis.ticks.x = element_blank(),
          axis.text.x=element_text(size=0),
          axis.title.x=element_text(size=0),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank())
  
  
  #Mean
  sd_sim <- matrix(unlist(sd_sim), ncol = n_comp, byrow = TRUE)
  Plt_dataset <- data.frame(Val = unlist(list(sd_sim)),
                            ID = rep(1:n_comp, each = n_sim))
  
  
  p2 <- ggplot(Plt_dataset)+
    geom_boxplot(mapping = aes(x = factor(ID), y=Val)) +
    geom_point(True_Values, 
               mapping = aes(x = factor(ID),y=SD), 
               color = 'red', size = 1)+
    ggtitle("Stan-Dev") + 
    scale_y_continuous(name = "CF") +
    theme_bw() +
    theme(axis.text=element_text(size=5),
          axis.title=element_text(size=5),
          plot.title = element_text(size=10),
          axis.ticks.x = element_blank(),
          axis.text.x=element_text(size=0),
          axis.title.x=element_text(size=0),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank())
  
  #Minimum
  min_sim <- matrix(unlist(min_sim), ncol = n_comp, byrow = TRUE)
  Plt_dataset <- data.frame(Val = unlist(list(min_sim)),
                            ID = rep(1:n_comp, each = n_sim))
  
  p3 <- ggplot(Plt_dataset)+
    geom_boxplot(mapping = aes(x = factor(ID), y=Val)) +
    geom_point(True_Values, 
               mapping = aes(x = factor(ID),y=Min), 
               color = 'red', size = 1)+
    ggtitle("Minimum") + 
    scale_y_continuous(name = "CF") +
    theme_bw() +
    theme(axis.text=element_text(size=5),
          axis.title=element_text(size=5),
          plot.title = element_text(size=10),
          axis.ticks.x = element_blank(),
          axis.text.x=element_text(size=0),
          axis.title.x=element_text(size=0),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank())
  
  
  #Maximum
  max_sim <- matrix(unlist(max_sim), ncol = n_comp, byrow = TRUE)
  Plt_dataset <- data.frame(Val = unlist(list(max_sim)),
                            ID = rep(1:n_comp, each = n_sim))
  
  p4 <- ggplot(Plt_dataset)+
    geom_boxplot(mapping = aes(x = factor(ID), y=Val)) +
    geom_point(True_Values, 
               mapping = aes(x = factor(ID),y=Max), 
               color = 'red', size = 1)+
    ggtitle("Maximum") + 
    scale_y_continuous(name = "CF") +
    theme_bw() +
    theme(axis.text=element_text(size=5),
          axis.title=element_text(size=5),
          plot.title = element_text(size=10),
          axis.ticks.x = element_blank(),
          axis.text.x=element_text(size=0),
          axis.title.x=element_text(size=0),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank())
    
  p <- plot_grid(p1,p2,p3,p4,
                 ncol = 2)
  
  return(p)
  
  
}


###Plotting the results
p1 <- get_base_moments(Fld = WP,
                       Fld_Sims = wind_ksts)

p2 <- get_base_moments(Fld = WP,
                       Fld_Sims = wind_knn)

p3 <- get_base_moments(Fld = ssrd,
                       Fld_Sims = solar_ksts)

p4 <- get_base_moments(Fld = ssrd,
                       Fld_Sims = solar_ksts)

plot_grid(p1,p2,p3,p4,
          ncol = 2,
          labels = c('A', 'B', 'C', 'D'), 
          label_size = 12)


#______________________________________________________________________________#
### - Quantiles
###Input
#1. True Data Field(Fld)
#2. Simulations(Fld_Sims)

###Output
#1. 2x2 Boxplot


get_quantile_plots <- function(Fld,Fld_Sims,Field_Name){
  
  #Hyper-parameters
  n_sim <- length(Fld_Sims)
  
  #Define quantiles needed
  qts <- c(0.01,0.05,0.1,0.25,0.75,0.9,0.95,0.99)
  ept <- c("st","th", "th", "th", "th","th", "th", "th")
  
  for(i in 1:length(qts)){
    
    #Legend
    lgd <- paste0(100*qts[i], ept[i], " quantile \n", Field_Name)
    
    #Create the plotting dataset
    #Set-up the storage
    Plt_dataset <- data.frame(Reanalysis = rep(NA, ncol(Fld)),
                              Median_Corr = rep(NA, ncol(Fld)),
                              upper = rep(NA, ncol(Fld)),
                              lower = rep(NA, ncol(Fld)),
                              Legend1 = lgd)
    #Compute data quantiles
    Plt_dataset$Reanalysis <- apply(Fld,2, quantile, probs = qts[i])
    
    #Compute the simulation quatiles
    temp <- matrix(NA, ncol = ncol(Fld), nrow = n_sim)
    for(k in 1:n_sim){
      temp[k,] <- apply(Fld_Sims[[k]], 2, quantile, probs = qts[i])
      }
    Plt_dataset$Median_Corr <- apply(temp,2, quantile, probs = 0.5)
    Plt_dataset$upper <- apply(temp,2, quantile, probs = 0.95)
    Plt_dataset$lower <- apply(temp,2, quantile, probs = 0.05)
    
    
    #Plotting the results
   p <- ggplot(Plt_dataset) +
      geom_abline(intercept = 0, slope = 1, size = 1) +
      geom_linerange(mapping = aes(x=Reanalysis, ymin = lower, 
                                   ymax = upper),size = 0.5, color = 'red') +
      geom_point(mapping = aes(x = Reanalysis, y = Median_Corr),
                 color ='blue', size = 0.75) +
      scale_x_continuous(name = "Data") +
      scale_y_continuous(name = "Simulation") +
     annotate(geom="text", label=lgd, size = 2,
              x = -Inf, y = Inf, 
              hjust = -0.1, vjust = 1.1) +
     theme_bw() +
     theme(axis.text=element_text(size=5),
            axis.title=element_text(size=8),
            plot.title = element_text(size=12),
            legend.position = c(0.2,0.9),
            legend.title=element_text(size=0),
            legend.text =element_text(size=7),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank())
    
   assign(paste0("p",i), p)
  }
  
  p_lower <- plot_grid(p1,p2,p5,p6,
            ncol = 2)
  p_upper <- plot_grid(p3,p4,p7,p8,
                       ncol = 2)
  
  out = list(p_lower = p_lower,
             p_upper = p_upper)
  
}

#Plotting the results
wind <- get_quantile_plots(Fld = WP,
                   Fld_Sims = wind_ksts,
                   Field_Name = "Wind - KSTS")

solar <- get_quantile_plots(Fld = ssrd,
                            Fld_Sims = solar_ksts,
                            Field_Name = "Solar - KSTS")

plot_grid(wind$p_lower,wind$p_upper,
          solar$p_lower,solar$p_upper,
          ncol = 2,
          labels = c('A', ' ', 'B', ' '), 
          label_size = 12)


wind <- get_quantile_plots(Fld = WP,
                           Fld_Sims = wind_knn,
                           Field_Name = "Wind - KNN")

solar <- get_quantile_plots(Fld = ssrd,
                            Fld_Sims = solar_knn,
                            Field_Name = "Solar - KNN")

plot_grid(wind$p_lower,wind$p_upper,
          solar$p_lower,solar$p_upper,
          ncol = 2,
          labels = c('A', ' ', 'B', ' '), 
          label_size = 12)


#______________________________________________________________________________#
### - Site Correlaion

#Load the function
source("plotting_functions/Get_Site_Correlation.R")

Get_Site_Correlation(Fld1 = WP,Fld2 = ssrd,
                     Fld1_Sims = wind_ksts,
                     Fld2_Sims = solar_ksts,
                     Grid = grid_locs)

Get_Site_Correlation(Fld1 = WP,Fld2 = ssrd,
                     Fld1_Sims = wind_knn,
                     Fld2_Sims = solar_knn,
                     Grid = grid_locs)



#______________________________________________________________________________#
### Get PDF plots for a single site
get_pdf <- function(Fld, Fld_Sim,Field_Name){
  
  #Hyper-Parameters
  nsim <- length(Fld_Sim)
  
  #Select a site at random
  sel_site <- sample(1:ncol(Fld),1)
  
  #Compute the KDE on the Reanalysis Data
  og_pdf <- density(Fld[,sel_site], from = 0, to = 1)
  
  #Compute the Simulation PDF
  sim_pdf <- matrix(NA, ncol = nsim, nrow = length(og_pdf$x))
  for(j in 1:nsim){
    #Computing each CDF
    sim <- Fld_Sim[[j]][,sel_site]
    pdf_sim <- density(sim, from = 0, to = 1)
    sim_pdf[,j] <- pdf_sim$y
  }
  lower_pdf <- apply(sim_pdf, 1, function(x) quantile(x, probs=.01))
  upper_pdf <- apply(sim_pdf, 1, function(x) quantile(x, probs=.99))
  median_pdf <- apply(sim_pdf, 1, median)
  
  #Plotting Dataset
  Plt_Data <- data.frame(X = rep(og_pdf$x,2),
                         lower = c(lower_pdf, lower_pdf),
                         upper = c(upper_pdf, upper_pdf))
  
  OG_PDF <- data.frame(X= rep(og_pdf$x,2),
                       Y=c(og_pdf$y,median_pdf),
                       Type = rep(c("Data", "Simulations"), each = length(og_pdf$x)))
  
  group.colors <- c(Data ="red", Simulations = "black")
  
  p <- ggplot(Plt_Data, aes(X)) +
    geom_ribbon(aes(ymin = lower, ymax =  upper), fill ='grey', alpha = 0.85) +
    geom_line(OG_PDF, mapping = aes(x = X, y = Y, color = Type), size = 0.9) +
    scale_x_continuous(name = "CF") +
    scale_y_continuous(name = "Density - f(x)") +
    labs(title = paste0("PDF - Single Site - ", Field_Name)) + 
    theme_bw() +
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=10),
          plot.title = element_text(size=12),
          legend.position = c(0.81,0.85),
          legend.title=element_text(size=0),
          legend.text =element_text(size=8),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank()) +
    scale_color_manual(values=group.colors)
  
  return(p)
  
}


p1 <- get_pdf(Fld = WP,
              Fld_Sim = wind_ksts,
              Field_Name = "Wind - KSTS")


p2 <- get_pdf(Fld = ssrd,
              Fld_Sim = solar_ksts,
              Field_Name = "Solar - KSTS")


p3 <- get_pdf(Fld = WP,
              Fld_Sim = wind_knn,
              Field_Name = "Wind - KNN")


p4 <- get_pdf(Fld = ssrd,
              Fld_Sim = solar_knn,
              Field_Name = "Solar - KNN")


plot_grid(p1,p3,p2,p4,
          ncol = 2,
          labels = c('A', 'B', 'C', 'D'), 
          label_size = 12)


#______________________________________________________________________________#
### Get Auto Correlation Factor plots

get_acf <- function(Fld,Fld_Sims,Field_Name){
  
  #Hyper-parameters
  n_sim <- length(Fld_Sims)
  
  #Define quantiles needed
  Lags <- 1:4
  
  for(i in 1:length(Lags)){
    
    #Create the plotting dataset
    #Set-up the storage
    Plt_dataset <- data.frame(Reanalysis = rep(NA, ncol(Fld)),
                              Median_Corr = rep(NA, ncol(Fld)),
                              upper = rep(NA, ncol(Fld)),
                              lower = rep(NA, ncol(Fld)))
    
    
    for(j in 1:ncol(Fld)){
      
      #Compute the data correlations
      tx <- acf(Fld[,j], plot = FALSE)
      Plt_dataset$Reanalysis[j] <- tx$acf[i+1]
    
      #Compute the simulation quatiles
      temp <- list()
      for(k in 1:n_sim){
        tx <- acf(Fld_Sims[[k]][,j], plot = FALSE)
        temp[[k]] <- tx$acf[i+1]
      }
      Plt_dataset$Median_Corr[j] <- quantile(unlist(temp), probs = 0.5)
      Plt_dataset$upper[j] <- quantile(unlist(temp), probs = 0.95)
      Plt_dataset$lower[j] <- quantile(unlist(temp), probs = 0.05)
    }
    
    
    #Plotting the results
    p <- ggplot(Plt_dataset) +
      geom_abline(intercept = 0, slope = 1, size = 1) +
      geom_linerange(mapping = aes(x=Reanalysis, ymin = lower, 
                                   ymax = upper),size = 0.5, color = 'red') +
      geom_point(mapping = aes(x = Reanalysis, y = Median_Corr),
                 color ='blue', size = 0.75) +
      scale_x_continuous(name = "Data ACF",
                         limits = c(min(Plt_dataset$lower), 
                                    max(Plt_dataset$upper))) +
      scale_y_continuous(name = "Simulation ACF",
                         limits = c(min(Plt_dataset$lower), 
                                    max(Plt_dataset$upper))) +
      ggtitle( paste0("Lag - ", Lags[i]," - ", Field_Name)) + 
      theme_bw() +
      theme(axis.text=element_text(size=5),
            axis.title=element_text(size=8),
            plot.title = element_text(size=8),
            legend.position = c(0.2,0.9),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank())
    
    assign(paste0("p",i), p)
  }
  

  p <- plot_grid(p1,p2,p3,p4,
                 ncol = 2)
 
  return(p)
  
}


p1 <- get_acf(Fld = WP,
              Fld_Sim = wind_ksts,
              Field_Name = "Wind - KSTS")


p2 <- get_acf(Fld = ssrd,
              Fld_Sim = solar_ksts,
              Field_Name = "Solar - KSTS")


p3 <- get_acf(Fld = WP,
              Fld_Sim = wind_knn,
              Field_Name = "Wind - KNN")


p4 <- get_acf(Fld = ssrd,
              Fld_Sim = solar_knn,
              Field_Name = "Solar - KNN")


plot_grid(p1,p3,p2,p4,
          ncol = 2,
          labels = c('A', 'B', 'C', 'D'), 
          label_size = 15)



#______________________________________________________________________________#
#Monthly Boxplots

#Load the function
source("plotting_functions/Get_Annual_Cycle.R")

p1 <- Get_Annual_Cycle(True_Data = WP, 
                 Simulations = wind_ksts,
                 Field_Name = "Wind - KSTS",
                 Start_Date = "01-01-1950")

p2 <- Get_Annual_Cycle(True_Data = WP, 
                              Simulations = wind_knn,
                              Field_Name = "Wind - KNN",
                              Start_Date = "01-01-1950")

p3 <- Get_Annual_Cycle(True_Data = ssrd, 
                              Simulations = solar_ksts,
                              Field_Name = "Solar - KSTS",
                              Start_Date = "01-01-1950")

p4 <- Get_Annual_Cycle(True_Data = ssrd, 
                               Simulations = solar_knn,
                               Field_Name = "Solar - KNN",
                               Start_Date = "01-01-1950")

p_monthly <- plot_grid(p1$p,p2$p,p3$p,p4$p,
          ncol = 2,
          labels = c('A', 'B', 'C', 'D'), 
          label_size = 15)
plot_grid(p_monthly, p1$legend, ncol = 1, rel_heights = c(1, .1))





#______________________________________________________________________________#
#Cross Correlation
get_cross_correlation <- function(Fld,Fld_Sims,Field_Name){
  
  #Hyper-parameters
  n_sim <- length(Fld_Sims)
  n_sites <- 40
  
  #Subset the sites
  sel_sites <- sample(1:ncol(Fld), n_sites, replace = FALSE)
  
  #Set up the plotting dataset
  Plt_dataset <- data.frame(Reanalysis = rep(NA, n_sites^2),
                            Median_Corr = rep(NA, n_sites^2),
                            upper = rep(NA, n_sites^2),
                            lower = rep(NA, n_sites^2))
  
  #Compute true value correlations
  Fld_Subset <- as.data.frame(Fld) %>% select(sel_sites)
  ct <- 1
  for(i in 1:n_sites){
    for(j in 1:n_sites){
      Plt_dataset$Reanalysis[ct] <- cor(Fld_Subset[,i],Fld_Subset[,j])
    ct = ct+1
    }
  }
  
  #Compute the Simulation Correlations
  ct <- 1
  pb = txtProgressBar(min = 1, max = n_sites, initial = 1) 
  for(i in 1:n_sites){
    setTxtProgressBar(pb,i)
    for(j in 1:n_sites){
      temp <- list()
      for(k in 1:n_sim){
        Sims_Subset <- as.data.frame(Fld_Sims[[k]]) %>% select(sel_sites)
        temp[[k]] <- cor(Sims_Subset[,i],Sims_Subset[,j])
      }
      Plt_dataset$Median_Corr[ct] <- quantile(unlist(temp), probs = 0.5)
      Plt_dataset$upper[ct] <- quantile(unlist(temp), probs = 0.95)
      Plt_dataset$lower[ct] <- quantile(unlist(temp), probs = 0.05)
      
      ct = ct+1
    }
  }
  
  
  
    
    #Plotting the results
    p <- ggplot(Plt_dataset) +
      geom_abline(intercept = 0, slope = 1, size = 1) +
      geom_linerange(mapping = aes(x=Reanalysis, ymin = lower, 
                                   ymax = upper),size = 0.5, color = 'red') +
      geom_point(mapping = aes(x = Reanalysis, y = Median_Corr),
                 color ='blue', size = 0.75) +
      scale_x_continuous(name = "Data Correlation",
                         limits = c(min(Plt_dataset$lower), 1)) +
      scale_y_continuous(name = "Simulation Correlation",
                         limits = c(min(Plt_dataset$lower), 1)) +
      annotate(geom="text", label=Field_Name,
               x = -Inf, y = Inf, 
               hjust = -0.15, vjust = 1.5) +
      ggtitle("Data vs Simulation - Cross Site Correlation") + 
      theme_bw() +
      theme(axis.text=element_text(size=5),
            axis.title=element_text(size=10),
            plot.title = element_text(size=10),
            legend.position = c(0.2,0.9),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank())

  
  return(p)
  
}

#Plotting the results
p1 <- get_cross_correlation(Fld <- WP,
                      Fld_Sims <- wind_ksts,
                      Field_Name <- "Wind - KSTS")

p2 <- get_cross_correlation(Fld <- WP,
                            Fld_Sims <- wind_knn,
                            Field_Name <- "Wind - KNN")

p3 <- get_cross_correlation(Fld <- ssrd,
                            Fld_Sims <- solar_ksts,
                            Field_Name <- "Solar - KSTS")

p4 <- get_cross_correlation(Fld <- ssrd,
                            Fld_Sims <- solar_knn,
                            Field_Name <- "Solar - KNN")



plot_grid(p1,p2,p3,p4,
          ncol = 2,
          labels = c('A', 'B', 'C', 'D'), 
          label_size = 15)


#______________________________________________________________________________#
#Seasonal Correlations

#Load the function
source("plotting_functions/Get_Seasonal_Correlation.R")


Get_Seasonal_Correlation(Fld1 = WP,Fld2= ssrd,
                         Fld1_Sims = wind_ksts,
                         Fld2_Sims = solar_ksts, 
                         Grid = grid_locs,
                         start_date = "01-01-1950",
                         col_hx = "#af8dc3")


Get_Seasonal_Correlation(Fld1 = WP,Fld2= ssrd,
                         Fld1_Sims = wind_knn,
                         Fld2_Sims = solar_knn, 
                         Grid = grid_locs,
                         start_date = "01-01-1950",
                         col_hx = "#7fbf7b")





dev.off()
