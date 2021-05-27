#______________________________________________________________________________#
###Code for Principal Component Analysis Using GGPLOT. 


###INPUT
#1. True Data Feild
#2. Simulated Data Field
#3. Lat-Lon Grids



###OUTPUT
#1. Spatial PC-1 and PC-2.
#2. Boxplots of the Eigenvectors
#3. Spatial of the ensemble Mean of PC-1 and PC-2.
#4. Field Name


#______________________________________________________________________________#
get_pca_plot <- function(X, Grid, Field, Sims){
  
  ###Load Dependencies
  library(ggplot2)
  library(gridExtra)
  
  
  ###PCA on the Data
  pcw <- prcomp(X, scale = TRUE, center = TRUE)
  var <- (pcw$sdev^2)
  var <- (var)/sum(var)
  
  
  ###Plotting the Data PC-1
  world <- map_data("world")
  us <- map_data("state")
  
  
  Grid$Eigenvectors <- scale(pcw$rotation[,1])
  p1 <- ggplot() +
    geom_map(data=world, map=world,
             aes(x=long, y=lat, map_id=region),
             fill="#D3D3D3", color="#000000", size=0.15) +
    geom_map(data=us, map=us,
             aes(x=long, y=lat, map_id=region),
             fill="#D3D3D3", color="#000000", size=0.15) +
    geom_tile(data = Grid, aes(x= lon-360, y = lat, fill = Eigenvectors)) + 
    scale_fill_gradient2(midpoint= 0,
                          low="red", mid="green",high="blue", 
                          limits = c(min(Grid$Eigenvectors), max(Grid$Eigenvectors))) +
    labs(title = paste0("Data PC-1 ", Field), y = " ", x = " ",
         subtitle = paste0("Fractional Variance - ", round(var[1]*100),"%")) + 
    scale_x_continuous(name = "lon", limits = c(-107, -92)) +
    scale_y_continuous(name = "lat", limits = c(25.5, 36.5)) +
    theme_bw() +
    theme(axis.text=element_text(size=0),
          axis.title=element_text(size=0),
          axis.ticks = element_blank())
  
  
  ###Plotting the Data PC-2
  Grid$Eigenvectors <- scale(pcw$rotation[,2])
  p2 <- ggplot() +
    geom_map(data=world, map=world,
             aes(x=long, y=lat, map_id=region),
             fill="#D3D3D3", color="#000000", size=0.15) +
    geom_map(data=us, map=us,
             aes(x=long, y=lat, map_id=region),
             fill="#D3D3D3", color="#000000", size=0.15) +
    geom_tile(data = Grid, aes(x= lon-360, y = lat, fill = Eigenvectors)) + 
    scale_fill_gradient2(midpoint= 0,
                         low="red", mid="green",high="blue", 
                         limits = c(min(Grid$Eigenvectors), max(Grid$Eigenvectors))) +
    labs(title = paste0("Data PC-1 ", Field), y = " ", x = " ",
         subtitle = paste0("Fractional Variance - ", round(var[2]*100),"%")) + 
    scale_x_continuous(name = "lon", limits = c(-107, -92)) +
    scale_y_continuous(name = "lat", limits = c(25.5, 36.5)) +
    theme_bw() +
    theme(axis.text=element_text(size=0),
          axis.title=element_text(size=0),
          axis.ticks = element_blank())
  
  
  ###Boxplot for PC-1 and PC-1
  pc1 <- pc2 <- list()
  n_sim <- length(Sims)
  
  for(i in 1:n_sim){
    temp <- as.data.frame(Sims[[i]])
    temp <- scale(temp)
    pcw_temp <-  prcomp(temp, scale = TRUE, center = TRUE)
    tvar <- (pcw_temp$sdev^2)
    tvar <- (tvar)/sum(tvar)
    pc1[[i]] <- tvar[1]
    pc2[[i]] <- tvar[2]
  }
  sims_pc <- data.frame(PC1 = unlist(pc1), PC2 = unlist(pc2))
  
  
  ###Ensemble Mean
  ens_avg <- as.data.frame(Sims[[1]])
  for(i in 2:n_sim){
    temp <- as.data.frame(Sims[[i]])
    ens_avg <- rbind(ens_avg,temp)
  }
  ens_avg <- scale(ens_avg)
  pcw_temp <-  prcomp(ens_avg, scale = TRUE, center = TRUE)
  tvar <- (pcw_temp$sdev^2)
  tvar <- (tvar)/sum(tvar)
  pc1_ens <- tvar[1]
  pc2_ens <- tvar[2]
  
  #Boxplot PC-1
  data_var <- data.frame(Variance = c(var[1],pc1_ens),
                         Type = c("Data", "Ensemble Mean"))
  
  p3 <- qplot(y=sims_pc$PC1*100, x= 1, geom = "boxplot") +
    geom_hline(aes(yintercept = 100*Variance, col = Type), data_var, size = 1.5) +
    labs(title = paste0("Simulated PC-1 Variance ", Field),
         y = "Fractional Variance (%)", x = " ", size = 1.5) + 
    theme_bw() + 
    theme(axis.ticks.x=element_blank(),
          axis.text.x=element_blank())
  
  #Boxplot PC-2
  data_var <- data.frame(Variance = c(var[2],pc2_ens),
                         Type = c("Data", "Ensemble Mean"))
  
  
  p4 <- qplot(y=sims_pc$PC2*100, x= 1, geom = "boxplot") +
    geom_hline(aes(yintercept = 100*Variance, col = Type), data_var, size = 1.5) +
    labs(title = paste0("Simulated PC-2 Variance ", Field),
         y = "Fractional Variance (%)", x = " ", size = 1.5) + 
    theme_bw() + 
    theme(axis.ticks.x=element_blank(),
          axis.text.x=element_blank())
  
  
  ###Plotting the Ensemble Mean PC-1
  Grid$Eigenvectors <- scale(pcw_temp$rotation[,1])
 
  p5 <- ggplot() +
    geom_map(data=world, map=world,
             aes(x=long, y=lat, map_id=region),
             fill="#D3D3D3", color="#000000", size=0.15) +
    geom_map(data=us, map=us,
             aes(x=long, y=lat, map_id=region),
             fill="#D3D3D3", color="#000000", size=0.15) +
    geom_tile(data = Grid, aes(x= lon-360, y = lat, fill = Eigenvectors)) + 
    scale_fill_gradient2(midpoint= 0,
                          low="red", mid="green",high="blue", 
                          limits = c(min(Grid$Eigenvectors), max(Grid$Eigenvectors))) +
    labs(title = paste0("Ensemble Mean PC-1 ", Field), y = " ", x = " ",
         subtitle = paste0("Fractional Variance - ", round(pc1_ens*100),"%")) + 
    scale_x_continuous(name = "lon", limits = c(-107, -92)) +
    scale_y_continuous(name = "lat", limits = c(25.5, 36.5)) +
    theme_bw() +
    theme(axis.text=element_text(size=0),
          axis.title=element_text(size=0),
          axis.ticks = element_blank())
  
  
  ###Plotting the Ensemble Mean PC-2
  Grid$Eigenvectors <- scale(pcw_temp$rotation[,2])
  
  p6 <- ggplot() +
    geom_map(data=world, map=world,
             aes(x=long, y=lat, map_id=region),
             fill="#D3D3D3", color="#000000", size=0.15) +
    geom_map(data=us, map=us,
             aes(x=long, y=lat, map_id=region),
             fill="#D3D3D3", color="#000000", size=0.15) +
    geom_tile(data = Grid, aes(x= lon-360, y = lat, fill = Eigenvectors)) + 
    scale_fill_gradient2(midpoint= 0,
                         low="red", mid="green",high="blue", 
                         limits = c(min(Grid$Eigenvectors), max(Grid$Eigenvectors))) +
    labs(title = paste0("Ensemble Mean PC-2 ", Field), y = " ", x = " ",
         subtitle = paste0("Fractional Variance - ", round(pc2_ens*100),"%")) +
    scale_x_continuous(name = "lon", limits = c(-107, -92)) +
    scale_y_continuous(name = "lat", limits = c(25.5, 36.5)) +
    theme_bw() +
    theme(axis.text=element_text(size=0),
          axis.title=element_text(size=0),
          axis.ticks = element_blank())
  
  
  grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 3)
}