#______________________________________________________________________________#
#Function to Visualize Simulation Skill for Multi-Site Data

#Input:- 
#True_Data - The Original Data Set. 
#          - Data Frame. (Required Format)
#          - Columns are Sites, Rows are observations at time t. 

#Simulations - The Generated Simulations
#            - List of Lists (Required Format)
#            - Each Simulation Run is in a format similar to True Data. 

#Field_Name - The name of the field. e.g "Wind"


#Output - Plots of Simulation Skill for different metrics. 

Get_Simulation_Skill <- function(True_Data, Simulations, Field_Name){

library(dplyr)
library(ggplot2)
library(maps)
  
#Variables Provided
X <- True_Data
Sims <- Simulations
Data_Type <- Field_Name


#Variables Derived
n_sim <- length(Sims)

#______________________________________________________________________________#
###Moments
n_comp <- 20
n_sel <- sample(1:ncol(X), n_comp, replace = FALSE) #Select n random locatins

#Moments
mean_sim <- sd_sim <- max_sim <- min_sim <- list()
for(i in 1:n_sim){
  temp <- as.data.frame(Sims[[i]]) 
  temp <- temp %>% select(n_sel)
  
  mean_sim[[i]] <- apply(temp,2,mean)
  sd_sim[[i]] <- apply(temp,2,sd)
  max_sim[[i]] <- apply(temp,2,max)
  min_sim[[i]] <- apply(temp, 2, min)
}

mean_sim <- matrix(unlist(mean_sim), ncol = n_comp, byrow = TRUE)
sd_sim <- matrix(unlist(sd_sim), ncol = n_comp, byrow = TRUE)
max_sim <- matrix(unlist(max_sim), ncol = n_comp, byrow = TRUE)
min_sim <- matrix(unlist(min_sim), ncol = n_comp, byrow = TRUE)
true_values <- as.data.frame(X) %>% select(n_sel)

ylab.text = expression(paste("W/m"^"2"))

par(mfrow = c(2,2), mar = c(3,5,4,2))
boxplot(mean_sim, use.cols = TRUE, 
        main = paste0('Mean - ', Data_Type ), 
        xaxt='n', ylab = ylab.text,
        cex.lab = 1.75)
points(1:n_comp, apply(true_values,2,mean), pch = 19, col = 'red', cex = 1)


boxplot(sd_sim, use.cols = TRUE, 
        main = paste0('Stan Dev - ', Data_Type ), 
        xaxt='n', ylab = ylab.text,
        cex.lab = 1.75)
points(1:n_comp, apply(true_values,2,sd), pch = 19, col = 'red', cex = 1)


boxplot(max_sim, use.cols = TRUE, 
        main = paste0('Maximum - ', Data_Type ), 
        xaxt='n', ylab = ylab.text,
        cex.lab = 1.75)
points(1:n_comp, apply(true_values,2,max), pch = 19, col = 'red', cex = 1)

boxplot(min_sim, use.cols = TRUE, 
        main = paste0('Minimum - ', Data_Type ), 
        xaxt='n', ylab = ylab.text,
        cex.lab = 1.75)
points(1:n_comp, apply(true_values,2,min), pch = 19, col = 'red', cex = 1)
par(mfrow = c(1,1), mar = c(3,3,4,2))

#______________________________________________________________________________#
###Code for Correlation Matrix
n_comp <- 40
n_sel <- sample(1:ncol(X), n_comp, replace = FALSE) #Select n random locatins

par(mfrow = c(2,2))
source("functions/Get_Cross_Correlation.R")
temp <- as.data.frame(X) %>% select(n_sel)
get_crosscorrelations(input_data = temp,
                      nam = " Data")


#Selecting a single Simulation
sel_sim <- sample(n_sim, 3)

temp <- Sims[[sel_sim[1]]]
temp <- as.data.frame(temp) %>% select(n_sel)
get_crosscorrelations(input_data = temp,
                      nam = paste0("Simulation Run - ", sel_sim[1]))

temp <- Sims[[sel_sim[2]]]
temp <- as.data.frame(temp) %>% select(n_sel)
get_crosscorrelations(input_data = temp,
                      nam = paste0("Simulation Run - ", sel_sim[2]))

temp <- Sims[[sel_sim[3]]]
temp <- as.data.frame(temp) %>% select(n_sel)
get_crosscorrelations(input_data = temp,
                      nam = paste0("Simulation Run - ", sel_sim[3]))

par(mfrow = c(1,1))
#####################

#------------------------------------------------------------------------------#
#Correlation for True Data
x <- as.data.frame(X) 
M <- cor(x)

#Correlations for Simulations
M_sim <- list()
for(i in 1:n_sim){
  sim <- Sims[[i]]
  sim <- as.data.frame(sim)
  M_sim[[i]] <- cor(sim)}
MM <- Reduce("+", M_sim) / length(M_sim)

par(mfrow=c(1,1),mar = c(4.5,4.5,4,2))
plot(M,MM, main = paste0("Cross Site Correlations - ", Data_Type),
     xlab = "Data (Reanalysis) Correlation", ylab = "Simulation Correlation",
     cex.lab = 1.25, cex.main=1.5,
     ylim = c((min(M)-0.05), (max(M)+0.05)),
     xlim = c((min(M)-0.05), (max(M)+0.05))
     )
abline(coef = c(0,1), lwd = 2, col='red')
##################

#------------------------------------------------------------------------------#
###Stick Plot for Correlation
#True Correlation
x <- as.data.frame(X) 
M <- cor(x)

#Correlations for Simulations
M_sim <- list()
for(i in 1:n_sim){
  sim <- Sims[[i]]
  sim <- as.data.frame(sim)
  M_sim[[i]] <- cor(sim)}

quant <- c(0.05, 0.25,0.75, 0.95) #Region of Interest

#Setting up the Plotting Region
plot(M,M_sim[[1]], main = paste0("Cross Site Correlations - ", Data_Type),
     xlab = "Reanalysis Data Correlation", ylab = "Simulation Correlation", type='n',
     cex.lab = 1.75, cex.main=1.5,
     ylim = c((min(M)-0.05), (max(M)+0.05)),
     xlim = c((min(M)-0.05), (max(M)+0.05)))
abline(coef = c(0,1), lwd = 3)

#Sub-Sampling for better graphical Optics
i_sub <- sample(ncol(M), 40)
j_sub <- sample(ncol(M), 40)


for(it in 1:length(i_sub)){
  i <- i_sub[it]
  for(jt in 1:length(j_sub)){
    j <- j_sub[jt]
    tcor <- list() 
    for(k in 1:n_sim){
      tcor[k] <- M_sim[[k]][i,j]
    }
    tcor <- unlist(tcor)
    med <- median(tcor)
    qt <- quantile(tcor, prob = quant)
    segments(M[i,j], qt[1], M[i,j], qt[4], col = "red")
    #segments(M[i,j], qt[2], M[i,j], qt[3], col = "black")
    points(M[i,j],med, pch = 19, col ='blue') 
  }
}

legend('topleft', legend = c("Median Simulation Value","5 - 95 Percentile Range"),
       col = c('blue','red'),
       lty = c(0,1), pch = c(19,NA), lwd = c(NA, 3), cex = 1.75)


#______________________________________________________________________________#
#Auto-Correlation
rand_count <- sample(1:ncol(X),1)
tx <- X[,rand_count]
og_acf <- acf(tx, plot = FALSE)
consolid_sims <- list()
consolid_points <- list()
acf_significance_level <- qnorm((1 + 0.95)/2)/sqrt(sum(!is.na(tx)))

for(j in 1:n_sim) {
  temp <- as.data.frame(Sims[[j]])
  temp <- temp[,rand_count]
  t_acf <- acf(temp, plot=FALSE)
  consolid_sims[[j]] <- t_acf$acf
  consolid_points[[j]] <- t_acf$lag
}
consolid_points <- unlist(consolid_points)
consolid_sims <- unlist(consolid_sims)
consolid_df <- data.frame(x=consolid_points,y=consolid_sims)
consolid_df$x <- cut(consolid_df$x, breaks = seq(-1,20,1))
mid_points <- cbind(lower = as.numeric( sub("\\((.+),.*", "\\1", consolid_df$x) ),
                    upper = as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", consolid_df$x) ))
consolid_df$x <- rowMeans(mid_points)
og_df <- data.frame(x1=og_acf$lag, y1= og_acf$acf)


consolid_df$Type <- "Simulations"
og_df$Type <- "Data"


plo <- ggplot()+
  scale_x_continuous(limits=c(-1,10.5)) +
  scale_y_continuous(limits=c(-0.1,1.05)) +
  geom_boxplot(consolid_df, mapping = aes(y = y, x = x+0.5,group = x, fill = Type),
               outlier.shape = NA, outlier.colour = NA)+ 
  ggtitle(paste0("ACF - Single Random Site - ", Data_Type )) +
  ylab("Auto-Correlation") +
  xlab(paste0("Lags (Days)")) +
  theme_bw(base_size = 15)


plo <- plo + 
  geom_point(og_df, mapping = aes(x=x1,y=y1), size = 3, color = "blue")+
  geom_hline(aes(yintercept= acf_significance_level),
             size=1.15, linetype = "dashed")+
  geom_hline(aes(yintercept= -acf_significance_level),
             size=1.15, linetype = "dashed")

print(plo)






#------------------------------------------------------------------------------#
#Auto-Correlation
temp <- X
data_acf <- matrix(NA, ncol = ncol(X), nrow = 10)

for(i in 1:length(n_sel)){
  t_acf <- acf(temp[,i], plot = FALSE)
  data_acf[,i] <- t_acf$acf[2:11]
}

sim_acf <- list()
for(j in 1:n_sim){
  t_acf <- matrix(NA, ncol = ncol(X), nrow = 10)
  for(i in 1:ncol(X)){
    sim <- as.data.frame(Sims[[j]])
    sim <- sim[,i]
    sim <- acf(sim, plot = FALSE)
    t_acf[,i] <- sim$acf[2:11]
  }
  sim_acf[[j]] <- t_acf
}

MM <- Reduce("+", sim_acf) / length(sim_acf)


par(mfrow=c(3,3))
for(i in 1:9){
  plot(data_acf[i,],MM[i,], main = paste0("Lag-",i," Correlations - ", Data_Type),
       xlab = "True Data Correlation", ylab = "Mean Simulated Correlation", 
       pch = 19, ylim = c(-1,1), xlim = c(-1,1))
  abline(coef = c(0,1), lwd = 2)
}

#------------------------------------------------------------------------------#
#ACF on True Data
temp <- X
data_acf <- matrix(NA, ncol = ncol(X), nrow = 10)

for(i in 1:ncol(X)){
  t_acf <- acf(temp[,i], plot = FALSE)
  data_acf[,i] <- t_acf$acf[2:11]
}

#ACF on Simulated Data
sim_acf <- list()
for(j in 1:n_sim){
  t_acf <- matrix(NA, ncol = ncol(X), nrow = 10)
  for(i in 1:ncol(X)){
    sim <- as.data.frame(Sims[[j]])
    sim <- sim[,i]
    sim <- acf(sim, plot = FALSE)
    t_acf[,i] <- sim$acf[2:11]
  }
  sim_acf[[j]] <- t_acf
}


###Plotting 
n_sel <- sample(ncol(X), replace = FALSE)
quant <- c(0.05, 0.25, 0.75,0.95)

par(mfrow=c(2,2))
for(i in 1:8){
  plot(data_acf[i,],sim_acf[[1]][i,], main = paste0("Lag-",i,"  - ", Data_Type),
       xlab = "Reanalysis Correlation", ylab = "Simulation Correlation", 
       ylim = c(min(sim_acf[[1]][i,])-0.1, max(sim_acf[[1]][i,])+0.1),
       xlim = c(min(sim_acf[[1]][i,])-0.1, max(sim_acf[[1]][i,])+0.1),
       type = 'n', cex.lab = 1.75, cex.main = 2)
  abline(coef = c(0,1), lwd = 2)
  
  for(j in 1:length(n_sel)){
    st <- n_sel[j]
    tacf <- list()
    for(k in 1:n_sim){
      tacf[k] <- sim_acf[[k]][i,st]
    }
    tacf <- unlist(tacf)
    med <- median(tacf)
    qt <- quantile(tacf, prob = quant)
    
    
    segments(data_acf[i,st], qt[1], data_acf[i,st], qt[4], col = "red", lwd = 3)
    #segments(data_acf[i,st], qt[2], data_acf[i,st], qt[3], col = "blue", lwd = 2)
    points(data_acf[i,st],med, col ="blue", pch = 19)
    #legend('topleft', legend = c("Median Simulation Value", "Interquartile Range (IQR) ", 
    #                             "5 - 95 Percentile Range"), col = c('black','blue','red'),
    #       lty = c(0,1,1), pch = c(19,NA,NA), lwd = c(NA, 3, 3), cex = 1)
  }
}



#______________________________________________________________________________#
#Cumulative Probability Density for a single site.
tx <- X[,rand_count]
cdf_og <- ts_eval <- seq(min(tx)-sd(tx),max(tx)+sd(tx),max(tx)/20)
og_cdf <- ecdf(tx)
for(j in 1:length(ts_eval)) cdf_og[j] <- og_cdf(ts_eval[j])
sim_cdf <- matrix(NA, ncol = n_sim, nrow = length(ts_eval)) #Storing the Simulated CDF's
for(j in 1:n_sim){
  #Computing each CDF
  sim <- as.data.frame(Sims[[j]])
  sim <- sim[,rand_count]
  cdf_sim <- ecdf(sim)
  for(i in 1:length(ts_eval)) {sim_cdf[i,j] <- cdf_sim(ts_eval[i])}
}
#Getting the percentiles
lower_percentile <- apply(sim_cdf, 1, function(x) quantile(x, probs=.05))
upper_percentile <- apply(sim_cdf, 1, function(x) quantile(x, probs=.95))
median_percentile <- apply(sim_cdf, 1, median)
par(mfrow=c(1,1), mar = c(4,4,3,1))
plot(ts_eval, cdf_og, type='l',col='red',
     lwd = 2, main = paste0("True vs Simulated CDF for Single Site - ", Data_Type), 
     xlab = "x", ylab = "F(x)")
polygon(c(ts_eval,rev(ts_eval)),c(lower_percentile,rev(upper_percentile)),col="gray")
lines(ts_eval, median_percentile, lwd = 2)
lines(ts_eval, cdf_og, col='red', lwd = 2)
legend('topleft', legend = c("Median Simulation Value", "Reanalysis Data","5 - 95 Percentile Range"),
       lty = 1, col = c('black','red','grey'), lwd = 3, cex = 1)

#______________________________________________________________________________#
#PDF for a single site
tx <- X[,rand_count]
og_pdf <- density(tx, from = 0, to = (max(tx) + 1.5*sd(tx)))

sim_pdf <- matrix(NA, ncol = n_sim, nrow = length(og_pdf$x)) #Storing the Simulated CDF's
for(j in 1:n_sim){
  #Computing each CDF
  sim <- as.data.frame(Sims[[j]])
  sim <- sim[,rand_count]
  pdf_sim <- density(sim, from = 0, to = (max(tx) + 1.5*sd(tx)))
  sim_pdf[,j] <- pdf_sim$y
}

#Getting the percentiles
lower_percentile <- apply(sim_pdf, 1, function(x) quantile(x, probs=.05))
upper_percentile <- apply(sim_pdf, 1, function(x) quantile(x, probs=.95))
median_percentile <- apply(sim_pdf, 1, median)
par(mfrow=c(1,1), mar = c(6,6,4,1))
plot(og_pdf$x, og_pdf$y, type='l',col='red',
     lwd = 2, main = paste0("PDF at a single site - ", Data_Type), 
     xlab = "Watts/sq-m", ylab = "Density - f(x) ",
     ylim = c(0, 1.05*max(og_pdf$y)),
     cex.lab = 1.5)
polygon(c(og_pdf$x,rev(og_pdf$x)),c(lower_percentile,rev(upper_percentile)),col="gray")
lines(og_pdf$x, median_percentile, lwd = 2)
lines(og_pdf$x, og_pdf$y, col='red', lwd = 2)
legend('topright', legend = c("Median Simulation Value", "Reanalysis Data", "5 - 95 Percentile Range"),
       lty = 1, col = c('black','red','grey'), lwd = 3, cex = 1.25)




#------------------------------------------------------------------------------#
###Joint Quantiles
quant <- c(0.01,0.05,.1,.25,.50,.75,.90,.95,.99)

#Compute Quantiles for True Data. 
data_quant <- matrix(NA, ncol = ncol(X), nrow = length(quant))

for(i in 1:ncol(X)){
  t <- quantile(X[,i], prob = quant)
  data_quant[,i] <- t
}

#Computing Quantiles for Simulations. 
sim_quant <- list()
for(j in 1:n_sim){
  t_quant <- matrix(NA, ncol = ncol(X), nrow = length(quant))
  for(i in 1:ncol(X)){
    sim <- as.data.frame(Sims[[j]])
    sim <- sim[,i]
    sim <- quantile(sim, prob = quant)
    t_quant[,i] <- sim
  }
  sim_quant[[j]] <- t_quant
}

MM <- Reduce("+", sim_quant) / length(sim_quant)


par(mfrow=c(3,3))
for(i in 1:9){
  plot(data_quant[i,],MM[i,], main = paste0("Quantile - ",quant[i]*100,"th - ", Data_Type),
       xlab = "True Data", ylab = "Mean Simulated ", 
       pch=19, col='red')
  abline(coef = c(0,1), lwd = 2)
}

#------------------------------------------------------------------------------#
##Quantiles and Uncertainty
quant <- c(0.01,0.05,.1,.25,.75,.90,.95,.99)

#Compute Quantiles for True Data. 
data_quant <- matrix(NA, ncol = ncol(X), nrow = length(quant))

for(i in 1:ncol(X)){
  t <- quantile(X[,i], prob = quant)
  data_quant[,i] <- t
}

#Computing Quantiles for Simulations. 
sim_quant <- list()
for(j in 1:n_sim){
  t_quant <- matrix(NA, ncol = ncol(X), nrow = length(quant))
  for(i in 1:ncol(X)){
    sim <- as.data.frame(Sims[[j]])
    sim <- sim[,i]
    sim <- quantile(sim, prob = quant)
    t_quant[,i] <- sim
  }
  sim_quant[[j]] <- t_quant
}


n_sel <- sample(ncol(X), replace = FALSE)

par(mfrow=c(2,2), mar = c(6,6,6,3))
for(i in 1:8){
  plot(data_quant[i,],sim_quant[[1]][i,], 
       main = paste0(quant[i]*100,"th Quantile - ",Data_Type),
       xlab = "Reanalysis Data", ylab = "Simulation", type='n',
       ylim = c(min(0.95*sim_quant[[1]][i,]),1.05*max(sim_quant[[1]][i,])),
       xlim = c(min(0.95*sim_quant[[1]][i,]),1.05*max(sim_quant[[1]][i,])),,
       pch=19, col='red', cex.lab = 1.75, cex.main = 2)
  abline(coef = c(0,1), lwd = 2)
  
  for(j in 1:length(n_sel)){
    st <- n_sel[j]
    tq <- list()
    for(k in 1:n_sim){
      tq[k] <- sim_quant[[k]][i,st]
    }
    tq <- unlist(tq)
    med <- median(tq)
    quants <- c(0.05,0.25,0.75,0.90)
    qt <- quantile(tq, prob = quants)
    
    segments(data_quant[i,st], qt[1], data_quant[i,st], qt[4], col = "red", lwd = 3)
    #segments(data_quant[i,st], qt[2], data_quant[i,st], qt[3], col = "blue", lwd = 2)
    points(data_quant[i,st],med, col ="blue", pch = 19)
  }
  #legend('topleft', legend = c("Median Simulation Value", "Interquartile Range (IQR) ", 
  #                             "5 - 95 Percentile Range"), col = c('black','blue','red'),
  #      lty = c(0,1,1), pch = c(19,NA,NA), lwd = c(NA, 3, 3), cex = 1)
  
}



#______________________________________________________________________________#
#Spatial Loadings - PCA True Data
par(mfrow=c(3,2))
library(plotrix)
sel_sim <- sample(1:n_sim,1)
pcw <- prcomp(X, scale = TRUE, center = TRUE)
ju <- pcw$rotation[,1]
var <- (pcw$sdev^2)
var <- (var)/sum(var)

map('state', 
    region = c("Texas","OKlahoma", "New Mexico", "Louisiana", "Arkansas"))
points(grid_locs$lon-360,grid_locs$lat,pch=19,cex=1.5, 
       col=color.scale(ju,c(1,0.5,0),c(0,0.5,0),c(0,0,1),color.spec="rgb"))
title(paste0("PC-1 - ", Data_Type,"-", round(var[1]*100),"%"))

ju <- pcw$rotation[,2]
map('state', 
    region = c("Texas","OKlahoma", "New Mexico", "Louisiana", "Arkansas"))
points(grid_locs$lon-360,grid_locs$lat,pch=19,cex=1.5, 
       col=color.scale(ju,c(1,0.5,0),c(0,0.5,0),c(0,0,1),color.spec="rgb"))
title(paste0("PC-2 - ", Data_Type,"-", round(var[2]*100),"%"))

#PCA on individual simulations
pc1 <- pc2 <- list()
for(i in 1:n_sim){
  temp <- as.data.frame(Sims[[i]])
  temp <- scale(temp)
  pcw_temp <-  prcomp(temp, scale = TRUE, center = TRUE)
  tvar <- (pcw_temp$sdev^2)
  tvar <- (tvar)/sum(tvar)
  pc1[[i]] <- tvar[1]
  pc2[[i]] <- tvar[2]
}


pc1 <- unlist(pc1)
pc2 <- unlist(pc2)

#PCA on the ensemble average
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
ju <- pcw_temp$rotation[,1]

boxplot(pc1, use.cols = TRUE, 
        main = paste0('Simulated PC-1 Variance ', Data_Type ))
abline(h = var[[1]], col = 'red', lwd = 2)
abline(h = pc1_ens, col = 'blue', lwd = 2)
legend('topright', legend = c('Data', "Ensemble Mean"),
       col = c('red','blue'), lwd = 2)

boxplot(pc2, use.cols = TRUE, 
        main = paste0('Simulated PC-2 Variance ', Data_Type ))
abline(h = var[[2]], col = 'red', lwd = 2)
abline(h = pc2_ens, col = 'blue', lwd = 2)
legend('topright', legend = c('Data', "Ensemble Mean"),
       col = c('red','blue'), lwd = 2)


map('state', 
    region = c("Texas","OKlahoma", "New Mexico", "Louisiana", "Arkansas"))
points(grid_locs$lon-360,grid_locs$lat,pch=19,cex=1.5, 
       col=color.scale(ju,c(1,0.5,0),c(0,0.5,0),c(0,0,1),color.spec="rgb"))
title(paste0("Ensemble Mean PC-1 - ", Data_Type,"-", round(pc1_ens*100),"%"))

ju <- pcw$rotation[,2]
map('state', 
    region = c("Texas","OKlahoma", "New Mexico", "Louisiana", "Arkansas"))
points(grid_locs$lon-360,grid_locs$lat,pch=19,cex=1.5, 
       col=color.scale(ju,c(1,0.5,0),c(0,0.5,0),c(0,0,1),color.spec="rgb"))
title(paste0("Ensemble Mean PC-2 - ", Data_Type,"-", round(pc2_ens*100),"%"))


#______________________________________________________________________________#
#CDF for th entire field. 
tx <- rowSums(X)
cdf_og <- ts_eval <- seq(min(tx)-sd(tx),max(tx)+sd(tx),max(tx)/20)
og_cdf <- ecdf(tx)
for(j in 1:length(ts_eval)) cdf_og[j] <- og_cdf(ts_eval[j])
sim_cdf <- matrix(NA, ncol = n_sim, nrow = length(ts_eval)) #Storing the Simulated CDF's
for(j in 1:n_sim){
  #Computing each CDF
  sim <- as.data.frame(Sims[[j]])
  sim <- rowSums(sim)
  cdf_sim <- ecdf(sim)
  for(i in 1:length(ts_eval)) {sim_cdf[i,j] <- cdf_sim(ts_eval[i])}
}
#Getting the percentiles
lower_percentile <- apply(sim_cdf, 1, function(x) quantile(x, probs=.05))
upper_percentile <- apply(sim_cdf, 1, function(x) quantile(x, probs=.95))
median_percentile <- apply(sim_cdf, 1, median)
par(mfrow=c(1,1), mar = c(4,4,3,1))
plot(ts_eval, cdf_og, type='l',col='red',
     lwd = 2, main = paste0("Simulated CDF for Entire Field - ", Data_Type), 
     xlab = "Aggregated Variable (X)", ylab = "CDF  -  F(X) ")
polygon(c(ts_eval,rev(ts_eval)),c(lower_percentile,rev(upper_percentile)),col="gray")
lines(ts_eval, median_percentile, lwd = 2)
lines(ts_eval, cdf_og, col='red', lwd = 2)
legend('topleft', legend = c("Median Simulation Value", "True Data CDF", "5 - 95 Percentile Range"),
       lty = 1, col = c('black','red','grey'), lwd = 3, cex = 1)


#______________________________________________________________________________#
#PDF for th entire field.
#COnvert W/sq-m to MWhr/sq-m
t_fac <- 24*10000/(10^6)   #hrs*kms/Mega

tx <- rowSums(X)*t_fac
og_pdf <- density(tx, from = 0, to = (max(tx) + 0.5*sd(tx)))

sim_pdf <- matrix(NA, ncol = n_sim, nrow = length(og_pdf$x)) #Storing the Simulated CDF's
for(j in 1:n_sim){
  #Computing each CDF
  sim <- as.data.frame(Sims[[j]])
  sim <- rowSums(sim)*t_fac
  pdf_sim <- density(sim, from = 0, to = (max(tx) + 0.5*sd(tx)))
  sim_pdf[,j] <- pdf_sim$y
}

#Getting the percentiles
lower_percentile <- apply(sim_pdf, 1, function(x) quantile(x, probs=.05))
upper_percentile <- apply(sim_pdf, 1, function(x) quantile(x, probs=.95))
median_percentile <- apply(sim_pdf, 1, median)
par(mfrow=c(1,1), mar = c(6,6,4,1))
plot(og_pdf$x, og_pdf$y, type='l',col='red',
     lwd = 2, main = paste0("Daily Generation Potential across ERCOT - ", Data_Type), 
     xlab = "MWhr", ylab = "Density - f(x) ",
     ylim = c(0, 1.05*max(og_pdf$y)),
     cex.lab = 1.5)
polygon(c(og_pdf$x,rev(og_pdf$x)),c(lower_percentile,rev(upper_percentile)),col="gray")
lines(og_pdf$x, median_percentile, lwd = 2)
lines(og_pdf$x, og_pdf$y, col='red', lwd = 2)
legend('topright', legend = c("Median Simulation Value", "True Data CDF", "5 - 95 Percentile Range"),
       lty = 1, col = c('black','red','grey'), lwd = 3, cex = 1.25)




#______________________________________________________________________________#
#Simulation Visualization
n_sel <- sample(ncol(X), 1) #Site
sim <- sample(6,1)
sel_sim <- sample(n_sim, 1) #Simulation Run
par(mfrow=c(3,2))


for(i in 1:6){
  str <- sample(nrow(X)-100,1)
  if(i == sim){
    temp <- as.data.frame(Sims[[sel_sim]])
    plot(temp[str:(str+90),n_sel],type='l',
         ylab = ylab.text, xlab = "Simulated")
  } else{
    
    plot(X[str:(str+90),n_sel], type='l',
         ylab = ylab.text, xlab = "Observed")
  }
}





}
