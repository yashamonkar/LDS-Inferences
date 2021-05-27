#______________________________________________________________________________#
#Code to get the PC-ARIMA Simulations


###Baseline PC-ARIMA Simulations###

###Method
#Step 1:- PCA on Domain
#Step 2:- ARIMA on leading PCs
#Step 3:- Gaussian Noise for other PCs
#Step 4:- Re-transform to other data. 

###Input
#1. Data_Field:- True Data to be Simulated
#2. Field_Name:- Name of the field. e.g. "Wind Power"
#3. PCsel :- Number of PC's which have ARIMA
#4. NSims :- Number of Ensemble Simulations 50

###Output
#1. List of List of Simulations 

#______________________________________________________________________________#
Get_PCAR_Simulations <- function(Data_Field, Field_Name, PCsel, NSims){
  
  library(MASS)
  library(forecast)
  
  #Principal Component Analysis
  pcw <- prcomp(Data_Field, scale = TRUE, center = TRUE)
  var <- pcw$sdev^2
  plot(1:10, cumsum(var[1:10]), pch =19, main = "Principal Component Analysis",
       xlab = "Leading PCs", ylab = "Variance")
  pcs <- pcw$x[,1:PCsel]
  
  #Re-Transforming the PCs
  pc_sd <- apply(pcs,2, sd)
  pcs_trans <- scale(pcs)
  
  ###Fitting ARIMA Processes 
  ar.fit <- list()
  pb = txtProgressBar(min = 1, max = ncol(pcs_trans), initial = 1) 
  for(i in 1:ncol(pcs_trans)){
    setTxtProgressBar(pb,i)
    ar.fit[[i]] <- auto.arima(pcs_trans[,i], seasonal = TRUE, allowmean = FALSE,
                              allowdrift = FALSE, start.p = 1, start.q=1)
  }
  
  
  ###Generate Emsemble Simulations
  Ensem_Sims <- list()
  for(j in 1:NSims){
    
    #Set-up storage for that run
    Sims <- matrix(NA, ncol = ncol(Data_Field), nrow = nrow(Data_Field)) 
    
    #Adding the ARIMA parts
    for(i in 1:ncol(pcs_trans)){
      #Get the Arima Order 
      ord = arimaorder(ar.fit[[i]])
      
      if(ord[1]==0 && ord[3] > 0){
        ma = ar.fit[[i]]$coef
        ts.sim <- arima.sim(list(order = ord, ma = ma),
                            n = nrow(Data_Field), 
                            sd = sqrt(ar.fit[[i]]$sigma2))
        
      } else if(ord[3]==0 && ord[1] > 0){
        ar = ar.fit[[i]]$coef
        ts.sim <- arima.sim(list(order = ord, ar = ar),
                            n = nrow(Data_Field), 
                            sd = sqrt(ar.fit[[i]]$sigma2))
        
      } else if(ord[3]==0 && ord[1]==0){
        
        ts.sim <- arima.sim(list(order = ord),
                            n = nrow(Data_Field), 
                            sd = sqrt(ar.fit[[i]]$sigma2))
        
      } else{
        ar = ar.fit[[i]]$coef[1:ord[1]]
        ma = ar.fit[[i]]$coef[(ord[1]+1):(ord[1]+ord[3])]
        ts.sim <- arima.sim(list(order = ord, ar = ar, ma = ma),
                            n = nrow(Data_Field), 
                            sd = sqrt(ar.fit[[i]]$sigma2))
      }
      
      Sims[,i] <- tail(ts.sim,nrow(Sims))
      
      #Renormalize the Data
      Sims[,i] <- Sims[,i]*pc_sd[i]
    }
    
    #Adding the Rest of the PCs as White Noise
    for(i in (PCsel+1):ncol(Data_Field)){
      sd_t <- sd(pcw$x[,i])
      Sims[,i] <- rnorm(nrow(Data_Field),0,sd_t)
    }
    
    #Re-Converting back to OG Spave from PC-Space
    Sims <- t(t(Sims %*% t(pcw$rotation)) * pcw$scale + pcw$center)
  
    Ensem_Sims[[j]] <- Sims
    
  }
  
  #Retun the Simulations
  return(Ensem_Sims)
  
}