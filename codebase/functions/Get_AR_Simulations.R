#______________________________________________________________________________#

###Baseline-Threshold Model###

#Method
#Step 1:- AR on Sites. 
#Step 2:- Simulate fitted AR model and to get simulations. 


#Input
#1. Data Field
#2. Data Field Name
#3. Number of Simulations
#4. ic - Information Criteria e.g. AIC or BIC. 

#Output
#1. Ensemble Simulations. 

Get_AR_Simulations <- function(Data_Field, Field_Name, NSims, ic){
  
  library(MASS)
  library(forecast)

  
  #Scale the transformed Data. 
  site_means <- apply(Data_Field, 2, mean)
  site_sd <- apply(Data_Field,2, sd)
  Data_Field_trans <- scale(Data_Field)
  
  
  ###Fitting ARIMA Processes 
  ar.fit <- list()
  pb = txtProgressBar(min = 1, max = ncol(Data_Field_trans), initial = 1) 
  for(i in 1:ncol(Data_Field_trans)){
    setTxtProgressBar(pb,i)
    ar.fit[[i]] <- auto.arima(Data_Field_trans[,i], seasonal = TRUE, allowmean = FALSE,
                              allowdrift = FALSE, start.p = 1, start.q=1, ic = ic)
  }
  
  ###Visualization of the fitted ARMA Process
  plot(0:5,0:5, type='n', xlab = "AR Model", ylab = "MA Model",
       main = "Site Fitted AR and MA Models")
  for(i in 1:length(ar.fit)){
    ord = arimaorder(ar.fit[[i]])
    points(jitter(ord[1]),jitter(ord[3]),pch=19)
  }
  
  
  ###Generate Emsemble Simulations
  Ensem_Sims <- list()
  for(j in 1:NSims){
    Sims <- matrix(NA, ncol = ncol(Data_Field), nrow = nrow(Data_Field))
    
    for(i in 1:ncol(Data_Field_trans)){
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
      Sims[,i] <- Sims[,i]*site_sd[i] + site_means[i]
    }
    
    Ensem_Sims[[j]] <- Sims
  }
  
  #Retun the Simulations
  return(Ensem_Sims)
}