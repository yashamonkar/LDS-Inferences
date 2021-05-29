#______________________________________________________________________________#
#Code for PC-KNN Simulations#


###Baseline PC-KNN Simulations###

###Method
#Step 1:- PCA on Domain
#Step 2:- KNN on leading PCs
#Step 3:- Gaussian Noise for other PCs
#Step 4:- Re-transform to other data. 

###Input
#1. Data_Field:- True Data to be Simulated
#2. Field_Name:- Name of the field. e.g. "Wind Power"
#3. PCsel :- Number of PC's which have KNN
#4. nneib :- Number of Nearest Neighbors
#4. NSims :- Number of Ensemble Simulations 50

###Output
#1. List of List of Simulations 

#______________________________________________________________________________#
###Functions
#KNN - non parametric nearest neighbor resampling based on training and test data.
knn.sim <- function(y,x,xtest,nneib,pj,nsample){
  
  #Compute Distances
  for(i in 1:ncol(x))
    d <- (x[,i]-xtest[,i])^2
  
  sorted.data = sort.int(d,method="quick",index.return=TRUE)
  
  sorted.neib = matrix(y[sorted.data$ix[1:nneib],])
  #  sorted.dates = sorted.data$ix[1:nneib]
  
  ynp = sample(1:nneib,nsample,replace=T,prob=pj[1:(nneib)])
  
  yknn = sorted.neib[ynp,]
  #  yknndates = sorted.dates[ynp]
  
  out = list(yknn=yknn
             #yknndates=yknndates
  )
}


#______________________________________________________________________________#
Get_PCKNN_Simulations <- function(Data_Field, Field_Name, PCsel, nneib, NSims){
  
  
  #Principal Component Analysis
  pcw <- prcomp(Data_Field, scale = TRUE, center = TRUE)
  var <- pcw$sdev^2
  plot(1:10, cumsum(var[1:10]), pch =19, main = "Principal Component Analysis",
       xlab = "Leading PCs", ylab = "Variance")
  pcs <- pcw$x[,1:PCsel]
 
  
  ###Generate Emsemble Simulations
  Ensem_Sims <- list()
  pb = txtProgressBar(min = 1, max = NSims, initial = 1) 
  for(j in 1:NSims){
    setTxtProgressBar(pb,j)
    
    #Set-up storage for that run
    Sims <- matrix(NA, ncol = ncol(Data_Field), nrow = nrow(Data_Field))
    
    for(i in 1:PCsel) {
      
      #Generate the feature space
      Y <- matrix(embed(pcs[,i],2)[,1], ncol=1)
      X <- matrix(embed(pcs[,i],2)[,2], ncol=1)
        
      #Initialize
      Sims[1,i] <- jitter(pcs[1,i]) 
      sj = matrix(NA,nneib,1);for (jk in 1:nneib) { sj[jk,1] = 1/jk }
      pj = matrix(NA,nneib,1);for (jk in 1:nneib) { pj[jk,1] = sj[jk,1]/sum(sj) }
        
      #Running the KNN Algorithm
      for(k in 2:nrow(Sims)){
        xtest <- matrix(Sims[(k-1),i], ncol = 1)
        Sims[k,i] <- knn.sim(Y,X,xtest,nneib,pj,1)$yknn
      }
    }

    boot_row <- sample(1:nrow(Sims),nrow(Sims),replace = TRUE)
    Sims[,(PCsel+1):ncol(Sims)] <- pcw$x[boot_row,(PCsel+1):ncol(pcw$x)]
      
    
    #Re-Converting back to OG Spave from PC-Space
    Sims <- t(t(Sims %*% t(pcw$rotation)) * pcw$scale + pcw$center)
    
    Ensem_Sims[[j]] <- Sims
  

  }
  
  #Retun the Simulations
  return(Ensem_Sims)
  
}