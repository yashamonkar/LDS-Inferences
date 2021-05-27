#______________________________________________________________________________#
#Code for generating KSTS without any diagnostics. 



###Input
#1. Data
#2. Number of Neighbors
#3. Weights
#4. Simulation Length
#5. Start_Date
#6. Maximum Embedding
#7. Selected Lags
#8. Total Lags

###Output
#1. Field KSTS Simulations




#______________________________________________________________________________#
#Functions
#Moving Window Indices - Returns the day/hour index of interest.
close_ind <- function(curr,max,window){
  if((curr-window) < 1){
    
    indx <- c(tail(1:max,(1+abs(curr-window))),
              1:(curr+window))
    
  } else if((curr+window) > max) {
    
    indx <- c(((curr-window):max),
              (1:((curr+window)-max)))
    
  } else {
    
    indx <- (curr-window):(curr+window)
    
  }
  return(indx)
}

#KNN - non parametric nearest neighbors resampling based on training and test data.
#Returns Index of Nearest neighbors
knn.sim.index <- function(x,xtest,nneib,weights){
  
  d <- matrix(NA, ncol=ncol(xtest), nrow=nrow(x))
  #Compute Distances
  for(i in 1:ncol(x))
    d[,i] <- weights[i]*(x[,i]-xtest[,i])^2
  d <- rowSums(d)
  sorted.data = sort.int(d,method="quick",index.return=TRUE)
  
  sorted.neib = matrix(sorted.data$ix[1:nneib])
  
  out = list(yknn=sorted.neib
  )
}


#KSTS Simulator
knngrids <- function(Fld,
                     nneib,weights, #KNN Parameter
                     sim_length, #Simulation length
                     start_date, day_mv,  #Date Seasonality Parameters
                     max_embd, sel_lags, n_lags) #Embedding Parameters
{ library(dplyr)
  library(lubridate)
  
  #Hyper-parameters
  N_valid <- nrow(Fld)
  ngrids <- ncol(Fld)
  
  #Day Indices
  st_date <- as.Date(start_date, format="%m-%d-%Y")
  time_stamps <- seq(st_date, by = "day", length.out = N_valid)
  day_index <- yday(time_stamps)
  
  #Simulation Day Indices
  time_stamps <- seq(st_date, by = "day", length.out = sim_length)
  sim_day_index <- yday(time_stamps)
  
  #Setting up Storage for Simulations
  Xnew = matrix(NA,nr=sim_length,nc=ngrids)
  for(i in 1:max_embd){
    Xnew[i,] = jitter(Fld[i,]) 
  }
  
  #Creating the feature Vector
  X = array(NA, c(N_valid-max_embd,n_lags,ngrids))
  Y = array(NA, c(N_valid-max_embd,1,ngrids))
  for (j in 1:ngrids) 
  {
    #Get Lagged Structure Upto Max Embedding
    x_fld <- embed(Fld[,j],max_embd+1) 
    
    X[,,j] <- x_fld[,sel_lags+1]
    Y[,,j] <- x_fld[,1]
  }
  
  
  pb = txtProgressBar(min = (max_embd+1), max = sim_length, initial = 1) 
  for (i in (max_embd+1):sim_length){
    
    setTxtProgressBar(pb,i)
    
    nn_index <- list() #Store all the nearest neighbours.
    day <- sim_day_index[i]
    sel_days <- close_ind(day,max = 366, window = day_mv)
    
    #Subset to the moving window
    indx <- day_index
    days <- tail(indx, -max_embd)
    X_t <- X[days %in% sel_days,,]
    Y_t <- Y[days %in% sel_days,,]
    
    for (j in 1:ngrids){
      
      #Setting the Test Parameters
      sel_pars <- i-sel_lags
      xtest <- matrix(Xnew[sel_pars,j], nrow =1)
      
      #Running the Simulator
      nn_index[[j]] <- knn.sim.index(X_t[,,j],xtest,nneib,weights)
    }
    
    #Computing the Resampling Probability
    un_index <- unique(unlist(nn_index))
    un_prob <- rep(NA, length(un_index))
    
    nn_index  <-  matrix(unlist(nn_index), nrow=nneib)
    
    for(k in 1:length(un_index)){
      temp <- which(nn_index == un_index[k]) %% nneib
      for(l in 1:length(temp)){
        if(temp[l]==0) {
          temp[l] = 1/nneib
        } else {
          temp[l]=1/temp[l]
        }
      }
      un_prob[k] <- sum(temp)
    }
    #Computing the Resampling Probablites 
    pj <- data.frame(cbind(un_prob,un_index)) 
    thresh <- apply(pj, 2, FUN = function(x) tail(sort(x), nneib+1)[1])[1]
    pj <- pj %>% filter(un_prob > thresh)  
    pj$un_prob <- pj$un_prob/sum(pj$un_prob)
    ns <- sample(pj$un_index,1,prob = pj$un_prob)
    
    Xnew[i,] <- Y_t[ns,]
    
  }
  
  out = list(Xnew=Xnew)
}