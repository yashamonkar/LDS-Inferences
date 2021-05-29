#______________________________________________________________________________#

#Function to get Visualize the Simulation Skill on the PC-Wavelet Domain. 
#Method. 
#1. Run PCA on the Domain
#2. Select the number of PC's to focus on. 
#3. Run Wavelets on each PC.
#4. Visualize the spread. 

#Inputs. 
#1. Original/True Data. - Dataframe 
#2. Simulations - List of Lists with each simulation as Dataframe. 
#3. Field Name

#Output:
#1. Plot showing the true and simulated Global Wavelet Spectrum.


Get_PCWavelet_Sim_Skill <- function(Data_Field, Field_Name, Sims, PCs){
  
  #Import Dependencies
  
  
  #Get Hyper-Parameters
  nsims <- length(Sims)
  
  #Import Functions
  ####Functions_Needed###########
  #----------------------------------------------------------------------------
  wavelet=function(Y,dj=0.025){
    
    #Y is time series to be analyzed
    DT=1# is timestep for annual data, 1
    pad=1
    #dj=0.025
    param=6
    #pad data ? 0=F, 1=T
    #dj= spacing between discrete scales (.025)
    #param = wavenumber (6)
    
    s0=2*DT
    
    n1 = length(Y)
    J1=floor((log2(n1*DT/s0))/dj)
    
    
    #....construct time series to analyze, pad if necessary
    x = Y - mean(Y)
    
    
    if (pad == 1){
      base2 = trunc(log(n1)/log(2) + 0.4999)   # power of 2 nearest to N
      x = c(x, rep(0, 2^(base2 + 1) - n1))
    }
    n = length(x)
    
    #....construct wavenumber array used in transform [Eqn(5)]
    k = (1:trunc(n/2))
    k = k*((2*pi)/(n*DT))
    k = c(0, k, -rev(k[1:floor((n-1)/2)]))
    
    #....compute FFT of the (padded) time series
    f = fft(x)    # [Eqn(3)]
    
    #....construct SCALE array & empty PERIOD & WAVE arrays
    scale = s0*2^((0:J1)*dj)
    period = scale;
    wave = matrix(data=0, ncol=n, nrow=J1+1)  # define the wavelet array
    wave = as.complex(wave)  # make it complex
    wave=matrix(data=wave, ncol=n, nrow=J1+1)
    
    # loop through all scales and compute transform
    for(a1 in 1:(J1+1)){
      scl=scale[a1]		
      
      nn = length(k);
      k0 = param
      expnt = -(scl*k - k0)^2/(2*(k > 0))
      norm = sqrt(scl*k[2])*(pi^(-0.25))*sqrt(nn)    # total energy=N   [Eqn(7)]
      daughter = norm*exp(expnt)
      daughter = daughter*(k > 0)    # Heaviside step function
      fourier_factor = (4*pi)/(k0 + sqrt(2 + k0^2)) # Scale-->Fourier [Sec.3h]
      coi = fourier_factor/sqrt(2)                  # Cone-of-influence [Sec.3g]
      dofmin = 2                                   # Degrees of freedom
      
      out <- list(daughter=daughter, fourier_factor=fourier_factor,coi=coi,dofmin=dofmin)
      
      daughter=out$daughter
      fourier_factor=out$fourier_factor
      coi=out$coi
      dofmin=out$dofmin	
      wave[a1,] = fft((f*daughter), inverse = TRUE)/(length(f*daughter))  # wavelet transform[Eqn(4)]
    }
    
    period = fourier_factor*scale
    
    coi = coi*c(1:(floor(n1 + 1)/2), rev(1:floor(n1/2))) * DT
    
    wave = wave[,1:n1]  # get rid of padding before returning
    power=abs(wave)^2
    ncol=length(power[1,])
    nrow=length(scale)
    avg.power=apply(power,1,mean)
    result=list(wave=wave, period=period, scale=scale, power=power, coi=coi,nc=ncol,nr=nrow,p.avg=avg.power)
    return(result)
  }
  
  ### Confidence level function
  CI=function(conf, dat,type){
    
    #### enter confidence as decimal 0-1
    #### two types of tests available 1) red noise enter: "r" , white noise enter: "w"
    # requires the wavelet function
    
    na=length(dat)
    wlt=wavelet(dat)
    
    if(type=="r"){
      
      zz=arima(dat/10^10, order = c(1, 0, 0))
      alpha=zz$coef[1]
      print(alpha)
    } else{
      alpha=0
    }
    
    ps=wlt$period
    LP= length(ps)
    
    freq = 1/ps
    
    CI=1:LP    ## confidence interval
    
    for(i in 1:LP){
      
      P=(1-(alpha^2))/(1+(alpha^2)-(2*alpha*cos(2*pi*freq[i])))    # descrete fourier power spectrum page 69 [qn 16] ( torrence and compo)... null hypothesis test
      df=2*sqrt(1+((na/(2.32*ps[i]))^2))
      CI[i] =P*(qchisq(conf, df)/df)*var(dat)          #divide by 2 removes power of 2.....for mean no chi dist.[ eqn 17]
    }
    
    
    list(sig=CI)
    
  }
  
  #Loop in for the Required PCs
  for(pc in 1:PCs){
  
    #Principal Component Analysis
    pcw <- prcomp(Data_Field, scale = TRUE, center = TRUE)
    tx <- pcw$x[,pc]
    wlt_og <- wavelet(tx)
  
    #Wavelets for Simulations
    sim_pow <- matrix(NA, nrow = nsims, ncol = length(wlt_og$p.avg))
      for(j in 1:nsims) {
        temp <- as.data.frame(Sims[[j]])
        pcw <- prcomp(temp, scale = TRUE, center = TRUE)
        wlt <- wavelet(pcw$x[,pc])
        sim_pow[j,] <- wlt$p.avg
      }
    
      Cw=CI(0.9,tx,"w");
      C=CI(0.9,tx,"r");
      par(mar=c(6,6,6,2),mfrow = c(1,1))
      plot(wlt_og$period,wlt_og$p.avg,xlim=c(0,length(tx)*0.25),
          main=paste0("Global Wavelet Spectrum PC - ",pc," ", Field_Name),
          xlab="Period (Days)",ylab="Variance", col ='red',
          cex.lab = 1.5, cex.main = 1.5)
    
      lower_percentile <- apply(sim_pow, 2, function(x) quantile(x, probs=.05))
      upper_percentile <- apply(sim_pow, 2, function(x) quantile(x, probs=.95))
      median_percentile <- apply(sim_pow, 2, median)
      low <- cbind(lower_percentile,wlt_og$period)
      high <- cbind(upper_percentile,wlt_og$period)
    
      polygon(c(wlt_og$period,rev(wlt_og$period)),
              c(lower_percentile,rev(upper_percentile)),col="gray")
      lines(wlt_og$period, lower_percentile)
      lines(wlt_og$period, upper_percentile)
      lines(wlt_og$period,C$sig, lty=2, col ='red')
      lines(wlt_og$period,Cw$sig, lty=2, col ='black')
      lines(wlt_og$period,wlt_og$p.avg, col ='red', lwd = 2)
      points(wlt_og$period,wlt_og$p.avg, col ='red')
      lines(wlt_og$period, median_percentile, lwd = 2)
      legend('topright', 
             legend = c("Median Simulation", "Data","5 - 95 Percentile",
                        "Red Noise Signif", "White Noise Signif"),
           lty = c(1,1,1,2,2), lwd = c(3,3,3,1,1), col = c('black','red','grey', 'red', 'black'),
           cex = 1.5)
      
    }


  
  
  
}
