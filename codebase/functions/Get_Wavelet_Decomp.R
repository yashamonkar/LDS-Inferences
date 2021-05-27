#Function to give out the 
#1. Global Wavelet Signal for a Time Series
#2. Portions which exceed a threshold over a confidence interval. 
#3. Variance explained by each period. 

#Inputs
#1. Signal itself. 
#2. Confidence Intervals. 
#3. If the variance should be a rolling mean.

#Outputs. 
#1. Plot with signal over the threshold. 
#2. Written output of signals over theshold. 
#3. Variance explanied by each signal. 

Get_Wavelet_Decomp <- function(signal, nam, signif){
  
  ### Wavelet Function
  
  #WAVELET  1D Wavelet transform with optional singificance testing
  #
  #   [WAVE,PERIOD,SCALE,COI] = wavelet(Y,DT,PAD,DJ,S0,J1,MOTHER,PARAM)
  #
  #   Computes the wavelet transform of the vector Y (length N),
  #   with sampling rate DT.
  #
  #   By default, the Morlet wavelet (k0=6) is used.
  #   The wavelet basis is normalized to have total energy=1 at all scales.
  #
  #
  # INPUTS:
  #
  #    Y = the time series of length N.
  #    DT = amount of time between each Y value, i.e. the sampling time.
  #
  # OUTPUTS:
  #
  #    WAVE is the WAVELET transform of Y. This is a complex array
  #    of dimensions (N,J1+1). FLOAT(WAVE) gives the WAVELET amplitude,
  #    ATAN(IMAGINARY(WAVE),FLOAT(WAVE) gives the WAVELET phase.
  #    The WAVELET power spectrum is ABS(WAVE)^2.
  #    Its units are sigma^2 (the time series variance).
  #
  #
  # OPTIONAL INPUTS:
  # 
  # *** Note *** setting any of the following to -1 will cause the default
  #               value to be used.
  #
  #    PAD = if set to 1 (default is 0), pad time series with enough zeroes to get
  #         N up to the next higher power of 2. This prevents wraparound
  #         from the end of the time series to the beginning, and also
  #         speeds up the FFT's used to do the wavelet transform.
  #         This will not eliminate all edge effects (see COI below).
  #
  #    DJ = the spacing between discrete scales. Default is 0.25.
  #         A smaller # will give better scale resolution, but be slower to plot.
  #
  #    S0 = the smallest scale of the wavelet.  Default is 2*DT.
  #
  #    J1 = the # of scales minus one. Scales range from S0 up to S0*2^(J1*DJ),
  #        to give a total of (J1+1) scales. Default is J1 = (LOG2(N DT/S0))/DJ.
  #
  #    MOTHER = the mother wavelet function.
  #             The choices are 'MORLET', 'PAUL', or 'DOG'
  #
  #    PARAM = the mother wavelet parameter.
  #            For 'MORLET' this is k0 (wavenumber), default is 6.
  #            For 'PAUL' this is m (order), default is 4.
  #            For 'DOG' this is m (m-th derivative), default is 2.
  #
  #
  # OPTIONAL OUTPUTS:
  #
  #    PERIOD = the vector of "Fourier" periods (in time units) that corresponds
  #           to the SCALEs.
  #
  #    SCALE = the vector of scale indices, given by S0*2^(j*DJ), j=0...J1
  #            where J1+1 is the total # of scales.
  #
  #    COI = if specified, then return the Cone-of-Influence, which is a vector
  #        of N points that contains the maximum period of useful information
  #        at that particular time.
  #        Periods greater than this are subject to edge effects.
  #        This can be used to plot COI lines on a contour plot by doing:
  #
  #              contour(time,log(period),log(power))
  #              plot(time,log(coi),'k')
  #
  #----------------------------------------------------------------------------
  #   Copyright (C) 1995-2004, Christopher Torrence and Gilbert P. Compo
  #
  #   This software may be used, copied, or redistributed as long as it is not
  #   sold and this copyright notice is reproduced on each copy made. This
  #   routine is provided as is without any express or implied warranties
  #   whatsoever.
  #
  # Notice: Please acknowledge the use of the above software in any publications:
  #    ``Wavelet software was provided by C. Torrence and G. Compo,
  #      and is available at URL: http://paos.colorado.edu/research/wavelets/''.
  #
  # Reference: Torrence, C. and G. P. Compo, 1998: A Practical Guide to
  #            Wavelet Analysis. <I>Bull. Amer. Meteor. Soc.</I>, 79, 61-78.
  #
  # Please send a copy of such publications to either C. Torrence or G. Compo:
  #  Dr. Christopher Torrence               Dr. Gilbert P. Compo
  #  Research Systems, Inc.                 Climate Diagnostics Center
  #  4990 Pearl East Circle                 325 Broadway R/CDC1
  #  Boulder, CO 80301, USA                 Boulder, CO 80305-3328, USA
  #  E-mail: chris[AT]rsinc[DOT]com         E-mail: compo[AT]colorado[DOT]edu
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
  
  #Reading the Data
  input_ts <- signal
  
  #Computing the wavelet signal. 
  wlt <- wavelet(input_ts)
  Cw=CI(signif,input_ts,"w")
  C=CI(signif,input_ts,"r")
  
  
  #Plotting the Wavelet Spectrum
  plot(wlt$period,wlt$p.avg,xlab="Period",ylab="Variance",
       main=paste0("Global Wavelet Spectrum - ",nam))
  lines(wlt$period,wlt$p.avg, lwd = 3)
  lines(wlt$period,Cw$sig, lwd = 3)
  lines(wlt$period,C$sig,col="red", lwd = 3)
  
  
  #Reconstruct the Signal
  Cd <- 0.776;psi0 <- pi^(-.25);dj=0.025 #From the Torrence and Compo
  reconst <- matrix(NA, ncol = length(input_ts), nrow = length(wlt$scale))
  for(j in 1:ncol(reconst)) {reconst[,j] <- Re(wlt$wave[,j])/(wlt$scale^0.5)}
  p_reconst <- colSums(reconst)*dj/(Cd*psi0)
  
  #Reconstruction Check
  plot(input_ts, type='l',
       main = "Reconstruction Check", lwd = 3, ylab = 'Signal')
  lines(1:length(input_ts), p_reconst, col='red', lwd = 3)
  legend('topright', c("Signal","Reconstructed"), col = c('black','red'),
         lty = 1, lwd= 3)
  
  if(length(input_ts) > 1000){
  #Reconstruction Check
  plot(input_ts[1:100], type='l',
       main = "Reconstruction Check", lwd = 3, ylab = 'Signal')
  lines(1:100, p_reconst[1:100], col='red', lwd = 3)
  legend('topright', c("Signal","Reconstructed"), col = c('black','red'),
         lty = 1, lwd= 3)}
  
  
  #Getting the Significant Scales
  wlt_dataset <- data.frame(Power = wlt$p.avg,
                            Red_Noise = C$sig)
  wlt_dataset$Signif <- wlt_dataset$Power-wlt_dataset$Red_Noise
  wlt_dataset$Signif <- ifelse(wlt_dataset$Signif < 0, 0, wlt_dataset$Signif)
  wlt_dataset$Signif <- ifelse(wlt_dataset$Signif > 0, 1, wlt_dataset$Signif)
  print(paste0("The wavelet periods over the threshold are ", wlt$period[which(wlt_dataset$Signif==1)]))
  
  #Associating Scales with each cluster
  clust <- which(wlt_dataset$Signif==1)
  n_clust <- length(which(diff(clust)>1))+1
  clust_list <- list(NA)
  if (n_clust > 1) { st <- 1
  for(j in 1:n_clust) {
    en <- c(which(diff(clust)>1),length(clust))
    en <- en[j]
    clust_list[[j]] <- clust[st:en] 
    st <- en+1
  } 
  }
  if(n_clust < 2) { clust_list <- list(clust)}
  n_clust <- length(rapply(clust_list, length))
  
  #Variance associated with Significant and Non-Significant Clusters
  if(length(clust) > 0){
  t <- reconst[unlist(clust_list),]
  t_reconst <- colSums(t)*dj/(Cd*psi0)
  t <- reconst[-unlist(clust_list),]
  s_reconst <- colSums(t)*dj/(Cd*psi0)  #Non-Sig
  sig_info <- 100*var(t_reconst)/var(p_reconst)

  #Variance Check
  plot(input_ts, type='l',
       main = paste0("Variance Check ", round(sig_info,2), " %"), 
       lwd = 3, ylab = 'Signal')
  lines(1:length(input_ts), t_reconst, col='red', lwd = 3)
  lines(1:length(input_ts), s_reconst, col='blue', lwd = 3)
  legend('topright', c("Signal","Signif", "Non-Signif"), col = c('black','red','blue'),
         lty = 1, lwd= 3)
  }
  
  
  #Variance explained at each scale. 
  scale_ts <- matrix(NA, ncol = length(input_ts), nrow = length(wlt$scale))
  for(i in 1:nrow(scale_ts)){
    scale_ts[i,] <-  reconst[i,]*dj/(Cd*psi0)
  }
  scale_var <- apply(scale_ts, 2, var)
  scale_var <- 2*scale_var/max(scale_var)
  
  
  #Plotting the Wavelet Spectrum
  plot(wlt$period,wlt$p.avg,xlab="Period",ylab="Variance",
       main=paste0("Scalewise Variance Wavelet Spectrum - ",nam), type = 'l')
  points(wlt$period,wlt$p.avg, pch = 19, cex = scale_var)
  lines(wlt$period,Cw$sig, lwd = 3)
  lines(wlt$period,C$sig,col="red", lwd = 3)
  
  
}