rm(list=ls())

noisetype <- 'ar' #change manually
phivalue <- 0.5 # ar correlation

timelength_sets <- c(1006,10054)
bd_sets <- c('33','31')
noiselevel_sets <- c(0.2,0.5)


tlag_sets <- c(0,2) 

henon_map <- function(b, d, mu, n) {
  x <- numeric()
  y <- numeric()
  x[1]<- runif(1)
  x[2] <- runif(1)
  y[1]<- runif(1)
  y[2] <- runif(1)
  for (k in 2:n) {
    x[k+1] <- 1.4+b*x[k-1]-(x[k])^2
  }
  
  for (k in 2:n) {
    y[k+1] <- 1.4+d*y[k-1]-(mu*x[k]+(1-mu)*y[k])*y[k]
  }
  return (list(x=x,y=y))
}

##############################################################################################
for (tlag in tlag_sets) {
  for (timelength in timelength_sets) {
    for (bd in bd_sets) {
      for (noiselevel in noiselevel_sets) {
        filtdir <- function(x,y) {
          #Parameters
          pi2 <- 2*pi
          npt <- length(x)
          ncut <- 1
          #Hilbert transform and computation of phases
          ht <- HilbertTransform(x)
          phi1 <- Im(log((ht/abs(ht))))
          phi1 <- phi1[ncut:(npt-ncut)]
          ht <- HilbertTransform(y)
          phi2 <- Im(log((ht/abs(ht))))
          phi2 <- phi2[ncut:(npt-ncut)]
          #Unwrap phases
          phi1 <- unwrap(phi1);phi2 <- unwrap(phi2)
          #%Check synchronization index
          synchroindex=(mean(cos(phi1-phi2)))^2 + (mean(sin(phi1-phi2)))^2
          
          return(synchroindex)
        }
        
        divergence_measures <- function(map1, map2, a_R, a_JR) {
          
          K12<- 0; K21 <- 0; K <- 0
          Da <- numeric(length(a_R))
          J <- 0
          JR <- numeric(length(a_JR))
          MI <- 0
          #MSE
          prod <- mean(colMeans(map1*map2))
          dif <- prod/(mean(colMeans(map1))*mean(colMeans(map2)))
          
          #KL 1->2
          K12 <- sum(map1*log(map1/map2))
          #K21
          K21 <- sum(map2*log(map2/map1))
          #average KL
          #if directly translated from MATLAB CODE: K <- sum(0.5*(map1*log(map1/map2)+map2*log(map2/map1)))
          #by the equation in the paper:
          K <-  (K21+K12)/2
          #Renyi
          for (cc in 1:length(a_R)) {
            Da[cc] <- 1/(a_R[cc]-1)*log(sum((map1^(a_R[cc]))*(map2^(1 - a_R[cc]))))
          }
          #Jensen-Shannon
          map <- 0.5*(map1 + map2)
          J <- sum(-map*log(map) + 0.5*(map1*log(map1)+map2*log(map2)))
          #Jensen-Renyi
          map <- sqrt(map1*map2)
          pmap <- pmap1 <- pmap2 <- c()
          for (cc in 1:length(a_JR)) {
            pmap[cc] <- sum(map^(a_JR[cc]))
            pmap1[cc] <- sum(map1^(a_JR[cc]))
            pmap2[cc] <- sum(map2^(a_JR[cc]))
            JR[cc] <- 1/(1-a_JR[cc])*(log(pmap[cc])-0.5*(log(pmap1[cc])+log(pmap2[cc])))
          }
          #Mutual information
          MI <- -sum(map*log(map))
          #output as a list
          out <- list(dif=dif, K12=K12, K21=K21, K=K, Da=Da, J=J, JR=JR, MI=MI)
          return(out)
        }
        ##########################################################
        #cwt
        library(wavScalogram)
        cwt_signal <- function(Sig, Cwt, Min, Max, FreqEch, Senuil, frqsmp, Aff, Fc, dt=1) {
          minscal <- 1/(Fc * FreqEch / Min)
          maxscal <- 1/(Fc * FreqEch / Max)
          step <- frqsmp*(maxscal - minscal)/(Max - Min)
          #scales
          s <- 1/seq(maxscal, minscal,-step)
          ccfs <- cwt_wst(signal=Sig, wname = Cwt, dt=dt, scales = s, makefigure=FALSE, powerscales=FALSE)$coefs
          ccfs <- t(ccfs)
          ccfs <- ccfs[rev(1:nrow(ccfs)), ]
          ascwt <- abs(ccfs)-Senuil
          ascwt[ascwt<0] <- 0
          ascwt <- ascwt[rev(1:nrow(ascwt)), ]
          out <- list(ccfs = ccfs, ascwt = ascwt, scales = s)
          return(out)
        }
        library(IRISSeismic)
        library(colorednoise)
        library(BAMBI) # the bivariate von Mises sine distribution
        library(hht)
        library(signal) # unwrap
        library(seewave) # hilbert transform
        set.seed(1)
        datasave <- datasavec <- c()
        datasave.all <- datasave.allc <- c()
        correlationsave <- phasesynsave <- misave <- kullbacksave <- jsave <- dasave <- jrsave <- coherencesave <- list()
        estimates <- list()
        true_synchrony_mu <- list()
        timestart <- Sys.time()
        ###########################################################################################      
        for (times in 1:20) {
          
          print(times)
          rlow <- runif(20, 0.1, 0.3)
          ravg <- runif(20, 0.4, 0.6)
          rhigh <- runif(20, 0.7, 0.9)
          r <- c(rlow, ravg, rhigh)
          computelambda <- function(r) log(1-r^2)/(-2) 
          x <- list();y <- list();xn <- list();xnoise <- list();ynoise <- list()
          xnolag <- ynolag <- list()
          correlation <- c()
          cr <- c()
          kullback <- c();j <- c();mi <- c()
          phasesyn <- c();da <- list();jr <- list()
          coherence <- c()
          
          for (i in 1:length(r)) {
            if (bd=='33') {
              b <- 0.3;d <- 0.3
            } else if (bd=='31') {
              b <- 0.3;d <- 0.1
            } else {
              b <- 0.1;d <- 0.3
            }
            # generate data
            henonmaps_data <- henon_map(b,d,r[[i]],timelength+10000)
            xnolag[[i]] <- henonmaps_data$x[10001:(timelength+10000)] # discard the first 10000 samples
            ynolag[[i]] <- henonmaps_data$y[(10001):(timelength+10000)]
            x[[i]] <- henonmaps_data$x[10001:(timelength+10000-tlag)] # discard the first 10000 samples
            y[[i]] <- henonmaps_data$y[(10001+tlag):(timelength+10000)]
            
            Varnoise <- function(signal, snr) varnoise <- mean(signal^2)/snr
            varofnoise1 <- Varnoise(x[[i]],((1-noiselevel)/noiselevel)) 
            varofnoise2 <- Varnoise(y[[i]],((1-noiselevel)/noiselevel)) 
            noise1 <-colored_noise(timesteps = length(x[[i]]), mean = 0, sd = sqrt(varofnoise1), phi = phivalue) 
            noise2 <-colored_noise(timesteps = length(y[[i]]), mean = 0, sd = sqrt(varofnoise2), phi = phivalue)
            xnoise[[i]] <- x[[i]] + noise1 #add noise to signal
            ynoise[[i]] <- y[[i]] + noise2
            
            xn[[i]] <- scale(cbind(xnolag[[i]],ynolag[[i]]))
            #compute correlation
            correlation[i] <- mean(xn[[i]][,1]*xn[[i]][,2])
            # compute synchrony estimates
            out <- cwt_signal(Sig=xnoise[[i]],Cwt='MORLET',Min=4,Max=30,
                              FreqEch=200,Senuil=0,frqsmp=1,Aff=0,
                              Fc=0.625)
            ascwt <- out$ascwt
            ccfs <- out$ccfs
            s <- out$scales
            
            #Normalization
            for (k in 1:nrow(ascwt)) {
              ascwt[k,] <- ascwt[k, ]/mean(ascwt[k, ])
            }
            ascwt1 <- ascwt/sum(ascwt) #divided by total sum
            
            out <- cwt_signal(Sig=ynoise[[i]],Cwt='MORLET',Min=4,Max=30,
                              FreqEch=200,Senuil=0,frqsmp=1,Aff=0,
                              Fc=0.625)
            
            ascwt <- out$ascwt
            ccfs <- out$ccfs
            s <- out$scales
            # Normalization
            for (k in 1:nrow(ascwt)) {
              ascwt[k, ] <- ascwt[k, ]/mean(ascwt[k, ])
            }
            ascwt2 <- ascwt/sum(ascwt)
            
            # part 3
            v1 <- seq(0.1, 0.9, 0.1)
            v2 <- seq(2, 10, 1)
            inftheo <- divergence_measures(ascwt1, ascwt2, v1, v2)
            print(i)
            
            kullback[i] <- inftheo$K
            j[i] <- inftheo$J
            mi[i] <- inftheo$MI
            da[[i]] <- inftheo$Da
            jr[[i]] <- inftheo$JR
            phasesyn[i] <- filtdir(xnoise[[i]],ynoise[[i]])
            fs <- 10
            ts1 <- ts(xnoise[[i]],frequency = fs)
            ts2 <- ts(ynoise[[i]], frequency = fs)
            DF <- crossSpectrum(ts.union(ts1,ts2), spans=c(5,7))
            coherence[i] <- mean(DF$coh)
            
          }
          
          
          
          #########################
          ##compute MSEs and correlations
          #use mu as the true synchrony
          MSEk_low <- mean((rlow - (1-kullback[1:20]))^2) 
          crk_low <- cor((1-kullback[1:20]),rlow)#-0.4
          MSEj_low <- mean((rlow - (1-j[1:20]))^2)
          crj_low <- cor(1-j[1:20],rlow)#-0.4
          MSEmi_low <- mean((rlow - mi[1:20])^2)
          crmi_low <- cor(mi[1:20],rlow)#0.000217
          MSEphasesyn_low <- mean((rlow - phasesyn[1:20])^2)
          crphasesyn_low <- cor(phasesyn[1:20],rlow)
          MSEcoherence_low <- mean((rlow - (coherence[1:20]))^2)
          crcoherence_low <- cor((coherence[1:20]),rlow)
          #use correlation as the true synchrony
          MSEk_lowc <- mean((correlation[1:20] - (1-kullback[1:20]))^2) 
          crk_lowc <- cor((1-kullback[1:20]),correlation[1:20])#-0.4
          MSEj_lowc <- mean((correlation[1:20] - (1-j[1:20]))^2)
          crj_lowc <- cor(1-j[1:20],correlation[1:20])#-0.4
          MSEmi_lowc <- mean((correlation[1:20] - mi[1:20])^2)
          crmi_lowc <- cor(mi[1:20],correlation[1:20])#0.000217
          MSEphasesyn_lowc <- mean((correlation[1:20] - phasesyn[1:20])^2)
          crphasesyn_lowc <- cor(phasesyn[1:20],correlation[1:20])
          MSEcoherence_lowc <- mean((correlation[1:20] - (coherence[1:20]))^2)
          crcoherence_lowc <- cor((coherence[1:20]),correlation[1:20])
          da01 <-  da05<- da09<- c()
          for (nr in 1:60) {
            da01[nr] <- da[[nr]][1]
            
            da05[nr] <- da[[nr]][5]
            
            da09[nr] <- da[[nr]][9]
          }
          jr2 <- jr6<- jr10<- c()
          for (nr in 1:60) {
            jr2[nr] <- jr[[nr]][1]
            jr6[nr] <- jr[[nr]][5]
            jr10[nr] <- jr[[nr]][9]
          }
          
          MSEda01_low <- mean((rlow - (1-da01[1:20]))^2)
          crda01_low <- cor((1-da01[1:20]),rlow)
          
          MSEda05_low <- mean((rlow - (1-da05[1:20]))^2)
          crda05_low <- cor(1-da05[1:20],rlow)
          
          MSEda09_low <- mean((rlow - (1-da09[1:20]))^2)
          crda09_low <- cor(1-da09[1:20],rlow)
          
          MSEjr2_low <- mean((rlow - (1-jr2[1:20]))^2)
          crjr2_low <- cor((1-jr2[1:20]),rlow)
          
          MSEjr6_low <- mean((rlow - (1-jr6[1:20]))^2)
          crjr6_low <- cor(1-jr6[1:20],rlow)
          
          MSEjr10_low <- mean((rlow - (1-jr10[1:20]))^2)
          crjr10_low <- cor((1-jr10[1:20]),rlow)
          
          MSEda01_lowc <- mean((correlation[1:20] - (1-da01[1:20]))^2)
          crda01_lowc <- cor((1-da01[1:20]),correlation[1:20])
          
          MSEda05_lowc <- mean((correlation[1:20] - (1-da05[1:20]))^2)
          crda05_lowc <- cor(1-da05[1:20],correlation[1:20])
          
          MSEda09_lowc <- mean((correlation[1:20] - (1-da09[1:20]))^2)
          crda09_lowc <- cor(1-da09[1:20],correlation[1:20])
          
          MSEjr2_lowc <- mean((correlation[1:20] - (1-jr2[1:20]))^2)
          crjr2_lowc <- cor((1-jr2[1:20]),correlation[1:20])
          
          MSEjr6_lowc <- mean((correlation[1:20] - (1-jr6[1:20]))^2)
          crjr6_lowc <- cor(1-jr6[1:20],correlation[1:20])
          
          MSEjr10_lowc <- mean((correlation[1:20] - (1-jr10[1:20]))^2)
          crjr10_lowc <- cor((1-jr10[1:20]),correlation[1:20])
          
          
          
          
          resultk_low <- c('Kullback','rlow',bd,tlag,timelength,noiselevel,noisetype,round(MSEk_low,2),round(crk_low,2))
          resultj_low <- c('JensenShannon','rlow',bd,tlag,timelength,noiselevel,noisetype,round(MSEj_low,2),round(crj_low,2))
          resultmi_low <- c('mutual_information','rlow',bd,tlag,timelength,noiselevel,noisetype,round(MSEmi_low,2),round(crmi_low,2))
          resultphasesyn_low <- c('phasesynchrony','rlow',bd,tlag,timelength,noiselevel,noisetype,round(MSEphasesyn_low,2),round(crphasesyn_low,2))
          resultda01_low <- c('Renyi01','rlow',bd,tlag,timelength,noiselevel,noisetype,round(MSEda01_low,2),round(crda01_low,2))
          resultda05_low <- c('Renyi05','rlow',bd,tlag,timelength,noiselevel,noisetype,round(MSEda05_low,2),round(crda05_low,2))
          resultda09_low <- c('Renyi09','rlow',bd,tlag,timelength,noiselevel,noisetype,round(MSEda09_low,2),round(crda09_low,2))
          
          resultjr2_low <- c('JensenRenyi2','rlow',bd,tlag,timelength,noiselevel,noisetype,round(MSEjr2_low,2),round(crjr2_low,2))
          resultjr6_low <- c('JensenRenyi6','rlow',bd,tlag,timelength,noiselevel,noisetype,round(MSEjr6_low,2),round(crjr6_low,2))
          resultjr10_low <- c('JensenRenyi10','rlow',bd,tlag,timelength,noiselevel,noisetype,round(MSEjr10_low,2),round(crjr10_low,2))
          
          
          resultk_lowc <- c('Kullback','rlow',bd,tlag,timelength,noiselevel,noisetype,round(MSEk_lowc,2),round(crk_lowc,2))
          resultj_lowc <- c('JensenShannon','rlow',bd,tlag,timelength,noiselevel,noisetype,round(MSEj_lowc,2),round(crj_lowc,2))
          resultmi_lowc <- c('mutual_information','rlow',bd,tlag,timelength,noiselevel,noisetype,round(MSEmi_lowc,2),round(crmi_lowc,2))
          resultphasesyn_lowc <- c('phasesynchrony','rlow',bd,tlag,timelength,noiselevel,noisetype,round(MSEphasesyn_lowc,2),round(crphasesyn_lowc,2))
          resultda01_lowc <- c('Renyi01','rlow',bd,tlag,timelength,noiselevel,noisetype,round(MSEda01_lowc,2),round(crda01_lowc,2))
          resultda05_lowc <- c('Renyi05','rlow',bd,tlag,timelength,noiselevel,noisetype,round(MSEda05_lowc,2),round(crda05_lowc,2))
          resultda09_lowc <- c('Renyi09','rlow',bd,tlag,timelength,noiselevel,noisetype,round(MSEda09_lowc,2),round(crda09_lowc,2))
          
          resultjr2_lowc <- c('JensenRenyi2','rlow',bd,tlag,timelength,noiselevel,noisetype,round(MSEjr2_lowc,2),round(crjr2_lowc,2))
          resultjr6_lowc <- c('JensenRenyi6','rlow',bd,tlag,timelength,noiselevel,noisetype,round(MSEjr6_lowc,2),round(crjr6_lowc,2))
          resultjr10_lowc <- c('JensenRenyi10','rlow',bd,tlag,timelength,noiselevel,noisetype,round(MSEjr10_lowc,2),round(crjr10_lowc,2))
          
          resultcoherence_low <- c('coherence','rlow',bd,tlag,timelength,noiselevel,noisetype,round(MSEcoherence_low,2),round(crcoherence_low,2))
          resultcoherence_lowc <- c('coherence','rlow',bd,tlag,timelength,noiselevel,noisetype,round(MSEcoherence_lowc,2),round(crcoherence_lowc,2))
          
          ####################################################################################################
          #avg
          MSEk_avg <- mean((ravg - (1-kullback[21:40]))^2) 
          crk_avg <- cor((1-kullback[21:40]),ravg)#-0.4
          MSEj_avg <- mean((ravg - (1-j[21:40]))^2)
          crj_avg <- cor(1-j[21:40],ravg)#-0.4
          MSEmi_avg <- mean((ravg - mi[21:40])^2)
          crmi_avg <- cor(mi[21:40],ravg)#0.000217
          MSEphasesyn_avg <- mean((ravg - phasesyn[21:40])^2)
          crphasesyn_avg <- cor(phasesyn[21:40],ravg)
          
          MSEda01_avg <- mean((ravg - (1-da01[21:40]))^2)
          crda01_avg <- cor((1-da01[21:40]),ravg)
          
          MSEda05_avg <- mean((ravg - (1-da05[21:40]))^2)
          crda05_avg <- cor(1-da05[21:40],ravg)
          
          MSEda09_avg <- mean((ravg - (1-da09[21:40]))^2)
          crda09_avg <- cor(1-da09[21:40],ravg)
          
          MSEjr2_avg <- mean((ravg - (1-jr2[21:40]))^2)
          crjr2_avg <- cor((1-jr2[21:40]),ravg)
          MSEjr6_avg <- mean((ravg - (1-jr6[21:40]))^2)
          crjr6_avg <- cor(1-jr6[21:40],ravg)
          
          MSEjr10_avg <- mean((ravg - (1-jr10[21:40]))^2)
          crjr10_avg <- cor((1-jr10[21:40]),ravg)
          MSEcoherence_avg <- mean((ravg - (coherence[21:40]))^2) 
          crcoherence_avg <- cor((coherence[21:40]),ravg)
          
          #
          MSEk_avgc <- mean((correlation[21:40] - (1-kullback[21:40]))^2) 
          crk_avgc <- cor((1-kullback[21:40]),correlation[21:40])#-0.4
          MSEj_avgc <- mean((correlation[21:40] - (1-j[21:40]))^2)
          crj_avgc <- cor(1-j[21:40],correlation[21:40])#-0.4
          MSEmi_avgc <- mean((correlation[21:40] - mi[21:40])^2)
          crmi_avgc <- cor(mi[21:40],correlation[21:40])#0.000217
          
          MSEphasesyn_avgc <- mean((correlation[21:40] - phasesyn[21:40])^2)
          crphasesyn_avgc <- cor(phasesyn[21:40],correlation[21:40])
          
          
          
          MSEda01_avgc <- mean((correlation[21:40] - (1-da01[21:40]))^2)
          crda01_avgc <- cor((1-da01[21:40]),correlation[21:40])
          
          MSEda05_avgc <- mean((correlation[21:40] - (1-da05[21:40]))^2)
          crda05_avgc <- cor(1-da05[21:40],correlation[21:40])
          
          MSEda09_avgc <- mean((correlation[21:40] - (1-da09[21:40]))^2)
          crda09_avgc <- cor(1-da09[21:40],correlation[21:40])
          
          MSEjr2_avgc <- mean((correlation[21:40] - (1-jr2[21:40]))^2)
          crjr2_avgc <- cor((1-jr2[21:40]),correlation[21:40])
          MSEjr6_avgc <- mean((correlation[21:40] - (1-jr6[21:40]))^2)
          crjr6_avgc <- cor(1-jr6[21:40],correlation[21:40])
          
          MSEjr10_avgc <- mean((correlation[21:40] - (1-jr10[21:40]))^2)
          crjr10_avgc <- cor((1-jr10[21:40]),correlation[21:40])
          MSEcoherence_avgc <- mean((correlation[21:40] - (coherence[21:40]))^2) 
          crcoherence_avgc <- cor((coherence[21:40]),correlation[21:40])
          
          resultk_avg <- c('Kullback','ravg',bd,tlag,timelength,noiselevel,noisetype,round(MSEk_avg,2),round(crk_avg,2))
          resultj_avg <- c('JensenShannon','ravg',bd,tlag,timelength,noiselevel,noisetype,round(MSEj_avg,2),round(crj_avg,2))
          resultmi_avg <- c('mutual_information','ravg',bd,tlag,timelength,noiselevel,noisetype,round(MSEmi_avg,2),round(crmi_avg,2))
          resultphasesyn_avg <- c('phasesynchrony','ravg',bd,tlag,timelength,noiselevel,noisetype,round(MSEphasesyn_avg,2),round(crphasesyn_avg,2))
          resultda01_avg <- c('Renyi01','ravg',bd,tlag,timelength,noiselevel,noisetype,round(MSEda01_avg,2),round(crda01_avg,2))
          resultda05_avg <- c('Renyi05','ravg',bd,tlag,timelength,noiselevel,noisetype,round(MSEda05_avg,2),round(crda05_avg,2))
          resultda09_avg <- c('Renyi09','ravg',bd,tlag,timelength,noiselevel,noisetype,round(MSEda09_avg,2),round(crda09_avg,2))
          
          resultjr2_avg <- c('JensenRenyi2','ravg',bd,tlag,timelength,noiselevel,noisetype,round(MSEjr2_avg,2),round(crjr2_avg,2))
          resultjr6_avg <- c('JensenRenyi6','ravg',bd,tlag,timelength,noiselevel,noisetype,round(MSEjr6_avg,2),round(crjr6_avg,2))
          resultjr10_avg <- c('JensenRenyi10','ravg',bd,tlag,timelength,noiselevel,noisetype,round(MSEjr10_avg,2),round(crjr10_avg,2))
          
          
          resultk_avgc <- c('Kullback','ravg',bd,tlag,timelength,noiselevel,noisetype,round(MSEk_avgc,2),round(crk_avgc,2))
          resultj_avgc <- c('JensenShannon','ravg',bd,tlag,timelength,noiselevel,noisetype,round(MSEj_avgc,2),round(crj_avgc,2))
          resultmi_avgc <- c('mutual_information','ravg',bd,tlag,timelength,noiselevel,noisetype,round(MSEmi_avgc,2),round(crmi_avgc,2))
          resultphasesyn_avgc <- c('phasesynchrony','ravg',bd,tlag,timelength,noiselevel,noisetype,round(MSEphasesyn_avgc,2),round(crphasesyn_avgc,2))
          resultda01_avgc <- c('Renyi01','ravg',bd,tlag,timelength,noiselevel,noisetype,round(MSEda01_avgc,2),round(crda01_avgc,2))
          resultda05_avgc <- c('Renyi05','ravg',bd,tlag,timelength,noiselevel,noisetype,round(MSEda05_avgc,2),round(crda05_avgc,2))
          resultda09_avgc <- c('Renyi09','ravg',bd,tlag,timelength,noiselevel,noisetype,round(MSEda09_avgc,2),round(crda09_avgc,2))
          
          resultjr2_avgc <- c('JensenRenyi2','ravg',bd,tlag,timelength,noiselevel,noisetype,round(MSEjr2_avgc,2),round(crjr2_avgc,2))
          resultjr6_avgc <- c('JensenRenyi6','ravg',bd,tlag,timelength,noiselevel,noisetype,round(MSEjr6_avgc,2),round(crjr6_avgc,2))
          resultjr10_avgc <- c('JensenRenyi10','ravg',bd,tlag,timelength,noiselevel,noisetype,round(MSEjr10_avgc,2),round(crjr10_avgc,2))
          resultcoherence_avg <- c('coherence','ravg',bd,tlag,timelength,noiselevel,noisetype,round(MSEcoherence_avg,2),round(crcoherence_avg,2))
          resultcoherence_avgc <- c('coherence','ravg',bd,tlag,timelength,noiselevel,noisetype,round(MSEcoherence_avgc,2),round(crcoherence_avgc,2))
          
          #################################
          ###################################
          #high
          
          MSEk_high <- mean((rhigh - (1-kullback[41:60]))^2) 
          crk_high <- cor((1-kullback[41:60]),rhigh)#-0.4
          MSEj_high <- mean((rhigh - (1-j[41:60]))^2)
          crj_high <- cor(1-j[41:60],rhigh)#-0.4
          MSEmi_high <- mean((rhigh - mi[41:60])^2)
          crmi_high <- cor(mi[41:60],rhigh)#0.000217
          
          MSEc_high <- mean((rhigh - correlation[41:60])^2)
          crc_high <- cor(correlation[41:60],rhigh) 
          MSEphasesyn_high <- mean((rhigh - phasesyn[41:60])^2)
          crphasesyn_high <- cor(phasesyn[41:60],rhigh)
          
          
          
          MSEda01_high <- mean((rhigh - (1-da01[41:60]))^2)
          crda01_high <- cor((1-da01[41:60]),rhigh)
          MSEda05_high <- mean((rhigh - (1-da05[41:60]))^2)
          crda05_high <- cor(1-da05[41:60],rhigh)
          MSEda09_high <- mean((rhigh - (1-da09[41:60]))^2)
          crda09_high <- cor(1-da09[41:60],rhigh)
          
          MSEjr2_high <- mean((rhigh - (1-jr2[41:60]))^2)
          crjr2_high <- cor((1-jr2[41:60]),rhigh)
          
          MSEjr6_high <- mean((rhigh - (1-jr6[41:60]))^2)
          crjr6_high <- cor(1-jr6[41:60],rhigh)
          
          MSEjr10_high <- mean((rhigh - (1-jr10[41:60]))^2)
          crjr10_high <- cor((1-jr10[41:60]),rhigh)
          MSEcoherence_high <- mean((rhigh - coherence[41:60])^2) 
          crcoherence_high <- cor(coherence[41:60],rhigh)
          
          #
          MSEk_highc <- mean((correlation[41:60] - (1-kullback[41:60]))^2) 
          crk_highc <- cor((1-kullback[41:60]),correlation[41:60])#-0.4
          MSEj_highc <- mean((correlation[41:60] - (1-j[41:60]))^2)
          crj_highc <- cor(1-j[41:60],correlation[41:60])#-0.4
          MSEmi_highc <- mean((correlation[41:60] - mi[41:60])^2)
          crmi_highc <- cor(mi[41:60],correlation[41:60])#0.000217
          
          MSEphasesyn_highc <- mean((correlation[41:60] - phasesyn[41:60])^2)
          crphasesyn_highc <- cor(phasesyn[41:60],correlation[41:60])
          
          
          
          MSEda01_highc <- mean((correlation[41:60] - (1-da01[41:60]))^2)
          crda01_highc <- cor((1-da01[41:60]),correlation[41:60])
          MSEda05_highc <- mean((correlation[41:60] - (1-da05[41:60]))^2)
          crda05_highc <- cor(1-da05[41:60],correlation[41:60])
          MSEda09_highc <- mean((correlation[41:60] - (1-da09[41:60]))^2)
          crda09_highc <- cor(1-da09[41:60],correlation[41:60])
          
          MSEjr2_highc <- mean((correlation[41:60] - (1-jr2[41:60]))^2)
          crjr2_highc <- cor((1-jr2[41:60]),correlation[41:60])
          
          MSEjr6_highc <- mean((correlation[41:60] - (1-jr6[41:60]))^2)
          crjr6_highc <- cor(1-jr6[41:60],correlation[41:60])
          
          MSEjr10_highc <- mean((correlation[41:60] - (1-jr10[41:60]))^2)
          crjr10_highc <- cor((1-jr10[41:60]),correlation[41:60])
          
          MSEcoherence_highc <- mean((correlation[41:60] - coherence[41:60])^2) 
          crcoherence_highc <- cor(coherence[41:60],correlation[41:60])
          
          resultk_high <- c('Kullback','rhigh',bd,tlag,timelength,noiselevel,noisetype,round(MSEk_high,2),round(crk_high,2))
          resultj_high <- c('JensenShannon','rhigh',bd,tlag,timelength,noiselevel,noisetype,round(MSEj_high,2),round(crj_high,2))
          resultmi_high <- c('mutual_information','rhigh',bd,tlag,timelength,noiselevel,noisetype,round(MSEmi_high,2),round(crmi_high,2))
          resultphasesyn_high <- c('phasesynchrony','rhigh',bd,tlag,timelength,noiselevel,noisetype,round(MSEphasesyn_high,2),round(crphasesyn_high,2))
          resultda01_high <- c('Renyi01','rhigh',bd,tlag,timelength,noiselevel,noisetype,round(MSEda01_high,2),round(crda01_high,2))
          resultda05_high <- c('Renyi05','rhigh',bd,tlag,timelength,noiselevel,noisetype,round(MSEda05_high,2),round(crda05_high,2))
          resultda09_high <- c('Renyi09','rhigh',bd,tlag,timelength,noiselevel,noisetype,round(MSEda09_high,2),round(crda09_high,2))
          
          resultjr2_high <- c('JensenRenyi2','rhigh',bd,tlag,timelength,noiselevel,noisetype,round(MSEjr2_high,2),round(crjr2_high,2))
          resultjr6_high <- c('JensenRenyi6','rhigh',bd,tlag,timelength,noiselevel,noisetype,round(MSEjr6_high,2),round(crjr6_high,2))
          resultjr10_high <- c('JensenRenyi10','rhigh',bd,tlag,timelength,noiselevel,noisetype,round(MSEjr10_high,2),round(crjr10_high,2))
          
          
          
          resultk_highc <- c('Kullback','rhigh',bd,tlag,timelength,noiselevel,noisetype,round(MSEk_highc,2),round(crk_highc,2))
          resultj_highc <- c('JensenShannon','rhigh',bd,tlag,timelength,noiselevel,noisetype,round(MSEj_highc,2),round(crj_highc,2))
          resultmi_highc <- c('mutual_information','rhigh',bd,tlag,timelength,noiselevel,noisetype,round(MSEmi_highc,2),round(crmi_highc,2))
          resultphasesyn_highc <- c('phasesynchrony','rhigh',bd,tlag,timelength,noiselevel,noisetype,round(MSEphasesyn_highc,2),round(crphasesyn_highc,2))
          resultda01_highc <- c('Renyi01','rhigh',bd,tlag,timelength,noiselevel,noisetype,round(MSEda01_highc,2),round(crda01_highc,2))
          resultda05_highc <- c('Renyi05','rhigh',bd,tlag,timelength,noiselevel,noisetype,round(MSEda05_highc,2),round(crda05_highc,2))
          resultda09_highc <- c('Renyi09','rhigh',bd,tlag,timelength,noiselevel,noisetype,round(MSEda09_highc,2),round(crda09_highc,2))
          
          resultjr2_highc <- c('JensenRenyi2','rhigh',bd,tlag,timelength,noiselevel,noisetype,round(MSEjr2_highc,2),round(crjr2_highc,2))
          resultjr6_highc <- c('JensenRenyi6','rhigh',bd,tlag,timelength,noiselevel,noisetype,round(MSEjr6_highc,2),round(crjr6_highc,2))
          resultjr10_highc <- c('JensenRenyi10','rhigh',bd,tlag,timelength,noiselevel,noisetype,round(MSEjr10_highc,2),round(crjr10_highc,2))
          resultcoherence_high <- c('coherence','rhigh',bd,tlag,timelength,noiselevel,noisetype,round(MSEcoherence_high,2),round(crcoherence_high,2))
          resultcoherence_highc <- c('coherence','rhigh',bd,tlag,timelength,noiselevel,noisetype,round(MSEcoherence_highc,2),round(crcoherence_highc,2))
          
          
          
          ###########################################################
          #overall
          MSEk <- mean((r - (1-kullback))^2) 
          crk <- cor((1-kullback),r)
          MSEj <- mean((r - (1-j))^2)
          crj <- cor(1-j,r)
          MSEmi <- mean((r - mi)^2)
          crmi <- cor(mi,r)
          
          MSEc <- mean((r - correlation)^2)
          crc <- cor(correlation,r) 
          MSEphasesyn <- mean((r - phasesyn)^2)
          crphasesyn <- cor(phasesyn,r)
          
          
          
          MSEda01 <- mean((r - (1-da01))^2)
          crda01 <- cor(1-da01,r)
          
          MSEda05 <- mean((r - (1-da05))^2)
          crda05 <- cor(1-da05,r)
          
          MSEda09 <- mean((r - (1-da09))^2)
          crda09 <- cor(1-da09,r)
          
          MSEjr2 <- mean((r - (1-jr2))^2)
          crjr2 <- cor((1-jr2),r)
          
          MSEjr6 <- mean((r - (1-jr6))^2)
          crjr6 <- cor(1-jr6,r)
          
          MSEjr10 <- mean((r - (1-jr10))^2)
          crjr10 <- cor((1-jr10),r)
          
          MSEcoherence <- mean((r - coherence)^2) 
          crcoherence <- cor(coherence,r)
          
          
          
          resultk <- c('Kullback','r',bd,tlag,timelength,noiselevel,noisetype,round(MSEk,2),round(crk,2))
          resultj <- c('JensenShannon','r',bd,tlag,timelength,noiselevel,noisetype,round(MSEj,2),round(crj,2))
          resultmi <- c('mutual_information','r',bd,tlag,timelength,noiselevel,noisetype,round(MSEmi,2),round(crmi,2))
          resultphasesyn <- c('phasesynchrony','r',bd,tlag,timelength,noiselevel,noisetype,round(MSEphasesyn,2),round(crphasesyn,2))
          resultda01 <- c('Renyi01','r',bd,tlag,timelength,noiselevel,noisetype,round(MSEda01,2),round(crda01,2))
          resultda05 <- c('Renyi05','r',bd,tlag,timelength,noiselevel,noisetype,round(MSEda05,2),round(crda05,2))
          resultda09 <- c('Renyi09','r',bd,tlag,timelength,noiselevel,noisetype,round(MSEda09,2),round(crda09,2))
          
          resultjr2 <- c('JensenRenyi2','r',bd,tlag,timelength,noiselevel,noisetype,round(MSEjr2,2),round(crjr2,2))
          resultjr6 <- c('JensenRenyi6','r',bd,tlag,timelength,noiselevel,noisetype,round(MSEjr6,2),round(crjr6,2))
          resultjr10 <- c('JensenRenyi10','r',bd,tlag,timelength,noiselevel,noisetype,round(MSEjr10,2),round(crjr10,2))
          
          resultcoherence <- c('coherence','r',bd,tlag,timelength,noiselevel,noisetype,round(MSEcoherence,2),round(crcoherence,2))
          
          
          MSEkc <- mean((correlation - (1-kullback))^2) 
          crkc <- cor((1-kullback),correlation)
          MSEjc <- mean((correlation - (1-j))^2)
          crjc <- cor(1-j,correlation)
          MSEmic <- mean((correlation - mi)^2)
          crmic <- cor(mi,correlation)
          
          
          MSEphasesync <- mean((correlation - phasesyn)^2)
          crphasesync <- cor(phasesyn,correlation)
          
          
          
          MSEda01c <- mean((r - (correlation-da01))^2)
          crda01c <- cor(1-da01,correlation)
          
          MSEda05c <- mean((correlation - (1-da05))^2)
          crda05c <- cor(1-da05,correlation)
          
          MSEda09c <- mean((correlation - (1-da09))^2)
          crda09c <- cor(1-da09,correlation)
          
          MSEjr2c <- mean((correlation - (1-jr2))^2)
          crjr2c <- cor((1-jr2),correlation)
          
          MSEjr6c <- mean((correlation - (1-jr6))^2)
          crjr6c <- cor(1-jr6,correlation)
          
          MSEjr10c <- mean((correlation - (1-jr10))^2)
          crjr10c <- cor((1-jr10),correlation)
          MSEcoherencec <- mean((r - correlation)^2) 
          crcoherencec <- cor(correlation,r)
          
          resultkc <- c('Kullback','r',bd,tlag,timelength,noiselevel,noisetype,round(MSEkc,2),round(crkc,2))
          resultjc <- c('JensenShannon','r',bd,tlag,timelength,noiselevel,noisetype,round(MSEjc,2),round(crjc,2))
          resultmic <- c('mutual_information','r',bd,tlag,timelength,noiselevel,noisetype,round(MSEmic,2),round(crmic,2))
          resultphasesync <- c('phasesynchrony','r',bd,tlag,timelength,noiselevel,noisetype,round(MSEphasesync,2),round(crphasesync,2))
          resultda01c <- c('Renyi01','r',bd,tlag,timelength,noiselevel,noisetype,round(MSEda01c,2),round(crda01c,2))
          resultda05c <- c('Renyi05','r',bd,tlag,timelength,noiselevel,noisetype,round(MSEda05c,2),round(crda05c,2))
          resultda09c <- c('Renyi09','r',bd,tlag,timelength,noiselevel,noisetype,round(MSEda09c,2),round(crda09c,2))
          
          resultjr2c <- c('JensenRenyi2','r',bd,tlag,timelength,noiselevel,noisetype,round(MSEjr2c,2),round(crjr2c,2))
          resultjr6c <- c('JensenRenyi6','r',bd,tlag,timelength,noiselevel,noisetype,round(MSEjr6c,2),round(crjr6c,2))
          resultjr10c <- c('JensenRenyi10','r',bd,tlag,timelength,noiselevel,noisetype,round(MSEjr10c,2),round(crjr10c,2))
          resultcoherencec <- c('coherence','r',bd,tlag,timelength,noiselevel,noisetype,round(MSEcoherencec,2),round(crcoherencec,2))
          
          
          
          #######################################
          datasave <- rbind(datasave,resultk_low,resultj_low,resultmi_low,
                            resultphasesyn_low,resultda01_low,resultda05_low,resultda09_low
                            ,resultjr2_low,resultjr6_low,resultjr10_low,resultcoherence_low,
                            resultk_avg,resultj_avg,resultmi_avg,
                            resultphasesyn_avg,resultda01_avg,resultda05_avg,resultda09_avg,resultjr2_avg,
                            resultjr6_avg,resultjr10_avg,resultcoherence_avg,
                            resultk_high,resultj_high,resultmi_high,
                            resultphasesyn_high,resultda01_high,resultda05_high,resultda09_high,resultjr2_high,
                            resultjr6_high,resultjr10_high,resultcoherence_high
          )
          
          datasave.all <- rbind(datasave.all,resultk,resultj,resultmi,resultphasesyn,resultda01,
                                resultda05,resultda09,resultjr2,resultjr6,resultjr10,resultcoherence)
          
          datasavec <- rbind(datasavec,resultk_lowc,resultj_lowc,resultmi_lowc,
                             resultphasesyn_lowc,resultda01_lowc,resultda05_lowc,resultda09_lowc
                             ,resultjr2_lowc,resultjr6_lowc,resultjr10_lowc,resultcoherence_lowc,
                             resultk_avgc,resultj_avgc,resultmi_avgc,
                             resultphasesyn_avgc,resultda01_avgc,resultda05_avgc,resultda09_avgc,resultjr2_avgc,
                             resultjr6_avgc,resultjr10_avgc,resultcoherence_avgc,
                             resultk_highc,resultj_highc,resultmi_highc,
                             resultphasesyn_highc,resultda01_highc,resultda05_highc,resultda09_highc,resultjr2_highc,
                             resultjr6_highc,resultjr10_highc,resultcoherence_highc
          )
          
          datasave.allc <- rbind(datasave.allc,resultkc,resultjc,resultmic,resultphasesync,resultda01c,
                                 resultda05c,resultda09c,resultjr2c,resultjr6c,resultjr10c,resultcoherencec)
          
          
          
          
          
          true_synchrony_mu[[times]] <- r
          correlationsave[[times]] <- correlation
          phasesynsave[[times]] <- phasesyn
          misave[[times]] <- mi
          jsave[[times]] <- j
          kullbacksave[[times]] <- kullback
          dasave[[times]] <- da
          jrsave[[times]] <- jr 
          coherencesave[[times]] <- coherence
          
          estimates$true_synchrony_mu <- true_synchrony_mu 
          estimates$correlation <- correlationsave
          estimates$phasesyn <- phasesynsave
          estimates$mi <- misave
          estimates$j <- jsave
          estimates$kullback <- kullbacksave
          estimates$da <- dasave
          estimates$jr <- jrsave
          estimates$coherence <- coherencesave
          
          
        }
        timeend <- Sys.time()
        
        save(datasave,file = paste(c('bd',bd,'_lag',tlag,'_timelength',timelength,'_noiselevel',noiselevel,'_noisetype',noisetype,'.RData'),sep = '',collapse = ''))
        save(datasave.all,file = paste(c('overall_','bd',bd,'_lag',tlag,'_timelength',timelength,'_noiselevel',noiselevel,'_noisetype',noisetype,'.RData'),sep = '',collapse = ''))
        
        save(datasavec,file = paste(c('c','bd' ,bd,'_lag',tlag,'_timelength',timelength,'_noiselevel',noiselevel,'_noisetype',noisetype,'.RData'),sep = '',collapse = ''))
        save(datasave.allc,file = paste(c('c','overall_','bd',bd,'_lag',tlag,'_timelength',timelength,'_noiselevel',noiselevel,'_noisetype',noisetype,'.RData'),sep = '',collapse = ''))
        
        save(estimates,file=paste(c('estimates','bd',bd,'_lag',tlag,'_timelength',timelength,'_noiselevel',noiselevel,'_noisetype',noisetype,'.RData'),sep = '',collapse = '')
        )
        timeend-timestart
      }
    }
  }}






















