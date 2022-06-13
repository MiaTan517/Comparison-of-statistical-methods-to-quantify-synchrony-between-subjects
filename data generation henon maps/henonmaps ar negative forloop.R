timelength_sets <- c(1006,10054)
bd_sets <- c('33','31','13')



noiselevel_sets <- c(0.2,0.5)

noisetype <- 'arn' #change manually
phivalue <- -0.5 # ar correlation

tlag <- 0 #change manually

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
      
      library(colorednoise)
      library(BAMBI) # the bivariate von Mises sine distribution
      library(hht)
      library(signal) # unwrap
      library(seewave) # hilbert transform
      set.seed(1)
      datasave <- c()
      datasave.all <- c()
      xsave <- list()
      ysave <- list()
      crosscorrelationsave <- phasesynsave <- misave <- kullbacksave <- jsave <- dasave <- jrsave <- list()
      estimates <- list()
      timestart <- Sys.time()
      
      for (times in 1:20) {
        
        print(times)
        rlow <- runif(20, 0.1, 0.3)
        ravg <- runif(20, 0.4, 0.6)
        rhigh <- runif(20, 0.7, 0.9)
        r <- c(rlow, ravg, rhigh)
        computelambda <- function(r) log(1-r^2)/(-2) 
        x <- list();y <- list();xn <- list();xnoise <- list();ynoise <- list()
        crosscorrelation <- c()
        cr <- c()
        kullback <- c();j <- c();mi <- c()
        phasesyn <- c();da <- list();jr <- list()
        
        for (i in 1:length(r)) {
          if (bd=='33') {
            b <- 3;d <- 3
          } else if (bd=='31') {
            b <- 3;d <- 1
          } else {
            b <- 1;d <- 3
          }
          henonmaps_data <- henon_map(b,d,r[[i]],timelength+10000)
          x[[i]] <- henonmaps_data$x[10001:(timelength+10000)] # discard the first 10000 samples
          y[[i]] <- henonmaps_data$y[10001:(timelength+10000)]
          
          Varnoise <- function(signal, snr) varnoise <- mean(signal^2)/snr
          varofnoise1 <- Varnoise(x[[i]],((1-noiselevel)/noiselevel)) 
          varofnoise2 <- Varnoise(y[[i]],((1-noiselevel)/noiselevel)) 
          noise1 <-colored_noise(timesteps = length(x[[i]]), mean = 0, sd = sqrt(varofnoise1), phi = phivalue) 
          noise2 <-colored_noise(timesteps = length(y[[i]]), mean = 0, sd = sqrt(varofnoise2), phi = phivalue)
          xnoise[[i]] <- x[[i]] + noise1 #add noise to signal
          ynoise[[i]] <- y[[i]] + noise2
          
          xn[[i]] <- scale(cbind(xnoise[[i]],ynoise[[i]]))
          crosscorrelation[i] <- mean(xn[[i]][,1]*xn[[i]][,2])
          out <- cwt_signal(Sig=xnoise[[i]],Cwt='MORLET',Min=4,Max=30,
                            FreqEch=200,Senuil=0,frqsmp=1,Aff=0,
                            Fc=0.8125)
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
                            Fc=0.8125)
          
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
        }
        #########################
        #compute MSEs and correlations
        MSEk_low <- mean((rlow - (1-kullback[1:20]))^2) 
        crk_low <- cor((1-kullback[1:20]),rlow)#-0.4
        MSEj_low <- mean((rlow - (1-j[1:20]))^2)
        crj_low <- cor(1-j[1:20],rlow)#-0.4
        MSEmi_low <- mean((rlow - mi[1:20])^2)
        crmi_low <- cor(mi[1:20],rlow)#0.000217
        
        MSEc_low <- mean((rlow - crosscorrelation[1:20])^2)
        crc_low <- cor(crosscorrelation[1:20],rlow) 
        MSEphasesyn_low <- mean((rlow - phasesyn[1:20])^2)
        crphasesyn_low <- cor(phasesyn[1:20],rlow)
        
        da01 <- da02 <- da03<- da04<- da05<- da06<- da07<- da08<- da09<- c()
        for (nr in 1:60) {
          da01[nr] <- da[[nr]][1]
          da02[nr] <- da[[nr]][2]
          da03[nr] <- da[[nr]][3]
          da04[nr] <- da[[nr]][4]
          da05[nr] <- da[[nr]][5]
          da06[nr] <- da[[nr]][6]
          da07[nr] <- da[[nr]][7]
          da08[nr] <- da[[nr]][8]
          da09[nr] <- da[[nr]][9]
        }
        jr2 <- jr3 <- jr4<- jr5<- jr6<- jr7<- jr8<- jr9<- jr10<- c()
        for (nr in 1:60) {
          jr2[nr] <- jr[[nr]][1]
          jr3[nr] <- jr[[nr]][2]
          jr4[nr] <- jr[[nr]][3]
          jr5[nr] <- jr[[nr]][4]
          jr6[nr] <- jr[[nr]][5]
          jr7[nr] <- jr[[nr]][6]
          jr8[nr] <- jr[[nr]][7]
          jr9[nr] <- jr[[nr]][8]
          jr10[nr] <- jr[[nr]][9]
        }
        
        MSEda01_low <- mean((rlow - (1-da01[1:20]))^2)
        crda01_low <- cor((1-da01[1:20]),rlow)
        MSEda02_low <- mean((rlow - (1-da02[1:20]))^2)
        crda02_low <- cor(1-da02[1:20],rlow)
        MSEda03_low <- mean((rlow - (1-da03[1:20]))^2)
        crda03_low <- cor(1-da03[1:20],rlow)
        MSEda04_low <- mean((rlow - (1-da04[1:20]))^2)
        crda04_low <- cor(1-da04[1:20],rlow)
        MSEda05_low <- mean((rlow - (1-da05[1:20]))^2)
        crda05_low <- cor(1-da05[1:20],rlow)
        MSEda06_low <- mean((rlow - (1-da06[1:20]))^2)
        crda06_low <- cor(1-da06[1:20],rlow)
        MSEda07_low <- mean((rlow - (1-da07[1:20]))^2)
        crda07_low <- cor(1-da07[1:20],rlow)
        MSEda08_low <- mean((rlow - (1-da08[1:20]))^2)
        crda08_low <- cor(1-da08[1:20],rlow)
        MSEda09_low <- mean((rlow - (1-da09[1:20]))^2)
        crda09_low <- cor(1-da09[1:20],rlow)
        
        MSEjr2_low <- mean((rlow - (1-jr2[1:20]))^2)
        crjr2_low <- cor((1-jr2[1:20]),rlow)
        MSEjr3_low <- mean((rlow - (1-jr3[1:20]))^2)
        crjr3_low <- cor(1-jr3[1:20],rlow)
        MSEjr4_low <- mean((rlow - (1-jr4[1:20]))^2)
        crjr4_low <- cor(1-jr4[1:20],rlow)
        MSEjr5_low <- mean((rlow - (1-jr5[1:20]))^2)
        crjr5_low <- cor(1-jr5[1:20],rlow)
        MSEjr6_low <- mean((rlow - (1-jr6[1:20]))^2)
        crjr6_low <- cor(1-jr6[1:20],rlow)
        MSEjr7_low <- mean((rlow - (1-jr7[1:20]))^2)
        crjr7_low <- cor(1-jr7[1:20],rlow)
        MSEjr8_low <- mean((rlow - (1-jr8[1:20]))^2)
        crjr8_low <- cor(1-jr8[1:20],rlow)
        MSEjr9_low <- mean((rlow - (1-jr9[1:20]))^2)
        crjr9_low <- cor(1-jr9[1:20],rlow)
        MSEjr10_low <- mean((rlow - (1-jr10[1:20]))^2)
        crjr10_low <- cor((1-jr10[1:20]),rlow)
        
        
        resultk_low <- c('Kullback','rlow',bd,tlag,timelength,noiselevel,noisetype,round(MSEk_low,2),round(crk_low,2))
        resultj_low <- c('JensenShannon','rlow',bd,tlag,timelength,noiselevel,noisetype,round(MSEj_low,2),round(crj_low,2))
        resultmi_low <- c('mutual_information','rlow',bd,tlag,timelength,noiselevel,noisetype,round(MSEmi_low,2),round(crmi_low,2))
        resultc_low <- c('cross_correlation','rlow',bd,tlag,timelength,noiselevel,noisetype,round(MSEc_low,2),round(crc_low,2))
        resultphasesyn_low <- c('phasesynchrony','rlow',bd,tlag,timelength,noiselevel,noisetype,round(MSEphasesyn_low,2),round(crphasesyn_low,2))
        resultda01_low <- c('Renyi01','rlow',bd,tlag,timelength,noiselevel,noisetype,round(MSEda01_low,2),round(crda01_low,2))
        resultda02_low <- c('Renyi02','rlow',bd,tlag,timelength,noiselevel,noisetype,round(MSEda02_low,2),round(crda02_low,2))
        resultda03_low <- c('Renyi03','rlow',bd,tlag,timelength,noiselevel,noisetype,round(MSEda03_low,2),round(crda03_low,2))
        resultda04_low <- c('Renyi04','rlow',bd,tlag,timelength,noiselevel,noisetype,round(MSEda04_low,2),round(crda04_low,2))
        resultda05_low <- c('Renyi05','rlow',bd,tlag,timelength,noiselevel,noisetype,round(MSEda05_low,2),round(crda05_low,2))
        resultda06_low <- c('Renyi06','rlow',bd,tlag,timelength,noiselevel,noisetype,round(MSEda06_low,2),round(crda06_low,2))
        resultda07_low <- c('Renyi07','rlow',bd,tlag,timelength,noiselevel,noisetype,round(MSEda07_low,2),round(crda07_low,2))
        resultda08_low <- c('Renyi08','rlow',bd,tlag,timelength,noiselevel,noisetype,round(MSEda08_low,2),round(crda08_low,2))
        resultda09_low <- c('Renyi09','rlow',bd,tlag,timelength,noiselevel,noisetype,round(MSEda09_low,2),round(crda09_low,2))
        
        resultjr2_low <- c('JensenRenyi2','rlow',bd,tlag,timelength,noiselevel,noisetype,round(MSEjr2_low,2),round(crjr2_low,2))
        resultjr3_low <- c('JensenRenyi3','rlow',bd,tlag,timelength,noiselevel,noisetype,round(MSEjr3_low,2),round(crjr3_low,2))
        resultjr4_low <- c('JensenRenyi4','rlow',bd,tlag,timelength,noiselevel,noisetype,round(MSEjr4_low,2),round(crjr4_low,2))
        resultjr5_low <- c('JensenRenyi5','rlow',bd,tlag,timelength,noiselevel,noisetype,round(MSEjr5_low,2),round(crjr5_low,2))
        resultjr6_low <- c('JensenRenyi6','rlow',bd,tlag,timelength,noiselevel,noisetype,round(MSEjr6_low,2),round(crjr6_low,2))
        resultjr7_low <- c('JensenRenyi7','rlow',bd,tlag,timelength,noiselevel,noisetype,round(MSEjr7_low,2),round(crjr7_low,2))
        resultjr8_low <- c('JensenRenyi8','rlow',bd,tlag,timelength,noiselevel,noisetype,round(MSEjr8_low,2),round(crjr8_low,2))
        resultjr9_low <- c('JensenRenyi9','rlow',bd,tlag,timelength,noiselevel,noisetype,round(MSEjr9_low,2),round(crjr9_low,2))
        resultjr10_low <- c('JensenRenyi10','rlow',bd,tlag,timelength,noiselevel,noisetype,round(MSEjr10_low,2),round(crjr10_low,2))
        #################################
        #avg
        MSEk_avg <- mean((ravg - (1-kullback[21:40]))^2) 
        crk_avg <- cor((1-kullback[21:40]),ravg)#-0.4
        MSEj_avg <- mean((ravg - (1-j[21:40]))^2)
        crj_avg <- cor(1-j[21:40],ravg)#-0.4
        MSEmi_avg <- mean((ravg - mi[21:40])^2)
        crmi_avg <- cor(mi[21:40],ravg)#0.000217
        
        MSEc_avg <- mean((ravg - crosscorrelation[21:40])^2)
        crc_avg <- cor(crosscorrelation[21:40],ravg) 
        MSEphasesyn_avg <- mean((ravg - phasesyn[21:40])^2)
        crphasesyn_avg <- cor(phasesyn[21:40],ravg)
        
        
        
        MSEda01_avg <- mean((ravg - (1-da01[21:40]))^2)
        crda01_avg <- cor((1-da01[21:40]),ravg)
        MSEda02_avg <- mean((ravg - (1-da02[21:40]))^2)
        crda02_avg <- cor(1-da02[21:40],ravg)
        MSEda03_avg <- mean((ravg - (1-da03[21:40]))^2)
        crda03_avg <- cor(1-da03[21:40],ravg)
        MSEda04_avg <- mean((ravg - (1-da04[21:40]))^2)
        crda04_avg <- cor(1-da04[21:40],ravg)
        MSEda05_avg <- mean((ravg - (1-da05[21:40]))^2)
        crda05_avg <- cor(1-da05[21:40],ravg)
        MSEda06_avg <- mean((ravg - (1-da06[21:40]))^2)
        crda06_avg <- cor(1-da06[21:40],ravg)
        MSEda07_avg <- mean((ravg - (1-da07[21:40]))^2)
        crda07_avg <- cor(1-da07[21:40],ravg)
        MSEda08_avg <- mean((ravg - (1-da08[21:40]))^2)
        crda08_avg <- cor(1-da08[21:40],ravg)
        MSEda09_avg <- mean((ravg - (1-da09[21:40]))^2)
        crda09_avg <- cor(1-da09[21:40],ravg)
        
        MSEjr2_avg <- mean((ravg - (1-jr2[21:40]))^2)
        crjr2_avg <- cor((1-jr2[21:40]),ravg)
        MSEjr3_avg <- mean((ravg - (1-jr3[21:40]))^2)
        crjr3_avg <- cor(1-jr3[21:40],ravg)
        MSEjr4_avg <- mean((ravg - (1-jr4[21:40]))^2)
        crjr4_avg <- cor(1-jr4[21:40],ravg)
        MSEjr5_avg <- mean((ravg - (1-jr5[21:40]))^2)
        crjr5_avg <- cor(1-jr5[21:40],ravg)
        MSEjr6_avg <- mean((ravg - (1-jr6[21:40]))^2)
        crjr6_avg <- cor(1-jr6[21:40],ravg)
        MSEjr7_avg <- mean((ravg - (1-jr7[21:40]))^2)
        crjr7_avg <- cor(1-jr7[21:40],ravg)
        MSEjr8_avg <- mean((ravg - (1-jr8[21:40]))^2)
        crjr8_avg <- cor(1-jr8[21:40],ravg)
        MSEjr9_avg <- mean((ravg - (1-jr9[21:40]))^2)
        crjr9_avg <- cor(1-jr9[21:40],ravg)
        MSEjr10_avg <- mean((ravg - (1-jr10[21:40]))^2)
        crjr10_avg <- cor((1-jr10[21:40]),ravg)
        
        
        
        
        resultk_avg <- c('Kullback','ravg',bd,tlag,timelength,noiselevel,noisetype,round(MSEk_avg,2),round(crk_avg,2))
        resultj_avg <- c('JensenShannon','ravg',bd,tlag,timelength,noiselevel,noisetype,round(MSEj_avg,2),round(crj_avg,2))
        resultmi_avg <- c('mutual_information','ravg',bd,tlag,timelength,noiselevel,noisetype,round(MSEmi_avg,2),round(crmi_avg,2))
        resultc_avg <- c('cross_correlation','ravg',bd,tlag,timelength,noiselevel,noisetype,round(MSEc_avg,2),round(crc_avg,2))
        resultphasesyn_avg <- c('phasesynchrony','ravg',bd,tlag,timelength,noiselevel,noisetype,round(MSEphasesyn_avg,2),round(crphasesyn_avg,2))
        resultda01_avg <- c('Renyi01','ravg',bd,tlag,timelength,noiselevel,noisetype,round(MSEda01_avg,2),round(crda01_avg,2))
        resultda02_avg <- c('Renyi02','ravg',bd,tlag,timelength,noiselevel,noisetype,round(MSEda02_avg,2),round(crda02_avg,2))
        resultda03_avg <- c('Renyi03','ravg',bd,tlag,timelength,noiselevel,noisetype,round(MSEda03_avg,2),round(crda03_avg,2))
        resultda04_avg <- c('Renyi04','ravg',bd,tlag,timelength,noiselevel,noisetype,round(MSEda04_avg,2),round(crda04_avg,2))
        resultda05_avg <- c('Renyi05','ravg',bd,tlag,timelength,noiselevel,noisetype,round(MSEda05_avg,2),round(crda05_avg,2))
        resultda06_avg <- c('Renyi06','ravg',bd,tlag,timelength,noiselevel,noisetype,round(MSEda06_avg,2),round(crda06_avg,2))
        resultda07_avg <- c('Renyi07','ravg',bd,tlag,timelength,noiselevel,noisetype,round(MSEda07_avg,2),round(crda07_avg,2))
        resultda08_avg <- c('Renyi08','ravg',bd,tlag,timelength,noiselevel,noisetype,round(MSEda08_avg,2),round(crda08_avg,2))
        resultda09_avg <- c('Renyi09','ravg',bd,tlag,timelength,noiselevel,noisetype,round(MSEda09_avg,2),round(crda09_avg,2))
        
        resultjr2_avg <- c('JensenRenyi2','ravg',bd,tlag,timelength,noiselevel,noisetype,round(MSEjr2_avg,2),round(crjr2_avg,2))
        resultjr3_avg <- c('JensenRenyi3','ravg',bd,tlag,timelength,noiselevel,noisetype,round(MSEjr3_avg,2),round(crjr3_avg,2))
        resultjr4_avg <- c('JensenRenyi4','ravg',bd,tlag,timelength,noiselevel,noisetype,round(MSEjr4_avg,2),round(crjr4_avg,2))
        resultjr5_avg <- c('JensenRenyi5','ravg',bd,tlag,timelength,noiselevel,noisetype,round(MSEjr5_avg,2),round(crjr5_avg,2))
        resultjr6_avg <- c('JensenRenyi6','ravg',bd,tlag,timelength,noiselevel,noisetype,round(MSEjr6_avg,2),round(crjr6_avg,2))
        resultjr7_avg <- c('JensenRenyi7','ravg',bd,tlag,timelength,noiselevel,noisetype,round(MSEjr7_avg,2),round(crjr7_avg,2))
        resultjr8_avg <- c('JensenRenyi8','ravg',bd,tlag,timelength,noiselevel,noisetype,round(MSEjr8_avg,2),round(crjr8_avg,2))
        resultjr9_avg <- c('JensenRenyi9','ravg',bd,tlag,timelength,noiselevel,noisetype,round(MSEjr9_avg,2),round(crjr9_avg,2))
        resultjr10_avg <- c('JensenRenyi10','ravg',bd,tlag,timelength,noiselevel,noisetype,round(MSEjr10_avg,2),round(crjr10_avg,2))
        #################################
        ###################################
        #high
        
        MSEk_high <- mean((rhigh - (1-kullback[41:60]))^2) 
        crk_high <- cor((1-kullback[41:60]),rhigh)#-0.4
        MSEj_high <- mean((rhigh - (1-j[41:60]))^2)
        crj_high <- cor(1-j[41:60],rhigh)#-0.4
        MSEmi_high <- mean((rhigh - mi[41:60])^2)
        crmi_high <- cor(mi[41:60],rhigh)#0.000217
        
        MSEc_high <- mean((rhigh - crosscorrelation[41:60])^2)
        crc_high <- cor(crosscorrelation[41:60],rhigh) 
        MSEphasesyn_high <- mean((rhigh - phasesyn[41:60])^2)
        crphasesyn_high <- cor(phasesyn[41:60],rhigh)
        
        
        
        MSEda01_high <- mean((rhigh - (1-da01[41:60]))^2)
        crda01_high <- cor((1-da01[41:60]),rhigh)
        MSEda02_high <- mean((rhigh - (1-da02[41:60]))^2)
        crda02_high <- cor(1-da02[41:60],rhigh)
        MSEda03_high <- mean((rhigh - (1-da03[41:60]))^2)
        crda03_high <- cor(1-da03[41:60],rhigh)
        MSEda04_high <- mean((rhigh - (1-da04[41:60]))^2)
        crda04_high <- cor(1-da04[41:60],rhigh)
        MSEda05_high <- mean((rhigh - (1-da05[41:60]))^2)
        crda05_high <- cor(1-da05[41:60],rhigh)
        MSEda06_high <- mean((rhigh - (1-da06[41:60]))^2)
        crda06_high <- cor(1-da06[41:60],rhigh)
        MSEda07_high <- mean((rhigh - (1-da07[41:60]))^2)
        crda07_high <- cor(1-da07[41:60],rhigh)
        MSEda08_high <- mean((rhigh - (1-da08[41:60]))^2)
        crda08_high <- cor(1-da08[41:60],rhigh)
        MSEda09_high <- mean((rhigh - (1-da09[41:60]))^2)
        crda09_high <- cor(1-da09[41:60],rhigh)
        
        MSEjr2_high <- mean((rhigh - (1-jr2[41:60]))^2)
        crjr2_high <- cor((1-jr2[41:60]),rhigh)
        MSEjr3_high <- mean((rhigh - (1-jr3[41:60]))^2)
        crjr3_high <- cor(1-jr3[41:60],rhigh)
        MSEjr4_high <- mean((rhigh - (1-jr4[41:60]))^2)
        crjr4_high <- cor(1-jr4[41:60],rhigh)
        MSEjr5_high <- mean((rhigh - (1-jr5[41:60]))^2)
        crjr5_high <- cor(1-jr5[41:60],rhigh)
        MSEjr6_high <- mean((rhigh - (1-jr6[41:60]))^2)
        crjr6_high <- cor(1-jr6[41:60],rhigh)
        MSEjr7_high <- mean((rhigh - (1-jr7[41:60]))^2)
        crjr7_high <- cor(1-jr7[41:60],rhigh)
        MSEjr8_high <- mean((rhigh - (1-jr8[41:60]))^2)
        crjr8_high <- cor(1-jr8[41:60],rhigh)
        MSEjr9_high <- mean((rhigh - (1-jr9[41:60]))^2)
        crjr9_high <- cor(1-jr9[41:60],rhigh)
        MSEjr10_high <- mean((rhigh - (1-jr10[41:60]))^2)
        crjr10_high <- cor((1-jr10[41:60]),rhigh)
        
        
        resultk_high <- c('Kullback','rhigh',bd,tlag,timelength,noiselevel,noisetype,round(MSEk_high,2),round(crk_high,2))
        resultj_high <- c('JensenShannon','rhigh',bd,tlag,timelength,noiselevel,noisetype,round(MSEj_high,2),round(crj_high,2))
        resultmi_high <- c('mutual_information','rhigh',bd,tlag,timelength,noiselevel,noisetype,round(MSEmi_high,2),round(crmi_high,2))
        resultc_high <- c('cross_correlation','rhigh',bd,tlag,timelength,noiselevel,noisetype,round(MSEc_high,2),round(crc_high,2))
        resultphasesyn_high <- c('phasesynchrony','rhigh',bd,tlag,timelength,noiselevel,noisetype,round(MSEphasesyn_high,2),round(crphasesyn_high,2))
        resultda01_high <- c('Renyi01','rhigh',bd,tlag,timelength,noiselevel,noisetype,round(MSEda01_high,2),round(crda01_high,2))
        resultda02_high <- c('Renyi02','rhigh',bd,tlag,timelength,noiselevel,noisetype,round(MSEda02_high,2),round(crda02_high,2))
        resultda03_high <- c('Renyi03','rhigh',bd,tlag,timelength,noiselevel,noisetype,round(MSEda03_high,2),round(crda03_high,2))
        resultda04_high <- c('Renyi04','rhigh',bd,tlag,timelength,noiselevel,noisetype,round(MSEda04_high,2),round(crda04_high,2))
        resultda05_high <- c('Renyi05','rhigh',bd,tlag,timelength,noiselevel,noisetype,round(MSEda05_high,2),round(crda05_high,2))
        resultda06_high <- c('Renyi06','rhigh',bd,tlag,timelength,noiselevel,noisetype,round(MSEda06_high,2),round(crda06_high,2))
        resultda07_high <- c('Renyi07','rhigh',bd,tlag,timelength,noiselevel,noisetype,round(MSEda07_high,2),round(crda07_high,2))
        resultda08_high <- c('Renyi08','rhigh',bd,tlag,timelength,noiselevel,noisetype,round(MSEda08_high,2),round(crda08_high,2))
        resultda09_high <- c('Renyi09','rhigh',bd,tlag,timelength,noiselevel,noisetype,round(MSEda09_high,2),round(crda09_high,2))
        
        resultjr2_high <- c('JensenRenyi2','rhigh',bd,tlag,timelength,noiselevel,noisetype,round(MSEjr2_high,2),round(crjr2_high,2))
        resultjr3_high <- c('JensenRenyi3','rhigh',bd,tlag,timelength,noiselevel,noisetype,round(MSEjr3_high,2),round(crjr3_high,2))
        resultjr4_high <- c('JensenRenyi4','rhigh',bd,tlag,timelength,noiselevel,noisetype,round(MSEjr4_high,2),round(crjr4_high,2))
        resultjr5_high <- c('JensenRenyi5','rhigh',bd,tlag,timelength,noiselevel,noisetype,round(MSEjr5_high,2),round(crjr5_high,2))
        resultjr6_high <- c('JensenRenyi6','rhigh',bd,tlag,timelength,noiselevel,noisetype,round(MSEjr6_high,2),round(crjr6_high,2))
        resultjr7_high <- c('JensenRenyi7','rhigh',bd,tlag,timelength,noiselevel,noisetype,round(MSEjr7_high,2),round(crjr7_high,2))
        resultjr8_high <- c('JensenRenyi8','rhigh',bd,tlag,timelength,noiselevel,noisetype,round(MSEjr8_high,2),round(crjr8_high,2))
        resultjr9_high <- c('JensenRenyi9','rhigh',bd,tlag,timelength,noiselevel,noisetype,round(MSEjr9_high,2),round(crjr9_high,2))
        resultjr10_high <- c('JensenRenyi10','rhigh',bd,tlag,timelength,noiselevel,noisetype,round(MSEjr10_high,2),round(crjr10_high,2))
        
        
        
        ###########################################################
        #overall
        MSEk <- mean((r - (1-kullback))^2) 
        crk <- cor((1-kullback),r)
        MSEj <- mean((r - (1-j))^2)
        crj <- cor(1-j,r)
        MSEmi <- mean((r - mi)^2)
        crmi <- cor(mi,r)
        
        MSEc <- mean((r - crosscorrelation)^2)
        crc <- cor(crosscorrelation,r) 
        MSEphasesyn <- mean((r - phasesyn)^2)
        crphasesyn <- cor(phasesyn,r)
        
        
        
        MSEda01 <- mean((r - (1-da01))^2)
        crda01 <- cor(1-da01,r)
        MSEda02 <- mean((r - (1-da02))^2)
        crda02 <- cor(1-da02,r)
        MSEda03 <- mean((r - (1-da03))^2)
        crda03 <- cor(1-da03,r)
        MSEda04 <- mean((r - (1-da04))^2)
        crda04 <- cor(1-da04,r)
        MSEda05 <- mean((r - (1-da05))^2)
        crda05 <- cor(1-da05,r)
        MSEda06 <- mean((r - (1-da06))^2)
        crda06 <- cor(1-da06,r)
        MSEda07 <- mean((r - (1-da07))^2)
        crda07 <- cor(1-da07,r)
        MSEda08 <- mean((r - (1-da08))^2)
        crda08 <- cor(1-da08,r)
        MSEda09 <- mean((r - (1-da09))^2)
        crda09 <- cor(1-da09,r)
        
        MSEjr2 <- mean((r - (1-jr2))^2)
        crjr2 <- cor((1-jr2),r)
        MSEjr3 <- mean((r - (1-jr3))^2)
        crjr3 <- cor(1-jr3,r)
        MSEjr4 <- mean((r - (1-jr4))^2)
        crjr4 <- cor(1-jr4,r)
        MSEjr5 <- mean((r - (1-jr5))^2)
        crjr5 <- cor(1-jr5,r)
        MSEjr6 <- mean((r - (1-jr6))^2)
        crjr6 <- cor(1-jr6,r)
        MSEjr7 <- mean((r - (1-jr7))^2)
        crjr7 <- cor(1-jr7,r)
        MSEjr8 <- mean((r - (1-jr8))^2)
        crjr8 <- cor(1-jr8,r)
        MSEjr9 <- mean((r - (1-jr9))^2)
        crjr9 <- cor(1-jr9,r)
        MSEjr10 <- mean((r - (1-jr10))^2)
        crjr10 <- cor((1-jr10),r)
        
        
        resultk <- c('Kullback','r',bd,tlag,timelength,noiselevel,noisetype,round(MSEk,2),round(crk,2))
        resultj <- c('JensenShannon','r',bd,tlag,timelength,noiselevel,noisetype,round(MSEj,2),round(crj,2))
        resultmi <- c('mutual_information','r',bd,tlag,timelength,noiselevel,noisetype,round(MSEmi,2),round(crmi,2))
        resultc <- c('cross_correlation','r',bd,tlag,timelength,noiselevel,noisetype,round(MSEc,2),round(crc,2))
        resultphasesyn <- c('phasesynchrony','r',bd,tlag,timelength,noiselevel,noisetype,round(MSEphasesyn,2),round(crphasesyn,2))
        resultda01 <- c('Renyi01','r',bd,tlag,timelength,noiselevel,noisetype,round(MSEda01,2),round(crda01,2))
        resultda02 <- c('Renyi02','r',bd,tlag,timelength,noiselevel,noisetype,round(MSEda02,2),round(crda02,2))
        resultda03 <- c('Renyi03','r',bd,tlag,timelength,noiselevel,noisetype,round(MSEda03,2),round(crda03,2))
        resultda04 <- c('Renyi04','r',bd,tlag,timelength,noiselevel,noisetype,round(MSEda04,2),round(crda04,2))
        resultda05 <- c('Renyi05','r',bd,tlag,timelength,noiselevel,noisetype,round(MSEda05,2),round(crda05,2))
        resultda06 <- c('Renyi06','r',bd,tlag,timelength,noiselevel,noisetype,round(MSEda06,2),round(crda06,2))
        resultda07 <- c('Renyi07','r',bd,tlag,timelength,noiselevel,noisetype,round(MSEda07,2),round(crda07,2))
        resultda08 <- c('Renyi08','r',bd,tlag,timelength,noiselevel,noisetype,round(MSEda08,2),round(crda08,2))
        resultda09 <- c('Renyi09','r',bd,tlag,timelength,noiselevel,noisetype,round(MSEda09,2),round(crda09,2))
        
        resultjr2 <- c('JensenRenyi2','r',bd,tlag,timelength,noiselevel,noisetype,round(MSEjr2,2),round(crjr2,2))
        resultjr3 <- c('JensenRenyi3','r',bd,tlag,timelength,noiselevel,noisetype,round(MSEjr3,2),round(crjr3,2))
        resultjr4 <- c('JensenRenyi4','r',bd,tlag,timelength,noiselevel,noisetype,round(MSEjr4,2),round(crjr4,2))
        resultjr5 <- c('JensenRenyi5','r',bd,tlag,timelength,noiselevel,noisetype,round(MSEjr5,2),round(crjr5,2))
        resultjr6 <- c('JensenRenyi6','r',bd,tlag,timelength,noiselevel,noisetype,round(MSEjr6,2),round(crjr6,2))
        resultjr7 <- c('JensenRenyi7','r',bd,tlag,timelength,noiselevel,noisetype,round(MSEjr7,2),round(crjr7,2))
        resultjr8 <- c('JensenRenyi8','r',bd,tlag,timelength,noiselevel,noisetype,round(MSEjr8,2),round(crjr8,2))
        resultjr9 <- c('JensenRenyi9','r',bd,tlag,timelength,noiselevel,noisetype,round(MSEjr9,2),round(crjr9,2))
        resultjr10 <- c('JensenRenyi10','r',bd,tlag,timelength,noiselevel,noisetype,round(MSEjr10,2),round(crjr10,2))
        
        
        
        
        #######################################
        datasave <- rbind(datasave,resultk_low,resultj_low,resultmi_low,resultc_low,
                          resultphasesyn_low,resultda01_low,resultda02_low,resultda03_low,
                          resultda04_low,resultda05_low,resultda06_low,resultda07_low,
                          resultda08_low,resultda09_low,resultjr2_low,resultjr3_low,
                          resultjr4_low,resultjr5_low,resultjr6_low,
                          resultjr7_low,resultjr8_low,resultjr9_low,resultjr10_low,
                          resultk_avg,resultj_avg,resultmi_avg,resultc_avg,
                          resultphasesyn_avg,resultda01_avg,resultda02_avg,resultda03_avg,
                          resultda04_avg,resultda05_avg,resultda06_avg,resultda07_avg,
                          resultda08_avg,resultda09_avg,resultjr2_avg,resultjr3_avg,
                          resultjr4_avg,resultjr5_avg,resultjr6_avg,
                          resultjr7_avg,resultjr8_avg,resultjr9_avg,resultjr10_avg,
                          resultk_high,resultj_high,resultmi_high,resultc_high,
                          resultphasesyn_high,resultda01_high,resultda02_high,resultda03_high,
                          resultda04_high,resultda05_high,resultda06_high,resultda07_high,
                          resultda08_high,resultda09_high,resultjr2_high,resultjr3_high,
                          resultjr4_high,resultjr5_high,resultjr6_high,
                          resultjr7_high,resultjr8_high,resultjr9_high,resultjr10_high
        )
        
        datasave.all <- rbind(datasave.all,resultk,resultj,resultmi,resultc,resultphasesyn,resultda01,
                              resultda02,resultda03,resultda04,resultda05,resultda06,resultda07,
                              resultda08,resultda09,resultjr2,resultjr3,resultjr4,resultjr5,resultjr6,
                              resultjr7,resultjr8,resultjr9,resultjr10)
        
        
        xsave[[times]] <- xnoise
        
        
        ysave[[times]] <- ynoise
        
        
        
        crosscorrelationsave[[times]] <- crosscorrelation
        phasesynsave[[times]] <- phasesyn
        misave[[times]] <- mi
        jsave[[times]] <- j
        kullbacksave[[times]] <- kullback
        dasave[[times]] <- da
        jrsave[[times]] <- jr 
        
        estimates$crosscorrelation <- crosscorrelationsave
        estimates$phasesyn <- phasesynsave
        estimates$mi <- misave
        estimates$j <- jsave
        estimates$kullback <- kullbacksave
        estimates$da <- dasave
        estimates$jr <- jrsave
        
      }
      timeend <- Sys.time()
      
      save(datasave,file = paste(c('bd',bd,'_lag',tlag,'_timelength',timelength,'_noiselevel',noiselevel,'_noisetype',noisetype,'.RData'),sep = '',collapse = ''))
      save(datasave.all,file = paste(c('overall_','bd',bd,'_lag',tlag,'_timelength',timelength,'_noiselevel',noiselevel,'_noisetype',noisetype,'.RData'),sep = '',collapse = ''))
      save(xsave,file = paste(c('timeseriesx','bd',bd,'_lag',tlag,'_timelength',timelength,'_noiselevel',noiselevel,'_noisetype',noisetype,'.RData'),sep = '',collapse = '')
      )
      save(ysave,file = paste(c('timeseriesy','bd',bd,'_lag',tlag,'_timelength',timelength,'_noiselevel',noiselevel,'_noisetype',noisetype,'.RData'),sep = '',collapse = '')
      )
      save(estimates,file=paste(c('estimates','bd',bd,'_lag',tlag,'_timelength',timelength,'_noiselevel',noiselevel,'_noisetype',noisetype,'.RData'),sep = '',collapse = '')
      )
      timeend-timestart
    }
  }
}
