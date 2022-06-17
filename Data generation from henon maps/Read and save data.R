timelength_sets <- c(1006,10054)
bd_sets <- c('33','31')
noiselevel_sets <- c(0.2,0.5)

noisetype_sets <- c('white','ar') #change manually
#phivalue <- 0.5 # ar correlation

tlag_sets <- c(0,2) #change manually
henondata <- henondatac <- c()
##############################################################################################

for (timelength in timelength_sets) {
  for (bd in bd_sets) {
    for (noiselevel in noiselevel_sets) {
      for (tlag in tlag_sets) {
        for (noisetype in noisetype_sets) {
          load(paste(c('bd',bd,'_lag',tlag,'_timelength',timelength,'_noiselevel',noiselevel,'_noisetype',noisetype,'.RData'),sep = '',collapse = ''))
          load(paste(c('c','bd',bd,'_lag',tlag,'_timelength',timelength,'_noiselevel',noiselevel,'_noisetype',noisetype,'.RData'),sep = '',collapse = ''))
          henondata <- rbind(henondata,datasave)
          henondatac <- rbind(henondatac,datasavec)
          
        }
      }
    }
  }
}
henondata <- cbind(henondata,henondatac[,8:9])
save(henondata,file = 'henondata.RData')
#####################################################
henondata.all <- c()
henondata.allc <- c()
##############################################################################################

for (timelength in timelength_sets) {
  for (bd in bd_sets) {
    for (noiselevel in noiselevel_sets) {
      for (tlag in tlag_sets) {
        for (noisetype in noisetype_sets) {
          load(paste(c('overall_','bd',bd,'_lag',tlag,'_timelength',timelength,'_noiselevel',noiselevel,'_noisetype',noisetype,'.RData'),sep = '',collapse = ''))
          load(paste(c('c','overall_','bd',bd,'_lag',tlag,'_timelength',timelength,'_noiselevel',noiselevel,'_noisetype',noisetype,'.RData'),sep = '',collapse = ''))
          
          henondata.all <- rbind(henondata.all,datasave.all)
          henondata.allc <- rbind(henondata.allc,datasave.allc)
          
        }
      }
    }
  }
}
henondata.all <- cbind(henondata.all,henondata.allc[,8:9])
save(henondata.all,file = 'henondata.all.RData')

####################################################################################################
load('henondata.RData')
datahenon <- data.frame(henondata)
colnames(datahenon) <- c("synchr_method","true_synchrony","bd",     
                         "tlag", "timelength", "noiselevel",
                         "noisetype", "MSE",'correlation','MSE_c','correlation_c')
datahenon$true_synchrony <-  factor(datahenon$true_synchrony)
datahenon$bd <- factor(datahenon$bd)
datahenon$tlag <- factor(datahenon$tlag)
datahenon$timelength <- factor(datahenon$timelength)
datahenon$noiselevel <- factor(datahenon$noiselevel)
datahenon$noisetype <- factor(datahenon$noisetype)
datahenon$MSE <- as.numeric(datahenon$MSE)
datahenon$correlation <- as.numeric(datahenon$correlation)
datahenon$MSE_c <- as.numeric(datahenon$MSE_c)
datahenon$correlation_c <- as.numeric(datahenon$correlation_c)
datahenon$synchr_method <- factor(datahenon$synchr_method)

id <- rep(1:1920,each=11)
datahenon <- cbind(id,datahenon)

save(datahenon,file = 'datahenon.RData')
load('datahenon.RData')

