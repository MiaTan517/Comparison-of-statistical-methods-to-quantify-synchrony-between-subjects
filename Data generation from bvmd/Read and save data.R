# Run this file after generating the data by running the other two R files in this folder

timelength_sets <- c(1006,10054)
kappavalue_sets <- c(0.5,1,2)
kappa_sets <- c(0.5,1,2)

noiselevel_sets <- c(0.2,0.5)

noisetype_sets <- c('white','ar0.5') 
#phivalue <- 0.5 # ar correlation

tlag_sets <- c(0,2) #change manually
bvmddata <- bvmddatac <-  c()
##############################################################################################

for (timelength in timelength_sets) {
  for (kappavalue in kappavalue_sets) {
    for (noiselevel in noiselevel_sets) {
      for (tlag in tlag_sets) {
        for (noisetype in noisetype_sets) {
          load(paste(c('kappa',kappavalue,'_lag',tlag,'_timelength',timelength,'_noiselevel',noiselevel,'_noisetype',noisetype,'.RData'),sep = '',collapse = ''))
          load(paste(c('ckappa',kappavalue,'_lag',tlag,'_timelength',timelength,'_noiselevel',noiselevel,'_noisetype',noisetype,'.RData'),sep = '',collapse = ''))
          bvmddata <- rbind(bvmddata,datasave)
          bvmddatac <- rbind(bvmddatac,datasavec)
          
        }
      }
    }
  }
}
bvmddata <- cbind(bvmddata,bvmddatac[,8:9])
save(bvmddata,file = 'bvmddata.RData')


load('bvmddata.RData')

###################################################################
bvmddata.all <- c()
bvmddata.allc <- c()
##############################################################################################

for (timelength in timelength_sets) {
  for (kappavalue in kappavalue_sets) {
    for (noiselevel in noiselevel_sets) {
      for (tlag in tlag_sets) {
        for (noisetype in noisetype_sets) {
          load(paste(c('overall_','kappa',kappavalue,'_lag',tlag,'_timelength',timelength,'_noiselevel',noiselevel,'_noisetype',noisetype,'.RData'),sep = '',collapse = ''))
          load(paste(c('coverall_','kappa',kappavalue,'_lag',tlag,'_timelength',timelength,'_noiselevel',noiselevel,'_noisetype',noisetype,'.RData'),sep = '',collapse = ''))
          bvmddata.all <- rbind(bvmddata.all,datasave.all)
          bvmddata.allc <- rbind(bvmddata.allc,datasave.allc)
          
        }
      }
    }
  }
}
bvmddata.all <- cbind(bvmddata.all,bvmddata.allc[,8:9])
save(bvmddata.all,file = 'bvmddata.all.RData')
##################################################################


bvmddata <- data.frame(bvmddata)
colnames(bvmddata) <- c("synchr_method","true_synchrony","kappa",     
                        "tlag", "timelength", "noiselevel",
                        "noisetype", "MSE",'correlation','MSE_c','correlation_c')
bvmddata$true_synchrony <-  factor(bvmddata$true_synchrony)
bvmddata$kappa <- factor(bvmddata$kappa)
bvmddata$tlag <- factor(bvmddata$tlag)
bvmddata$timelength <- factor(bvmddata$timelength)
bvmddata$noiselevel <- factor(bvmddata$noiselevel)
bvmddata$noisetype <- factor(bvmddata$noisetype)
bvmddata$MSE <- as.numeric(bvmddata$MSE)
bvmddata$correlation <- as.numeric(bvmddata$correlation)
bvmddata$MSE_c <- as.numeric(bvmddata$MSE_c)
bvmddata$correlation_c <- as.numeric(bvmddata$correlation_c)
bvmddata$synchr_method <- factor(bvmddata$synchr_method)

id <- rep(1:2880,each=11)
bvmddata <- cbind(id,bvmddata)
save(bvmddata,file='bvmddata.RData')
load('bvmddata.RData')
