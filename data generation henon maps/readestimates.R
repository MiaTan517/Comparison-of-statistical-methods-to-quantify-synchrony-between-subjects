


timelength_sets <- c(1006,10054)
bd_sets <- c('33','31')
noiselevel_sets <- c(0.2,0.5)

noisetype_sets <- c('white','ar') #change manually
#phivalue <- 0.5 # ar correlation

tlag_sets <- c(0,2) #change manually

#########################################################
#phasesyn
datalow <- NULL


for (timelength in timelength_sets) {
  for (bd in bd_sets) {
    for (noiselevel in noiselevel_sets) {
      for (tlag in tlag_sets) {
        for (noisetype in noisetype_sets) {
          load(paste(c('estimates','bd',bd,'_lag',tlag,'_timelength',timelength,'_noiselevel',noiselevel,'_noisetype',noisetype,'.RData'),sep = '',collapse = ''))
          for (m in 1:20) {
            for (n in 1:20) {
              datalow <- rbind(datalow,cbind('phasesynchrony','rlow',bd,tlag,timelength,noiselevel,noisetype,round(estimates$phasesyn[[m]][n],2),round(estimates$true_synchrony_mu[[m]][n],2),round(estimates$crosscorrelation[[m]][n],2)))
            }
          }
        }
      }
    }
  }
}
dataavg <- NULL

for (timelength in timelength_sets) {
  for (bd in bd_sets) {
    for (noiselevel in noiselevel_sets) {
      for (tlag in tlag_sets) {
        for (noisetype in noisetype_sets) {
          load(paste(c('estimates','bd',bd,'_lag',tlag,'_timelength',timelength,'_noiselevel',noiselevel,'_noisetype',noisetype,'.RData'),sep = '',collapse = ''))
          for (m in 1:20) {
            for (n in 21:40) {
              dataavg <- rbind(dataavg,cbind('phasesynchrony','ravg',bd,tlag,timelength,noiselevel,noisetype,round(estimates$phasesyn[[m]][n],2),round(estimates$true_synchrony_mu[[m]][n],2),round(estimates$crosscorrelation[[m]][n],2)))
            }
          }
        }
      }
    }
  }
}

datahigh <- NULL
for (timelength in timelength_sets) {
  for (bd in bd_sets) {
    for (noiselevel in noiselevel_sets) {
      for (tlag in tlag_sets) {
        for (noisetype in noisetype_sets) {
          load(paste(c('estimates','bd',bd,'_lag',tlag,'_timelength',timelength,'_noiselevel',noiselevel,'_noisetype',noisetype,'.RData'),sep = '',collapse = ''))
          for (m in 1:20) {
            for (n in 41:60) {
              datahigh <- rbind(datahigh,cbind('phasesynchrony','rhigh',bd,tlag,timelength,noiselevel,noisetype,round(estimates$phasesyn[[m]][n],2),round(estimates$true_synchrony_mu[[m]][n],2),round(estimates$crosscorrelation[[m]][n],2)))
            }
          }
        }
      }
    }
  }
}

phasesyn <- rbind(datalow,dataavg,datahigh)
##################################
#kullback
datalow <- NULL


for (timelength in timelength_sets) {
  for (bd in bd_sets) {
    for (noiselevel in noiselevel_sets) {
      for (tlag in tlag_sets) {
        for (noisetype in noisetype_sets) {
          load(paste(c('estimates','bd',bd,'_lag',tlag,'_timelength',timelength,'_noiselevel',noiselevel,'_noisetype',noisetype,'.RData'),sep = '',collapse = ''))
          for (m in 1:20) {
            for (n in 1:20) {
              datalow <- rbind(datalow,cbind('kullback','rlow',bd,tlag,timelength,noiselevel,noisetype,round(estimates$kullback[[m]][n],2),round(estimates$true_synchrony_mu[[m]][n],2),round(estimates$crosscorrelation[[m]][n],2)))
            }
          }
        }
      }
    }
  }
}
dataavg <- NULL

for (timelength in timelength_sets) {
  for (bd in bd_sets) {
    for (noiselevel in noiselevel_sets) {
      for (tlag in tlag_sets) {
        for (noisetype in noisetype_sets) {
          load(paste(c('estimates','bd',bd,'_lag',tlag,'_timelength',timelength,'_noiselevel',noiselevel,'_noisetype',noisetype,'.RData'),sep = '',collapse = ''))
          for (m in 1:20) {
            for (n in 21:40) {
              dataavg <- rbind(dataavg,cbind('kullback','ravg',bd,tlag,timelength,noiselevel,noisetype,round(estimates$kullback[[m]][n],2),round(estimates$true_synchrony_mu[[m]][n],2),round(estimates$crosscorrelation[[m]][n],2)))
            }
          }
        }
      }
    }
  }
}

datahigh <- NULL
for (timelength in timelength_sets) {
  for (bd in bd_sets) {
    for (noiselevel in noiselevel_sets) {
      for (tlag in tlag_sets) {
        for (noisetype in noisetype_sets) {
          load(paste(c('estimates','bd',bd,'_lag',tlag,'_timelength',timelength,'_noiselevel',noiselevel,'_noisetype',noisetype,'.RData'),sep = '',collapse = ''))
          for (m in 1:20) {
            for (n in 41:60) {
              datahigh <- rbind(datahigh,cbind('kullback','rhigh',bd,tlag,timelength,noiselevel,noisetype,round(estimates$kullback[[m]][n],2),round(estimates$true_synchrony_mu[[m]][n],2),round(estimates$crosscorrelation[[m]][n],2)))
            }
          }
        }
      }
    }
  }
}

Kullback <- rbind(datalow,dataavg,datahigh)
#######################################################
#mutualinformation
datalow <- NULL


for (timelength in timelength_sets) {
  for (bd in bd_sets) {
    for (noiselevel in noiselevel_sets) {
      for (tlag in tlag_sets) {
        for (noisetype in noisetype_sets) {
          load(paste(c('estimates','bd',bd,'_lag',tlag,'_timelength',timelength,'_noiselevel',noiselevel,'_noisetype',noisetype,'.RData'),sep = '',collapse = ''))
          for (m in 1:20) {
            for (n in 1:20) {
              datalow <- rbind(datalow,cbind('mutual_information','rlow',bd,tlag,timelength,noiselevel,noisetype,round(estimates$mi[[m]][n],2),round(estimates$true_synchrony_mu[[m]][n],2),round(estimates$crosscorrelation[[m]][n],2)))
            }
          }
        }
      }
    }
  }
}
dataavg <- NULL

for (timelength in timelength_sets) {
  for (bd in bd_sets) {
    for (noiselevel in noiselevel_sets) {
      for (tlag in tlag_sets) {
        for (noisetype in noisetype_sets) {
          load(paste(c('estimates','bd',bd,'_lag',tlag,'_timelength',timelength,'_noiselevel',noiselevel,'_noisetype',noisetype,'.RData'),sep = '',collapse = ''))
          for (m in 1:20) {
            for (n in 21:40) {
              dataavg <- rbind(dataavg,cbind('mutual_information','ravg',bd,tlag,timelength,noiselevel,noisetype,round(estimates$mi[[m]][n],2),round(estimates$true_synchrony_mu[[m]][n],2),round(estimates$crosscorrelation[[m]][n],2)))
            }
          }
        }
      }
    }
  }
}

datahigh <- NULL
for (timelength in timelength_sets) {
  for (bd in bd_sets) {
    for (noiselevel in noiselevel_sets) {
      for (tlag in tlag_sets) {
        for (noisetype in noisetype_sets) {
          load(paste(c('estimates','bd',bd,'_lag',tlag,'_timelength',timelength,'_noiselevel',noiselevel,'_noisetype',noisetype,'.RData'),sep = '',collapse = ''))
          for (m in 1:20) {
            for (n in 41:60) {
              datahigh <- rbind(datahigh,cbind('mutual_information','rhigh',bd,tlag,timelength,noiselevel,noisetype,round(estimates$mi[[m]][n],2),round(estimates$true_synchrony_mu[[m]][n],2),round(estimates$crosscorrelation[[m]][n],2)))
            }
          }
        }
      }
    }
  }
}

mutual_information <- rbind(datalow,dataavg,datahigh)
#####################################################

##########################################
# JensenShannon
estimates$j
datalow <- NULL


for (timelength in timelength_sets) {
  for (bd in bd_sets) {
    for (noiselevel in noiselevel_sets) {
      for (tlag in tlag_sets) {
        for (noisetype in noisetype_sets) {
          load(paste(c('estimates','bd',bd,'_lag',tlag,'_timelength',timelength,'_noiselevel',noiselevel,'_noisetype',noisetype,'.RData'),sep = '',collapse = ''))
          for (m in 1:20) {
            for (n in 1:20) {
              datalow <- rbind(datalow,cbind('JensenShannon','rlow',bd,tlag,timelength,noiselevel,noisetype,round(estimates$j[[m]][n],2),round(estimates$true_synchrony_mu[[m]][n],2),round(estimates$crosscorrelation[[m]][n],2)))
            }
          }
        }
      }
    }
  }
}
dataavg <- NULL

for (timelength in timelength_sets) {
  for (bd in bd_sets) {
    for (noiselevel in noiselevel_sets) {
      for (tlag in tlag_sets) {
        for (noisetype in noisetype_sets) {
          load(paste(c('estimates','bd',bd,'_lag',tlag,'_timelength',timelength,'_noiselevel',noiselevel,'_noisetype',noisetype,'.RData'),sep = '',collapse = ''))
          for (m in 1:20) {
            for (n in 21:40) {
              dataavg <- rbind(dataavg,cbind('JensenShannon','ravg',bd,tlag,timelength,noiselevel,noisetype,round(estimates$j[[m]][n],2),round(estimates$true_synchrony_mu[[m]][n],2),round(estimates$crosscorrelation[[m]][n],2)))
            }
          }
        }
      }
    }
  }
}

datahigh <- NULL
for (timelength in timelength_sets) {
  for (bd in bd_sets) {
    for (noiselevel in noiselevel_sets) {
      for (tlag in tlag_sets) {
        for (noisetype in noisetype_sets) {
          load(paste(c('estimates','bd',bd,'_lag',tlag,'_timelength',timelength,'_noiselevel',noiselevel,'_noisetype',noisetype,'.RData'),sep = '',collapse = ''))
          for (m in 1:20) {
            for (n in 41:60) {
              datahigh <- rbind(datahigh,cbind('JensenShannon','rhigh',bd,tlag,timelength,noiselevel,noisetype,round(estimates$j[[m]][n],2),round(estimates$true_synchrony_mu[[m]][n],2),round(estimates$crosscorrelation[[m]][n],2)))
            }
          }
        }
      }
    }
  }
}

JensenShannon <- rbind(datalow,dataavg,datahigh)
################################################################
datalow <- NULL


for (timelength in timelength_sets) {
  for (bd in bd_sets) {
    for (noiselevel in noiselevel_sets) {
      for (tlag in tlag_sets) {
        for (noisetype in noisetype_sets) {
          load(paste(c('estimates','bd',bd,'_lag',tlag,'_timelength',timelength,'_noiselevel',noiselevel,'_noisetype',noisetype,'.RData'),sep = '',collapse = ''))
          for (m in 1:20) {
            for (n in 1:20) {
              datalow <- rbind(datalow,cbind('Renyi01','rlow',bd,tlag,timelength,noiselevel,noisetype,round(estimates$da[[m]][[n]][1],2),round(estimates$true_synchrony_mu[[m]][n],2),round(estimates$crosscorrelation[[m]][n],2)))
            }
          }
        }
      }
    }
  }
}
dataavg <- NULL

for (timelength in timelength_sets) {
  for (bd in bd_sets) {
    for (noiselevel in noiselevel_sets) {
      for (tlag in tlag_sets) {
        for (noisetype in noisetype_sets) {
          load(paste(c('estimates','bd',bd,'_lag',tlag,'_timelength',timelength,'_noiselevel',noiselevel,'_noisetype',noisetype,'.RData'),sep = '',collapse = ''))
          for (m in 1:20) {
            for (n in 21:40) {
              dataavg <- rbind(dataavg,cbind('Renyi01','ravg',bd,tlag,timelength,noiselevel,noisetype,round(estimates$da[[m]][[n]][1],2),round(estimates$true_synchrony_mu[[m]][n],2),round(estimates$crosscorrelation[[m]][n],2)))
            }
          }
        }
      }
    }
  }
}

datahigh <- NULL
for (timelength in timelength_sets) {
  for (bd in bd_sets) {
    for (noiselevel in noiselevel_sets) {
      for (tlag in tlag_sets) {
        for (noisetype in noisetype_sets) {
          load(paste(c('estimates','bd',bd,'_lag',tlag,'_timelength',timelength,'_noiselevel',noiselevel,'_noisetype',noisetype,'.RData'),sep = '',collapse = ''))
          for (m in 1:20) {
            for (n in 41:60) {
              datahigh <- rbind(datahigh,cbind('Renyi01','rhigh',bd,tlag,timelength,noiselevel,noisetype,round(estimates$da[[m]][[n]][1],2),round(estimates$true_synchrony_mu[[m]][n],2),round(estimates$crosscorrelation[[m]][n],2)))
            }
          }
        }
      }
    }
  }
}

Renyi01 <- rbind(datalow,dataavg,datahigh)
#################################
datalow <- NULL


for (timelength in timelength_sets) {
  for (bd in bd_sets) {
    for (noiselevel in noiselevel_sets) {
      for (tlag in tlag_sets) {
        for (noisetype in noisetype_sets) {
          load(paste(c('estimates','bd',bd,'_lag',tlag,'_timelength',timelength,'_noiselevel',noiselevel,'_noisetype',noisetype,'.RData'),sep = '',collapse = ''))
          for (m in 1:20) {
            for (n in 1:20) {
              datalow <- rbind(datalow,cbind('Renyi05','rlow',bd,tlag,timelength,noiselevel,noisetype,round(estimates$da[[m]][[n]][5],2),round(estimates$true_synchrony_mu[[m]][n],2),round(estimates$crosscorrelation[[m]][n],2)))
            }
          }
        }
      }
    }
  }
}
dataavg <- NULL

for (timelength in timelength_sets) {
  for (bd in bd_sets) {
    for (noiselevel in noiselevel_sets) {
      for (tlag in tlag_sets) {
        for (noisetype in noisetype_sets) {
          load(paste(c('estimates','bd',bd,'_lag',tlag,'_timelength',timelength,'_noiselevel',noiselevel,'_noisetype',noisetype,'.RData'),sep = '',collapse = ''))
          for (m in 1:20) {
            for (n in 21:40) {
              dataavg <- rbind(dataavg,cbind('Renyi05','ravg',bd,tlag,timelength,noiselevel,noisetype,round(estimates$da[[m]][[n]][5],2),round(estimates$true_synchrony_mu[[m]][n],2),round(estimates$crosscorrelation[[m]][n],2)))
            }
          }
        }
      }
    }
  }
}

datahigh <- NULL
for (timelength in timelength_sets) {
  for (bd in bd_sets) {
    for (noiselevel in noiselevel_sets) {
      for (tlag in tlag_sets) {
        for (noisetype in noisetype_sets) {
          load(paste(c('estimates','bd',bd,'_lag',tlag,'_timelength',timelength,'_noiselevel',noiselevel,'_noisetype',noisetype,'.RData'),sep = '',collapse = ''))
          for (m in 1:20) {
            for (n in 41:60) {
              datahigh <- rbind(datahigh,cbind('Renyi05','rhigh',bd,tlag,timelength,noiselevel,noisetype,round(estimates$da[[m]][[n]][5],2),round(estimates$true_synchrony_mu[[m]][n],2),round(estimates$crosscorrelation[[m]][n],2)))
            }
          }
        }
      }
    }
  }
}

Renyi05 <- rbind(datalow,dataavg,datahigh)
#################################################################
datalow <- NULL


for (timelength in timelength_sets) {
  for (bd in bd_sets) {
    for (noiselevel in noiselevel_sets) {
      for (tlag in tlag_sets) {
        for (noisetype in noisetype_sets) {
          load(paste(c('estimates','bd',bd,'_lag',tlag,'_timelength',timelength,'_noiselevel',noiselevel,'_noisetype',noisetype,'.RData'),sep = '',collapse = ''))
          for (m in 1:20) {
            for (n in 1:20) {
              datalow <- rbind(datalow,cbind('Renyi09','rlow',bd,tlag,timelength,noiselevel,noisetype,round(estimates$da[[m]][[n]][9],2),round(estimates$true_synchrony_mu[[m]][n],2),round(estimates$crosscorrelation[[m]][n],2)))
            }
          }
        }
      }
    }
  }
}
dataavg <- NULL

for (timelength in timelength_sets) {
  for (bd in bd_sets) {
    for (noiselevel in noiselevel_sets) {
      for (tlag in tlag_sets) {
        for (noisetype in noisetype_sets) {
          load(paste(c('estimates','bd',bd,'_lag',tlag,'_timelength',timelength,'_noiselevel',noiselevel,'_noisetype',noisetype,'.RData'),sep = '',collapse = ''))
          for (m in 1:20) {
            for (n in 21:40) {
              dataavg <- rbind(dataavg,cbind('Renyi09','ravg',bd,tlag,timelength,noiselevel,noisetype,round(estimates$da[[m]][[n]][9],2),round(estimates$true_synchrony_mu[[m]][n],2),round(estimates$crosscorrelation[[m]][n],2)))
            }
          }
        }
      }
    }
  }
}

datahigh <- NULL
for (timelength in timelength_sets) {
  for (bd in bd_sets) {
    for (noiselevel in noiselevel_sets) {
      for (tlag in tlag_sets) {
        for (noisetype in noisetype_sets) {
          load(paste(c('estimates','bd',bd,'_lag',tlag,'_timelength',timelength,'_noiselevel',noiselevel,'_noisetype',noisetype,'.RData'),sep = '',collapse = ''))
          for (m in 1:20) {
            for (n in 41:60) {
              datahigh <- rbind(datahigh,cbind('Renyi09','rhigh',bd,tlag,timelength,noiselevel,noisetype,round(estimates$da[[m]][[n]][9],2),round(estimates$true_synchrony_mu[[m]][n],2),round(estimates$crosscorrelation[[m]][n],2)))
            }
          }
        }
      }
    }
  }
}

Renyi09 <- rbind(datalow,dataavg,datahigh)
##############################################################
datalow <- NULL


for (timelength in timelength_sets) {
  for (bd in bd_sets) {
    for (noiselevel in noiselevel_sets) {
      for (tlag in tlag_sets) {
        for (noisetype in noisetype_sets) {
          load(paste(c('estimates','bd',bd,'_lag',tlag,'_timelength',timelength,'_noiselevel',noiselevel,'_noisetype',noisetype,'.RData'),sep = '',collapse = ''))
          for (m in 1:20) {
            for (n in 1:20) {
              datalow <- rbind(datalow,cbind('JensenRenyi2','rlow',bd,tlag,timelength,noiselevel,noisetype,round(estimates$da[[m]][[n]][1],2),round(estimates$true_synchrony_mu[[m]][n],2),round(estimates$crosscorrelation[[m]][n],2)))
            }
          }
        }
      }
    }
  }
}
dataavg <- NULL

for (timelength in timelength_sets) {
  for (bd in bd_sets) {
    for (noiselevel in noiselevel_sets) {
      for (tlag in tlag_sets) {
        for (noisetype in noisetype_sets) {
          load(paste(c('estimates','bd',bd,'_lag',tlag,'_timelength',timelength,'_noiselevel',noiselevel,'_noisetype',noisetype,'.RData'),sep = '',collapse = ''))
          for (m in 1:20) {
            for (n in 21:40) {
              dataavg <- rbind(dataavg,cbind('JensenRenyi2','ravg',bd,tlag,timelength,noiselevel,noisetype,round(estimates$da[[m]][[n]][1],2),round(estimates$true_synchrony_mu[[m]][n],2),round(estimates$crosscorrelation[[m]][n],2)))
            }
          }
        }
      }
    }
  }
}

datahigh <- NULL
for (timelength in timelength_sets) {
  for (bd in bd_sets) {
    for (noiselevel in noiselevel_sets) {
      for (tlag in tlag_sets) {
        for (noisetype in noisetype_sets) {
          load(paste(c('estimates','bd',bd,'_lag',tlag,'_timelength',timelength,'_noiselevel',noiselevel,'_noisetype',noisetype,'.RData'),sep = '',collapse = ''))
          for (m in 1:20) {
            for (n in 41:60) {
              datahigh <- rbind(datahigh,cbind('JensenRenyi2','rhigh',bd,tlag,timelength,noiselevel,noisetype,round(estimates$da[[m]][[n]][1],2),round(estimates$true_synchrony_mu[[m]][n],2),round(estimates$crosscorrelation[[m]][n],2)))
            }
          }
        }
      }
    }
  }
}

JensenRenyi2 <- rbind(datalow,dataavg,datahigh)
#######################################################
datalow <- NULL


for (timelength in timelength_sets) {
  for (bd in bd_sets) {
    for (noiselevel in noiselevel_sets) {
      for (tlag in tlag_sets) {
        for (noisetype in noisetype_sets) {
          load(paste(c('estimates','bd',bd,'_lag',tlag,'_timelength',timelength,'_noiselevel',noiselevel,'_noisetype',noisetype,'.RData'),sep = '',collapse = ''))
          for (m in 1:20) {
            for (n in 1:20) {
              datalow <- rbind(datalow,cbind('JensenRenyi6','rlow',bd,tlag,timelength,noiselevel,noisetype,round(estimates$da[[m]][[n]][5],2),round(estimates$true_synchrony_mu[[m]][n],2),round(estimates$crosscorrelation[[m]][n],2)))
            }
          }
        }
      }
    }
  }
}
dataavg <- NULL

for (timelength in timelength_sets) {
  for (bd in bd_sets) {
    for (noiselevel in noiselevel_sets) {
      for (tlag in tlag_sets) {
        for (noisetype in noisetype_sets) {
          load(paste(c('estimates','bd',bd,'_lag',tlag,'_timelength',timelength,'_noiselevel',noiselevel,'_noisetype',noisetype,'.RData'),sep = '',collapse = ''))
          for (m in 1:20) {
            for (n in 21:40) {
              dataavg <- rbind(dataavg,cbind('JensenRenyi6','ravg',bd,tlag,timelength,noiselevel,noisetype,round(estimates$da[[m]][[n]][5],2),round(estimates$true_synchrony_mu[[m]][n],2),round(estimates$crosscorrelation[[m]][n],2)))
            }
          }
        }
      }
    }
  }
}

datahigh <- NULL
for (timelength in timelength_sets) {
  for (bd in bd_sets) {
    for (noiselevel in noiselevel_sets) {
      for (tlag in tlag_sets) {
        for (noisetype in noisetype_sets) {
          load(paste(c('estimates','bd',bd,'_lag',tlag,'_timelength',timelength,'_noiselevel',noiselevel,'_noisetype',noisetype,'.RData'),sep = '',collapse = ''))
          for (m in 1:20) {
            for (n in 41:60) {
              datahigh <- rbind(datahigh,cbind('JensenRenyi6','rhigh',bd,tlag,timelength,noiselevel,noisetype,round(estimates$da[[m]][[n]][5],2),round(estimates$true_synchrony_mu[[m]][n],2),round(estimates$crosscorrelation[[m]][n],2)))
            }
          }
        }
      }
    }
  }
}

JensenRenyi6 <- rbind(datalow,dataavg,datahigh)
#################################################################
datalow <- NULL


for (timelength in timelength_sets) {
  for (bd in bd_sets) {
    for (noiselevel in noiselevel_sets) {
      for (tlag in tlag_sets) {
        for (noisetype in noisetype_sets) {
          load(paste(c('estimates','bd',bd,'_lag',tlag,'_timelength',timelength,'_noiselevel',noiselevel,'_noisetype',noisetype,'.RData'),sep = '',collapse = ''))
          for (m in 1:20) {
            for (n in 1:20) {
              datalow <- rbind(datalow,cbind('JensenRenyi10','rlow',bd,tlag,timelength,noiselevel,noisetype,round(estimates$da[[m]][[n]][9],2),round(estimates$true_synchrony_mu[[m]][n],2),round(estimates$crosscorrelation[[m]][n],2)))
            }
          }
        }
      }
    }
  }
}
dataavg <- NULL

for (timelength in timelength_sets) {
  for (bd in bd_sets) {
    for (noiselevel in noiselevel_sets) {
      for (tlag in tlag_sets) {
        for (noisetype in noisetype_sets) {
          load(paste(c('estimates','bd',bd,'_lag',tlag,'_timelength',timelength,'_noiselevel',noiselevel,'_noisetype',noisetype,'.RData'),sep = '',collapse = ''))
          for (m in 1:20) {
            for (n in 21:40) {
              dataavg <- rbind(dataavg,cbind('JensenRenyi10','ravg',bd,tlag,timelength,noiselevel,noisetype,round(estimates$da[[m]][[n]][9],2),round(estimates$true_synchrony_mu[[m]][n],2),round(estimates$crosscorrelation[[m]][n],2)))
            }
          }
        }
      }
    }
  }
}

datahigh <- NULL
for (timelength in timelength_sets) {
  for (bd in bd_sets) {
    for (noiselevel in noiselevel_sets) {
      for (tlag in tlag_sets) {
        for (noisetype in noisetype_sets) {
          load(paste(c('estimates','bd',bd,'_lag',tlag,'_timelength',timelength,'_noiselevel',noiselevel,'_noisetype',noisetype,'.RData'),sep = '',collapse = ''))
          for (m in 1:20) {
            for (n in 41:60) {
              datahigh <- rbind(datahigh,cbind('JensenRenyi10','rhigh',bd,tlag,timelength,noiselevel,noisetype,round(estimates$da[[m]][[n]][9],2),round(estimates$true_synchrony_mu[[m]][n],2),round(estimates$crosscorrelation[[m]][n],2)))
            }
          }
        }
      }
    }
  }
}

JensenRenyi10 <- rbind(datalow,dataavg,datahigh)
load('phasecoherence.RData')
dataestimates <- rbind(Kullback,phasesyn,Renyi01,Renyi05,Renyi09,
                       JensenShannon,JensenRenyi2,JensenRenyi6,JensenRenyi10,
                       mutual_information,phasecoherence)
dataestimates <- data.frame(dataestimates)
colnames(dataestimates) <- c("synchr_method","true_synchrony","bd",     
                        "tlag", "timelength", "noiselevel",
                        "noisetype",'estimates','true_synchrony_mu','true_synchrony_c')
dataestimates$true_synchrony <-  factor(dataestimates$true_synchrony)
dataestimates$bd <- factor(dataestimates$bd)
dataestimates$tlag <- factor(dataestimates$tlag)
dataestimates$timelength <- factor(dataestimates$timelength)
dataestimates$noiselevel <- factor(dataestimates$noiselevel)
dataestimates$noisetype <- factor(dataestimates$noisetype)
dataestimates$estimates <- as.numeric(dataestimates$estimates)
dataestimates$true_synchrony_mu <- as.numeric(dataestimates$true_synchrony_mu)
dataestimates$true_synchrony_c <- as.numeric(dataestimates$true_synchrony_c)
dataestimates$synchr_method <- factor(dataestimates$synchr_method)

save(dataestimates,file='estimates.RData')
View(dataestimates)
summary(dataestimates)
