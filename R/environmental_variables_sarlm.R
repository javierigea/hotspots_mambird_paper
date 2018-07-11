library(Hmisc) 
library(corrplot)
library(MASS)

source('./R/plot_boxplot_predictedtable_map.R')

corr_variables<-function(table.bioclim,table.hab.elevation,table.CCV,table.tectonic){
  #calculate COVs
  table.bioclim$cov.AR<-sqrt(table.bioclim$var.AR)/table.bioclim$mean.AR
  table.bioclim$cov.ATR<-sqrt(table.bioclim$var.ATR)/table.bioclim$mean.ATR
  table.bioclim$cov.MAT<-sqrt(table.bioclim$var.MAT)/table.bioclim$mean.MAT
  table.hab.elevation$cov.TRI<-sqrt(table.hab.elevation$var.TRI)/table.hab.elevation$mean.TRI
  table.CCV$cov.AR.CCV<-sqrt(table.CCV$varAR.CCV)/table.CCV$mean.AR.CCV
  table.CCV$cov.MAT.CCV<-sqrt(table.CCV$varMAT.CCV)/table.CCV$mean.MAT.CCV
  #drop unwanted columns
  table.bioclim<-table.bioclim[,-which(colnames(table.bioclim)%in%c('number.of.species','number.of.species.wend'))]
  table.hab.elevation<-table.hab.elevation[,-which(colnames(table.hab.elevation)%in%c('number.of.species','number.of.species.wend','hotspot'))]
  table.CCV<-table.CCV[,-which(colnames(table.CCV)%in%c('number.of.species','number.of.species.wend','hotspot'))]
  table.tectonic<-table.tectonic[,which(colnames(table.tectonic)%in%c('cells','tectonic.movement'))]
  #merge all tables
  table.environment<-merge(table.bioclim,table.hab.elevation)
  table.environment<-merge(table.environment,table.CCV)
  table.environment<-merge(table.environment,table.tectonic)
  #drop cells column
  table.environment<-table.environment[-which(colnames(table.environment)%in%c('cells'))]
  #center and scale variables
  scaled.table.environment<-as.data.frame(scale(table.environment),stringsAsFactors=F)
  #check correlations of predictor variables
  corr<-rcorr(as.matrix(scaled.table.environment)) # compute Pearson's (or spearman's corr) with rcorr from Hmisc package. I like rcorr as it allows to separately access the correlations, the # or observations and the p-value. ?rcorr is worth a read.
  corr_r<-as.matrix(corr[[1]])# Access the correlation matrix. 
  corr_r[,1]# subset the correlation of "a" (=var1 ) with the rest if you want.
  pval<-as.matrix(corr[[3]])# get the p-values
  corrplot(corr_r,method="number",type="lower",diag=FALSE,tl.col="black",tl.cex=.5,tl.offset=0.1,tl.srt=45,number.cex=.7)# plot all pairs
  corrplot(corr_r,p.mat = pval,sig.level=0.05,insig = "blank",method="circle",type="lower",diag=FALSE,tl.col="black",tl.cex=1,tl.offset=0.1,tl.srt=45)# plot pairs with significance cutoff defined by "p.mat"
}


corr_variablesNPP<-function(table.NPP,table.hab.elevation,table.CCV,table.tectonic){
  #calculate COVs
  table.NPP$cov.NPP<-sqrt(table.NPP$var.NPP)/table.NPP$mean.NPP
  table.hab.elevation$cov.TRI<-sqrt(table.hab.elevation$var.TRI)/table.hab.elevation$mean.TRI
  table.CCV$cov.AR.CCV<-sqrt(table.CCV$varAR.CCV)/table.CCV$mean.AR.CCV
  table.CCV$cov.MAT.CCV<-sqrt(table.CCV$varMAT.CCV)/table.CCV$mean.MAT.CCV
  #drop unwanted columns
  table.NPP<-table.NPP[,-which(colnames(table.NPP)%in%c('number.of.species','number.of.species.wend'))]
  table.hab.elevation<-table.hab.elevation[,-which(colnames(table.hab.elevation)%in%c('number.of.species','number.of.species.wend','hotspot'))]
  table.CCV<-table.CCV[,-which(colnames(table.CCV)%in%c('number.of.species','number.of.species.wend','hotspot'))]
  table.tectonic<-table.tectonic[,which(colnames(table.tectonic)%in%c('cells','tectonic.movement'))]
  #merge all tables
  table.environment<-merge(table.NPP,table.hab.elevation)
  table.environment<-merge(table.environment,table.CCV)
  table.environment<-merge(table.environment,table.tectonic)
  #drop cells column
  table.environment<-table.environment[-which(colnames(table.environment)%in%c('cells'))]
  #center and scale variables
  scaled.table.environment<-as.data.frame(scale(table.environment),stringsAsFactors=F)
  #check correlations of predictor variables
  corr<-rcorr(as.matrix(scaled.table.environment)) # compute Pearson's (or spearman's corr) with rcorr from Hmisc package. I like rcorr as it allows to separately access the correlations, the # or observations and the p-value. ?rcorr is worth a read.
  corr_r<-as.matrix(corr[[1]])# Access the correlation matrix. 
  corr_r[,1]# subset the correlation of "a" (=var1 ) with the rest if you want.
  pval<-as.matrix(corr[[3]])# get the p-values
  corrplot(corr_r,method="number",type="lower",diag=FALSE,tl.col="black",tl.cex=.5,tl.offset=0.1,tl.srt=45,number.cex=.7)# plot all pairs
  corrplot(corr_r,p.mat = pval,sig.level=0.05,insig = "blank",method="circle",type="lower",diag=FALSE,tl.col="black",tl.cex=1,tl.offset=0.1,tl.srt=45)# plot pairs with significance cutoff defined by "p.mat"
}

corr_variablesNPP_sarlmlat<-function(table.NPP,table.hab.elevation,table.CCV,table.tectonic){
  #calculate COVs
  table.NPP$cov.NPP<-sqrt(table.NPP$var.NPP)/table.NPP$mean.NPP
  table.hab.elevation$cov.TRI<-sqrt(table.hab.elevation$var.TRI)/table.hab.elevation$mean.TRI
  table.CCV$cov.AR.CCV<-sqrt(table.CCV$varAR.CCV)/table.CCV$mean.AR.CCV
  table.CCV$cov.MAT.CCV<-sqrt(table.CCV$varMAT.CCV)/table.CCV$mean.MAT.CCV
  #drop unwanted columns
  #table.NPP<-table.NPP[,-which(colnames(table.NPP)%in%c('number.of.species','number.of.species.wend'))]
  #table.hab.elevation<-table.hab.elevation[,-which(colnames(table.hab.elevation)%in%c('number.of.species','number.of.species.wend','hotspot'))]
  #table.CCV<-table.CCV[,-which(colnames(table.CCV)%in%c('number.of.species','number.of.species.wend','hotspot'))]
  table.tectonic<-table.tectonic[,which(colnames(table.tectonic)%in%c('cells','tectonic.movement'))]
  #merge all tables
  table.environment<-merge(table.NPP,table.hab.elevation)
  table.environment<-merge(table.environment,table.CCV)
  table.environment<-merge(table.environment,table.tectonic)
  #drop cells column
  table.environment<-table.environment[-which(colnames(table.environment)%in%c('cells'))]
  #center and scale variables
  scaled.table.environment<-as.data.frame(scale(table.environment),stringsAsFactors=F)
  #check correlations of predictor variables
  corr<-rcorr(as.matrix(scaled.table.environment)) # compute Pearson's (or spearman's corr) with rcorr from Hmisc package. I like rcorr as it allows to separately access the correlations, the # or observations and the p-value. ?rcorr is worth a read.
  corr_r<-as.matrix(corr[[1]])# Access the correlation matrix. 
  corr_r[,1]# subset the correlation of "a" (=var1 ) with the rest if you want.
  pval<-as.matrix(corr[[3]])# get the p-values
  corrplot(corr_r,method="number",type="lower",diag=FALSE,tl.col="black",tl.cex=.5,tl.offset=0.1,tl.srt=45,number.cex=.7)# plot all pairs
  corrplot(corr_r,p.mat = pval,sig.level=0.05,insig = "blank",method="circle",type="lower",diag=FALSE,tl.col="black",tl.cex=1,tl.offset=0.1,tl.srt=45)# plot pairs with significance cutoff defined by "p.mat"
}
environmentalmodels_sarlm_variables<-function(table.bioclim,table.hab.elevation,table.CCV,table.tectonic,variable.vector){
  #calculate COVs
  table.bioclim$cov.AR<-sqrt(table.bioclim$var.AR)/table.bioclim$mean.AR
  table.bioclim$cov.ATR<-sqrt(table.bioclim$var.ATR)/table.bioclim$mean.ATR
  table.bioclim$cov.MAT<-sqrt(table.bioclim$var.MAT)/table.bioclim$mean.MAT
  table.hab.elevation$cov.TRI<-sqrt(table.hab.elevation$var.TRI)/table.hab.elevation$mean.TRI
  table.CCV$cov.AR.CCV<-sqrt(table.CCV$varAR.CCV)/table.CCV$mean.AR.CCV
  table.CCV$cov.MAT.CCV<-sqrt(table.CCV$varMAT.CCV)/table.CCV$mean.MAT.CCV
  #drop unwanted columns
  table.bioclim<-table.bioclim[,which(colnames(table.bioclim)%in%c('cells',variable.vector,'hotspot'))]
  table.hab.elevation<-table.hab.elevation[,which(colnames(table.hab.elevation)%in%c('cells',variable.vector))]
  table.CCV<-table.CCV[,which(colnames(table.CCV)%in%c('cells',variable.vector))]
  table.tectonic<-table.tectonic[,which(colnames(table.tectonic)%in%c('cells',variable.vector))]
  #merge all tables
  table.environment<-merge(table.bioclim,table.hab.elevation)
  table.environment<-merge(table.environment,table.CCV)
  table.environment<-merge(table.environment,table.tectonic)
  #table.environment.scaled<-as.data.frame(cbind(table.environment$cells,table.environment$hotspot,scale(table.environment[,-c(which(colnames(table.environment)=='cells'),which(colnames(table.environment)=='hotspot'))])),stringsAsFactors=F)
  #colnames(table.environment.scaled)[c(1,2)]<-c('cells','hotspot')
  environmental.predicted.scaled<-predict_sarlm_table(table=table.environment.scaled,predictor.variable='hotspot',response.variables = variable.vector,mode = 'realms')
  environmental.predicted<-predict_sarlm_table(table=table.environment,predictor.variable='hotspot',response.variables = variable.vector,mode = 'realms')
  lapply(variable.vector,function(x)plot_propdifferencesmeans_plus_vioplots(table.object =environmental.predicted,variable = x))
}

environmentalmodels_scaled_sarlm_variables<-function(table.bioclim,table.hab.elevation,table.CCV,table.tectonic,variable.vector){
  #calculate COVs
  table.bioclim$cov.AR<-sqrt(table.bioclim$var.AR)/table.bioclim$mean.AR
  table.bioclim$cov.ATR<-sqrt(table.bioclim$var.ATR)/table.bioclim$mean.ATR
  table.bioclim$cov.MAT<-sqrt(table.bioclim$var.MAT)/table.bioclim$mean.MAT
  table.hab.elevation$cov.TRI<-sqrt(table.hab.elevation$var.TRI)/table.hab.elevation$mean.TRI
  table.CCV$cov.AR.CCV<-sqrt(table.CCV$varAR.CCV)/table.CCV$mean.AR.CCV
  table.CCV$cov.MAT.CCV<-sqrt(table.CCV$varMAT.CCV)/table.CCV$mean.MAT.CCV
  #drop unwanted columns
  table.bioclim<-table.bioclim[,which(colnames(table.bioclim)%in%c('cells',variable.vector,'hotspot'))]
  table.hab.elevation<-table.hab.elevation[,which(colnames(table.hab.elevation)%in%c('cells',variable.vector))]
  table.CCV<-table.CCV[,which(colnames(table.CCV)%in%c('cells',variable.vector))]
  table.tectonic<-table.tectonic[,which(colnames(table.tectonic)%in%c('cells',variable.vector))]
  #merge all tables
  table.environment<-merge(table.bioclim,table.hab.elevation)
  table.environment<-merge(table.environment,table.CCV)
  table.environment<-merge(table.environment,table.tectonic)
  table.environment.scaled<-as.data.frame(cbind(table.environment$cells,table.environment$hotspot,scale(table.environment[,-c(which(colnames(table.environment)=='cells'),which(colnames(table.environment)=='hotspot'))])),stringsAsFactors=F)
  colnames(table.environment.scaled)[c(1,2)]<-c('cells','hotspot')
  environmental.predicted.scaled<-predict_sarlm_table(table=table.environment.scaled,predictor.variable='hotspot',response.variables = variable.vector,mode = 'realms')
  environmental.predicted<-predict_sarlm_table(table=table.environment,predictor.variable='hotspot',response.variables = variable.vector,mode = 'realms')
  lapply(variable.vector,function(x)plot_propdifferencesmeans_plus_vioplots(table.object =environmental.predicted,variable = x))
}

environmentalmodels_lm_variables<-function(table.bioclim,table.hab.elevation,table.CCV,table.tectonic,variable.vector){
  #calculate COVs
  table.bioclim$cov.AR<-sqrt(table.bioclim$var.AR)/table.bioclim$mean.AR
  table.bioclim$cov.ATR<-sqrt(table.bioclim$var.ATR)/table.bioclim$mean.ATR
  table.bioclim$cov.MAT<-sqrt(table.bioclim$var.MAT)/table.bioclim$mean.MAT
  table.hab.elevation$cov.TRI<-sqrt(table.hab.elevation$var.TRI)/table.hab.elevation$mean.TRI
  table.CCV$cov.AR.CCV<-sqrt(table.CCV$varAR.CCV)/table.CCV$mean.AR.CCV
  table.CCV$cov.MAT.CCV<-sqrt(table.CCV$varMAT.CCV)/table.CCV$mean.MAT.CCV
  #drop unwanted columns
  table.bioclim<-table.bioclim[,which(colnames(table.bioclim)%in%c('cells',variable.vector,'hotspot'))]
  table.hab.elevation<-table.hab.elevation[,which(colnames(table.hab.elevation)%in%c('cells',variable.vector))]
  table.CCV<-table.CCV[,which(colnames(table.CCV)%in%c('cells',variable.vector))]
  table.tectonic<-table.tectonic[,which(colnames(table.tectonic)%in%c('cells',variable.vector))]
  #merge all tables
  table.environment<-merge(table.bioclim,table.hab.elevation)
  table.environment<-merge(table.environment,table.CCV)
  table.environment<-merge(table.environment,table.tectonic)
  #table.environment.scaled<-as.data.frame(cbind(table.environment$cells,table.environment$hotspot,scale(table.environment[,-c(which(colnames(table.environment)=='cells'),which(colnames(table.environment)=='hotspot'))])),stringsAsFactors=F)
  #colnames(table.environment.scaled)[c(1,2)]<-c('cells','hotspot')
  environmental.predicted.scaled<-predict_sarlm_table(table=table.environment.scaled,predictor.variable='hotspot',response.variables = variable.vector,mode = 'realms')
  environmental.predicted<-predict_sarlm_table(table=table.environment,predictor.variable='hotspot',response.variables = variable.vector,mode = 'realms')
  lapply(variable.vector,function(x)plot_propdifferencesmeans_plus_vioplots(table.object =environmental.predicted,variable = x))
}










