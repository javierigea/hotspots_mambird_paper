library(RColorBrewer)
library(colorspace)
library(gplots)
environmentalmodels_lm_variables_heatmap<-function(table.bioclim,table.hab.elevation,table.CCV,table.tectonic,table.distance,variable.vector){
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
  table.distance<-table.distance[,which(colnames(table.distance)%in%c('cells',variable.vector))]
  #merge all tables
  table.environment<-merge(table.bioclim,table.hab.elevation)
  table.environment<-merge(table.environment,table.CCV)
  table.environment<-merge(table.environment,table.tectonic)
  table.environment<-merge(table.environment,table.distance)
  table.environment.scaled<-as.data.frame(cbind(table.environment$cells,table.environment$hotspot,scale(table.environment[,-c(which(colnames(table.environment)=='cells'),which(colnames(table.environment)=='hotspot'))])),stringsAsFactors=F)
  colnames(table.environment.scaled)[c(1,2)]<-c('cells','hotspot')
  #analyse global
  lm.environmental.global<-lapply(colnames(table.environment.scaled)[-c(1,2)],function(x) lm(get(x)~hotspot,data=table.environment.scaled))
  lm.environmental.global.coef<-lapply(lm.environmental.global,function(x) c(round(summary(x)$coefficients[2,1],5),round(summary(x)$coefficients[2,4],5)))
  names(lm.environmental.global.coef)<-colnames(table.environment.scaled)[-c(1,2)]
  #analyse by realms
  grids<-list.files(path='./output/grids/',pattern='grid_.*_100.rds')
  world.grid<-grep('World_RealmsMerged',grids)
  grid.world<-grids[grep('World_RealmsMerged',grids)]
  grid.realms.names<-grids[-grep('World_RealmsMerged',grids)]
  grid.world<-readRDS(paste('./output/grids/',grid.world,sep=''))
  grid.realms<-lapply(grid.realms.names,function(x) readRDS(paste('./output/grids/',x,sep='')))
  grid.realms.names<-sub(grid.realms.names,pattern='grid_',replacement='')
  grid.realms.names<-sub(grid.realms.names,pattern='_100.rds',replacement='')
  grid.realms.ncells<-unlist(lapply(grid.realms,function(x)length(x)))
  grid.realms.start.cell<-numeric(length(grid.realms.ncells))
  grid.realms.end.cell<-numeric(length(grid.realms.ncells))
  for(i in 1:length(grid.realms.ncells)){
    if(i==1){
      grid.realms.start.cell[i]<-1
      grid.realms.end.cell[i]<-grid.realms.ncells[i]
      next
    }else{
      grid.realms.start.cell[i]<-grid.realms.end.cell[i-1]+1
      grid.realms.end.cell[i]<-grid.realms.start.cell[i]+grid.realms.ncells[i]-1
      next
    }
  }
  grid.cells.df<-as.data.frame(cbind(grid.realms.names,grid.realms.ncells,grid.realms.start.cell,grid.realms.end.cell),stringsAsFactors = F)
  colnames(grid.cells.df)<-c('realms.names','ncells','start.cell','end.cell')
  grid.cells.df$ncells<-as.numeric(grid.cells.df$ncells)
  grid.cells.df$start.cell<-as.numeric(grid.cells.df$start.cell)
  grid.cells.df$end.cell<-as.numeric(grid.cells.df$end.cell)
  #remove Oceanic grid.cells.df (just 79 cells and most of them hotspot)
  grid.cells.df<-grid.cells.df[-grep('Oceanic',grid.cells.df$realms.names),]
  response.hotspot.realm<-list()
  response.hotspot.coef.realm<-list()
  for (i in 1:nrow(grid.cells.df)){
    table.environment.scaled.realm<-table.environment.scaled[table.environment.scaled$cells%in%c(grid.cells.df[i,'start.cell']:grid.cells.df[i,'end.cell']),]
    response.hotspot.realm[[i]]<-lapply(colnames(table.environment.scaled.realm)[-c(1,2)],function(x) lm(get(x)~hotspot,data=table.environment.scaled.realm))
    response.hotspot.coef.realm[[i]]<-lapply(response.hotspot.realm[[i]],function(x) c(round(summary(x)$coefficients[2,1],5),round(summary(x)$coefficients[2,4],5)))
    names(response.hotspot.coef.realm[[i]])<-colnames(table.environment.scaled.realm)[-c(1,2)]
    
  }
  names(response.hotspot.realm)<-grid.cells.df$realms.names
  names(response.hotspot.coef.realm)<-grid.cells.df$realms.names
  #build matrix for heatmap
  #rows are variable, cols are global + realms
  heatmap.matrix<-matrix(c(unlist(lapply(lm.environmental.global.coef,function(x)x[1])),unlist(lapply(response.hotspot.coef.realm,function(x)lapply(x,function(x)x[1])))),ncol=length(response.hotspot.coef.realm)+1,nrow=length(response.hotspot.coef.realm[[1]]),byrow = F)
  #significance matrix
  pvalue.matrix<-matrix(c(unlist(lapply(lm.environmental.global.coef,function(x)x[2])),unlist(lapply(response.hotspot.coef.realm,function(x)lapply(x,function(x)x[2])))),ncol=length(response.hotspot.coef.realm)+1,nrow=length(response.hotspot.coef.realm[[1]]),byrow = F)
  colnames(heatmap.matrix)<-c('global',names(response.hotspot.coef.realm))
  rownames(heatmap.matrix)<-names(response.hotspot.coef.realm[[1]])
  colnames(pvalue.matrix)<-c('global',names(response.hotspot.coef.realm))
  rownames(pvalue.matrix)<-names(response.hotspot.coef.realm[[1]])
  pvalue.adjustmatrix<-matrix(p.adjust(pvalue.matrix),byrow = FALSE,ncol=7)
  colnames(pvalue.adjustmatrix)<-c('global',names(response.hotspot.coef.realm))
  rownames(pvalue.adjustmatrix)<-names(response.hotspot.coef.realm[[1]])
  
  #this plots the heatmap
  heatmap.ordered.matrix<-heatmap.matrix
  heatmap.ordered.matrix[which(pvalue.adjustmatrix>0.05)]<-NA
  heatmap.ordered.matrix<-heatmap.ordered.matrix[c('mean.AR','mean.MAT','mean.AR.CCV','mean.MAT.CCV','tectonic.movement','mean.TRI','n.habitats','average.distance.to.class.1000neighbours'),]
  #heatmap.ordered.matrix<-heatmap.ordered.matrix[c(1,2,3,4,8,9,10,11,12,5,6,7),]
  #this heatmap has to be modified: the paleyellow overlapping 0 will be grey, remove the borders,etc
  heatmap.2(heatmap.ordered.matrix, Rowv=NA, Colv=NA, dendrogram='none',col = rev(brewer.pal(11,"RdYlBu")), scale="none", trace='none',margins=c(5,10),keysize = 1,denscol = NA,key.title=NA,key.xlab='Z-score',key.ylab=NA,sepwidth=c(0.02,0.02),sepcolor="black",colsep=1:ncol(heatmap.ordered.matrix),rowsep=1:nrow(heatmap.ordered.matrix))
  
  #this plots a heatmap where the black boxes indicate p.value<0.05
  heatmap.2(pvalue.matrix, Rowv=NA, Colv=NA, dendrogram='none',col = c('black',rep('white',19)), scale="none", trace='none',margins=c(5,10),keysize = 1,denscol = NA,key.title=NA,key.xlab='Z-score',key.ylab=NA,sepwidth=c(0.01,0.01),sepcolor="grey",colsep=1:ncol(heatmap.matrix),rowsep=1:nrow(heatmap.matrix))
  
  
  
}


environmentalmodelsNPP_lm_variables_heatmap<-function(table.NPP,table.hab.elevation,table.CCV,table.tectonic,table.distance,variable.vector){
  #calculate COVs
  table.NPP$cov.NPP<-sqrt(table.NPP$var.NPP)/table.NPP$mean.NPP
  table.hab.elevation$cov.TRI<-sqrt(table.hab.elevation$var.TRI)/table.hab.elevation$mean.TRI
  table.CCV$cov.AR.CCV<-sqrt(table.CCV$varAR.CCV)/table.CCV$mean.AR.CCV
  table.CCV$cov.MAT.CCV<-sqrt(table.CCV$varMAT.CCV)/table.CCV$mean.MAT.CCV
  #drop unwanted columns
  table.NPP<-table.NPP[,which(colnames(table.NPP)%in%c('cells',variable.vector,'hotspot'))]
  table.hab.elevation<-table.hab.elevation[,which(colnames(table.hab.elevation)%in%c('cells',variable.vector))]
  table.CCV<-table.CCV[,which(colnames(table.CCV)%in%c('cells',variable.vector))]
  table.tectonic<-table.tectonic[,which(colnames(table.tectonic)%in%c('cells',variable.vector))]
  table.distance<-table.distance[,which(colnames(table.distance)%in%c('cells',variable.vector))]
  #merge all tables
  table.environment<-merge(table.NPP,table.hab.elevation)
  table.environment<-merge(table.environment,table.CCV)
  table.environment<-merge(table.environment,table.tectonic)
  table.environment<-merge(table.environment,table.distance)
  table.environment.scaled<-as.data.frame(cbind(table.environment$cells,table.environment$hotspot,scale(table.environment[,-c(which(colnames(table.environment)=='cells'),which(colnames(table.environment)=='hotspot'))])),stringsAsFactors=F)
  colnames(table.environment.scaled)[c(1,2)]<-c('cells','hotspot')
  #analyse global
  lm.environmental.global<-lapply(colnames(table.environment.scaled)[-c(1,2)],function(x) lm(get(x)~hotspot,data=table.environment.scaled))
  lm.environmental.global.coef<-lapply(lm.environmental.global,function(x) c(round(summary(x)$coefficients[2,1],5),round(summary(x)$coefficients[2,4],5)))
  names(lm.environmental.global.coef)<-colnames(table.environment.scaled)[-c(1,2)]
  #analyse by realms
  grids<-list.files(path='./output/grids/',pattern='grid_.*_100.rds')
  world.grid<-grep('World_RealmsMerged',grids)
  grid.world<-grids[grep('World_RealmsMerged',grids)]
  grid.realms.names<-grids[-grep('World_RealmsMerged',grids)]
  grid.world<-readRDS(paste('./output/grids/',grid.world,sep=''))
  grid.realms<-lapply(grid.realms.names,function(x) readRDS(paste('./output/grids/',x,sep='')))
  grid.realms.names<-sub(grid.realms.names,pattern='grid_',replacement='')
  grid.realms.names<-sub(grid.realms.names,pattern='_100.rds',replacement='')
  grid.realms.ncells<-unlist(lapply(grid.realms,function(x)length(x)))
  grid.realms.start.cell<-numeric(length(grid.realms.ncells))
  grid.realms.end.cell<-numeric(length(grid.realms.ncells))
  for(i in 1:length(grid.realms.ncells)){
    if(i==1){
      grid.realms.start.cell[i]<-1
      grid.realms.end.cell[i]<-grid.realms.ncells[i]
      next
    }else{
      grid.realms.start.cell[i]<-grid.realms.end.cell[i-1]+1
      grid.realms.end.cell[i]<-grid.realms.start.cell[i]+grid.realms.ncells[i]-1
      next
    }
  }
  grid.cells.df<-as.data.frame(cbind(grid.realms.names,grid.realms.ncells,grid.realms.start.cell,grid.realms.end.cell),stringsAsFactors = F)
  colnames(grid.cells.df)<-c('realms.names','ncells','start.cell','end.cell')
  grid.cells.df$ncells<-as.numeric(grid.cells.df$ncells)
  grid.cells.df$start.cell<-as.numeric(grid.cells.df$start.cell)
  grid.cells.df$end.cell<-as.numeric(grid.cells.df$end.cell)
  #remove Oceanic grid.cells.df (just 79 cells and most of them hotspot)
  grid.cells.df<-grid.cells.df[-grep('Oceanic',grid.cells.df$realms.names),]
  response.hotspot.realm<-list()
  response.hotspot.coef.realm<-list()
  for (i in 1:nrow(grid.cells.df)){
    table.environment.scaled.realm<-table.environment.scaled[table.environment.scaled$cells%in%c(grid.cells.df[i,'start.cell']:grid.cells.df[i,'end.cell']),]
    response.hotspot.realm[[i]]<-lapply(colnames(table.environment.scaled.realm)[-c(1,2)],function(x) lm(get(x)~hotspot,data=table.environment.scaled.realm))
    response.hotspot.coef.realm[[i]]<-lapply(response.hotspot.realm[[i]],function(x) c(round(summary(x)$coefficients[2,1],5),round(summary(x)$coefficients[2,4],5)))
    names(response.hotspot.coef.realm[[i]])<-colnames(table.environment.scaled.realm)[-c(1,2)]
    
  }
  names(response.hotspot.realm)<-grid.cells.df$realms.names
  names(response.hotspot.coef.realm)<-grid.cells.df$realms.names
  #build matrix for heatmap
  #rows are variable, cols are global + realms
  heatmap.matrix<-matrix(c(unlist(lapply(lm.environmental.global.coef,function(x)x[1])),unlist(lapply(response.hotspot.coef.realm,function(x)lapply(x,function(x)x[1])))),ncol=length(response.hotspot.coef.realm)+1,nrow=length(response.hotspot.coef.realm[[1]]),byrow = F)
  #significance matrix
  pvalue.matrix<-matrix(c(unlist(lapply(lm.environmental.global.coef,function(x)x[2])),unlist(lapply(response.hotspot.coef.realm,function(x)lapply(x,function(x)x[2])))),ncol=length(response.hotspot.coef.realm)+1,nrow=length(response.hotspot.coef.realm[[1]]),byrow = F)
  colnames(heatmap.matrix)<-c('global',names(response.hotspot.coef.realm))
  rownames(heatmap.matrix)<-names(response.hotspot.coef.realm[[1]])
  colnames(pvalue.matrix)<-c('global',names(response.hotspot.coef.realm))
  rownames(pvalue.matrix)<-names(response.hotspot.coef.realm[[1]])
  pvalue.adjustmatrix<-matrix(p.adjust(pvalue.matrix),byrow = FALSE,ncol=7)
  colnames(pvalue.adjustmatrix)<-c('global',names(response.hotspot.coef.realm))
  rownames(pvalue.adjustmatrix)<-names(response.hotspot.coef.realm[[1]])
  
  #this plots the heatmap
  heatmap.ordered.matrix<-heatmap.matrix
  heatmap.ordered.matrix[which(pvalue.adjustmatrix>0.05)]<-NA
  heatmap.ordered.matrix<-heatmap.ordered.matrix[c('mean.NPP','mean.AR.CCV','mean.MAT.CCV','tectonic.movement','mean.TRI','n.habitats','average.distance.to.class.1000neighbours'),]
  #heatmap.ordered.matrix<-heatmap.ordered.matrix[c(1,2,3,4,8,9,10,11,12,5,6,7),]
  #this heatmap has to be modified: the paleyellow overlapping 0 will be grey, remove the borders,etc
  heatmap.2(heatmap.ordered.matrix, Rowv=NA, Colv=NA, dendrogram='none',col = rev(brewer.pal(11,"RdYlBu")), scale="none", trace='none',margins=c(5,10),keysize = 1,denscol = NA,key.title=NA,key.xlab='Z-score',key.ylab=NA,sepwidth=c(0.02,0.02),sepcolor="black",colsep=1:ncol(heatmap.ordered.matrix),rowsep=1:nrow(heatmap.ordered.matrix))
  
  #this plots a heatmap where the black boxes indicate p.value<0.05
  heatmap.2(pvalue.matrix, Rowv=NA, Colv=NA, dendrogram='none',col = c('black',rep('white',19)), scale="none", trace='none',margins=c(5,10),keysize = 1,denscol = NA,key.title=NA,key.xlab='Z-score',key.ylab=NA,sepwidth=c(0.01,0.01),sepcolor="grey",colsep=1:ncol(heatmap.matrix),rowsep=1:nrow(heatmap.matrix))
  
  
  
}
environmentalmodels_sarlm_variables_heatmap<-function(table.bioclim,table.hab.elevation,table.CCV,table.tectonic,variable.vector){
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
  table.distance<-table.distance[,which(colnames(table.distance)%in%c('cells',variable.vector))]
  #merge all tables
  table.environment<-merge(table.bioclim,table.hab.elevation)
  table.environment<-merge(table.environment,table.CCV)
  table.environment<-merge(table.environment,table.tectonic)
  table.environment<-merge(table.environment,table.distance)
  table.environment.scaled<-as.data.frame(cbind(table.environment$cells,table.environment$hotspot,scale(table.environment[,-c(which(colnames(table.environment)=='cells'),which(colnames(table.environment)=='hotspot'))])),stringsAsFactors=F)
  colnames(table.environment.scaled)[c(1,2)]<-c('cells','hotspot')
  #analyse global
  grids<-list.files(path='./output/grids/',pattern='grid_.*_100.rds')
  world.grid<-grep('World_RealmsMerged',grids)
  grid.world<-grids[grep('World_RealmsMerged',grids)]
  grid.realms.names<-grids[-grep('World_RealmsMerged',grids)]
  grid.world<-readRDS(paste('./output/grids/',grid.world,sep=''))
  grid.coordinates<-lapply(grid.world[table.environment.scaled$cells],function(x) coordinates(x))
  grid.coordinates<-do.call("rbind", grid.coordinates)
  neighbours.1000<-dnearneigh(grid.coordinates,d1=0,d2=1000)
  neighbours.1000.w<-nb2listw(neighbours.1000,style="W",zero.policy = TRUE)
  #response variables vs hotspot in global
  response.hotspot.global.sarlm<-lapply(colnames(table.environment.scaled)[-c(1,2)],function(x) errorsarlm(get(x)~hotspot,data=table.environment.scaled,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'))
  response.hotspot.global.sarlm.coef<-lapply(response.hotspot.global.sarlm,function(x) c(round(summary.sarlm(x)$Coef[2,1],3),round(summary.sarlm(x)$Coef[2,4],3)))
  names(response.hotspot.global.sarlm.coef)<-colnames(table.environment.scaled)[-c(1,2)]
  #analyse by realms
  grids<-list.files(path='./output/grids/',pattern='grid_.*_100.rds')
  world.grid<-grep('World_RealmsMerged',grids)
  grid.world<-grids[grep('World_RealmsMerged',grids)]
  grid.realms.names<-grids[-grep('World_RealmsMerged',grids)]
  grid.world<-readRDS(paste('./output/grids/',grid.world,sep=''))
  grid.realms<-lapply(grid.realms.names,function(x) readRDS(paste('./output/grids/',x,sep='')))
  grid.realms.names<-sub(grid.realms.names,pattern='grid_',replacement='')
  grid.realms.names<-sub(grid.realms.names,pattern='_100.rds',replacement='')
  grid.realms.ncells<-unlist(lapply(grid.realms,function(x)length(x)))
  grid.realms.start.cell<-numeric(length(grid.realms.ncells))
  grid.realms.end.cell<-numeric(length(grid.realms.ncells))
  for(i in 1:length(grid.realms.ncells)){
    if(i==1){
      grid.realms.start.cell[i]<-1
      grid.realms.end.cell[i]<-grid.realms.ncells[i]
      next
    }else{
      grid.realms.start.cell[i]<-grid.realms.end.cell[i-1]+1
      grid.realms.end.cell[i]<-grid.realms.start.cell[i]+grid.realms.ncells[i]-1
      next
    }
  }
  grid.cells.df<-as.data.frame(cbind(grid.realms.names,grid.realms.ncells,grid.realms.start.cell,grid.realms.end.cell),stringsAsFactors = F)
  colnames(grid.cells.df)<-c('realms.names','ncells','start.cell','end.cell')
  grid.cells.df$ncells<-as.numeric(grid.cells.df$ncells)
  grid.cells.df$start.cell<-as.numeric(grid.cells.df$start.cell)
  grid.cells.df$end.cell<-as.numeric(grid.cells.df$end.cell)
  #remove Oceanic grid.cells.df (just 79 cells and most of them hotspot)
  grid.cells.df<-grid.cells.df[-grep('Oceanic',grid.cells.df$realms.names),]
  response.hotspot.realm.sarlm<-list()
  response.hotspot.realm.sarlm.coef<-list()
  for (i in 1:nrow(grid.cells.df)){
    table.environment.scaled.realm<-table.environment.scaled[table.environment.scaled$cells%in%c(grid.cells.df[i,'start.cell']:grid.cells.df[i,'end.cell']),]
    grid.coordinates<-lapply(grid.world[table.environment.scaled.realm$cells],function(x) coordinates(x))
    grid.coordinates<-do.call("rbind", grid.coordinates)
    neighbours.1000<-dnearneigh(grid.coordinates,d1=0,d2=1000)
    neighbours.1000.w<-nb2listw(neighbours.1000,style="W",zero.policy = TRUE)
    response.hotspot.realm.sarlm[[i]]<-lapply(colnames(table.environment.scaled.realm)[-c(1,2)],function(x) errorsarlm(get(x)~hotspot,data=table.environment.scaled.realm,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'))
    response.hotspot.realm.sarlm.coef[[i]]<-lapply(response.hotspot.realm.sarlm[[i]],function(x) c(round(summary.sarlm(x)$Coef[2,1],3),round(summary.sarlm(x)$Coef[2,4],3)))
    names(response.hotspot.realm.sarlm.coef[[i]])<-colnames(table.environment.scaled.realm)[-c(1,2)]
    
  }
  names(response.hotspot.realm.sarlm)<-grid.cells.df$realms.names
  names(response.hotspot.realm.sarlm.coef)<-grid.cells.df$realms.names
  #build matrix for heatmap
  #rows are variable, cols are global + realms
  heatmap.matrix.sarlm<-matrix(c(unlist(lapply(response.hotspot.global.sarlm.coef,function(x)x[1])),unlist(lapply(response.hotspot.realm.sarlm.coef,function(x)lapply(x,function(x)x[1])))),ncol=length(response.hotspot.realm.sarlm.coef)+1,nrow=length(response.hotspot.realm.sarlm.coef[[1]]),byrow = F)
  #significance matrix
  pvalue.matrix.sarlm<-matrix(c(unlist(lapply(response.hotspot.global.sarlm.coef,function(x)x[2])),unlist(lapply(response.hotspot.realm.sarlm.coef,function(x)lapply(x,function(x)x[2])))),ncol=length(response.hotspot.realm.sarlm.coef)+1,nrow=length(response.hotspot.realm.sarlm.coef[[1]]),byrow = F)
  colnames(heatmap.matrix.sarlm)<-c('global',names(response.hotspot.realm.sarlm.coef))
  rownames(heatmap.matrix.sarlm)<-names(response.hotspot.realm.sarlm.coef[[1]])
  colnames(pvalue.matrix.sarlm)<-c('global',names(response.hotspot.realm.sarlm.coef))
  rownames(pvalue.matrix.sarlm)<-names(response.hotspot.realm.sarlm.coef[[1]])
  pvalue.adjustmatrix.sarlm<-matrix(p.adjust(pvalue.matrix.sarlm),byrow = FALSE,ncol=7)
  colnames(pvalue.adjustmatrix.sarlm)<-c('global',names(response.hotspot.realm.sarlm.coef))
  rownames(pvalue.adjustmatrix.sarlm)<-names(response.hotspot.realm.sarlm.coef[[1]])
  
  #this plots the heatmap
  heatmap.ordered.matrix.sarlm<-heatmap.matrix.sarlm
  heatmap.ordered.matrix.sarlm[which(pvalue.adjustmatrix.sarlm>0.05)]<-NA
  heatmap.ordered.matrix.sarlm<-heatmap.ordered.matrix.sarlm[c('mean.AR','mean.MAT','mean.AR.CCV','mean.MAT.CCV','tectonic.movement','mean.TRI','n.habitats','average.distance.to.class.1000neighbours'),]
  #heatmap.ordered.matrix.sarlm<-heatmap.ordered.matrix.sarlm[c(1,2,3,4,8,9,10,11,12,5,6,7),]
  #this heatmap has to be modified: the paleyellow overlapping 0 will be grey, remove the borders,etc
  heatmap.2(heatmap.ordered.matrix.sarlm, Rowv=NA, Colv=NA, dendrogram='none',col = rev(brewer.pal(11,"RdYlBu")), scale="none", trace='none',margins=c(5,10),keysize = 1,denscol = NA,key.title=NA,key.xlab='Z-score',key.ylab=NA,sepwidth=c(0.02,0.02),sepcolor="black",colsep=1:ncol(heatmap.ordered.matrix.sarlm),rowsep=1:nrow(heatmap.ordered.matrix.sarlm))
  #this plots a heatmap where the black boxes indicate p.value<0.05
  heatmap.2(pvalue.matrix.sarlm, Rowv=NA, Colv=NA, dendrogram='none',col = c('black',rep('white',19)), scale="none", trace='none',margins=c(5,10),keysize = 1,denscol = NA,key.title=NA,key.xlab='Z-score',key.ylab=NA,sepwidth=c(0.01,0.01),sepcolor="grey",colsep=1:ncol(heatmap.matrix),rowsep=1:nrow(heatmap.matrix))
  
  
  
  
}

environmentalmodelsNPP_sarlm_variables_heatmap<-function(table.NPP,table.hab.elevation,table.CCV,table.tectonic,variable.vector){
  #calculate COVs
  table.NPP$cov.NPP<-sqrt(table.NPP$var.NPP)/table.NPP$mean.NPP
  table.hab.elevation$cov.TRI<-sqrt(table.hab.elevation$var.TRI)/table.hab.elevation$mean.TRI
  table.CCV$cov.AR.CCV<-sqrt(table.CCV$varAR.CCV)/table.CCV$mean.AR.CCV
  table.CCV$cov.MAT.CCV<-sqrt(table.CCV$varMAT.CCV)/table.CCV$mean.MAT.CCV
  #drop unwanted columns
  table.NPP<-table.NPP[,which(colnames(table.NPP)%in%c('cells',variable.vector,'hotspot'))]
  table.hab.elevation<-table.hab.elevation[,which(colnames(table.hab.elevation)%in%c('cells',variable.vector))]
  table.CCV<-table.CCV[,which(colnames(table.CCV)%in%c('cells',variable.vector))]
  table.tectonic<-table.tectonic[,which(colnames(table.tectonic)%in%c('cells',variable.vector))]
  #merge all tables
  table.environment<-merge(table.NPP,table.hab.elevation)
  table.environment<-merge(table.environment,table.CCV)
  table.environment<-merge(table.environment,table.tectonic)
  table.environment.scaled<-as.data.frame(cbind(table.environment$cells,table.environment$hotspot,scale(table.environment[,-c(which(colnames(table.environment)=='cells'),which(colnames(table.environment)=='hotspot'))])),stringsAsFactors=F)
  colnames(table.environment.scaled)[c(1,2)]<-c('cells','hotspot')
  #analyse global
  grids<-list.files(path='./output/grids/',pattern='grid_.*_100.rds')
  world.grid<-grep('World_RealmsMerged',grids)
  grid.world<-grids[grep('World_RealmsMerged',grids)]
  grid.realms.names<-grids[-grep('World_RealmsMerged',grids)]
  grid.world<-readRDS(paste('./output/grids/',grid.world,sep=''))
  grid.coordinates<-lapply(grid.world[table.environment.scaled$cells],function(x) coordinates(x))
  grid.coordinates<-do.call("rbind", grid.coordinates)
  neighbours.1000<-dnearneigh(grid.coordinates,d1=0,d2=1000)
  neighbours.1000.w<-nb2listw(neighbours.1000,style="W",zero.policy = TRUE)
  #response variables vs hotspot in global
  response.hotspot.global.sarlm<-lapply(colnames(table.environment.scaled)[-c(1,2)],function(x) errorsarlm(get(x)~hotspot,data=table.environment.scaled,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'))
  response.hotspot.global.sarlm.coef<-lapply(response.hotspot.global.sarlm,function(x) c(round(summary.sarlm(x)$Coef[2,1],3),round(summary.sarlm(x)$Coef[2,4],3)))
  names(response.hotspot.global.sarlm.coef)<-colnames(table.environment.scaled)[-c(1,2)]
  #analyse by realms
  grids<-list.files(path='./output/grids/',pattern='grid_.*_100.rds')
  world.grid<-grep('World_RealmsMerged',grids)
  grid.world<-grids[grep('World_RealmsMerged',grids)]
  grid.realms.names<-grids[-grep('World_RealmsMerged',grids)]
  grid.world<-readRDS(paste('./output/grids/',grid.world,sep=''))
  grid.realms<-lapply(grid.realms.names,function(x) readRDS(paste('./output/grids/',x,sep='')))
  grid.realms.names<-sub(grid.realms.names,pattern='grid_',replacement='')
  grid.realms.names<-sub(grid.realms.names,pattern='_100.rds',replacement='')
  grid.realms.ncells<-unlist(lapply(grid.realms,function(x)length(x)))
  grid.realms.start.cell<-numeric(length(grid.realms.ncells))
  grid.realms.end.cell<-numeric(length(grid.realms.ncells))
  for(i in 1:length(grid.realms.ncells)){
    if(i==1){
      grid.realms.start.cell[i]<-1
      grid.realms.end.cell[i]<-grid.realms.ncells[i]
      next
    }else{
      grid.realms.start.cell[i]<-grid.realms.end.cell[i-1]+1
      grid.realms.end.cell[i]<-grid.realms.start.cell[i]+grid.realms.ncells[i]-1
      next
    }
  }
  grid.cells.df<-as.data.frame(cbind(grid.realms.names,grid.realms.ncells,grid.realms.start.cell,grid.realms.end.cell),stringsAsFactors = F)
  colnames(grid.cells.df)<-c('realms.names','ncells','start.cell','end.cell')
  grid.cells.df$ncells<-as.numeric(grid.cells.df$ncells)
  grid.cells.df$start.cell<-as.numeric(grid.cells.df$start.cell)
  grid.cells.df$end.cell<-as.numeric(grid.cells.df$end.cell)
  #remove Oceanic grid.cells.df (just 79 cells and most of them hotspot)
  grid.cells.df<-grid.cells.df[-grep('Oceanic',grid.cells.df$realms.names),]
  response.hotspot.realm.sarlm<-list()
  response.hotspot.realm.sarlm.coef<-list()
  for (i in 1:nrow(grid.cells.df)){
    table.environment.scaled.realm<-table.environment.scaled[table.environment.scaled$cells%in%c(grid.cells.df[i,'start.cell']:grid.cells.df[i,'end.cell']),]
    grid.coordinates<-lapply(grid.world[table.environment.scaled.realm$cells],function(x) coordinates(x))
    grid.coordinates<-do.call("rbind", grid.coordinates)
    neighbours.1000<-dnearneigh(grid.coordinates,d1=0,d2=1000)
    neighbours.1000.w<-nb2listw(neighbours.1000,style="W",zero.policy = TRUE)
    response.hotspot.realm.sarlm[[i]]<-lapply(colnames(table.environment.scaled.realm)[-c(1,2)],function(x) errorsarlm(get(x)~hotspot,data=table.environment.scaled.realm,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'))
    response.hotspot.realm.sarlm.coef[[i]]<-lapply(response.hotspot.realm.sarlm[[i]],function(x) c(round(summary.sarlm(x)$Coef[2,1],3),round(summary.sarlm(x)$Coef[2,4],3)))
    names(response.hotspot.realm.sarlm.coef[[i]])<-colnames(table.environment.scaled.realm)[-c(1,2)]
    
  }
  names(response.hotspot.realm.sarlm)<-grid.cells.df$realms.names
  names(response.hotspot.realm.sarlm.coef)<-grid.cells.df$realms.names
  #build matrix for heatmap
  #rows are variable, cols are global + realms
  #heatmap.matrix.sarlm<-matrix(c(unlist(lapply(response.hotspot.global.sarlm.coef,function(x)x[1])),unlist(lapply(response.hotspot.realm.sarlm.coef,function(x)lapply(x,function(x)x[1])))),ncol=length(response.hotspot.realm.sarlm.coef)+1,nrow=length(response.hotspot.realm.sarlm.coef[[1]]),byrow = F)
  heatmap.matrix.sarlm<-matrix(c(unlist(lapply(response.hotspot.realm.sarlm.coef,function(x)lapply(x,function(x)x[1])))),ncol=length(response.hotspot.realm.sarlm.coef),nrow=length(response.hotspot.realm.sarlm.coef[[1]]),byrow = F)
  #significance matrix
  #pvalue.matrix.sarlm<-matrix(c(unlist(lapply(response.hotspot.global.sarlm.coef,function(x)x[2])),unlist(lapply(response.hotspot.realm.sarlm.coef,function(x)lapply(x,function(x)x[2])))),ncol=length(response.hotspot.realm.sarlm.coef)+1,nrow=length(response.hotspot.realm.sarlm.coef[[1]]),byrow = F)
  pvalue.matrix.sarlm<-matrix(c(unlist(lapply(response.hotspot.realm.sarlm.coef,function(x)lapply(x,function(x)x[2])))),ncol=length(response.hotspot.realm.sarlm.coef),nrow=length(response.hotspot.realm.sarlm.coef[[1]]),byrow = F)
  #colnames(heatmap.matrix.sarlm)<-c('global',names(response.hotspot.realm.sarlm.coef))
  #colnames(heatmap.matrix.sarlm)<-c(names(response.hotspot.realm.sarlm.coef))
  colnames(heatmap.matrix.sarlm)<-c(names(response.hotspot.realm.sarlm.coef))
  rownames(heatmap.matrix.sarlm)<-names(response.hotspot.realm.sarlm.coef[[1]])
  #colnames(pvalue.matrix.sarlm)<-c('global',names(response.hotspot.realm.sarlm.coef))
  colnames(pvalue.matrix.sarlm)<-c(names(response.hotspot.realm.sarlm.coef))
  rownames(pvalue.matrix.sarlm)<-names(response.hotspot.realm.sarlm.coef[[1]])
  pvalue.adjustmatrix.sarlm<-matrix(p.adjust(pvalue.matrix.sarlm),byrow = FALSE,ncol=6)
  #colnames(pvalue.adjustmatrix.sarlm)<-c('global',names(response.hotspot.realm.sarlm.coef))
  colnames(pvalue.adjustmatrix.sarlm)<-c(names(response.hotspot.realm.sarlm.coef))
  rownames(pvalue.adjustmatrix.sarlm)<-names(response.hotspot.realm.sarlm.coef[[1]])
  
  #this plots the heatmap
  heatmap.ordered.matrix.sarlm<-heatmap.matrix.sarlm
  heatmap.ordered.matrix.sarlm[which(pvalue.adjustmatrix.sarlm>0.05)]<-NA
  heatmap.ordered.matrix.sarlm<-heatmap.ordered.matrix.sarlm[c('mean.NPP','mean.AR.CCV','mean.MAT.CCV','tectonic.movement','mean.TRI','n.habitats'),]
  #heatmap.ordered.matrix.sarlm<-heatmap.ordered.matrix.sarlm[c(1,2,3,4,8,9,10,11,12,5,6,7),]
  #this heatmap has to be modified: the paleyellow overlapping 0 will be grey, remove the borders,etc
  heatmap.2(heatmap.ordered.matrix.sarlm, Rowv=NA, Colv=NA, dendrogram='none',col = rev(brewer.pal(11,"RdYlBu")), scale="none", trace='none',margins=c(5,10),keysize = 1,denscol = NA,key.title=NA,key.xlab='Z-score',key.ylab=NA,sepwidth=c(0.02,0.02),sepcolor="black",colsep=1:ncol(heatmap.ordered.matrix.sarlm),rowsep=1:nrow(heatmap.ordered.matrix.sarlm))
  #this plots a heatmap where the black boxes indicate p.value<0.05
  heatmap.2(pvalue.matrix.sarlm, Rowv=NA, Colv=NA, dendrogram='none',col = c('black',rep('white',19)), scale="none", trace='none',margins=c(5,10),keysize = 1,denscol = NA,key.title=NA,key.xlab='Z-score',key.ylab=NA,sepwidth=c(0.01,0.01),sepcolor="grey",colsep=1:ncol(heatmap.matrix.sarlm),rowsep=1:nrow(heatmap.matrix.sarlm))
  
  
  
  
}




#get table of environment

#scale all variables

#run lms of all variables by hotspot, get variables coefficient + pvalue

#split by realms

#run lms of all variables by hotspot, get variables coefficient + pvalue

#plot all in heat map, columns global + 6 realms, rows variables
