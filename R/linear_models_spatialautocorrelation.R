library(colorspace)
library(spdep)
library(gridExtra)

####this divides a global table into realms
divide_table_realms<-function(table,predictor.variable){
  table[,predictor.variable]<-as.factor(table[,predictor.variable])
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
  variable.hotspot.realm<-list()
  for (i in 1:nrow(grid.cells.df)){
    table.realm<-table[table$cells%in%c(grid.cells.df[i,'start.cell']:grid.cells.df[i,'end.cell']),]
    variable.hotspot.realm[[i]]<-table.realm
  }
  names(variable.hotspot.realm)<-grid.cells.df$realms.names
  return(variable.hotspot.realm)
}



prepare_linearmodel_table<-function(vector.of.tables){
  tables<-lapply(vector.of.tables,function(x) read.table(x,header=T,sep='\t',stringsAsFactors = F))
  tables<-lapply(tables,function(x) {colnames(x)[1]<-'cells';return(x)})
  table<-Reduce(function(...) merge(..., by='cells'), tables)
  table<-na.omit(table)
  return(table)
}

####this runs a lm and and lm with moran.test (to check for spatial autocorrelation in the residuals)
run_lm_univariate_moran_output<-function(data,response,predictor,neighbours.list.w,name){
  cat('running lm','\n')
  lm.uni<-lm(get(response)~get(predictor),data=data)
  cat('running lm.morantest','\n')
  lm.uni.moran<-lm.morantest(lm.uni,listw=neighbours.list.w,zero.policy = TRUE)
  result<-cbind(response,predictor,round(summary(lm.uni)$coefficients[1,1],5),round(summary(lm.uni)$coefficients[2,1],5),round(summary(lm.uni)$coefficients[2,4],5),round(summary(lm.uni)$r.squared,5),round(summary(lm.uni)$adj.r.squared,5),round(AIC(lm.uni),5),round(unname(lm.uni.moran$estimate[1]),5),round(unname(lm.uni.moran$estimate[2]),5),round(unname(lm.uni.moran$p.value[1]),5))
  colnames(result)<-c('response','predictor','intercept','estimate','Pr(>|t|)','r.squared','adj.r.squared','AIC','obs.Moran.I','expected.Moran.I','p-value.Moran.I')
  #pdf(paste('./',response,'_',predictor,'_',name,'_LM_mammals.pdf',sep=''))
  #par(mfrow=c(2,2))
  #boxplot of variable
  #boxplot(get(response)~get(predictor),data=data,main=name,ylab=response,col=c('lightcoral','lightcyan'))
  #pie chart of hot vs not hot numbers
  #slices <- c(nrow(data[predictor==1,]),nrow(data[data$MamHotspot!='Hot',]))
  #lbls <- c(nrow(data[data$MamHotspot=='Hot',]), nrow(data[data$MamHotspot!='Hot',]))
  #pie(slices, labels = lbls, main=paste(name,'_hotspots',sep=''),col=c('lightcoral','lightcyan'))
  #tt<-ttheme_default(base_size = 6)
  #grid.table(as.data.frame(result,stringsAsFactors = F), theme=tt)
  #dev.off()
  return(result)
}

####this runs a sarlm 
run_errorsarlm_univariate<-function(data,response,predictor,neighbours.list.w,name){
  cat('running errorsarlm','\n')
  errorsarlm.uni<-errorsarlm(get(response)~get(predictor),data=data,listw = neighbours.list.w,zero.policy = TRUE,quiet=FALSE,method='spam')
  AIC.sarlm<-2*summary.sarlm(errorsarlm.uni,Nagelkerke = TRUE)$parameters-2*summary.sarlm(errorsarlm.uni,Nagelkerke = TRUE)$LL
  result<-cbind(response,predictor,round(summary.sarlm(errorsarlm.uni,Nagelkerke = TRUE)$Coef[1,1],5),round(summary.sarlm(errorsarlm.uni,Nagelkerke = TRUE)$Coef[2,1],5),round(summary.sarlm(errorsarlm.uni,Nagelkerke = TRUE)$Coef[2,4],5),round(summary(errorsarlm.uni,Nagelkerke = TRUE)$NK,5),round(AIC.sarlm,5))
  colnames(result)<-c('response','predictor','intercept','estimate','Pr(>|t|)','r.squared.Nagelkerke','AIC')
  #pdf(paste('./',response,'_',predictor,'_',name,'_SARLM_mammals.pdf',sep=''))
  #par(mfrow=c(2,2))
  #boxplot of variable
  #boxplot(get(response)~get(predictor),data=data,main=name,ylab=response,col=c('lightcoral','lightcyan'))
  #pie chart of hot vs not hot numbers
  #slices <- c(nrow(data[data$MamHotspot=='Hot',]),nrow(data[data$MamHotspot!='Hot',]))
  #lbls <- c(nrow(data[data$MamHotspot=='Hot',]), nrow(data[data$MamHotspot!='Hot',]))
  #pie(slices, labels = lbls, main=paste(name,'_hotspots',sep=''),col=c('lightcoral','lightcyan'))
  #tt<-ttheme_default(base_size = 6)
  #grid.table(as.data.frame(result,stringsAsFactors = F), theme=tt)
  #dev.off()
  return(result)
}

predict_errorsarlm_univariate<-function(data,response,predictor,neighbours.list.w,name){
  cat('running errorsarlm','\n')
  row.names(data)<-data$cells
  errorsarlm.uni<-errorsarlm(get(response)~get(predictor),data=data,listw = neighbours.list.w,zero.policy = TRUE,quiet=FALSE,method='spam')
  predict.sarlm<-predict.sarlm(errorsarlm.uni,list.w=neighbours.1000.w,zero.policy=TRUE,spChk=TRUE)
  predict.sarlm.df<-as.data.frame(cbind(data$cells,predict.sarlm),stringsAsFactors = F)
  colnames(predict.sarlm.df)<-c('cells','value')
  predict.sarlm.df$value<-as.numeric(predict.sarlm.df$value)
  table.sarlm.df<-merge(data,predict.sarlm.df)
  return(table.sarlm.df)
}

predict_errorsarlm_univariate2<-function(data,response,predictor,neighbours.list.w,name){
  cat('running errorsarlm','\n')
  row.names(data)<-data$cells
  errorsarlm.uni<-errorsarlm(get(response)~get(predictor),data=data,listw = neighbours.list.w,zero.policy = TRUE,quiet=FALSE,method='spam',tol.solve = 1.0e-12)
  predict.sarlm<-predict.sarlm(errorsarlm.uni,list.w=neighbours.1000.w,zero.policy=TRUE,spChk=TRUE)
  predict.sarlm.df<-as.data.frame(cbind(data$cells,predict.sarlm),stringsAsFactors = F)
  colnames(predict.sarlm.df)<-c('cells','value')
  predict.sarlm.df$value<-as.numeric(predict.sarlm.df$value)
  table.sarlm.df<-merge(data,predict.sarlm.df)
  return(table.sarlm.df)
}


#this is a wrapper to run a lm and a moran.test for spatial autocorrelation in a table
#table is the table with the grid data
#predictor.variable is 'hotspot'
#response.variables is a vector with responses to test c('Q1.residuals','Q4.residuals')
#mode is either 'global' or 'realms'
run_lm.morantest_table<-function(table,predictor.variable,response.variables,mode){
  table[,predictor.variable]<-as.factor(table[,predictor.variable])
  grids<-list.files(path='./output/grids/',pattern='grid_.*_100.rds')
  world.grid<-grep('World_RealmsMerged',grids)
  grid.world<-grids[grep('World_RealmsMerged',grids)]
  grid.realms.names<-grids[-grep('World_RealmsMerged',grids)]
  grid.world<-readRDS(paste('./output/grids/',grid.world,sep=''))
  if(mode=='global'){
    grid.coordinates<-lapply(grid.world[table$cells],function(x) coordinates(x))
    grid.coordinates<-do.call("rbind", grid.coordinates)
    neighbours.1000<-dnearneigh(grid.coordinates,d1=0,d2=1000)
    neighbours.1000.w<-nb2listw(neighbours.1000,style="W",zero.policy = TRUE)
    #response variables vs hotspot in global
    response.hotspot.global<-lapply(response.variables,function(x) run_lm_univariate_moran_output(data=table,response=x,predictor="hotspot",neighbours.list.w = neighbours.1000.w,name='global'))
    return(response.hotspot.global)
  }else if(mode=='realms'){
    ####response variables vs hotspot in realms
    #######divide dataset by realms
    #load realm grids
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
    for (i in 1:nrow(grid.cells.df)){
      table.realm<-table[table$cells%in%c(grid.cells.df[i,'start.cell']:grid.cells.df[i,'end.cell']),]
      grid.coordinates<-lapply(grid.world[table.realm$cells],function(x) coordinates(x))
      grid.coordinates<-do.call("rbind", grid.coordinates)
      neighbours.1000<-dnearneigh(grid.coordinates,d1=0,d2=1000)
      neighbours.1000.w<-nb2listw(neighbours.1000,style="W",zero.policy = TRUE)
      response.hotspot.realm[[i]]<-lapply(response.variables,function(x) run_lm_univariate_moran_output(data=table.realm,response=x,predictor='hotspot',neighbours.list.w = neighbours.1000.w,name=paste(grid.cells.df[i,'realms.names'],sep='')))
    }
    names(response.hotspot.realm)<-grid.cells.df[,'realms.names']
    return(response.hotspot.realm)
  }else{return('mode has to be one of global/realms')}
}

#this is a wrapper to run a lm and a moran.test for spatial autocorrelation in a table
#table is the table with the grid data
#predictor.variable is 'hotspot'
#response.variables is a vector with responses to test c('Q1.residuals','Q4.residuals')
#mode is either 'global' or 'realms'
run_sarlm_table<-function(table,predictor.variable,response.variables,mode){
  table[,predictor.variable]<-as.factor(table[,predictor.variable])
  grids<-list.files(path='./output/grids/',pattern='grid_.*_100.rds')
  world.grid<-grep('World_RealmsMerged',grids)
  grid.world<-grids[grep('World_RealmsMerged',grids)]
  grid.realms.names<-grids[-grep('World_RealmsMerged',grids)]
  grid.world<-readRDS(paste('./output/grids/',grid.world,sep=''))
  if(mode=='global'){
    grid.coordinates<-lapply(grid.world[table$cells],function(x) coordinates(x))
    grid.coordinates<-do.call("rbind", grid.coordinates)
    neighbours.1000<-dnearneigh(grid.coordinates,d1=0,d2=1000)
    neighbours.1000.w<-nb2listw(neighbours.1000,style="W",zero.policy = TRUE)
    #response variables vs hotspot in global
    response.hotspot.global<-lapply(response.variables,function(x) run_errorsarlm_univariate(data=table,response=x,predictor="hotspot",neighbours.list.w = neighbours.1000.w,name='global'))
    return(response.hotspot.global)
  }else if(mode=='realms'){
    ####response variables vs hotspot in realms
    #######divide dataset by realms
    #load realm grids
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
    for (i in 1:nrow(grid.cells.df)){
      table.realm<-table[table$cells%in%c(grid.cells.df[i,'start.cell']:grid.cells.df[i,'end.cell']),]
      grid.coordinates<-lapply(grid.world[table.realm$cells],function(x) coordinates(x))
      grid.coordinates<-do.call("rbind", grid.coordinates)
      neighbours.1000<-dnearneigh(grid.coordinates,d1=0,d2=1000)
      neighbours.1000.w<-nb2listw(neighbours.1000,style="W",zero.policy = TRUE)
      response.hotspot.realm[[i]]<-lapply(response.variables,function(x) run_errorsarlm_univariate(data=table.realm,response=x,predictor='hotspot',neighbours.list.w = neighbours.1000.w,name=paste(grid.cells.df[i,'realms.names'],sep='')))
    }
    names(response.hotspot.realm)<-grid.cells.df[,'realms.names']
    return(response.hotspot.realm)
  }else{return('mode has to be one of global/realms')}
}

#this is a wrapper to run a sarlm and correct the values accounting for spatial autocorrelation
#table is the table with the grid data
#predictor.variable is 'hotspot'
#response.variables is a vector with responses to test c('Q1.residuals','Q4.residuals')
#mode is either 'global' or 'realms'
predict_sarlm_table<-function(table,predictor.variable,response.variables,mode){
  table[,predictor.variable]<-as.factor(table[,predictor.variable])
  grids<-list.files(path='./output/grids/',pattern='grid_.*_100.rds')
  world.grid<-grep('World_RealmsMerged',grids)
  grid.world<-grids[grep('World_RealmsMerged',grids)]
  grid.realms.names<-grids[-grep('World_RealmsMerged',grids)]
  grid.world<-readRDS(paste('./output/grids/',grid.world,sep=''))
  if(mode=='global'){
    grid.coordinates<-lapply(grid.world[table$cells],function(x) coordinates(x))
    grid.coordinates<-do.call("rbind", grid.coordinates)
    neighbours.1000<-dnearneigh(grid.coordinates,d1=0,d2=1000)
    neighbours.1000.w<-nb2listw(neighbours.1000,style="W",zero.policy = TRUE)
    #response variables vs hotspot in global
    response.hotspot.global<-lapply(response.variables,function(x) predict_errorsarlm_univariate(data=table,response=x,predictor="hotspot",neighbours.list.w = neighbours.1000.w,name='global'))
    names(response.hotspot.global)<-response.variables
    return(response.hotspot.global)
  }else if(mode=='realms'){
    ####response variables vs hotspot in realms
    #######divide dataset by realms
    #load realm grids
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
    for (i in 1:nrow(grid.cells.df)){
      table.realm<-table[table$cells%in%c(grid.cells.df[i,'start.cell']:grid.cells.df[i,'end.cell']),]
      grid.coordinates<-lapply(grid.world[table.realm$cells],function(x) coordinates(x))
      grid.coordinates<-do.call("rbind", grid.coordinates)
      neighbours.1000<-dnearneigh(grid.coordinates,d1=0,d2=1000)
      neighbours.1000.w<-nb2listw(neighbours.1000,style="W",zero.policy = TRUE)
      response.hotspot.realm[[i]]<-lapply(response.variables,function(x) predict_errorsarlm_univariate(data=table.realm,response=x,predictor='hotspot',neighbours.list.w = neighbours.1000.w,name=paste(grid.cells.df[i,'realms.names'],sep='')))
      names(response.hotspot.realm[[i]])<-response.variables
    }
    names(response.hotspot.realm)<-grid.cells.df[,'realms.names']
    return(response.hotspot.realm)
  }else{return('mode has to be one of global/realms')}
}

predict_sarlm_table2<-function(table,predictor.variable,response.variables,mode){
  table[,predictor.variable]<-as.factor(table[,predictor.variable])
  grids<-list.files(path='./output/grids/',pattern='grid_.*_100.rds')
  world.grid<-grep('World_RealmsMerged',grids)
  grid.world<-grids[grep('World_RealmsMerged',grids)]
  grid.realms.names<-grids[-grep('World_RealmsMerged',grids)]
  grid.world<-readRDS(paste('./output/grids/',grid.world,sep=''))
  if(mode=='global'){
    grid.coordinates<-lapply(grid.world[table$cells],function(x) coordinates(x))
    grid.coordinates<-do.call("rbind", grid.coordinates)
    neighbours.1000<-dnearneigh(grid.coordinates,d1=0,d2=1000)
    neighbours.1000.w<-nb2listw(neighbours.1000,style="W",zero.policy = TRUE)
    #response variables vs hotspot in global
    response.hotspot.global<-lapply(response.variables,function(x) predict_errorsarlm_univariate2(data=table,response=x,predictor="hotspot",neighbours.list.w = neighbours.1000.w,name='global'))
    names(response.hotspot.global)<-response.variables
    return(response.hotspot.global)
  }else if(mode=='realms'){
    ####response variables vs hotspot in realms
    #######divide dataset by realms
    #load realm grids
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
    for (i in 1:nrow(grid.cells.df)){
      table.realm<-table[table$cells%in%c(grid.cells.df[i,'start.cell']:grid.cells.df[i,'end.cell']),]
      grid.coordinates<-lapply(grid.world[table.realm$cells],function(x) coordinates(x))
      grid.coordinates<-do.call("rbind", grid.coordinates)
      neighbours.1000<-dnearneigh(grid.coordinates,d1=0,d2=1000)
      neighbours.1000.w<-nb2listw(neighbours.1000,style="W",zero.policy = TRUE)
      response.hotspot.realm[[i]]<-lapply(response.variables,function(x) predict_errorsarlm_univariate2(data=table.realm,response=x,predictor='hotspot',neighbours.list.w = neighbours.1000.w,name=paste(grid.cells.df[i,'realms.names'],sep='')))
      names(response.hotspot.realm[[i]])<-response.variables
    }
    names(response.hotspot.realm)<-grid.cells.df[,'realms.names']
    return(response.hotspot.realm)
  }else{return('mode has to be one of global/realms')}
}
