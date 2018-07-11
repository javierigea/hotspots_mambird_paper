library(vioplot)
library(spdep)
get_averagedistanceneighboursclass<-function(hotspots.table.file,outputpath,name){
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
  #get the hotspots in
  hotspots.table<-read.table(hotspots.table.file,header=T,sep='\t',stringsAsFactors = F)
  hotspots.table<-hotspots.table[,c('cells','hotspot')]
  #for each realm, get coordinates of centroid
  median.distance.class.realm<-numeric()
  summary.hotspots.neighbour.class.distances<-list()
  summary.nonhotspots.neighbour.class.distances<-list()
  pdf(paste(outputpath,'/median_distance_to_class_neighbours_',name,'.pdf',sep=''),width=10,height=7)
  hotspot.table.realm.all<-vector()
  par(mfrow=c(2,3))
  hotspot.distances.realms<-list()
  non.hotspot.distances.realms<-list()
  for (i in 1:nrow(grid.cells.df)){
    cat(i,'\n')
    #get coordinates of centroid of each cell
    grid.coordinates<-lapply(grid.world,function(x) coordinates(x))
    grid.coordinates<-do.call("rbind", grid.coordinates)
    #get table for realm
    hotspot.table.realm<-hotspots.table[(hotspots.table$cells%in%c(grid.cells.df$start.cell[i]:grid.cells.df$end.cell[i])),]
    #get distances to neigbouring cells in a 1000 km radius
    neighbours.1000<-dnearneigh(grid.coordinates[grid.cells.df$start.cell[i]:grid.cells.df$end.cell[i],],d1=0,d2=1000)
    #get a matrix of distances of each cell in realm with neighbours
    dist.realm.matrix<-as.matrix(dist(grid.coordinates[hotspot.table.realm$cells,],upper=TRUE))
    #for each cell, check if hotspot or non hotspot, get neighbours that are the same class and average distance
    hotspot.vector<-hotspot.table.realm$cells
    names(hotspot.vector)<-hotspot.table.realm$hotspot
    median.distance.class<-numeric()
    for(a in 1:length(hotspot.vector)){
      if(names(hotspot.vector[a])=='1'){
        neighbours.1000.subseted<-neighbours.1000[[a]][neighbours.1000[[a]]%in%which(names(hotspot.vector)==1)]
        median.distance.class[a]<-mean(dist.realm.matrix[a,neighbours.1000.subseted])
      }else if(names(hotspot.vector[a])=='0'){
        neighbours.1000.subseted<-neighbours.1000[[a]][neighbours.1000[[a]]%in%which(names(hotspot.vector)==0)]
        median.distance.class[a]<-mean(dist.realm.matrix[a,neighbours.1000.subseted])
      }
    }
    
    #another way, storing all distances to plot them later
    cell.distance.to.neighbours<-list()
    for(a in 1:length(hotspot.vector)){
      if(names(hotspot.vector[a])=='1'){
        neighbours.1000.subseted<-neighbours.1000[[a]][neighbours.1000[[a]]%in%which(names(hotspot.vector)==1)]
        cell.distance.to.neighbours[[a]]<-median(dist.realm.matrix[a,neighbours.1000.subseted])
      }else if(names(hotspot.vector[a])=='0'){
        neighbours.1000.subseted<-neighbours.1000[[a]][neighbours.1000[[a]]%in%which(names(hotspot.vector)==0)]
        cell.distance.to.neighbours[[a]]<-median(dist.realm.matrix[a,neighbours.1000.subseted])
      }
    }
    
    
    hotspot.table.realm$median.distance.to.class.neighbours<-unlist(cell.distance.to.neighbours)
    names(median.distance.class)<-hotspot.table.realm$cells
    cat(nrow(hotspot.table.realm),'\n')
    hotspot.table.realm.all<-rbind(hotspot.table.realm.all,hotspot.table.realm)
    #dist.realm.hotspot.matrix<-as.matrix(dist(grid.coordinates[hotspot.table.realm[hotspot.table.realm$hotspot==1,'cells'],],upper=TRUE))
    #median.dist.realm.hotspot<-apply(dist.realm.hotspot.matrix,1,function(x) median(x))
    #mean.dist.realm.hotspot<-apply(dist.realm.hotspot.matrix,1,function(x) mean(x))
    #dist.realm.nonhotspot.matrix<-as.matrix(dist(grid.coordinates[hotspot.table.realm[hotspot.table.realm$hotspot==0,'cells'],],upper=TRUE))
    #median.dist.realm.nonhotspot<-apply(dist.realm.nonhotspot.matrix,1,function(x) median(x))
    #mean.dist.realm.nonhotspot<-apply(dist.realm.nonhotspot.matrix,1,function(x) mean(x))
    #hotspot.table.realm$dist.to.class<-NA
    #hotspot.table.realm[hotspot.table.realm$hotspot==1,'dist.to.class']<-median.dist.realm.hotspot
    #hotspot.table.realm[hotspot.table.realm$hotspot==0,'dist.to.class']<-median.dist.realm.nonhotspot
    
    #boxplot(median.distance.to.class.neighbours~hotspot,data=hotspot.table.realm,main=grid.cells.df[i,'realms.names'])
    plot(c(1,1),x=c(0,2),ylim=c(0,1000),type='n',xaxt='n',ylab='median.distance.to.class.neighbours',xlab='',main=grid.cells.df[i,'realms.names'])
    vioplot(hotspot.table.realm[hotspot.table.realm$hotspot==0,'median.distance.to.class.neighbours'][!is.na(hotspot.table.realm[hotspot.table.realm$hotspot==0,'median.distance.to.class.neighbours'])],at=0.5,col='blue',add=T)
    vioplot(hotspot.table.realm[hotspot.table.realm$hotspot==1,'median.distance.to.class.neighbours'][!is.na(hotspot.table.realm[hotspot.table.realm$hotspot==1,'median.distance.to.class.neighbours'])],at=1.5,col='red',add=T)
    summary.hotspots.neighbour.class.distances[[i]]<-summary(hotspot.table.realm[hotspot.table.realm$hotspot==1,'median.distance.to.class.neighbours'][!is.na(hotspot.table.realm[hotspot.table.realm$hotspot==1,'median.distance.to.class.neighbours'])])
    summary.nonhotspots.neighbour.class.distances[[i]]<-summary(hotspot.table.realm[hotspot.table.realm$hotspot==0,'median.distance.to.class.neighbours'][!is.na(hotspot.table.realm[hotspot.table.realm$hotspot==0,'median.distance.to.class.neighbours'])])
    t.test<-t.test(hotspot.table.realm[hotspot.table.realm$hotspot==0,'median.distance.to.class.neighbours'][!is.na(hotspot.table.realm[hotspot.table.realm$hotspot==0,'median.distance.to.class.neighbours'])],hotspot.table.realm[hotspot.table.realm$hotspot==1,'median.distance.to.class.neighbours'][!is.na(hotspot.table.realm[hotspot.table.realm$hotspot==1,'median.distance.to.class.neighbours'])])
    cat(t.test$p.value,'\n')
    hotspot.distances.realms[[i]]<-hotspot.table.realm[hotspot.table.realm$hotspot==1,'median.distance.to.class.neighbours'][!is.na(hotspot.table.realm[hotspot.table.realm$hotspot==1,'median.distance.to.class.neighbours'])]
    non.hotspot.distances.realms[[i]]<-hotspot.table.realm[hotspot.table.realm$hotspot==0,'median.distance.to.class.neighbours'][!is.na(hotspot.table.realm[hotspot.table.realm$hotspot==0,'median.distance.to.class.neighbours'])]
    cat(mean(hotspot.table.realm[hotspot.table.realm$hotspot==0,'median.distance.to.class.neighbours'][!is.na(hotspot.table.realm[hotspot.table.realm$hotspot==0,'median.distance.to.class.neighbours'])]),'\n')
    cat(mean(hotspot.table.realm[hotspot.table.realm$hotspot==1,'median.distance.to.class.neighbours'][!is.na(hotspot.table.realm[hotspot.table.realm$hotspot==1,'median.distance.to.class.neighbours'])]),'\n')
    axis(1,at=c(0.5,1.5),labels=c('NH','H'))
  }
  dev.off()
  summary.hotspots.neighbour.class.distances<-as.data.frame(do.call('rbind',summary.hotspots.neighbour.class.distances))
  summary.hotspots.neighbour.class.distances$realm<-c('Afro','Austral','Indo','Nearctic','Neotrop','Palearctic')
  write.table(summary.hotspots.neighbour.class.distances,file=paste(outputpath,'/mediandistancetoclass_1000neighbours_hotspots_',name,'.txt',sep=''),sep='\t',quote=F,row.names=F)
    
  summary.nonhotspots.neighbour.class.distances<-as.data.frame(do.call('rbind',summary.nonhotspots.neighbour.class.distances))
  summary.nonhotspots.neighbour.class.distances$realm<-c('Afro','Austral','Indo','Nearctic','Neotrop','Palearctic')
  write.table(summary.nonhotspots.neighbour.class.distances,file=paste(outputpath,'/mediandistancetoclass_1000neighbours_nonhotspots_',name,'.txt',sep=''),sep='\t',quote=F,row.names=F)
  
}

################for mammals
##
##library(spdep)
###analyse by realms
##grids<-list.files(path='./output/grids/',pattern='grid_.*_100.rds')
##world.grid<-grep('World_RealmsMerged',grids)
##grid.world<-grids[grep('World_RealmsMerged',grids)]
##grid.realms.names<-grids[-grep('World_RealmsMerged',grids)]
##grid.world<-readRDS(paste('./output/grids/',grid.world,sep=''))
##grid.realms<-lapply(grid.realms.names,function(x) readRDS(paste('./output/grids/',x,sep='')))
##grid.realms.names<-sub(grid.realms.names,pattern='grid_',replacement='')
##grid.realms.names<-sub(grid.realms.names,pattern='_100.rds',replacement='')
##grid.realms.ncells<-unlist(lapply(grid.realms,function(x)length(x)))
##grid.realms.start.cell<-numeric(length(grid.realms.ncells))
##grid.realms.end.cell<-numeric(length(grid.realms.ncells))
##for(i in 1:length(grid.realms.ncells)){
##  if(i==1){
##    grid.realms.start.cell[i]<-1
##    grid.realms.end.cell[i]<-grid.realms.ncells[i]
##    next
##  }else{
##    grid.realms.start.cell[i]<-grid.realms.end.cell[i-1]+1
##    grid.realms.end.cell[i]<-grid.realms.start.cell[i]+grid.realms.ncells[i]-1
##    next
##  }
##}
##grid.cells.df<-as.data.frame(cbind(grid.realms.names,grid.realms.ncells,grid.realms.start.cell,grid.realms.end.cell),stringsAsFactors = F)
##colnames(grid.cells.df)<-c('realms.names','ncells','start.cell','end.cell')
##grid.cells.df$ncells<-as.numeric(grid.cells.df$ncells)
##grid.cells.df$start.cell<-as.numeric(grid.cells.df$start.cell)
##grid.cells.df$end.cell<-as.numeric(grid.cells.df$end.cell)
###remove Oceanic grid.cells.df (just 79 cells and most of them hotspot)
##grid.cells.df<-grid.cells.df[-grep('Oceanic',grid.cells.df$realms.names),]
###get the hotspots in
##hotspots.table<-read.table('./output/mammals/tables/100_all_realms_number.of.species.wend_0.8.txt',header=T,sep='\t',stringsAsFactors = F)
##hotspots.table<-hotspots.table[,c('cells','hotspot')]
###for each realm, get coordinates of centroid
##median.distance.class.realm<-numeric()
##summary.hotspots.neighbour.class.distances<-list()
##summary.nonhotspots.neighbour.class.distances<-list()
##pdf('./median_distance_to_class_neighbours_mammals_realms.pdf',width=10,height=7)
##hotspot.table.realm.all<-vector()
##par(mfrow=c(2,3))
##for (i in 1:nrow(grid.cells.df)){
##  cat(i,'\n')
##  #get coordinates of centroid of each cell
##  grid.coordinates<-lapply(grid.world,function(x) coordinates(x))
##  grid.coordinates<-do.call("rbind", grid.coordinates)
##  #get table for realm
##  hotspot.table.realm<-hotspots.table[(hotspots.table$cells%in%c(grid.cells.df$start.cell[i]:grid.cells.df$end.cell[i])),]
##  #get distances to neigbouring cells in a 1000 km radius
##  neighbours.1000<-dnearneigh(grid.coordinates[grid.cells.df$start.cell[i]:grid.cells.df$end.cell[i],],d1=0,d2=1000)
##  #get a matrix of distances of each cell in realm with neighbours
##  dist.realm.matrix<-as.matrix(dist(grid.coordinates[hotspot.table.realm$cells,],upper=TRUE))
##  #for each cell, check if hotspot or non hotspot, get neighbours that are the same class and average distance
##  hotspot.vector<-hotspot.table.realm$cells
##  names(hotspot.vector)<-hotspot.table.realm$hotspot
##  median.distance.class<-numeric()
##  for(a in 1:length(hotspot.vector)){
##    if(names(hotspot.vector[a])=='1'){
##      neighbours.1000.subseted<-neighbours.1000[[a]][neighbours.1000[[a]]%in%which(names(hotspot.vector)==1)]
##      median.distance.class[a]<-mean(dist.realm.matrix[a,neighbours.1000.subseted])
##    }else if(names(hotspot.vector[a])=='0'){
##      neighbours.1000.subseted<-neighbours.1000[[a]][neighbours.1000[[a]]%in%which(names(hotspot.vector)==0)]
##      median.distance.class[a]<-mean(dist.realm.matrix[a,neighbours.1000.subseted])
##    }
##  }
##  
##  #another way, storing all distances to plot them later
##  cell.distance.to.neighbours<-list()
##  for(a in 1:length(hotspot.vector)){
##    if(names(hotspot.vector[a])=='1'){
##      neighbours.1000.subseted<-neighbours.1000[[a]][neighbours.1000[[a]]%in%which(names(hotspot.vector)==1)]
##      cell.distance.to.neighbours[[a]]<-median(dist.realm.matrix[a,neighbours.1000.subseted])
##    }else if(names(hotspot.vector[a])=='0'){
##      neighbours.1000.subseted<-neighbours.1000[[a]][neighbours.1000[[a]]%in%which(names(hotspot.vector)==0)]
##      cell.distance.to.neighbours[[a]]<-median(dist.realm.matrix[a,neighbours.1000.subseted])
##    }
##  }
##  
##
##  hotspot.table.realm$median.distance.to.class.neighbours<-unlist(cell.distance.to.neighbours)
##  names(median.distance.class)<-hotspot.table.realm$cells
##  cat(nrow(hotspot.table.realm),'\n')
##  hotspot.table.realm.all<-rbind(hotspot.table.realm.all,hotspot.table.realm)
##  #dist.realm.hotspot.matrix<-as.matrix(dist(grid.coordinates[hotspot.table.realm[hotspot.table.realm$hotspot==1,'cells'],],upper=TRUE))
##  #median.dist.realm.hotspot<-apply(dist.realm.hotspot.matrix,1,function(x) median(x))
##  #mean.dist.realm.hotspot<-apply(dist.realm.hotspot.matrix,1,function(x) mean(x))
##  #dist.realm.nonhotspot.matrix<-as.matrix(dist(grid.coordinates[hotspot.table.realm[hotspot.table.realm$hotspot==0,'cells'],],upper=TRUE))
##  #median.dist.realm.nonhotspot<-apply(dist.realm.nonhotspot.matrix,1,function(x) median(x))
##  #mean.dist.realm.nonhotspot<-apply(dist.realm.nonhotspot.matrix,1,function(x) mean(x))
##  #hotspot.table.realm$dist.to.class<-NA
##  #hotspot.table.realm[hotspot.table.realm$hotspot==1,'dist.to.class']<-median.dist.realm.hotspot
##  #hotspot.table.realm[hotspot.table.realm$hotspot==0,'dist.to.class']<-median.dist.realm.nonhotspot
##  
##  #boxplot(median.distance.to.class.neighbours~hotspot,data=hotspot.table.realm,main=grid.cells.df[i,'realms.names'])
##  plot(c(1,1),x=c(0,2),ylim=c(0,1000),type='n',xaxt='n',ylab='median.distance.to.class.neighbours',xlab='',main=grid.cells.df[i,'realms.names'])
##  vioplot(hotspot.table.realm[hotspot.table.realm$hotspot==0,'median.distance.to.class.neighbours'][!is.na(hotspot.table.realm[hotspot.table.realm$hotspot==0,'median.distance.to.class.neighbours'])],at=0.5,col='blue',add=T)
##  vioplot(hotspot.table.realm[hotspot.table.realm$hotspot==1,'median.distance.to.class.neighbours'][!is.na(hotspot.table.realm[hotspot.table.realm$hotspot==1,'median.distance.to.class.neighbours'])],at=1.5,col='red',add=T)
##  summary.hotspots.neighbour.class.distances[[i]]<-summary(hotspot.table.realm[hotspot.table.realm$hotspot==1,'median.distance.to.class.neighbours'][!is.na(hotspot.table.realm[hotspot.table.realm$hotspot==1,'median.distance.to.class.neighbours'])])
##  summary.nonhotspots.neighbour.class.distances[[i]]<-summary(hotspot.table.realm[hotspot.table.realm$hotspot==0,'median.distance.to.class.neighbours'][!is.na(hotspot.table.realm[hotspot.table.realm$hotspot==0,'median.distance.to.class.neighbours'])])
##  
##  axis(1,at=c(0.5,1.5),labels=c('NH','H'))
##}
##dev.off()
##summary.hotspots.neighbour.class.distances<-as.data.frame(do.call('rbind',summary.hotspots.neighbour.class.distances))
##summary.hotspots.neighbour.class.distances$realm<-c('Afro','Austral','Indo','Nearctic','Neotrop','Palearctic')
##write.table(summary.hotspots.neighbour.class.distances,file='./output/mammals/tables/mediandistancetoclass_1000neighbours_hotspots.txt',sep='\t',quote=F,row.names=F)
##write.table(hotspot.table.realm.all,file='./output/mammals/tables/mediandistancetoclass_1000neighbours_hotspots_table.txt',sep='\t',quote=F,row.names=F)
##summary.nonhotspots.neighbour.class.distances<-as.data.frame(do.call('rbind',summary.nonhotspots.neighbour.class.distances))
##summary.nonhotspots.neighbour.class.distances$realm<-c('Afro','Austral','Indo','Nearctic','Neotrop','Palearctic')
##write.table(summary.nonhotspots.neighbour.class.distances,file='./output/mammals/tables/mediandistancetoclass_1000neighbours_nonhotspots.txt',sep='\t',quote=F,row.names=F)
##
##
##
#########for birds
##
##library(spdep)
###analyse by realms
##grids<-list.files(path='./output/grids/',pattern='grid_.*_100.rds')
##world.grid<-grep('World_RealmsMerged',grids)
##grid.world<-grids[grep('World_RealmsMerged',grids)]
##grid.realms.names<-grids[-grep('World_RealmsMerged',grids)]
##grid.world<-readRDS(paste('./output/grids/',grid.world,sep=''))
##grid.realms<-lapply(grid.realms.names,function(x) readRDS(paste('./output/grids/',x,sep='')))
##grid.realms.names<-sub(grid.realms.names,pattern='grid_',replacement='')
##grid.realms.names<-sub(grid.realms.names,pattern='_100.rds',replacement='')
##grid.realms.ncells<-unlist(lapply(grid.realms,function(x)length(x)))
##grid.realms.start.cell<-numeric(length(grid.realms.ncells))
##grid.realms.end.cell<-numeric(length(grid.realms.ncells))
##for(i in 1:length(grid.realms.ncells)){
##  if(i==1){
##    grid.realms.start.cell[i]<-1
##    grid.realms.end.cell[i]<-grid.realms.ncells[i]
##    next
##  }else{
##    grid.realms.start.cell[i]<-grid.realms.end.cell[i-1]+1
##    grid.realms.end.cell[i]<-grid.realms.start.cell[i]+grid.realms.ncells[i]-1
##    next
##  }
##}
##grid.cells.df<-as.data.frame(cbind(grid.realms.names,grid.realms.ncells,grid.realms.start.cell,grid.realms.end.cell),stringsAsFactors = F)
##colnames(grid.cells.df)<-c('realms.names','ncells','start.cell','end.cell')
##grid.cells.df$ncells<-as.numeric(grid.cells.df$ncells)
##grid.cells.df$start.cell<-as.numeric(grid.cells.df$start.cell)
##grid.cells.df$end.cell<-as.numeric(grid.cells.df$end.cell)
###remove Oceanic grid.cells.df (just 79 cells and most of them hotspot)
##grid.cells.df<-grid.cells.df[-grep('Oceanic',grid.cells.df$realms.names),]
###get the hotspots in
##hotspots.table<-read.table('./output/birds/tables/100_all_realms_number.of.species.wend_0.8.txt',header=T,sep='\t',stringsAsFactors = F)
##hotspots.table<-hotspots.table[,c('cells','hotspot')]
###for each realm, get coordinates of centroid
##median.distance.class.realm<-numeric()
##summary.hotspots.neighbour.class.distances<-list()
##summary.nonhotspots.neighbour.class.distances<-list()
##pdf('./median_distance_to_class_neighbours_birds_realms.pdf',width=10,height=7)
##par(mfrow=c(2,3))
##for (i in 1:nrow(grid.cells.df)){
##  cat(i,'\n')
##  #get coordinates of centroid of each cell
##  grid.coordinates<-lapply(grid.world,function(x) coordinates(x))
##  grid.coordinates<-do.call("rbind", grid.coordinates)
##  #get table for realm
##  hotspot.table.realm<-hotspots.table[(hotspots.table$cells%in%c(grid.cells.df$start.cell[i]:grid.cells.df$end.cell[i])),]
##  #get distances to neigbouring cells in a 1000 km radius
##  neighbours.1000<-dnearneigh(grid.coordinates[grid.cells.df$start.cell[i]:grid.cells.df$end.cell[i],],d1=0,d2=1000)
##  #get a matrix of distances of each cell in realm with neighbours
##  dist.realm.matrix<-as.matrix(dist(grid.coordinates[hotspot.table.realm$cells,],upper=TRUE))
##  #for each cell, check if hotspot or non hotspot, get neighbours that are the same class and average distance
##  hotspot.vector<-hotspot.table.realm$cells
##  names(hotspot.vector)<-hotspot.table.realm$hotspot
##  median.distance.class<-numeric()
##  for(a in 1:length(hotspot.vector)){
##    if(names(hotspot.vector[a])=='1'){
##      neighbours.1000.subseted<-neighbours.1000[[a]][neighbours.1000[[a]]%in%which(names(hotspot.vector)==1)]
##      median.distance.class[a]<-mean(dist.realm.matrix[a,neighbours.1000.subseted])
##    }else if(names(hotspot.vector[a])=='0'){
##      neighbours.1000.subseted<-neighbours.1000[[a]][neighbours.1000[[a]]%in%which(names(hotspot.vector)==0)]
##      median.distance.class[a]<-mean(dist.realm.matrix[a,neighbours.1000.subseted])
##    }
##  }
##  
##  #another way, storing all distances to plot them later
##  cell.distance.to.neighbours<-list()
##  for(a in 1:length(hotspot.vector)){
##    if(names(hotspot.vector[a])=='1'){
##      neighbours.1000.subseted<-neighbours.1000[[a]][neighbours.1000[[a]]%in%which(names(hotspot.vector)==1)]
##      cell.distance.to.neighbours[[a]]<-median(dist.realm.matrix[a,neighbours.1000.subseted])
##    }else if(names(hotspot.vector[a])=='0'){
##      neighbours.1000.subseted<-neighbours.1000[[a]][neighbours.1000[[a]]%in%which(names(hotspot.vector)==0)]
##      cell.distance.to.neighbours[[a]]<-median(dist.realm.matrix[a,neighbours.1000.subseted])
##    }
##  }
##  
##  
##  hotspot.table.realm$median.distance.to.class.neighbours<-unlist(cell.distance.to.neighbours)
##  names(median.distance.class)<-hotspot.table.realm$cells
##  cat(nrow(hotspot.table.realm),'\n')
##  
##  #dist.realm.hotspot.matrix<-as.matrix(dist(grid.coordinates[hotspot.table.realm[hotspot.table.realm$hotspot==1,'cells'],],upper=TRUE))
##  #median.dist.realm.hotspot<-apply(dist.realm.hotspot.matrix,1,function(x) median(x))
##  #mean.dist.realm.hotspot<-apply(dist.realm.hotspot.matrix,1,function(x) mean(x))
##  #dist.realm.nonhotspot.matrix<-as.matrix(dist(grid.coordinates[hotspot.table.realm[hotspot.table.realm$hotspot==0,'cells'],],upper=TRUE))
##  #median.dist.realm.nonhotspot<-apply(dist.realm.nonhotspot.matrix,1,function(x) median(x))
##  #mean.dist.realm.nonhotspot<-apply(dist.realm.nonhotspot.matrix,1,function(x) mean(x))
##  #hotspot.table.realm$dist.to.class<-NA
##  #hotspot.table.realm[hotspot.table.realm$hotspot==1,'dist.to.class']<-median.dist.realm.hotspot
##  #hotspot.table.realm[hotspot.table.realm$hotspot==0,'dist.to.class']<-median.dist.realm.nonhotspot
##  
##  #boxplot(median.distance.to.class.neighbours~hotspot,data=hotspot.table.realm,main=grid.cells.df[i,'realms.names'])
##  plot(c(1,1),x=c(0,2),ylim=c(0,1000),type='n',xaxt='n',ylab='median.distance.to.class.neighbours',xlab='',main=grid.cells.df[i,'realms.names'])
##  vioplot(hotspot.table.realm[hotspot.table.realm$hotspot==0,'median.distance.to.class.neighbours'][!is.na(hotspot.table.realm[hotspot.table.realm$hotspot==0,'median.distance.to.class.neighbours'])],at=0.5,col='blue',add=T)
##  vioplot(hotspot.table.realm[hotspot.table.realm$hotspot==1,'median.distance.to.class.neighbours'][!is.na(hotspot.table.realm[hotspot.table.realm$hotspot==1,'median.distance.to.class.neighbours'])],at=1.5,col='red',add=T)
##  summary.hotspots.neighbour.class.distances[[i]]<-summary(hotspot.table.realm[hotspot.table.realm$hotspot==1,'median.distance.to.class.neighbours'][!is.na(hotspot.table.realm[hotspot.table.realm$hotspot==1,'median.distance.to.class.neighbours'])])
##  summary.nonhotspots.neighbour.class.distances[[i]]<-summary(hotspot.table.realm[hotspot.table.realm$hotspot==0,'median.distance.to.class.neighbours'][!is.na(hotspot.table.realm[hotspot.table.realm$hotspot==0,'median.distance.to.class.neighbours'])])
##  
##  axis(1,at=c(0.5,1.5),labels=c('NH','H'))
##}
##dev.off()
##
##
##