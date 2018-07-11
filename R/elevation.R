library(raster)
library(parallel)

####################################################################################
#this function builds a dem object with a list of objects for the whole globe
#get info for elevation data, dem will be stored in 'world_dem.Rsave'
#all world dems are stored in ./raw_data/gtopo30_dems/
build_world_dem<-function(path){
  demfiles<-list.files(path,pattern='\\.dem')
  #not loading 1 which is Antarctica
  dem.list<-lapply(demfiles[2:length(demfiles)],function(x) raster(paste('./raw_data/gtopo30_dems/',x,sep='')))
  #merge all dems into one
  new.dem<-dem.list[[1]]
  for (i in 2:length(dem.list)){
    new.dem<-merge(new.dem,dem.list[[i]])
  }
  writeRaster(new.dem, file="./output/grids/world_dem")
}

####################################################################
#this function extracts data from a raster for a grid
#layer is the raster that contains the data
#x is the number of polygon in the grid
#grid is the grid object
extract_data_grid<-function(layer,x,grid){
  cat('getting values for cell ',x,'\n')
  data<-as.data.frame(extract(layer,spTransform(grid[[x]],CRS=proj4string(layer))))
  data<-na.omit(data)
  data<-data[,1]
  return(data)
}
##########################################



get_elevation_variables_grid<-function(gridfile,dem,path){
  grid<-readRDS(gridfile)
  cat('extracting mean elevation...','\n')
  mean.elevation<-unlist(lapply(c(1:length(grid)),function(x) mean(extract_data_grid(dem,x,grid))))
  cat('extracting range of elevation...','\n')
  range.elevation<-unlist(lapply(c(1:length(grid)),function(x) {data<-extract_data_grid(dem,x,grid);if(length(data)==0){range<-NaN}else{range<-range(data);range<-range[2]-range[1]};return(range)}))
  cat('extracting mean TRI','\n')  
  mean.TRI<-unlist(lapply(c(1:length(grid)),function(x) mean(extract(terrain(crop(dem, extent(spTransform(grid[[x]],CRS=proj4string(dem)))),opt='TRI'),extent(spTransform(grid[[x]],CRS=proj4string(dem)))),na.rm=TRUE)))
  cat('extracting var TRI','\n')
  var.TRI<-unlist(lapply(c(1:length(grid)),function(x) var(extract(terrain(crop(dem, extent(spTransform(grid[[x]],CRS=proj4string(dem)))),opt='TRI'),extent(spTransform(grid[[x]],CRS=proj4string(dem)))),na.rm=TRUE)))
  elevation.table<-as.data.frame(cbind(c(1:length(grid)),mean.elevation,range.elevation,mean.TRI,var.TRI),stringAsFactors=F)
  colnames(elevation.table)[1]<-'cells'
  write.table(elevation.table,file=paste(path,'100_all_realms_elevation_table.txt',sep=''),sep='\t',row.names=F,quote=F)
}

####################################################################################
##plot meanTRI
##non hotspot vs hotspot
#hotspot0<-subset(grid.table.subset,grid.table.subset$hotspots=='0')
#hotspot1<-subset(grid.table.subset,grid.table.subset$hotspots!='0')
#
#plot(density(hotspot0$mean.TRI,na.rm=TRUE),xaxs='i',yaxs='i',main='mean.TRI',col='blue',xlab='mean.TRI')
#lines(density(hotspot1$mean.TRI,na.rm=TRUE),col='red')
#legend ('topright',c('nonhotspot','hotspot'),col = c('blue','red'), lty=1,lwd=2, bty="n", cex=.6)
#
##different means
#wilcox.test(hotspot0$mean.TRI,hotspot1$mean.TRI)
##ranks
#hotspot0$mean.TRI.01<-rank(hotspot0$mean.TRI)/length(hotspot0$mean.TRI)
#hotspot1$mean.TRI.01<-rank(hotspot1$mean.TRI)/length(hotspot1$mean.TRI)
#plot(hotspot0$mean.TRI,hotspot0$mean.TRI.01,pch=16,col='blue',main='rank of mean.TRI',xlab='mean.TRI',ylab='rank')
#points(hotspot1$mean.TRI,hotspot1$mean.TRI.01,pch=16,col='red')
#legend ('bottomright',c('nonhotspot','hotspot'),col = c('blue','red'), lty=1,lwd=2, bty="n", cex=.6)
#
##plot varTRI
#plot(density(hotspot0$var.TRI,na.rm=TRUE),xaxs='i',yaxs='i',main='var.TRI',col='blue',xlab='var.TRI')
#lines(density(hotspot1$var.TRI,na.rm=TRUE),col='red')
#legend ('topright',c('nonhotspot','hotspot'),col = c('blue','red'), lty=1,lwd=2, bty="n", cex=.6)
#
##different means
#wilcox.test(hotspot0$var.TRI,hotspot1$var.TRI)
##ranks
#hotspot0$var.TRI.01<-rank(hotspot0$var.TRI)/length(hotspot0$var.TRI)
#hotspot1$var.TRI.01<-rank(hotspot1$var.TRI)/length(hotspot1$var.TRI)
#plot(hotspot0$var.TRI,hotspot0$var.TRI.01,pch=16,col='blue',main='rank of var.TRI',xlab='var.TRI',ylab='rank')
#points(hotspot1$var.TRI,hotspot1$var.TRI.01,pch=16,col='red')
#legend ('bottomright',c('nonhotspot','hotspot'),col = c('blue','red'), lty=1,lwd=2, bty="n", cex=.6)
#
##############do permutations for the meanTRI
##omit one cell that has no elevation data (it's just off the coast)
#grid.table.subset<-na.omit(grid.table.subset)
##also omit the row and the column for that cell
#distances.subset.corrected<-distances.subset.corrected[-27,]
#distances.subset.corrected<-distances.subset.corrected[,-27]
#
#results.mean<-matrix(NA,nrow=1000,ncol=2)
#for(i in 1:1000){
#  new.rates<-vector()
#  #for every cell get a random rate from the rate vector (weighted by the probability of distance of cell to rest of the cells)
#  for (a in 1:length(grid.table.subset$mean.TRI)){
#    new.rates<-c(new.rates,sample(grid.table.subset$mean.TRI,1,prob=distances.subset.corrected[a,]))
#  }
#  #keep the same hotspots tag
#  new.rates.df<-data.frame(new.rates,grid.table.subset$hotspots)
#  colnames(new.rates.df)<-c('mean.TRI','hotspots')
#  cat(i,'\n')
#  results.mean[i,]<-c(i,mean(new.rates.df[!(new.rates.df$hotspots=='0'),]$mean.TRI)-mean(new.rates.df[new.rates.df$hotspots=='0',]$mean.TRI))
#}
#results.mean<-as.data.frame(results.mean)
#colnames(results.mean)<-c('rep','diff.mean.TRI.hotspots.non')
#hist(results.mean$diff.mean.TRI.hotspots.non,main='diff(hot vs non) in mean of mean.TRI',yaxs='i',xaxs='i',xlab='mean hot - mean nonhot',xlim=c(-20,20))
#abline(v=mean(grid.table.subset[!(grid.table.subset$hotspots=='0'),]$mean.TRI)-mean(grid.table.subset[grid.table.subset$hotspots=='0',]$mean.TRI),col='red')
#
#
##############do permutations for the varTRI
#results.mean<-matrix(NA,nrow=1000,ncol=2)
#for(i in 1:1000){
#  new.rates<-vector()
#  #for every cell get a random rate from the rate vector (weighted by the probability of distance of cell to rest of the cells)
#  for (a in 1:length(grid.table.subset$mean.TRI)){
#    new.rates<-c(new.rates,sample(grid.table.subset$var.TRI,1,prob=distances.subset.corrected[a,]))
#  }
#  #keep the same hotspots tag
#  new.rates.df<-data.frame(new.rates,grid.table.subset$hotspots)
#  colnames(new.rates.df)<-c('var.TRI','hotspots')
#  cat(i,'\n')
#  results.mean[i,]<-c(i,mean(new.rates.df[!(new.rates.df$hotspots=='0'),]$var.TRI)-mean(new.rates.df[new.rates.df$hotspots=='0',]$var.TRI))
#}
#results.mean<-as.data.frame(results.mean)
#colnames(results.mean)<-c('rep','diff.var.TRI.hotspots.non')
#hist(results.mean$diff.var.TRI.hotspots.non,main='diff(hot vs non) in mean of var.TRI',xlim=c(-1210,1210),yaxs='i',xaxs='i',xlab='mean hot - mean nonhot')
#abline(v=mean(grid.table.subset[!(grid.table.subset$hotspots=='0'),]$var.TRI)-mean(grid.table.subset[grid.table.subset$hotspots=='0',]$var.TRI),col='red')
#
#hotspots.afro<-c("afromontane","cape","coasteastafrica","guinea","horn","karoo","maputaland")
#
#############buildling "control sets of polygons" (i.e, non hotspot regions of same size as hotspots)
####subset the distance matrix (for only the cells with species in the grid)
#distances.subset<-distances[grid.table.subset$cell,grid.table.subset$cell]
####subset the matrix for non hotspot cells only
#distances.subset.non<-distances[grid.table.subset[(grid.table.subset$hotspots=='0'),]$cell,grid.table.subset[(grid.table.subset$hotspots=='0'),]$cell]
####get the total number of cells in the target hotspot (7 in the example)
#desired.size<-length(grep(7,grid.table.subset$hotspots))
##x is a row (the distance of one cell to the rest of the cells)
#control.cells<-list()
#for (row in 1:nrow(distances.subset.non)){
#  x<-distances.subset.non[row,]
#  distance.vector<-vector()
#  position.vector<-vector()
#  for (i in 1:desired.size){
#    distances<-vector()
#    positions<-vector()
#    distances<-x[as.numeric(which(x==i))]
#    if (length(distances)==0){
#      cat('unable to make a continuous polygon with distance ',i,'\n')
#      break
#    }
#    positions<-as.numeric(names(distances))
#    distances<-unname(unlist(distances))
#    distance.vector<-c(distance.vector,distances)
#    position.vector<-c(position.vector,positions)
#    if(length(position.vector)>=desired.size){
#      distance.vector<-distance.vector[1:desired.size]
#      position.vector<-position.vector[1:desired.size]
#      control.cells[[row]]<-position.vector
#      break
#    }
#    
#  }
#  
#}
#
#grid.table.subset.hotspot<-grid.table.subset[grep(7,grid.table.subset$hotspots),]
#mean(grid.table.subset.hotspot$var.TRI)
#control.results<-data.frame(matrix(NA,nrow=length(control.cells),ncol=2))
#for (i in 1:length(control.cells)){
#  control.subset<-grid.table.subset[which(grid.table.subset$cell %in% control.cells[[i]]),]
#  control.results[i,]<-c(i,mean(control.subset$var.TRI))
#}
#colnames(control.results)<-c('number.control','mean.variable')
#plot(density(control.results$mean.variable,na.rm=TRUE),main='var.TRI',xaxs='i',yaxs='i',xlab='rate',xlim=c(0,3500))
#abline(v=mean(grid.table.subset.hotspot$var.TRI),col='red')
#legend ('topright',c(hotspots.afro[7]),col = c('red'), lty=1,lwd=2, bty="n", cex=.6)
#
#
#