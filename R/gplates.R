library(raster)
library(maptools)
library(colorspace)
library(fields)
library(rgdal)
library(sp)
library(rgeos)
library(maptools)
library(parallel)
library(PBSmapping)
library(ape)
library(geiger)
library(phytools)
library(phangorn)
library(diversitree)
library(cleangeo)
#read in the 0Ma_xy
list.xy.files<-list.files(path='/Users/javier/Documents/',pattern='*.xy')
#sort in correct order
list.xy.files<-list.xy.files[c(c(1,2,13,24,35,46,57,64,65,66),c(1:length(list.xy.files))[-c(1,2,13,24,35,46,57,64,65,66)])]
centroids.time<-list()
for(i in 1:length(list.xy.files)){
  cat(i,'\n')
  lines<-readLines(paste('/Users/javier/Documents/',list.xy.files[i],sep=''))
  lines<-lines[grep('^  ',lines)]
  lines<-sapply(lines,function(x) gsub(x,pattern=' +',replacement=' '))
  lines<-sapply(lines,function(x) gsub(x,pattern='^ ',replacement=''))
  lines<-lapply(lines,function(x)as.numeric(unlist(strsplit(x,split=' '))))
  centroids.time[[i]]<-as.data.frame(matrix(unlist(lapply(c(seq(from=1,to=length(lines),by=5)),function(x) c(mean(c(unlist(lines[x])[1],unlist(lines[x+1])[1],unlist(lines[x+2])[1],unlist(lines[x+3])[1])),mean(c(unlist(lines[x])[2],unlist(lines[x+1])[2],unlist(lines[x+2])[2],unlist(lines[x+3])[2]))))),ncol=2,byrow=TRUE),stringsAsFactors=TRUE)
}

present.points<-centroids.time[[1]]
colnames(present.points)<-c('x','y')
plot(present.points,pch=16,cex=.3,col='blue')
past.points<-centroids.time[[66]]
plot(past.points,pch=16,cex=.3,col='red')
#get coordinates of centroids of cells in grid
grid<-readRDS('./output/grids/grid_World_RealmsMerged_100.rds')
grid<-lapply(grid,function(x) spTransform(x,CRS("+proj=longlat")))
#save to raster
grid.df<-lapply(grid,function(x)coordinates(x))
grid.df<-as.data.frame(do.call('rbind',grid.df),stringsAsFactors=F)
colnames(grid.df)<-c('x','y')
#calculate distances between 65 and 0 myr
dist.matrix<-as.data.frame(matrix(NA,ncol=1,nrow=nrow(past.points)),stringsAsFactors=F)
for(i in 1:nrow(past.points)){
  cat(i,'\n')
  dist.matrix[i,]<-dist(rbind(present.points[i,],past.points[i,]))[1]
}
hist(dist.matrix$V1)
plot(present.points,pch=16,cex=.1)
#check the distance between the grid cell centroid and the present point assigned
system.time(cell.distance.to.points<-unlist(sapply(c(1:nrow(grid.df)),function(x)dist(rbind(grid.df[x,],present.points[x,]))[1])))
hist(cell.distance.to.points)
cell.distance.to.points.colour<-create_colours_vector(ncategories = 10,vector=cell.distance.to.points)
#no obvious bias anywhere
plot(grid.df[,'x'],grid.df[,'y'],pch=16,cex=.5,col=cell.distance.to.points.colour)
library(maptools)
library(spdep)
grid<-readRDS('./output/grids/grid_World_RealmsMerged_100.rds')
grid<-lapply(grid,function(x) spTransform(x,CRS("+proj=longlat")))
grid.list<-list()
for (i in 1:length(grid)){
  cat(i,'\n')
  if(i==1){
    grid.list<-grid[[i]]
  }else{
    grid.list<-rbind(grid.list,grid[[i]],makeUniqueIDs = TRUE) 
  }
  
}
centroids.time<-lapply(centroids.time,function(x){colnames(x)<-c('x','y');return(x)})
#get queen neighbours
grid.queen<-poly2nb(grid.list, queen=TRUE)
#plot number of neighbours
n.of.neighbours<-unlist(lapply(grid.queen,function(x)length(x)))
n.of.neighbours.colour<-sapply(n.of.neighbours,function(x) diverge_hsv(n=8)[x])
plot(grid.df[,'x'],grid.df[,'y'],pch=16,cex=.5,col=n.of.neighbours.colour)

grid.queen.distances<-list()
grid.queen.distances.stdev<-list()
for (i in 1:length(grid.queen)){
  cat(i,'\n')
  grid.queen.distances[[i]]<-lapply(grid.queen[[i]],function(x) lapply(centroids.time,function(y) dist(rbind(y[i,],y[x,]))[1]))
  grid.queen.distances.stdev[[i]]<-lapply(grid.queen.distances[[i]],function(x)sd(unlist(x)))
}
grid.queen.distances.mean.stdev<-unlist(lapply(grid.queen.distances.stdev,function(x) mean(unlist(x))))
grid.queen.distances.mean.stdev.colour<-create_colours_vector(ncategories = 10,vector=grid.queen.distances.mean.stdev)
grid.df$tectonic.movement<-grid.queen.distances.mean.stdev/100
hist(grid.df$tectonic.movement)
grid.queen.distances.mean.stdev.colour<-create_colours_vector(ncategories = 10,vector=grid.df$tectonic.movement)
plot(grid.df[,'x'],grid.df[,'y'],pch=16,cex=.4,col=grid.queen.distances.mean.stdev.colour)
#output to table
tectonic.output<-as.data.frame(cbind(c(1:nrow(grid.df)),grid.df$tectonic.movement),stringsAsFactors=F)
colnames(tectonic.output)<-c('cells','tectonic.movement')
write.table(tectonic.output,'./output/grids/tables/100_all_realms_tectonicmovement.txt',row.names=F,quote=F,sep='\t')
