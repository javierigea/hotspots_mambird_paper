library(spdep)
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
######just for plotting purposes, divides a vector into quartiles and assigns colours
create_colours_vector<-function(ncategories,vector){
  colours.vector<-diverge_hsv(n=ncategories)
  quantiles.colour<-quantile(vector, seq(from=0,to=1,by=(1/ncategories)),na.rm=TRUE)
  grid.colours<-vector()
  #assign colours to each element of vector 
  for (i in 1:length(vector)){
    if(is.na(vector[i])){
      grid.colours[i]<-'grey'
      next;
    }
    new.vector<-sort(c(quantiles.colour,vector[i]))
    interval<-which(new.vector==vector[i])-1
    if(length(interval)>1){
      interval<-interval[1]
    }
    if(interval==0){
      interval<-1
    }
    grid.colours[i]<-colours.vector[interval]
  }
  return(grid.colours)
}

###########this functions creates a shp from the world grid file to be read in by Gplates
create_shp_from_grid<-function(world.gridfile,path){
  grid<-readRDS(world.gridfile)
  grid<-lapply(grid,function(x) spTransform(x,CRS("+proj=longlat")))
  #save to shp
  grid.list<-list()
  cat('binding grid cells','\n')
  for (i in 1:length(grid)){
    if(i==1){
      grid.list<-grid[[i]]
    }else{
      grid.list<-rbind(grid.list,grid[[i]],makeUniqueIDs = TRUE) 
    }
    
  }
  currentwd<-getwd()
  setwd(path)
  writeSpatialShape(x = as(grid.list, "SpatialPolygonsDataFrame" ),fn = 'all_worldgrid')
  setwd(currentwd)
}
#this takes a path with reconstructed positions (x,y) of a set of points (that determine the cells in a grid)
#and outputs a table with the average standard deviations of the distance of a given cell with all its queen neighbours
#across time. It first discards any cell that hasn't moved in 65 My (= errors when mapping to Gplates)
#path is the path where the output table will be saved
get_average_stdev_tectonicmovements<-function(reconstruction.path,world.gridfile,path){
  #read in the 0Ma_xy
  list.xy.files<-list.files(path=reconstruction.path,pattern='*.xy')
  #sort in correct order
  list.xy.files<-list.xy.files[c(c(1,2,13,24,35,46,57,64,65,66),c(1:length(list.xy.files))[-c(1,2,13,24,35,46,57,64,65,66)])]
  centroids.time<-list()
  #read in all reconstruction files, clean the format and get the coordinates of the centroids of each cell
  cat('process reconstruction files ')
  for(i in 1:length(list.xy.files)){
    cat(i,', ')
    lines<-readLines(paste(reconstruction.path,list.xy.files[i],sep=''))
    #clean to get just coordinates (one point; x and y columns; per line)
    lines<-lines[grep('^  ',lines)]
    lines<-sapply(lines,function(x) gsub(x,pattern=' +',replacement=' '))
    lines<-sapply(lines,function(x) gsub(x,pattern='^ ',replacement=''))
    lines<-lapply(lines,function(x)as.numeric(unlist(strsplit(x,split=' '))))
    #there are 5 coordinates per cell (the last one is the same as the 4th one) so I read batches of 5 skipping the 5th one
    #to get the centroid, I take the mean of all x and all y coordinates in the vertices
    centroids.time[[i]]<-as.data.frame(matrix(unlist(lapply(c(seq(from=1,to=length(lines),by=5)),function(x) c(mean(c(unlist(lines[x])[1],unlist(lines[x+1])[1],unlist(lines[x+2])[1],unlist(lines[x+3])[1])),mean(c(unlist(lines[x])[2],unlist(lines[x+1])[2],unlist(lines[x+2])[2],unlist(lines[x+3])[2]))))),ncol=2,byrow=TRUE),stringsAsFactors=TRUE)
  }
  present.points<-centroids.time[[1]]
  colnames(present.points)<-c('x','y')
  plot(present.points,pch=16,cex=.3,col='blue')
  past.points<-centroids.time[[66]]
  colnames(past.points)<-c('x','y')
  plot(past.points,pch=16,cex=.3,col='red')
  #calculate distances between 65 and 0 myr
  dist.matrix<-as.data.frame(matrix(NA,ncol=1,nrow=nrow(past.points)),stringsAsFactors=F)
  cat('calculating distance betweeen present and 65 Ma','\n')
  for(i in 1:nrow(past.points)){
    dist.matrix[i,]<-dist(rbind(present.points[i,],past.points[i,]))[1]
  }
  hist(dist.matrix$V1)
  plot(present.points,pch=16,cex=.2)
  #these are the points that didn't move between 0 and 65 Ma (errors of Gplates mapping)
  points(present.points[which(dist.matrix$V1==0),],pch=16,cex=.3,col='red')
  #store points that haven't moved
  static.points<-which(dist.matrix$V1==0)
  #get coordinates of centroids of cells in grid
  grid<-readRDS(world.gridfile)
  grid<-lapply(grid,function(x) spTransform(x,CRS("+proj=longlat")))
  grid.df<-lapply(grid,function(x)coordinates(x))
  grid.df<-as.data.frame(do.call('rbind',grid.df),stringsAsFactors=F)
  colnames(grid.df)<-c('x','y')
  #check the distance between the grid cell centroid and the present point assigned
  system.time(cell.distance.to.points<-unlist(sapply(c(1:nrow(grid.df)),function(x)dist(rbind(grid.df[x,],present.points[x,]))[1])))
  hist(cell.distance.to.points)
  cell.distance.to.points.colour<-create_colours_vector(ncategories = 10,vector=cell.distance.to.points)
  #no obvious bias anywhere
  plot(grid.df[,'x'],grid.df[,'y'],pch=16,cex=.5,col=cell.distance.to.points.colour)
  #add column names to centroids.time
  centroids.time<-lapply(centroids.time,function(x){colnames(x)<-c('x','y');return(x)})
  #drop static points in grid.df and in centroids.time
  grid.df<-grid.df[-static.points,]
  centroids.time<-lapply(centroids.time,function(x){return(x[-static.points,])})
  #to get queen neighbours I need all the polygons in a list object
  grid<-readRDS('./output/grids/grid_World_RealmsMerged_100.rds')
  grid<-grid[-static.points]
  grid<-lapply(grid,function(x) spTransform(x,CRS("+proj=longlat")))
  grid.list<-list()
  cat('aggregating cells into a polygon','\n')
  for (i in 1:length(grid)){
   cat(i,'\n')
    if(i==1){
      grid.list<-grid[[i]]
    }else{
      grid.list<-rbind(grid.list,grid[[i]],makeUniqueIDs = TRUE) 
    }
    
  }
  #get queen neighbours of each cell
  grid.queen<-poly2nb(grid.list, queen=TRUE)
  #plot number of neighbours
  n.of.neighbours<-unlist(lapply(grid.queen,function(x)length(x)))
  n.of.neighbours.colour<-sapply(n.of.neighbours,function(x) diverge_hsv(n=8)[x])
  plot(grid.df[,'x'],grid.df[,'y'],pch=16,cex=.5,col=n.of.neighbours.colour)
  #get the standard deviation of the distances to all the neighbouring cells for each cell across time
  grid.queen.distances<-list()
  grid.queen.distances.stdev<-list()
  for (i in 1:length(grid.queen)){
    cat(i,'\n')
    #for each cell in the grid, get the distances through time for each neighbour
    grid.queen.distances[[i]]<-lapply(grid.queen[[i]],function(x) lapply(centroids.time,function(y) dist(rbind(y[i,],y[x,]))[1]))
    #get the stdev of the distances
    grid.queen.distances.stdev[[i]]<-lapply(grid.queen.distances[[i]],function(x)sd(unlist(x)))
  }
  grid.queen.distances.mean.stdev<-unlist(lapply(grid.queen.distances.stdev,function(x) mean(unlist(x))))
  grid.df$tectonic.movement<-grid.queen.distances.mean.stdev/100
  hist(grid.df$tectonic.movement)
  grid.queen.distances.mean.stdev.colour<-create_colours_vector(ncategories = 10,vector=grid.df$tectonic.movement)
  plot(grid.df[,'x'],grid.df[,'y'],pch=16,cex=.4,col=grid.queen.distances.mean.stdev.colour)
  #output to table
  tectonic.output<-as.data.frame(cbind(c(1:nrow(grid.df)),grid.df$tectonic.movement),stringsAsFactors=F)
  colnames(tectonic.output)<-c('cells','tectonic.movement')
  write.table(tectonic.output,file=paste(path,'100_all_realms_tectonicmovement.txt',sep=''),row.names=F,quote=F,sep='\t')
}




