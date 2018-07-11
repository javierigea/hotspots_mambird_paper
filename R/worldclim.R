library(raster)
library(rgdal)
library(sp)
library(rgeos)
library(maptools)
library(multicore)
library(PBSmapping)
library(plyr)
library(ade4)
library(factoextra)
library(parallel)
library(ape)
library(BAMMtools)
library(cleangeo)
#####################################################################
#function to get worlclim data at 2.5 and store it in raw_data folder
#get the worldclim data
fetch_worldclim<-function(){
  climate<-getData('worldclim', var='bio', res=2.5,path='./raw_data/')
}


####################################################################
#this function extracts all variables of data from a raster for a grid
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
get_variable_grid<-function(gridfile,raster,name,path){
  grid<-readRDS(gridfile)
  cat('extracting mean and variance of variable','\n')  
  mean.var<-lapply(c(1:length(grid)),function(x) {data.grid<-extract_data_grid(raster,x,grid);mean<-mean(data.grid);var<-var(data.grid);return(c(mean,var))})
  #cat('extracting var variable','\n')  
  #var<-unlist(lapply(c(1:length(grid)),function(x) var(extract_data_grid(raster,x,grid))))
  mean<-unlist(lapply(mean.var,function(x)x[1]))
  var<-unlist(lapply(mean.var,function(x)x[2]))
  grid.table<-as.data.frame(cbind(c(1:length(grid)),mean,var),stringAsFactors=F)
  colnames(grid.table)[1]<-'cells'
  write.table(grid.table,file=paste(path,name,'_table.txt',sep=''),quote=F,sep='\t',row.names=F)
}
######################################
get_worldclim_variables_grid<-function(gridfile,gridtablefile,layer){
  load(gridfile)
  grid.table<-read.table(gridtablefile,header=T,sep='\t')
  colnames(grid.table)[1]<-'cells'
  means<-unlist(lapply(grid.table$cells,function(x) apply(as.data.frame(extract_data_grid(layer,x,grid)),2,function(x)mean(x))))
  means<-as.data.frame(matrix(means,ncol=19,byrow=TRUE))
  colnames(means)<-paste(rep('mean.bio',19),c(1:19),sep='')
  vars<-unlist(lapply(grid.table$cells,function(x) apply(as.data.frame(extract_data_grid(layer,x,grid)),2,function(x)var(x))))
  vars<-as.data.frame(matrix(vars,ncol=19,byrow=TRUE))
  colnames(vars)<-paste(rep('var.bio',19),c(1:19),sep='')
  grid.table<-cbind(grid.table,means,vars)
  write.table(grid.table,file=sub(gridtablefile,pattern='.txt',replacement='_worldclim_meansvars.txt'),sep='\t',row.names=F,quote=F)
}

########################################
#this function transforms a row of bioclim data (19 variables) to its principal components using a dudi.pca.object

get_pca_from_object<-function(pca.object,points){
  results<-apply(pca.object$c1,2,function(x) sum(x*((points-pca.object$cent)/pca.object$norm)))
  return(results)
}

########################################
#this function transforms a row of bioclim data (19 variables) to its principal components using a dudi.pca.object

merge_tables_bioclim<-function(vector.tablefiles,path,name){
  names.variables<-sapply(vector.tablefiles,function(x) {y<-gsub(x,pattern='.*/',replacement='');y<-gsub(y,pattern='_table.txt',replacement='');return(y)})
  tables<-lapply(vector.tablefiles,function(x) read.table(x,header=T,sep='\t',stringsAsFactors = F))
  for (i in 1:length(tables)){
    colnames(tables[[i]])<-c('cells',paste(colnames(tables[[i]][2]),names.variables[i],sep='.'),paste(colnames(tables[[i]][3]),names.variables[i],sep='.'))
  }
  merged.data.frame<-Reduce(function(...) merge(..., all=T), tables)
  write.table(merged.data.frame,file=paste(path,name,'.txt',sep=''),sep='\t',quote=F,row.names=F)
}

