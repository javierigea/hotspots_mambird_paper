library(SDMTools)
library(raster)
library(rgdal)

#path_current_raster
#path_LGM_raster
#type: if 'temperature' divide raster by 10 (to get degrees C)
#name: for output file
generate_LGMCCV_rasters<-function(path_current_raster,path_LGM_raster,type,name){
  cat('reading rasters','\n')
  current<-raster(path_current_raster)
  LGM<-raster(path_LGM_raster)
  if (type =='temperature'){
    current<-current/10
    LGM<-LGM/10
  }
  current<-projectRaster(current,crs = '+proj=cea +units=km +ellps=WGS84')
  LGM<-projectRaster(LGM,crs = '+proj=cea +units=km +ellps=WGS84')
  #stack both layers 
  cat('calculating temporal gradient','\n')
  stack<-stack(current,LGM)
  #substract current - LGM to get difference in degrees ("temperature gradient")
  overlay<-overlay(stack[[1]],stack[[2]],fun=function(x,y){x-y})
  #calculate slope ("spatial gradient")
  cat('calculating spatial gradient','\n')
  current.slope<-slope(current, latlon = F)
  #remove values < 0.01
  current.slope[current.slope<0.01]<-0.01
  #stack the temperature and spatial gradients and divide them
  stack.gradients<-stack(overlay,current.slope)
  overlay.gradients<-overlay(stack.gradients[[1]],stack.gradients[[2]],fun=function(x,y){x/y})
  #this gives raster units as m/year
  overlay.gradients<-(overlay.gradients*1000)/21000
  writeRaster(overlay.gradients,filename=paste('./output/rasters/',name,sep=''))
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

get_variable_grid<-function(gridfile,raster,name,path){
  grid<-readRDS(gridfile)
  cat('extracting mean and variance of variable','\n')  
  mean.var<-lapply(c(1:length(grid)),function(x) {data.grid<-extract_data_grid(raster,x,grid);mean<-mean(data.grid);var<-var(data.grid);return(c(mean,var))})
  #cat('extracting var variable','\n')  
  #var<-unlist(lapply(c(1:length(grid)),function(x) var(extract_data_grid(raster,x,grid))))
  mean<-unlist(lapply(mean.var,function(x)x[1]))
  var<-unlist(lapply(mean.var,function(x)x[2]))
  grid.table<-as.data.frame(cbind(c(1:length(grid)),mean,var),stringAsFactors=F)
  colnames(grid.table)<-c('cells',paste('mean.',name,sep=''),paste('var',name,sep=''))
  write.table(grid.table,file=paste(path,name,'_table.txt',sep=''),quote=F,sep='\t',row.names=F)
}


#get_variable_grid<-function(gridfile,raster,path,name){
#  grid<-readRDS(gridfile)
#  grid.table<-as.data.frame(matrix(NA,nrow=length(grid),ncol=2),stringsAsFactors=F)
#  grid.table$mean<-unlist(lapply(1:length(grid),function(x) mean(extract_data_grid(raster,x,grid))))
#  #grid.table$var<-unlist(lapply(1:length(grid),function(x) var(extract_data_grid(raster,x,grid))))
#  write.table(grid.table,file=paste(path,name,'_table.txt',sep=''),quote=F,sep='\t',row.names=F)
#}


