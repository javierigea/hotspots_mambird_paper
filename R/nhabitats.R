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
library(cleangeo)

### number of habitats data
get_habitat_count_grid<-function(gridfile,raster,name,path){
  grid<-readRDS(gridfile)
  n.habitats<-unlist(lapply(c(1:length(grid)),function(x) length(unique(extract_data_grid(raster,x,grid)))))
  grid.table<-as.data.frame(cbind(c(1:length(grid)),n.habitats),stringAsFactors=F)
  write.table(grid.table,file=paste(path,name,'_table.txt',sep=''),quote=F,sep='\t',row.names=F)
}
