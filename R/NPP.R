library(raster)

####################################################################################
#this function builds a raster object with the NPP data from Cramer 1999 (stored in ave2npp.zo1)
#raster
build_raster_npp<-function(tablefile){
  npp.table<-read.csv(tablefile,header=F)
  colnames(npp.table)<-c('lon','lat','value')
  #three cells have missing values, represented by -9999.9 in the dataset
  npp.table[npp.table$value==-9999.9,]$value<-NA
  NPP.raster<-rasterFromXYZ(npp.table)
  #assign coord.ref
  proj4string(NPP.raster)<-"+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"
  save(NPP.raster, file="./output/world_NPP.Rsave")
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

get_NPP_variables_grid<-function(gridfile,gridtablefile,raster){
  load(gridfile)
  grid.table<-read.table(gridtablefile,header=T,sep='\t')
  colnames(grid.table)[1]<-'cells'
  grid.table$mean.NPP<-unlist(lapply(grid.table$cells,function(x) mean(extract_data_grid(raster,x,grid))))
  grid.table$var.NPP<-unlist(lapply(grid.table$cells,function(x) var(extract_data_grid(raster,x,grid))))
  write.table(grid.table,file=sub(gridtablefile,pattern='.txt',replacement='_NPP.txt'),sep='\t',row.names=F,quote=F)
}

####################################################################################
