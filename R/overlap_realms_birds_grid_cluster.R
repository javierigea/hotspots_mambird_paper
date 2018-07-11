#run model in hydrogen cluster
.libPaths('/home/ecosys/ji247/R/x86_64-pc-linux-gnu-library/3.3')
setwd('/home/ecosys/ji247/hotspots_vertebrates/')

library(rgdal)
library(sp)
library(rgeos)
library(maptools)
library(parallel)
library(PBSmapping)
library(cleangeo)
library(R.utils)
library(geiger)
library(rgdal)
library(sp)
library(rgeos)
library(maptools)
library(parallel)
library(PBSmapping)
library(plyr)
library(ape)
library(cleangeo)
library(R.utils)
library(raster)

generate_birds_spp.ranges<-function(ncores,timeout.seconds){
  ##load spp ranges and process
  cat('opening layers','\n')
  shp.list<-list.files(path='./raw_data/birds/birds_ranges_2016/',pattern='*.shp',recursive = TRUE)
  for (i in seq(from=1, to=length(shp.list), by=500)){
    start<-i
    end<-i+499
    if (end > length(shp.list)){
      end<-length(shp.list)
    } 
    shp.list.small<-shp.list[start:end]
    message<-paste('reading from ',start,' to ',end,'\n',sep='')
    write(message,file=paste('./birds_control_',timeout.seconds,'.txt',sep=''),append = T)
    time<-system.time(spp.ranges<-mclapply(mc.cores=ncores,shp.list.small, function(x) try(withTimeout(readOGR(sub(strsplit(x,'/')[[1]][1],pattern='^',replacement='./raw_data/birds/birds_ranges_2016/'),sub(strsplit(x,'/')[[1]][2],pattern='.shp',replacement='')),timeout=timeout.seconds,onTimeout='error'))))
    write(time,file=paste('./birds_control_',timeout.seconds,'.txt',sep=''),append = T)
    #selecting breeding regions only
    breeding.regions<-c(1,2,5)
    #clean errors
    errors<-vector('numeric')
    for (f in 1:length(spp.ranges)){
      if (class(spp.ranges[[f]])!='SpatialPolygonsDataFrame'){
        errors<-c(errors,f)
      }
    }
    
    if(length(errors)>0){
      spp.ranges<-spp.ranges[-errors]
      cat(paste('some spp layers were not read: ',errors,'\n'))
      error.species<-shp.list.small[errors]
      write(error.species,file=paste('./birds_errors_',timeout.seconds,'.txt',sep=''),append = T)
      
    }
    spp.ranges<-mclapply(mc.cores=ncores,spp.ranges, function(x) try(x[x@data$SEASONAL%in%breeding.regions,]))
    spp<-unique(unlist(lapply(spp.ranges,function(x) as.character(x@data$SCINAME[1]))))
    spp.ranges.separated<-mclapply(mc.cores=ncores,spp, function(x){grep.spp<-grep(x,unlist(lapply(spp.ranges,function(x)as.character(x@data$SCINAME[1]))));spp.ranges[grep.spp]})
    
    saveRDS(spp.ranges.separated,file=paste('spp.ranges.',end,'.rds',sep=''))
    
    
  }
}
generate_birds_spp_errors.ranges<-function(ncores,timeout.seconds,spp.list){
  ##load spp ranges and process
  cat('opening layers','\n')
  spp.list<-spp.list
  for (i in seq(from=1, to=length(spp.list), by=500)){
    start<-i
    end<-i+499
    if (end > length(spp.list)){
      end<-length(spp.list)
    } 
    spp.list.small<-spp.list[start:end]
    ptm <- proc.time()
    message<-paste('reading from ',start,' to ',end,'\n',sep='')
    write(message,file=paste('./birds_control_errors_',timeout.seconds,'.txt',sep=''),append = T)
    time<-system.time(spp.ranges<-mclapply(mc.cores=ncores,spp.list.small, function(x) try(withTimeout(readOGR(sub(strsplit(x,'/')[[1]][1],pattern='^',replacement='./raw_data/birds/birds_ranges_2016/'),sub(strsplit(x,'/')[[1]][2],pattern='.shp',replacement='')),timeout=timeout.seconds,onTimeout='error'))))
    write(time,file=paste('./birds_control_errors_',timeout.seconds,'.txt',sep=''),append = T)
    #selecting breeding regions only
    breeding.regions<-c(1,2,5)
    #clean errors
    errors<-vector('numeric')
    for (f in 1:length(spp.ranges)){
      if (class(spp.ranges[[f]])!='SpatialPolygonsDataFrame'){
        errors<-c(errors,f)
      }
    }
    
    if(length(errors)>0){
      spp.ranges<-spp.ranges[-errors]
      cat(paste('some spp layers were not read: ',errors,'\n'))
      error.species<-spp[errors]
      write(error.species,file=paste('./birds_errors_',timeout.seconds,'.txt',sep=''),append = T)
      
    }
    spp.ranges<-mclapply(mc.cores=ncores,spp.ranges, function(x) try(x[x@data$SEASONAL%in%breeding.regions,]))
    spp<-unique(spp.ranges@data$SCINAME)
    spp.ranges.separated<-mclapply(mc.cores=ncores,spp, function(x) spp.ranges[spp.ranges@data$SCINAME==x,])
    
    
    saveRDS(spp.ranges.separated,file=paste('spp.ranges.',timeout.seconds,'.second.errors.',end,'.rds',sep=''))

    
  }
}
generate_birds_spp.ranges.separated.clean<-function(ncores,spp.ranges.separated){
  ##get species names
  spp<-unlist(lapply(spp.ranges.separated,function(x) as.character(x$SCINAME[1])))
  spp<-unlist(lapply(spp,function(x) sub(x,pattern=' ',replacement='_')))
  ##load spp ranges and process
  write('cleaning spp polygons',file='./birds_control_pipeline.txt',append=T)
  #version 0.2 of cleangeo does not work, use 0.1.1
  spp.ranges.separated<-mclapply(mc.cores=ncores,spp.ranges.separated, function(x) try(clgeo_Clean(x,print.log = TRUE)))
  #this is for version 0.2 of cleangeo but it takes too long
  #spp.ranges.separated<-mclapply(mc.cores=ncores,spp.ranges.separated, function(x) try(withTimeout(clgeo_Clean(x),timeout=6,onTimeout = 'error')))
  #merge spp polygons
  write ('merging spp polygons',file='./birds_control_pipeline.txt',append=T)
  spp.ranges.separated<-mclapply(mc.cores=ncores,spp.ranges.separated,function (x) try(gUnaryUnion(x)))
  #fix polygons
  write ('fixing spp polygons',file='./birds_control_pipeline.txt',append=T)
  spp.ranges.separated<-mclapply(mc.cores=ncores,spp.ranges.separated, function(x) try(gBuffer(x, width=0)))
  #check they're valid
  write ('validity of spp polygons',file='./birds_control_pipeline.txt',append=T)
  valid<-mclapply(mc.cores=ncores,spp.ranges.separated, function(x) try(gIsValid(x)))
  write ('number of non valid spp polygons',file='./birds_control_pipeline.txt',append=T)
  nonvalid<-length(valid[valid==F])
  write (nonvalid,file='./birds_control_pipeline.txt',append=T)
  write ('deleting non valid polygons',file='./birds_control_pipeline.txt',append=T)
  #delete non valid species
  if (length(valid[valid==F])>0){
    print (length(valid[valid==F]))
    not.valid<-vector('numeric')
    for (g in 1:length(valid)){
      if(valid[[g]]!=TRUE){
        print (g)
        not.valid<-c(g,not.valid)    
      }
    }
    write(spp[not.valid],file='./birds_errors_cleaned_separated.txt')
    spp<-spp[-not.valid]
    spp.ranges.separated<-spp.ranges.separated[-not.valid]
    write(spp,file='./birds_cleaned_separated.txt')
  }
  saveRDS(spp.ranges.separated,file='./output/birds/birds.spp.ranges.cleaned.separated.rds')
}

#################################################
######this loads the WWF realms & splits, cleans, merges, fixes and validates them
######and saves result into an realm.world.proj object
######it also creates a grid in each realm (size determined by grid.size) and saves them to output/birds
get_gridded_realms<-function(grid.size){
  #load WWF Biogeographic Realms 2004
  realms.layer<-readOGR('../raw_data/GIS/Generalised_Biogeographic_Realms_2004/','brpol_fin_dd')
  #get region names
  realms.names<-unique(realms.layer@data$REALM)
  realms.names<-as.character(realms.names[realms.names!='Antarctic'])
  realms<-unique(realms.layer@data$REALM)
  cat('getting realms ready','\n')
  #drop Antarctica
  realms<-realms[realms!='Antarctic']
  #separate by spp name
  realms.separated<-mclapply(mc.cores=ncores,realms, function(x) realms.layer[realms.layer@data$REALM==x,])
  #merge spp polygons
  realms.separated<-mclapply(mc.cores=ncores,realms.separated, gUnaryUnion)
  #fix polygons
  realms.separated<-mclapply(mc.cores=ncores,realms.separated, function(x) gBuffer(x, width=0))
  valid<-mclapply(mc.cores=ncores,realms.separated, gIsValid)
  cat('non valid realms','\n')
  length(valid[valid==F])
  #read the world
  world<-readOGR('../raw_data/GIS/realms/','newRealms')
  #merge the world
  world.merge<-gUnaryUnion(world)
  #fix the world
  cat ('fixing layer','\n')
  world.merge<-gBuffer(world.merge, width=0)
  #check validity of world
  valid<-gIsValid(world.merge)
  length(valid[valid==F])
  #intersect the world with realms
  cat('intersecting world with realms','\n')
  realm.world<-list()
  for (i in 1:length(realms.separated)){
    realm.world[[i]]<-try(gIntersection(realms.separated[[i]],world.merge,drop_lower_td=T))
    #plot(realm.world[[i]])
  }
  realm.world.proj<-lapply(realm.world,function(x) spTransform(x, CRS("+proj=cea +units=km")))
  #create grid (or load it)
  for (i in 1:length(realm.world.proj)){
    create_grid_layer(realm.world.proj[[i]],grid.size,realms.names[i])  
  }
  saveRDS(realm.world.proj,file='./output/realm.world.proj.rds')
}



####################################################
#########function to create a grid in a layer
##########layer is a realm/world layer, merged, fixed,validated and projected
##########cell size is the size of the cells (units depend on the projection of the layer)
###########name
create_grid_layer<-function(layer,cell.size,name){
  ##################build grid
  #get the limits of the layer
  bb<-bbox(layer)
  #define cell size (cs)
  cs<-c(cell.size,cell.size)
  cat('Cell size ',cs,'\n')
  #get the number of cells per direction
  cd <- ceiling(diff(t(bb))/cs)
  #get the number of cells in x and y axes
  x.width<-(bb['x','max']-bb['x','min'])/cd['max','x']
  y.width<-(bb['y','max']-bb['y','min'])/cd['max','y']
  #create an empty list for the grid
  grid<-list()
  #this string stores the number of cycles of the loop
  cont<-0
  #for each cell in x axis
  cat ('creating grid','\n')
  for (a in 1:(cd['max','x'])){
    #for each cell in y axis
    for (b in 1:(cd['max','y'])){
      cont<-cont+1
      #create the polygon
      grid[[cont]]<-Polygon(matrix(c(bb['x','min']+((a-1)*x.width),bb['x','min']+((a*x.width)),bb['x','min']+((a*x.width)),bb['x','min']+((a-1)*x.width),bb['x','min']+((a-1)*x.width),bb['y','max']-((b-1)*y.width),bb['y','max']-((b-1)*y.width),bb['y','max']-((b*y.width)),bb['y','max']-((b*y.width)),bb['y','max']-((b-1)*y.width)),ncol=2),hole=F)
      grid[[cont]]<-Polygons(list(grid[[cont]]),1)
      grid[[cont]]<-SpatialPolygons(list(grid[[cont]]))
      #adding correct projection
      proj4string(grid[[cont]])<-proj4string(layer)
      
    }
  }
  #save 
  grid<-lapply(grid,function(x){if(gIntersects(x,layer)){return(x)}})
  grid<-grid[!sapply(grid, is.null)]
  saveRDS(grid,file=paste('./output/birds/grid','_',name,'_',cell.size,'.rds',sep=''))
  return(grid)
}

####################################################################
######this generates the spp.ranges.separated.realm for a given realm
######it overlaps all spp.ranges.separated with a realm.world.proj
overlap_realm_species<-function(spp.ranges.separated,position,realm.names.vector,ncores){
  write (paste('analysing ',realm.names.vector[position],sep=''),file='./birds_control_pipeline.txt',append=T)
  cat('project spp','\n')
  spp.ranges.separated.proj<-mclapply(mc.cores=ncores,spp.ranges.separated, function(x) if (!is.null(x)) try(spTransform(x, CRS("+proj=cea +units=km"))))
  #this cleaning step solves problems with intersections later on, DO NOT REMOVE
  cat('cleaning spp ranges','\n')
  spp.ranges.separated.proj.clean<-mclapply(mc.cores=ncores,spp.ranges.separated.proj, function(x) try(clgeo_Clean(x, print.log = TRUE)))
  cat('intersecting spp with geography','\n')
  spp.ranges.separated.realm<-mclapply(mc.cores=ncores,spp.ranges.separated.proj.clean, function(x) try(gIntersection(x,realm.world.proj[[position]],drop_lower_td=T)))
  cat('still intersecting spp with geography','\n')
  spp.ranges.separated.realm.no.drop<-mclapply(mc.cores=ncores,spp.ranges.separated.proj.clean, function(x) try(gIntersection(x,realm.world.proj[[position]])))
  a<-vector('numeric')
  for (i in 1:length(spp.ranges.separated.realm)){
    if(identical(spp.ranges.separated.realm[[i]],spp.ranges.separated.realm.no.drop[[i]]) ==FALSE){
      a<-c(a,i)
    }
  }
  if(length(a)>0){
    for (b in 1:length(a)){
      spp.ranges.separated.realm[[a[b]]]<-spp.ranges.separated.realm.no.drop[[a[b]]]@polyobj
    }
    
  }
  #Fix some objects (have lines apart from polygons; they appear as SpatialCollections rgeos objects instead of SpatialPolygons sp objects)
  cat ('fixing intersections','\n')
  spp.ranges.separated.realm<-mclapply(mc.cores=ncores,spp.ranges.separated.realm, function(x) if (!is.null(x) && class(x)!="SpatialPolygons") try(gPolygonize(x)) else x)
  saveRDS(spp.ranges.separated.realm,file=paste('./output/',realm.names.vector[position],'.birds.spp.ranges.separated.realm.rds',sep=''))
}

#####################################################
#function to overlap a grid with species intersected with a layer
#species.object: spp object intersected with layer
#grid: grid object
#path: is the path that will be used to save the output
#cell.size: is the cell size of the grid, will be used just for naming the output
#name: name to append to the output (i.e,name of the realm)
grid_overlap_species<-function(grid,species.object,species.list,path,cell.size,name){
  write (paste('intersecting grid with spp in realm','\n',sep=''),file=paste('./output/birds_intersecting_control',name,'.txt',sep=''),append=T)
  spp.ranges.grid<-list()
  #  spp.ranges.grid.no.drop<-list()
  for (i in 1:length(species.object)){
    write(paste('analysing species ',i,'\n',sep=''),file=paste('./output/birds_intersecting_control',name,'.txt',sep=''),append=T)
    if(is.null(species.object[[i]])){
      spp.ranges.grid[[i]]<-list(NULL)
      next;
    } 
    else{
      spp.ranges.grid[[i]]<-mclapply(mc.cores=ncores,grid,function(x) if (!is.null(x)){if(gIntersects(x,species.object[[i]])){try(gIntersection(x,species.object[[i]],drop_lower_td = T))}else{return(NULL)}} else{x})
      
    }
    
    #THIS MIGHT BE IMPORTANT FOR BIRDS
    # if (!is.null(species.object[[i]])){
    #    spp.ranges.grid[[i]]<-mclapply(mc.cores=ncores,grid,function(x) try(gIntersection(x,species.object[[i]],drop_lower_td = T)))
    #    numbers.ok<-c(numbers.ok,i)
    #spp.ranges.grid.no.drop[[i]]<-mclapply(mc.cores=ncores,grid,function(x) try(gIntersection(x,species.object[[i]])))
    #if(identical(spp.ranges.grid[[i]],spp.ranges.grid.no.drop[[i]]) ==FALSE){
    #cat(i,'\n')
    #  spp.ranges.grid[[i]]<-spp.ranges.grid.no.drop[[i]]@polyobj
    #}
    #  }  
    #  else{
    #    spp.ranges.grid[[i]]<-NULL
    #   numbers.null<-c(numbers.null,i)
    #  }
  }
  spp.ranges.grid.count<-vector("list", length(spp.ranges.grid))
  #for each species get the numbers of the grid where they occur
  for (i in 1:length(spp.ranges.grid)){
    if(length(spp.ranges.grid[[i]])==1){
      spp.ranges.grid.count[[i]]<-list()
      spp.ranges.grid.count[[i]]<-list(NULL)
      next;
    }
    else if (!is.null(spp.ranges.grid[[i]])){
      spp.ranges.grid.count[[i]]<-list()
      spp.ranges.grid.count[[i]]<-list(which(unlist(lapply(spp.ranges.grid[[i]], function (x) !is.null(x)))))
      next;
    }
    
  }
  write (paste('analysing spp occurrence in grid','\n',sep=''),file=paste('./output/birds_intersecting_control',name,'.txt',sep=''),append=T)
  #for each grid, get how many species hit
  grid.occurrences<-data.frame(table(unlist(spp.ranges.grid.count)))
  colnames(grid.occurrences)<-c('grid','number.of.species')
  grid.occurrences$grid<-as.numeric(as.character(grid.occurrences$grid))
  #fill the grid table with grids with no hit
  grid.occurrences<-merge(data.frame(1:length(grid)),grid.occurrences,by.x="X1.length.grid.",by.y='grid',all.x=TRUE)
  colnames(grid.occurrences)<-c('cell','number.of.species')
  #for each species, get the cells where it occurs and store it in species.grid
  species.grid<-data.frame()
  for (i in 1:length(spp.ranges.grid.count)){
    species.grid[i,1]<-species.list[i]
    species.grid[i,2]<-paste(unlist(spp.ranges.grid.count[[i]]),collapse=' ')
  }
  colnames(species.grid)<-c('spp','cells')
  species.grid$spp<-gsub(species.grid$spp,pattern=' ',replacement='_')
  species.grid$range.cells<-apply(species.grid,1,function(x) length(strsplit(as.character(x['cells']),' ')[1][[1]]))
  #write outputs
  write.table(grid.occurrences,file=paste(path,cell.size,'_',name,'_richness_grid_table.txt',sep=''),quote=F,row.names=F,sep='\t')
  write.table(species.grid,file=paste(path,cell.size,'_',name,'_species_gridoccurrence_table.txt',sep=''),quote=F,row.names=F,sep='\t')
}


args<-commandArgs(trailingOnly = TRUE)
cat(args,'\n',file='./output/birds_intersecting_control.txt')
ncores<-args[1]

generate_birds_spp.ranges(ncores,timeout.seconds=200)
errors.spp.list<-read.table('./birds_errors_200.txt',stringsAsFactors = F)
errors.spp.list<-errors.spp.list$V1
errors.spp.list<-unname(sapply(errors.spp.list,function(x)paste(x,'.shp',sep='')))
generate_birds_spp_errors.ranges(ncores=ncores,timeout.seconds=10000,spp.list=errors.spp.list)

###concatenate all spp.ranges.rds
list.rds<-list.files(path='./',pattern='spp.ranges..+.rds')
list.rds.objects<-lapply(list.rds,function(x) readRDS(x))
list.rds.objects<-unlist(list.rds.objects)
species.names<-lapply(list.rds.objects,function(x) as.character(x@data$SCINAME))
spp<-unique(unlist(species.names))
spp.ranges.separated<-mclapply(mc.cores=ncores,spp, function(x) {grep.spp<-grep(paste('^',x,'$',sep=''),species.names);if(length(grep.spp)==1){x<-list.rds.objects[[grep.spp]]}else{x<-do.call('rbind',list.rds.objects[grep.spp])};return(x)})
saveRDS(spp.ranges.separated,file='./birds.spp.ranges.separated.rds')

spp.ranges.separated<-readRDS('./birds.spp.ranges.separated.rds')
generate_birds_spp.ranges.separated.clean(ncores=ncores,spp.ranges.separated=spp.ranges.separated)

####this generates a realm.world.proj object and grids of 100x100 in each realm
get_gridded_realms(100)

####we have to load the realm.world.proj object now for the rest to work
realm.world.proj<-readRDS('./output/realm.world.proj.rds')

#####we get the names of the realms here
#load WWF Biogeographic Realms 2004
realms.layer<-readOGR('../raw_data/GIS/Generalised_Biogeographic_Realms_2004/','brpol_fin_dd')
#get region names
realms.names<-unique(realms.layer@data$REALM)
realms.names<-as.character(realms.names[realms.names!='Antarctic'])

spp.ranges.separated<-readRDS('./output/birds/birds.spp.ranges.cleaned.separated.rds')

#this generates a "REALMNAME.spp.ranges.separated.realm.Rsave that contains a spp.ranges.separated.realm object
#this spp.ranges.separated.realm is to be used with the corresponding grid to intersect and find the species ocurring in each cell in the grid
overlap_realm_species(spp.ranges.separated,position=1,realms.names,ncores=ncores)
overlap_realm_species(spp.ranges.separated,position=2,realms.names,ncores=ncores)
overlap_realm_species(spp.ranges.separated,position=3,realms.names,ncores=ncores)
overlap_realm_species(spp.ranges.separated,position=4,realms.names,ncores=ncores)
overlap_realm_species(spp.ranges.separated,position=5,realms.names,ncores=ncores)
overlap_realm_species(spp.ranges.separated,position=6,realms.names,ncores=ncores)
overlap_realm_species(spp.ranges.separated,position=7,realms.names,ncores=ncores)

####WARNING: this commented bit below should be split by realms or combining a couple of realms (e.g, run Australasian + Neotropical)#####
###to speed things up a bit on hydrogen (using ~24 cores works fine)
###important: Palearctic and Nearctic take the longest, so do not run them in the same analysis!

#####get spp names
####spp<-read.table('./birds_cleaned_separated.txt',stringsAsFactors=F)
####spp<-spp$V1
####grid<-readRDS('./output/birds/grid_Australasian_100.rds')
####spp.ranges.separated.realm<-readRDS('./output/Australasian.birds.spp.ranges.separated.realm.rds')
####grid_overlap_species(grid,spp.ranges.separated.realm,spp,'./output/birds/',100,'Australasian_realms')

####
#####get spp names
####spp<-read.table('./birds_cleaned_separated.txt',stringsAsFactors=F)
####spp<-spp$V1
####grid<-readRDS('./output/birds/grid_Palearctic_100.rds')
####spp.ranges.separated.realm<-readRDS('./output/Palearctic.birds.spp.ranges.separated.realm.rds')
####grid_overlap_species(grid,spp.ranges.separated.realm,spp,'./output/birds/',100,'Palearctic_realms')
####

#####get spp names
####spp<-read.table('./birds_cleaned_separated.txt',stringsAsFactors=F)
####spp<-spp$V1
####grid<-readRDS('./output/birds/grid_Nearctic_100.rds')
####spp.ranges.separated.realm<-readRDS('./output/Nearctic.birds.spp.ranges.separated.realm.rds')
####grid_overlap_species(grid,spp.ranges.separated.realm,spp,'./output/birds/',100,'Nearctic_realms_new')

#####get spp names
####spp<-read.table('./birds_cleaned_separated.txt',stringsAsFactors=F)
####spp<-spp$V1
####grid<-readRDS('./output/birds/grid_Oceanic_100.rds')
####spp.ranges.separated.realm<-readRDS('./output/Oceanic.birds.spp.ranges.separated.realm.rds')
####grid_overlap_species(grid,spp.ranges.separated.realm,spp,'./output/birds/',100,'Oceanic_realms')

#####get spp names
####spp<-read.table('./birds_cleaned_separated.txt',stringsAsFactors=F)
####spp<-spp$V1
####grid<-readRDS('./output/birds/grid_Indo-Malay_100.rds')
####spp.ranges.separated.realm<-readRDS('./output/Indo-Malay.birds.spp.ranges.separated.realm.rds')
####grid_overlap_species(grid,spp.ranges.separated.realm,spp,'./output/birds/',100,'Indo-Malay_realms')
####

#####get spp names
####spp<-read.table('./birds_cleaned_separated.txt',stringsAsFactors=F)
####spp<-spp$V1
####grid<-readRDS('./output/birds/grid_Afrotropical_100.rds')
####spp.ranges.separated.realm<-readRDS('./output/Afrotropical.birds.spp.ranges.separated.realm.rds')
####grid_overlap_species(grid,spp.ranges.separated.realm,spp,'./output/birds/',100,'Afrotropical_realms')
####

#####get spp names
####spp<-read.table('./birds_cleaned_separated.txt',stringsAsFactors=F)
####spp<-spp$V1
####grid<-readRDS('./output/birds/grid_Neotropical_100.rds')
####spp.ranges.separated.realm<-readRDS('./output/Neotropical.birds.spp.ranges.separated.realm.rds')
####grid_overlap_species(grid,spp.ranges.separated.realm,spp,'./output/birds/',100,'Neotropical_realms')


