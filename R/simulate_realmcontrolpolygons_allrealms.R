library(raster)
library(spatstat)
library(colorspace)
library(spdep)
library(gridExtra)
library(sp)
library(fields)
########################################################################################################################
########################################################################################################################
#ppp is a ppp object resulting from simulate.kppm
#grid coordinats is a df with x and y coordinates
get_closest_cells_from_sim<-function(ppp,grid.coordinates){
    #store coordinates of sim
    points.df<-as.data.frame(cbind(ppp$x,ppp$y),stringsAsFactors=F)
    #distance between all points in sim vs all points in grid (grid.coordinates)
    rdist<-rdist(points.df,grid.coordinates)
    #get the minimum distance by rows
    rdist.df<-as.data.frame(rdist,stringsAsFactors=F)
    #have to get the check in a loop because of multiple hits to the same cell
    list.of.cells<-numeric()
    cat('fitting points to the grid','\n')
    for (i in 1:nrow(rdist.df)){
      row<-numeric()
      hit<-numeric()
      row<-rdist.df[i,]
      hit<-which.min(row)
      tried.hits<-numeric()
      while(hit%in%list.of.cells){
          #cat(i, ' checkin','\n')
          row<-row[-hit]
          hit<-which.min(row)
          tried.hits<-c(tried.hits,hit)
          if (length(tried.hits)+1==ncol(rdist.df)){
            break
          }
      }
      list.of.cells[i]<-hit
    }
    return(unique(list.of.cells))
}


#table.realm is a table with $cells
#hotspot.cells is a vector with the hotspot cells to create the control for
#grid.file is the path to the grid.world
#nreplicates is the number of control sets to generate
simulate_control_regions_realm<-function(table.realm,hotspot.cells,grid.file,nreplicates){
  grid.world<-readRDS(grid.file)
  #get coordinates of cells in a realm
  grid.coordinates<-lapply(grid.world[table.realm$cells],function(x) coordinates(x))
  grid.coordinates<-do.call("rbind", grid.coordinates)
  grid.coordinates<-as.data.frame(grid.coordinates,stringsAsFactors = F)
  #add a dummy z column for the raster
  grid.coordinates[,3]<-1
  colnames(grid.coordinates)<-c('x','y','z')
  rownames(grid.coordinates)<-c(1:nrow(grid.coordinates))
  #convert to SpatialPointsDataFrame
  #help here http://chris35wills.github.io/gridding_data/
  #and more help here https://gis.stackexchange.com/questions/187798/create-polygons-of-the-extents-of-a-given-raster-in-r
  coordinates(grid.coordinates) = ~x+y 									   
  #rasterize
  rast <- raster(ext=extent(grid.coordinates), resolution=100)
  # rasterize your irregular points 
  rasOut<-rasterize(grid.coordinates, rast, grid.coordinates$z, fun = mean) # we use a mean function here to regularly grid the irregular input points
  #plot the raster to check
  #plot(rasOut)
  #To get the rectangular extent
  e <- extent(rasOut)
  p <- as(e, 'SpatialPolygons') 
  #To get a polygon that surrounds cells that are not NA
  # make all values the same. Either do
  r <- rasOut > -Inf
  # or alternatively
  # r <- reclassify(x, cbind(-Inf, Inf, 1))
  # convert to polygons (you need to have package 'rgeos' installed for this to work)
  pp <- rasterToPolygons(r, dissolve=TRUE)
  #pp is a SpatialPolygons dataframe with the polygons we want to use to build the window 
  #extract coordinates of all Polygon objects, put them in a list
  list.polygon.coordinates<-list()
  for (i in 1:length(pp@polygons[[1]]@Polygons)){
    list.polygon.coordinates[[i]]<-pp@polygons[[1]]@Polygons[[i]]@coords
    
  }
  #reverse all (owin needs coordinates anticlockwise)
  list.polygon.coordinates<-lapply(list.polygon.coordinates,function(x)cbind(rev(x[,'x']),rev(x[,'y'])))
  owin.realm<-owin(poly=list.polygon.coordinates)
  #plot to check it's ok
  #plot(owin.realm)
  #get the coordinates of the hotspot cells
  hotspot.coordinates<-lapply(grid.world[hotspot.cells],function(x) coordinates(x))
  hotspot.coordinates<-do.call("rbind", hotspot.coordinates)
  hotspot.coordinates<-as.data.frame(hotspot.coordinates,stringsAsFactors = F)
  #get the hotspot points in
  spat.hotspot<-as.ppp(hotspot.coordinates,owin.realm)
  #plot(spat.hotspot,pch=16,cex=.5)
  #fit kppm to fit cluster process to the hotspot point patter
  fitted<-kppm(spat.hotspot)
  #Generates simulated realisations from a fitted cluster point process model.
  kfitted<-simulate.kppm(fitted,nsim=nreplicates)
  #plot(kfitted,pch=16,cex=.5)
  #plot(grid.coordinates)
  #points(kfitted[[1]]$x,kfitted[[1]]$y,col='red',pch=16)
  grid.coordinates<-lapply(grid.world[table.realm$cells],function(x) coordinates(x))
  grid.coordinates<-do.call("rbind", grid.coordinates)
  grid.coordinates<-as.data.frame(grid.coordinates,stringsAsFactors = F)
  colnames(grid.coordinates)<-c('x','y')
  #this controls for polygons that have more points than cells in the grid and deletes them
  #can use this filter later if needed to constrain size of simulated polygons
  kfitted.size<-unlist(lapply(kfitted,function(x)length(x$x)))
  kfitted.too.big<-which(kfitted.size>=nrow(grid.coordinates))
  while(length(kfitted.too.big)>0){
    kfitted<-kfitted[-kfitted.too.big]
    #add new sims to replace the big ones
    kfitted.new<-simulate.kppm(fitted,nsim=nreplicates)
    kfitted.new.size<-unlist(lapply(kfitted.new,function(x)length(x$x)))
    kfitted.new.ok<-which(kfitted.new.size<nrow(grid.coordinates))
    kfitted.new<-kfitted.new[kfitted.new.ok[1:length(kfitted.too.big)]]
    kfitted<-c(kfitted,kfitted.new)
    kfitted.size<-unlist(lapply(kfitted,function(x)length(x$x)))
    kfitted.too.big<-which(kfitted.size>=nrow(grid.coordinates))
  }
  kfitted.size<-unlist(lapply(kfitted,function(x)length(x$x)))
  #here I only select polygons that are at least half of the size of the real polygons
  kfitted.too.small<-which(kfitted.size<(nrow(grid.coordinates)/2))
  while(length(kfitted.too.small)>0){
    kfitted<-kfitted[-kfitted.too.small]
    #add new sims to replace the big ones
    kfitted.new<-simulate.kppm(fitted,nsim=nreplicates)
    kfitted.new.size<-unlist(lapply(kfitted.new,function(x)length(x$x)))
    kfitted.new.ok<-which(kfitted.new.size>0)
    kfitted.new<-kfitted.new[kfitted.new.ok[1:length(kfitted.too.small)]]
    kfitted<-c(kfitted,kfitted.new)
    kfitted.size<-unlist(lapply(kfitted,function(x)length(x$x)))
    kfitted.too.big<-which(kfitted.size>=nrow(grid.coordinates))
    kfitted.too.small<-which(kfitted.size==0)
    
  }
  points<-lapply(kfitted,function(x)get_closest_cells_from_sim(ppp=x,grid.coordinates = grid.coordinates))
  return(points)
}

simulate_realms_control_polygons<-function(hotspots.tablefile,grid.file,path,nreplicates){
  table<-read.table(hotspots.tablefile,header=T,sep='\t',stringsAsFactors = F)
  table$hotspot<-as.factor(table$hotspot)
  #######divide dataset by realms
  #load realm grids
  grids<-list.files(path='./output/grids/',pattern='grid_.*_100.rds')
  grid.world<-grids[grep('World_RealmsMerged',grids)]
  if(length(grid.world)>0){
    grid.world<-grids[grep('World_RealmsMerged',grids)]
    grid.realms.names<-grids[-grep('World_RealmsMerged',grids)]
  }
  grid.world<-readRDS(paste('./output/grids/',grid.world,sep=''))
  grid.realms<-lapply(grid.realms.names,function(x) readRDS(paste('./output/grids/',x,sep='')))
  
  #get the limits of each realm in the global grid
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
  #grid.cells.df has the n of cells of each realm and the starting and ending cell
  
  #subset it to get the table.realm for each realm
  table.afro<-table[table$cells%in%c(grid.cells.df[1,'start.cell']:grid.cells.df[1,'end.cell']),]
  table.austral<-table[table$cells%in%c(grid.cells.df[2,'start.cell']:grid.cells.df[2,'end.cell']),]
  table.indo<-table[table$cells%in%c(grid.cells.df[3,'start.cell']:grid.cells.df[3,'end.cell']),]
  table.nearctic<-table[table$cells%in%c(grid.cells.df[4,'start.cell']:grid.cells.df[4,'end.cell']),]
  table.neotropical<-table[table$cells%in%c(grid.cells.df[5,'start.cell']:grid.cells.df[5,'end.cell']),]
  table.palearctic<-table[table$cells%in%c(grid.cells.df[6,'start.cell']:grid.cells.df[6,'end.cell']),]
  
  
  afrotropical.hotspot<-table.afro[table.afro$hotspot==1,]$cells
  austral.hotspot<-table.austral[table.austral$hotspot==1,]$cells
  indo.hotspot<-table.indo[table.indo$hotspot==1,]$cells
  nearctic.hotspot<-table.nearctic[table.nearctic$hotspot==1,]$cells
  neotropical.hotspot<-table.neotropical[table.neotropical$hotspot==1,]$cells
  palearctic.hotspot<-table.palearctic[table.palearctic$hotspot==1,]$cells
  
  #simulating nreplicates control polygons for realm with a lower size limit
  cat(paste('simulating afrotrop control - ',nreplicates,sep=''),'\n')
  points.afro.replicates<-simulate_control_regions_realm(table.realm=table.afro,hotspot.cells = afrotropical.hotspot,grid.file=grid.file,nreplicates=nreplicates)
  saveRDS(points.afro.replicates,file=paste(path,'points.afrotropical.',nreplicates,'.RDS',sep=''))
  cat(paste('simulating austral control - ',nreplicates,sep=''),'\n')
  points.austral.replicates<-simulate_control_regions_realm(table.realm=table.austral,hotspot.cells = austral.hotspot,grid.file=grid.file,nreplicates=nreplicates)
  saveRDS(points.austral.replicates,file=paste(path,'points.austral.',nreplicates,'.RDS',sep=''))
  cat(paste('simulating indo control - ',nreplicates,sep=''),'\n')
  points.indo.replicates<-simulate_control_regions_realm(table.realm=table.indo,hotspot.cells = indo.hotspot,grid.file=grid.file,nreplicates=nreplicates)
  saveRDS(points.indo.replicates,file=paste(path,'points.indo.',nreplicates,'.RDS',sep=''))
  cat(paste('simulating nearctic control - ',nreplicates,sep=''),'\n')
  points.nearctic.replicates<-simulate_control_regions_realm(table.realm=table.nearctic,hotspot.cells = nearctic.hotspot,grid.file=grid.file,nreplicates=nreplicates)
  saveRDS(points.nearctic.replicates,file=paste(path,'points.nearctic.',nreplicates,'.RDS',sep=''))
  cat(paste('simulating neotropical control - ',nreplicates,sep=''),'\n')
  points.neotropical.replicates<-simulate_control_regions_realm(table.realm=table.neotropical,hotspot.cells = neotropical.hotspot,grid.file=grid.file,nreplicates=nreplicates)
  saveRDS(points.neotropical.replicates,file=paste(path,'points.neotropical.',nreplicates,'.RDS',sep=''))
  cat(paste('simulating palearctic control - ',nreplicates,sep=''),'\n')
  points.palearctic.replicates<-simulate_control_regions_realm(table.realm=table.palearctic,hotspot.cells = palearctic.hotspot,grid.file=grid.file,nreplicates=nreplicates)
  saveRDS(points.palearctic.replicates,file=paste(path,'points.palearctic.',nreplicates,'.RDS',sep=''))
  
}






#########plotting
##par(mfrow=c(4,5))
##grid.coordinates<-lapply(grid.world[table.afro$cells],function(x) coordinates(x))
##grid.coordinates<-do.call("rbind", grid.coordinates)
##grid.coordinates<-as.data.frame(grid.coordinates,stringsAsFactors = F)
##plot(grid.coordinates[,1],grid.coordinates[,2],pch=16,cex=.5);points(grid.coordinates[afrotropical.hotspot,1],grid.coordinates[afrotropical.hotspot,2],pch=16,col='red',cex=.5)
##lapply(points.afro.50,function(x){plot(grid.coordinates[,1],grid.coordinates[,2],pch=16,cex=.5);points(grid.coordinates[x,1],grid.coordinates[x,2],pch=16,col='red',cex=.5)})
##dev.off()
##
##par(mfrow=c(4,5))
##grid.coordinates<-lapply(grid.world[table.austral$cells],function(x) coordinates(x))
##grid.coordinates<-do.call("rbind", grid.coordinates)
##grid.coordinates<-as.data.frame(grid.coordinates,stringsAsFactors = F)
##austral.hotspot.ref<-austral.hotspot-grid.cells.df[grid.cells.df$realms.names=='Australasian',]$start.cell
##plot(grid.coordinates[,1],grid.coordinates[,2],pch=16,cex=.5);points(grid.coordinates[austral.hotspot.ref,1],grid.coordinates[austral.hotspot.ref,2],pch=16,col='red',cex=.5)
##lapply(points.austral.50,function(x){plot(grid.coordinates[,1],grid.coordinates[,2],pch=16,cex=.5);points(grid.coordinates[x,1],grid.coordinates[x,2],pch=16,col='red',cex=.5)})
##dev.off()
##
##par(mfrow=c(4,5))
##grid.coordinates<-lapply(grid.world[table.indo$cells],function(x) coordinates(x))
##grid.coordinates<-do.call("rbind", grid.coordinates)
##grid.coordinates<-as.data.frame(grid.coordinates,stringsAsFactors = F)
##indo.hotspot.ref<-indo.hotspot-grid.cells.df[grid.cells.df$realms.names=='Indo-Malay',]$start.cell
##plot(grid.coordinates[,1],grid.coordinates[,2],pch=16,cex=.5);points(grid.coordinates[indo.hotspot.ref,1],grid.coordinates[indo.hotspot.ref,2],pch=16,col='red',cex=.5)
##lapply(points.indo.50,function(x){plot(grid.coordinates[,1],grid.coordinates[,2],pch=16,cex=.5);points(grid.coordinates[x,1],grid.coordinates[x,2],pch=16,col='red',cex=.5)})
##dev.off()
##
##par(mfrow=c(4,5))
##grid.coordinates<-lapply(grid.world[table.nearctic$cells],function(x) coordinates(x))
##grid.coordinates<-do.call("rbind", grid.coordinates)
##grid.coordinates<-as.data.frame(grid.coordinates,stringsAsFactors = F)
##nearctic.hotspot.ref<-nearctic.hotspot-grid.cells.df[grid.cells.df$realms.names=='Nearctic',]$start.cell
##plot(grid.coordinates[,1],grid.coordinates[,2],pch=16,cex=.5);points(grid.coordinates[nearctic.hotspot.ref,1],grid.coordinates[nearctic.hotspot.ref,2],pch=16,col='red',cex=.5)
##lapply(points.nearctic.50,function(x){plot(grid.coordinates[,1],grid.coordinates[,2],pch=16,cex=.5);points(grid.coordinates[x,1],grid.coordinates[x,2],pch=16,col='red',cex=.5)})
##dev.off()
##
##par(mfrow=c(4,5))
##grid.coordinates<-lapply(grid.world[table.neotropical$cells],function(x) coordinates(x))
##grid.coordinates<-do.call("rbind", grid.coordinates)
##grid.coordinates<-as.data.frame(grid.coordinates,stringsAsFactors = F)
##neotropical.hotspot.ref<-neotropical.hotspot-grid.cells.df[grid.cells.df$realms.names=='Neotropical',]$start.cell
##plot(grid.coordinates[,1],grid.coordinates[,2],pch=16,cex=.5);points(grid.coordinates[neotropical.hotspot.ref,1],grid.coordinates[neotropical.hotspot.ref,2],pch=16,col='red',cex=.5)
##lapply(points.neotropical.50,function(x){plot(grid.coordinates[,1],grid.coordinates[,2],pch=16,cex=.5);points(grid.coordinates[x,1],grid.coordinates[x,2],pch=16,col='red',cex=.5)})
##dev.off()
##
##par(mfrow=c(4,5))
##grid.coordinates<-lapply(grid.world[table.palearctic$cells],function(x) coordinates(x))
##grid.coordinates<-do.call("rbind", grid.coordinates)
##grid.coordinates<-as.data.frame(grid.coordinates,stringsAsFactors = F)
##palearctic.hotspot.ref<-palearctic.hotspot-grid.cells.df[grid.cells.df$realms.names=='Palearctic',]$start.cell
##plot(grid.coordinates[,1],grid.coordinates[,2],pch=16,cex=.5);points(grid.coordinates[palearctic.hotspot.ref,1],grid.coordinates[palearctic.hotspot.ref,2],pch=16,col='red',cex=.5)
##lapply(points.palearctic.50,function(x){plot(grid.coordinates[,1],grid.coordinates[,2],pch=16,cex=.5);points(grid.coordinates[x,1],grid.coordinates[x,2],pch=16,col='red',cex=.5)})
##dev.off()
##
##############histogram of sizes of hotspots vs control regions in each realm
##points.afro.50<-readRDS(file='./output/mammals/new_tables_hotspots/immigration/controls/points.afro.50.RDS')
##points.austral.50<-readRDS(file='./output/mammals/new_tables_hotspots/immigration/controls/points.austral.50.RDS')
##points.indo.50<-readRDS(file='./output/mammals/new_tables_hotspots/immigration/controls/points.indo.50.RDS')
##points.nearctic.50<-readRDS(file='./output/mammals/new_tables_hotspots/immigration/controls/points.nearctic.50.RDS')
##points.neotropical.50<-readRDS(file='./output/mammals/new_tables_hotspots/immigration/controls/points.neotropical.50.RDS')
##points.palearctic.50<-readRDS(file='./output/mammals/new_tables_hotspots/immigration/controls/points.palearctic.50.RDS')
##
##
##pdf('./size_comparison_realmhotspots_vs_control.pdf',paper='a4r')
##par(mfrow=c(2,3))
##plot(density(unlist(lapply(points.afro.50,function(x)length(x)))),yaxt='n',ylab='',main='afrotropical control polygons size',xlab='number of grids')
##abline(v=length(afrotropical.hotspot),col='red')
##
##plot(density(unlist(lapply(points.austral.50,function(x)length(x)))),yaxt='n',ylab='',main='austral control polygons size',xlab='number of grids')
##abline(v=length(austral.hotspot),col='red')
##
##plot(density(unlist(lapply(points.indo.50,function(x)length(x)))),yaxt='n',ylab='',main='indo control polygons size',xlab='number of grids')
##abline(v=length(indo.hotspot),col='red')
##
##plot(density(unlist(lapply(points.nearctic.50,function(x)length(x)))),yaxt='n',ylab='',main='nearctic control polygons size',xlab='number of grids')
##abline(v=length(nearctic.hotspot),col='red')
##
##plot(density(unlist(lapply(points.neotropical.50,function(x)length(x)))),yaxt='n',ylab='',main='neotropical control polygons size',xlab='number of grids')
##abline(v=length(neotropical.hotspot),col='red')
##
##plot(density(unlist(lapply(points.palearctic.50,function(x)length(x)))),yaxt='n',ylab='',main='palearctic control polygons size',xlab='number of grids')
##abline(v=length(palearctic.hotspot),col='red')
##dev.off()
############they are all right in the middle of the distribution (as expected)
######