#run model in hydrogen cluster
.libPaths('/home/ecosys/ji247/R/x86_64-pc-linux-gnu-library/3.3')
setwd('/home/ecosys/ji247/hotspots_vertebrates/')

library(ape)
library(plyr)
library(phangorn)
library(geiger)
library(parallel)
library(diversitree)
###############################################
###########################


#world.table.file is 100_all_realms_ranges_plus_hotspots.txt
#treefile is the tree
#name has to be one of (afrotrop,austral,indo,nearctic,neotrop,palearctic)
#overlap is 0.8
prepare_realm_input_7areas<-function(world.table.file,treefile,name,overlap){
  tree<-read.tree(treefile)
  world.table<-read.table(world.table.file,header=T,stringsAsFactors = F,sep='\t')
  #collate tree and table
  tree<-drop.tip(tree,setdiff(tree$tip.label,world.table$spp))
  world.table<-world.table[world.table$spp%in%tree$tip.label,]
  #keep columns with proportion of range inside realm and proportion of range within hotspot
  drop.columns<-grep('cells',colnames(world.table))
  world.table<-world.table[,-drop.columns]
  colnames(world.table)<-sub(colnames(world.table),pattern='_WWF',replacement='')
  prepare_realm_input_7areas_plus_inhotoutrealm(table=world.table,name=name,overlap=0.80,tree=tree)
}

prepare_realm_input_7areas_plus_inhotoutrealm<-function(table,name,overlap,tree,replicate){
  table<-table[,c('spp',colnames(table)[grep('_realm.area',colnames(table))],paste(name,'_hotspot.area',sep=''))]
  table[is.na(table)]<-0
  #get a table without the hotspot area first
  realm.table<-table[,-grep('hotspot',colnames(table))]
  #process realm.table
  realm.table$model.state<-0
  #if overlap in any column is >0.8 assign to that region
  realm.exclusive.species<-apply(realm.table,1,function(x) unname(which(x[c(2:7)]>=overlap)))
  realm.exclusive.species<-lapply(realm.exclusive.species,function(x)if(length(x)==0){x<-0}else{return(x)})
  realm.table$model.state<-unlist(realm.exclusive.species)
  #sort out the multirealm species
  multirrealm.species<-realm.table[realm.table$model.state==0,]
  multirrealm.species.regions<-apply(multirrealm.species,1,function(x)unname(which(x[c(2:7)]>(1-overlap))))
  multirrealm.species.regions<-unlist(lapply(multirrealm.species.regions,function(x) if(length(x)==2){x<-paste(x,collapse=',')}else{x<-0}))
  realm.table[realm.table$model.state==0,'model.state']<-multirrealm.species.regions
  #drop species with 0 = species occurring in Oceanic realm
  realm.table<-realm.table[!realm.table$model.state==0,]
  #select column of realm to be analysed
  #e.g: 1 for afrotrop, 2 for austral,3 for indo, etc
  realm.character<-as.character(grep(name,colnames(realm.table)[-c(1,8)]))
  #create df with realm and hotspot overlaps
  realm.hotspot.table<-table[,grep(name,colnames(table))]
  #select species that with overlap > 0.8 with realm
  realm.hotspot.table<-realm.hotspot.table[which(realm.hotspot.table[,1]>=overlap),]
  realm.hotspot.table$model.state<-0
  #species occurring in hotspot and realm
  realm.hotspot.table[realm.hotspot.table[,2]<(1-overlap),'model.state']<-realm.character
  #species occurring in hotspot
  realm.hotspot.table[realm.hotspot.table[,2]>=overlap,'model.state']<-7
  #species occurring in realm
  realm.hotspot.table[realm.hotspot.table[,2]>(1-overlap)&realm.hotspot.table[,2]<overlap,'model.state']<-paste(realm.character,7,sep=',')
  realm.table[realm.table$model.state==realm.character,'model.state']<-realm.hotspot.table$model.state
  #select species occurring in target realm plus another realm
  target.states<-names(table(realm.table$model.state))[grep(realm.character,names(table(realm.table$model.state)))]
  if(length(grep('7',target.states))>0){
    target.states<-target.states[-grep('7',target.states)]
  }
  target.states<-target.states[-match(realm.character,target.states)]
  if(!is.na(target.states)){
    realm.plus.outside.table<-realm.table[realm.table$model.state%in%target.states,]
    #add hotspot overlap within realm to realm.plus.outside table
    realm.plus.outside.table.hotspots<-as.data.frame(cbind(table[table$spp%in%realm.plus.outside.table$spp,'spp'],table[table$spp%in%realm.plus.outside.table$spp,paste(name,'_hotspot.area',sep='')]),stringsAsFactors = F)
    colnames(realm.plus.outside.table.hotspots)<-c('spp',paste(name,'_hotspot.area',sep=''))
    realm.plus.outside.table.hotspots$model.state<-realm.table[realm.table$model.state%in%target.states,]$model.state
    realm.plus.outside.table.hotspots[realm.plus.outside.table.hotspots[,2]>=overlap,]$model.state<-sub(realm.plus.outside.table.hotspots[realm.plus.outside.table.hotspots[,2]>=overlap,]$model.state,pattern=realm.character,replacement='7')
    realm.plus.outside.table.hotspots[realm.plus.outside.table.hotspots$model.state==paste('7,','1',sep=''),'model.state']<-'1,7'
    realm.plus.outside.table.hotspots[realm.plus.outside.table.hotspots$model.state==paste('7,','2',sep=''),'model.state']<-'2,7'
    realm.plus.outside.table.hotspots[realm.plus.outside.table.hotspots$model.state==paste('7,','3',sep=''),'model.state']<-'3,7'
    realm.plus.outside.table.hotspots[realm.plus.outside.table.hotspots$model.state==paste('7,','4',sep=''),'model.state']<-'4,7'
    realm.plus.outside.table.hotspots[realm.plus.outside.table.hotspots$model.state==paste('7,','5',sep=''),'model.state']<-'5,7'
    realm.plus.outside.table.hotspots[realm.plus.outside.table.hotspots$model.state==paste('7,','6',sep=''),'model.state']<-'6,7'
    realm.table[realm.table$model.state%in%target.states,'model.state']<-realm.plus.outside.table.hotspots$model.state
  }
  
  table<-realm.table
  cat('state distribution','\n')
  cat(table(table$model.state),'\n')
  results.table<-table
  model5<-results.table$model.state
  table.model.state<-table(table$model.state)
  names(model5)<-results.table$spp
  name.check<-name.check(tree,model5)
  if(name.check=='OK'){
    tree.model<-tree
  }else{
    #drop tips in tree without data
    tree.model<-drop.tip(tree,name.check$tree_not_data) 
    #drop species not in tree from trait vector
    model5<-model5[!(names(model5) %in% name.check$data_not_tree)]
  }
  cat(name.check(tree.model,model5),'namecheck','\n')
  write.tree(tree.model,paste('./',name,'_7areas_inhotoutrealm.tree',sep=''))
  #create the table for BioGeoBEARS
  #1-> 100; 2-> 110; 3-> 010; 4-> 011; 5-> 001; 6->101
  
  geographytable<-data.frame(names(model5),unname(model5))
  colnames(geographytable)<-c('spp','model.state')
  geographytable$code<-'0000000'
  geographytable[geographytable$model.state=='1','code']<-'1000000'
  geographytable[geographytable$model.state=='2','code']<-'0100000'
  geographytable[geographytable$model.state=='3','code']<-'0010000'
  geographytable[geographytable$model.state=='4','code']<-'0001000'
  geographytable[geographytable$model.state=='5','code']<-'0000100'
  geographytable[geographytable$model.state=='6','code']<-'0000010'
  geographytable[geographytable$model.state=='7','code']<-'0000001'
  geographytable[geographytable$model.state=='1,7','code']<-'1000001'
  geographytable[geographytable$model.state=='2,7','code']<-'0100001'
  geographytable[geographytable$model.state=='3,7','code']<-'0010001'
  geographytable[geographytable$model.state=='4,7','code']<-'0001001'
  geographytable[geographytable$model.state=='5,7','code']<-'0000101'
  geographytable[geographytable$model.state=='6,7','code']<-'0000011'
  geographytable[geographytable$model.state=='1,2','code']<-'1100000'
  geographytable[geographytable$model.state=='1,3','code']<-'1010000'
  geographytable[geographytable$model.state=='1,4','code']<-'1001000'
  geographytable[geographytable$model.state=='1,5','code']<-'1000100'
  geographytable[geographytable$model.state=='1,6','code']<-'1000010'
  geographytable[geographytable$model.state=='2,3','code']<-'0110000'
  geographytable[geographytable$model.state=='2,4','code']<-'0101000'
  geographytable[geographytable$model.state=='2,5','code']<-'0100100'
  geographytable[geographytable$model.state=='2,6','code']<-'0100010'
  geographytable[geographytable$model.state=='3,4','code']<-'0011000'
  geographytable[geographytable$model.state=='3,5','code']<-'0010100'
  geographytable[geographytable$model.state=='3,6','code']<-'0010010'
  geographytable[geographytable$model.state=='4,5','code']<-'0001100'
  geographytable[geographytable$model.state=='4,6','code']<-'0001010'
  geographytable[geographytable$model.state=='5,6','code']<-'0000110'
  geographytable$model.state<-NULL
  header<-cbind(nrow(geographytable),7)
  write.table(header,paste('./',name,'_',replicate,'_7areas_inhotoutrealm_geographyfile.txt',sep=''),sep='\t',quote=F,row.names=F,col.names = F)
  write.table(geographytable,paste('./',name,'_',replicate,'_7areas_inhotoutrealm_geographyfile.txt',sep=''),append=T,sep=' ',quote=F,row.names=F,col.names = F)
  return(table.model.state)
  
}
#points is a vector with the points to get the overlap for
#returns a table with overlap with realm and hotspot for the set of points
#table is world.table
#name is the name of the realm (Afrotropical,Australasian,IndoMalay,Nearctic,Neotropical,Palearctic)
get_rangeoverlap_with_points<-function(points,name,table){
  if(!name%in%c('Afrotropical','Australasian','Indo-Malay','Nearctic','Neotropical','Palearctic')){
    cat('incorrect realm name','\n')
    break
  }
  world.table<-table
  #calculate proportion of species range inside hotspots
  grids<-list.files(path='./output/grids/',pattern='grid_.*_100.rds')
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
  
  #correct realm points (points are in reference to the realm but the table is worldwide)
  points<-grid.cells.df[grid.cells.df$realms.names==name,]$start.cell+points
  hotspots.cells.realm<-sort(points)
  world.table.cells<-as.character(world.table$cells)
  str(unlist(strsplit(world.table.cells[2],' ')))
  world.table.cells<-sapply(world.table.cells,function(x) unlist(strsplit(x,' ')))
  world.table.cells<-sapply(world.table.cells,function(x) x[x != ""])
  world.table.cells<-sapply(world.table.cells,function(x)as.numeric(x))
  #intersect each entry in world.table.cells with all elements in hotspots.cells.realm and record
  world.table.cells.hotspots<-lapply(world.table.cells,function(x) length(intersect(x,hotspots.cells.realm)))
  if(name=='Afrotropical'){
    world.table$range.cells.Afrotropical.hotspot<-unlist(world.table.cells.hotspots)
    world.table$afrotrop_WWF_hotspot.area<-world.table$range.cells.Afrotropical.hotspot/world.table$range.cells.Afrotropical
    world.table[is.na(world.table$afrotrop_WWF_hotspot.area),]$afrotrop_WWF_hotspot.area<-0
  }else if (name=='Australasian'){
    world.table$range.cells.Australasian.hotspot<-unlist(world.table.cells.hotspots)
    world.table$austral_WWF_hotspot.area<-world.table$range.cells.Australasian.hotspot/world.table$range.cells.Australasian
    world.table[is.na(world.table$austral_WWF_hotspot.area),]$austral_WWF_hotspot.area<-0
  }else if (name=='Indo-Malay'){
    world.table$range.cells.IndoMalay.hotspot<-unlist(world.table.cells.hotspots)
    world.table$indo_WWF_hotspot.area<-world.table$range.cells.IndoMalay.hotspot/world.table$range.cells.IndoMalay
    world.table[is.na(world.table$indo_WWF_hotspot.area),]$indo_WWF_hotspot.area<-0
  }else if (name=='Nearctic'){
    world.table$range.cells.Nearctic.hotspot<-unlist(world.table.cells.hotspots)
    world.table$nearctic_WWF_hotspot.area<-world.table$range.cells.Nearctic.hotspot/world.table$range.cells.Nearctic
    world.table[is.na(world.table$nearctic_WWF_hotspot.area),]$nearctic_WWF_hotspot.area<-0
  }else if (name=='Neotropical'){
    world.table$range.cells.Neotropical.hotspot<-unlist(world.table.cells.hotspots)
    world.table$neotrop_WWF_hotspot.area<-world.table$range.cells.Neotropical.hotspot/world.table$range.cells.Neotropical
    world.table[is.na(world.table$neotrop_WWF_hotspot.area),]$neotrop_WWF_hotspot.area<-0
  }else if (name=='Palearctic'){
    world.table$range.cells.Palearctic.hotspot<-unlist(world.table.cells.hotspots)
    world.table$palearctic_WWF_hotspot.area<-world.table$range.cells.Palearctic.hotspot/world.table$range.cells.Palearctic
    world.table[is.na(world.table$palearctic_WWF_hotspot.area),]$palearctic_WWF_hotspot.area<-0
  }
  return(world.table)  
}
####################################################################################################
####################################################################################################

####
#points.objectfile<-'./output/birds/immigration/controls/points.afrotropical.50.RDS'
#world.tablefile<-'./output/birds/tables/100_all_realms_ranges_plus_hotspots.txt'
#region.name is one of 'Afrotropical','Australasian','Indo-Malay','Neotropical','Nearctic','Palearctic'
#treefile<-'./output/birds/trees/birds_tree_IUCN.tree'
#name is one of 'afrotrop','austral','indo','nearctic','neotrop','palearctic'
#path<-'./output/birds/immigration/controls/'
#path is the relative path to store the control files ('./output/birds/immigration/controls/')

get_table_from_points<-function(points.objectfile,world.tablefile,region.name,treefile,name,path){
  dir.create(paste(path,name,sep=''))
  points<-readRDS(points.objectfile)
  world.table<-read.table(world.tablefile,header=T,sep='\t',stringsAsFactors = F)
  #drop columns with hotspot info in world.tablefile
  world.table<-world.table[,-grep('hotspot',colnames(world.table))]
  points.table<-lapply(points,function(x) get_rangeoverlap_with_points(points=x,name=region.name,table=world.table))
  #collate tree and table for austral, check species distributions
  tree<-read.tree(treefile)
  tree<-drop.tip(tree,setdiff(tree$tip.label,points.table[[1]]$spp))
  points.table<-lapply(points.table,function(x) x[x$spp%in%tree$tip.label,])
  #keep columns with proportion of range inside realm and proportion of range within hotspot
  points.table<-lapply(points.table,function(x){drop.columns<-grep('cells',colnames(x));x<-x[,-drop.columns];colnames(x)<-sub(colnames(x),pattern='_WWF',replacement='');return(x)})
  points.control.input<-list()
  for(i in 1:length(points.table)){
    name<-name
    overlap<-0.8
    folder<-paste(path,name,'/',i,'/',sep='')
    dir.create(folder)
    setwd(folder)
    points.control.input[[i]]<-prepare_realm_input_7areas_plus_inhotoutrealm(table=points.table[[i]],name=name,overlap=0.80,replicate=i,tree=tree)
    setwd('/home/ji247/hotspots_vertebrates/')
  }
  saveRDS(points.control.input,file=paste(path,name,'points.control.input.rds',sep=''))
}

get_table_from_points(points.objectfile='./output/birds/immigration/controls/points.afrotropical.50.RDS',world.tablefile='./output/birds/tables/100_all_realms_ranges_plus_hotspots.txt',region.name='Afrotropical',treefile='./output/birds/trees/birds_tree_IUCN.tree',name='afrotrop',path='./output/birds/immigration/controls/')
get_table_from_points(points.objectfile='./output/birds/immigration/controls/points.austral.50.RDS',world.tablefile='./output/birds/tables/100_all_realms_ranges_plus_hotspots.txt',region.name='Australasian',treefile='./output/birds/trees/birds_tree_IUCN.tree',name='austral',path='./output/birds/immigration/controls/')
get_table_from_points(points.objectfile='./output/birds/immigration/controls/points.indo.50.RDS',world.tablefile='./output/birds/tables/100_all_realms_ranges_plus_hotspots.txt',region.name='Indo-Malay',treefile='./output/birds/trees/birds_tree_IUCN.tree',name='indo',path='./output/birds/immigration/controls/')
get_table_from_points(points.objectfile='./output/birds/immigration/controls/points.nearctic.50.RDS',world.tablefile='./output/birds/tables/100_all_realms_ranges_plus_hotspots.txt',region.name='Nearctic',treefile='./output/birds/trees/birds_tree_IUCN.tree',name='nearctic',path='./output/birds/immigration/controls/')
get_table_from_points(points.objectfile='./output/birds/immigration/controls/points.neotropical.50.RDS',world.tablefile='./output/birds/tables/100_all_realms_ranges_plus_hotspots.txt',region.name='Neotropical',treefile='./output/birds/trees/birds_tree_IUCN.tree',name='neotrop',path='./output/birds/immigration/controls/')
get_table_from_points(points.objectfile='./output/birds/immigration/controls/points.palearctic.50.RDS',world.tablefile='./output/birds/tables/100_all_realms_ranges_plus_hotspots.txt',region.name='Palearctic',treefile='./output/birds/trees/birds_tree_IUCN.tree',name='palearctic',path='./output/birds/immigration/controls/')

########plot number of species in each category for afrotrop
####afrotrop.realm<-read.table('./output/birds/new_tables_hotspots/whole_realms/afrotrop/afrotrop_geographyfile.txt',header=F,sep=' ',skip = 1)
####afrotrop.realm.table<-table(afrotrop.realm$V2)
####afrotrop.control.input.points<-list()
####for(i in 1:5){
####  afrotrop.control.input.points[[i]]<-unlist(lapply(afrotrop.control.input,function(x)x[i]))
####}
####par(mfrow=c(2,2))
####for(i in 1:4){
####    if(i==1){
####    hist(afrotrop.control.input.points[[i]],xlim=c(0,max(c(afrotrop.control.input.points[[i]],afrotrop.realm.table['100']))),main='afrotrop hot endemics')
####    abline(v=afrotrop.realm.table['100'],col='red')
####  }else if(i==2){
####    hist(afrotrop.control.input.points[[i]],xlim=c(0,max(c(afrotrop.control.input.points[[i]],afrotrop.realm.table['110']))),main='afrotrop hot/nonhot')
####    abline(v=afrotrop.realm.table['110'],col='red')
####  }else if(i==3){
####    hist(afrotrop.control.input.points[[i]],xlim=c(0,max(c(afrotrop.control.input.points[[i]],afrotrop.realm.table['10']))),main='afrotrop nonhot')
####    abline(v=afrotrop.realm.table['10'],col='red')
####  }else if(i==4){
####    hist(afrotrop.control.input.points[[i]],xlim=c(0,max(c(afrotrop.control.input.points[[i]],afrotrop.realm.table['11']))),main='afrotrop nonhot/world')
####    abline(v=afrotrop.realm.table['11'],col='red')
####  }
####}
####dev.off()
####
##########plot number of species in each category for austral
####austral.realm<-read.table('./output/birds/new_tables_hotspots/whole_realms/austral/austral_geographyfile.txt',header=F,sep=' ',skip = 1)
####austral.realm.table<-table(austral.realm$V2)
####austral.control.input.points<-list()
####for(i in 1:5){
####  austral.control.input.points[[i]]<-unlist(lapply(austral.control.input,function(x)x[i]))
####}
####par(mfrow=c(2,2))
####for(i in 1:4){
####  if(i==1){
####    hist(austral.control.input.points[[i]],xlim=c(0,max(c(austral.control.input.points[[i]],austral.realm.table['100']))),main='austral hot endemics')
####    abline(v=austral.realm.table['100'],col='red')
####  }else if(i==2){
####    hist(austral.control.input.points[[i]],xlim=c(0,max(c(austral.control.input.points[[i]],austral.realm.table['110']))),main='austral hot/nonhot')
####    abline(v=austral.realm.table['110'],col='red')
####  }else if(i==3){
####    hist(austral.control.input.points[[i]],xlim=c(0,max(c(austral.control.input.points[[i]],austral.realm.table['10']))),main='austral nonhot')
####    abline(v=austral.realm.table['10'],col='red')
####  }else if(i==4){
####    hist(austral.control.input.points[[i]],xlim=c(0,max(c(austral.control.input.points[[i]],austral.realm.table['11']))),main='austral nonhot/world')
####    abline(v=austral.realm.table['11'],col='red')
####  }
####}
####dev.off()
####
##########plot number of species in each category for indo
#####no whole realm analysis for indo, so skipping
####
##########plot number of species in each category for nearctic
####nearctic.realm<-read.table('./output/birds/new_tables_hotspots/whole_realms/nearctic/nearctic_geographyfile.txt',header=F,sep=' ',skip = 1)
####nearctic.realm.table<-table(nearctic.realm$V2)
####nearctic.control.input.points<-list()
####for(i in 1:5){
####  nearctic.control.input.points[[i]]<-unlist(lapply(nearctic.control.input,function(x)x[i]))
####}
####par(mfrow=c(2,2))
####for(i in 1:4){
####  if(i==1){
####    hist(nearctic.control.input.points[[i]],xlim=c(0,max(c(nearctic.control.input.points[[i]],nearctic.realm.table['100']))),main='nearctic hot endemics')
####    abline(v=nearctic.realm.table['100'],col='red')
####  }else if(i==2){
####    hist(nearctic.control.input.points[[i]],xlim=c(0,max(c(nearctic.control.input.points[[i]],nearctic.realm.table['110']))),main='nearctic hot/nonhot')
####    abline(v=nearctic.realm.table['110'],col='red')
####  }else if(i==3){
####    hist(nearctic.control.input.points[[i]],xlim=c(0,max(c(nearctic.control.input.points[[i]],nearctic.realm.table['10']))),main='nearctic nonhot')
####    abline(v=nearctic.realm.table['10'],col='red')
####  }else if(i==4){
####    hist(nearctic.control.input.points[[i]],xlim=c(0,max(c(nearctic.control.input.points[[i]],nearctic.realm.table['11']))),main='nearctic nonhot/world')
####    abline(v=nearctic.realm.table['11'],col='red')
####  }
####}
####dev.off()
####
####
##########plot number of species in each category for neotrop
####neotrop.realm<-read.table('./output/birds/new_tables_hotspots/whole_realms/neotrop/neotrop_geographyfile.txt',header=F,sep=' ',skip = 1)
####neotrop.realm.table<-table(neotrop.realm$V2)
####neotrop.control.input.points<-list()
####for(i in 1:5){
####  neotrop.control.input.points[[i]]<-unlist(lapply(neotrop.control.input,function(x)x[i]))
####}
####par(mfrow=c(2,2))
####for(i in 1:4){
####  if(i==1){
####    hist(neotrop.control.input.points[[i]],xlim=c(0,max(c(neotrop.control.input.points[[i]],neotrop.realm.table['100']))),main='neotrop hot endemics')
####    abline(v=neotrop.realm.table['100'],col='red')
####  }else if(i==2){
####    hist(neotrop.control.input.points[[i]],xlim=c(0,max(c(neotrop.control.input.points[[i]],neotrop.realm.table['110']))),main='neotrop hot/nonhot')
####    abline(v=neotrop.realm.table['110'],col='red')
####  }else if(i==3){
####    hist(neotrop.control.input.points[[i]],xlim=c(0,max(c(neotrop.control.input.points[[i]],neotrop.realm.table['10']))),main='neotrop nonhot')
####    abline(v=neotrop.realm.table['10'],col='red')
####  }else if(i==4){
####    hist(neotrop.control.input.points[[i]],xlim=c(0,max(c(neotrop.control.input.points[[i]],neotrop.realm.table['11']))),main='neotrop nonhot/world')
####    abline(v=neotrop.realm.table['11'],col='red')
####  }
####}
####dev.off()
####
##########plot number of species in each category for palearctic
####palearctic.realm<-read.table('./output/birds/new_tables_hotspots/whole_realms/palearctic/palearctic_geographyfile.txt',header=F,sep=' ',skip = 1)
####palearctic.realm.table<-table(palearctic.realm$V2)
####palearctic.control.input.points<-list()
####for(i in 1:5){
####  palearctic.control.input.points[[i]]<-unlist(lapply(palearctic.control.input,function(x)x[i]))
####}
####par(mfrow=c(2,2))
####for(i in 1:4){
####  if(i==1){
####    hist(palearctic.control.input.points[[i]],xlim=c(0,max(c(palearctic.control.input.points[[i]],palearctic.realm.table['100']))),main='palearctic hot endemics')
####    abline(v=palearctic.realm.table['100'],col='red')
####  }else if(i==2){
####    hist(palearctic.control.input.points[[i]],xlim=c(0,max(c(palearctic.control.input.points[[i]],palearctic.realm.table['110']))),main='palearctic hot/nonhot')
####    abline(v=palearctic.realm.table['110'],col='red')
####  }else if(i==3){
####    hist(palearctic.control.input.points[[i]],xlim=c(0,max(c(palearctic.control.input.points[[i]],palearctic.realm.table['10']))),main='palearctic nonhot')
####    abline(v=palearctic.realm.table['10'],col='red')
####  }else if(i==4){
####    hist(palearctic.control.input.points[[i]],xlim=c(0,max(c(palearctic.control.input.points[[i]],palearctic.realm.table['11']))),main='palearctic nonhot/world')
####    abline(v=palearctic.realm.table['11'],col='red')
####  }
####}
####dev.off()
##
##realm.name<-'afrotrop'
##number<-'1'
##run_whole_realm_BioGeoBEARS_plusBSM<-function(realm.name,number){
##  setwd('~/hotspots_vertebrates/')
##  source('./BioGeoBEARS_run_models.R')
##  
##  path<-paste('./output/birds/new_tables_hotspots/immigration/controls/',realm.name,'/',number,'/',sep='')
##  names.trees<-list.files(path,pattern='.tree')
##  names.geography<-list.files(path,pattern='_geographyfile.txt')
##  if(length(names.trees)==0){
##    return('no trees in this folder')
##    
##  }
##  if(length(names.trees)!=length(names.geography)){
##    return('different lenghts, double check folder')
##    
##  }
##  for(i in 1:length(names.trees)){
##    run_BioGeoBEARS_models(treefile = names.trees[i],geographyfile =names.geography[i],path = path,name=sub(names.trees[i],pattern=paste('_',realm.name,'.tree',sep=''),replacement='') )
##    setwd('~/Desktop/HOTSPOTS/repositories/hotspots_vertebrates/')
##  }
##  list.aicc.files<-list.files(path,pattern='_AICc_rellike_formatted.txt')
##  #read AICC tables, select the best fitting model (sort by AICc weight), store the best fitting models in a list
##  model.list<-list(0)
##  for (i in 1:length(list.aicc.files)){
##    table<-read.table(file=paste(path,list.aicc.files[i],sep=''),header=T,sep='\t')
##    table<-table[order(-table$AICc_wt),]
##    cat(i,'   ',table[1,'AICc_wt'],'    ',row.names(table)[1],'\n')
##    model.list[[i]]<-c(row.names(table)[1],table[1,'AICc_wt'])
##  }
##  results.BSM<-list(0)
##  for (i in 1:length(model.list)){
##    results.BSM[[i]]<-run_BioGeoBEARS_selectedmodel_BSM(treefile=names.trees[i],geographyfile=names.geography[i],path=path,name=sub(names.trees[i],pattern='_neotropical.tree',replacement=''),model_name=model.list[[i]][1]) 
##    setwd('/Users/javier/Desktop/HOTSPOTS/repositories/hotspots_vertebrates')
##  }
##  saveRDS(results.BSM,file=paste(path,realm.name,'_',number,'_BSM.RDS',sep=''))
##  
##}
##

##
##
##