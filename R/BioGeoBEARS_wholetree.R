library(ape)
library(plyr)
library(phangorn)
library(geiger)
library(parallel)
library(diversitree)
###############################################
#####TO DOS####
#####remove prepare_realm_simple_input_5states
#####remove prepare_realm_simple_input_6states
#####remove prepare_realm_input_7areas
#####remove run_whole_realm_BioGeoBEARS_plusBSM
#####run_whole_realm_BioGeoBEARS_plusBSM_7areas
prepare_realm_simple_input_5states<-function(table,name,overlap){
  
  #states defined as:
  #1- species occurring only in hotspot in region
  #2- species occurring both in and outside hotspot in region
  #3- species occurring only outside hotspot in region
  #4- species occurring both outside hotspot in region and outside region
  #5- species occurring only outside region
  table<-table[,c('spp',paste(name,'_realm.area',sep=''),paste(name,'_hotspot.area',sep=''))]
  colnames(table)<-c('spp','realm.area','hotspot.area')
  table[is.na(table)]<-0
  table$model.state<-0
  table[table$realm.area<(1-overlap),'model.state']<-5
  table[(table$realm.area>(1-overlap))&(table$realm.area<overlap),'model.state']<-4
  table[(table$hotspot.area>=overlap)&(table$realm.area>=overlap),'model.state']<-1
  table[(table$hotspot.area<overlap)&(table$hotspot.area>(1-overlap))&(table$realm.area>=overlap),'model.state']<-2
  table[(table$hotspot.area<(1-overlap))&(table$realm.area>=overlap),'model.state']<-3
  if(length(table(table$model.state))<5){
    cat('not all states are present','\n')
    results.transitions.df<-NA
    return(results.transitions.df)
  }
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
  name.check(tree.model,model5)
  write.tree(tree.model,paste('./mammals_Rolland_terrestrial_IUCN_',name,'.tree',sep=''))
  #create the table for BioGeoBEARS
  #1-> 100; 2-> 110; 3-> 010; 4-> 011; 5-> 001
  
  geographytable<-data.frame(names(model5),unname(model5))
  colnames(geographytable)<-c('spp','code')
  geographytable[geographytable$code==1,]$code<-'100'
  geographytable[geographytable$code==2,]$code<-'110'
  geographytable[geographytable$code==3,]$code<-'010'
  geographytable[geographytable$code==4,]$code<-'011'
  geographytable[geographytable$code==5,]$code<-'001'
  header<-cbind(nrow(geographytable),3)
  write.table(header,paste('./',name,'_geographyfile.txt',sep=''),sep='\t',quote=F,row.names=F,col.names = F)
  write.table(geographytable,paste('./',name,'_geographyfile.txt',sep=''),append=T,sep=' ',quote=F,row.names=F,col.names = F)
  return(table.model.state)
  
}

prepare_realm_simple_input_6states<-function(table,name,overlap){
  
  #states defined as:
  #1- species occurring only in hotspot in region
  #2- species occurring both in and outside hotspot in region
  #3- species occurring only outside hotspot in region
  #4- species occurring both outside hotspot in region and outside region
  #6- species occurring only outside region
  #5- species occurring only outside region + inside hotspot
  table<-table[,c('spp',paste(name,'_realm.area',sep=''),paste(name,'_hotspot.area',sep=''))]
  colnames(table)<-c('spp','realm.area','hotspot.area')
  table[is.na(table)]<-0
  table$model.state<-0
  table[(table$realm.area<(1-overlap)),'model.state']<-6
  table[(table$realm.area>(1-overlap))&(table$realm.area<overlap)&(table$hotspot.area<overlap),'model.state']<-4
  table[(table$realm.area>(1-overlap))&(table$realm.area<overlap)&(table$hotspot.area>=overlap),'model.state']<-5
  table[(table$hotspot.area>=overlap)&(table$realm.area>=overlap),'model.state']<-1
  table[(table$hotspot.area<overlap)&(table$hotspot.area>(1-overlap))&(table$realm.area>=overlap),'model.state']<-2
  table[(table$hotspot.area<(1-overlap))&(table$realm.area>=overlap),'model.state']<-3
  if(length(table(table$model.state))<6){
    cat('not all states are present','\n')
   # results.transitions.df<-NA
  #  return(results.transitions.df)
  }
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
  name.check(tree.model,model5)
  write.tree(tree.model,paste('./mammals_Rolland_terrestrial_IUCN_',name,'_6states.tree',sep=''))
  #create the table for BioGeoBEARS
  #1-> 100; 2-> 110; 3-> 010; 4-> 011; 5-> 001; 6->101
  
  geographytable<-data.frame(names(model5),unname(model5))
  colnames(geographytable)<-c('spp','code')
  geographytable[geographytable$code==1,'code']<-'100'
  geographytable[geographytable$code==2,'code']<-'110'
  geographytable[geographytable$code==3,'code']<-'010'
  geographytable[geographytable$code==4,'code']<-'011'
  geographytable[geographytable$code==5,'code']<-'101'
  geographytable[geographytable$code==6,'code']<-'001'
  header<-cbind(nrow(geographytable),3)
  write.table(header,paste('./',name,'_6states_geographyfile.txt',sep=''),sep='\t',quote=F,row.names=F,col.names = F)
  write.table(geographytable,paste('./',name,'_6states_geographyfile.txt',sep=''),append=T,sep=' ',quote=F,row.names=F,col.names = F)
  return(table.model.state)
  
}




prepare_realm_input_7areas_plus_inhotoutrealm<-function(table,name,overlap,tree){
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
  realm.character<-as.character(grep(name,colnames(realm.table)[-c(1,8)]))
  realm.hotspot.table<-table[,grep(name,colnames(table))]
  realm.hotspot.table<-realm.hotspot.table[which(realm.hotspot.table[,1]>=overlap),]
  realm.hotspot.table$model.state<-0
  realm.hotspot.table[realm.hotspot.table[,2]<(1-overlap),'model.state']<-realm.character
  realm.hotspot.table[realm.hotspot.table[,2]>=overlap,'model.state']<-7
  realm.hotspot.table[realm.hotspot.table[,2]>(1-overlap)&realm.hotspot.table[,2]<overlap,'model.state']<-paste(realm.character,7,sep=',')
  realm.table[realm.table$model.state==realm.character,'model.state']<-realm.hotspot.table$model.state
  #select species occurring in target realm plus another realm to check if  
  target.states<-names(table(realm.table$model.state))[grep(realm.character,names(table(realm.table$model.state)))]
  target.states<-target.states[-grep('7',target.states)]
  target.states<-target.states[-match(realm.character,target.states)]
  realm.plus.outside.table<-realm.table[realm.table$model.state%in%target.states,]
  #add hotspot overlap within realm to realm.plus.outside table
  realm.plus.outside.table.hotspots<-as.data.frame(cbind(table[table$spp%in%realm.plus.outside.table$spp,'spp'],table[table$spp%in%realm.plus.outside.table$spp,paste(name,'_hotspot.area',sep='')]),stringsAsFactors = F)
  colnames(realm.plus.outside.table.hotspots)<-c('spp',paste(name,'_hotspot.area',sep=''))
  realm.plus.outside.table.hotspots$model.state<-realm.table[realm.table$model.state%in%target.states,]$model.state
  realm.plus.outside.table.hotspots[realm.plus.outside.table.hotspots[,2]>=overlap,'model.state']<-sub(realm.plus.outside.table.hotspots[realm.plus.outside.table.hotspots[,2]>=overlap,'model.state'],pattern=realm.character,replacement='7')
  realm.plus.outside.table.hotspots[realm.plus.outside.table.hotspots$model.state==paste('7,','1',sep=''),'model.state']<-'1,7'
  realm.plus.outside.table.hotspots[realm.plus.outside.table.hotspots$model.state==paste('7,','2',sep=''),'model.state']<-'2,7'
  realm.plus.outside.table.hotspots[realm.plus.outside.table.hotspots$model.state==paste('7,','3',sep=''),'model.state']<-'3,7'
  realm.plus.outside.table.hotspots[realm.plus.outside.table.hotspots$model.state==paste('7,','4',sep=''),'model.state']<-'4,7'
  realm.plus.outside.table.hotspots[realm.plus.outside.table.hotspots$model.state==paste('7,','5',sep=''),'model.state']<-'5,7'
  realm.plus.outside.table.hotspots[realm.plus.outside.table.hotspots$model.state==paste('7,','6',sep=''),'model.state']<-'6,7'
  realm.table[realm.table$model.state%in%target.states,]$model.state<-realm.plus.outside.table.hotspots$model.state
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
  write.table(header,paste('./',name,'_7areas_inhotoutrealm_geographyfile.txt',sep=''),sep='\t',quote=F,row.names=F,col.names = F)
  write.table(geographytable,paste('./',name,'_7areas_inhotoutrealm_geographyfile.txt',sep=''),append=T,sep=' ',quote=F,row.names=F,col.names = F)
  return(table.model.state)
  
}

###########################
#this function generates a table with the range of each species (in cells) in the world + in each realm and also adds ranges in predefined hotspot cells (and proportion of ranges in and outside hotspots)
#path<-'./output/mammals/new_tables_hotspots/'
#world.table.file<-'./output/mammals/new_tables_hotspots/100_all_realms_species_gridoccurrence_table.txt'
#weighted.endemism.file<-'./output/mammals/new_tables_hotspots/100_all_realms_realms_richness_wend_grid_table.txt'
#quantile.hotspot<-0.80
#path.grids<-'./output/mammals/'
build_world_range_hotspots_table<-function(path,world.table.file,weighted.endemism.file,quantile.hotspot,path.grids){
  list<-list.files(path=path,pattern='100_.*_realms_species_gridoccurrence_table.txt')
  #removes all_realm_file
  if(length(grep('all_realms',list))>0){
    list<-list[-grep('all_realms',list)]
  }
  #remove Oceanic
  if(length(grep('Oceanic',list))>0){
    list<-list[-grep('Oceanic',list)]
  }
  #there should be 6 files (Afrotropical,Australasian,Indo-Malay,Nearctic,Neotropical,Palearctic)
  #if not, exit
  if(!length(list)==6){
    return('not 6 files for 6 realms in this folder!')
  }
  #read region tables and paste them to data frames in a list
  all.regions<-lapply(list, function(x) read.table(paste(path,x,sep=''),header=T,sep='\t'))
  #reorder dataframes by world.table
  world.table<-read.table(world.table.file,header=T,sep='\t')
  all.regions<-lapply(all.regions,function(x) x[match(world.table$spp,x$spp),])
  all.regions.ranges<-lapply(all.regions,function(x)x[,3])
  world.table$range.cells.Afrotropical<-all.regions.ranges[[1]]
  world.table$range.cells.Australasian<-all.regions.ranges[[2]]
  world.table$range.cells.IndoMalay<-all.regions.ranges[[3]]
  world.table$range.cells.Nearctic<-all.regions.ranges[[4]]
  world.table$range.cells.Neotropical<-all.regions.ranges[[5]]
  world.table$range.cells.Palearctic<-all.regions.ranges[[6]]
  #get proportions of range in each realm
  world.table$afrotrop_WWF_realm.area<-world.table$range.cells.Afrotropical/world.table$range.cells
  world.table$austral_WWF_realm.area<-world.table$range.cells.Australasian/world.table$range.cells
  world.table$indo_WWF_realm.area<-world.table$range.cells.IndoMalay/world.table$range.cells
  world.table$nearctic_WWF_realm.area<-world.table$range.cells.Nearctic/world.table$range.cells
  world.table$neotrop_WWF_realm.area<-world.table$range.cells.Neotropical/world.table$range.cells
  world.table$palearctic_WWF_realm.area<-world.table$range.cells.Palearctic/world.table$range.cells
  #define the hotspots by weighted.endemism
  wend<-read.table(weighted.endemism.file,header=T,sep='\t',stringsAsFactors = F)
  wend<-na.omit(wend)
  colnames(wend)[1]<-'cells'
  quantile<-quantile.hotspot
  hotspots<-wend[wend$number.of.species.wend>quantile(wend$number.of.species.wend,quantile),]$cells
  #calculate proportion of species range inside hotspots
  grids<-list.files(path=path.grids,pattern='grid_.*_100.rds')
  grid.world<-grids[grep('World_RealmsMerged',grids)]
  grid.realms.names<-grids[-grep('World_RealmsMerged',grids)]
  grid.world<-readRDS(paste(path.grids,'/',grid.world,sep=''))
  grid.realms<-lapply(grid.realms.names,function(x) readRDS(paste(path.grids,'/',x,sep='')))
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
  #make a hotspot list per realm
  hotspots.cells.realm<-list()
  for (i in 1:nrow(grid.cells.df)){
    hotspots.cells.realm[[i]]<-hotspots[intersect(which(hotspots>=grid.cells.df[i,'start.cell']),which(hotspots<=grid.cells.df[i,'end.cell']))]
  }
  #get all the worldwide ranges in a list
  world.table.cells<-as.character(world.table$cells)
  str(unlist(strsplit(world.table.cells[2],' ')))
  world.table.cells<-sapply(world.table.cells,function(x) unlist(strsplit(x,' ')))
  world.table.cells<-sapply(world.table.cells,function(x) x[x != ""])
  world.table.cells<-sapply(world.table.cells,function(x)as.numeric(x))
  #intersect each entry in world.table.cells with all elements in hotspots.cells.realm and record
  world.table.cells.hotspots<-list()
  for (i in 1:length(world.table.cells)){
    world.table.cells.hotspots[[i]]<-unlist(lapply(hotspots.cells.realm,function(x) length(intersect(x,world.table.cells[[i]]))))
  }
  world.table$range.cells.Afrotropical.hotspot<-unlist(lapply(world.table.cells.hotspots,function(x)x[1]))
  world.table$range.cells.Australasian.hotspot<-unlist(lapply(world.table.cells.hotspots,function(x)x[2]))
  world.table$range.cells.IndoMalay.hotspot<-unlist(lapply(world.table.cells.hotspots,function(x)x[3]))
  world.table$range.cells.Nearctic.hotspot<-unlist(lapply(world.table.cells.hotspots,function(x)x[4]))
  world.table$range.cells.Neotropical.hotspot<-unlist(lapply(world.table.cells.hotspots,function(x)x[5]))
  world.table$range.cells.Palearctic.hotspot<-unlist(lapply(world.table.cells.hotspots,function(x)x[6]))
  
  world.table$afrotrop_WWF_hotspot.area<-world.table$range.cells.Afrotropical.hotspot/world.table$range.cells.Afrotropical
  world.table$austral_WWF_hotspot.area<-world.table$range.cells.Australasian.hotspot/world.table$range.cells.Australasian
  world.table$indo_WWF_hotspot.area<-world.table$range.cells.IndoMalay.hotspot/world.table$range.cells.IndoMalay
  world.table$nearctic_WWF_hotspot.area<-world.table$range.cells.Nearctic.hotspot/world.table$range.cells.Nearctic
  world.table$neotrop_WWF_hotspot.area<-world.table$range.cells.Neotropical.hotspot/world.table$range.cells.Neotropical
  world.table$palearctic_WWF_hotspot.area<-world.table$range.cells.Palearctic.hotspot/world.table$range.cells.Palearctic
  world.table[is.na(world.table)] <- 0
  write.table(world.table,file=paste(path,'/100_all_realms_ranges_plus_hotspots.txt',sep=''),sep='\t',quote=F,row.names=F)
}

#this function generates a table with the range of each species (in cells) in the world + in each realm and also adds ranges in predefined hotspot cells (and proportion of ranges in and outside hotspots)
#path<-'./output/mammals/new_tables_hotspots/'
#world.table.file<-'./output/mammals/new_tables_hotspots/100_all_realms_species_gridoccurrence_table.txt'
#hotspots.file<-'./output/mammals/new_tables_hotspots/100_all_realms_realms_richness_wend_grid_table.txt'
#quantile.hotspot<-0.80
#path.grids<-'./output/mammals/'

build_world_range_hotspots_table_predefined<-function(path,world.table.file,hotspots.file,quantile.hotspot,path.grids){
  list<-list.files(path=path,pattern='100_.*_realms_species_gridoccurrence_table.txt')
  #removes all_realm_file
  if(length(grep('all_realms',list))>0){
    list<-list[-grep('all_realms',list)]
  }
  #remove Oceanic
  if(length(grep('Oceanic',list))>0){
    list<-list[-grep('Oceanic',list)]
  }
  #there should be 6 files (Afrotropical,Australasian,Indo-Malay,Nearctic,Neotropical,Palearctic)
  #if not, exit
  if(!length(list)==6){
    return('not 6 files for 6 realms in this folder!')
  }
  #read region tables and paste them to data frames in a list
  all.regions<-lapply(list, function(x) read.table(paste(path,x,sep=''),header=T,sep='\t'))
  #reorder dataframes by world.table
  world.table<-read.table(world.table.file,header=T,sep='\t')
  all.regions<-lapply(all.regions,function(x) x[match(world.table$spp,x$spp),])
  all.regions.ranges<-lapply(all.regions,function(x)x[,3])
  world.table$range.cells.Afrotropical<-all.regions.ranges[[1]]
  world.table$range.cells.Australasian<-all.regions.ranges[[2]]
  world.table$range.cells.IndoMalay<-all.regions.ranges[[3]]
  world.table$range.cells.Nearctic<-all.regions.ranges[[4]]
  world.table$range.cells.Neotropical<-all.regions.ranges[[5]]
  world.table$range.cells.Palearctic<-all.regions.ranges[[6]]
  #get proportions of range in each realm
  world.table$afrotrop_WWF_realm.area<-world.table$range.cells.Afrotropical/world.table$range.cells
  world.table$austral_WWF_realm.area<-world.table$range.cells.Australasian/world.table$range.cells
  world.table$indo_WWF_realm.area<-world.table$range.cells.IndoMalay/world.table$range.cells
  world.table$nearctic_WWF_realm.area<-world.table$range.cells.Nearctic/world.table$range.cells
  world.table$neotrop_WWF_realm.area<-world.table$range.cells.Neotropical/world.table$range.cells
  world.table$palearctic_WWF_realm.area<-world.table$range.cells.Palearctic/world.table$range.cells
  #define the hotspots by weighted.endemism
  hotspots.table<-read.table(hotspots.file,header=T,sep='\t',stringsAsFactors = F)
  colnames(hotspots.table)<-c('cells','hotspot')
  hotspots<-as.numeric(hotspots.table[hotspots.table$hotspot==1,]$cells)
  #calculate proportion of species range inside hotspots
  grids<-list.files(path=path.grids,pattern='grid_.*_100.rds')
  grid.world<-grids[grep('World_RealmsMerged',grids)]
  grid.realms.names<-grids[-grep('World_RealmsMerged',grids)]
  grid.world<-readRDS(paste(path.grids,'/',grid.world,sep=''))
  grid.realms<-lapply(grid.realms.names,function(x) readRDS(paste(path.grids,'/',x,sep='')))
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
  #make a hotspot list per realm
  hotspots.cells.realm<-list()
  for (i in 1:nrow(grid.cells.df)){
    hotspots.cells.realm[[i]]<-hotspots[intersect(which(hotspots>=grid.cells.df[i,'start.cell']),which(hotspots<=grid.cells.df[i,'end.cell']))]
  }
  #get all the worldwide ranges in a list
  world.table.cells<-as.character(world.table$cells)
  str(unlist(strsplit(world.table.cells[2],' ')))
  world.table.cells<-sapply(world.table.cells,function(x) unlist(strsplit(x,' ')))
  world.table.cells<-sapply(world.table.cells,function(x) x[x != ""])
  world.table.cells<-sapply(world.table.cells,function(x)as.numeric(x))
  #intersect each entry in world.table.cells with all elements in hotspots.cells.realm and record
  world.table.cells.hotspots<-list()
  for (i in 1:length(world.table.cells)){
    world.table.cells.hotspots[[i]]<-unlist(lapply(hotspots.cells.realm,function(x) length(intersect(x,world.table.cells[[i]]))))
  }
  world.table$range.cells.Afrotropical.hotspot<-unlist(lapply(world.table.cells.hotspots,function(x)x[1]))
  world.table$range.cells.Australasian.hotspot<-unlist(lapply(world.table.cells.hotspots,function(x)x[2]))
  world.table$range.cells.IndoMalay.hotspot<-unlist(lapply(world.table.cells.hotspots,function(x)x[3]))
  world.table$range.cells.Nearctic.hotspot<-unlist(lapply(world.table.cells.hotspots,function(x)x[4]))
  world.table$range.cells.Neotropical.hotspot<-unlist(lapply(world.table.cells.hotspots,function(x)x[5]))
  world.table$range.cells.Palearctic.hotspot<-unlist(lapply(world.table.cells.hotspots,function(x)x[6]))
  
  world.table$afrotrop_WWF_hotspot.area<-world.table$range.cells.Afrotropical.hotspot/world.table$range.cells.Afrotropical
  world.table$austral_WWF_hotspot.area<-world.table$range.cells.Australasian.hotspot/world.table$range.cells.Australasian
  world.table$indo_WWF_hotspot.area<-world.table$range.cells.IndoMalay.hotspot/world.table$range.cells.IndoMalay
  world.table$nearctic_WWF_hotspot.area<-world.table$range.cells.Nearctic.hotspot/world.table$range.cells.Nearctic
  world.table$neotrop_WWF_hotspot.area<-world.table$range.cells.Neotropical.hotspot/world.table$range.cells.Neotropical
  world.table$palearctic_WWF_hotspot.area<-world.table$range.cells.Palearctic.hotspot/world.table$range.cells.Palearctic
  world.table[is.na(world.table)] <- 0
  write.table(world.table,file=paste(path,'/100_all_realms_ranges_plus_hotspots_SRlatSARLM.txt',sep=''),sep='\t',quote=F,row.names=F)
}

build_world_range_hotspots_table_predefined_lmSRlat<-function(path,world.table.file,hotspots.file,quantile.hotspot,path.grids){
  list<-list.files(path=path,pattern='100_.*_realms_species_gridoccurrence_table.txt')
  #removes all_realm_file
  if(length(grep('all_realms',list))>0){
    list<-list[-grep('all_realms',list)]
  }
  #remove Oceanic
  if(length(grep('Oceanic',list))>0){
    list<-list[-grep('Oceanic',list)]
  }
  #there should be 6 files (Afrotropical,Australasian,Indo-Malay,Nearctic,Neotropical,Palearctic)
  #if not, exit
  if(!length(list)==6){
    return('not 6 files for 6 realms in this folder!')
  }
  #read region tables and paste them to data frames in a list
  all.regions<-lapply(list, function(x) read.table(paste(path,x,sep=''),header=T,sep='\t'))
  #reorder dataframes by world.table
  world.table<-read.table(world.table.file,header=T,sep='\t')
  all.regions<-lapply(all.regions,function(x) x[match(world.table$spp,x$spp),])
  all.regions.ranges<-lapply(all.regions,function(x)x[,3])
  world.table$range.cells.Afrotropical<-all.regions.ranges[[1]]
  world.table$range.cells.Australasian<-all.regions.ranges[[2]]
  world.table$range.cells.IndoMalay<-all.regions.ranges[[3]]
  world.table$range.cells.Nearctic<-all.regions.ranges[[4]]
  world.table$range.cells.Neotropical<-all.regions.ranges[[5]]
  world.table$range.cells.Palearctic<-all.regions.ranges[[6]]
  #get proportions of range in each realm
  world.table$afrotrop_WWF_realm.area<-world.table$range.cells.Afrotropical/world.table$range.cells
  world.table$austral_WWF_realm.area<-world.table$range.cells.Australasian/world.table$range.cells
  world.table$indo_WWF_realm.area<-world.table$range.cells.IndoMalay/world.table$range.cells
  world.table$nearctic_WWF_realm.area<-world.table$range.cells.Nearctic/world.table$range.cells
  world.table$neotrop_WWF_realm.area<-world.table$range.cells.Neotropical/world.table$range.cells
  world.table$palearctic_WWF_realm.area<-world.table$range.cells.Palearctic/world.table$range.cells
  #define the hotspots by weighted.endemism
  hotspots.table<-read.table(hotspots.file,header=T,sep='\t',stringsAsFactors = F)
  colnames(hotspots.table)<-c('cells','hotspot')
  hotspots<-as.numeric(hotspots.table[hotspots.table$hotspot==1,]$cells)
  #calculate proportion of species range inside hotspots
  grids<-list.files(path=path.grids,pattern='grid_.*_100.rds')
  grid.world<-grids[grep('World_RealmsMerged',grids)]
  grid.realms.names<-grids[-grep('World_RealmsMerged',grids)]
  grid.world<-readRDS(paste(path.grids,'/',grid.world,sep=''))
  grid.realms<-lapply(grid.realms.names,function(x) readRDS(paste(path.grids,'/',x,sep='')))
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
  #make a hotspot list per realm
  hotspots.cells.realm<-list()
  for (i in 1:nrow(grid.cells.df)){
    hotspots.cells.realm[[i]]<-hotspots[intersect(which(hotspots>=grid.cells.df[i,'start.cell']),which(hotspots<=grid.cells.df[i,'end.cell']))]
  }
  #get all the worldwide ranges in a list
  world.table.cells<-as.character(world.table$cells)
  str(unlist(strsplit(world.table.cells[2],' ')))
  world.table.cells<-sapply(world.table.cells,function(x) unlist(strsplit(x,' ')))
  world.table.cells<-sapply(world.table.cells,function(x) x[x != ""])
  world.table.cells<-sapply(world.table.cells,function(x)as.numeric(x))
  #intersect each entry in world.table.cells with all elements in hotspots.cells.realm and record
  world.table.cells.hotspots<-list()
  for (i in 1:length(world.table.cells)){
    world.table.cells.hotspots[[i]]<-unlist(lapply(hotspots.cells.realm,function(x) length(intersect(x,world.table.cells[[i]]))))
  }
  world.table$range.cells.Afrotropical.hotspot<-unlist(lapply(world.table.cells.hotspots,function(x)x[1]))
  world.table$range.cells.Australasian.hotspot<-unlist(lapply(world.table.cells.hotspots,function(x)x[2]))
  world.table$range.cells.IndoMalay.hotspot<-unlist(lapply(world.table.cells.hotspots,function(x)x[3]))
  world.table$range.cells.Nearctic.hotspot<-unlist(lapply(world.table.cells.hotspots,function(x)x[4]))
  world.table$range.cells.Neotropical.hotspot<-unlist(lapply(world.table.cells.hotspots,function(x)x[5]))
  world.table$range.cells.Palearctic.hotspot<-unlist(lapply(world.table.cells.hotspots,function(x)x[6]))
  
  world.table$afrotrop_WWF_hotspot.area<-world.table$range.cells.Afrotropical.hotspot/world.table$range.cells.Afrotropical
  world.table$austral_WWF_hotspot.area<-world.table$range.cells.Australasian.hotspot/world.table$range.cells.Australasian
  world.table$indo_WWF_hotspot.area<-world.table$range.cells.IndoMalay.hotspot/world.table$range.cells.IndoMalay
  world.table$nearctic_WWF_hotspot.area<-world.table$range.cells.Nearctic.hotspot/world.table$range.cells.Nearctic
  world.table$neotrop_WWF_hotspot.area<-world.table$range.cells.Neotropical.hotspot/world.table$range.cells.Neotropical
  world.table$palearctic_WWF_hotspot.area<-world.table$range.cells.Palearctic.hotspot/world.table$range.cells.Palearctic
  world.table[is.na(world.table)] <- 0
  write.table(world.table,file=paste(path,'/100_all_realms_ranges_plus_hotspots_SRlatLM.txt',sep=''),sep='\t',quote=F,row.names=F)
}

build_world_range_hotspots_table_predefined_SR<-function(path,world.table.file,hotspots.file,quantile.hotspot,path.grids){
  list<-list.files(path=path,pattern='100_.*_realms_species_gridoccurrence_table.txt')
  #removes all_realm_file
  if(length(grep('all_realms',list))>0){
    list<-list[-grep('all_realms',list)]
  }
  #remove Oceanic
  if(length(grep('Oceanic',list))>0){
    list<-list[-grep('Oceanic',list)]
  }
  #there should be 6 files (Afrotropical,Australasian,Indo-Malay,Nearctic,Neotropical,Palearctic)
  #if not, exit
  if(!length(list)==6){
    return('not 6 files for 6 realms in this folder!')
  }
  #read region tables and paste them to data frames in a list
  all.regions<-lapply(list, function(x) read.table(paste(path,x,sep=''),header=T,sep='\t'))
  #reorder dataframes by world.table
  world.table<-read.table(world.table.file,header=T,sep='\t')
  all.regions<-lapply(all.regions,function(x) x[match(world.table$spp,x$spp),])
  all.regions.ranges<-lapply(all.regions,function(x)x[,3])
  world.table$range.cells.Afrotropical<-all.regions.ranges[[1]]
  world.table$range.cells.Australasian<-all.regions.ranges[[2]]
  world.table$range.cells.IndoMalay<-all.regions.ranges[[3]]
  world.table$range.cells.Nearctic<-all.regions.ranges[[4]]
  world.table$range.cells.Neotropical<-all.regions.ranges[[5]]
  world.table$range.cells.Palearctic<-all.regions.ranges[[6]]
  #get proportions of range in each realm
  world.table$afrotrop_WWF_realm.area<-world.table$range.cells.Afrotropical/world.table$range.cells
  world.table$austral_WWF_realm.area<-world.table$range.cells.Australasian/world.table$range.cells
  world.table$indo_WWF_realm.area<-world.table$range.cells.IndoMalay/world.table$range.cells
  world.table$nearctic_WWF_realm.area<-world.table$range.cells.Nearctic/world.table$range.cells
  world.table$neotrop_WWF_realm.area<-world.table$range.cells.Neotropical/world.table$range.cells
  world.table$palearctic_WWF_realm.area<-world.table$range.cells.Palearctic/world.table$range.cells
  #define the hotspots by weighted.endemism
  hotspots.table<-read.table(hotspots.file,header=T,sep='\t',stringsAsFactors = F)
  colnames(hotspots.table)<-c('cells','hotspot')
  hotspots<-as.numeric(hotspots.table[hotspots.table$hotspot==1,]$cells)
  #calculate proportion of species range inside hotspots
  grids<-list.files(path=path.grids,pattern='grid_.*_100.rds')
  grid.world<-grids[grep('World_RealmsMerged',grids)]
  grid.realms.names<-grids[-grep('World_RealmsMerged',grids)]
  grid.world<-readRDS(paste(path.grids,'/',grid.world,sep=''))
  grid.realms<-lapply(grid.realms.names,function(x) readRDS(paste(path.grids,'/',x,sep='')))
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
  #make a hotspot list per realm
  hotspots.cells.realm<-list()
  for (i in 1:nrow(grid.cells.df)){
    hotspots.cells.realm[[i]]<-hotspots[intersect(which(hotspots>=grid.cells.df[i,'start.cell']),which(hotspots<=grid.cells.df[i,'end.cell']))]
  }
  #get all the worldwide ranges in a list
  world.table.cells<-as.character(world.table$cells)
  str(unlist(strsplit(world.table.cells[2],' ')))
  world.table.cells<-sapply(world.table.cells,function(x) unlist(strsplit(x,' ')))
  world.table.cells<-sapply(world.table.cells,function(x) x[x != ""])
  world.table.cells<-sapply(world.table.cells,function(x)as.numeric(x))
  #intersect each entry in world.table.cells with all elements in hotspots.cells.realm and record
  world.table.cells.hotspots<-list()
  for (i in 1:length(world.table.cells)){
    world.table.cells.hotspots[[i]]<-unlist(lapply(hotspots.cells.realm,function(x) length(intersect(x,world.table.cells[[i]]))))
  }
  world.table$range.cells.Afrotropical.hotspot<-unlist(lapply(world.table.cells.hotspots,function(x)x[1]))
  world.table$range.cells.Australasian.hotspot<-unlist(lapply(world.table.cells.hotspots,function(x)x[2]))
  world.table$range.cells.IndoMalay.hotspot<-unlist(lapply(world.table.cells.hotspots,function(x)x[3]))
  world.table$range.cells.Nearctic.hotspot<-unlist(lapply(world.table.cells.hotspots,function(x)x[4]))
  world.table$range.cells.Neotropical.hotspot<-unlist(lapply(world.table.cells.hotspots,function(x)x[5]))
  world.table$range.cells.Palearctic.hotspot<-unlist(lapply(world.table.cells.hotspots,function(x)x[6]))
  
  world.table$afrotrop_WWF_hotspot.area<-world.table$range.cells.Afrotropical.hotspot/world.table$range.cells.Afrotropical
  world.table$austral_WWF_hotspot.area<-world.table$range.cells.Australasian.hotspot/world.table$range.cells.Australasian
  world.table$indo_WWF_hotspot.area<-world.table$range.cells.IndoMalay.hotspot/world.table$range.cells.IndoMalay
  world.table$nearctic_WWF_hotspot.area<-world.table$range.cells.Nearctic.hotspot/world.table$range.cells.Nearctic
  world.table$neotrop_WWF_hotspot.area<-world.table$range.cells.Neotropical.hotspot/world.table$range.cells.Neotropical
  world.table$palearctic_WWF_hotspot.area<-world.table$range.cells.Palearctic.hotspot/world.table$range.cells.Palearctic
  world.table[is.na(world.table)] <- 0
  write.table(world.table,file=paste(path,'/100_all_realms_ranges_plus_hotspots_SR.txt',sep=''),sep='\t',quote=F,row.names=F)
}
build_world_range_hotspots_table_predefined_narrow<-function(path,world.table.file,hotspots.file,quantile.hotspot,path.grids){
  list<-list.files(path=path,pattern='100_.*_realms_species_gridoccurrence_table.txt')
  #removes all_realm_file
  if(length(grep('all_realms',list))>0){
    list<-list[-grep('all_realms',list)]
  }
  #remove Oceanic
  if(length(grep('Oceanic',list))>0){
    list<-list[-grep('Oceanic',list)]
  }
  #there should be 6 files (Afrotropical,Australasian,Indo-Malay,Nearctic,Neotropical,Palearctic)
  #if not, exit
  if(!length(list)==6){
    return('not 6 files for 6 realms in this folder!')
  }
  #read region tables and paste them to data frames in a list
  all.regions<-lapply(list, function(x) read.table(paste(path,x,sep=''),header=T,sep='\t'))
  #reorder dataframes by world.table
  world.table<-read.table(world.table.file,header=T,sep='\t')
  all.regions<-lapply(all.regions,function(x) x[match(world.table$spp,x$spp),])
  all.regions.ranges<-lapply(all.regions,function(x)x[,3])
  world.table$range.cells.Afrotropical<-all.regions.ranges[[1]]
  world.table$range.cells.Australasian<-all.regions.ranges[[2]]
  world.table$range.cells.IndoMalay<-all.regions.ranges[[3]]
  world.table$range.cells.Nearctic<-all.regions.ranges[[4]]
  world.table$range.cells.Neotropical<-all.regions.ranges[[5]]
  world.table$range.cells.Palearctic<-all.regions.ranges[[6]]
  #get proportions of range in each realm
  world.table$afrotrop_WWF_realm.area<-world.table$range.cells.Afrotropical/world.table$range.cells
  world.table$austral_WWF_realm.area<-world.table$range.cells.Australasian/world.table$range.cells
  world.table$indo_WWF_realm.area<-world.table$range.cells.IndoMalay/world.table$range.cells
  world.table$nearctic_WWF_realm.area<-world.table$range.cells.Nearctic/world.table$range.cells
  world.table$neotrop_WWF_realm.area<-world.table$range.cells.Neotropical/world.table$range.cells
  world.table$palearctic_WWF_realm.area<-world.table$range.cells.Palearctic/world.table$range.cells
  #define the hotspots by weighted.endemism
  hotspots.table<-read.table(hotspots.file,header=T,sep='\t',stringsAsFactors = F)
  colnames(hotspots.table)<-c('cells','hotspot')
  hotspots<-as.numeric(hotspots.table[hotspots.table$hotspot==1,]$cells)
  #calculate proportion of species range inside hotspots
  grids<-list.files(path=path.grids,pattern='grid_.*_100.rds')
  grid.world<-grids[grep('World_RealmsMerged',grids)]
  grid.realms.names<-grids[-grep('World_RealmsMerged',grids)]
  grid.world<-readRDS(paste(path.grids,'/',grid.world,sep=''))
  grid.realms<-lapply(grid.realms.names,function(x) readRDS(paste(path.grids,'/',x,sep='')))
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
  #make a hotspot list per realm
  hotspots.cells.realm<-list()
  for (i in 1:nrow(grid.cells.df)){
    hotspots.cells.realm[[i]]<-hotspots[intersect(which(hotspots>=grid.cells.df[i,'start.cell']),which(hotspots<=grid.cells.df[i,'end.cell']))]
  }
  #get all the worldwide ranges in a list
  world.table.cells<-as.character(world.table$cells)
  str(unlist(strsplit(world.table.cells[2],' ')))
  world.table.cells<-sapply(world.table.cells,function(x) unlist(strsplit(x,' ')))
  world.table.cells<-sapply(world.table.cells,function(x) x[x != ""])
  world.table.cells<-sapply(world.table.cells,function(x)as.numeric(x))
  #intersect each entry in world.table.cells with all elements in hotspots.cells.realm and record
  world.table.cells.hotspots<-list()
  for (i in 1:length(world.table.cells)){
    world.table.cells.hotspots[[i]]<-unlist(lapply(hotspots.cells.realm,function(x) length(intersect(x,world.table.cells[[i]]))))
  }
  world.table$range.cells.Afrotropical.hotspot<-unlist(lapply(world.table.cells.hotspots,function(x)x[1]))
  world.table$range.cells.Australasian.hotspot<-unlist(lapply(world.table.cells.hotspots,function(x)x[2]))
  world.table$range.cells.IndoMalay.hotspot<-unlist(lapply(world.table.cells.hotspots,function(x)x[3]))
  world.table$range.cells.Nearctic.hotspot<-unlist(lapply(world.table.cells.hotspots,function(x)x[4]))
  world.table$range.cells.Neotropical.hotspot<-unlist(lapply(world.table.cells.hotspots,function(x)x[5]))
  world.table$range.cells.Palearctic.hotspot<-unlist(lapply(world.table.cells.hotspots,function(x)x[6]))
  
  world.table$afrotrop_WWF_hotspot.area<-world.table$range.cells.Afrotropical.hotspot/world.table$range.cells.Afrotropical
  world.table$austral_WWF_hotspot.area<-world.table$range.cells.Australasian.hotspot/world.table$range.cells.Australasian
  world.table$indo_WWF_hotspot.area<-world.table$range.cells.IndoMalay.hotspot/world.table$range.cells.IndoMalay
  world.table$nearctic_WWF_hotspot.area<-world.table$range.cells.Nearctic.hotspot/world.table$range.cells.Nearctic
  world.table$neotrop_WWF_hotspot.area<-world.table$range.cells.Neotropical.hotspot/world.table$range.cells.Neotropical
  world.table$palearctic_WWF_hotspot.area<-world.table$range.cells.Palearctic.hotspot/world.table$range.cells.Palearctic
  world.table[is.na(world.table)] <- 0
  write.table(world.table,file=paste(path,'/100_all_realms_ranges_plus_hotspots_narrowrangesp.txt',sep=''),sep='\t',quote=F,row.names=F)
}


build_world_range_hotspots_table_predefined_SR_WEsize<-function(path,world.table.file,hotspots.file,quantile.hotspot,path.grids){
  list<-list.files(path=path,pattern='100_.*_realms_species_gridoccurrence_table.txt')
  #removes all_realm_file
  if(length(grep('all_realms',list))>0){
    list<-list[-grep('all_realms',list)]
  }
  #remove Oceanic
  if(length(grep('Oceanic',list))>0){
    list<-list[-grep('Oceanic',list)]
  }
  #there should be 6 files (Afrotropical,Australasian,Indo-Malay,Nearctic,Neotropical,Palearctic)
  #if not, exit
  if(!length(list)==6){
    return('not 6 files for 6 realms in this folder!')
  }
  #read region tables and paste them to data frames in a list
  all.regions<-lapply(list, function(x) read.table(paste(path,x,sep=''),header=T,sep='\t'))
  #reorder dataframes by world.table
  world.table<-read.table(world.table.file,header=T,sep='\t')
  all.regions<-lapply(all.regions,function(x) x[match(world.table$spp,x$spp),])
  all.regions.ranges<-lapply(all.regions,function(x)x[,3])
  world.table$range.cells.Afrotropical<-all.regions.ranges[[1]]
  world.table$range.cells.Australasian<-all.regions.ranges[[2]]
  world.table$range.cells.IndoMalay<-all.regions.ranges[[3]]
  world.table$range.cells.Nearctic<-all.regions.ranges[[4]]
  world.table$range.cells.Neotropical<-all.regions.ranges[[5]]
  world.table$range.cells.Palearctic<-all.regions.ranges[[6]]
  #get proportions of range in each realm
  world.table$afrotrop_WWF_realm.area<-world.table$range.cells.Afrotropical/world.table$range.cells
  world.table$austral_WWF_realm.area<-world.table$range.cells.Australasian/world.table$range.cells
  world.table$indo_WWF_realm.area<-world.table$range.cells.IndoMalay/world.table$range.cells
  world.table$nearctic_WWF_realm.area<-world.table$range.cells.Nearctic/world.table$range.cells
  world.table$neotrop_WWF_realm.area<-world.table$range.cells.Neotropical/world.table$range.cells
  world.table$palearctic_WWF_realm.area<-world.table$range.cells.Palearctic/world.table$range.cells
  #define the hotspots by weighted.endemism
  hotspots.table<-read.table(hotspots.file,header=T,sep='\t',stringsAsFactors = F)
  colnames(hotspots.table)<-c('cells','hotspot')
  hotspots<-as.numeric(hotspots.table[hotspots.table$hotspot==1,]$cells)
  #calculate proportion of species range inside hotspots
  grids<-list.files(path=path.grids,pattern='grid_.*_100.rds')
  grid.world<-grids[grep('World_RealmsMerged',grids)]
  grid.realms.names<-grids[-grep('World_RealmsMerged',grids)]
  grid.world<-readRDS(paste(path.grids,'/',grid.world,sep=''))
  grid.realms<-lapply(grid.realms.names,function(x) readRDS(paste(path.grids,'/',x,sep='')))
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
  #make a hotspot list per realm
  hotspots.cells.realm<-list()
  for (i in 1:nrow(grid.cells.df)){
    hotspots.cells.realm[[i]]<-hotspots[intersect(which(hotspots>=grid.cells.df[i,'start.cell']),which(hotspots<=grid.cells.df[i,'end.cell']))]
  }
  #get all the worldwide ranges in a list
  world.table.cells<-as.character(world.table$cells)
  str(unlist(strsplit(world.table.cells[2],' ')))
  world.table.cells<-sapply(world.table.cells,function(x) unlist(strsplit(x,' ')))
  world.table.cells<-sapply(world.table.cells,function(x) x[x != ""])
  world.table.cells<-sapply(world.table.cells,function(x)as.numeric(x))
  #intersect each entry in world.table.cells with all elements in hotspots.cells.realm and record
  world.table.cells.hotspots<-list()
  for (i in 1:length(world.table.cells)){
    world.table.cells.hotspots[[i]]<-unlist(lapply(hotspots.cells.realm,function(x) length(intersect(x,world.table.cells[[i]]))))
  }
  world.table$range.cells.Afrotropical.hotspot<-unlist(lapply(world.table.cells.hotspots,function(x)x[1]))
  world.table$range.cells.Australasian.hotspot<-unlist(lapply(world.table.cells.hotspots,function(x)x[2]))
  world.table$range.cells.IndoMalay.hotspot<-unlist(lapply(world.table.cells.hotspots,function(x)x[3]))
  world.table$range.cells.Nearctic.hotspot<-unlist(lapply(world.table.cells.hotspots,function(x)x[4]))
  world.table$range.cells.Neotropical.hotspot<-unlist(lapply(world.table.cells.hotspots,function(x)x[5]))
  world.table$range.cells.Palearctic.hotspot<-unlist(lapply(world.table.cells.hotspots,function(x)x[6]))
  
  world.table$afrotrop_WWF_hotspot.area<-world.table$range.cells.Afrotropical.hotspot/world.table$range.cells.Afrotropical
  world.table$austral_WWF_hotspot.area<-world.table$range.cells.Australasian.hotspot/world.table$range.cells.Australasian
  world.table$indo_WWF_hotspot.area<-world.table$range.cells.IndoMalay.hotspot/world.table$range.cells.IndoMalay
  world.table$nearctic_WWF_hotspot.area<-world.table$range.cells.Nearctic.hotspot/world.table$range.cells.Nearctic
  world.table$neotrop_WWF_hotspot.area<-world.table$range.cells.Neotropical.hotspot/world.table$range.cells.Neotropical
  world.table$palearctic_WWF_hotspot.area<-world.table$range.cells.Palearctic.hotspot/world.table$range.cells.Palearctic
  world.table[is.na(world.table)] <- 0
  write.table(world.table,file=paste(path,'/100_all_realms_ranges_plus_hotspots_SR_WEsize.txt',sep=''),sep='\t',quote=F,row.names=F)
}
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

build_world_range_hotspots_table_SR<-function(path,world.table.file,weighted.endemism.file,quantile.hotspot,path.grids){
  list<-list.files(path=path,pattern='100_.*_realms_species_gridoccurrence_table.txt')
  #removes all_realm_file
  if(length(grep('all_realms',list))>0){
    list<-list[-grep('all_realms',list)]
  }
  #remove Oceanic
  if(length(grep('Oceanic',list))>0){
    list<-list[-grep('Oceanic',list)]
  }
  #there should be 6 files (Afrotropical,Australasian,Indo-Malay,Nearctic,Neotropical,Palearctic)
  #if not, exit
  if(!length(list)==6){
    return('not 6 files for 6 realms in this folder!')
  }
  #read region tables and paste them to data frames in a list
  all.regions<-lapply(list, function(x) read.table(paste(path,x,sep=''),header=T,sep='\t'))
  #reorder dataframes by world.table
  world.table<-read.table(world.table.file,header=T,sep='\t')
  all.regions<-lapply(all.regions,function(x) x[match(world.table$spp,x$spp),])
  all.regions.ranges<-lapply(all.regions,function(x)x[,3])
  world.table$range.cells.Afrotropical<-all.regions.ranges[[1]]
  world.table$range.cells.Australasian<-all.regions.ranges[[2]]
  world.table$range.cells.IndoMalay<-all.regions.ranges[[3]]
  world.table$range.cells.Nearctic<-all.regions.ranges[[4]]
  world.table$range.cells.Neotropical<-all.regions.ranges[[5]]
  world.table$range.cells.Palearctic<-all.regions.ranges[[6]]
  #get proportions of range in each realm
  world.table$afrotrop_WWF_realm.area<-world.table$range.cells.Afrotropical/world.table$range.cells
  world.table$austral_WWF_realm.area<-world.table$range.cells.Australasian/world.table$range.cells
  world.table$indo_WWF_realm.area<-world.table$range.cells.IndoMalay/world.table$range.cells
  world.table$nearctic_WWF_realm.area<-world.table$range.cells.Nearctic/world.table$range.cells
  world.table$neotrop_WWF_realm.area<-world.table$range.cells.Neotropical/world.table$range.cells
  world.table$palearctic_WWF_realm.area<-world.table$range.cells.Palearctic/world.table$range.cells
  #define the hotspots by weighted.endemism
  wend<-read.table(weighted.endemism.file,header=T,sep='\t',stringsAsFactors = F)
  wend<-na.omit(wend)
  colnames(wend)[1]<-'cells'
  quantile<-quantile.hotspot
  hotspots<-wend[wend$number.of.species>quantile(wend$number.of.species,quantile),]$cells
  #calculate proportion of species range inside hotspots
  grids<-list.files(path=path.grids,pattern='grid_.*_100.rds')
  grid.world<-grids[grep('World_RealmsMerged',grids)]
  grid.realms.names<-grids[-grep('World_RealmsMerged',grids)]
  grid.world<-readRDS(paste(path.grids,'/',grid.world,sep=''))
  grid.realms<-lapply(grid.realms.names,function(x) readRDS(paste(path.grids,'/',x,sep='')))
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
  #make a hotspot list per realm
  hotspots.cells.realm<-list()
  for (i in 1:nrow(grid.cells.df)){
    hotspots.cells.realm[[i]]<-hotspots[intersect(which(hotspots>=grid.cells.df[i,'start.cell']),which(hotspots<=grid.cells.df[i,'end.cell']))]
  }
  #get all the worldwide ranges in a list
  world.table.cells<-as.character(world.table$cells)
  str(unlist(strsplit(world.table.cells[2],' ')))
  world.table.cells<-sapply(world.table.cells,function(x) unlist(strsplit(x,' ')))
  world.table.cells<-sapply(world.table.cells,function(x) x[x != ""])
  world.table.cells<-sapply(world.table.cells,function(x)as.numeric(x))
  #intersect each entry in world.table.cells with all elements in hotspots.cells.realm and record
  world.table.cells.hotspots<-list()
  for (i in 1:length(world.table.cells)){
    world.table.cells.hotspots[[i]]<-unlist(lapply(hotspots.cells.realm,function(x) length(intersect(x,world.table.cells[[i]]))))
  }
  world.table$range.cells.Afrotropical.hotspot<-unlist(lapply(world.table.cells.hotspots,function(x)x[1]))
  world.table$range.cells.Australasian.hotspot<-unlist(lapply(world.table.cells.hotspots,function(x)x[2]))
  world.table$range.cells.IndoMalay.hotspot<-unlist(lapply(world.table.cells.hotspots,function(x)x[3]))
  world.table$range.cells.Nearctic.hotspot<-unlist(lapply(world.table.cells.hotspots,function(x)x[4]))
  world.table$range.cells.Neotropical.hotspot<-unlist(lapply(world.table.cells.hotspots,function(x)x[5]))
  world.table$range.cells.Palearctic.hotspot<-unlist(lapply(world.table.cells.hotspots,function(x)x[6]))
  
  world.table$afrotrop_WWF_hotspot.area<-world.table$range.cells.Afrotropical.hotspot/world.table$range.cells.Afrotropical
  world.table$austral_WWF_hotspot.area<-world.table$range.cells.Australasian.hotspot/world.table$range.cells.Australasian
  world.table$indo_WWF_hotspot.area<-world.table$range.cells.IndoMalay.hotspot/world.table$range.cells.IndoMalay
  world.table$nearctic_WWF_hotspot.area<-world.table$range.cells.Nearctic.hotspot/world.table$range.cells.Nearctic
  world.table$neotrop_WWF_hotspot.area<-world.table$range.cells.Neotropical.hotspot/world.table$range.cells.Neotropical
  world.table$palearctic_WWF_hotspot.area<-world.table$range.cells.Palearctic.hotspot/world.table$range.cells.Palearctic
  world.table[is.na(world.table)] <- 0
  write.table(world.table,file=paste(path,'/100_all_realms_ranges_plus_speciesrichnesshotspots.txt',sep=''),sep='\t',quote=F,row.names=F)
}

prepare_realm_input_7areas_plus_inhotoutrealm_SR<-function(table,name,overlap,tree){
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
  realm.character<-as.character(grep(name,colnames(realm.table)[-c(1,8)]))
  realm.hotspot.table<-table[,grep(name,colnames(table))]
  realm.hotspot.table<-realm.hotspot.table[which(realm.hotspot.table[,1]>=overlap),]
  realm.hotspot.table$model.state<-0
  realm.hotspot.table[realm.hotspot.table[,2]<(1-overlap),'model.state']<-realm.character
  realm.hotspot.table[realm.hotspot.table[,2]>=overlap,'model.state']<-7
  realm.hotspot.table[realm.hotspot.table[,2]>(1-overlap)&realm.hotspot.table[,2]<overlap,'model.state']<-paste(realm.character,7,sep=',')
  realm.table[realm.table$model.state==realm.character,'model.state']<-realm.hotspot.table$model.state
  #select species occurring in target realm plus another realm to check if  
  target.states<-names(table(realm.table$model.state))[grep(realm.character,names(table(realm.table$model.state)))]
  target.states<-target.states[-grep('7',target.states)]
  target.states<-target.states[-match(realm.character,target.states)]
  realm.plus.outside.table<-realm.table[realm.table$model.state%in%target.states,]
  #add hotspot overlap within realm to realm.plus.outside table
  realm.plus.outside.table.hotspots<-as.data.frame(cbind(table[table$spp%in%realm.plus.outside.table$spp,'spp'],table[table$spp%in%realm.plus.outside.table$spp,paste(name,'_hotspot.area',sep='')]),stringsAsFactors = F)
  colnames(realm.plus.outside.table.hotspots)<-c('spp',paste(name,'_hotspot.area',sep=''))
  realm.plus.outside.table.hotspots$model.state<-realm.table[realm.table$model.state%in%target.states,'model.state']
  realm.plus.outside.table.hotspots[realm.plus.outside.table.hotspots[,2]>=overlap,'model.state']<-sub(realm.plus.outside.table.hotspots[realm.plus.outside.table.hotspots[,2]>=overlap,'model.state'],pattern=realm.character,replacement='7')
  realm.plus.outside.table.hotspots[realm.plus.outside.table.hotspots$model.state==paste('7,','1',sep=''),'model.state']<-'1,7'
  realm.plus.outside.table.hotspots[realm.plus.outside.table.hotspots$model.state==paste('7,','2',sep=''),'model.state']<-'2,7'
  realm.plus.outside.table.hotspots[realm.plus.outside.table.hotspots$model.state==paste('7,','3',sep=''),'model.state']<-'3,7'
  realm.plus.outside.table.hotspots[realm.plus.outside.table.hotspots$model.state==paste('7,','4',sep=''),'model.state']<-'4,7'
  realm.plus.outside.table.hotspots[realm.plus.outside.table.hotspots$model.state==paste('7,','5',sep=''),'model.state']<-'5,7'
  realm.plus.outside.table.hotspots[realm.plus.outside.table.hotspots$model.state==paste('7,','6',sep=''),'model.state']<-'6,7'
  realm.table[realm.table$model.state%in%target.states,'model.state']<-realm.plus.outside.table.hotspots$model.state
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
  write.table(header,paste('./',name,'_7areas_inhotoutrealm_geographyfile_SRhotspot.txt',sep=''),sep='\t',quote=F,row.names=F,col.names = F)
  write.table(geographytable,paste('./',name,'_7areas_inhotoutrealm_geographyfile_SRhotspot.txt',sep=''),append=T,sep=' ',quote=F,row.names=F,col.names = F)
  return(table.model.state)
  
}

prepare_realm_input_7areas_plus_inhotoutrealm_lmSRlat<-function(table,name,overlap,tree){
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
  realm.character<-as.character(grep(name,colnames(realm.table)[-c(1,8)]))
  realm.hotspot.table<-table[,grep(name,colnames(table))]
  realm.hotspot.table<-realm.hotspot.table[which(realm.hotspot.table[,1]>=overlap),]
  realm.hotspot.table$model.state<-0
  realm.hotspot.table[realm.hotspot.table[,2]<(1-overlap),'model.state']<-realm.character
  realm.hotspot.table[realm.hotspot.table[,2]>=overlap,'model.state']<-7
  realm.hotspot.table[realm.hotspot.table[,2]>(1-overlap)&realm.hotspot.table[,2]<overlap,'model.state']<-paste(realm.character,7,sep=',')
  realm.table[realm.table$model.state==realm.character,'model.state']<-realm.hotspot.table$model.state
  #select species occurring in target realm plus another realm to check if  
  target.states<-names(table(realm.table$model.state))[grep(realm.character,names(table(realm.table$model.state)))]
  target.states<-target.states[-grep('7',target.states)]
  target.states<-target.states[-match(realm.character,target.states)]
  realm.plus.outside.table<-realm.table[realm.table$model.state%in%target.states,]
  #add hotspot overlap within realm to realm.plus.outside table
  realm.plus.outside.table.hotspots<-as.data.frame(cbind(table[table$spp%in%realm.plus.outside.table$spp,'spp'],table[table$spp%in%realm.plus.outside.table$spp,paste(name,'_hotspot.area',sep='')]),stringsAsFactors = F)
  colnames(realm.plus.outside.table.hotspots)<-c('spp',paste(name,'_hotspot.area',sep=''))
  realm.plus.outside.table.hotspots$model.state<-realm.table[realm.table$model.state%in%target.states,'model.state']
  realm.plus.outside.table.hotspots[realm.plus.outside.table.hotspots[,2]>=overlap,'model.state']<-sub(realm.plus.outside.table.hotspots[realm.plus.outside.table.hotspots[,2]>=overlap,'model.state'],pattern=realm.character,replacement='7')
  realm.plus.outside.table.hotspots[realm.plus.outside.table.hotspots$model.state==paste('7,','1',sep=''),'model.state']<-'1,7'
  realm.plus.outside.table.hotspots[realm.plus.outside.table.hotspots$model.state==paste('7,','2',sep=''),'model.state']<-'2,7'
  realm.plus.outside.table.hotspots[realm.plus.outside.table.hotspots$model.state==paste('7,','3',sep=''),'model.state']<-'3,7'
  realm.plus.outside.table.hotspots[realm.plus.outside.table.hotspots$model.state==paste('7,','4',sep=''),'model.state']<-'4,7'
  realm.plus.outside.table.hotspots[realm.plus.outside.table.hotspots$model.state==paste('7,','5',sep=''),'model.state']<-'5,7'
  realm.plus.outside.table.hotspots[realm.plus.outside.table.hotspots$model.state==paste('7,','6',sep=''),'model.state']<-'6,7'
  realm.table[realm.table$model.state%in%target.states,'model.state']<-realm.plus.outside.table.hotspots$model.state
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
  write.table(header,paste('./',name,'_7areas_inhotoutrealm_geographyfile_lmSRlat_hotspot.txt',sep=''),sep='\t',quote=F,row.names=F,col.names = F)
  write.table(geographytable,paste('./',name,'_7areas_inhotoutrealm_geographyfile_lmSRlat_hotspot.txt',sep=''),append=T,sep=' ',quote=F,row.names=F,col.names = F)
  return(table.model.state)
  
}
prepare_realm_input_7areas_plus_inhotoutrealm_narrow<-function(table,name,overlap,tree){
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
  realm.character<-as.character(grep(name,colnames(realm.table)[-c(1,8)]))
  realm.hotspot.table<-table[,grep(name,colnames(table))]
  realm.hotspot.table<-realm.hotspot.table[which(realm.hotspot.table[,1]>=overlap),]
  realm.hotspot.table$model.state<-0
  realm.hotspot.table[realm.hotspot.table[,2]<(1-overlap),'model.state']<-realm.character
  realm.hotspot.table[realm.hotspot.table[,2]>=overlap,'model.state']<-7
  realm.hotspot.table[realm.hotspot.table[,2]>(1-overlap)&realm.hotspot.table[,2]<overlap,'model.state']<-paste(realm.character,7,sep=',')
  realm.table[realm.table$model.state==realm.character,'model.state']<-realm.hotspot.table$model.state
  #select species occurring in target realm plus another realm to check if  
  target.states<-names(table(realm.table$model.state))[grep(realm.character,names(table(realm.table$model.state)))]
  target.states<-target.states[-grep('7',target.states)]
  target.states<-target.states[-match(realm.character,target.states)]
  realm.plus.outside.table<-realm.table[realm.table$model.state%in%target.states,]
  #add hotspot overlap within realm to realm.plus.outside table
  realm.plus.outside.table.hotspots<-as.data.frame(cbind(table[table$spp%in%realm.plus.outside.table$spp,'spp'],table[table$spp%in%realm.plus.outside.table$spp,paste(name,'_hotspot.area',sep='')]),stringsAsFactors = F)
  colnames(realm.plus.outside.table.hotspots)<-c('spp',paste(name,'_hotspot.area',sep=''))
  realm.plus.outside.table.hotspots$model.state<-realm.table[realm.table$model.state%in%target.states,'model.state']
  realm.plus.outside.table.hotspots[realm.plus.outside.table.hotspots[,2]>=overlap,'model.state']<-sub(realm.plus.outside.table.hotspots[realm.plus.outside.table.hotspots[,2]>=overlap,'model.state'],pattern=realm.character,replacement='7')
  realm.plus.outside.table.hotspots[realm.plus.outside.table.hotspots$model.state==paste('7,','1',sep=''),'model.state']<-'1,7'
  realm.plus.outside.table.hotspots[realm.plus.outside.table.hotspots$model.state==paste('7,','2',sep=''),'model.state']<-'2,7'
  realm.plus.outside.table.hotspots[realm.plus.outside.table.hotspots$model.state==paste('7,','3',sep=''),'model.state']<-'3,7'
  realm.plus.outside.table.hotspots[realm.plus.outside.table.hotspots$model.state==paste('7,','4',sep=''),'model.state']<-'4,7'
  realm.plus.outside.table.hotspots[realm.plus.outside.table.hotspots$model.state==paste('7,','5',sep=''),'model.state']<-'5,7'
  realm.plus.outside.table.hotspots[realm.plus.outside.table.hotspots$model.state==paste('7,','6',sep=''),'model.state']<-'6,7'
  realm.table[realm.table$model.state%in%target.states,'model.state']<-realm.plus.outside.table.hotspots$model.state
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
  write.table(header,paste('./',name,'_7areas_inhotoutrealm_geographyfile_narrowhotspot.txt',sep=''),sep='\t',quote=F,row.names=F,col.names = F)
  write.table(geographytable,paste('./',name,'_7areas_inhotoutrealm_geographyfile_narrowhotspot.txt',sep=''),append=T,sep=' ',quote=F,row.names=F,col.names = F)
  return(table.model.state)
  
}

prepare_realm_input_7areas_plus_inhotoutrealm_SR_WEhotsize<-function(table,name,overlap,tree){
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
  realm.table[realm.table$model.state==0,]$model.state<-multirrealm.species.regions
  #drop species with 0 = species occurring in Oceanic realm
  realm.table<-realm.table[!realm.table$model.state==0,]
  #select column of realm to be analysed
  realm.character<-as.character(grep(name,colnames(realm.table)[-c(1,8)]))
  realm.hotspot.table<-table[,grep(name,colnames(table))]
  realm.hotspot.table<-realm.hotspot.table[which(realm.hotspot.table[,1]>=overlap),]
  realm.hotspot.table$model.state<-0
  realm.hotspot.table[realm.hotspot.table[,2]<(1-overlap),'model.state']<-realm.character
  realm.hotspot.table[realm.hotspot.table[,2]>=overlap,'model.state']<-7
  realm.hotspot.table[realm.hotspot.table[,2]>(1-overlap)&realm.hotspot.table[,2]<overlap,'model.state']<-paste(realm.character,7,sep=',')
  realm.table[realm.table$model.state==realm.character,'model.state']<-realm.hotspot.table$model.state
  #select species occurring in target realm plus another realm to check if  
  target.states<-names(table(realm.table$model.state))[grep(realm.character,names(table(realm.table$model.state)))]
  target.states<-target.states[-grep('7',target.states)]
  target.states<-target.states[-match(realm.character,target.states)]
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
  realm.table[realm.table$model.state%in%target.states,]$model.state<-realm.plus.outside.table.hotspots$model.state
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
  write.table(header,paste('./',name,'_7areas_inhotoutrealm_geographyfile_SR_WEhotsize_hotspot.txt',sep=''),sep='\t',quote=F,row.names=F,col.names = F)
  write.table(geographytable,paste('./',name,'_7areas_inhotoutrealm_geographyfile_SR_WEhotsize_hotspot.txt',sep=''),append=T,sep=' ',quote=F,row.names=F,col.names = F)
  return(table.model.state)
  
}

prepare_realm_input_7areas_plus_inhotoutrealm_SRsarlm<-function(table,name,overlap,tree){
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
  realm.table[realm.table$model.state==0,]$model.state<-multirrealm.species.regions
  #drop species with 0 = species occurring in Oceanic realm
  realm.table<-realm.table[!realm.table$model.state==0,]
  #select column of realm to be analysed
  realm.character<-as.character(grep(name,colnames(realm.table)[-c(1,8)]))
  realm.hotspot.table<-table[,grep(name,colnames(table))]
  realm.hotspot.table<-realm.hotspot.table[which(realm.hotspot.table[,1]>=overlap),]
  realm.hotspot.table$model.state<-0
  realm.hotspot.table[realm.hotspot.table[,2]<(1-overlap),'model.state']<-realm.character
  realm.hotspot.table[realm.hotspot.table[,2]>=overlap,'model.state']<-7
  realm.hotspot.table[realm.hotspot.table[,2]>(1-overlap)&realm.hotspot.table[,2]<overlap,'model.state']<-paste(realm.character,7,sep=',')
  realm.table[realm.table$model.state==realm.character,'model.state']<-realm.hotspot.table$model.state
  #select species occurring in target realm plus another realm to check if  
  target.states<-names(table(realm.table$model.state))[grep(realm.character,names(table(realm.table$model.state)))]
  target.states<-target.states[-grep('7',target.states)]
  target.states<-target.states[-match(realm.character,target.states)]
  realm.plus.outside.table<-realm.table[realm.table$model.state%in%target.states,]
  #add hotspot overlap within realm to realm.plus.outside table
  realm.plus.outside.table.hotspots<-as.data.frame(cbind(table[table$spp%in%realm.plus.outside.table$spp,'spp'],table[table$spp%in%realm.plus.outside.table$spp,paste(name,'_hotspot.area',sep='')]),stringsAsFactors = F)
  colnames(realm.plus.outside.table.hotspots)<-c('spp',paste(name,'_hotspot.area',sep=''))
  realm.plus.outside.table.hotspots$model.state<-realm.table[realm.table$model.state%in%target.states,'model.state']
  realm.plus.outside.table.hotspots[realm.plus.outside.table.hotspots[,2]>=overlap,'model.state']<-sub(realm.plus.outside.table.hotspots[realm.plus.outside.table.hotspots[,2]>=overlap,'model.state'],pattern=realm.character,replacement='7')
  realm.plus.outside.table.hotspots[realm.plus.outside.table.hotspots$model.state==paste('7,','1',sep=''),'model.state']<-'1,7'
  realm.plus.outside.table.hotspots[realm.plus.outside.table.hotspots$model.state==paste('7,','2',sep=''),'model.state']<-'2,7'
  realm.plus.outside.table.hotspots[realm.plus.outside.table.hotspots$model.state==paste('7,','3',sep=''),'model.state']<-'3,7'
  realm.plus.outside.table.hotspots[realm.plus.outside.table.hotspots$model.state==paste('7,','4',sep=''),'model.state']<-'4,7'
  realm.plus.outside.table.hotspots[realm.plus.outside.table.hotspots$model.state==paste('7,','5',sep=''),'model.state']<-'5,7'
  realm.plus.outside.table.hotspots[realm.plus.outside.table.hotspots$model.state==paste('7,','6',sep=''),'model.state']<-'6,7'
  realm.table[realm.table$model.state%in%target.states,'model.state']<-realm.plus.outside.table.hotspots$model.state
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
  write.table(header,paste('./',name,'_7areas_inhotoutrealm_geographyfile_SRsarlmhotspot.txt',sep=''),sep='\t',quote=F,row.names=F,col.names = F)
  write.table(geographytable,paste('./',name,'_7areas_inhotoutrealm_geographyfile_SRsarlmhotspot.txt',sep=''),append=T,sep=' ',quote=F,row.names=F,col.names = F)
  return(table.model.state)
  
}


#world.table.file is 100_all_realms_ranges_plus_hotspots.txt
#treefile is the tree
#name has to be one of (afrotrop,austral,indo,nearctic,neotrop,palearctic)
#overlap is 0.8
prepare_realm_input_7areas_SR<-function(world.table.file,treefile,name,overlap){
  tree<-read.tree(treefile)
  world.table<-read.table(world.table.file,header=T,stringsAsFactors = F,sep='\t')
  #collate tree and table
  tree<-drop.tip(tree,setdiff(tree$tip.label,world.table$spp))
  world.table<-world.table[world.table$spp%in%tree$tip.label,]
  #keep columns with proportion of range inside realm and proportion of range within hotspot
  drop.columns<-grep('cells',colnames(world.table))
  world.table<-world.table[,-drop.columns]
  colnames(world.table)<-sub(colnames(world.table),pattern='_WWF',replacement='')
  prepare_realm_input_7areas_plus_inhotoutrealm_SR(table=world.table,name=name,overlap=0.80,tree=tree)
}
prepare_realm_input_7areas_narrow<-function(world.table.file,treefile,name,overlap){
  tree<-read.tree(treefile)
  world.table<-read.table(world.table.file,header=T,stringsAsFactors = F,sep='\t')
  #collate tree and table
  tree<-drop.tip(tree,setdiff(tree$tip.label,world.table$spp))
  world.table<-world.table[world.table$spp%in%tree$tip.label,]
  #keep columns with proportion of range inside realm and proportion of range within hotspot
  drop.columns<-grep('cells',colnames(world.table))
  world.table<-world.table[,-drop.columns]
  colnames(world.table)<-sub(colnames(world.table),pattern='_WWF',replacement='')
  prepare_realm_input_7areas_plus_inhotoutrealm_narrow(table=world.table,name=name,overlap=0.80,tree=tree)
}

prepare_realm_input_7areas_SRsarlm<-function(world.table.file,treefile,name,overlap){
  tree<-read.tree(treefile)
  world.table<-read.table(world.table.file,header=T,stringsAsFactors = F,sep='\t')
  #collate tree and table
  tree<-drop.tip(tree,setdiff(tree$tip.label,world.table$spp))
  world.table<-world.table[world.table$spp%in%tree$tip.label,]
  #keep columns with proportion of range inside realm and proportion of range within hotspot
  drop.columns<-grep('cells',colnames(world.table))
  world.table<-world.table[,-drop.columns]
  colnames(world.table)<-sub(colnames(world.table),pattern='_WWF',replacement='')
  prepare_realm_input_7areas_plus_inhotoutrealm_SRsarlm(table=world.table,name=name,overlap=0.80,tree=tree)
}

prepare_realm_input_7areas_lmSRlat<-function(world.table.file,treefile,name,overlap){
  tree<-read.tree(treefile)
  world.table<-read.table(world.table.file,header=T,stringsAsFactors = F,sep='\t')
  #collate tree and table
  tree<-drop.tip(tree,setdiff(tree$tip.label,world.table$spp))
  world.table<-world.table[world.table$spp%in%tree$tip.label,]
  #keep columns with proportion of range inside realm and proportion of range within hotspot
  drop.columns<-grep('cells',colnames(world.table))
  world.table<-world.table[,-drop.columns]
  colnames(world.table)<-sub(colnames(world.table),pattern='_WWF',replacement='')
  prepare_realm_input_7areas_plus_inhotoutrealm_lmSRlat(table=world.table,name=name,overlap=0.80,tree=tree)
}
prepare_realm_input_7areas_SR_WEhotsize<-function(world.table.file,treefile,name,overlap){
  tree<-read.tree(treefile)
  world.table<-read.table(world.table.file,header=T,stringsAsFactors = F,sep='\t')
  #collate tree and table
  tree<-drop.tip(tree,setdiff(tree$tip.label,world.table$spp))
  world.table<-world.table[world.table$spp%in%tree$tip.label,]
  #keep columns with proportion of range inside realm and proportion of range within hotspot
  drop.columns<-grep('cells',colnames(world.table))
  world.table<-world.table[,-drop.columns]
  colnames(world.table)<-sub(colnames(world.table),pattern='_WWF',replacement='')
  prepare_realm_input_7areas_plus_inhotoutrealm_SR_WEhotsize(table=world.table,name=name,overlap=0.80,tree=tree)
}




##this is for 5 states
#afro.real<-prepare_realm_simple_input_5states(table=world.table,name='afrotrop',overlap=0.80)
#austral.real<-prepare_realm_simple_input_5states(table=world.table,name='austral',overlap=0.80)
#indo.real<-prepare_realm_simple_input_5states(table=world.table,name='indo',overlap=0.80)
#nearctic.real<-prepare_realm_simple_input_5states(table=world.table,name='nearctic',overlap=0.80)
#neotrop.real<-prepare_realm_simple_input_5states(table=world.table,name='neotrop',overlap=0.80)
#palearctic.real<-prepare_realm_simple_input_5states(table=world.table,name='palearctic',overlap=0.80)
#
##this is for 6 states (state numbers are different from 5state runs!!)
#afro.real6<-prepare_realm_simple_input_6states(table=world.table,name='afrotrop',overlap=0.80)
#austral.real6<-prepare_realm_simple_input_6states(table=world.table,name='austral',overlap=0.80)
#indo.real6<-prepare_realm_simple_input_6states(table=world.table,name='indo',overlap=0.80)
#nearctic.real6<-prepare_realm_simple_input_6states(table=world.table,name='nearctic',overlap=0.80)
#neotrop.real6<-prepare_realm_simple_input_6states(table=world.table,name='neotrop',overlap=0.80)
#palearctic.real6<-prepare_realm_simple_input_6states(table=world.table,name='palearctic',overlap=0.80)
#
##this is for 7 areas (all realms + plus hotspot in each case)
#afro.7areas<-prepare_realm_input_7areas(table=world.table,name='afrotrop',overlap=0.80)
#austral.7areas<-prepare_realm_input_7areas(table=world.table,name='austral',overlap=0.80)
#indo.7areas<-prepare_realm_input_7areas(table=world.table,name='indo',overlap=0.80)
#nearctic.7areas<-prepare_realm_input_7areas(table=world.table,name='nearctic',overlap=0.80)
#neotrop.7areas<-prepare_realm_input_7areas(table=world.table,name='neotrop',overlap=0.80)
#palearctic.7areas<-prepare_realm_input_7areas(table=world.table,name='palearctic',overlap=0.80)
#
##this is for 7 areas (all realms + plus hotspot in each case)
#afro.7areas.inhot_outrealm<-prepare_realm_input_7areas_plus_inhotoutrealm(table=world.table,name='afrotrop',overlap=0.80)
#austral.7areas.inhot_outrealm<-prepare_realm_input_7areas_plus_inhotoutrealm(table=world.table,name='austral',overlap=0.80)
#indo.7areas.inhot_outrealm<-prepare_realm_input_7areas_plus_inhotoutrealm(table=world.table,name='indo',overlap=0.80)
#nearctic.7areas.inhot_outrealm<-prepare_realm_input_7areas_plus_inhotoutrealm(table=world.table,name='nearctic',overlap=0.80)
#neotrop.7areas.inhot_outrealm<-prepare_realm_input_7areas_plus_inhotoutrealm(table=world.table,name='neotrop',overlap=0.80)
#palearctic.7areas.inhot_outrealm<-prepare_realm_input_7areas_plus_inhotoutrealm(table=world.table,name='palearctic',overlap=0.80)
#
####have to move outputs from prepare_realm_simple_input_6states to correct folder (not in /hotspots_vertebrates/)

run_whole_realm_BioGeoBEARS_plusBSM<-function(name){
  setwd('~/Desktop/HOTSPOTS/repositories/hotspots_vertebrates/')
  source('./R/BioGeoBEARS_run_models.R')
  
  path<-paste('./output/mammals/new_tables_hotspots/whole_realms_6states/',name,'/',sep='')
  names.trees<-list.files(path,pattern='.tree')
  names.geography<-list.files(path,pattern='_geographyfile.txt')
  if(length(names.trees)!=length(names.geography)){
    cat('different lenghts, double check folder','\n')
  }
  for(i in 1:length(names.trees)){
      run_BioGeoBEARS_models(treefile = names.trees[i],geographyfile =names.geography[i],path = path,name=sub(names.trees[i],pattern=paste('_',name,'.tree',sep=''),replacement=''))
      setwd('~/Desktop/HOTSPOTS/repositories/hotspots_vertebrates/')
  }
  list.aicc.files<-list.files(path,pattern='_AICc_rellike_formatted.txt')
  #read AICC tables, select the best fitting model (sort by AICc weight), store the best fitting models in a list
  model.list<-list(0)
  for (i in 1:length(list.aicc.files)){
    table<-read.table(file=paste(path,list.aicc.files[i],sep=''),header=T,sep='\t')
    table<-table[order(-table$AICc_wt),]
    cat(i,'   ',table[1,'AICc_wt'],'    ',row.names(table)[1],'\n')
    model.list[[i]]<-c(row.names(table)[1],table[1,'AICc_wt'])
  }
  results.BSM<-list(0)
  for (i in 1:length(model.list)){
    results.BSM[[i]]<-run_BioGeoBEARS_selectedmodel_BSM(treefile=names.trees[i],geographyfile=names.geography[i],path=path,name=sub(names.trees[i],pattern='_neotropical.tree',replacement=''),model_name=model.list[[i]][1]) 
    setwd('/Users/javier/Desktop/HOTSPOTS/repositories/hotspots_vertebrates')
  }
  setwd('/Users/javier/Desktop/HOTSPOTS/repositories/hotspots_vertebrates')
  assign(paste(name,'.whole_realms_BSM',sep=''),results.BSM)
  saveRDS(get(paste(name,'.whole_realms_BSM',sep='')),file=paste(name,'.whole_realms_BSM.RDS',sep=''))
}

run_whole_realm_BioGeoBEARS_plusBSM_7areas<-function(name){
  setwd('~/Desktop/HOTSPOTS/repositories/hotspots_vertebrates/')
  source('./R/BioGeoBEARS_run_models.R')
  
  path<-paste('./output/mammals/new_tables_hotspots/whole_realms_7areas/',name,'/',sep='')
  names.trees<-list.files(path,pattern='.tree')
  names.geography<-list.files(path,pattern='_geographyfile.txt')
  if(length(names.trees)!=length(names.geography)){
    cat('different lenghts, double check folder','\n')
  }
  for(i in 1:length(names.trees)){
    run_BioGeoBEARS_models(treefile = names.trees[i],geographyfile =names.geography[i],path = path,name=sub(names.trees[i],pattern=paste('_',name,'.tree',sep=''),replacement=''))
    setwd('~/Desktop/HOTSPOTS/repositories/hotspots_vertebrates/')
  }
  list.aicc.files<-list.files(path,pattern='_AICc_rellike_formatted.txt')
  #read AICC tables, select the best fitting model (sort by AICc weight), store the best fitting models in a list
  model.list<-list(0)
  for (i in 1:length(list.aicc.files)){
    table<-read.table(file=paste(path,list.aicc.files[i],sep=''),header=T,sep='\t')
    table<-table[order(-table$AICc_wt),]
    cat(i,'   ',table[1,'AICc_wt'],'    ',row.names(table)[1],'\n')
    model.list[[i]]<-c(row.names(table)[1],table[1,'AICc_wt'])
  }
  results.BSM<-list(0)
  for (i in 1:length(model.list)){
    results.BSM[[i]]<-run_BioGeoBEARS_selectedmodel_BSM(treefile=names.trees[i],geographyfile=names.geography[i],path=path,name=sub(names.trees[i],pattern='_neotropical.tree',replacement=''),model_name=model.list[[i]][1]) 
    setwd('/Users/javier/Desktop/HOTSPOTS/repositories/hotspots_vertebrates')
  }
  setwd('/Users/javier/Desktop/HOTSPOTS/repositories/hotspots_vertebrates')
  assign(paste(name,'.whole_realms_7areas_BSM',sep=''),results.BSM)
  saveRDS(get(paste(name,'.whole_realms_7areas_BSM',sep='')),file=paste(name,'.whole_realms_7areas_BSM.RDS',sep=''))
}

run_whole_realm_BioGeoBEARS_plusBSM_7areas_inhotoutrealm<-function(name,path){
  setwd('/Users/javier/Documents/Work/HOTSPOTS/repositories/hotspots_vertebrates/')
  source('./R/BioGeoBEARS_run_models.R')
  
  path<-path
  names.trees<-list.files(path,pattern='.tree$')
  names.geography<-list.files(path,pattern='_geographyfile.txt')
  if(length(names.trees)!=length(names.geography)){
    cat('different lenghts, double check folder','\n')
  }
  for(i in 1:length(names.trees)){
    run_BioGeoBEARS_models(treefile = names.trees[i],geographyfile =names.geography[i],path = path,name=sub(names.trees[i],pattern=paste('_',name,'.tree',sep=''),replacement=''))
    setwd('/Users/javier/Documents/Work/HOTSPOTS/repositories/hotspots_vertebrates/')
  }
  list.aicc.files<-list.files(path,pattern='_AICc_rellike_formatted.txt')
  #read AICC tables, select the best fitting model (sort by AICc weight), store the best fitting models in a list
  model.list<-list(0)
  for (i in 1:length(list.aicc.files)){
    table<-read.table(file=paste(path,list.aicc.files[i],sep=''),header=T,sep='\t')
    table<-table[order(-table$AICc_wt),]
    cat(i,'   ',table[1,'AICc_wt'],'    ',row.names(table)[1],'\n')
    model.list[[i]]<-c(row.names(table)[1],table[1,'AICc_wt'])
  }
  results.BSM<-list(0)
  results.BSM.summary<-list(0)
  for (i in 1:length(model.list)){
    results.BSM[[i]]<-run_BioGeoBEARS_selectedmodel_BSM_object(treefile=names.trees[i],geographyfile=names.geography[i],path=path,name=sub(names.trees[i],pattern=paste('_',name,'.tree',sep=''),replacement=''),model_name=model.list[[i]][1],nreplicates=50) 
    setwd('/Users/javier/Documents/Work/HOTSPOTS/repositories/hotspots_vertebrates/')
    assign(paste(name,'.whole_realms_7areas_inhotoutrealm_BSM',sep=''),results.BSM[[i]])
    saveRDS(get(paste(name,'.whole_realms_7areas_inhotoutrealm_BSM',sep='')),file=paste(path,name,'.whole_realms_7areas_inhotoutrealm_BSM.RDS',sep=''))
    setwd('/Users/javier/Documents/Work/HOTSPOTS/repositories/hotspots_vertebrates/')
    results.BSM.summary[[i]]<-summarise_BSMobject(geographyfile=names.geography[i],path=path,name=sub(names.trees[i],pattern=paste('_',name,'.tree',sep=''),replacement=''),BSM.objectfile=paste(path,name,'.whole_realms_7areas_inhotoutrealm_BSM.RDS',sep=''))
    setwd('/Users/javier/Documents/Work/HOTSPOTS/repositories/hotspots_vertebrates/')
    assign(paste(name,'.whole_realms_7areas_inhotoutrealm_BSM_summary',sep=''),results.BSM.summary[[i]])
    saveRDS(get(paste(name,'.whole_realms_7areas_inhotoutrealm_BSM_summary',sep='')),file=paste(path,name,'.whole_realms_7areas_inhotoutrealm_BSM_summary.RDS',sep=''))
    
  }  
  setwd('/Users/javier/Documents/Work/HOTSPOTS/repositories/hotspots_vertebrates/')
  
}


run_whole_realm_BioGeoBEARS_plusBSM_7areas_inhotoutrealm_cluster<-function(name,path){
  setwd('/home/ji247/hotspots_vertebrates/')
  source('./R/BioGeoBEARS_run_models.R')
  
  path<-path
  names.trees<-list.files(path,pattern='.tree$')
  names.geography<-list.files(path,pattern='_geographyfile.txt')
  if(length(names.trees)!=length(names.geography)){
    cat('different lenghts, double check folder','\n')
  }
  for(i in 1:length(names.trees)){
    run_BioGeoBEARS_models(treefile = names.trees[i],geographyfile =names.geography[i],path = path,name=sub(names.trees[i],pattern=paste('_',name,'.tree',sep=''),replacement=''))
    setwd('/home/ji247/hotspots_vertebrates/')
  }
  list.aicc.files<-list.files(path,pattern='_AICc_rellike_formatted.txt')
  #read AICC tables, select the best fitting model (sort by AICc weight), store the best fitting models in a list
  model.list<-list(0)
  for (i in 1:length(list.aicc.files)){
    table<-read.table(file=paste(path,list.aicc.files[i],sep=''),header=T,sep='\t')
    table<-table[order(-table$AICc_wt),]
    cat(i,'   ',table[1,'AICc_wt'],'    ',row.names(table)[1],'\n')
    model.list[[i]]<-c(row.names(table)[1],table[1,'AICc_wt'])
  }
  results.BSM<-list(0)
  results.BSM.summary<-list(0)
  for (i in 1:length(model.list)){
    results.BSM[[i]]<-run_BioGeoBEARS_selectedmodel_BSM_object(treefile=names.trees[i],geographyfile=names.geography[i],path=path,name=sub(names.trees[i],pattern=paste('_',name,'.tree',sep=''),replacement=''),model_name=model.list[[i]][1],nreplicates=50) 
    setwd('/home/ji247/hotspots_vertebrates/')
    assign(paste(name,'.whole_realms_7areas_inhotoutrealm_BSM',sep=''),results.BSM[[i]])
    saveRDS(get(paste(name,'.whole_realms_7areas_inhotoutrealm_BSM',sep='')),file=paste(path,name,'.whole_realms_7areas_inhotoutrealm_BSM.RDS',sep=''))
    setwd('/home/ji247/hotspots_vertebrates/')
    results.BSM.summary[[i]]<-summarise_BSMobject(geographyfile=names.geography[i],path=path,name=sub(names.trees[i],pattern=paste('_',name,'.tree',sep=''),replacement=''),BSM.objectfile=paste(path,name,'.whole_realms_7areas_inhotoutrealm_BSM.RDS',sep=''))
    setwd('/home/ji247/hotspots_vertebrates/')
    assign(paste(name,'.whole_realms_7areas_inhotoutrealm_BSM_summary',sep=''),results.BSM.summary[[i]])
    saveRDS(get(paste(name,'.whole_realms_7areas_inhotoutrealm_BSM_summary',sep='')),file=paste(path,name,'.whole_realms_7areas_inhotoutrealm_BSM_summary.RDS',sep=''))
    
  }  
  setwd('/home/ji247/hotspots_vertebrates/')
  
}

#####for 6 states
####source('./BioGeoBEARS_run_models.R')
####setwd('~/Desktop/HOTSPOTS/repositories/hotspots_vertebrates/')
####run_whole_realm_BioGeoBEARS_plusBSM('afrotrop')
####setwd('~/Desktop/HOTSPOTS/repositories/hotspots_vertebrates/')
####run_whole_realm_BioGeoBEARS_plusBSM('austral')
####setwd('~/Desktop/HOTSPOTS/repositories/hotspots_vertebrates/')
####run_whole_realm_BioGeoBEARS_plusBSM('indo')
####setwd('~/Desktop/HOTSPOTS/repositories/hotspots_vertebrates/')
####run_whole_realm_BioGeoBEARS_plusBSM('nearctic')
####setwd('~/Desktop/HOTSPOTS/repositories/hotspots_vertebrates/')
####run_whole_realm_BioGeoBEARS_plusBSM('neotrop')
####setwd('~/Desktop/HOTSPOTS/repositories/hotspots_vertebrates/')
####run_whole_realm_BioGeoBEARS_plusBSM('palearctic')
####setwd('~/Desktop/HOTSPOTS/repositories/hotspots_vertebrates/')
####
#####for 7 areas
####source('./BioGeoBEARS_run_models.R')
####setwd('~/Desktop/HOTSPOTS/repositories/hotspots_vertebrates/')
####run_whole_realm_BioGeoBEARS_plusBSM_7areas('afrotrop')
####setwd('~/Desktop/HOTSPOTS/repositories/hotspots_vertebrates/')
####run_whole_realm_BioGeoBEARS_plusBSM_7areas('austral')
####setwd('~/Desktop/HOTSPOTS/repositories/hotspots_vertebrates/')
####run_whole_realm_BioGeoBEARS_plusBSM_7areas('indo')
####setwd('~/Desktop/HOTSPOTS/repositories/hotspots_vertebrates/')
####run_whole_realm_BioGeoBEARS_plusBSM_7areas('nearctic')
####setwd('~/Desktop/HOTSPOTS/repositories/hotspots_vertebrates/')
####run_whole_realm_BioGeoBEARS_plusBSM_7areas('neotrop')
####setwd('~/Desktop/HOTSPOTS/repositories/hotspots_vertebrates/')
####run_whole_realm_BioGeoBEARS_plusBSM_7areas('palearctic')
####setwd('~/Desktop/HOTSPOTS/repositories/hotspots_vertebrates/')
####
#####for 7 areas inhotoutrealm
####source('./BioGeoBEARS_run_models.R')
####setwd('~/Desktop/HOTSPOTS/repositories/hotspots_vertebrates/')
####run_whole_realm_BioGeoBEARS_plusBSM_7areas_inhotoutrealm('afrotrop')
####setwd('~/Desktop/HOTSPOTS/repositories/hotspots_vertebrates/')
####run_whole_realm_BioGeoBEARS_plusBSM_7areas_inhotoutrealm('austral')
####setwd('~/Desktop/HOTSPOTS/repositories/hotspots_vertebrates/')
####run_whole_realm_BioGeoBEARS_plusBSM_7areas_inhotoutrealm('indo')
####setwd('~/Desktop/HOTSPOTS/repositories/hotspots_vertebrates/')
####run_whole_realm_BioGeoBEARS_plusBSM_7areas_inhotoutrealm('nearctic')
####setwd('~/Desktop/HOTSPOTS/repositories/hotspots_vertebrates/')
####run_whole_realm_BioGeoBEARS_plusBSM_7areas_inhotoutrealm('neotrop')
####setwd('~/Desktop/HOTSPOTS/repositories/hotspots_vertebrates/')
####run_whole_realm_BioGeoBEARS_plusBSM_7areas_inhotoutrealm('palearctic')
####setwd('~/Desktop/HOTSPOTS/repositories/hotspots_vertebrates/')
####
####
