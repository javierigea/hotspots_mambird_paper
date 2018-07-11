library(ape)
library(picante)
library(geiger)
######################################################################################
#this generates a weighted endemism measure in a grid file
#species.grid.table: _realms_species_gridoccurrence_table.txt
#cell.grid.table: _realms_richness_grid_table.txt
#path: to store the outputfile
#name:to append to outputfile name
generate_grid_weightedendemism<-function(speciesgrid.tablefile,cellgrid.tablefile,path,name){
  species.grid.table<-read.table(speciesgrid.tablefile,header=T,sep='\t')
  cell.grid.table<-read.table(cellgrid.tablefile,header=T,sep='\t')
  #add a space at the beginning and end of the cell string for grep to work
  species.grid.table$cells<-paste(' ',species.grid.table$cells,' ')
  WE<-lapply(c(1:nrow(cell.grid.table)),function(x) sum(1/as.numeric(species.grid.table[grep(paste(' ',x,' ',sep=''),species.grid.table$cells),'range.cells']),na.rm=TRUE))
  cell.species.df<-cell.grid.table
  cell.species.df$number.of.species.wend<-unlist(WE)
  cell.species.df[cell.species.df$number.of.species.wend==0,'number.of.species.wend']<-NA
  write.table(cell.species.df,file=paste(path,'/100_',name,'_realms_richness_wend_grid_table.txt',sep=''),sep='\t',row.names=F,quote=F)
}
##########################################################################################

######################################################################################
#this generates a list with all cells in the grid, the species richness and names of the species in each cell
#species.grid.table: _realms_species_gridoccurrence_table.txt
#cell.grid.table: _realms_richness_grid_table.txt
#path: to store the outputfile
#name:to append to outputfile name
get_speciesnames_grid<-function(species.grid.table,cell.grid.table,path,name){
  species.grid.table<-read.table(species.grid.table,header=T,sep='\t')
  cell.grid.table<-read.table(cell.grid.table,header=T,sep='\t')
  species.grid.table$cells<-paste(' ',species.grid.table$cells,' ')
  cell.species.matrix<-matrix(NA,ncol=2,nrow=nrow(cell.grid.table))
  for (i in 1:nrow(cell.grid.table)){
    cat(i,'\n')
    species.in.cell<-as.character(species.grid.table[grep(paste(' ',i,' ',sep=''),species.grid.table$cells),]$spp)
    species.in.cell<-paste(species.in.cell, collapse = " ")
    cell.species.matrix[i,]<-cbind(i,species.in.cell)
  }
  cell.species.df<-as.data.frame(cell.species.matrix,stringsAsFactors = F)
  cell.species.df$number.of.species<-cell.grid.table[,2]
  colnames(cell.species.df)<-c('cell','species','number.of.species')
  write.table(cell.species.df,file=paste(path,name,'_richness_grid_speciesnames_table.txt',sep=''),sep='\t',row.names=F,quote=F)
}
################################################################
###this function defines hotspot cells in a grid table
###hotspots are defined using a quantile value (q80 = hotspots are the top 20%) + a variable.name
#grid.tablefile is the tablefile to be used, e.g. _realms_richness_wend_grid_table.txt
#quantile defines the top % that will be hotspots
#variable.name is the column name
#path to store the outputfile
#name to append to outputfile
define_hotspots<-function(grid.tablefile,quantile,variable.name,path,name){
  table<-read.table(grid.tablefile,header=T,sep='\t')
  colnames(table)[1]<-'cells'
  #check presence of column
  if(!variable.name%in%colnames(table)){
    return('cannot find variable.name in tablefile')
  }
  table$hotspot<-0
  threshold<-quantile(table[,variable.name],quantile,na.rm=TRUE)
  table[which(table[,variable.name]>threshold),'hotspot']<-1
  table$hotspot<-as.factor(table$hotspot)
  write.table(table,file=paste(path,'/100_',name,'_',variable.name,'_',quantile,'.txt',sep=''),sep='\t',row.names=F,quote=F)
}

#this outputs a hotspot file using ranges, it selects the narrow ranged species (range <= threshold) and hotspots are cells where these occur
define_narrowrangedspecies_hotspots<-function(species.ranges.tablefile,threshold,path){
  #first get narrow ranged species (i.e. species present in 10 or less cells)
  species.ranges<-read.table(species.ranges.tablefile,header=T,sep='\t',stringsAsFactors = F)
  cat('number of narrow ranged species','\n')
  cat(nrow(species.ranges[species.ranges$range.cells<=threshold,]),'\n')
  #there are 1438 species with narrow ranges.
  cells.narrow.ranged.species<-species.ranges[species.ranges$range.cells<=threshold,'cells']
  cells.narrow.ranged.species<-strsplit(cells.narrow.ranged.species,' ')
  cells.narrow.ranged.species<-lapply(cells.narrow.ranged.species,function(x)x[x!=''])
  cells.narrow.ranged.species<-unique(unlist(cells.narrow.ranged.species))
  #create table
  table.narrow.ranged.hotspots<-as.data.frame(cbind(c(1:16517)),stringsAsFactors = F)
  colnames(table.narrow.ranged.hotspots)<-'cells'
  table.narrow.ranged.hotspots$hotspot<-as.numeric(table.narrow.ranged.hotspots$cells%in%cells.narrow.ranged.species)
  write.table(table.narrow.ranged.hotspots,file=paste(path,'/100_all_realms_narrow.ranged.species_hotspots.txt',sep=''),quote=F,sep='\t',row.names=F)
}
######################################################################################
#this first divides species in a DR_table into 4 quartiles
#and then counts the species richness in each quartile in the cells of a grid
#dr.table: DR_*.txt
#species.grid.table: _realms_species_gridoccurrence_table.txt
#cell.grid.table: _realms_richness_grid_table.txt
#path: to store the outputfile
#name:to append to outputfile name
generate_grid_richness_quartilesDR<-function(dr.table,species.grid.tablefile,cellgrid.tablefile,path,name){
  species.grid.table<-read.table(species.grid.tablefile,header=T,sep='\t')
  cell.grid.table<-read.table(cellgrid.tablefile,header=T,sep='\t')
  #add a space at the beginning and end of the cell string for grep to work
  species.grid.table$cells<-paste(' ',species.grid.table$cells,' ')
  dr.table<-read.table(dr.table,header=T,sep='\t')
  colnames(dr.table)[2]<-'DR'
  #getting the quartiles of DR into the DR table
  dr.table[dr.table$DR<=quantile(dr.table$DR,0.25),'quartile.DR']<-'Q1'
  dr.table[dr.table$DR>quantile(dr.table$DR,0.25)&dr.table$DR<=quantile(dr.table$DR,0.5),'quartile.DR']<-'Q2'
  dr.table[dr.table$DR>quantile(dr.table$DR,0.5)&dr.table$DR<=quantile(dr.table$DR,0.75),'quartile.DR']<-'Q3'
  dr.table[dr.table$DR>quantile(dr.table$DR,0.75),'quartile.DR']<-'Q4'
  #add a space at the beginning and end of the cell string for grep to work
  species.grid.table$cells<-paste(' ',species.grid.table$cells,' ')
  species.grid.table<-merge(species.grid.table,dr.table,by.x='spp',by.y='Species',all.x=TRUE)
  cell.species.matrix<-matrix(NA,ncol=6,nrow=nrow(cell.grid.table))
  for (i in 1:nrow(cell.grid.table)){
    cat(i,'\n')
    species.in.cell<-as.character(species.grid.table[grep(paste(' ',i,' ',sep=''),species.grid.table$cells),]$spp)
    quartiles.in.cell.vector<-as.character(species.grid.table[grep(paste(' ',i,' ',sep=''),species.grid.table$cells),]$quartile.DR)
    quartiles.in.cell<-unname(sapply(c('Q1','Q2','Q3','Q4'),function(x) length(which(quartiles.in.cell.vector==x))))
    species.in.cell<-paste(species.in.cell, collapse = " ")
    cell.species.matrix[i,]<-cbind(i,species.in.cell,rbind(quartiles.in.cell))
  }
  cell.species.df<-as.data.frame(cell.species.matrix,stringsAsFactors = F)
  cell.species.df$number.of.species<-cell.grid.table[,2]
  colnames(cell.species.df)<-c('cell','species','Q1','Q2','Q3','Q4','number.of.species')
  cell.species.df[3:7] <- lapply(cell.species.df[3:7], as.numeric)
  cell.species.df$proportion.Q1<-cell.species.df$Q1/cell.species.df$number.of.species
  cell.species.df$proportion.Q2<-cell.species.df$Q2/cell.species.df$number.of.species
  cell.species.df$proportion.Q3<-cell.species.df$Q3/cell.species.df$number.of.species
  cell.species.df$proportion.Q4<-cell.species.df$Q4/cell.species.df$number.of.species
  write.table(cell.species.df,file=paste(path,'/',name,'_quartilesDR_grid_table.txt',sep=''),sep='\t',row.names=F,quote=F)
}

#this first divides species in a DR_table into 4 quartiles
#and then counts the species richness in each quartile in the cells of a grid
#dr.table: DR_*.txt
#species.grid.table: _realms_species_gridoccurrence_table.txt
#cell.grid.table: _realms_richness_grid_table.txt
#path: to store the outputfile
#name:to append to outputfile name
generate_grid_richness_quartilesage<-function(age.table,species.grid.tablefile,cellgrid.tablefile,path,name){
  species.grid.table<-read.table(species.grid.tablefile,header=T,sep='\t')
  cell.grid.table<-read.table(cellgrid.tablefile,header=T,sep='\t')
  #add a space at the beginning and end of the cell string for grep to work
  species.grid.table$cells<-paste(' ',species.grid.table$cells,' ')
  age.table<-read.table(age.table,header=T,sep='\t')
  colnames(age.table)[2]<-'species.age'
  #getting the quartiles of species.age into the species.age table
  age.table[age.table$species.age<=quantile(age.table$species.age,0.25),'quartile.species.age']<-'Q1'
  age.table[age.table$species.age>quantile(age.table$species.age,0.25)&age.table$species.age<=quantile(age.table$species.age,0.5),'quartile.species.age']<-'Q2'
  age.table[age.table$species.age>quantile(age.table$species.age,0.5)&age.table$species.age<=quantile(age.table$species.age,0.75),'quartile.species.age']<-'Q3'
  age.table[age.table$species.age>quantile(age.table$species.age,0.75),'quartile.species.age']<-'Q4'
  #add a space at the beginning and end of the cell string for grep to work
  species.grid.table$cells<-paste(' ',species.grid.table$cells,' ')
  species.grid.table<-merge(species.grid.table,age.table,by.x='spp',by.y='species',all.x=TRUE)
  cell.species.matrix<-matrix(NA,ncol=6,nrow=nrow(cell.grid.table))
  for (i in 1:nrow(cell.grid.table)){
    cat(i,'\n')
    species.in.cell<-as.character(species.grid.table[grep(paste(' ',i,' ',sep=''),species.grid.table$cells),]$spp)
    quartiles.in.cell.vector<-as.character(species.grid.table[grep(paste(' ',i,' ',sep=''),species.grid.table$cells),]$quartile.species.age)
    quartiles.in.cell<-unname(sapply(c('Q1','Q2','Q3','Q4'),function(x) length(which(quartiles.in.cell.vector==x))))
    species.in.cell<-paste(species.in.cell, collapse = " ")
    cell.species.matrix[i,]<-cbind(i,species.in.cell,rbind(quartiles.in.cell))
  }
  cell.species.df<-as.data.frame(cell.species.matrix,stringsAsFactors = F)
  cell.species.df$number.of.species<-cell.grid.table[,2]
  colnames(cell.species.df)<-c('cell','species','Q1','Q2','Q3','Q4','number.of.species')
  cell.species.df[3:7] <- lapply(cell.species.df[3:7], as.numeric)
  cell.species.df$proportion.Q1<-cell.species.df$Q1/cell.species.df$number.of.species
  cell.species.df$proportion.Q2<-cell.species.df$Q2/cell.species.df$number.of.species
  cell.species.df$proportion.Q3<-cell.species.df$Q3/cell.species.df$number.of.species
  cell.species.df$proportion.Q4<-cell.species.df$Q4/cell.species.df$number.of.species
  write.table(cell.species.df,file=paste(path,'/',name,'_quartilesspecies.age_grid_table.txt',sep=''),sep='\t',row.names=F,quote=F)
}

######################################################################################
#this first divides species in a DR_table into 4 quartiles
#and then counts the species richness in each quartile in the cells of a grid
#BAMM.table: *_BAMMTipRates.txt
#species.grid.table: _realms_species_gridoccurrence_table.txt
#cell.grid.table: _realms_richness_grid_table.txt
#path: to store the outputfile
#name:to append to outputfile name
generate_grid_richness_quartilesBAMM<-function(BAMM.table,species.grid.tablefile,cellgrid.tablefile,path,name){
  species.grid.table<-read.table(species.grid.tablefile,header=T,sep='\t')
  cell.grid.table<-read.table(cellgrid.tablefile,header=T,sep='\t')
  #add a space at the beginning and end of the cell string for grep to work
  species.grid.table$cells<-paste(' ',species.grid.table$cells,' ')
  BAMM.table<-read.table(BAMM.table,header=T,sep='\t')
  #getting the quartiles of lambda.avg into the BAMM table
  BAMM.table[BAMM.table$lambda.avg<=quantile(BAMM.table$lambda.avg,0.25),'quartile.lambda.avg']<-'Q1.lambda.avg'
  BAMM.table[BAMM.table$lambda.avg>quantile(BAMM.table$lambda.avg,0.25)&BAMM.table$lambda.avg<=quantile(BAMM.table$lambda.avg,0.5),'quartile.lambda.avg']<-'Q2.lambda.avg'
  BAMM.table[BAMM.table$lambda.avg>quantile(BAMM.table$lambda.avg,0.5)&BAMM.table$lambda.avg<=quantile(BAMM.table$lambda.avg,0.75),'quartile.lambda.avg']<-'Q3.lambda.avg'
  BAMM.table[BAMM.table$lambda.avg>quantile(BAMM.table$lambda.avg,0.75),'quartile.lambda.avg']<-'Q4.lambda.avg'
  BAMM.table[BAMM.table$mu.avg<=quantile(BAMM.table$mu.avg,0.25),'quartile.mu.avg']<-'Q1.mu.avg'
  BAMM.table[BAMM.table$mu.avg>quantile(BAMM.table$mu.avg,0.25)&BAMM.table$mu.avg<=quantile(BAMM.table$mu.avg,0.5),'quartile.mu.avg']<-'Q2.mu.avg'
  BAMM.table[BAMM.table$mu.avg>quantile(BAMM.table$mu.avg,0.5)&BAMM.table$mu.avg<=quantile(BAMM.table$mu.avg,0.75),'quartile.mu.avg']<-'Q3.mu.avg'
  BAMM.table[BAMM.table$mu.avg>quantile(BAMM.table$mu.avg,0.75),'quartile.mu.avg']<-'Q4.mu.avg'
  BAMM.table[BAMM.table$netdiv.avg<=quantile(BAMM.table$netdiv.avg,0.25),'quartile.netdiv.avg']<-'Q1.netdiv.avg'
  BAMM.table[BAMM.table$netdiv.avg>quantile(BAMM.table$netdiv.avg,0.25)&BAMM.table$netdiv.avg<=quantile(BAMM.table$netdiv.avg,0.5),'quartile.netdiv.avg']<-'Q2.netdiv.avg'
  BAMM.table[BAMM.table$netdiv.avg>quantile(BAMM.table$netdiv.avg,0.5)&BAMM.table$netdiv.avg<=quantile(BAMM.table$netdiv.avg,0.75),'quartile.netdiv.avg']<-'Q3.netdiv.avg'
  BAMM.table[BAMM.table$netdiv.avg>quantile(BAMM.table$netdiv.avg,0.75),'quartile.netdiv.avg']<-'Q4.netdiv.avg'
  #add a space at the beginning and end of the cell string for grep to work
  species.grid.table$cells<-paste(' ',species.grid.table$cells,' ')
  species.grid.table<-merge(species.grid.table,BAMM.table,by.x='spp',by.y='spp',all.x=TRUE)
  cell.species.matrix<-matrix(NA,ncol=14,nrow=nrow(cell.grid.table))
  for (i in 1:nrow(cell.grid.table)){
    cat(i,'\n')
    species.in.cell<-as.character(species.grid.table[grep(paste(' ',i,' ',sep=''),species.grid.table$cells),]$spp)
    quartiles.in.cell.vector<-as.character(species.grid.table[grep(paste(' ',i,' ',sep=''),species.grid.table$cells),c('quartile.lambda.avg','quartile.mu.avg','quartile.netdiv.avg')])
    quartiles.in.cell.vector<-as.character(c(species.grid.table[grep(paste(' ',i,' ',sep=''),species.grid.table$cells),c('quartile.lambda.avg')],species.grid.table[grep(paste(' ',i,' ',sep=''),species.grid.table$cells),c('quartile.mu.avg')],species.grid.table[grep(paste(' ',i,' ',sep=''),species.grid.table$cells),c('quartile.netdiv.avg')]))
    quartiles.in.cell<-unname(sapply(c('Q1.lambda.avg','Q2.lambda.avg','Q3.lambda.avg','Q4.lambda.avg','Q1.mu.avg','Q2.mu.avg','Q3.mu.avg','Q4.mu.avg','Q1.netdiv.avg','Q2.netdiv.avg','Q3.netdiv.avg','Q4.netdiv.avg'),function(x) length(which(quartiles.in.cell.vector==x))))
    species.in.cell<-paste(species.in.cell, collapse = " ")
    cell.species.matrix[i,]<-cbind(i,species.in.cell,rbind(quartiles.in.cell))
    
  }
  cell.species.df<-as.data.frame(cell.species.matrix,stringsAsFactors = F)
  cell.species.df$number.of.species<-cell.grid.table[,2]
  colnames(cell.species.df)<-c('cell','species','Q1.lambda.avg','Q2.lambda.avg','Q3.lambda.avg','Q4.lambda.avg','Q1.mu.avg','Q2.mu.avg','Q3.mu.avg','Q4.mu.avg','Q1.netdiv.avg','Q2.netdiv.avg','Q3.netdiv.avg','Q4.netdiv.avg','number.of.species')
  cell.species.df[3:15] <- lapply(cell.species.df[3:15], as.numeric)
  cell.species.df$proportion.Q1.lambda.avg<-cell.species.df$Q1.lambda.avg/cell.species.df$number.of.species
  cell.species.df$proportion.Q2.lambda.avg<-cell.species.df$Q2.lambda.avg/cell.species.df$number.of.species
  cell.species.df$proportion.Q3.lambda.avg<-cell.species.df$Q3.lambda.avg/cell.species.df$number.of.species
  cell.species.df$proportion.Q4.lambda.avg<-cell.species.df$Q4.lambda.avg/cell.species.df$number.of.species
  write.table(cell.species.df,file=paste(path,'/',name,'_quartilesBAMMTipRates_grid_table.txt',sep=''),sep='\t',row.names=F,quote=F)
}

######################################################################################
#this first divides species in a DR_table into 4 quartiles
#and then counts the species richness in each quartile in the cells of a grid
#dr.table: DR_*.txt
#species.grid.table: _realms_species_gridoccurrence_table.txt
#cell.grid.table: _realms_richness_grid_table.txt
#path: to store the outputfile
#name:to append to outputfile name
generate_grid_richness_quartilesDR_order<-function(dr.table,species.grid.tablefile,cellgrid.tablefile,path,name){
  species.grid.table<-read.table(species.grid.tablefile,header=T,sep='\t')
  cell.grid.table<-read.table(cellgrid.tablefile,header=T,sep='\t')
  #add a space at the beginning and end of the cell string for grep to work
  species.grid.table$cells<-paste(' ',species.grid.table$cells,' ')
  dr.table<-read.table(dr.table,header=T,sep='\t')
  #subset DR.table to species in order
  dr.table<-dr.table[dr.table$Species%in%species.grid.table$spp,]
  #getting the quartiles of DR into the DR table
  dr.table[dr.table$DR<=quantile(dr.table$DR,0.25),'quartile.DR']<-'Q1'
  dr.table[dr.table$DR>quantile(dr.table$DR,0.25)&dr.table$DR<=quantile(dr.table$DR,0.5),'quartile.DR']<-'Q2'
  dr.table[dr.table$DR>quantile(dr.table$DR,0.5)&dr.table$DR<=quantile(dr.table$DR,0.75),'quartile.DR']<-'Q3'
  dr.table[dr.table$DR>quantile(dr.table$DR,0.75),'quartile.DR']<-'Q4'
  #add a space at the beginning and end of the cell string for grep to work
  species.grid.table$cells<-paste(' ',species.grid.table$cells,' ')
  species.grid.table<-merge(species.grid.table,dr.table,by.x='spp',by.y='Species',all.x=TRUE)
  cell.species.matrix<-matrix(NA,ncol=6,nrow=nrow(cell.grid.table))
  for (i in 1:nrow(cell.grid.table)){
    cat(i,'\n')
    species.in.cell<-as.character(species.grid.table[grep(paste(' ',i,' ',sep=''),species.grid.table$cells),]$spp)
    quartiles.in.cell.vector<-as.character(species.grid.table[grep(paste(' ',i,' ',sep=''),species.grid.table$cells),]$quartile.DR)
    quartiles.in.cell<-unname(sapply(c('Q1','Q2','Q3','Q4'),function(x) length(which(quartiles.in.cell.vector==x))))
    species.in.cell<-paste(species.in.cell, collapse = " ")
    cell.species.matrix[i,]<-cbind(i,species.in.cell,rbind(quartiles.in.cell))
  }
  cell.species.df<-as.data.frame(cell.species.matrix,stringsAsFactors = F)
  cell.species.df$number.of.species<-cell.grid.table[,2]
  colnames(cell.species.df)<-c('cell','species','Q1','Q2','Q3','Q4','number.of.species')
  cell.species.df[3:7] <- lapply(cell.species.df[3:7], as.numeric)
  cell.species.df$proportion.Q1<-cell.species.df$Q1/cell.species.df$number.of.species
  cell.species.df$proportion.Q2<-cell.species.df$Q2/cell.species.df$number.of.species
  cell.species.df$proportion.Q3<-cell.species.df$Q3/cell.species.df$number.of.species
  cell.species.df$proportion.Q4<-cell.species.df$Q4/cell.species.df$number.of.species
  write.table(cell.species.df,file=paste(path,'/',name,'_quartilesDR_grid_table.txt',sep=''),sep='\t',row.names=F,quote=F)
}

######################################################################################
#runs a linear model QuartileRichness ~ total richness and returns residual of each model
####dr.quartiles.gridtable: _quartilesDR_grid_table.txt
build_lm_speciesagequartiles<-function(species.age.quartiles.table,path,name){
  species.age.quartiles.table<-read.table(species.age.quartiles.table,header=T,sep='\t')
  species.age.quartiles.table<-species.age.quartiles.table[,-2]
  Q1.lm<-lm(species.age.quartiles.table$Q1~species.age.quartiles.table$number.of.species)
  Q2.lm<-lm(species.age.quartiles.table$Q2~species.age.quartiles.table$number.of.species)
  Q3.lm<-lm(species.age.quartiles.table$Q3~species.age.quartiles.table$number.of.species)
  Q4.lm<-lm(species.age.quartiles.table$Q4~species.age.quartiles.table$number.of.species)
  residuals.lm<-as.data.frame(cbind(Q1.lm$residuals,Q2.lm$residuals,Q3.lm$residuals,Q4.lm$residuals))
  residuals.lm$cell<-row.names(residuals.lm)
  row.names(residuals.lm)<-NULL
  residuals.lm<-residuals.lm[,c(5,1,2,3,4)]
  colnames(residuals.lm)[c(2,3,4,5)]<-c('Q1.residuals','Q2.residuals','Q3.residuals','Q4.residuals')
  species.age.gridtable.cells<-as.data.frame(species.age.quartiles.table[,1])
  colnames(species.age.gridtable.cells)<-'cells'
  colnames(residuals.lm)[1]<-'cells'
  residuals.lm$cells<-as.numeric(residuals.lm$cells)
  residuals.lm<-merge(residuals.lm,species.age.gridtable.cells,all.y=TRUE)
  residuals.lm<-residuals.lm[match(seq(1:nrow(residuals.lm)),residuals.lm$cell),]
  write.table(residuals.lm,file=paste(path,name,'_quartilesspeciesage_grid_residualstable.txt',sep=''),sep='\t',row.names=F,quote=F)
}



#runs a linear model QuartileRichness ~ total richness and returns residual of each model
####BAMM.quartiles.gridtable: _quartilesDR_grid_table.txt
build_lm_DRquartiles<-function(dr.quartiles.gridtable,path,name){
  dr.quartiles.gridtable<-read.table(dr.quartiles.gridtable,header=T,sep='\t')
  dr.quartiles.gridtable<-dr.quartiles.gridtable[,-2]
  Q1.lm<-lm(dr.quartiles.gridtable$Q1~dr.quartiles.gridtable$number.of.species)
  Q2.lm<-lm(dr.quartiles.gridtable$Q2~dr.quartiles.gridtable$number.of.species)
  Q3.lm<-lm(dr.quartiles.gridtable$Q3~dr.quartiles.gridtable$number.of.species)
  Q4.lm<-lm(dr.quartiles.gridtable$Q4~dr.quartiles.gridtable$number.of.species)
  residuals.lm<-as.data.frame(cbind(Q1.lm$residuals,Q2.lm$residuals,Q3.lm$residuals,Q4.lm$residuals))
  residuals.lm$cell<-row.names(residuals.lm)
  row.names(residuals.lm)<-NULL
  residuals.lm<-residuals.lm[,c(5,1,2,3,4)]
  colnames(residuals.lm)[c(2,3,4,5)]<-c('Q1.residuals','Q2.residuals','Q3.residuals','Q4.residuals')
  dr.gridtable.cells<-as.data.frame(dr.quartiles.gridtable[,1])
  colnames(dr.gridtable.cells)<-'cells'
  colnames(residuals.lm)[1]<-'cells'
  residuals.lm$cells<-as.numeric(residuals.lm$cells)
  residuals.lm<-merge(residuals.lm,dr.gridtable.cells,all.y=TRUE)
  residuals.lm<-residuals.lm[match(seq(1:nrow(residuals.lm)),residuals.lm$cell),]
  write.table(residuals.lm,file=paste(path,name,'_quartilesDR_grid_residualstable.txt',sep=''),sep='\t',row.names=F,quote=F)
}
#runs a linear model QuartileRichness ~ total richness and returns residual of each model
##BAMM.quartiles.gridtable: _quartilesBAMMTipRates_grid_table.txt
##
build_lm_BAMMquartiles<-function(BAMM.quartiles.gridtable,variable.name,path,name){
  BAMM.quartiles.gridtable<-read.table(BAMM.quartiles.gridtable,header=T,sep='\t')
  BAMM.quartiles.gridtable<-BAMM.quartiles.gridtable[,-2]
  Q1.lm<-lm(BAMM.quartiles.gridtable[,paste('Q1.',variable.name,sep='')]~BAMM.quartiles.gridtable[,'number.of.species'])
  Q2.lm<-lm(BAMM.quartiles.gridtable[,paste('Q2.',variable.name,sep='')]~BAMM.quartiles.gridtable[,'number.of.species'])
  Q3.lm<-lm(BAMM.quartiles.gridtable[,paste('Q3.',variable.name,sep='')]~BAMM.quartiles.gridtable[,'number.of.species'])
  Q4.lm<-lm(BAMM.quartiles.gridtable[,paste('Q4.',variable.name,sep='')]~BAMM.quartiles.gridtable[,'number.of.species'])
  residuals.lm<-as.data.frame(cbind(Q1.lm$residuals,Q2.lm$residuals,Q3.lm$residuals,Q4.lm$residuals))
  residuals.lm$cell<-row.names(residuals.lm)
  row.names(residuals.lm)<-NULL
  residuals.lm<-residuals.lm[,c(5,1,2,3,4)]
  colnames(residuals.lm)[c(2,3,4,5)]<-c('Q1.residuals','Q2.residuals','Q3.residuals','Q4.residuals')
  dr.gridtable.cells<-as.data.frame(BAMM.quartiles.gridtable[,1])
  colnames(dr.gridtable.cells)<-'cells'
  colnames(residuals.lm)[1]<-'cells'
  residuals.lm$cells<-as.numeric(residuals.lm$cells)
  residuals.lm<-merge(residuals.lm,dr.gridtable.cells,all.y=TRUE)
  residuals.lm<-residuals.lm[match(seq(1:nrow(residuals.lm)),residuals.lm$cell),]
  write.table(residuals.lm,file=paste(path,variable.name,'_',name,'_quartilesBAMMTipRates_grid_residualstable.txt',sep=''),sep='\t',row.names=F,quote=F)
}

########################################################################################
#this subsets a total species richness grid to a particular Order
#species.grid.table: _realms_species_gridoccurrence_table.txt
#dictionary file is the path to iucnmammaltaxonomy.csv (downloaded from IUCN)
#order name is the order to select
#path is to output the file
#name to append to output file
subset_speciesgridoccurrence_to_Order<-function(species.grid.tablefile,dictionaryfile,order.name,path,name){
  species.grid.table<-read.table(species.grid.tablefile,header=T,sep='\t',stringsAsFactors = F)
  spp<-species.grid.table$spp
  table.synonyms<-read.csv(dictionaryfile,header=T)
  table.synonyms<-data.frame(table.synonyms$Genus,table.synonyms$Species,paste(table.synonyms$Genus,table.synonyms$Species,sep=' '),table.synonyms$Order,stringsAsFactors = FALSE)
  colnames(table.synonyms)<-c('Genus','Species','Binomial','Order')
  table.synonyms$Binomial<-gsub(table.synonyms$Binomial,pattern=' ',replacement='_')
  #select order
  table.synonyms<-table.synonyms[table.synonyms$Order==toupper(order.name),]
  species.grid.table<-species.grid.table[(species.grid.table$spp%in%table.synonyms$Binomial),]
  #write table
  write.table(species.grid.table,file=paste(path,name,'_',order.name,'_gridoccurrence_table.txt',sep=''),sep='\t',quote=F,row.names=F)
  #calculate cell species richness
  cells<-species.grid.table$cells
  cells<-sapply(cells,function(x)as.numeric(unlist(strsplit(x,split=' '))))
  cells<-sapply(cells,function(x)x[!is.na(x)])
  cells<-unlist(cells)
  cell.counts<-as.data.frame(table(cells),stringsAsFactors = F)
  colnames(cell.counts)<-c('cells','number.of.species')
  cell.counts$cells<-as.numeric(cell.counts$cells)
  ncells<-as.data.frame(c(1:16517),stringsAsFactors=F)
  colnames(ncells)<-'cells'
  species.cell.table<-merge(ncells,cell.counts,all.x=TRUE,by='cells')
  write.table(species.cell.table,file=paste(path,name,'_',order.name,'_richness_grid_table.txt',sep=''),sep='\t',quote=F,row.names=F)
  
}

calculate_species_ages<-function(treefile,path,name){
  tree<-read.tree(treefile)
  #calculating species ages
  cat('calculating species ages','\n')
  node.age(tree)->phy.age
  cbind(phy.age$edge,phy.age$age, tree$edge.length)->BL.position
  max(phy.age$age)-BL.position[,3]->dist.tip
  cbind(BL.position,dist.tip)->BL.positions
  BL.positions[,5]+BL.positions[,4]->ages
  cbind(BL.positions,ages)->BL.positions
  as.data.frame(BL.positions)->node.ages
  names(node.ages)<-c("parental.node","daughter.node","dist.root","BL","dist.tip","mrca.age")
  ## node.ages is a data frame listing as variables the identity of parental and
  #daughter nodes, the distance from the root and from the present of each node,
  #the branch lenght and the age of the most recent common ancestor
  node.ages[node.ages[,2]<length(tree$tip)+1,]->species.ages
  row.names(species.ages)<-tree$tip
  ## species ages is node.ages data frame reduced to the tips (species)
  species.ages<-species.ages[order(row.names(species.ages)),]
  output.table<-as.data.frame(cbind(row.names(species.ages),species.ages$mrca.age))
  colnames(output.table)<-c('species','species.age')
  cat('writing output table','\n')
  write.table(output.table,file=paste(path,name,'species_ages.txt',sep=''),row.names=F,quote=F,sep='\t')
}


generate_grid_medianDR<-function(dr.table,species.grid.tablefile,cellgrid.tablefile,path,name){
  species.grid.table<-read.table(species.grid.tablefile,header=T,sep='\t')
  cell.grid.table<-read.table(cellgrid.tablefile,header=T,sep='\t')
  #add a space at the beginning and end of the cell string for grep to work
  species.grid.table$cells<-paste(' ',species.grid.table$cells,' ')
  dr.table<-read.table(dr.table,header=T,sep='\t')
  colnames(dr.table)[2]<-'DR'
  #add a space at the beginning and end of the cell string for grep to work
  species.grid.table<-merge(species.grid.table,dr.table,by.x='spp',by.y='Species',all.x=TRUE)
  cell.species.matrix<-matrix(NA,ncol=2,nrow=nrow(cell.grid.table))
  for (i in 1:nrow(cell.grid.table)){
    cat(i,'\n')
    
    species.in.cell<-as.character(species.grid.table[grep(paste(' ',i,' ',sep=''),species.grid.table$cells),]$spp)
    cell.species.matrix[i,]<-cbind(i,median(species.grid.table[species.grid.table$spp%in%species.in.cell,'DR'],na.rm=TRUE))
    
  }
  cell.species.df<-as.data.frame(cell.species.matrix,stringsAsFactors = F)
  colnames(cell.species.df)<-c('cells','median.DR')
  write.table(cell.species.df,file=paste(path,'/',name,'_medianDR_grid_table.txt',sep=''),sep='\t',row.names=F,quote=F)
}

generate_grid_medianlambda<-function(BAMM.table,species.grid.tablefile,cellgrid.tablefile,path,name){
  species.grid.table<-read.table(species.grid.tablefile,header=T,sep='\t')
  cell.grid.table<-read.table(cellgrid.tablefile,header=T,sep='\t')
  #add a space at the beginning and end of the cell string for grep to work
  species.grid.table$cells<-paste(' ',species.grid.table$cells,' ')
  BAMM.table<-read.table(BAMM.table,header=T,sep='\t')
  colnames(BAMM.table)[2]<-'lambda.avg'
  #add a space at the beginning and end of the cell string for grep to work
  species.grid.table<-merge(species.grid.table,BAMM.table,by.x='spp',by.y='spp',all.x=TRUE)
  cell.species.matrix<-matrix(NA,ncol=2,nrow=nrow(cell.grid.table))
  for (i in 1:nrow(cell.grid.table)){
    cat(i,'\n')
    
    species.in.cell<-as.character(species.grid.table[grep(paste(' ',i,' ',sep=''),species.grid.table$cells),]$spp)
    cell.species.matrix[i,]<-cbind(i,median(species.grid.table[species.grid.table$spp%in%species.in.cell,'lambda.avg'],na.rm=TRUE))
    
  }
  cell.species.df<-as.data.frame(cell.species.matrix,stringsAsFactors = F)
  colnames(cell.species.df)<-c('cells','median.lambda')
  write.table(cell.species.df,file=paste(path,'/',name,'_medianlambda_grid_table.txt',sep=''),sep='\t',row.names=F,quote=F)
}


generate_grid_medianage<-function(age.table,species.grid.tablefile,cellgrid.tablefile,path,name){
  species.grid.table<-read.table(species.grid.tablefile,header=T,sep='\t')
  cell.grid.table<-read.table(cellgrid.tablefile,header=T,sep='\t')
  #add a space at the beginning and end of the cell string for grep to work
  species.grid.table$cells<-paste(' ',species.grid.table$cells,' ')
  age.table<-read.table(age.table,header=T,sep='\t')
  colnames(age.table)[2]<-'species.age'
  #add a space at the beginning and end of the cell string for grep to work
  species.grid.table<-merge(species.grid.table,age.table,by.x='spp',by.y='species',all.x=TRUE)
  cell.species.matrix<-matrix(NA,ncol=2,nrow=nrow(cell.grid.table))
  for (i in 1:nrow(cell.grid.table)){
    cat(i,'\n')
    
    species.in.cell<-as.character(species.grid.table[grep(paste(' ',i,' ',sep=''),species.grid.table$cells),]$spp)
    cell.species.matrix[i,]<-cbind(i,median(species.grid.table[species.grid.table$spp%in%species.in.cell,'species.age'],na.rm=TRUE))
    
  }
  cell.species.df<-as.data.frame(cell.species.matrix,stringsAsFactors = F)
  colnames(cell.species.df)<-c('cells','median.age')
  write.table(cell.species.df,file=paste(path,'/',name,'_medianage_grid_table.txt',sep=''),sep='\t',row.names=F,quote=F)
}


############
correlation_residualsDR_BAMM<-function(DRresidualstable.file,BAMMresidualstable.file){
  DRresidualstable<-read.table(DRresidualstable.file,header=T,sep='\t',stringsAsFactors = F)
  colnames(DRresidualstable)<-c('cells','Q1.DR','Q2.DR','Q3.DR','Q4.DR')
  BAMMresidualstable<-read.table(BAMMresidualstable.file,header=T,sep='\t',stringsAsFactors = F)
  colnames(BAMMresidualstable)<-c('cells','Q1.BAMM','Q2.BAMM','Q3.BAMM','Q4.BAMM')
  DRBAMMresiduals<-merge(DRresidualstable,BAMMresidualstable,by='cells')
  corrQ1<-cor.test(DRBAMMresiduals$Q1.DR,DRBAMMresiduals$Q1.BAMM)
  #plot(DRBAMMresiduals$Q1.DR,DRBAMMresiduals$Q1.BAMM,xlab='Q1.DRresiduals',ylab='Q1.BAMMresiduals',pch=16,cex=.5,main=paste('BAMM-DR residuals correlation (r=',round(corrQ1$estimate,3),')',sep=''),cex.main=.8)
  plot(DRBAMMresiduals$Q1.DR,DRBAMMresiduals$Q1.BAMM,xlab='Q1.DRresiduals',ylab='Q1.BAMMresiduals',pch=16,cex=.5,main=paste('BAMM-DR residuals correlation (r=',round(corrQ1$estimate,3),')',sep=''),cex.main=.8,col="#00000033")
  #smoothScatter(DRBAMMresiduals$Q1.DR,DRBAMMresiduals$Q1.BAMM,xlab='Q1.DRresiduals',ylab='Q1.BAMMresiduals',pch=16,cex=.5,main=paste('BAMM-DR residuals correlation (r=',round(corrQ1$estimate,3),')',sep=''),cex.main=.8)
  corrQ4<-cor.test(DRBAMMresiduals$Q4.DR,DRBAMMresiduals$Q4.BAMM)
  #plot(DRBAMMresiduals$Q4.DR,DRBAMMresiduals$Q4.BAMM,xlab='Q4.DRresiduals',ylab='Q4.BAMMresiduals',pch=16,cex=.5,main=paste('BAMM-DR residuals correlation (r=',round(corrQ4$estimate,3),')',sep=''),cex.main=.8)
  plot(DRBAMMresiduals$Q4.DR,DRBAMMresiduals$Q4.BAMM,xlab='Q4.DRresiduals',ylab='Q4.BAMMresiduals',pch=16,cex=.5,main=paste('BAMM-DR residuals correlation (r=',round(corrQ4$estimate,3),')',sep=''),cex.main=.8,col="#00000033")
  #smoothScatter(DRBAMMresiduals$Q4.DR,DRBAMMresiduals$Q4.BAMM,xlab='Q4.DRresiduals',ylab='Q4.BAMMresiduals',pch=16,cex=.5,main=paste('BAMM-DR residuals correlation (r=',round(corrQ4$estimate,3),')',sep=''),cex.main=.8)
  #plot(log10(DRBAMMrates$DR),log10(DRBAMMrates$netdiv.avg),xlab='log10.DR',ylab='log10.BAMM.netdiv',pch=16,cex=.5,main=paste('BAMM-DR correlation (r=',round(corr$estimate,3),')',sep=''),cex.main=.8)
}