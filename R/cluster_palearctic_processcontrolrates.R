#run model in hydrogen cluster
.libPaths('/home/ecosys/ji247/R/x86_64-pc-linux-gnu-library/3.3')
setwd('/home/ecosys/ji247/hotspots_vertebrates/')

source('./R/process_BSM_events.R')
#control files rates for 2my
speciation.rate.treefile='./output/mammals/trees/palearctic_speciationrates.tree'
speciation.rate.tree<-read.tree(speciation.rate.treefile)
table.speciation.rate.palearctic<-prt(speciation.rate.tree,printflag = FALSE)
control.polygons.BSM<-list.files('./output/mammals/immigration/controls/palearctic/palearctic.BSM/',pattern='_BSM.RDS')
results.BSM.control.list<-lapply(control.polygons.BSM,function(x) readRDS(paste('./output/mammals/immigration/controls/palearctic/palearctic.BSM/',x,sep='')))
BSMoutput.control.list<-lapply(results.BSM.control.list,function(x) get_BSM_BAMM_fullrates_events_2my_object(results.BSM=x[[2]],name='palearctic',speciation.rate.treefile=table.speciation.rate.palearctic))
saveRDS(BSMoutput.control.list,file='./output/mammals/immigration/controls/palearctic/palearctic_BSMcontrol_output_2my_full.RDS')

source('./R/process_BSM_events.R')
#control files rates for 2my
speciation.rate.treefile='./output/birds/trees/palearctic_speciationrates.tree'
speciation.rate.tree<-read.tree(speciation.rate.treefile)
table.speciation.rate.palearctic<-prt(speciation.rate.tree,printflag = FALSE)
control.polygons.BSM<-list.files('./output/birds/immigration/controls/palearctic/palearctic.BSM/',pattern='_BSM.RDS')
results.BSM.control.list<-lapply(control.polygons.BSM,function(x) readRDS(paste('./output/birds/immigration/controls/palearctic/palearctic.BSM/',x,sep='')))
BSMoutput.control.list<-lapply(results.BSM.control.list,function(x) get_BSM_BAMM_fullrates_events_2my_object(results.BSM=x[[2]],name='palearctic',speciation.rate.treefile=table.speciation.rate.palearctic))
saveRDS(BSMoutput.control.list,file='./output/birds/immigration/controls/palearctic/palearctic_BSMcontrol_output_2my_full.RDS')
