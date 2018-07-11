#run model in hydrogen cluster
.libPaths('/home/ecosys/ji247/R/x86_64-pc-linux-gnu-library/3.3')
setwd('/home/ecosys/ji247/hotspots_vertebrates/')

source('./R/process_BSM_events.R')
#control files rates for 2my
speciation.rate.treefile='./output/mammals/trees/neotrop_speciationrates.tree'
speciation.rate.tree<-read.tree(speciation.rate.treefile)
table.speciation.rate.neotrop<-prt(speciation.rate.tree,printflag = FALSE)
control.polygons.BSM<-list.files('./output/mammals/immigration/controls/neotrop/neotrop.BSM/',pattern='_BSM.RDS')
results.BSM.control.list<-lapply(control.polygons.BSM,function(x) readRDS(paste('./output/mammals/immigration/controls/neotrop/neotrop.BSM/',x,sep='')))
BSMoutput.control.list<-lapply(results.BSM.control.list,function(x) get_BSM_BAMM_fullrates_events_2my_object(results.BSM=x[[2]],name='neotrop',speciation.rate.treefile=table.speciation.rate.neotrop))
saveRDS(BSMoutput.control.list,file='./output/mammals/immigration/controls/neotrop/neotrop_BSMcontrol_output_2my_full.RDS')

source('./R/process_BSM_events.R')
#control files rates for 2my
speciation.rate.treefile='./output/birds/trees/neotrop_speciationrates.tree'
speciation.rate.tree<-read.tree(speciation.rate.treefile)
table.speciation.rate.neotrop<-prt(speciation.rate.tree,printflag = FALSE)
control.polygons.BSM<-list.files('./output/birds/immigration/controls/neotrop/',pattern='_BSM.RDS')
results.BSM.control.list<-lapply(control.polygons.BSM,function(x) readRDS(paste('./output/birds/immigration/controls/neotrop/',x,sep='')))
BSMoutput.control.list<-lapply(results.BSM.control.list,function(x) get_BSM_BAMM_fullrates_events_2my_object(results.BSM=x[[2]],name='neotrop',speciation.rate.treefile=table.speciation.rate.neotrop))
saveRDS(BSMoutput.control.list,file='./output/birds/immigration/controls/neotrop/neotrop_BSMcontrol_output_2my_full.RDS')
