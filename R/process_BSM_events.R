library(grDevices)
library(optimx)         # You need to have some version of optimx available
# as it is a BioGeoBEARS dependency; however, if you
# don't want to use optimx, and use optim() (from R core) 
# you can set:
# BioGeoBEARS_run_object$use_optimx = FALSE
# ...everything should work either way -- NJM 2014-01-08
library(FD)       # for FD::maxent() (make sure this is up-to-date)
library(snow)     # (if you want to use multicore functionality; some systems/R versions prefer library(parallel), try either)
library(parallel)
library(BioGeoBEARS)
source("http://phylo.wdfiles.com/local--files/biogeobears/cladoRcpp.R") # (needed now that traits model added; source FIRST!)
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_add_fossils_randomly_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_basics_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_calc_transition_matrices_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_classes_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_detection_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_DNA_cladogenesis_sim_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_extract_Qmat_COOmat_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_generics_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_models_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_on_multiple_trees_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_plots_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_readwrite_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_simulate_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_SSEsim_makePlots_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_SSEsim_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_stochastic_mapping_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_stratified_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_univ_model_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/calc_uppass_probs_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/calc_loglike_sp_v01.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/get_stratified_subbranch_top_downpass_likelihoods_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/runBSM_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/stochastic_map_given_inputs.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/summarize_BSM_tables_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_traits_v1.R") # added traits model
calc_loglike_sp = compiler::cmpfun(calc_loglike_sp_prebyte)    # crucial to fix bug in uppass calculations
calc_independent_likelihoods_on_each_branch = compiler::cmpfun(calc_independent_likelihoods_on_each_branch_prebyte)
# slight speedup hopefully
library(ape)
####TO DO: delete "old functions"####
####TO DO: double check notation of events is correct in scripts####





plot_all_cumulative_events_BSM_logCI<-function(results.BSM,name,age){
  events.list<-list()
  for (a in 1:length(results.BSM$RES_clado_events_tables)){
      clado.df<-results.BSM$RES_clado_events_tables[[a]]
      clado.df<-clado.df[,c('time_bp','clado_event_txt')]
      clado.df<-clado.df[clado.df$clado_event_txt!='',]
      ana.df<-results.BSM$RES_ana_events_tables[[a]]
      ana.df<-ana.df[,c('abs_event_time','event_txt')]
      ana.df<-ana.df[ana.df$abs_event_time!='',]
      colnames(clado.df)<-c('time_bp','event_txt')
      colnames(ana.df)<-c('time_bp','event_txt')
      events.list[[a]]<-rbind(clado.df,ana.df)
    
  }
  
  event.type.table<-lapply(events.list,function(x)table(x$event_txt))
  
  
  if(name=='afrotrop'){
    realm.code<-'A'
  }else if(name=='austral'){
    realm.code<-'B'
  }else if(name=='indo'){
    realm.code<-'C'
  }else if(name=='nearctic'){
    realm.code<-'D'
  }else if(name=='neotrop'){
    realm.code<-'E'
  }else if(name=='palearctic'){
    realm.code<-'F'
  }
  hot.code<-'G'
  event.names<-unlist(lapply(event.type.table,function(x)unlist(x)))
  event.names<-names(table(names(event.names)))
  #anything that involves the realm
  realm.event.names<-event.names[grep(realm.code,event.names)]
  #anything that involves the hotspot
  hot.event.names<-event.names[grep(hot.code,event.names)]
  #get insitu cladogenesis in realm
  insitu.realm.event.names<-realm.event.names[grep(paste('->.*',realm.code,'.*,*',sep=''),realm.event.names)]
  insitu.realm.event.names<-insitu.realm.event.names[grep(paste(realm.code,'.*->',sep=''),insitu.realm.event.names)]
  insitu.realm.event.names<-insitu.realm.event.names[grep(paste('.*->.*',realm.code,'.*,.*',realm.code,sep=''),insitu.realm.event.names)]
  #get strict insitu cladogenesis in realm
  #insitu.strict.realm.event.names<-insitu.realm.event.names[grep(paste('^',realm.code,'->',realm.code,',',realm.code,'$',sep=''),insitu.realm.event.names)]
  #get insitu cladogenesis in hotspot
  insitu.hot.event.names<-hot.event.names[grep(paste('->.*',hot.code,'.*,*',sep=''),hot.event.names)]
  insitu.hot.event.names<-insitu.hot.event.names[grep(paste(hot.code,'.*->',sep=''),insitu.hot.event.names)]
  insitu.hot.event.names<-insitu.hot.event.names[grep(paste('.*->.*',hot.code,'.*,.*',hot.code,sep=''),insitu.hot.event.names)]
  #get strict insitu cladogenesis in realm
  #insitu.strict.hot.event.names<-insitu.hot.event.names[grep(paste('^',hot.code,'->',hot.code,',',hot.code,'$',sep=''),insitu.hot.event.names)]
  
  #get dispersal from hot into realm
  hot.into.realm.event.names<-realm.event.names[-grep(paste(realm.code,'.*->',sep=''),realm.event.names)]
  hot.into.realm.event.names<-hot.into.realm.event.names[grep(paste(hot.code,'.*->.*',sep=''),hot.into.realm.event.names)]
  #get dispersal from realm into hot
  realm.into.hot.event.names<-hot.event.names[-grep(paste(hot.code,'.*->',sep=''),hot.event.names)]
  realm.into.hot.event.names<-realm.into.hot.event.names[grep(paste(realm.code,'.*->.*',sep=''),realm.into.hot.event.names)]
  #get local extinctions in hot
  hot.local.extinction.event.names<-event.names[grep(paste('.*',hot.code,'.*->*',sep=''),event.names)]
  hot.local.extinction.event.names<-hot.local.extinction.event.names[-grep(paste('.*->.*',hot.code,'.*',sep=''),hot.local.extinction.event.names)]
  #get local extinctions in realm
  realm.local.extinction.event.names<-event.names[grep(paste('.*',realm.code,'.*->.*',sep=''),event.names)]
  realm.local.extinction.event.names<-realm.local.extinction.event.names[-grep(paste('.*->.*',realm.code,'.*',sep=''),realm.local.extinction.event.names)]
  #get dispersal into realm
  dispersal.into.realm.event.names<-realm.event.names[-grep(paste(realm.code,'.*->',sep=''),realm.event.names)]
  #get dispersal into hotspot
  dispersal.into.hot.event.names<-hot.event.names[-grep(paste(hot.code,'.*->',sep=''),hot.event.names)]
  #get dispersal from hotspot
  dispersal.from.hot.event.names<-event.names[grep(paste('.*',hot.code,'.*->*',sep=''),event.names)]
  dispersal.from.hot.event.names<-dispersal.from.hot.event.names[!(dispersal.from.hot.event.names%in%insitu.hot.event.names)]
  dispersal.from.hot.event.names<-dispersal.from.hot.event.names[!(dispersal.from.hot.event.names%in%hot.local.extinction.event.names)]
  if(length(grep(paste('->',hot.code,'$',sep=''),dispersal.from.hot.event.names))>0){
    dispersal.from.hot.event.names<-dispersal.from.hot.event.names[-grep(paste('->',hot.code,'$',sep=''),dispersal.from.hot.event.names)]
  }
  
  #get dispersal from realm
  dispersal.from.realm.event.names<-event.names[grep(paste('.*',realm.code,'.*->*',sep=''),event.names)]
  dispersal.from.realm.event.names<-dispersal.from.realm.event.names[!(dispersal.from.realm.event.names%in%insitu.realm.event.names)]
  dispersal.from.realm.event.names<-dispersal.from.realm.event.names[!(dispersal.from.realm.event.names%in%realm.local.extinction.event.names)]
  if(length(grep(paste('->',realm.code,'$',sep=''),dispersal.from.realm.event.names))>0){
    dispersal.from.realm.event.names<-dispersal.from.realm.event.names[-grep(paste('->',realm.code,'$',sep=''),dispersal.from.realm.event.names)]
  }
  
  
  
  dispersal.into.realm<-unlist(lapply(event.type.table,function(x)sum(x[dispersal.into.realm.event.names],na.rm=TRUE)))
  dispersal.into.hot<-unlist(lapply(event.type.table,function(x)sum(x[dispersal.into.hot.event.names],na.rm=TRUE)))
  dispersal.from.realm<-unlist(lapply(event.type.table,function(x)sum(x[dispersal.from.realm.event.names],na.rm=TRUE)))
  dispersal.from.hot<-unlist(lapply(event.type.table,function(x)sum(x[dispersal.from.hot.event.names],na.rm=TRUE)))
  dispersal.hot.into.realm<-unlist(lapply(event.type.table,function(x)sum(x[hot.into.realm.event.names],na.rm=TRUE)))
  dispersal.realm.into.hot<-unlist(lapply(event.type.table,function(x)sum(x[realm.into.hot.event.names],na.rm=TRUE)))
  #insitu.strict.realm<-unlist(lapply(event.type.table,function(x)sum(x[insitu.strict.realm.event.names],na.rm=TRUE)))
  #insitu.strict.hot<-unlist(lapply(event.type.table,function(x)sum(x[insitu.strict.hot.event.names],na.rm=TRUE)))
  insitu.realm<-unlist(lapply(event.type.table,function(x)sum(x[insitu.realm.event.names],na.rm=TRUE)))
  insitu.hot<-unlist(lapply(event.type.table,function(x)sum(x[insitu.hot.event.names],na.rm=TRUE)))
  hot.local.extinction<-unlist(lapply(event.type.table,function(x)sum(x[hot.local.extinction.event.names],na.rm=TRUE)))
  realm.local.extinction<-unlist(lapply(event.type.table,function(x)sum(x[realm.local.extinction.event.names],na.rm=TRUE)))
  
  dispersal.from.realm<-lapply(dispersal.from.realm,function(x) replace(x,is.na(x),0))
  dispersal.from.hot<-lapply(dispersal.from.hot,function(x) replace(x,is.na(x),0))
  dispersal.into.realm<-lapply(dispersal.into.realm,function(x) replace(x,is.na(x),0))
  dispersal.into.hot<-lapply(dispersal.into.hot,function(x) replace(x,is.na(x),0))
  dispersal.hot.into.realm<-lapply(dispersal.hot.into.realm,function(x) replace(x,is.na(x),0))
  dispersal.realm.into.hot<-lapply(dispersal.realm.into.hot,function(x) replace(x,is.na(x),0))
  #insitu.strict.realm<-lapply(insitu.strict.realm,function(x) replace(x,is.na(x),0))
  #insitu.strict.hot<-lapply(insitu.strict.hot,function(x) replace(x,is.na(x),0))
  insitu.realm<-lapply(insitu.realm,function(x) replace(x,is.na(x),0))
  insitu.hot<-lapply(insitu.hot,function(x) replace(x,is.na(x),0))
  hot.local.extinction<-lapply(hot.local.extinction,function(x) replace(x,is.na(x),0))
  realm.local.extinction<-lapply(realm.local.extinction,function(x) replace(x,is.na(x),0))
  
  ###########cumulative sum of events plot
  #process a results.BSM to extract events and times in a list of lists (each element of the list is one clade with a list of 50 replicate BSMs)
  events.list<-list()
  for (a in 1:length(results.BSM$RES_clado_events_tables)){
    clado.df<-results.BSM$RES_clado_events_tables[[a]]
    clado.df<-clado.df[,c('time_bp','clado_event_txt')]
    clado.df<-clado.df[clado.df$clado_event_txt!='',]
    ana.df<-results.BSM$RES_ana_events_tables[[a]]
    ana.df<-ana.df[,c('abs_event_time','event_txt')]
    ana.df<-ana.df[ana.df$abs_event_time!='',]
    colnames(clado.df)<-c('time_bp','event_txt')
    colnames(ana.df)<-c('time_bp','event_txt')
    events.list[[a]]<-rbind(clado.df,ana.df)
  }
  
  #create a df that bins events into 1 myr slots
  #get the roots of all clades
  clade.roots<-results.BSM$RES_clado_events_tables[[1]][results.BSM$RES_clado_events_tables[[1]]$node.type=='root',]$time_bp
  #another way of binning events, a list of lists (binning events for each BSM replicate)
  binned.events.list<-vector('list',length=length(events.list))
  cat('binning events','\n')
  for (i in 1:length(events.list)){
    # cat(i,'\n')
    binned.events.list[[i]]<-vector('list',length=length(seq(0,floor(clade.roots)))+1)
    for (b in 1:nrow(events.list[[i]])){
      # cat(b,'\n')
      time_bp<-floor(events.list[[i]][b,'time_bp'])
      event<-events.list[[i]][b,'event_txt']
      names(event)<-row.names(events.list[[i]][b,])
      if(is.null((binned.events.list[[i]][[time_bp+1]]))){
        binned.events.list[[i]][[time_bp+1]]<-c('0',event)
      }else{
        binned.events.list[[i]][[time_bp+1]]<-c(binned.events.list[[i]][[time_bp+1]],event)  
      }
    }
  }
  binned.events.list.from.realm<-list()
  binned.events.list.from.hot<-list()
  binned.events.list.into.realm<-list()
  binned.events.list.into.hot<-list()
  binned.events.list.hot.into.realm<-list()
  binned.events.list.realm.into.hot<-list()
  binned.events.list.insitu.realm<-list()
  binned.events.list.insitu.hot<-list()
  binned.events.list.hot.local.extinction<-list()
  binned.events.list.realm.local.extinction<-list()
  
  
  
  cat('subsetting binned events to realm','\n')
  for(i in 1:length(binned.events.list)){
    #cat(i,'\n')
    binned.events.list.from.realm[[i]]<-vector('list',length=length(seq(0,floor(clade.roots)))+1)
    binned.events.list.from.hot[[i]]<-vector('list',length=length(seq(0,floor(clade.roots)))+1)
    binned.events.list.into.realm[[i]]<-vector('list',length=length(seq(0,floor(clade.roots)))+1)
    binned.events.list.into.hot[[i]]<-vector('list',length=length(seq(0,floor(clade.roots)))+1)
    binned.events.list.hot.into.realm[[i]]<-vector('list',length=length(seq(0,floor(clade.roots)))+1)
    binned.events.list.realm.into.hot[[i]]<-vector('list',length=length(seq(0,floor(clade.roots)))+1)
    binned.events.list.insitu.realm[[i]]<-vector('list',length=length(seq(0,floor(clade.roots)))+1)
    binned.events.list.insitu.hot[[i]]<-vector('list',length=length(seq(0,floor(clade.roots)))+1)
    binned.events.list.hot.local.extinction[[i]]<-vector('list',length=length(seq(0,floor(clade.roots)))+1)
    binned.events.list.realm.local.extinction[[i]]<-vector('list',length=length(seq(0,floor(clade.roots)))+1)
    
    
    for (a in 1:length(binned.events.list[[i]])){
      # cat(a,'\n')
      binned.events.list.from.realm[[i]][[a]]<-binned.events.list[[i]][[a]][unlist(sapply(dispersal.from.realm.event.names,function(x) grep(paste('^',x,'$',sep=''),binned.events.list[[i]][[a]])))]
      binned.events.list.from.hot[[i]][[a]]<-binned.events.list[[i]][[a]][unlist(sapply(dispersal.from.hot.event.names,function(x)grep(paste('^',x,'$',sep=''),binned.events.list[[i]][[a]])))]
      binned.events.list.into.realm[[i]][[a]]<-binned.events.list[[i]][[a]][unlist(sapply(dispersal.into.realm.event.names,function(x) grep(paste('^',x,'$',sep=''),binned.events.list[[i]][[a]])))]
      binned.events.list.into.hot[[i]][[a]]<-binned.events.list[[i]][[a]][unlist(sapply(dispersal.into.hot.event.names,function(x)grep(paste('^',x,'$',sep=''),binned.events.list[[i]][[a]])))]
      binned.events.list.hot.into.realm[[i]][[a]]<-binned.events.list[[i]][[a]][unlist(sapply(hot.into.realm.event.names,function(x) grep(paste('^',x,'$',sep=''),binned.events.list[[i]][[a]])))]
      binned.events.list.realm.into.hot[[i]][[a]]<-binned.events.list[[i]][[a]][unlist(sapply(realm.into.hot.event.names,function(x)grep(paste('^',x,'$',sep=''),binned.events.list[[i]][[a]])))]
      binned.events.list.insitu.realm[[i]][[a]]<-binned.events.list[[i]][[a]][unlist(sapply(insitu.realm.event.names,function(x)grep(paste('^',x,'$',sep=''),binned.events.list[[i]][[a]])))]
      binned.events.list.insitu.hot[[i]][[a]]<-binned.events.list[[i]][[a]][unlist(sapply(insitu.hot.event.names,function(x)grep(paste('^',x,'$',sep=''),binned.events.list[[i]][[a]])))]
      binned.events.list.hot.local.extinction[[i]][[a]]<-binned.events.list[[i]][[a]][unlist(sapply(hot.local.extinction.event.names,function(x) grep(paste('^',x,'$',sep=''),binned.events.list[[i]][[a]])))]
      binned.events.list.realm.local.extinction[[i]][[a]]<-binned.events.list[[i]][[a]][unlist(sapply(realm.local.extinction.event.names,function(x)grep(paste('^',x,'$',sep=''),binned.events.list[[i]][[a]])))]
      
    }
    
    
    
    
  }
  binned.events.list.from.realm.lengths<-lapply(binned.events.list.from.realm,function(x)unlist(lapply(x,function(x) length(x))))
  binned.events.list.from.hot.lengths<-lapply(binned.events.list.from.hot,function(x)unlist(lapply(x,function(x) length(x))))
  binned.events.list.into.realm.lengths<-lapply(binned.events.list.into.realm,function(x)unlist(lapply(x,function(x) length(x))))
  binned.events.list.into.hot.lengths<-lapply(binned.events.list.into.hot,function(x)unlist(lapply(x,function(x) length(x))))
  binned.events.list.hot.into.realm.lengths<-lapply(binned.events.list.hot.into.realm,function(x)unlist(lapply(x,function(x) length(x))))
  binned.events.list.realm.into.hot.lengths<-lapply(binned.events.list.realm.into.hot,function(x)unlist(lapply(x,function(x) length(x))))
  binned.events.list.insitu.realm.lengths<-lapply(binned.events.list.insitu.realm,function(x)unlist(lapply(x,function(x) length(x))))
  binned.events.list.insitu.hot.lengths<-lapply(binned.events.list.insitu.hot,function(x)unlist(lapply(x,function(x) length(x))))
  binned.events.list.hot.local.extinction.lengths<-lapply(binned.events.list.hot.local.extinction,function(x)unlist(lapply(x,function(x) length(x))))
  binned.events.list.realm.local.extinction.lengths<-lapply(binned.events.list.realm.local.extinction,function(x)unlist(lapply(x,function(x) length(x))))
  
  
  binned.events.list.from.realm.lengths.rev<-lapply(binned.events.list.from.realm.lengths,function(x) rev(x))
  binned.events.list.from.hot.lengths.rev<-lapply(binned.events.list.from.hot.lengths,function(x) rev(x))
  binned.events.list.into.realm.lengths.rev<-lapply(binned.events.list.into.realm.lengths,function(x) rev(x))
  binned.events.list.into.hot.lengths.rev<-lapply(binned.events.list.into.hot.lengths,function(x) rev(x))
  binned.events.list.into.realm.lengths.rev<-lapply(binned.events.list.into.realm.lengths,function(x) rev(x))
  binned.events.list.into.hot.lengths.rev<-lapply(binned.events.list.into.hot.lengths,function(x) rev(x))
  binned.events.list.hot.into.realm.lengths.rev<-lapply(binned.events.list.hot.into.realm.lengths,function(x) rev(x))
  binned.events.list.realm.into.hot.lengths.rev<-lapply(binned.events.list.realm.into.hot.lengths,function(x) rev(x))
  binned.events.list.insitu.realm.lengths.rev<-lapply(binned.events.list.insitu.realm.lengths,function(x) rev(x))
  binned.events.list.insitu.hot.lengths.rev<-lapply(binned.events.list.insitu.hot.lengths,function(x) rev(x))
  binned.events.list.hot.local.extinction.lengths.rev<-lapply(binned.events.list.hot.local.extinction.lengths,function(x) rev(x))
  binned.events.list.realm.local.extinction.lengths.rev<-lapply(binned.events.list.realm.local.extinction.lengths,function(x) rev(x))
  
  
  
  binned.events.list.from.realm.lengths.rev.cumsum<-lapply(binned.events.list.from.realm.lengths.rev,function(x) cumsum(x))
  binned.events.list.from.hot.lengths.rev.cumsum<-lapply(binned.events.list.from.hot.lengths.rev,function(x) cumsum(x))
  binned.events.list.into.realm.lengths.rev.cumsum<-lapply(binned.events.list.into.realm.lengths.rev,function(x) cumsum(x))
  binned.events.list.into.hot.lengths.rev.cumsum<-lapply(binned.events.list.into.hot.lengths.rev,function(x) cumsum(x))
  binned.events.list.hot.into.realm.lengths.rev.cumsum<-lapply(binned.events.list.hot.into.realm.lengths.rev,function(x) cumsum(x))
  binned.events.list.realm.into.hot.lengths.rev.cumsum<-lapply(binned.events.list.realm.into.hot.lengths.rev,function(x) cumsum(x))
  binned.events.list.insitu.realm.lengths.rev.cumsum<-lapply(binned.events.list.insitu.realm.lengths.rev,function(x) cumsum(x))
  binned.events.list.insitu.hot.lengths.rev.cumsum<-lapply(binned.events.list.insitu.hot.lengths.rev,function(x) cumsum(x))
  binned.events.list.hot.local.extinction.lengths.rev.cumsum<-lapply(binned.events.list.hot.local.extinction.lengths.rev,function(x) cumsum(x))
  binned.events.list.realm.local.extinction.lengths.rev.cumsum<-lapply(binned.events.list.realm.local.extinction.lengths.rev,function(x) cumsum(x))
  
  
  binned.events.list.from.realm.lengths.rev.cumsum.CI<-list()
  binned.events.list.from.hot.lengths.rev.cumsum.CI<-list()
  binned.events.list.into.realm.lengths.rev.cumsum.CI<-list()
  binned.events.list.into.hot.lengths.rev.cumsum.CI<-list()
  binned.events.list.hot.into.realm.lengths.rev.cumsum.CI<-list()
  binned.events.list.realm.into.hot.lengths.rev.cumsum.CI<-list()
  binned.events.list.insitu.realm.lengths.rev.cumsum.CI<-list()
  binned.events.list.insitu.hot.lengths.rev.cumsum.CI<-list()
  binned.events.list.hot.local.extinction.lengths.rev.cumsum.CI<-list()
  binned.events.list.realm.local.extinction.lengths.rev.cumsum.CI<-list()
  
  
  for (i in 1:length(binned.events.list.into.realm.lengths.rev.cumsum[[1]])){
    binned.events.list.from.realm.lengths.rev.cumsum.CI[[i]]<-sapply(binned.events.list.from.realm.lengths.rev.cumsum,"[[",i)
    binned.events.list.from.hot.lengths.rev.cumsum.CI[[i]]<-sapply(binned.events.list.from.hot.lengths.rev.cumsum,"[[",i)
    binned.events.list.into.realm.lengths.rev.cumsum.CI[[i]]<-sapply(binned.events.list.into.realm.lengths.rev.cumsum,"[[",i)
    binned.events.list.into.hot.lengths.rev.cumsum.CI[[i]]<-sapply(binned.events.list.into.hot.lengths.rev.cumsum,"[[",i)
    binned.events.list.hot.into.realm.lengths.rev.cumsum.CI[[i]]<-sapply(binned.events.list.hot.into.realm.lengths.rev.cumsum,"[[",i)
    binned.events.list.realm.into.hot.lengths.rev.cumsum.CI[[i]]<-sapply(binned.events.list.realm.into.hot.lengths.rev.cumsum,"[[",i)
    binned.events.list.insitu.realm.lengths.rev.cumsum.CI[[i]]<-sapply(binned.events.list.insitu.realm.lengths.rev.cumsum,"[[",i)
    binned.events.list.insitu.hot.lengths.rev.cumsum.CI[[i]]<-sapply(binned.events.list.insitu.hot.lengths.rev.cumsum,"[[",i)
    binned.events.list.hot.local.extinction.lengths.rev.cumsum.CI[[i]]<-sapply(binned.events.list.hot.local.extinction.lengths.rev.cumsum,"[[",i)
    binned.events.list.realm.local.extinction.lengths.rev.cumsum.CI[[i]]<-sapply(binned.events.list.realm.local.extinction.lengths.rev.cumsum,"[[",i)
  }
  
  binned.events.list.from.realm.lengths.rev.cumsum.CI<-lapply(binned.events.list.from.realm.lengths.rev.cumsum.CI,function(x)quantile(x,c(0.05,0.95,0.50)))
  binned.events.list.from.hot.lengths.rev.cumsum.CI<-lapply(binned.events.list.from.hot.lengths.rev.cumsum.CI,function(x)quantile(x,c(0.05,0.95,0.50)))
  binned.events.list.into.realm.lengths.rev.cumsum.CI<-lapply(binned.events.list.into.realm.lengths.rev.cumsum.CI,function(x)quantile(x,c(0.05,0.95,0.50)))
  binned.events.list.into.hot.lengths.rev.cumsum.CI<-lapply(binned.events.list.into.hot.lengths.rev.cumsum.CI,function(x)quantile(x,c(0.05,0.95,0.50)))
  binned.events.list.hot.into.realm.lengths.rev.cumsum.CI<-lapply(binned.events.list.hot.into.realm.lengths.rev.cumsum.CI,function(x)quantile(x,c(0.05,0.95,0.50)))
  binned.events.list.realm.into.hot.lengths.rev.cumsum.CI<-lapply(binned.events.list.realm.into.hot.lengths.rev.cumsum.CI,function(x)quantile(x,c(0.05,0.95,0.50)))
  binned.events.list.insitu.realm.lengths.rev.cumsum.CI<-lapply(binned.events.list.insitu.realm.lengths.rev.cumsum.CI,function(x)quantile(x,c(0.05,0.95,0.50)))
  binned.events.list.insitu.hot.lengths.rev.cumsum.CI<-lapply(binned.events.list.insitu.hot.lengths.rev.cumsum.CI,function(x)quantile(x,c(0.05,0.95,0.50)))
  binned.events.list.hot.local.extinction.lengths.rev.cumsum.CI<-lapply(binned.events.list.hot.local.extinction.lengths.rev.cumsum.CI,function(x)quantile(x,c(0.05,0.95,0.50)))
  binned.events.list.realm.local.extinction.lengths.rev.cumsum.CI<-lapply(binned.events.list.realm.local.extinction.lengths.rev.cumsum.CI,function(x)quantile(x,c(0.05,0.95,0.50)))
  
  
  
  
  binned.events.list.from.realm.lengths.rev.cumsum.CI<-binned.events.list.from.realm.lengths.rev.cumsum.CI[c((length(binned.events.list.from.realm.lengths.rev.cumsum.CI)-age):length(binned.events.list.from.realm.lengths.rev.cumsum.CI))]
  binned.events.list.from.hot.lengths.rev.cumsum.CI<-binned.events.list.from.hot.lengths.rev.cumsum.CI[c((length(binned.events.list.from.hot.lengths.rev.cumsum.CI)-age):length(binned.events.list.from.hot.lengths.rev.cumsum.CI))]
  binned.events.list.into.realm.lengths.rev.cumsum.CI<-binned.events.list.into.realm.lengths.rev.cumsum.CI[c((length(binned.events.list.into.realm.lengths.rev.cumsum.CI)-age):length(binned.events.list.into.realm.lengths.rev.cumsum.CI))]
  binned.events.list.into.hot.lengths.rev.cumsum.CI<-binned.events.list.into.hot.lengths.rev.cumsum.CI[c((length(binned.events.list.into.hot.lengths.rev.cumsum.CI)-age):length(binned.events.list.into.hot.lengths.rev.cumsum.CI))]
  binned.events.list.hot.into.realm.lengths.rev.cumsum.CI<-binned.events.list.hot.into.realm.lengths.rev.cumsum.CI[c((length(binned.events.list.hot.into.realm.lengths.rev.cumsum.CI)-age):length(binned.events.list.hot.into.realm.lengths.rev.cumsum.CI))]
  binned.events.list.realm.into.hot.lengths.rev.cumsum.CI<-binned.events.list.realm.into.hot.lengths.rev.cumsum.CI[c((length(binned.events.list.realm.into.hot.lengths.rev.cumsum.CI)-age):length(binned.events.list.realm.into.hot.lengths.rev.cumsum.CI))]
  binned.events.list.insitu.realm.lengths.rev.cumsum.CI<-binned.events.list.insitu.realm.lengths.rev.cumsum.CI[c((length(binned.events.list.insitu.realm.lengths.rev.cumsum.CI)-age):length(binned.events.list.insitu.realm.lengths.rev.cumsum.CI))]
  binned.events.list.insitu.hot.lengths.rev.cumsum.CI<-binned.events.list.insitu.hot.lengths.rev.cumsum.CI[c((length(binned.events.list.insitu.hot.lengths.rev.cumsum.CI)-age):length(binned.events.list.insitu.hot.lengths.rev.cumsum.CI))]
  binned.events.list.hot.local.extinction.lengths.rev.cumsum.CI<-binned.events.list.hot.local.extinction.lengths.rev.cumsum.CI[c((length(binned.events.list.hot.local.extinction.lengths.rev.cumsum.CI)-age):length(binned.events.list.hot.local.extinction.lengths.rev.cumsum.CI))]
  binned.events.list.realm.local.extinction.lengths.rev.cumsum.CI<-binned.events.list.realm.local.extinction.lengths.rev.cumsum.CI[c((length(binned.events.list.realm.local.extinction.lengths.rev.cumsum.CI)-age):length(binned.events.list.realm.local.extinction.lengths.rev.cumsum.CI))]
  
  
  
  binned.events.list.from.realm.lengths.rev.cumsum.CI.dw<-unlist(sapply(binned.events.list.from.realm.lengths.rev.cumsum.CI,function(x) x[1]))
  binned.events.list.from.realm.lengths.rev.cumsum.CI.up<-unlist(sapply(binned.events.list.from.realm.lengths.rev.cumsum.CI,function(x) x[2]))
  binned.events.list.from.hot.lengths.rev.cumsum.CI.dw<-sapply(binned.events.list.from.hot.lengths.rev.cumsum.CI,function(x) x[1])
  binned.events.list.from.hot.lengths.rev.cumsum.CI.up<-sapply(binned.events.list.from.hot.lengths.rev.cumsum.CI,function(x) x[2])
  binned.events.list.into.realm.lengths.rev.cumsum.CI.dw<-unlist(sapply(binned.events.list.into.realm.lengths.rev.cumsum.CI,function(x) x[1]))
  binned.events.list.into.realm.lengths.rev.cumsum.CI.up<-unlist(sapply(binned.events.list.into.realm.lengths.rev.cumsum.CI,function(x) x[2]))
  binned.events.list.into.hot.lengths.rev.cumsum.CI.dw<-sapply(binned.events.list.into.hot.lengths.rev.cumsum.CI,function(x) x[1])
  binned.events.list.into.hot.lengths.rev.cumsum.CI.up<-sapply(binned.events.list.into.hot.lengths.rev.cumsum.CI,function(x) x[2])
  binned.events.list.hot.into.realm.lengths.rev.cumsum.CI.dw<-unlist(sapply(binned.events.list.hot.into.realm.lengths.rev.cumsum.CI,function(x) x[1]))
  binned.events.list.hot.into.realm.lengths.rev.cumsum.CI.up<-unlist(sapply(binned.events.list.hot.into.realm.lengths.rev.cumsum.CI,function(x) x[2]))
  binned.events.list.realm.into.hot.lengths.rev.cumsum.CI.dw<-sapply(binned.events.list.realm.into.hot.lengths.rev.cumsum.CI,function(x) x[1])
  binned.events.list.realm.into.hot.lengths.rev.cumsum.CI.up<-sapply(binned.events.list.realm.into.hot.lengths.rev.cumsum.CI,function(x) x[2])
  binned.events.list.insitu.realm.lengths.rev.cumsum.CI.dw<-sapply(binned.events.list.insitu.realm.lengths.rev.cumsum.CI,function(x) x[1])
  binned.events.list.insitu.realm.lengths.rev.cumsum.CI.up<-sapply(binned.events.list.insitu.realm.lengths.rev.cumsum.CI,function(x) x[2])
  binned.events.list.insitu.hot.lengths.rev.cumsum.CI.dw<-sapply(binned.events.list.insitu.hot.lengths.rev.cumsum.CI,function(x) x[1])
  binned.events.list.insitu.hot.lengths.rev.cumsum.CI.up<-sapply(binned.events.list.insitu.hot.lengths.rev.cumsum.CI,function(x) x[2])
  binned.events.list.hot.local.extinction.lengths.rev.cumsum.CI.dw<-sapply(binned.events.list.hot.local.extinction.lengths.rev.cumsum.CI,function(x) x[1])
  binned.events.list.hot.local.extinction.lengths.rev.cumsum.CI.up<-sapply(binned.events.list.hot.local.extinction.lengths.rev.cumsum.CI,function(x) x[2])
  binned.events.list.realm.local.extinction.lengths.rev.cumsum.CI.dw<-sapply(binned.events.list.realm.local.extinction.lengths.rev.cumsum.CI,function(x) x[1])
  binned.events.list.realm.local.extinction.lengths.rev.cumsum.CI.up<-sapply(binned.events.list.realm.local.extinction.lengths.rev.cumsum.CI,function(x) x[2])
  
  
  binned.events.list.from.realm.lengths.rev.cumsum.median<-unlist(sapply(binned.events.list.from.realm.lengths.rev.cumsum.CI,function(x) x[3]))
  binned.events.list.from.hot.lengths.rev.cumsum.median<-sapply(binned.events.list.from.hot.lengths.rev.cumsum.CI,function(x) x[3])
  binned.events.list.into.realm.lengths.rev.cumsum.median<-unlist(sapply(binned.events.list.into.realm.lengths.rev.cumsum.CI,function(x) x[3]))
  binned.events.list.into.hot.lengths.rev.cumsum.median<-sapply(binned.events.list.into.hot.lengths.rev.cumsum.CI,function(x) x[3])
  binned.events.list.hot.into.realm.lengths.rev.cumsum.median<-unlist(sapply(binned.events.list.hot.into.realm.lengths.rev.cumsum.CI,function(x) x[3]))
  binned.events.list.realm.into.hot.lengths.rev.cumsum.median<-sapply(binned.events.list.realm.into.hot.lengths.rev.cumsum.CI,function(x) x[3])
  binned.events.list.insitu.realm.lengths.rev.cumsum.median<-sapply(binned.events.list.insitu.realm.lengths.rev.cumsum.CI,function(x) x[3])
  binned.events.list.insitu.hot.lengths.rev.cumsum.median<-sapply(binned.events.list.insitu.hot.lengths.rev.cumsum.CI,function(x) x[3])
  binned.events.list.hot.local.extinction.lengths.rev.cumsum.median<-sapply(binned.events.list.hot.local.extinction.lengths.rev.cumsum.CI,function(x) x[3])
  binned.events.list.realm.local.extinction.lengths.rev.cumsum.median<-sapply(binned.events.list.realm.local.extinction.lengths.rev.cumsum.CI,function(x) x[3])
  
  
  binned.events.list.from.realm.lengths.rev.cumsum.CI.dw[binned.events.list.from.realm.lengths.rev.cumsum.CI.dw==0]<-0.1
  binned.events.list.from.realm.lengths.rev.cumsum.CI.up[binned.events.list.from.realm.lengths.rev.cumsum.CI.up==0]<-0.1
  binned.events.list.from.hot.lengths.rev.cumsum.CI.dw[binned.events.list.from.hot.lengths.rev.cumsum.CI.dw==0]<-0.1
  binned.events.list.from.hot.lengths.rev.cumsum.CI.up[binned.events.list.from.hot.lengths.rev.cumsum.CI.up==0]<-0.1
  binned.events.list.into.realm.lengths.rev.cumsum.CI.dw[binned.events.list.into.realm.lengths.rev.cumsum.CI.dw==0]<-0.1
  binned.events.list.into.realm.lengths.rev.cumsum.CI.up[binned.events.list.into.realm.lengths.rev.cumsum.CI.up==0]<-0.1
  binned.events.list.into.hot.lengths.rev.cumsum.CI.dw[binned.events.list.into.hot.lengths.rev.cumsum.CI.dw==0]<-0.1
  binned.events.list.into.hot.lengths.rev.cumsum.CI.up[binned.events.list.into.hot.lengths.rev.cumsum.CI.up==0]<-0.1
  binned.events.list.hot.into.realm.lengths.rev.cumsum.CI.dw[binned.events.list.hot.into.realm.lengths.rev.cumsum.CI.dw==0]<-0.1
  binned.events.list.hot.into.realm.lengths.rev.cumsum.CI.up[binned.events.list.hot.into.realm.lengths.rev.cumsum.CI.up==0]<-0.1
  binned.events.list.realm.into.hot.lengths.rev.cumsum.CI.dw[binned.events.list.realm.into.hot.lengths.rev.cumsum.CI.dw==0]<-0.1
  binned.events.list.realm.into.hot.lengths.rev.cumsum.CI.up[binned.events.list.realm.into.hot.lengths.rev.cumsum.CI.up==0]<-0.1
  binned.events.list.insitu.realm.lengths.rev.cumsum.CI.dw[binned.events.list.insitu.realm.lengths.rev.cumsum.CI.dw==0]<-0.1
  binned.events.list.insitu.realm.lengths.rev.cumsum.CI.up[binned.events.list.insitu.realm.lengths.rev.cumsum.CI.up==0]<-0.1
  binned.events.list.insitu.hot.lengths.rev.cumsum.CI.dw[binned.events.list.insitu.hot.lengths.rev.cumsum.CI.dw==0]<-0.1
  binned.events.list.insitu.hot.lengths.rev.cumsum.CI.up[binned.events.list.insitu.hot.lengths.rev.cumsum.CI.up==0]<-0.1
  binned.events.list.hot.local.extinction.lengths.rev.cumsum.CI.dw[binned.events.list.hot.local.extinction.lengths.rev.cumsum.CI.dw==0]<-0.1
  binned.events.list.hot.local.extinction.lengths.rev.cumsum.CI.up[binned.events.list.hot.local.extinction.lengths.rev.cumsum.CI.up==0]<-0.1
  binned.events.list.realm.local.extinction.lengths.rev.cumsum.CI.dw[binned.events.list.realm.local.extinction.lengths.rev.cumsum.CI.dw==0]<-0.1
  binned.events.list.realm.local.extinction.lengths.rev.cumsum.CI.up[binned.events.list.realm.local.extinction.lengths.rev.cumsum.CI.up==0]<-0.1
  
 ##plot without logging
 #plot(c(4,4),xlim=c(0,age),ylim=c(0,max(c(unlist(binned.events.list.into.realm.lengths.rev.cumsum),unlist(binned.events.list.into.hot.lengths.rev.cumsum),unlist(binned.events.list.insitu.realm.lengths.rev.cumsum),unlist(binned.events.list.insitu.hot.lengths.rev.cumsum)))),type='n',xaxt='n',xlab='age',ylab='cumulative number of events',main=paste(name,' all events',sep=''))
 ##plot(c(4,4),xlim=c(0,217),ylim=c(0,log10(max(c(unlist(binned.events.list.into.realm.lengths.rev.cumsum),unlist(binned.events.list.into.hot.lengths.rev.cumsum),unlist(binned.events.list.insitu.realm.lengths.rev.cumsum),unlist(binned.events.list.insitu.hot.lengths.rev.cumsum))))),type='n',xaxt='n',xlab='age',ylab='cumulative number of events',main=paste(name,' all events',sep=''))
 #polygon(x=c(c(0:age),c(age:0)),y=c(binned.events.list.into.realm.lengths.rev.cumsum.CI.dw,rev(binned.events.list.into.realm.lengths.rev.cumsum.CI.up)),col=adjustcolor( "blue", alpha.f = 0.2))
 #lines(x=c(0:age),y=binned.events.list.into.realm.lengths.rev.cumsum.median,col=adjustcolor( "blue", alpha.f = 0.9))
 #text(x=age,y=binned.events.list.into.realm.lengths.rev.cumsum.median[length(binned.events.list.into.realm.lengths.rev.cumsum.median)],binned.events.list.into.realm.lengths.rev.cumsum.median[length(binned.events.list.into.realm.lengths.rev.cumsum.median)],cex=.5,col='blue')
 #polygon(x=c(c(0:age),c(age:0)),y=c(binned.events.list.into.hot.lengths.rev.cumsum.CI.dw,rev(binned.events.list.into.hot.lengths.rev.cumsum.CI.up)),col=adjustcolor( "red", alpha.f = 0.2))
 #lines(x=c(0:age),y=binned.events.list.into.hot.lengths.rev.cumsum.median,col=adjustcolor( "red", alpha.f = 0.9))
 #text(x=age,y=binned.events.list.into.hot.lengths.rev.cumsum.median[length(binned.events.list.into.hot.lengths.rev.cumsum.median)],binned.events.list.into.hot.lengths.rev.cumsum.median[length(binned.events.list.into.hot.lengths.rev.cumsum.median)],cex=.5,col='red')
 #polygon(x=c(c(0:age),c(age:0)),y=c(binned.events.list.hot.into.realm.lengths.rev.cumsum.CI.dw,rev(binned.events.list.hot.into.realm.lengths.rev.cumsum.CI.up)),col=adjustcolor( "purple", alpha.f = 0.2))
 #lines(x=c(0:age),y=binned.events.list.hot.into.realm.lengths.rev.cumsum.median,col=adjustcolor( "purple", alpha.f = 0.9))
 #text(x=age,y=binned.events.list.hot.into.realm.lengths.rev.cumsum.median[length(binned.events.list.hot.into.realm.lengths.rev.cumsum.median)],binned.events.list.hot.into.realm.lengths.rev.cumsum.median[length(binned.events.list.hot.into.realm.lengths.rev.cumsum.median)],cex=.5,col='purple')
 #polygon(x=c(c(0:age),c(age:0)),y=c(binned.events.list.realm.into.hot.lengths.rev.cumsum.CI.dw,rev(binned.events.list.realm.into.hot.lengths.rev.cumsum.CI.up)),col=adjustcolor( "yellow", alpha.f = 0.2))
 #lines(x=c(0:age),y=binned.events.list.realm.into.hot.lengths.rev.cumsum.median,col=adjustcolor( "yellow", alpha.f = 0.9))
 #text(x=age,y=binned.events.list.realm.into.hot.lengths.rev.cumsum.median[length(binned.events.list.realm.into.hot.lengths.rev.cumsum.median)],binned.events.list.realm.into.hot.lengths.rev.cumsum.median[length(binned.events.list.realm.into.hot.lengths.rev.cumsum.median)],cex=.5,col='yellow')
 #polygon(x=c(c(0:age),c(age:0)),y=c(binned.events.list.insitu.realm.lengths.rev.cumsum.CI.dw,rev(binned.events.list.insitu.realm.lengths.rev.cumsum.CI.up)),col=adjustcolor( "green", alpha.f = 0.2))
 #lines(x=c(0:age),y=binned.events.list.insitu.realm.lengths.rev.cumsum.median,col=adjustcolor( "green", alpha.f = 0.9))
 #text(x=age,y=binned.events.list.insitu.realm.lengths.rev.cumsum.median[length(binned.events.list.insitu.realm.lengths.rev.cumsum.median)],binned.events.list.insitu.realm.lengths.rev.cumsum.median[length(binned.events.list.insitu.realm.lengths.rev.cumsum.median)],cex=.5,col='green')
 #polygon(x=c(c(0:age),c(age:0)),y=c(binned.events.list.insitu.hot.lengths.rev.cumsum.CI.dw,rev(binned.events.list.insitu.hot.lengths.rev.cumsum.CI.up)),col=adjustcolor( "orange", alpha.f = 0.2))
 #lines(x=c(0:age),y=binned.events.list.insitu.hot.lengths.rev.cumsum.median,col=adjustcolor( "orange", alpha.f = 0.9))
 #text(x=age,y=binned.events.list.insitu.hot.lengths.rev.cumsum.median[length(binned.events.list.insitu.hot.lengths.rev.cumsum.median)],binned.events.list.insitu.hot.lengths.rev.cumsum.median[length(binned.events.list.insitu.hot.lengths.rev.cumsum.median)],cex=.5,col='orange')
 #legend('topleft',c('into.hot','into.non.hot.realm','insitu.hot','insitu.non.hot.realm','hot.into.realm','realm.into.hot'),col=c('red','blue','orange','green','purple','yellow'),lty=1,cex=.7,bty='n')
 #axis(1,at=c(0:age),labels=c(-age:0))
 #
 ##plot  logging
 #plot(c(4,4),xlim=c(0,age),ylim=c(0,log10(max(c(unlist(binned.events.list.into.realm.lengths.rev.cumsum),unlist(binned.events.list.into.hot.lengths.rev.cumsum),unlist(binned.events.list.insitu.realm.lengths.rev.cumsum),unlist(binned.events.list.insitu.hot.lengths.rev.cumsum))))),type='n',xaxt='n',yaxt='n',xlab='age',ylab='cumulative number of events',main=paste(name,' in situ events',sep=''))
 ##plot(c(4,4),xlim=c(0,217),ylim=c(0,log10(max(c(unlist(binned.events.list.into.realm.lengths.rev.cumsum),unlist(binned.events.list.into.hot.lengths.rev.cumsum),unlist(binned.events.list.insitu.realm.lengths.rev.cumsum),unlist(binned.events.list.insitu.hot.lengths.rev.cumsum))))),type='n',xaxt='n',xlab='age',ylab='cumulative number of events',main=paste(name,' all events',sep=''))
 #polygon(x=c(c(0:age),c(age:0)),y=log10(c(binned.events.list.into.realm.lengths.rev.cumsum.CI.dw,rev(binned.events.list.into.realm.lengths.rev.cumsum.CI.up))),col=adjustcolor( "blue", alpha.f = 0.2))
 #lines(x=c(0:age),y=log10(binned.events.list.into.realm.lengths.rev.cumsum.median),col=adjustcolor( "blue", alpha.f = 0.9))
 #text(x=age,y=log10(binned.events.list.into.realm.lengths.rev.cumsum.median[length(binned.events.list.into.realm.lengths.rev.cumsum.median)]),binned.events.list.into.realm.lengths.rev.cumsum.median[length(binned.events.list.into.realm.lengths.rev.cumsum.median)],cex=.5,col='blue')
 #polygon(x=c(c(0:age),c(age:0)),y=log10(c(binned.events.list.into.hot.lengths.rev.cumsum.CI.dw,rev(binned.events.list.into.hot.lengths.rev.cumsum.CI.up))),col=adjustcolor( "red", alpha.f = 0.2))
 #lines(x=c(0:age),y=log10(binned.events.list.into.hot.lengths.rev.cumsum.median),col=adjustcolor( "red", alpha.f = 0.9))
 #text(x=age,y=log10(binned.events.list.into.hot.lengths.rev.cumsum.median[length(binned.events.list.into.hot.lengths.rev.cumsum.median)]),binned.events.list.into.hot.lengths.rev.cumsum.median[length(binned.events.list.into.hot.lengths.rev.cumsum.median)],cex=.5,col='red')
 #polygon(x=c(c(0:age),c(age:0)),y=log10(c(binned.events.list.hot.into.realm.lengths.rev.cumsum.CI.dw,rev(binned.events.list.hot.into.realm.lengths.rev.cumsum.CI.up))),col=adjustcolor( "purple", alpha.f = 0.2))
 #lines(x=c(0:age),y=log10(binned.events.list.hot.into.realm.lengths.rev.cumsum.median),col=adjustcolor( "purple", alpha.f = 0.9))
 #text(x=age,y=log10(binned.events.list.hot.into.realm.lengths.rev.cumsum.median[length(binned.events.list.hot.into.realm.lengths.rev.cumsum.median)]),binned.events.list.hot.into.realm.lengths.rev.cumsum.median[length(binned.events.list.hot.into.realm.lengths.rev.cumsum.median)],cex=.5,col='purple')
 #polygon(x=c(c(0:age),c(age:0)),y=log10(c(binned.events.list.realm.into.hot.lengths.rev.cumsum.CI.dw,rev(binned.events.list.realm.into.hot.lengths.rev.cumsum.CI.up))),col=adjustcolor( "yellow", alpha.f = 0.2))
 #lines(x=c(0:age),y=log10(binned.events.list.realm.into.hot.lengths.rev.cumsum.median),col=adjustcolor( "yellow", alpha.f = 0.9))
 #text(x=age,y=log10(binned.events.list.realm.into.hot.lengths.rev.cumsum.median[length(binned.events.list.realm.into.hot.lengths.rev.cumsum.median)]),binned.events.list.realm.into.hot.lengths.rev.cumsum.median[length(binned.events.list.realm.into.hot.lengths.rev.cumsum.median)],cex=.5,col='yellow')
 #polygon(x=c(c(0:age),c(age:0)),y=log10(c(binned.events.list.insitu.realm.lengths.rev.cumsum.CI.dw,rev(binned.events.list.insitu.realm.lengths.rev.cumsum.CI.up))),col=adjustcolor( "green", alpha.f = 0.2))
 #lines(x=c(0:age),y=log10(binned.events.list.insitu.realm.lengths.rev.cumsum.median),col=adjustcolor( "green", alpha.f = 0.9))
 #text(x=age,y=log10(binned.events.list.insitu.realm.lengths.rev.cumsum.median[length(binned.events.list.insitu.realm.lengths.rev.cumsum.median)]),binned.events.list.insitu.realm.lengths.rev.cumsum.median[length(binned.events.list.insitu.realm.lengths.rev.cumsum.median)],cex=.5,col='green')
 #polygon(x=c(c(0:age),c(age:0)),y=log10(c(binned.events.list.insitu.hot.lengths.rev.cumsum.CI.dw,rev(binned.events.list.insitu.hot.lengths.rev.cumsum.CI.up))),col=adjustcolor( "orange", alpha.f = 0.2))
 #lines(x=c(0:age),y=log10(binned.events.list.insitu.hot.lengths.rev.cumsum.median),col=adjustcolor( "orange", alpha.f = 0.9))
 #text(x=age,y=log10(binned.events.list.insitu.hot.lengths.rev.cumsum.median[length(binned.events.list.insitu.hot.lengths.rev.cumsum.median)]),binned.events.list.insitu.hot.lengths.rev.cumsum.median[length(binned.events.list.insitu.hot.lengths.rev.cumsum.median)],cex=.5,col='orange')
 #legend('topleft',c('into.hot','into.non.hot.realm','insitu.hot','insitu.non.hot.realm','hot.into.realm','realm.into.hot'),col=c('red','blue','orange','green','purple','yellow'),lty=1,cex=.7,bty='n')
 #axis(1,at=c(0:age),labels=c(-age:0))
  #axis(2,at=log10(c(1,5,25,100,500)),labels=c('1','5','25','100','500'),las=2)
  
  #plot without logging
  #plot into and from events
  plot(c(4,4),xlim=c(0,age),ylim=c(0,max(c(unlist(binned.events.list.into.realm.lengths.rev.cumsum),unlist(binned.events.list.into.hot.lengths.rev.cumsum),unlist(binned.events.list.from.realm.lengths.rev.cumsum),unlist(binned.events.list.from.hot.lengths.rev.cumsum)))),type='n',xaxt='n',xlab='age',ylab='cumulative number of events',main=paste(name,' all events',sep=''))
  #plot(c(4,4),xlim=c(0,217),ylim=c(0,log10(max(c(unlist(binned.events.list.into.realm.lengths.rev.cumsum),unlist(binned.events.list.into.hot.lengths.rev.cumsum),unlist(binned.events.list.insitu.realm.lengths.rev.cumsum),unlist(binned.events.list.insitu.hot.lengths.rev.cumsum))))),type='n',xaxt='n',xlab='age',ylab='cumulative number of events',main=paste(name,' all events',sep=''))
  polygon(x=c(c(0:age),c(age:0)),y=c(binned.events.list.into.realm.lengths.rev.cumsum.CI.dw,rev(binned.events.list.into.realm.lengths.rev.cumsum.CI.up)),col=adjustcolor( "green", alpha.f = 0.2))
  lines(x=c(0:age),y=binned.events.list.into.realm.lengths.rev.cumsum.median,col=adjustcolor( "green", alpha.f = 0.9))
  text(x=age,y=binned.events.list.into.realm.lengths.rev.cumsum.median[length(binned.events.list.into.realm.lengths.rev.cumsum.median)],binned.events.list.into.realm.lengths.rev.cumsum.median[length(binned.events.list.into.realm.lengths.rev.cumsum.median)],cex=.5,col='green')
  polygon(x=c(c(0:age),c(age:0)),y=c(binned.events.list.into.hot.lengths.rev.cumsum.CI.dw,rev(binned.events.list.into.hot.lengths.rev.cumsum.CI.up)),col=adjustcolor( "orange", alpha.f = 0.2))
  lines(x=c(0:age),y=binned.events.list.into.hot.lengths.rev.cumsum.median,col=adjustcolor( "orange", alpha.f = 0.9))
  text(x=age,y=binned.events.list.into.hot.lengths.rev.cumsum.median[length(binned.events.list.into.hot.lengths.rev.cumsum.median)],binned.events.list.into.hot.lengths.rev.cumsum.median[length(binned.events.list.into.hot.lengths.rev.cumsum.median)],cex=.5,col='orange')
  polygon(x=c(c(0:age),c(age:0)),y=c(binned.events.list.from.realm.lengths.rev.cumsum.CI.dw,rev(binned.events.list.from.realm.lengths.rev.cumsum.CI.up)),col=adjustcolor( "blue", alpha.f = 0.2))
  lines(x=c(0:age),y=binned.events.list.from.realm.lengths.rev.cumsum.median,col=adjustcolor( "blue", alpha.f = 0.9))
  text(x=age,y=binned.events.list.from.realm.lengths.rev.cumsum.median[length(binned.events.list.from.realm.lengths.rev.cumsum.median)],binned.events.list.from.realm.lengths.rev.cumsum.median[length(binned.events.list.from.realm.lengths.rev.cumsum.median)],cex=.5,col='green')
  polygon(x=c(c(0:age),c(age:0)),y=c(binned.events.list.from.hot.lengths.rev.cumsum.CI.dw,rev(binned.events.list.from.hot.lengths.rev.cumsum.CI.up)),col=adjustcolor( "red", alpha.f = 0.2))
  lines(x=c(0:age),y=binned.events.list.from.hot.lengths.rev.cumsum.median,col=adjustcolor( "red", alpha.f = 0.9))
  text(x=age,y=binned.events.list.from.hot.lengths.rev.cumsum.median[length(binned.events.list.from.hot.lengths.rev.cumsum.median)],binned.events.list.from.hot.lengths.rev.cumsum.median[length(binned.events.list.from.hot.lengths.rev.cumsum.median)],cex=.5,col='orange')
  legend('topleft',c('into.hot','into.non.hot.realm','from.hot','from.realm'),col=c('orange','green','red','blue'),lty=1,cex=.7,bty='n')
  axis(1,at=c(0:age),labels=c(-age:0))
  
  #plot into and from events hot-into
  plot(c(4,4),xlim=c(0,age),ylim=c(0,max(c(unlist(binned.events.list.hot.into.realm.lengths.rev.cumsum),unlist(binned.events.list.realm.into.hot.lengths.rev.cumsum)))),type='n',xaxt='n',xlab='age',ylab='cumulative number of events',main=paste(name,' all events',sep=''))
  #plot(c(4,4),xlim=c(0,217),ylim=c(0,log10(max(c(unlist(binned.events.list.hot.into.realm.lengths.rev.cumsum),unlist(binned.events.list.realm.into.hot.lengths.rev.cumsum),unlist(binned.events.list.insitu.realm.lengths.rev.cumsum),unlist(binned.events.list.insitu.hot.lengths.rev.cumsum))))),type='n',xaxt='n',xlab='age',ylab='cumulative number of events',main=paste(name,' all events',sep=''))
  polygon(x=c(c(0:age),c(age:0)),y=c(binned.events.list.hot.into.realm.lengths.rev.cumsum.CI.dw,rev(binned.events.list.hot.into.realm.lengths.rev.cumsum.CI.up)),col=adjustcolor( "red", alpha.f = 0.2))
  lines(x=c(0:age),y=binned.events.list.hot.into.realm.lengths.rev.cumsum.median,col=adjustcolor( "red", alpha.f = 0.9))
  text(x=age,y=binned.events.list.hot.into.realm.lengths.rev.cumsum.median[length(binned.events.list.hot.into.realm.lengths.rev.cumsum.median)],binned.events.list.hot.into.realm.lengths.rev.cumsum.median[length(binned.events.list.hot.into.realm.lengths.rev.cumsum.median)],cex=.5,col='red')
  polygon(x=c(c(0:age),c(age:0)),y=c(binned.events.list.realm.into.hot.lengths.rev.cumsum.CI.dw,rev(binned.events.list.realm.into.hot.lengths.rev.cumsum.CI.up)),col=adjustcolor( "blue", alpha.f = 0.2))
  lines(x=c(0:age),y=binned.events.list.realm.into.hot.lengths.rev.cumsum.median,col=adjustcolor( "blue", alpha.f = 0.9))
  text(x=age,y=binned.events.list.realm.into.hot.lengths.rev.cumsum.median[length(binned.events.list.realm.into.hot.lengths.rev.cumsum.median)],binned.events.list.realm.into.hot.lengths.rev.cumsum.median[length(binned.events.list.realm.into.hot.lengths.rev.cumsum.median)],cex=.5,col='blue')
  legend('topleft',c('realm.into.hot','hot.into.realm'),col=c('blue','red'),lty=1,cex=.7,bty='n')
  axis(1,at=c(0:age),labels=c(-age:0))
  
  #plot in situ events
  plot(c(4,4),xlim=c(0,age),ylim=c(0,max(c(unlist(binned.events.list.insitu.hot.lengths.rev.cumsum),unlist(binned.events.list.insitu.realm.lengths.rev.cumsum)))),type='n',xaxt='n',xlab='age',ylab='cumulative number of events',main=paste(name,' all events',sep=''))
  #plot(c(4,4),xlim=c(0,217),ylim=c(0,log10(max(c(unlist(binned.events.list.insitu.hot.lengths.rev.cumsum),unlist(binned.events.list.insitu.realm.lengths.rev.cumsum),unlist(binned.events.list.insitu.realm.lengths.rev.cumsum),unlist(binned.events.list.insitu.hot.lengths.rev.cumsum))))),type='n',xaxt='n',xlab='age',ylab='cumulative number of events',main=paste(name,' all events',sep=''))
  polygon(x=c(c(0:age),c(age:0)),y=c(binned.events.list.insitu.hot.lengths.rev.cumsum.CI.dw,rev(binned.events.list.insitu.hot.lengths.rev.cumsum.CI.up)),col=adjustcolor( "red", alpha.f = 0.2))
  lines(x=c(0:age),y=binned.events.list.insitu.hot.lengths.rev.cumsum.median,col=adjustcolor( "red", alpha.f = 0.9))
  text(x=age,y=binned.events.list.insitu.hot.lengths.rev.cumsum.median[length(binned.events.list.insitu.hot.lengths.rev.cumsum.median)],binned.events.list.insitu.hot.lengths.rev.cumsum.median[length(binned.events.list.insitu.hot.lengths.rev.cumsum.median)],cex=.5,col='red')
  polygon(x=c(c(0:age),c(age:0)),y=c(binned.events.list.insitu.realm.lengths.rev.cumsum.CI.dw,rev(binned.events.list.insitu.realm.lengths.rev.cumsum.CI.up)),col=adjustcolor( "blue", alpha.f = 0.2))
  lines(x=c(0:age),y=binned.events.list.insitu.realm.lengths.rev.cumsum.median,col=adjustcolor( "blue", alpha.f = 0.9))
  text(x=age,y=binned.events.list.insitu.realm.lengths.rev.cumsum.median[length(binned.events.list.insitu.realm.lengths.rev.cumsum.median)],binned.events.list.insitu.realm.lengths.rev.cumsum.median[length(binned.events.list.insitu.realm.lengths.rev.cumsum.median)],cex=.5,col='blue')
  legend('topleft',c('insitu.realm','insitu.hot'),col=c('blue','red'),lty=1,cex=.7,bty='n')
  axis(1,at=c(0:age),labels=c(-age:0))
  
  #plot logging
  #plot into and from events
  plot(c(4,4),xlim=c(0,age),ylim=c(0,log10(max(c(unlist(binned.events.list.into.realm.lengths.rev.cumsum),unlist(binned.events.list.into.hot.lengths.rev.cumsum),unlist(binned.events.list.from.realm.lengths.rev.cumsum),unlist(binned.events.list.from.hot.lengths.rev.cumsum))))),type='n',xaxt='n',yaxt='n',xlab='age',ylab='cumulative number of events',main=paste(name,' all events',sep=''))
  #plot(c(4,4),xlim=c(0,217),ylim=c(0,log10(max(c(unlist(binned.events.list.into.realm.lengths.rev.cumsum),unlist(binned.events.list.into.hot.lengths.rev.cumsum),unlist(binned.events.list.insitu.realm.lengths.rev.cumsum),unlist(binned.events.list.insitu.hot.lengths.rev.cumsum))))),type='n',xaxt='n',xlab='age',ylab='cumulative number of events',main=paste(name,' all events',sep=''))
  polygon(x=c(c(0:age),c(age:0)),y=log10(c(binned.events.list.into.realm.lengths.rev.cumsum.CI.dw,rev(binned.events.list.into.realm.lengths.rev.cumsum.CI.up))),col=adjustcolor( "green", alpha.f = 0.2))
  lines(x=c(0:age),y=log10(binned.events.list.into.realm.lengths.rev.cumsum.median),col=adjustcolor( "green", alpha.f = 0.9))
  polygon(x=c(c(0:age),c(age:0)),y=log10(c(binned.events.list.into.hot.lengths.rev.cumsum.CI.dw,rev(binned.events.list.into.hot.lengths.rev.cumsum.CI.up))),col=adjustcolor( "orange", alpha.f = 0.2))
  lines(x=c(0:age),y=log10(binned.events.list.into.hot.lengths.rev.cumsum.median),col=adjustcolor( "orange", alpha.f = 0.9))
  polygon(x=c(c(0:age),c(age:0)),y=log10(c(binned.events.list.from.realm.lengths.rev.cumsum.CI.dw,rev(binned.events.list.from.realm.lengths.rev.cumsum.CI.up))),col=adjustcolor( "blue", alpha.f = 0.2))
  lines(x=c(0:age),y=log10(binned.events.list.from.realm.lengths.rev.cumsum.median),col=adjustcolor( "blue", alpha.f = 0.9))
  polygon(x=c(c(0:age),c(age:0)),y=log10(c(binned.events.list.from.hot.lengths.rev.cumsum.CI.dw,rev(binned.events.list.from.hot.lengths.rev.cumsum.CI.up))),col=adjustcolor( "red", alpha.f = 0.2))
  lines(x=c(0:age),y=log10(binned.events.list.from.hot.lengths.rev.cumsum.median),col=adjustcolor( "red", alpha.f = 0.9))
  legend('topleft',c('into.hot','into.non.hot.realm','from.hot','from.realm'),col=c('orange','green','red','blue'),lty=1,cex=.7,bty='n')
  axis(1,at=c(0:age),labels=c(-age:0))
  axis(2,at=log10(c(1,5,25,100,500)),labels=c('1','5','25','100','500'),las=2)
  
  #plot into and from events hot-into
  plot(c(4,4),xlim=c(0,age),ylim=c(0,log10(max(c(unlist(binned.events.list.hot.into.realm.lengths.rev.cumsum),unlist(binned.events.list.realm.into.hot.lengths.rev.cumsum))))),type='n',xaxt='n',yaxt='n',xlab='age',ylab='cumulative number of events',main=paste(name,' all events',sep=''))
  #plot(c(4,4),xlim=c(0,217),ylim=c(0,log10(max(c(unlist(binned.events.list.hot.into.realm.lengths.rev.cumsum),unlist(binned.events.list.realm.into.hot.lengths.rev.cumsum),unlist(binned.events.list.insitu.realm.lengths.rev.cumsum),unlist(binned.events.list.insitu.hot.lengths.rev.cumsum))))),type='n',xaxt='n',xlab='age',ylab='cumulative number of events',main=paste(name,' all events',sep=''))
  polygon(x=c(c(0:age),c(age:0)),y=log10(c(binned.events.list.hot.into.realm.lengths.rev.cumsum.CI.dw,rev(binned.events.list.hot.into.realm.lengths.rev.cumsum.CI.up))),col=adjustcolor( "red", alpha.f = 0.2))
  lines(x=c(0:age),y=log10(binned.events.list.hot.into.realm.lengths.rev.cumsum.median),col=adjustcolor( "red", alpha.f = 0.9))
  polygon(x=c(c(0:age),c(age:0)),y=log10(c(binned.events.list.realm.into.hot.lengths.rev.cumsum.CI.dw,rev(binned.events.list.realm.into.hot.lengths.rev.cumsum.CI.up))),col=adjustcolor( "blue", alpha.f = 0.2))
  lines(x=c(0:age),y=log10(binned.events.list.realm.into.hot.lengths.rev.cumsum.median),col=adjustcolor( "blue", alpha.f = 0.9))
  legend('topleft',c('realm.into.hot','hot.into.realm'),col=c('blue','red'),lty=1,cex=.7,bty='n')
  axis(1,at=c(0:age),labels=c(-age:0))
  axis(2,at=log10(c(1,5,25,100,500)),labels=c('1','5','25','100','500'),las=2)
  
  #plot in situ events
  plot(c(4,4),xlim=c(0,age),ylim=c(0,log10(max(c(unlist(binned.events.list.insitu.hot.lengths.rev.cumsum),unlist(binned.events.list.insitu.realm.lengths.rev.cumsum))))),type='n',xaxt='n',yaxt='n',xlab='age',ylab='cumulative number of events',main=paste(name,' all events',sep=''))
  #plot(c(4,4),xlim=c(0,217),ylim=c(0,log10(max(c(unlist(binned.events.list.insitu.hot.lengths.rev.cumsum),unlist(binned.events.list.insitu.realm.lengths.rev.cumsum),unlist(binned.events.list.insitu.realm.lengths.rev.cumsum),unlist(binned.events.list.insitu.hot.lengths.rev.cumsum))))),type='n',xaxt='n',xlab='age',ylab='cumulative number of events',main=paste(name,' all events',sep=''))
  polygon(x=c(c(0:age),c(age:0)),y=log10(c(binned.events.list.insitu.hot.lengths.rev.cumsum.CI.dw,rev(binned.events.list.insitu.hot.lengths.rev.cumsum.CI.up))),col=adjustcolor( "red", alpha.f = 0.2))
  lines(x=c(0:age),y=log10(binned.events.list.insitu.hot.lengths.rev.cumsum.median),col=adjustcolor( "red", alpha.f = 0.9))
  polygon(x=c(c(0:age),c(age:0)),y=log10(c(binned.events.list.insitu.realm.lengths.rev.cumsum.CI.dw,rev(binned.events.list.insitu.realm.lengths.rev.cumsum.CI.up))),col=adjustcolor( "blue", alpha.f = 0.2))
  lines(x=c(0:age),y=log10(binned.events.list.insitu.realm.lengths.rev.cumsum.median),col=adjustcolor( "blue", alpha.f = 0.9))
  legend('topleft',c('insitu.realm','insitu.hot'),col=c('blue','red'),lty=1,cex=.7,bty='n')
  axis(1,at=c(0:age),labels=c(-age:0))
  axis(2,at=log10(c(1,5,25,100,500)),labels=c('1','5','25','100','500'),las=2)
  
  
}

plot_rates_cumulative_events_BSM_logCI<-function(results.BSM,name,age){
  events.list<-list()
  for (a in 1:length(results.BSM$RES_clado_events_tables)){
    clado.df<-results.BSM$RES_clado_events_tables[[a]]
    clado.df<-clado.df[,c('time_bp','clado_event_txt')]
    clado.df<-clado.df[clado.df$clado_event_txt!='',]
    ana.df<-results.BSM$RES_ana_events_tables[[a]]
    ana.df<-ana.df[,c('abs_event_time','event_txt')]
    ana.df<-ana.df[ana.df$abs_event_time!='',]
    colnames(clado.df)<-c('time_bp','event_txt')
    colnames(ana.df)<-c('time_bp','event_txt')
    events.list[[a]]<-rbind(clado.df,ana.df)
    
  }
  
  event.type.table<-lapply(events.list,function(x)table(x$event_txt))
  
  
  if(name=='afrotrop'){
    realm.code<-'A'
  }else if(name=='austral'){
    realm.code<-'B'
  }else if(name=='indo'){
    realm.code<-'C'
  }else if(name=='nearctic'){
    realm.code<-'D'
  }else if(name=='neotrop'){
    realm.code<-'E'
  }else if(name=='palearctic'){
    realm.code<-'F'
  }
  hot.code<-'G'
  event.names<-unlist(lapply(event.type.table,function(x)unlist(x)))
  event.names<-names(table(names(event.names)))
  #anything that involves the realm
  realm.event.names<-event.names[grep(realm.code,event.names)]
  #anything that involves the hotspot
  hot.event.names<-event.names[grep(hot.code,event.names)]
  #get insitu cladogenesis in realm
  insitu.realm.event.names<-realm.event.names[grep(paste('->.*',realm.code,'.*,*',sep=''),realm.event.names)]
  insitu.realm.event.names<-insitu.realm.event.names[grep(paste(realm.code,'.*->',sep=''),insitu.realm.event.names)]
  insitu.realm.event.names<-insitu.realm.event.names[grep(paste('.*->.*',realm.code,'.*,.*',realm.code,sep=''),insitu.realm.event.names)]
  #get strict insitu cladogenesis in realm
  #insitu.strict.realm.event.names<-insitu.realm.event.names[grep(paste('^',realm.code,'->',realm.code,',',realm.code,'$',sep=''),insitu.realm.event.names)]
  #get insitu cladogenesis in hotspot
  insitu.hot.event.names<-hot.event.names[grep(paste('->.*',hot.code,'.*,*',sep=''),hot.event.names)]
  insitu.hot.event.names<-insitu.hot.event.names[grep(paste(hot.code,'.*->',sep=''),insitu.hot.event.names)]
  insitu.hot.event.names<-insitu.hot.event.names[grep(paste('.*->.*',hot.code,'.*,.*',hot.code,sep=''),insitu.hot.event.names)]
  #get strict insitu cladogenesis in realm
  #insitu.strict.hot.event.names<-insitu.hot.event.names[grep(paste('^',hot.code,'->',hot.code,',',hot.code,'$',sep=''),insitu.hot.event.names)]
  
  #get dispersal from hot into realm
  hot.into.realm.event.names<-realm.event.names[-grep(paste(realm.code,'.*->',sep=''),realm.event.names)]
  hot.into.realm.event.names<-hot.into.realm.event.names[grep(paste(hot.code,'.*->.*',sep=''),hot.into.realm.event.names)]
  #get dispersal from realm into hot
  realm.into.hot.event.names<-hot.event.names[-grep(paste(hot.code,'.*->',sep=''),hot.event.names)]
  realm.into.hot.event.names<-realm.into.hot.event.names[grep(paste(realm.code,'.*->.*',sep=''),realm.into.hot.event.names)]
  #get local extinctions in hot
  hot.local.extinction.event.names<-event.names[grep(paste('.*',hot.code,'.*->*',sep=''),event.names)]
  hot.local.extinction.event.names<-hot.local.extinction.event.names[-grep(paste('.*->.*',hot.code,'.*',sep=''),hot.local.extinction.event.names)]
  #get local extinctions in realm
  realm.local.extinction.event.names<-event.names[grep(paste('.*',realm.code,'.*->.*',sep=''),event.names)]
  realm.local.extinction.event.names<-realm.local.extinction.event.names[-grep(paste('.*->.*',realm.code,'.*',sep=''),realm.local.extinction.event.names)]
  #get dispersal into realm
  dispersal.into.realm.event.names<-realm.event.names[-grep(paste(realm.code,'.*->',sep=''),realm.event.names)]
  #get dispersal into hotspot
  dispersal.into.hot.event.names<-hot.event.names[-grep(paste(hot.code,'.*->',sep=''),hot.event.names)]
  #get dispersal from hotspot
  dispersal.from.hot.event.names<-event.names[grep(paste('.*',hot.code,'.*->*',sep=''),event.names)]
  dispersal.from.hot.event.names<-dispersal.from.hot.event.names[!(dispersal.from.hot.event.names%in%insitu.hot.event.names)]
  dispersal.from.hot.event.names<-dispersal.from.hot.event.names[!(dispersal.from.hot.event.names%in%hot.local.extinction.event.names)]
  dispersal.from.hot.event.names<-dispersal.from.hot.event.names[-grep(paste('->',hot.code,'$',sep=''),dispersal.from.hot.event.names)]
  #get dispersal from realm
  dispersal.from.realm.event.names<-event.names[grep(paste('.*',realm.code,'.*->*',sep=''),event.names)]
  dispersal.from.realm.event.names<-dispersal.from.realm.event.names[!(dispersal.from.realm.event.names%in%insitu.realm.event.names)]
  dispersal.from.realm.event.names<-dispersal.from.realm.event.names[!(dispersal.from.realm.event.names%in%realm.local.extinction.event.names)]
  dispersal.from.realm.event.names<-dispersal.from.realm.event.names[-grep(paste('->',realm.code,'$',sep=''),dispersal.from.realm.event.names)]
  
  dispersal.into.realm<-unlist(lapply(event.type.table,function(x)sum(x[dispersal.into.realm.event.names],na.rm=TRUE)))
  dispersal.into.hot<-unlist(lapply(event.type.table,function(x)sum(x[dispersal.into.hot.event.names],na.rm=TRUE)))
  dispersal.from.realm<-unlist(lapply(event.type.table,function(x)sum(x[dispersal.from.realm.event.names],na.rm=TRUE)))
  dispersal.from.hot<-unlist(lapply(event.type.table,function(x)sum(x[dispersal.from.hot.event.names],na.rm=TRUE)))
  dispersal.hot.into.realm<-unlist(lapply(event.type.table,function(x)sum(x[hot.into.realm.event.names],na.rm=TRUE)))
  dispersal.realm.into.hot<-unlist(lapply(event.type.table,function(x)sum(x[realm.into.hot.event.names],na.rm=TRUE)))
  #insitu.strict.realm<-unlist(lapply(event.type.table,function(x)sum(x[insitu.strict.realm.event.names],na.rm=TRUE)))
  #insitu.strict.hot<-unlist(lapply(event.type.table,function(x)sum(x[insitu.strict.hot.event.names],na.rm=TRUE)))
  insitu.realm<-unlist(lapply(event.type.table,function(x)sum(x[insitu.realm.event.names],na.rm=TRUE)))
  insitu.hot<-unlist(lapply(event.type.table,function(x)sum(x[insitu.hot.event.names],na.rm=TRUE)))
  hot.local.extinction<-unlist(lapply(event.type.table,function(x)sum(x[hot.local.extinction.event.names],na.rm=TRUE)))
  realm.local.extinction<-unlist(lapply(event.type.table,function(x)sum(x[realm.local.extinction.event.names],na.rm=TRUE)))
  
  dispersal.from.realm<-lapply(dispersal.from.realm,function(x) replace(x,is.na(x),0))
  dispersal.from.hot<-lapply(dispersal.from.hot,function(x) replace(x,is.na(x),0))
  dispersal.into.realm<-lapply(dispersal.into.realm,function(x) replace(x,is.na(x),0))
  dispersal.into.hot<-lapply(dispersal.into.hot,function(x) replace(x,is.na(x),0))
  dispersal.hot.into.realm<-lapply(dispersal.hot.into.realm,function(x) replace(x,is.na(x),0))
  dispersal.realm.into.hot<-lapply(dispersal.realm.into.hot,function(x) replace(x,is.na(x),0))
  #insitu.strict.realm<-lapply(insitu.strict.realm,function(x) replace(x,is.na(x),0))
  #insitu.strict.hot<-lapply(insitu.strict.hot,function(x) replace(x,is.na(x),0))
  insitu.realm<-lapply(insitu.realm,function(x) replace(x,is.na(x),0))
  insitu.hot<-lapply(insitu.hot,function(x) replace(x,is.na(x),0))
  hot.local.extinction<-lapply(hot.local.extinction,function(x) replace(x,is.na(x),0))
  realm.local.extinction<-lapply(realm.local.extinction,function(x) replace(x,is.na(x),0))
  
  ###########cumulative sum of events plot
  #process a results.BSM to extract events and times in a list of lists (each element of the list is one clade with a list of 50 replicate BSMs)
  events.list<-list()
  for (a in 1:length(results.BSM$RES_clado_events_tables)){
    clado.df<-results.BSM$RES_clado_events_tables[[a]]
    clado.df<-clado.df[,c('time_bp','clado_event_txt')]
    clado.df<-clado.df[clado.df$clado_event_txt!='',]
    ana.df<-results.BSM$RES_ana_events_tables[[a]]
    ana.df<-ana.df[,c('abs_event_time','event_txt')]
    ana.df<-ana.df[ana.df$abs_event_time!='',]
    colnames(clado.df)<-c('time_bp','event_txt')
    colnames(ana.df)<-c('time_bp','event_txt')
    events.list[[a]]<-rbind(clado.df,ana.df)
  }
  
  #create a df that bins events into 1 myr slots
  #get the roots of all clades
  clade.roots<-results.BSM$RES_clado_events_tables[[1]][results.BSM$RES_clado_events_tables[[1]]$node.type=='root',]$time_bp
  #another way of binning events, a list of lists (binning events for each BSM replicate)
  binned.events.list<-vector('list',length=length(events.list))
  cat('binning events','\n')
  for (i in 1:length(events.list)){
    # cat(i,'\n')
    binned.events.list[[i]]<-vector('list',length=length(seq(0,floor(clade.roots)))+1)
    for (b in 1:nrow(events.list[[i]])){
      # cat(b,'\n')
      time_bp<-floor(events.list[[i]][b,'time_bp'])
      event<-events.list[[i]][b,'event_txt']
      names(event)<-row.names(events.list[[i]][b,])
      if(is.null((binned.events.list[[i]][[time_bp+1]]))){
        binned.events.list[[i]][[time_bp+1]]<-c('0',event)
      }else{
        binned.events.list[[i]][[time_bp+1]]<-c(binned.events.list[[i]][[time_bp+1]],event)  
      }
    }
  }
  binned.events.list.from.realm<-list()
  binned.events.list.from.hot<-list()
  binned.events.list.into.realm<-list()
  binned.events.list.into.hot<-list()
  binned.events.list.hot.into.realm<-list()
  binned.events.list.realm.into.hot<-list()
  binned.events.list.insitu.realm<-list()
  binned.events.list.insitu.hot<-list()
  binned.events.list.hot.local.extinction<-list()
  binned.events.list.realm.local.extinction<-list()
  
  
  
  cat('subsetting binned events to realm','\n')
  for(i in 1:length(binned.events.list)){
    #cat(i,'\n')
    binned.events.list.from.realm[[i]]<-vector('list',length=length(seq(0,floor(clade.roots)))+1)
    binned.events.list.from.hot[[i]]<-vector('list',length=length(seq(0,floor(clade.roots)))+1)
    binned.events.list.into.realm[[i]]<-vector('list',length=length(seq(0,floor(clade.roots)))+1)
    binned.events.list.into.hot[[i]]<-vector('list',length=length(seq(0,floor(clade.roots)))+1)
    binned.events.list.hot.into.realm[[i]]<-vector('list',length=length(seq(0,floor(clade.roots)))+1)
    binned.events.list.realm.into.hot[[i]]<-vector('list',length=length(seq(0,floor(clade.roots)))+1)
    binned.events.list.insitu.realm[[i]]<-vector('list',length=length(seq(0,floor(clade.roots)))+1)
    binned.events.list.insitu.hot[[i]]<-vector('list',length=length(seq(0,floor(clade.roots)))+1)
    binned.events.list.hot.local.extinction[[i]]<-vector('list',length=length(seq(0,floor(clade.roots)))+1)
    binned.events.list.realm.local.extinction[[i]]<-vector('list',length=length(seq(0,floor(clade.roots)))+1)
    
    
    for (a in 1:length(binned.events.list[[i]])){
      # cat(a,'\n')
      binned.events.list.from.realm[[i]][[a]]<-binned.events.list[[i]][[a]][unlist(sapply(dispersal.from.realm.event.names,function(x) grep(paste('^',x,'$',sep=''),binned.events.list[[i]][[a]])))]
      binned.events.list.from.hot[[i]][[a]]<-binned.events.list[[i]][[a]][unlist(sapply(dispersal.from.hot.event.names,function(x)grep(paste('^',x,'$',sep=''),binned.events.list[[i]][[a]])))]
      binned.events.list.into.realm[[i]][[a]]<-binned.events.list[[i]][[a]][unlist(sapply(dispersal.into.realm.event.names,function(x) grep(paste('^',x,'$',sep=''),binned.events.list[[i]][[a]])))]
      binned.events.list.into.hot[[i]][[a]]<-binned.events.list[[i]][[a]][unlist(sapply(dispersal.into.hot.event.names,function(x)grep(paste('^',x,'$',sep=''),binned.events.list[[i]][[a]])))]
      binned.events.list.hot.into.realm[[i]][[a]]<-binned.events.list[[i]][[a]][unlist(sapply(hot.into.realm.event.names,function(x) grep(paste('^',x,'$',sep=''),binned.events.list[[i]][[a]])))]
      binned.events.list.realm.into.hot[[i]][[a]]<-binned.events.list[[i]][[a]][unlist(sapply(realm.into.hot.event.names,function(x)grep(paste('^',x,'$',sep=''),binned.events.list[[i]][[a]])))]
      binned.events.list.insitu.realm[[i]][[a]]<-binned.events.list[[i]][[a]][unlist(sapply(insitu.realm.event.names,function(x)grep(paste('^',x,'$',sep=''),binned.events.list[[i]][[a]])))]
      binned.events.list.insitu.hot[[i]][[a]]<-binned.events.list[[i]][[a]][unlist(sapply(insitu.hot.event.names,function(x)grep(paste('^',x,'$',sep=''),binned.events.list[[i]][[a]])))]
      binned.events.list.hot.local.extinction[[i]][[a]]<-binned.events.list[[i]][[a]][unlist(sapply(hot.local.extinction.event.names,function(x) grep(paste('^',x,'$',sep=''),binned.events.list[[i]][[a]])))]
      binned.events.list.realm.local.extinction[[i]][[a]]<-binned.events.list[[i]][[a]][unlist(sapply(realm.local.extinction.event.names,function(x)grep(paste('^',x,'$',sep=''),binned.events.list[[i]][[a]])))]
      
    }
    
    
    
    
  }
  binned.events.list.from.realm.lengths<-lapply(binned.events.list.from.realm,function(x)unlist(lapply(x,function(x) length(x))))
  binned.events.list.from.hot.lengths<-lapply(binned.events.list.from.hot,function(x)unlist(lapply(x,function(x) length(x))))
  binned.events.list.into.realm.lengths<-lapply(binned.events.list.into.realm,function(x)unlist(lapply(x,function(x) length(x))))
  binned.events.list.into.hot.lengths<-lapply(binned.events.list.into.hot,function(x)unlist(lapply(x,function(x) length(x))))
  binned.events.list.hot.into.realm.lengths<-lapply(binned.events.list.hot.into.realm,function(x)unlist(lapply(x,function(x) length(x))))
  binned.events.list.realm.into.hot.lengths<-lapply(binned.events.list.realm.into.hot,function(x)unlist(lapply(x,function(x) length(x))))
  binned.events.list.insitu.realm.lengths<-lapply(binned.events.list.insitu.realm,function(x)unlist(lapply(x,function(x) length(x))))
  binned.events.list.insitu.hot.lengths<-lapply(binned.events.list.insitu.hot,function(x)unlist(lapply(x,function(x) length(x))))
  binned.events.list.hot.local.extinction.lengths<-lapply(binned.events.list.hot.local.extinction,function(x)unlist(lapply(x,function(x) length(x))))
  binned.events.list.realm.local.extinction.lengths<-lapply(binned.events.list.realm.local.extinction,function(x)unlist(lapply(x,function(x) length(x))))
  
  
  binned.events.list.from.realm.lengths.rev<-lapply(binned.events.list.from.realm.lengths,function(x) rev(x))
  binned.events.list.from.hot.lengths.rev<-lapply(binned.events.list.from.hot.lengths,function(x) rev(x))
  binned.events.list.into.realm.lengths.rev<-lapply(binned.events.list.into.realm.lengths,function(x) rev(x))
  binned.events.list.into.hot.lengths.rev<-lapply(binned.events.list.into.hot.lengths,function(x) rev(x))
  binned.events.list.into.realm.lengths.rev<-lapply(binned.events.list.into.realm.lengths,function(x) rev(x))
  binned.events.list.into.hot.lengths.rev<-lapply(binned.events.list.into.hot.lengths,function(x) rev(x))
  binned.events.list.hot.into.realm.lengths.rev<-lapply(binned.events.list.hot.into.realm.lengths,function(x) rev(x))
  binned.events.list.realm.into.hot.lengths.rev<-lapply(binned.events.list.realm.into.hot.lengths,function(x) rev(x))
  binned.events.list.insitu.realm.lengths.rev<-lapply(binned.events.list.insitu.realm.lengths,function(x) rev(x))
  binned.events.list.insitu.hot.lengths.rev<-lapply(binned.events.list.insitu.hot.lengths,function(x) rev(x))
  binned.events.list.hot.local.extinction.lengths.rev<-lapply(binned.events.list.hot.local.extinction.lengths,function(x) rev(x))
  binned.events.list.realm.local.extinction.lengths.rev<-lapply(binned.events.list.realm.local.extinction.lengths,function(x) rev(x))
  
  
  
  binned.events.list.from.realm.lengths.rev.cumsum<-lapply(binned.events.list.from.realm.lengths.rev,function(x) cumsum(x))
  binned.events.list.from.hot.lengths.rev.cumsum<-lapply(binned.events.list.from.hot.lengths.rev,function(x) cumsum(x))
  binned.events.list.into.realm.lengths.rev.cumsum<-lapply(binned.events.list.into.realm.lengths.rev,function(x) cumsum(x))
  binned.events.list.into.hot.lengths.rev.cumsum<-lapply(binned.events.list.into.hot.lengths.rev,function(x) cumsum(x))
  binned.events.list.hot.into.realm.lengths.rev.cumsum<-lapply(binned.events.list.hot.into.realm.lengths.rev,function(x) cumsum(x))
  binned.events.list.realm.into.hot.lengths.rev.cumsum<-lapply(binned.events.list.realm.into.hot.lengths.rev,function(x) cumsum(x))
  binned.events.list.insitu.realm.lengths.rev.cumsum<-lapply(binned.events.list.insitu.realm.lengths.rev,function(x) cumsum(x))
  binned.events.list.insitu.hot.lengths.rev.cumsum<-lapply(binned.events.list.insitu.hot.lengths.rev,function(x) cumsum(x))
  binned.events.list.hot.local.extinction.lengths.rev.cumsum<-lapply(binned.events.list.hot.local.extinction.lengths.rev,function(x) cumsum(x))
  binned.events.list.realm.local.extinction.lengths.rev.cumsum<-lapply(binned.events.list.realm.local.extinction.lengths.rev,function(x) cumsum(x))
  
  
  binned.events.list.from.realm.lengths.rev.cumsum.CI<-list()
  binned.events.list.from.hot.lengths.rev.cumsum.CI<-list()
  binned.events.list.into.realm.lengths.rev.cumsum.CI<-list()
  binned.events.list.into.hot.lengths.rev.cumsum.CI<-list()
  binned.events.list.hot.into.realm.lengths.rev.cumsum.CI<-list()
  binned.events.list.realm.into.hot.lengths.rev.cumsum.CI<-list()
  binned.events.list.insitu.realm.lengths.rev.cumsum.CI<-list()
  binned.events.list.insitu.hot.lengths.rev.cumsum.CI<-list()
  binned.events.list.hot.local.extinction.lengths.rev.cumsum.CI<-list()
  binned.events.list.realm.local.extinction.lengths.rev.cumsum.CI<-list()
  
  
  for (i in 1:length(binned.events.list.into.realm.lengths.rev.cumsum[[1]])){
    binned.events.list.from.realm.lengths.rev.cumsum.CI[[i]]<-sapply(binned.events.list.from.realm.lengths.rev.cumsum,"[[",i)
    binned.events.list.from.hot.lengths.rev.cumsum.CI[[i]]<-sapply(binned.events.list.from.hot.lengths.rev.cumsum,"[[",i)
    binned.events.list.into.realm.lengths.rev.cumsum.CI[[i]]<-sapply(binned.events.list.into.realm.lengths.rev.cumsum,"[[",i)
    binned.events.list.into.hot.lengths.rev.cumsum.CI[[i]]<-sapply(binned.events.list.into.hot.lengths.rev.cumsum,"[[",i)
    binned.events.list.hot.into.realm.lengths.rev.cumsum.CI[[i]]<-sapply(binned.events.list.hot.into.realm.lengths.rev.cumsum,"[[",i)
    binned.events.list.realm.into.hot.lengths.rev.cumsum.CI[[i]]<-sapply(binned.events.list.realm.into.hot.lengths.rev.cumsum,"[[",i)
    binned.events.list.insitu.realm.lengths.rev.cumsum.CI[[i]]<-sapply(binned.events.list.insitu.realm.lengths.rev.cumsum,"[[",i)
    binned.events.list.insitu.hot.lengths.rev.cumsum.CI[[i]]<-sapply(binned.events.list.insitu.hot.lengths.rev.cumsum,"[[",i)
    binned.events.list.hot.local.extinction.lengths.rev.cumsum.CI[[i]]<-sapply(binned.events.list.hot.local.extinction.lengths.rev.cumsum,"[[",i)
    binned.events.list.realm.local.extinction.lengths.rev.cumsum.CI[[i]]<-sapply(binned.events.list.realm.local.extinction.lengths.rev.cumsum,"[[",i)
  }
  
  binned.events.list.from.realm.lengths.rev.cumsum.CI<-lapply(binned.events.list.from.realm.lengths.rev.cumsum.CI,function(x)quantile(x,c(0.05,0.95,0.50)))
  binned.events.list.from.hot.lengths.rev.cumsum.CI<-lapply(binned.events.list.from.hot.lengths.rev.cumsum.CI,function(x)quantile(x,c(0.05,0.95,0.50)))
  binned.events.list.into.realm.lengths.rev.cumsum.CI<-lapply(binned.events.list.into.realm.lengths.rev.cumsum.CI,function(x)quantile(x,c(0.05,0.95,0.50)))
  binned.events.list.into.hot.lengths.rev.cumsum.CI<-lapply(binned.events.list.into.hot.lengths.rev.cumsum.CI,function(x)quantile(x,c(0.05,0.95,0.50)))
  binned.events.list.hot.into.realm.lengths.rev.cumsum.CI<-lapply(binned.events.list.hot.into.realm.lengths.rev.cumsum.CI,function(x)quantile(x,c(0.05,0.95,0.50)))
  binned.events.list.realm.into.hot.lengths.rev.cumsum.CI<-lapply(binned.events.list.realm.into.hot.lengths.rev.cumsum.CI,function(x)quantile(x,c(0.05,0.95,0.50)))
  binned.events.list.insitu.realm.lengths.rev.cumsum.CI<-lapply(binned.events.list.insitu.realm.lengths.rev.cumsum.CI,function(x)quantile(x,c(0.05,0.95,0.50)))
  binned.events.list.insitu.hot.lengths.rev.cumsum.CI<-lapply(binned.events.list.insitu.hot.lengths.rev.cumsum.CI,function(x)quantile(x,c(0.05,0.95,0.50)))
  binned.events.list.hot.local.extinction.lengths.rev.cumsum.CI<-lapply(binned.events.list.hot.local.extinction.lengths.rev.cumsum.CI,function(x)quantile(x,c(0.05,0.95,0.50)))
  binned.events.list.realm.local.extinction.lengths.rev.cumsum.CI<-lapply(binned.events.list.realm.local.extinction.lengths.rev.cumsum.CI,function(x)quantile(x,c(0.05,0.95,0.50)))
  
  
  
  
  binned.events.list.from.realm.lengths.rev.cumsum.CI<-binned.events.list.from.realm.lengths.rev.cumsum.CI[c((length(binned.events.list.from.realm.lengths.rev.cumsum.CI)-age):length(binned.events.list.from.realm.lengths.rev.cumsum.CI))]
  binned.events.list.from.hot.lengths.rev.cumsum.CI<-binned.events.list.from.hot.lengths.rev.cumsum.CI[c((length(binned.events.list.from.hot.lengths.rev.cumsum.CI)-age):length(binned.events.list.from.hot.lengths.rev.cumsum.CI))]
  binned.events.list.into.realm.lengths.rev.cumsum.CI<-binned.events.list.into.realm.lengths.rev.cumsum.CI[c((length(binned.events.list.into.realm.lengths.rev.cumsum.CI)-age):length(binned.events.list.into.realm.lengths.rev.cumsum.CI))]
  binned.events.list.into.hot.lengths.rev.cumsum.CI<-binned.events.list.into.hot.lengths.rev.cumsum.CI[c((length(binned.events.list.into.hot.lengths.rev.cumsum.CI)-age):length(binned.events.list.into.hot.lengths.rev.cumsum.CI))]
  binned.events.list.hot.into.realm.lengths.rev.cumsum.CI<-binned.events.list.hot.into.realm.lengths.rev.cumsum.CI[c((length(binned.events.list.hot.into.realm.lengths.rev.cumsum.CI)-age):length(binned.events.list.hot.into.realm.lengths.rev.cumsum.CI))]
  binned.events.list.realm.into.hot.lengths.rev.cumsum.CI<-binned.events.list.realm.into.hot.lengths.rev.cumsum.CI[c((length(binned.events.list.realm.into.hot.lengths.rev.cumsum.CI)-age):length(binned.events.list.realm.into.hot.lengths.rev.cumsum.CI))]
  binned.events.list.insitu.realm.lengths.rev.cumsum.CI<-binned.events.list.insitu.realm.lengths.rev.cumsum.CI[c((length(binned.events.list.insitu.realm.lengths.rev.cumsum.CI)-age):length(binned.events.list.insitu.realm.lengths.rev.cumsum.CI))]
  binned.events.list.insitu.hot.lengths.rev.cumsum.CI<-binned.events.list.insitu.hot.lengths.rev.cumsum.CI[c((length(binned.events.list.insitu.hot.lengths.rev.cumsum.CI)-age):length(binned.events.list.insitu.hot.lengths.rev.cumsum.CI))]
  binned.events.list.hot.local.extinction.lengths.rev.cumsum.CI<-binned.events.list.hot.local.extinction.lengths.rev.cumsum.CI[c((length(binned.events.list.hot.local.extinction.lengths.rev.cumsum.CI)-age):length(binned.events.list.hot.local.extinction.lengths.rev.cumsum.CI))]
  binned.events.list.realm.local.extinction.lengths.rev.cumsum.CI<-binned.events.list.realm.local.extinction.lengths.rev.cumsum.CI[c((length(binned.events.list.realm.local.extinction.lengths.rev.cumsum.CI)-age):length(binned.events.list.realm.local.extinction.lengths.rev.cumsum.CI))]
  
  
  
  binned.events.list.from.realm.lengths.rev.cumsum.CI.dw<-unlist(sapply(binned.events.list.from.realm.lengths.rev.cumsum.CI,function(x) x[1]))
  binned.events.list.from.realm.lengths.rev.cumsum.CI.up<-unlist(sapply(binned.events.list.from.realm.lengths.rev.cumsum.CI,function(x) x[2]))
  binned.events.list.from.hot.lengths.rev.cumsum.CI.dw<-sapply(binned.events.list.from.hot.lengths.rev.cumsum.CI,function(x) x[1])
  binned.events.list.from.hot.lengths.rev.cumsum.CI.up<-sapply(binned.events.list.from.hot.lengths.rev.cumsum.CI,function(x) x[2])
  binned.events.list.into.realm.lengths.rev.cumsum.CI.dw<-unlist(sapply(binned.events.list.into.realm.lengths.rev.cumsum.CI,function(x) x[1]))
  binned.events.list.into.realm.lengths.rev.cumsum.CI.up<-unlist(sapply(binned.events.list.into.realm.lengths.rev.cumsum.CI,function(x) x[2]))
  binned.events.list.into.hot.lengths.rev.cumsum.CI.dw<-sapply(binned.events.list.into.hot.lengths.rev.cumsum.CI,function(x) x[1])
  binned.events.list.into.hot.lengths.rev.cumsum.CI.up<-sapply(binned.events.list.into.hot.lengths.rev.cumsum.CI,function(x) x[2])
  binned.events.list.hot.into.realm.lengths.rev.cumsum.CI.dw<-unlist(sapply(binned.events.list.hot.into.realm.lengths.rev.cumsum.CI,function(x) x[1]))
  binned.events.list.hot.into.realm.lengths.rev.cumsum.CI.up<-unlist(sapply(binned.events.list.hot.into.realm.lengths.rev.cumsum.CI,function(x) x[2]))
  binned.events.list.realm.into.hot.lengths.rev.cumsum.CI.dw<-sapply(binned.events.list.realm.into.hot.lengths.rev.cumsum.CI,function(x) x[1])
  binned.events.list.realm.into.hot.lengths.rev.cumsum.CI.up<-sapply(binned.events.list.realm.into.hot.lengths.rev.cumsum.CI,function(x) x[2])
  binned.events.list.insitu.realm.lengths.rev.cumsum.CI.dw<-sapply(binned.events.list.insitu.realm.lengths.rev.cumsum.CI,function(x) x[1])
  binned.events.list.insitu.realm.lengths.rev.cumsum.CI.up<-sapply(binned.events.list.insitu.realm.lengths.rev.cumsum.CI,function(x) x[2])
  binned.events.list.insitu.hot.lengths.rev.cumsum.CI.dw<-sapply(binned.events.list.insitu.hot.lengths.rev.cumsum.CI,function(x) x[1])
  binned.events.list.insitu.hot.lengths.rev.cumsum.CI.up<-sapply(binned.events.list.insitu.hot.lengths.rev.cumsum.CI,function(x) x[2])
  binned.events.list.hot.local.extinction.lengths.rev.cumsum.CI.dw<-sapply(binned.events.list.hot.local.extinction.lengths.rev.cumsum.CI,function(x) x[1])
  binned.events.list.hot.local.extinction.lengths.rev.cumsum.CI.up<-sapply(binned.events.list.hot.local.extinction.lengths.rev.cumsum.CI,function(x) x[2])
  binned.events.list.realm.local.extinction.lengths.rev.cumsum.CI.dw<-sapply(binned.events.list.realm.local.extinction.lengths.rev.cumsum.CI,function(x) x[1])
  binned.events.list.realm.local.extinction.lengths.rev.cumsum.CI.up<-sapply(binned.events.list.realm.local.extinction.lengths.rev.cumsum.CI,function(x) x[2])
  
  
  binned.events.list.from.realm.lengths.rev.cumsum.median<-unlist(sapply(binned.events.list.from.realm.lengths.rev.cumsum.CI,function(x) x[3]))
  binned.events.list.from.hot.lengths.rev.cumsum.median<-sapply(binned.events.list.from.hot.lengths.rev.cumsum.CI,function(x) x[3])
  binned.events.list.into.realm.lengths.rev.cumsum.median<-unlist(sapply(binned.events.list.into.realm.lengths.rev.cumsum.CI,function(x) x[3]))
  binned.events.list.into.hot.lengths.rev.cumsum.median<-sapply(binned.events.list.into.hot.lengths.rev.cumsum.CI,function(x) x[3])
  binned.events.list.hot.into.realm.lengths.rev.cumsum.median<-unlist(sapply(binned.events.list.hot.into.realm.lengths.rev.cumsum.CI,function(x) x[3]))
  binned.events.list.realm.into.hot.lengths.rev.cumsum.median<-sapply(binned.events.list.realm.into.hot.lengths.rev.cumsum.CI,function(x) x[3])
  binned.events.list.insitu.realm.lengths.rev.cumsum.median<-sapply(binned.events.list.insitu.realm.lengths.rev.cumsum.CI,function(x) x[3])
  binned.events.list.insitu.hot.lengths.rev.cumsum.median<-sapply(binned.events.list.insitu.hot.lengths.rev.cumsum.CI,function(x) x[3])
  binned.events.list.hot.local.extinction.lengths.rev.cumsum.median<-sapply(binned.events.list.hot.local.extinction.lengths.rev.cumsum.CI,function(x) x[3])
  binned.events.list.realm.local.extinction.lengths.rev.cumsum.median<-sapply(binned.events.list.realm.local.extinction.lengths.rev.cumsum.CI,function(x) x[3])
  
  
  binned.events.list.from.realm.lengths.rev.cumsum.CI.dw[binned.events.list.from.realm.lengths.rev.cumsum.CI.dw==0]<-0.1
  binned.events.list.from.realm.lengths.rev.cumsum.CI.up[binned.events.list.from.realm.lengths.rev.cumsum.CI.up==0]<-0.1
  binned.events.list.from.hot.lengths.rev.cumsum.CI.dw[binned.events.list.from.hot.lengths.rev.cumsum.CI.dw==0]<-0.1
  binned.events.list.from.hot.lengths.rev.cumsum.CI.up[binned.events.list.from.hot.lengths.rev.cumsum.CI.up==0]<-0.1
  binned.events.list.into.realm.lengths.rev.cumsum.CI.dw[binned.events.list.into.realm.lengths.rev.cumsum.CI.dw==0]<-0.1
  binned.events.list.into.realm.lengths.rev.cumsum.CI.up[binned.events.list.into.realm.lengths.rev.cumsum.CI.up==0]<-0.1
  binned.events.list.into.hot.lengths.rev.cumsum.CI.dw[binned.events.list.into.hot.lengths.rev.cumsum.CI.dw==0]<-0.1
  binned.events.list.into.hot.lengths.rev.cumsum.CI.up[binned.events.list.into.hot.lengths.rev.cumsum.CI.up==0]<-0.1
  binned.events.list.hot.into.realm.lengths.rev.cumsum.CI.dw[binned.events.list.hot.into.realm.lengths.rev.cumsum.CI.dw==0]<-0.1
  binned.events.list.hot.into.realm.lengths.rev.cumsum.CI.up[binned.events.list.hot.into.realm.lengths.rev.cumsum.CI.up==0]<-0.1
  binned.events.list.realm.into.hot.lengths.rev.cumsum.CI.dw[binned.events.list.realm.into.hot.lengths.rev.cumsum.CI.dw==0]<-0.1
  binned.events.list.realm.into.hot.lengths.rev.cumsum.CI.up[binned.events.list.realm.into.hot.lengths.rev.cumsum.CI.up==0]<-0.1
  binned.events.list.insitu.realm.lengths.rev.cumsum.CI.dw[binned.events.list.insitu.realm.lengths.rev.cumsum.CI.dw==0]<-0.1
  binned.events.list.insitu.realm.lengths.rev.cumsum.CI.up[binned.events.list.insitu.realm.lengths.rev.cumsum.CI.up==0]<-0.1
  binned.events.list.insitu.hot.lengths.rev.cumsum.CI.dw[binned.events.list.insitu.hot.lengths.rev.cumsum.CI.dw==0]<-0.1
  binned.events.list.insitu.hot.lengths.rev.cumsum.CI.up[binned.events.list.insitu.hot.lengths.rev.cumsum.CI.up==0]<-0.1
  binned.events.list.hot.local.extinction.lengths.rev.cumsum.CI.dw[binned.events.list.hot.local.extinction.lengths.rev.cumsum.CI.dw==0]<-0.1
  binned.events.list.hot.local.extinction.lengths.rev.cumsum.CI.up[binned.events.list.hot.local.extinction.lengths.rev.cumsum.CI.up==0]<-0.1
  binned.events.list.realm.local.extinction.lengths.rev.cumsum.CI.dw[binned.events.list.realm.local.extinction.lengths.rev.cumsum.CI.dw==0]<-0.1
  binned.events.list.realm.local.extinction.lengths.rev.cumsum.CI.up[binned.events.list.realm.local.extinction.lengths.rev.cumsum.CI.up==0]<-0.1
  
  #in situ rates in a 1-my period as number of events in a period / number of lineages (=cumulative number of in situ speciation + colonization into - local extinction)
  binned.events.list.from.realm.lengths.rev.median<-list()
  binned.events.list.from.hot.lengths.rev.median<-list()
  binned.events.list.into.realm.lengths.rev.median<-list()
  binned.events.list.into.hot.lengths.rev.median<-list()  
  binned.events.list.hot.into.realm.lengths.rev.median<-list()
  binned.events.list.realm.into.hot.lengths.rev.median<-list()
  binned.events.list.insitu.realm.lengths.rev.median<-list()
  binned.events.list.insitu.hot.lengths.rev.median<-list()
  
  for (i in 1:length(binned.events.list.into.realm.lengths.rev.cumsum[[1]])){
    binned.events.list.from.realm.lengths.rev.median[[i]]<-sapply(binned.events.list.from.realm.lengths.rev.cumsum,"[[",i)
    binned.events.list.from.hot.lengths.rev.median[[i]]<-sapply(binned.events.list.from.hot.lengths.rev.cumsum,"[[",i)
    binned.events.list.into.realm.lengths.rev.median[[i]]<-sapply(binned.events.list.into.realm.lengths.rev.cumsum,"[[",i)
    binned.events.list.into.hot.lengths.rev.median[[i]]<-sapply(binned.events.list.into.hot.lengths.rev.cumsum,"[[",i)
    binned.events.list.hot.into.realm.lengths.rev.median[[i]]<-sapply(binned.events.list.hot.into.realm.lengths.rev,"[[",i)
    binned.events.list.realm.into.hot.lengths.rev.median[[i]]<-sapply(binned.events.list.realm.into.hot.lengths.rev,"[[",i)
    binned.events.list.insitu.realm.lengths.rev.median[[i]]<-sapply(binned.events.list.insitu.realm.lengths.rev,"[[",i)
    binned.events.list.insitu.hot.lengths.rev.median[[i]]<-sapply(binned.events.list.insitu.hot.lengths.rev,"[[",i)
    #binned.events.list.hot.local.extinction.lengths.rev.cumsum.CI[[i]]<-sapply(binned.events.list.hot.local.extinction.lengths.rev.cumsum,"[[",i)
    #binned.events.list.realm.local.extinction.lengths.rev.cumsum.CI[[i]]<-sapply(binned.events.list.realm.local.extinction.lengths.rev.cumsum,"[[",i)
  }
  
  binned.events.list.from.realm.lengths.rev.median<-lapply(binned.events.list.from.realm.lengths.rev.median,function(x)quantile(x,c(0.05,0.95,0.50)))
  binned.events.list.from.hot.lengths.rev.median<-lapply(binned.events.list.from.hot.lengths.rev.median,function(x)quantile(x,c(0.05,0.95,0.50)))
  binned.events.list.into.realm.lengths.rev.median<-lapply(binned.events.list.into.realm.lengths.rev.median,function(x)quantile(x,c(0.05,0.95,0.50)))
  binned.events.list.into.hot.lengths.rev.median<-lapply(binned.events.list.into.hot.lengths.rev.median,function(x)quantile(x,c(0.05,0.95,0.50)))
  binned.events.list.hot.into.realm.lengths.rev.median<-lapply(binned.events.list.hot.into.realm.lengths.rev.median,function(x)quantile(x,c(0.05,0.95,0.50)))
  binned.events.list.realm.into.hot.lengths.rev.median<-lapply(binned.events.list.realm.into.hot.lengths.rev.median,function(x)quantile(x,c(0.05,0.95,0.50)))
  binned.events.list.insitu.realm.lengths.rev.median<-lapply(binned.events.list.insitu.realm.lengths.rev.median,function(x)quantile(x,c(0.05,0.95,0.50)))
  binned.events.list.insitu.hot.lengths.rev.median<-lapply(binned.events.list.insitu.hot.lengths.rev.median,function(x)quantile(x,c(0.05,0.95,0.50)))
  
  
  
  #binned.events.list.insitu.hot.lengths.rev.median<-unlist(lapply(binned.events.list.insitu.hot.lengths.rev.median,function(x) median(x)))
  #binned.events.list.insitu.hot.lengths.rev.median<-binned.events.list.insitu.hot.lengths.rev.median[c((length(binned.events.list.insitu.hot.lengths.rev.median)-age):length(binned.events.list.insitu.hot.lengths.rev.median))]
  #binned.events.list.insitu.realm.lengths.rev.median<-unlist(lapply(binned.events.list.insitu.realm.lengths.rev.median,function(x) median(x)))
  #binned.events.list.insitu.realm.lengths.rev.median<-binned.events.list.insitu.realm.lengths.rev.median[c((length(binned.events.list.insitu.realm.lengths.rev.median)-age):length(binned.events.list.insitu.realm.lengths.rev.median))]
  #binned.events.list.hot.into.realm.lengths.rev.median<-unlist(lapply(binned.events.list.hot.into.realm.lengths.rev.median,function(x) median(x)))
  #binned.events.list.hot.into.realm.lengths.rev.median<-binned.events.list.hot.into.realm.lengths.rev.median[c((length(binned.events.list.hot.into.realm.lengths.rev.median)-age):length(binned.events.list.hot.into.realm.lengths.rev.median))]
  #binned.events.list.realm.into.hot.lengths.rev.median<-unlist(lapply(binned.events.list.realm.into.hot.lengths.rev.median,function(x) median(x)))
  #binned.events.list.realm.into.hot.lengths.rev.median<-binned.events.list.realm.into.hot.lengths.rev.median[c((length(binned.events.list.realm.into.hot.lengths.rev.median)-age):length(binned.events.list.realm.into.hot.lengths.rev.median))]
  
  binned.events.list.from.realm.lengths.rev.median<-binned.events.list.from.realm.lengths.rev.median[c((length(binned.events.list.from.realm.lengths.rev.median)-age):length(binned.events.list.from.realm.lengths.rev.median))]
  binned.events.list.from.hot.lengths.rev.median<-binned.events.list.from.hot.lengths.rev.median[c((length(binned.events.list.from.hot.lengths.rev.median)-age):length(binned.events.list.from.hot.lengths.rev.median))]
  binned.events.list.into.realm.lengths.rev.median<-binned.events.list.into.realm.lengths.rev.median[c((length(binned.events.list.into.realm.lengths.rev.median)-age):length(binned.events.list.into.realm.lengths.rev.median))]
  binned.events.list.into.hot.lengths.rev.median<-binned.events.list.into.hot.lengths.rev.median[c((length(binned.events.list.into.hot.lengths.rev.median)-age):length(binned.events.list.into.hot.lengths.rev.median))]
  binned.events.list.hot.into.realm.lengths.rev.median<-binned.events.list.hot.into.realm.lengths.rev.median[c((length(binned.events.list.hot.into.realm.lengths.rev.median)-age):length(binned.events.list.hot.into.realm.lengths.rev.median))]
  binned.events.list.realm.into.hot.lengths.rev.median<-binned.events.list.realm.into.hot.lengths.rev.median[c((length(binned.events.list.realm.into.hot.lengths.rev.median)-age):length(binned.events.list.realm.into.hot.lengths.rev.median))]
  binned.events.list.insitu.realm.lengths.rev.median<-binned.events.list.insitu.realm.lengths.rev.median[c((length(binned.events.list.insitu.realm.lengths.rev.median)-age):length(binned.events.list.insitu.realm.lengths.rev.median))]
  binned.events.list.insitu.hot.lengths.rev.median<-binned.events.list.insitu.hot.lengths.rev.median[c((length(binned.events.list.insitu.hot.lengths.rev.median)-age):length(binned.events.list.insitu.hot.lengths.rev.median))]
  
  binned.events.list.hot.lineages.rev.cumsum.median<-binned.events.list.insitu.hot.lengths.rev.cumsum.median+binned.events.list.into.hot.lengths.rev.cumsum.median-binned.events.list.hot.local.extinction.lengths.rev.cumsum.median
  binned.events.list.hot.lineages.rev.cumsum.CI.up<-binned.events.list.insitu.hot.lengths.rev.cumsum.CI.up+binned.events.list.into.hot.lengths.rev.cumsum.CI.up- binned.events.list.hot.local.extinction.lengths.rev.cumsum.CI.up
  binned.events.list.hot.lineages.rev.cumsum.CI.dw<-binned.events.list.insitu.hot.lengths.rev.cumsum.CI.dw+binned.events.list.into.hot.lengths.rev.cumsum.CI.dw- binned.events.list.hot.local.extinction.lengths.rev.cumsum.CI.dw
  binned.events.list.realm.lineages.rev.cumsum.median<-binned.events.list.insitu.realm.lengths.rev.cumsum.median+binned.events.list.into.realm.lengths.rev.cumsum.median-binned.events.list.realm.local.extinction.lengths.rev.cumsum.median
  binned.events.list.realm.lineages.rev.cumsum.CI.up<-binned.events.list.insitu.realm.lengths.rev.cumsum.CI.up+binned.events.list.into.realm.lengths.rev.cumsum.CI.up- binned.events.list.realm.local.extinction.lengths.rev.cumsum.CI.up
  binned.events.list.realm.lineages.rev.cumsum.CI.dw<-binned.events.list.insitu.realm.lengths.rev.cumsum.CI.dw+binned.events.list.into.realm.lengths.rev.cumsum.CI.dw- binned.events.list.realm.local.extinction.lengths.rev.cumsum.CI.dw
  
  
  rates.into.hot.median<-numeric()
  rates.into.realm.median<-numeric()
  rates.insitu.hot.median<-numeric()
  rates.insitu.realm.median<-numeric()
  rates.hot.into.realm.median<-numeric()
  rates.realm.into.hot.median<-numeric()
  
  rates.into.hot.CI.up<-numeric()
  rates.into.realm.CI.up<-numeric()
  rates.insitu.hot.CI.up<-numeric()
  rates.insitu.realm.CI.up<-numeric()
  rates.hot.into.realm.CI.up<-numeric()
  rates.realm.into.hot.CI.up<-numeric()
  
  rates.into.realm.CI.dw<-numeric()
  rates.into.hot.CI.dw<-numeric()
  rates.insitu.hot.CI.dw<-numeric()
  rates.insitu.realm.CI.dw<-numeric()
  rates.hot.into.realm.CI.dw<-numeric()
  rates.realm.into.hot.CI.dw<-numeric()
  
  
  for(i in 2:length(binned.events.list.hot.lineages.rev.cumsum.median)){
    rates.into.hot.median[i]<-binned.events.list.into.hot.lengths.rev.median[[i]][3]/binned.events.list.hot.lineages.rev.cumsum.median[i-1]
    rates.into.realm.median[i]<-binned.events.list.into.realm.lengths.rev.median[[i]][3]/binned.events.list.realm.lineages.rev.cumsum.median[i-1]
    rates.insitu.hot.median[i]<-binned.events.list.insitu.hot.lengths.rev.median[[i]][3]/binned.events.list.hot.lineages.rev.cumsum.median[i-1]
    rates.insitu.realm.median[i]<-binned.events.list.insitu.realm.lengths.rev.median[[i]][3]/binned.events.list.realm.lineages.rev.cumsum.median[i-1]
    rates.hot.into.realm.median[i]<-binned.events.list.hot.into.realm.lengths.rev.median[[i]][3]/binned.events.list.realm.lineages.rev.cumsum.median[i-1]
    rates.realm.into.hot.median[i]<-binned.events.list.realm.into.hot.lengths.rev.median[[i]][3]/binned.events.list.hot.lineages.rev.cumsum.median[i-1]
    
    rates.into.hot.CI.up[i]<-binned.events.list.into.hot.lengths.rev.median[[i]][2]/binned.events.list.hot.lineages.rev.cumsum.median[i-1]
    rates.into.realm.CI.up[i]<-binned.events.list.into.realm.lengths.rev.median[[i]][2]/binned.events.list.realm.lineages.rev.cumsum.median[i-1]
    rates.insitu.hot.CI.up[i]<-binned.events.list.insitu.hot.lengths.rev.median[[i]][2]/binned.events.list.hot.lineages.rev.cumsum.median[i-1]
    rates.insitu.realm.CI.up[i]<-binned.events.list.insitu.realm.lengths.rev.median[[i]][2]/binned.events.list.realm.lineages.rev.cumsum.median[i-1]
    rates.hot.into.realm.CI.up[i]<-binned.events.list.hot.into.realm.lengths.rev.median[[i]][2]/binned.events.list.realm.lineages.rev.cumsum.median[i-1]
    rates.realm.into.hot.CI.up[i]<-binned.events.list.realm.into.hot.lengths.rev.median[[i]][2]/binned.events.list.hot.lineages.rev.cumsum.median[i-1]
    
    rates.into.hot.CI.dw[i]<-binned.events.list.into.hot.lengths.rev.median[[i]][1]/binned.events.list.hot.lineages.rev.cumsum.median[i-1]
    rates.into.realm.CI.dw[i]<-binned.events.list.into.realm.lengths.rev.median[[i]][1]/binned.events.list.realm.lineages.rev.cumsum.median[i-1]
    rates.insitu.hot.CI.dw[i]<-binned.events.list.insitu.hot.lengths.rev.median[[i]][1]/binned.events.list.hot.lineages.rev.cumsum.median[i-1]
    rates.insitu.realm.CI.dw[i]<-binned.events.list.insitu.realm.lengths.rev.median[[i]][1]/binned.events.list.realm.lineages.rev.cumsum.median[i-1]
    rates.hot.into.realm.CI.dw[i]<-binned.events.list.hot.into.realm.lengths.rev.median[[i]][1]/binned.events.list.realm.lineages.rev.cumsum.median[i-1]
    rates.realm.into.hot.CI.dw[i]<-binned.events.list.realm.into.hot.lengths.rev.median[[i]][1]/binned.events.list.hot.lineages.rev.cumsum.median[i-1]
  }
  rates.into.hot.median<-as.numeric(gsub(rates.into.hot.median,pattern='Inf',replacement=''))
  rates.into.realm.median<-as.numeric(gsub(rates.into.realm.median,pattern='Inf',replacement=''))
  rates.insitu.hot.median<-as.numeric(gsub(rates.insitu.hot.median,pattern='Inf',replacement=''))
  rates.insitu.realm.median<-as.numeric(gsub(rates.insitu.realm.median,pattern='Inf',replacement=''))
  rates.hot.into.realm.median<-as.numeric(gsub(rates.hot.into.realm.median,pattern='Inf',replacement=''))
  rates.realm.into.hot.median<-as.numeric(gsub(rates.realm.into.hot.median,pattern='Inf',replacement=''))
  
  rates.into.hot.CI.up<-as.numeric(gsub(rates.into.hot.CI.up,pattern='Inf',replacement=''))
  rates.into.realm.CI.up<-as.numeric(gsub(rates.into.realm.CI.up,pattern='Inf',replacement=''))
  rates.insitu.hot.CI.up<-as.numeric(gsub(rates.insitu.hot.CI.up,pattern='Inf',replacement=''))
  rates.insitu.realm.CI.up<-as.numeric(gsub(rates.insitu.realm.CI.up,pattern='Inf',replacement=''))
  rates.hot.into.realm.CI.up<-as.numeric(gsub(rates.hot.into.realm.CI.up,pattern='Inf',replacement=''))
  rates.realm.into.hot.CI.up<-as.numeric(gsub(rates.realm.into.hot.CI.up,pattern='Inf',replacement=''))
  
  rates.into.hot.CI.dw<-as.numeric(gsub(rates.into.hot.CI.dw,pattern='Inf',replacement=''))
  rates.into.realm.CI.dw<-as.numeric(gsub(rates.into.realm.CI.dw,pattern='Inf',replacement=''))
  rates.insitu.hot.CI.dw<-as.numeric(gsub(rates.insitu.hot.CI.dw,pattern='Inf',replacement=''))
  rates.insitu.realm.CI.dw<-as.numeric(gsub(rates.insitu.realm.CI.dw,pattern='Inf',replacement=''))
  rates.hot.into.realm.CI.dw<-as.numeric(gsub(rates.hot.into.realm.CI.dw,pattern='Inf',replacement=''))
  rates.realm.into.hot.CI.dw<-as.numeric(gsub(rates.realm.into.hot.CI.dw,pattern='Inf',replacement=''))
  
  rates.into.hot.median[is.na(rates.into.hot.median)]<-0
  rates.into.realm.median[is.na(rates.into.realm.median)]<-0
  rates.insitu.hot.median[is.na(rates.insitu.hot.median)]<-0
  rates.insitu.realm.median[is.na(rates.insitu.realm.median)]<-0
  rates.hot.into.realm.median[is.na(rates.hot.into.realm.median)]<-0
  rates.realm.into.hot.median[is.na(rates.realm.into.hot.median)]<-0
  
  rates.into.hot.CI.up[is.na(rates.into.hot.CI.up)]<-0
  rates.into.realm.CI.up[is.na(rates.into.realm.CI.up)]<-0
  rates.insitu.hot.CI.up[is.na(rates.insitu.hot.CI.up)]<-0
  rates.insitu.realm.CI.up[is.na(rates.insitu.realm.CI.up)]<-0
  rates.hot.into.realm.CI.up[is.na(rates.hot.into.realm.CI.up)]<-0
  rates.realm.into.hot.CI.up[is.na(rates.realm.into.hot.CI.up)]<-0
  
  rates.into.hot.CI.dw[is.na(rates.into.hot.CI.dw)]<-0
  rates.into.realm.CI.dw[is.na(rates.into.realm.CI.dw)]<-0
  rates.insitu.hot.CI.dw[is.na(rates.insitu.hot.CI.dw)]<-0
  rates.insitu.realm.CI.dw[is.na(rates.insitu.realm.CI.dw)]<-0
  rates.hot.into.realm.CI.dw[is.na(rates.hot.into.realm.CI.dw)]<-0
  rates.realm.into.hot.CI.dw[is.na(rates.realm.into.hot.CI.dw)]<-0
  
  #plot without logging
  plot(c(4,4),xlim=c(0,age),ylim=c(0,max(c(rates.insitu.hot.CI.up,rates.insitu.realm.CI.up),na.rm=TRUE)),type='n',xaxt='n',xlab='age',ylab='speciation rates in situ',main=paste(name,' rates',sep=''))
  #plot(c(4,4),xlim=c(0,217),ylim=c(0,log10(max(c(unlist(binned.events.list.into.realm.lengths.rev.cumsum),unlist(binned.events.list.into.hot.lengths.rev.cumsum),unlist(binned.events.list.insitu.realm.lengths.rev.cumsum),unlist(binned.events.list.insitu.hot.lengths.rev.cumsum))))),type='n',xaxt='n',xlab='age',ylab='cumulative number of events',main=paste(name,' all events',sep=''))
  polygon(x=c(c(0:age),c(age:0)),y=c(rates.insitu.hot.CI.dw,rev(rates.insitu.hot.CI.up)),col=adjustcolor( "red", alpha.f = 0.2))
  lines(x=c(0:age),y=rates.insitu.hot.median,col=adjustcolor( "red", alpha.f = 0.9))
  polygon(x=c(c(0:age),c(age:0)),y=c(rates.insitu.realm.CI.dw,rev(rates.insitu.realm.CI.up)),col=adjustcolor( "blue", alpha.f = 0.2))
  lines(x=c(0:age),y=rates.insitu.realm.median,col=adjustcolor( "blue", alpha.f = 0.9))
  legend('topleft',c('rate.insitu.hot','rate.insitu.non.hot.realm'),col=c('red','blue'),lty=1,cex=.7,bty='n')
  axis(1,at=c(0:age),labels=c(-age:0))
  
  plot(c(4,4),xlim=c(0,age),ylim=c(0,max(c(rates.hot.into.realm.CI.up,rates.realm.into.hot.CI.up),na.rm=TRUE)),type='n',xaxt='n',xlab='age',ylab='colonisation rates',main=paste(name,' rates',sep=''))
  #plot(c(4,4),xlim=c(0,217),ylim=c(0,log10(max(c(unlist(binned.events.list.into.realm.lengths.rev.cumsum),unlist(binned.events.list.into.hot.lengths.rev.cumsum),unlist(binned.events.list.insitu.realm.lengths.rev.cumsum),unlist(binned.events.list.insitu.hot.lengths.rev.cumsum))))),type='n',xaxt='n',xlab='age',ylab='cumulative number of events',main=paste(name,' all events',sep=''))
  polygon(x=c(c(0:age),c(age:0)),y=c(rates.hot.into.realm.CI.dw,rev(rates.hot.into.realm.CI.up)),col=adjustcolor( "red", alpha.f = 0.2))
  lines(x=c(0:age),y=rates.hot.into.realm.median,col=adjustcolor( "red", alpha.f = 0.9))
  polygon(x=c(c(0:age),c(age:0)),y=c(rates.realm.into.hot.CI.dw,rev(rates.realm.into.hot.CI.up)),col=adjustcolor( "blue", alpha.f = 0.2))
  lines(x=c(0:age),y=rates.realm.into.hot.median,col=adjustcolor( "blue", alpha.f = 0.9))
  legend('topleft',c('rate.hot.into.realm','rate.realm.into.hot'),col=c('red','blue'),lty=1,cex=.7,bty='n')
  axis(1,at=c(0:age),labels=c(-age:0))
  
  plot(c(4,4),xlim=c(0,age),ylim=c(0,max(c(rates.into.realm.CI.up,rates.into.hot.CI.up),na.rm=TRUE)),type='n',xaxt='n',xlab='age',ylab='colonisation rates',main=paste(name,' rates',sep=''))
  #plot(c(4,4),xlim=c(0,217),ylim=c(0,log10(max(c(unlist(binned.events.list.into.realm.lengths.rev.cumsum),unlist(binned.events.list.into.hot.lengths.rev.cumsum),unlist(binned.events.list.insitu.realm.lengths.rev.cumsum),unlist(binned.events.list.insitu.hot.lengths.rev.cumsum))))),type='n',xaxt='n',xlab='age',ylab='cumulative number of events',main=paste(name,' all events',sep=''))
  polygon(x=c(c(0:age),c(age:0)),y=c(rates.into.realm.CI.dw,rev(rates.into.realm.CI.up)),col=adjustcolor( "purple", alpha.f = 0.2))
  lines(x=c(0:age),y=rates.into.realm.median,col=adjustcolor( "purple", alpha.f = 0.9))
  polygon(x=c(c(0:age),c(age:0)),y=c(rates.into.hot.CI.dw,rev(rates.into.hot.CI.up)),col=adjustcolor( "orange", alpha.f = 0.2))
  lines(x=c(0:age),y=rates.into.hot.median,col=adjustcolor( "orange", alpha.f = 0.9))
  legend('topleft',c('rate.into.realm','rate.into.hot'),col=c('purple','orange'),lty=1,cex=.7,bty='n')
  axis(1,at=c(0:age),labels=c(-age:0))
  
  
  
}

plot_insitu_cumulative_events_BSM_old<-function(results.BSM,name){
  results.BSM<-results.BSM[[2]][[2]]
  events.list<-list()
  for (i in 1:length(results.BSM)){
    events.list[[i]]<-list()
    for (a in 1:length(results.BSM[[i]]$RES_clado_events_tables)){
      clado.df<-results.BSM[[i]]$RES_clado_events_tables[[a]]
      clado.df<-clado.df[,c('time_bp','clado_event_txt')]
      clado.df<-clado.df[clado.df$clado_event_txt!='',]
      ana.df<-results.BSM[[i]]$RES_ana_events_tables[[a]]
      ana.df<-ana.df[,c('abs_event_time','event_txt')]
      ana.df<-ana.df[ana.df$abs_event_time!='',]
      colnames(clado.df)<-c('time_bp','event_txt')
      colnames(ana.df)<-c('time_bp','event_txt')
      events.list[[i]][[a]]<-rbind(clado.df,ana.df)
    }
  }
  event.type.table<-list()
  for (i in 1:length(events.list)){
    event.type.table[[i]]<-lapply(events.list[[i]],function(x)table(x$event_txt))
  }
  insitu.realm<-list()
  insitu.hot<-list()
  if(name=='afrotrop'){
    realm.code<-'A'
  }else if(name=='austral'){
    realm.code<-'B'
  }else if(name=='indo'){
    realm.code<-'C'
  }else if(name=='nearctic'){
    realm.code<-'D'
  }else if(name=='neotrop'){
    realm.code<-'E'
  }else if(name=='palearctic'){
    realm.code<-'F'
  }
  hot.code<-'G'
  for (i in 1:length(event.type.table)){
    insitu.realm[[i]]<-unlist(lapply(event.type.table[[i]],function(x)x[paste(realm.code,'->',realm.code,',',realm.code,sep='')]))
    insitu.hot[[i]]<-unlist(lapply(event.type.table[[i]],function(x)x[paste(hot.code,'->',hot.code,',',hot.code,sep='')]))
    
  }
  insitu.realm<-lapply(insitu.realm,function(x) replace(x,is.na(x),0))
  insitu.hot<-lapply(insitu.hot,function(x) replace(x,is.na(x),0))
  max<-max(unlist(insitu.realm),unlist(insitu.hot))
  for (i in 1:length(insitu.hot)){
    hist(insitu.hot[[i]],xlim=c(0,max),col=adjustcolor('red', alpha.f = 0.3),ylim=c(0,30),main=paste(name,'.In.Situ.Speciation.BSM',sep=''),xlab='number of events')
    par(new=TRUE)
    hist(insitu.realm[[i]],xlim=c(0,max),col=adjustcolor('blue', alpha.f = 0.3),ylim=c(0,30),main='',xlab='',ylab='')
    legend('topleft',c('hot','non.hot.realm'),col=c('red','blue'),lty=1,cex=.7,bty='n')
  }
  #A->AB; B->AB
  realm.to.hot<-list()
  hot.to.realm<-list()
  for (i in 1:length(event.type.table)){
    ####CHECK THAT THESE ARE THE ONLY TWO NOTATIONS POSSIBLE FOR DISPERSAL
    realm.to.hot[[i]]<-unlist(lapply(event.type.table[[i]],function(x)x[paste(realm.code,'->',realm.code,hot.code,sep='')]))
    hot.to.realm[[i]]<-unlist(lapply(event.type.table[[i]],function(x)x[paste(hot.code,'->',realm.code,hot.code,sep='')]))
    
  }
  realm.to.hot<-lapply(realm.to.hot,function(x) replace(x,is.na(x),0))
  hot.to.realm<-lapply(hot.to.realm,function(x) replace(x,is.na(x),0))
  max<-max(unlist(realm.to.hot),unlist(hot.to.realm))
  for (i in 1:length(hot.to.realm)){
    hist(hot.to.realm[[i]],xlim=c(0,max),col=adjustcolor('red', alpha.f = 0.3),ylim=c(0,30),main=paste(name,'.Dispersal.events',sep=''),xlab='number of events')
    par(new=TRUE)
    hist(realm.to.hot[[i]],xlim=c(0,max),col=adjustcolor('blue', alpha.f = 0.3),ylim=c(0,30),main='',xlab='',ylab='')
    legend('topleft',c('hot.to.non.hot.realm','non.hot.realm.to.hot'),col=c('red','blue'),lty=1,cex=.7,bty='n')
    
  }
  ###########cumulative sum of events plot
  #process a results.BSM to extract events and times in a list of lists (each element of the list is one clade with a list of 50 replicate BSMs)
  events.list<-list()
  for (i in 1:length(results.BSM)){
    events.list[[i]]<-list()
    for (a in 1:length(results.BSM[[i]]$RES_clado_events_tables)){
      clado.df<-results.BSM[[i]]$RES_clado_events_tables[[a]]
      clado.df<-clado.df[,c('time_bp','clado_event_txt')]
      clado.df<-clado.df[clado.df$clado_event_txt!='',]
      ana.df<-results.BSM[[i]]$RES_ana_events_tables[[a]]
      ana.df<-ana.df[,c('abs_event_time','event_txt')]
      ana.df<-ana.df[ana.df$abs_event_time!='',]
      colnames(clado.df)<-c('time_bp','event_txt')
      colnames(ana.df)<-c('time_bp','event_txt')
      events.list[[i]][[a]]<-rbind(clado.df,ana.df)
    }
  }
  #create a df that bins events into 1 myr slots
  #get the roots of all clades
  clade.roots<-vector()
  for (i in 1:length(results.BSM)){
    clade.roots[i]<-results.BSM[[i]]$RES_clado_events_tables[[1]][results.BSM[[i]]$RES_clado_events_tables[[1]]$node.type=='root',]$time_bp
  }
  #another way of binning events, a list of lists (binning events for each BSM replicate)
  binned.events.list<-vector('list',length=length(events.list[[1]]))
  cat('binning events','\n')
  for (i in 1:length(events.list[[1]])){
   # cat(i,'\n')
    binned.events.list[[i]]<-vector('list',length=length(seq(0,floor(clade.roots)))+1)
    for (b in 1:nrow(events.list[[1]][[i]])){
     # cat(b,'\n')
      time_bp<-floor(events.list[[1]][[i]][b,'time_bp'])
      event<-events.list[[1]][[i]][b,'event_txt']
      names(event)<-row.names(events.list[[1]][[i]][b,])
      if(is.null((binned.events.list[[i]][[time_bp+1]]))){
        binned.events.list[[i]][[time_bp+1]]<-c('0',event)
      }else{
        binned.events.list[[i]][[time_bp+1]]<-c(binned.events.list[[i]][[time_bp+1]],event)  
      }
    }
  }
  binned.events.list.insitu.realm<-list()
  binned.events.list.insitu.hot<-list()
  cat('subsetting binned events to realm','\n')
  for(i in 1:length(binned.events.list)){
    #cat(i,'\n')
    binned.events.list.insitu.realm[[i]]<-vector('list',length=length(seq(0,floor(clade.roots)))+1)
    binned.events.list.insitu.hot[[i]]<-vector('list',length=length(seq(0,floor(clade.roots)))+1)
    for (a in 1:length(binned.events.list[[i]])){
     # cat(a,'\n')
      binned.events.list.insitu.realm[[i]][[a]]<-binned.events.list[[i]][[a]][grep(paste(realm.code,'->',realm.code,',',realm.code,sep=''),binned.events.list[[i]][[a]])]
      binned.events.list.insitu.hot[[i]][[a]]<-binned.events.list[[i]][[a]][grep(paste(hot.code,'->',hot.code,',',hot.code,sep=''),binned.events.list[[i]][[a]])]
    }
    
  }
  binned.events.list.insitu.realm.lengths<-lapply(binned.events.list.insitu.realm,function(x)unlist(lapply(x,function(x) length(x))))
  binned.events.list.insitu.hot.lengths<-lapply(binned.events.list.insitu.hot,function(x)unlist(lapply(x,function(x) length(x))))
  
  binned.events.list.insitu.realm.lengths.rev<-lapply(binned.events.list.insitu.realm.lengths,function(x) rev(x))
  binned.events.list.insitu.hot.lengths.rev<-lapply(binned.events.list.insitu.hot.lengths,function(x) rev(x))
  
  binned.events.list.insitu.realm.lengths.rev.20myr<-lapply(binned.events.list.insitu.realm.lengths.rev,function(x) x[c((length(x)-20):length(x))])
  binned.events.list.insitu.hot.lengths.rev.20myr<-lapply(binned.events.list.insitu.hot.lengths.rev,function(x) x[c((length(x)-20):length(x))])
  
  binned.events.list.insitu.realm.lengths.rev.20myr.cumsum<-lapply(binned.events.list.insitu.realm.lengths.rev.20myr,function(x) cumsum(x))
  binned.events.list.insitu.hot.lengths.rev.20myr.cumsum<-lapply(binned.events.list.insitu.hot.lengths.rev.20myr,function(x) cumsum(x))
  
  
  
  
  plot(c(4,4),xlim=c(0,21),ylim=c(0,max(c(unlist(binned.events.list.insitu.realm.lengths.rev.20myr.cumsum),unlist(binned.events.list.insitu.hot.lengths.rev.20myr.cumsum)))),type='n',xaxt='n',xlab='age',ylab='cumulative number of events',main=paste(name,' in situ speciation',sep=''))
  for (i in 1:length(binned.events.list.insitu.realm.lengths.rev.20myr.cumsum)){
    lines(x=c(1:21),binned.events.list.insitu.realm.lengths.rev.20myr.cumsum[[i]],col= adjustcolor( "blue", alpha.f = 0.2))
  }
  axis(1,at=c(1:20),labels=c(-20:-1))
  for (i in 1:length(binned.events.list.insitu.hot.lengths.rev.20myr.cumsum)){
    lines(x=c(1:21),binned.events.list.insitu.hot.lengths.rev.20myr.cumsum[[i]],col= adjustcolor( "red", alpha.f = 0.2))
  }
  legend('topleft',c('hot','non.hot.realm'),col=c('red','blue'),lty=1,cex=.7,bty='n')
  plot(c(4,4),xlim=c(0,21),ylim=c(0,max(c(unlist(binned.events.list.insitu.realm.lengths.rev.20myr.cumsum)))),type='n',xaxt='n',xlab='age',ylab='cumulative number of events',main=paste(name,' in situ non.hot.realm speciation',sep=''))
  for (i in 1:length(binned.events.list.insitu.realm.lengths.rev.20myr.cumsum)){
    lines(x=c(1:21),binned.events.list.insitu.realm.lengths.rev.20myr.cumsum[[i]],col= adjustcolor( "blue", alpha.f = 0.2))
  }
  axis(1,at=c(1:20),labels=c(-20:-1))
  
  plot(c(4,4),xlim=c(0,21),ylim=c(0,max(c(unlist(binned.events.list.insitu.hot.lengths.rev.20myr.cumsum)))),type='n',xaxt='n',xlab='age',ylab='cumulative number of events',main=paste(name,' in situ hot speciation',sep=''))
  for (i in 1:length(binned.events.list.insitu.hot.lengths.rev.20myr.cumsum)){
    lines(x=c(1:21),binned.events.list.insitu.hot.lengths.rev.20myr.cumsum[[i]],col= adjustcolor( "red", alpha.f = 0.2))
  }
  axis(1,at=c(1:20),labels=c(-20:-1))
  
  
}


get_BSM_BAMM_rates_events_object<-function(results.BSM,name,speciation.rate.treefile){
  cat('getting speciation rates','\n')
  speciation.rate.tree<-read.tree(speciation.rate.treefile)  
  table.speciation.rate<-prt(speciation.rate.tree,printflag = FALSE)
  area.codes<-c('O','A','B','C','D','E','F','G','AB','AC','AD','AE','AF','AG','BC','BD','BE','BF','BG','CD','CE','CF','CG','DE','DF','DG','EF','EG','FG')
  #get event names from BSM
  if(name=='afrotrop'){
    realm.code<-'A'
  }else if(name=='austral'){
    realm.code<-'B'
  }else if(name=='indo'){
    realm.code<-'C'
  }else if(name=='nearctic'){
    realm.code<-'D'
  }else if(name=='neotrop'){
    realm.code<-'E'
  }else if(name=='palearctic'){
    realm.code<-'F'
  }
  hot.code<-'G'
  #read clado and ana events table of each BSM replicate
  insitu.realm.rates.BSM<-list()
  insitu.hot.rates.BSM<-list()
  insitu.realm.rates.BAMM.strict<-list()
  insitu.hot.rates.BAMM.strict<-list()
  insitu.realm.rates.BAMM<-list()
  insitu.hot.rates.BAMM<-list()
  dispersal.into.realm.rates.BSM<-list()
  dispersal.into.hot.rates.BSM<-list()
  hot.into.realm.rates.BSM<-list()
  realm.into.hot.rates.BSM<-list()
  
  dispersal.into.hot.time.bins.counts.rev.cumsum<-list()
  dispersal.into.realm.time.bins.counts.rev.cumsum<-list()
  hot.into.realm.time.bins.counts.rev.cumsum<-list()
  realm.into.hot.time.bins.counts.rev.cumsum<-list()
  insitu.hot.time.bins.counts.rev.cumsum<-list()
  insitu.realm.time.bins.counts.rev.cumsum<-list()
  dispersal.from.hot.time.bins.counts.rev.cumsum<-list()
  dispersal.from.realm.time.bins.counts.rev.cumsum<-list()
  realm.into.hot.time.bins.counts.rev.cumsum<-list()
  into.into.realm.time.bins.counts.rev.cumsum<-list()
  realm.lineages.counts.BSM<-list()
  hot.lineages.counts.BSM<-list()
  cat('analysing BSMs','\n')
  for (i in 1:length(results.BSM$RES_clado_events_tables)){
    ###BioGeoBEARS
    #cat(i,'\n')
    a1<-results.BSM$RES_clado_events_tables[[i]]
    a2<-results.BSM$RES_ana_events_tables[[i]]
    #assign events to time bin
    a1$time.bin<-ceiling(a1$time_bp)
    a1[a1$time.bin==0,'time.bin']<-1
    a2$time.bin<-ceiling(a2$time_bp)
    a2[a2$time.bin==0,'time.bin']<-1
    max.age<-max(a1$time.bin)
    #get time bins of each event type
    #dispersal into hot
    dispersal.into.hot.time.bins<-sort(c(a1[a1$clado_dispersal_to==hot.code,'time.bin'],a2[a2$dispersal_to==hot.code,'time.bin']))
    dispersal.into.hot.time.bins.counts<-sapply(c(1:max.age),function(x)length(dispersal.into.hot.time.bins[dispersal.into.hot.time.bins%in%x]))
    dispersal.into.hot.time.bins.counts.rev<-rev(dispersal.into.hot.time.bins.counts)
    dispersal.into.hot.time.bins.counts.rev.cumsum[[i]]<-cumsum(dispersal.into.hot.time.bins.counts.rev)
    #dispersal into realm
    dispersal.into.realm.time.bins<-sort(c(a1[a1$clado_dispersal_to==realm.code,'time.bin'],a2[a2$dispersal_to==realm.code,'time.bin']))
    dispersal.into.realm.time.bins.counts<-sapply(c(1:max.age),function(x)length(dispersal.into.realm.time.bins[dispersal.into.realm.time.bins%in%x]))
    dispersal.into.realm.time.bins.counts.rev<-rev(dispersal.into.realm.time.bins.counts)
    dispersal.into.realm.time.bins.counts.rev.cumsum[[i]]<-cumsum(dispersal.into.realm.time.bins.counts.rev)
    #hot local extinctions
    hot.local.extinctions.time.bins<-sort(a2[a2$extirpation_from==hot.code,'time.bin'])
    hot.local.extinctions.time.bins.counts<-sapply(c(1:max.age),function(x)length(hot.local.extinctions.time.bins[hot.local.extinctions.time.bins%in%x]))
    hot.local.extinctions.time.bins.counts.rev<-rev(hot.local.extinctions.time.bins.counts)
    hot.local.extinctions.time.bins.counts.rev.cumsum<-cumsum(hot.local.extinctions.time.bins.counts.rev)
    #realm local extinctions
    realm.local.extinctions.time.bins<-sort(a2[a2$extirpation_from==realm.code,'time.bin'])
    realm.local.extinctions.time.bins.counts<-sapply(c(1:max.age),function(x)length(realm.local.extinctions.time.bins[realm.local.extinctions.time.bins%in%x]))
    realm.local.extinctions.time.bins.counts.rev<-rev(realm.local.extinctions.time.bins.counts)
    realm.local.extinctions.time.bins.counts.rev.cumsum<-cumsum(realm.local.extinctions.time.bins.counts.rev)
    #hot into realm dispersal
    hot.into.realm.time.bins<-sort(c(a2[intersect(intersect(grep(hot.code,a2$current_rangetxt),grep(realm.code,a2$current_rangetxt,invert=TRUE)),grep(realm.code,a2$new_rangetxt)),'time.bin'],a1[a1$sampled_states_AT_nodes%in%intersect(grep(hot.code,area.codes),grep(realm.code,area.codes,invert=TRUE))&a1$clado_dispersal_to==realm.code,'time.bin']))
    hot.into.realm.time.bins.counts<-sapply(c(1:max.age),function(x)length(hot.into.realm.time.bins[hot.into.realm.time.bins%in%x]))
    hot.into.realm.time.bins.counts.rev<-rev(hot.into.realm.time.bins.counts)
    hot.into.realm.time.bins.counts.rev.cumsum[[i]]<-cumsum(hot.into.realm.time.bins.counts.rev)
    #realm into hot dispersal
    realm.into.hot.time.bins<-sort(c(a2[intersect(intersect(grep(realm.code,a2$current_rangetxt),grep(hot.code,a2$current_rangetxt,invert=TRUE)),grep(hot.code,a2$new_rangetxt)),'time.bin'],a1[a1$sampled_states_AT_nodes%in%intersect(grep(realm.code,area.codes),grep(hot.code,area.codes,invert=TRUE))&a1$clado_dispersal_to==hot.code,'time.bin']))
    realm.into.hot.time.bins.counts<-sapply(c(1:max.age),function(x)length(realm.into.hot.time.bins[realm.into.hot.time.bins%in%x]))
    realm.into.hot.time.bins.counts.rev<-rev(realm.into.hot.time.bins.counts)
    realm.into.hot.time.bins.counts.rev.cumsum[[i]]<-cumsum(realm.into.hot.time.bins.counts.rev)
    #in situ hot cladogenesis
    insitu.hot.time.bins<-sort(a1[a1$clado_dispersal_to==''&(a1$sampled_states_AT_nodes%in%grep(hot.code,area.codes))&a1$clado_event_txt!='','time.bin'])
    insitu.hot.time.bins.counts<-sapply(c(1:max.age),function(x)length(insitu.hot.time.bins[insitu.hot.time.bins%in%x]))
    insitu.hot.time.bins.counts.rev<-rev(insitu.hot.time.bins.counts)
    insitu.hot.time.bins.counts.rev.cumsum[[i]]<-cumsum(insitu.hot.time.bins.counts.rev)
    #in situ realm cladogenesis
    insitu.realm.time.bins<-sort(a1[a1$clado_dispersal_to==''&(a1$sampled_states_AT_nodes%in%grep(realm.code,area.codes))&a1$clado_event_txt!='','time.bin'])
    insitu.realm.time.bins.counts<-sapply(c(1:max.age),function(x)length(insitu.realm.time.bins[insitu.realm.time.bins%in%x]))
    insitu.realm.time.bins.counts.rev<-rev(insitu.realm.time.bins.counts)
    insitu.realm.time.bins.counts.rev.cumsum[[i]]<-cumsum(insitu.realm.time.bins.counts.rev)
    #dispersal from hot
    dispersal.from.hot.time.bins<-sort(c(a2[a2$extirpation_from=='-'&(a2$current_rangetxt%in%area.codes[grep(hot.code,area.codes)]),'time.bin'],a1[(a1$clado_event_txt!='')&(a1$sampled_states_AT_nodes%in%grep(hot.code,area.codes))&(!(a1$samp_LEFT_dcorner%in%grep(hot.code,area.codes))|!(a1$samp_RIGHT_dcorner%in%grep(hot.code,area.codes))),'time.bin']))
    dispersal.from.hot.time.bins.counts<-sapply(c(1:max.age),function(x)length(dispersal.from.hot.time.bins[dispersal.from.hot.time.bins%in%x]))
    dispersal.from.hot.time.bins.counts.rev<-rev(dispersal.from.hot.time.bins.counts)
    dispersal.from.hot.time.bins.counts.rev.cumsum[[i]]<-cumsum(dispersal.from.hot.time.bins.counts.rev)
    #dispersal from realm
    dispersal.from.realm.time.bins<-sort(c(a2[a2$extirpation_from=='-'&(a2$current_rangetxt%in%area.codes[grep(realm.code,area.codes)]),'time.bin'],a1[(a1$clado_event_txt!='')&(a1$sampled_states_AT_nodes%in%grep(realm.code,area.codes))&(!(a1$samp_LEFT_dcorner%in%grep(realm.code,area.codes))|!(a1$samp_RIGHT_dcorner%in%grep(realm.code,area.codes))),'time.bin']))
    dispersal.from.realm.time.bins.counts<-sapply(c(1:max.age),function(x)length(dispersal.from.realm.time.bins[dispersal.from.realm.time.bins%in%x]))
    dispersal.from.realm.time.bins.counts.rev<-rev(dispersal.from.realm.time.bins.counts)
    dispersal.from.realm.time.bins.counts.rev.cumsum[[i]]<-cumsum(dispersal.from.realm.time.bins.counts.rev)
    #calculate number of lineages in hot and realm at each time bin
    realm.lineages.counts<-insitu.realm.time.bins.counts.rev.cumsum[[i]]+dispersal.into.realm.time.bins.counts.rev.cumsum[[i]]-realm.local.extinctions.time.bins.counts.rev.cumsum
    hot.lineages.counts<-insitu.hot.time.bins.counts.rev.cumsum[[i]]+dispersal.into.hot.time.bins.counts.rev.cumsum[[i]]-hot.local.extinctions.time.bins.counts.rev.cumsum
    realm.lineages.counts.BSM[[i]]<-list()
    hot.lineages.counts.BSM[[i]]<-list()
    realm.lineages.counts.BSM[[i]]<-realm.lineages.counts
    hot.lineages.counts.BSM[[i]]<-hot.lineages.counts
    #calculate in.situ.rates.BSM
    insitu.realm.rates.BSM[[i]]<-unlist(lapply(2:length(realm.lineages.counts),function(x)insitu.realm.time.bins.counts.rev[x]/realm.lineages.counts[x-1]))
    insitu.hot.rates.BSM[[i]]<-unlist(lapply(2:length(hot.lineages.counts),function(x)insitu.hot.time.bins.counts.rev[x]/hot.lineages.counts[x-1]))
    #calculate colonisation rates BSM
    dispersal.into.realm.rates.BSM[[i]]<-unlist(lapply(2:length(realm.lineages.counts),function(x)dispersal.into.realm.time.bins.counts.rev[x]/realm.lineages.counts[x-1]))
    dispersal.into.hot.rates.BSM[[i]]<-unlist(lapply(2:length(hot.lineages.counts),function(x)dispersal.into.hot.time.bins.counts.rev[x]/hot.lineages.counts[x-1]))
    #calculate dispersal rates hot->realm and viceversa BSM
    hot.into.realm.rates.BSM[[i]]<-unlist(lapply(2:length(realm.lineages.counts),function(x)hot.into.realm.time.bins.counts.rev[x]/realm.lineages.counts[x-1]))
    realm.into.hot.rates.BSM[[i]]<-unlist(lapply(2:length(hot.lineages.counts),function(x)realm.into.hot.time.bins.counts.rev[x]/hot.lineages.counts[x-1]))
    
    ###BAMM
    #assign regions to each node
    regions<-vector()
    for(a in 1:nrow(a1)){
      if(a1[a,'node.type']=='root'){
        regions[a]<-area.codes[a1[a,'sampled_states_AT_nodes']]
      }else if((a1[a,'sampled_states_AT_nodes']==a1[a,'sampled_states_AT_brbots'])&(a1[a,'anagenetic_events_txt_below_node']=='none')){
        regions[a]<-area.codes[a1[a,'sampled_states_AT_nodes']]
      }else{
        regions[a]<-NA
      }#this ignores anagenetic changes, only assigns rates to branches that do not change state
      #else{
      #   new.range<-strsplit(a1[a,'anagenetic_events_txt_below_node'],split=',')[[1]][grep('new_rangetxt',strsplit(a1[a,'anagenetic_events_txt_below_node'],split=',')[[1]])[length(grep('new_rangetxt',strsplit(a1[a,'anagenetic_events_txt_below_node'],split=',')[[1]]))]]
      #   regions[a]<-gsub(new.range,pattern='new_rangetxt:',replacement='')
      #  }
    }
    a1$speciation.rate<-table.speciation.rate$edge.length
    a1$region<-regions
    #assign events to time bin
    a1$time.bin<-ceiling(a1$time_bp)
    a1[a1$time.bin==0,'time.bin']<-1
    #get mean rates per region (strict = species only can occur in a single region)
    if(length(a1$region[a1$region%in%c(realm.code)])>1){
      insitu.realm.rates.BAMM.strict[[i]]<-aggregate(speciation.rate~time.bin,data=a1[a1$region%in%c(realm.code),],function(x)mean(x))
    }else if (length(a1$region[a1$region%in%c(realm.code)])==1){
      insitu.realm.rates.BAMM.strict[[i]]<-as.data.frame(cbind(a1[a1$region%in%c(realm.code),'time.bin'],a1[a1$region%in%c(realm.code),'speciation.rate']),stringsAsFactors = F)
      colnames(insitu.realm.rates.BAMM.strict[[i]])<-c('time.bin','speciation.rate')
    }else{
      insitu.realm.rates.BAMM.strict[[i]]<-as.data.frame(cbind(1,0),stringsAsFactors = F)
      colnames(insitu.realm.rates.BAMM.strict[[i]])<-c('time.bin','speciation.rate')
    }
    if(length(a1$region[a1$region%in%c(hot.code)])>1){
      insitu.hot.rates.BAMM.strict[[i]]<-aggregate(speciation.rate~time.bin,data=a1[a1$region%in%c(hot.code),],function(x)mean(x))
    }else if (length(a1$region[a1$region%in%c(hot.code)])==1){
      insitu.hot.rates.BAMM.strict[[i]]<-as.data.frame(cbind(a1[a1$region%in%c(hot.code),'time.bin'],a1[a1$region%in%c(hot.code),'speciation.rate']),stringsAsFactors = F)
      colnames(insitu.hot.rates.BAMM.strict[[i]])<-c('time.bin','speciation.rate')
    }else{
      insitu.hot.rates.BAMM.strict[[i]]<-as.data.frame(cbind(1,0),stringsAsFactors = F)
      colnames(insitu.hot.rates.BAMM.strict[[i]])<-c('time.bin','speciation.rate')
    }
    if(length(a1$region[a1$region%in%c(area.codes[grep(realm.code,area.codes)])])>1){
      insitu.realm.rates.BAMM[[i]]<-aggregate(speciation.rate~time.bin,data=a1[a1$region%in%c(area.codes[grep(realm.code,area.codes)]),],function(x)mean(x))
    }else if (length(a1$region[a1$region%in%c(area.codes[grep(realm.code,area.codes)])])==1){
      insitu.realm.rates.BAMM[[i]]<-as.data.frame(cbind(a1[a1$region%in%c(area.codes[grep(realm.code,area.codes)]),'time.bin'],a1[a1$region%in%c(area.codes[grep(realm.code,area.codes)]),'speciation.rate']),stringsAsFactors = F)
      colnames(insitu.realm.rates.BAMM[[i]])<-c('time.bin','speciation.rate')
    }else{
      insitu.realm.rates.BAMM[[i]]<-as.data.frame(cbind(1,0),stringsAsFactors = F)
      colnames(insitu.realm.rates.BAMM[[i]])<-c('time.bin','speciation.rate')
    }
    if(length(a1$region[a1$region%in%c(area.codes[grep(hot.code,area.codes)])])>1){
      insitu.hot.rates.BAMM[[i]]<-aggregate(speciation.rate~time.bin,data=a1[a1$region%in%c(area.codes[grep(hot.code,area.codes)]),],function(x)mean(x))
    }else if (length(a1$region[a1$region%in%c(area.codes[grep(hot.code,area.codes)])])==1){
      insitu.hot.rates.BAMM[[i]]<-as.data.frame(cbind(a1[a1$region%in%c(area.codes[grep(hot.code,area.codes)]),'time.bin'],a1[a1$region%in%c(area.codes[grep(hot.code,area.codes)]),'speciation.rate']),stringsAsFactors = F)
      colnames(insitu.hot.rates.BAMM[[i]])<-c('time.bin','speciation.rate')
    }else{
      insitu.hot.rates.BAMM[[i]]<-as.data.frame(cbind(1,0),stringsAsFactors = F)
      colnames(insitu.hot.rates.BAMM[[i]])<-c('time.bin','speciation.rate')
    }
    
    
  }
  #add 0s to empty bins
  insitu.realm.rates.BAMM<-lapply(insitu.realm.rates.BAMM,function(x) {missing.time.bins<-c(1:length(insitu.hot.rates.BSM[[1]]))[!c(1:length(insitu.hot.rates.BSM[[1]])%in%x$time.bin)];missing.df<-as.data.frame(cbind(missing.time.bins,rep(0,length(missing.time.bins))));colnames(missing.df)<-c('time.bin','speciation.rate');new.df<-rbind(x,missing.df);new.df<-new.df[order(new.df$time.bin),];return(rev(new.df$speciation.rate))})
  insitu.hot.rates.BAMM<-lapply(insitu.hot.rates.BAMM,function(x) {missing.time.bins<-c(1:length(insitu.hot.rates.BSM[[1]]))[!c(1:length(insitu.hot.rates.BSM[[1]])%in%x$time.bin)];missing.df<-as.data.frame(cbind(missing.time.bins,rep(0,length(missing.time.bins))));colnames(missing.df)<-c('time.bin','speciation.rate');new.df<-rbind(x,missing.df);new.df<-new.df[order(new.df$time.bin),];return(rev(new.df$speciation.rate))})
  #add 0s to empty bins
  insitu.realm.rates.BAMM.strict<-lapply(insitu.realm.rates.BAMM.strict,function(x) {missing.time.bins<-c(1:length(insitu.hot.rates.BSM[[1]]))[!c(1:length(insitu.hot.rates.BSM[[1]])%in%x$time.bin)];missing.df<-as.data.frame(cbind(missing.time.bins,rep(0,length(missing.time.bins))));colnames(missing.df)<-c('time.bin','speciation.rate');new.df<-rbind(x,missing.df);new.df<-new.df[order(new.df$time.bin),];return(rev(new.df$speciation.rate))})
  insitu.hot.rates.BAMM.strict<-lapply(insitu.hot.rates.BAMM.strict,function(x) {missing.time.bins<-c(1:length(insitu.hot.rates.BSM[[1]]))[!c(1:length(insitu.hot.rates.BSM[[1]])%in%x$time.bin)];missing.df<-as.data.frame(cbind(missing.time.bins,rep(0,length(missing.time.bins))));colnames(missing.df)<-c('time.bin','speciation.rate');new.df<-rbind(x,missing.df);new.df<-new.df[order(new.df$time.bin),];return(rev(new.df$speciation.rate))})
  
  #get the BSM rates CIs to draw on plot
  insitu.realm.rates.BSM.CI<-list()
  insitu.hot.rates.BSM.CI<-list()
  for (i in 1:length(insitu.realm.rates.BSM[[1]])){
    insitu.realm.rates.BSM.CI[[i]]<-sapply(insitu.realm.rates.BSM,"[[",i)
    insitu.hot.rates.BSM.CI[[i]]<-sapply(insitu.hot.rates.BSM,"[[",i)
  }
  insitu.realm.rates.BSM.CI<-lapply(insitu.realm.rates.BSM.CI,function(x)quantile(x,c(0.05,0.95,0.50),na.rm=TRUE))
  insitu.hot.rates.BSM.CI<-lapply(insitu.hot.rates.BSM.CI,function(x)quantile(x,c(0.05,0.95,0.50),na.rm=TRUE))
  #get the BSM dispersal rates CIs to draw on plot
  hot.into.realm.rates.BSM.CI<-list()
  realm.into.hot.rates.BSM.CI<-list()
  for (i in 1:length(hot.into.realm.rates.BSM[[1]])){
    hot.into.realm.rates.BSM.CI[[i]]<-sapply(hot.into.realm.rates.BSM,"[[",i)
    realm.into.hot.rates.BSM.CI[[i]]<-sapply(realm.into.hot.rates.BSM,"[[",i)
  }
  hot.into.realm.rates.BSM.CI<-lapply(hot.into.realm.rates.BSM.CI,function(x)quantile(x,c(0.05,0.95,0.50),na.rm=TRUE))
  realm.into.hot.rates.BSM.CI<-lapply(realm.into.hot.rates.BSM.CI,function(x)quantile(x,c(0.05,0.95,0.50),na.rm=TRUE))
  dispersal.into.realm.rates.BSM.CI<-list()
  dispersal.into.hot.rates.BSM.CI<-list()
  for (i in 1:length(dispersal.into.realm.rates.BSM[[1]])){
    dispersal.into.realm.rates.BSM.CI[[i]]<-sapply(dispersal.into.realm.rates.BSM,"[[",i)
    dispersal.into.hot.rates.BSM.CI[[i]]<-sapply(dispersal.into.hot.rates.BSM,"[[",i)
  }
  dispersal.into.realm.rates.BSM.CI<-lapply(dispersal.into.realm.rates.BSM.CI,function(x)quantile(x,c(0.05,0.95,0.50),na.rm=TRUE))
  dispersal.into.hot.rates.BSM.CI<-lapply(dispersal.into.hot.rates.BSM.CI,function(x)quantile(x,c(0.05,0.95,0.50),na.rm=TRUE))
  #get the BAMM strict rates CIs to draw on plot
  
  insitu.realm.rates.strict.BAMM.CI<-list()
  insitu.hot.rates.strict.BAMM.CI<-list()
  for (i in 1:length(insitu.realm.rates.BAMM.strict[[1]])){
    insitu.realm.rates.strict.BAMM.CI[[i]]<-sapply(insitu.realm.rates.BAMM.strict,"[[",i)
    insitu.hot.rates.strict.BAMM.CI[[i]]<-sapply(insitu.hot.rates.BAMM.strict,"[[",i)
  }
  insitu.realm.rates.strict.BAMM.CI<-lapply(insitu.realm.rates.strict.BAMM.CI,function(x)quantile(x,c(0.05,0.95,0.50),na.rm=TRUE))
  insitu.hot.rates.strict.BAMM.CI<-lapply(insitu.hot.rates.strict.BAMM.CI,function(x)quantile(x,c(0.05,0.95,0.50),na.rm=TRUE))
  #get the BAMM rates CIs to draw on plot
  insitu.realm.rates.BAMM.CI<-list()
  insitu.hot.rates.BAMM.CI<-list()
  for (i in 1:length(insitu.realm.rates.BAMM[[1]])){
    insitu.realm.rates.BAMM.CI[[i]]<-sapply(insitu.realm.rates.BAMM,"[[",i)
    insitu.hot.rates.BAMM.CI[[i]]<-sapply(insitu.hot.rates.BAMM,"[[",i)
  }
  insitu.realm.rates.BAMM.CI<-lapply(insitu.realm.rates.BAMM.CI,function(x)quantile(x,c(0.05,0.95,0.50),na.rm=TRUE))
  insitu.hot.rates.BAMM.CI<-lapply(insitu.hot.rates.BAMM.CI,function(x)quantile(x,c(0.05,0.95,0.50),na.rm=TRUE))
  
  hot.into.realm.time.bins.counts.rev.cumsum.CI<-list()
  realm.into.hot.time.bins.counts.rev.cumsum.CI<-list()
  insitu.hot.time.bins.counts.rev.cumsum.CI<-list()
  insitu.realm.time.bins.counts.rev.cumsum.CI<-list()
  dispersal.into.hot.time.bins.counts.rev.cumsum.CI<-list()
  dispersal.into.realm.time.bins.counts.rev.cumsum.CI<-list()
  dispersal.from.realm.time.bins.counts.rev.cumsum.CI<-list()
  dispersal.from.hot.time.bins.counts.rev.cumsum.CI<-list()
  
  #get cumulative events counts
  for (i in 1:length(dispersal.from.realm.time.bins.counts.rev.cumsum[[1]])){
    dispersal.from.realm.time.bins.counts.rev.cumsum.CI[[i]]<-sapply(dispersal.from.realm.time.bins.counts.rev.cumsum,"[[",i)
    dispersal.from.hot.time.bins.counts.rev.cumsum.CI[[i]]<-sapply(dispersal.from.hot.time.bins.counts.rev.cumsum,"[[",i)
    dispersal.into.realm.time.bins.counts.rev.cumsum.CI[[i]]<-sapply(dispersal.into.realm.time.bins.counts.rev.cumsum,"[[",i)
    dispersal.into.hot.time.bins.counts.rev.cumsum.CI[[i]]<-sapply(dispersal.into.hot.time.bins.counts.rev.cumsum,"[[",i)
    hot.into.realm.time.bins.counts.rev.cumsum.CI[[i]]<-sapply(hot.into.realm.time.bins.counts.rev.cumsum,"[[",i)
    realm.into.hot.time.bins.counts.rev.cumsum.CI[[i]]<-sapply(realm.into.hot.time.bins.counts.rev.cumsum,"[[",i)
    insitu.realm.time.bins.counts.rev.cumsum.CI[[i]]<-sapply(insitu.realm.time.bins.counts.rev.cumsum,"[[",i)
    insitu.hot.time.bins.counts.rev.cumsum.CI[[i]]<-sapply(insitu.hot.time.bins.counts.rev.cumsum,"[[",i)
  }
  
  
  dispersal.from.realm.time.bins.counts.rev.cumsum.CI<-lapply(dispersal.from.realm.time.bins.counts.rev.cumsum.CI,function(x)quantile(x,c(0.05,0.95,0.50)))
  dispersal.from.hot.time.bins.counts.rev.cumsum.CI<-lapply(dispersal.from.hot.time.bins.counts.rev.cumsum.CI,function(x)quantile(x,c(0.05,0.95,0.50)))
  dispersal.into.realm.time.bins.counts.rev.cumsum.CI<-lapply(dispersal.into.realm.time.bins.counts.rev.cumsum.CI,function(x)quantile(x,c(0.05,0.95,0.50)))
  dispersal.into.hot.time.bins.counts.rev.cumsum.CI<-lapply(dispersal.into.hot.time.bins.counts.rev.cumsum.CI,function(x)quantile(x,c(0.05,0.95,0.50)))
  hot.into.realm.time.bins.counts.rev.cumsum.CI<-lapply(hot.into.realm.time.bins.counts.rev.cumsum.CI,function(x)quantile(x,c(0.05,0.95,0.50)))
  realm.into.hot.time.bins.counts.rev.cumsum.CI<-lapply(realm.into.hot.time.bins.counts.rev.cumsum.CI,function(x)quantile(x,c(0.05,0.95,0.50)))
  insitu.realm.time.bins.counts.rev.cumsum.CI<-lapply(insitu.realm.time.bins.counts.rev.cumsum.CI,function(x)quantile(x,c(0.05,0.95,0.50)))
  insitu.hot.time.bins.counts.rev.cumsum.CI<-lapply(insitu.hot.time.bins.counts.rev.cumsum.CI,function(x)quantile(x,c(0.05,0.95,0.50)))
  
  rates.list<-list(insitu.realm.rates.BSM.CI,insitu.hot.rates.BSM.CI,hot.into.realm.rates.BSM.CI,realm.into.hot.rates.BSM.CI,dispersal.into.realm.rates.BSM.CI,dispersal.into.hot.rates.BSM.CI,insitu.realm.rates.strict.BAMM.CI,insitu.hot.rates.strict.BAMM.CI,insitu.realm.rates.BAMM.CI,insitu.hot.rates.BAMM.CI)
  names(rates.list)<-c("insitu.realm.rates.BSM.CI","insitu.hot.rates.BSM.CI","hot.into.realm.rates.BSM.CI","realm.into.hot.rates.BSM.CI","dispersal.into.realm.rates.BSM.CI","dispersal.into.hot.rates.BSM.CI","insitu.realm.rates.strict.BAMM.CI","insitu.hot.rates.strict.BAMM.CI","insitu.realm.rates.BAMM.CI","insitu.hot.rates.BAMM.CI")
  events.list<-list(hot.into.realm.time.bins.counts.rev.cumsum.CI,realm.into.hot.time.bins.counts.rev.cumsum.CI,insitu.hot.time.bins.counts.rev.cumsum.CI,insitu.realm.time.bins.counts.rev.cumsum.CI,dispersal.into.hot.time.bins.counts.rev.cumsum.CI,dispersal.into.realm.time.bins.counts.rev.cumsum.CI,dispersal.from.realm.time.bins.counts.rev.cumsum.CI,dispersal.from.hot.time.bins.counts.rev.cumsum.CI)
  names(events.list)<-c("hot.into.realm.time.bins.counts.rev.cumsum.CI","realm.into.hot.time.bins.counts.rev.cumsum.CI","insitu.hot.time.bins.counts.rev.cumsum.CI","insitu.realm.time.bins.counts.rev.cumsum.CI","dispersal.into.hot.time.bins.counts.rev.cumsum.CI","dispersal.into.realm.time.bins.counts.rev.cumsum.CI","dispersal.from.realm.time.bins.counts.rev.cumsum.CI","dispersal.from.hot.time.bins.counts.rev.cumsum.CI")
  output.list<-list(rates.list,events.list,realm.lineages.counts.BSM,hot.lineages.counts.BSM)
  return(output.list)
}

get_BSM_BAMM_rates_events_2my_object<-function(results.BSM,name,speciation.rate.treefile){
  if(class(speciation.rate.treefile)=='character'){
    cat('getting speciation rates','\n')
    speciation.rate.tree<-read.tree(speciation.rate.treefile)  
    table.speciation.rate<-prt(speciation.rate.tree,printflag = FALSE)
  }else if (class(speciation.rate.treefile)=='data.frame'){
    table.speciation.rate<-speciation.rate.treefile
  }
  area.codes<-c('O','A','B','C','D','E','F','G','AB','AC','AD','AE','AF','AG','BC','BD','BE','BF','BG','CD','CE','CF','CG','DE','DF','DG','EF','EG','FG')
  #get event names from BSM
  if(name=='afrotrop'){
    realm.code<-'A'
  }else if(name=='austral'){
    realm.code<-'B'
  }else if(name=='indo'){
    realm.code<-'C'
  }else if(name=='nearctic'){
    realm.code<-'D'
  }else if(name=='neotrop'){
    realm.code<-'E'
  }else if(name=='palearctic'){
    realm.code<-'F'
  }
  hot.code<-'G'
  #read clado and ana events table of each BSM replicate
  insitu.realm.rates.BSM<-list()
  insitu.hot.rates.BSM<-list()
  insitu.realm.rates.BAMM.strict<-list()
  insitu.hot.rates.BAMM.strict<-list()
  insitu.realm.rates.BAMM<-list()
  insitu.hot.rates.BAMM<-list()
  dispersal.into.realm.rates.BSM<-list()
  dispersal.into.hot.rates.BSM<-list()
  hot.into.realm.rates.BSM<-list()
  realm.into.hot.rates.BSM<-list()
  
  dispersal.into.hot.time.bins.counts.rev.cumsum<-list()
  dispersal.into.realm.time.bins.counts.rev.cumsum<-list()
  hot.into.realm.time.bins.counts.rev.cumsum<-list()
  realm.into.hot.time.bins.counts.rev.cumsum<-list()
  insitu.hot.time.bins.counts.rev.cumsum<-list()
  insitu.realm.time.bins.counts.rev.cumsum<-list()
  dispersal.from.hot.time.bins.counts.rev.cumsum<-list()
  dispersal.from.realm.time.bins.counts.rev.cumsum<-list()
  realm.into.hot.time.bins.counts.rev.cumsum<-list()
  into.into.realm.time.bins.counts.rev.cumsum<-list()
  realm.lineages.counts.BSM<-list()
  hot.lineages.counts.BSM<-list()
  cat('analysing BSMs','\n')
  for (i in 1:length(results.BSM$RES_clado_events_tables)){
    ###BioGeoBEARS
    #cat(i,'\n')
    a1<-results.BSM$RES_clado_events_tables[[i]]
    a2<-results.BSM$RES_ana_events_tables[[i]]
    #assign events to time bin
    a1$time.bin<-floor(a1$time_bp/2)*2
    a2$time.bin<-floor(a2$time_bp/2)*2
    max.age<-max(a1$time.bin)
    #get time bins of each event type
    #dispersal into hot
    dispersal.into.hot.time.bins<-sort(c(a1[a1$clado_dispersal_to==hot.code,'time.bin'],a2[a2$dispersal_to==hot.code,'time.bin']))
    dispersal.into.hot.time.bins.counts<-sapply(seq(from=0,to=max.age,by=2),function(x)length(dispersal.into.hot.time.bins[dispersal.into.hot.time.bins%in%x]))
    dispersal.into.hot.time.bins.counts.rev<-rev(dispersal.into.hot.time.bins.counts)
    dispersal.into.hot.time.bins.counts.rev.cumsum[[i]]<-cumsum(dispersal.into.hot.time.bins.counts.rev)
    #dispersal into realm
    dispersal.into.realm.time.bins<-sort(c(a1[a1$clado_dispersal_to==realm.code,'time.bin'],a2[a2$dispersal_to==realm.code,'time.bin']))
    dispersal.into.realm.time.bins.counts<-sapply(seq(from=0,to=max.age,by=2),function(x)length(dispersal.into.realm.time.bins[dispersal.into.realm.time.bins%in%x]))
    dispersal.into.realm.time.bins.counts.rev<-rev(dispersal.into.realm.time.bins.counts)
    dispersal.into.realm.time.bins.counts.rev.cumsum[[i]]<-cumsum(dispersal.into.realm.time.bins.counts.rev)
    #hot local extinctions
    hot.local.extinctions.time.bins<-sort(a2[a2$extirpation_from==hot.code,'time.bin'])
    hot.local.extinctions.time.bins.counts<-sapply(seq(from=0,to=max.age,by=2),function(x)length(hot.local.extinctions.time.bins[hot.local.extinctions.time.bins%in%x]))
    hot.local.extinctions.time.bins.counts.rev<-rev(hot.local.extinctions.time.bins.counts)
    hot.local.extinctions.time.bins.counts.rev.cumsum<-cumsum(hot.local.extinctions.time.bins.counts.rev)
    #realm local extinctions
    realm.local.extinctions.time.bins<-sort(a2[a2$extirpation_from==realm.code,'time.bin'])
    realm.local.extinctions.time.bins.counts<-sapply(seq(from=0,to=max.age,by=2),function(x)length(realm.local.extinctions.time.bins[realm.local.extinctions.time.bins%in%x]))
    realm.local.extinctions.time.bins.counts.rev<-rev(realm.local.extinctions.time.bins.counts)
    realm.local.extinctions.time.bins.counts.rev.cumsum<-cumsum(realm.local.extinctions.time.bins.counts.rev)
    #hot into realm dispersal
    hot.into.realm.time.bins<-sort(c(a2[intersect(intersect(grep(hot.code,a2$current_rangetxt),grep(realm.code,a2$current_rangetxt,invert=TRUE)),grep(realm.code,a2$new_rangetxt)),'time.bin'],a1[a1$sampled_states_AT_nodes%in%intersect(grep(hot.code,area.codes),grep(realm.code,area.codes,invert=TRUE))&a1$clado_dispersal_to==realm.code,'time.bin']))
    hot.into.realm.time.bins.counts<-sapply(seq(from=0,to=max.age,by=2),function(x)length(hot.into.realm.time.bins[hot.into.realm.time.bins%in%x]))
    hot.into.realm.time.bins.counts.rev<-rev(hot.into.realm.time.bins.counts)
    hot.into.realm.time.bins.counts.rev.cumsum[[i]]<-cumsum(hot.into.realm.time.bins.counts.rev)
    #realm into hot dispersal
    realm.into.hot.time.bins<-sort(c(a2[intersect(intersect(grep(realm.code,a2$current_rangetxt),grep(hot.code,a2$current_rangetxt,invert=TRUE)),grep(hot.code,a2$new_rangetxt)),'time.bin'],a1[a1$sampled_states_AT_nodes%in%intersect(grep(realm.code,area.codes),grep(hot.code,area.codes,invert=TRUE))&a1$clado_dispersal_to==hot.code,'time.bin']))
    realm.into.hot.time.bins.counts<-sapply(seq(from=0,to=max.age,by=2),function(x)length(realm.into.hot.time.bins[realm.into.hot.time.bins%in%x]))
    realm.into.hot.time.bins.counts.rev<-rev(realm.into.hot.time.bins.counts)
    realm.into.hot.time.bins.counts.rev.cumsum[[i]]<-cumsum(realm.into.hot.time.bins.counts.rev)
    #in situ hot cladogenesis
    insitu.hot.time.bins<-sort(a1[a1$clado_dispersal_to==''&(a1$sampled_states_AT_nodes%in%grep(hot.code,area.codes))&a1$clado_event_txt!='','time.bin'])
    insitu.hot.time.bins.counts<-sapply(seq(from=0,to=max.age,by=2),function(x)length(insitu.hot.time.bins[insitu.hot.time.bins%in%x]))
    insitu.hot.time.bins.counts.rev<-rev(insitu.hot.time.bins.counts)
    insitu.hot.time.bins.counts.rev.cumsum[[i]]<-cumsum(insitu.hot.time.bins.counts.rev)
    #in situ realm cladogenesis
    insitu.realm.time.bins<-sort(a1[a1$clado_dispersal_to==''&(a1$sampled_states_AT_nodes%in%grep(realm.code,area.codes))&a1$clado_event_txt!='','time.bin'])
    insitu.realm.time.bins.counts<-sapply(seq(from=0,to=max.age,by=2),function(x)length(insitu.realm.time.bins[insitu.realm.time.bins%in%x]))
    insitu.realm.time.bins.counts.rev<-rev(insitu.realm.time.bins.counts)
    insitu.realm.time.bins.counts.rev.cumsum[[i]]<-cumsum(insitu.realm.time.bins.counts.rev)
    #dispersal from hot
    dispersal.from.hot.time.bins<-sort(c(a2[a2$extirpation_from=='-'&(a2$current_rangetxt%in%area.codes[grep(hot.code,area.codes)]),'time.bin'],a1[(a1$clado_event_txt!='')&(a1$sampled_states_AT_nodes%in%grep(hot.code,area.codes))&(!(a1$samp_LEFT_dcorner%in%grep(hot.code,area.codes))|!(a1$samp_RIGHT_dcorner%in%grep(hot.code,area.codes))),'time.bin']))
    dispersal.from.hot.time.bins.counts<-sapply(seq(from=0,to=max.age,by=2),function(x)length(dispersal.from.hot.time.bins[dispersal.from.hot.time.bins%in%x]))
    dispersal.from.hot.time.bins.counts.rev<-rev(dispersal.from.hot.time.bins.counts)
    dispersal.from.hot.time.bins.counts.rev.cumsum[[i]]<-cumsum(dispersal.from.hot.time.bins.counts.rev)
    #dispersal from realm
    dispersal.from.realm.time.bins<-sort(c(a2[a2$extirpation_from=='-'&(a2$current_rangetxt%in%area.codes[grep(realm.code,area.codes)]),'time.bin'],a1[(a1$clado_event_txt!='')&(a1$sampled_states_AT_nodes%in%grep(realm.code,area.codes))&(!(a1$samp_LEFT_dcorner%in%grep(realm.code,area.codes))|!(a1$samp_RIGHT_dcorner%in%grep(realm.code,area.codes))),'time.bin']))
    dispersal.from.realm.time.bins.counts<-sapply(seq(from=0,to=max.age,by=2),function(x)length(dispersal.from.realm.time.bins[dispersal.from.realm.time.bins%in%x]))
    dispersal.from.realm.time.bins.counts.rev<-rev(dispersal.from.realm.time.bins.counts)
    dispersal.from.realm.time.bins.counts.rev.cumsum[[i]]<-cumsum(dispersal.from.realm.time.bins.counts.rev)
    #calculate number of lineages in hot and realm at each time bin
    realm.lineages.counts<-insitu.realm.time.bins.counts.rev.cumsum[[i]]+dispersal.into.realm.time.bins.counts.rev.cumsum[[i]]-realm.local.extinctions.time.bins.counts.rev.cumsum
    hot.lineages.counts<-insitu.hot.time.bins.counts.rev.cumsum[[i]]+dispersal.into.hot.time.bins.counts.rev.cumsum[[i]]-hot.local.extinctions.time.bins.counts.rev.cumsum
    realm.lineages.counts.BSM[[i]]<-list()
    hot.lineages.counts.BSM[[i]]<-list()
    realm.lineages.counts.BSM[[i]]<-realm.lineages.counts
    hot.lineages.counts.BSM[[i]]<-hot.lineages.counts
    #calculate in.situ.rates.BSM
    insitu.realm.rates.BSM[[i]]<-unlist(lapply(2:length(realm.lineages.counts),function(x)insitu.realm.time.bins.counts.rev[x]/realm.lineages.counts[x-1]))
    insitu.hot.rates.BSM[[i]]<-unlist(lapply(2:length(hot.lineages.counts),function(x)insitu.hot.time.bins.counts.rev[x]/hot.lineages.counts[x-1]))
    #calculate colonisation rates BSM
    dispersal.into.realm.rates.BSM[[i]]<-unlist(lapply(2:length(realm.lineages.counts),function(x)dispersal.into.realm.time.bins.counts.rev[x]/realm.lineages.counts[x-1]))
    dispersal.into.hot.rates.BSM[[i]]<-unlist(lapply(2:length(hot.lineages.counts),function(x)dispersal.into.hot.time.bins.counts.rev[x]/hot.lineages.counts[x-1]))
    #calculate dispersal rates hot->realm and viceversa BSM
    hot.into.realm.rates.BSM[[i]]<-unlist(lapply(2:length(realm.lineages.counts),function(x)hot.into.realm.time.bins.counts.rev[x]/realm.lineages.counts[x-1]))
    realm.into.hot.rates.BSM[[i]]<-unlist(lapply(2:length(hot.lineages.counts),function(x)realm.into.hot.time.bins.counts.rev[x]/hot.lineages.counts[x-1]))
    
    ###BAMM
    #assign regions to each node
    regions<-vector()
    for(a in 1:nrow(a1)){
      if(a1[a,'node.type']=='root'){
        regions[a]<-area.codes[a1[a,'sampled_states_AT_nodes']]
      }else if((a1[a,'sampled_states_AT_nodes']==a1[a,'sampled_states_AT_brbots'])&(a1[a,'anagenetic_events_txt_below_node']=='none')){
        regions[a]<-area.codes[a1[a,'sampled_states_AT_nodes']]
      }else{
        regions[a]<-NA
      }#this ignores anagenetic changes, only assigns rates to branches that do not change state
      #else{
      #   new.range<-strsplit(a1[a,'anagenetic_events_txt_below_node'],split=',')[[1]][grep('new_rangetxt',strsplit(a1[a,'anagenetic_events_txt_below_node'],split=',')[[1]])[length(grep('new_rangetxt',strsplit(a1[a,'anagenetic_events_txt_below_node'],split=',')[[1]]))]]
      #   regions[a]<-gsub(new.range,pattern='new_rangetxt:',replacement='')
      #  }
    }
    a1$speciation.rate<-table.speciation.rate$edge.length
    a1$region<-regions
    #assign events to time bin
    a1$time.bin<-floor(a1$time_bp/2)*2
    #get mean rates per region (strict = species only can occur in a single region)
    if(length(a1$region[a1$region%in%c(realm.code)])>1){
      insitu.realm.rates.BAMM.strict[[i]]<-aggregate(speciation.rate~time.bin,data=a1[a1$region%in%c(realm.code),],function(x)mean(x))
    }else if (length(a1$region[a1$region%in%c(realm.code)])==1){
      insitu.realm.rates.BAMM.strict[[i]]<-as.data.frame(cbind(a1[a1$region%in%c(realm.code),'time.bin'],a1[a1$region%in%c(realm.code),'speciation.rate']),stringsAsFactors = F)
      colnames(insitu.realm.rates.BAMM.strict[[i]])<-c('time.bin','speciation.rate')
    }else{
      insitu.realm.rates.BAMM.strict[[i]]<-as.data.frame(cbind(1,0),stringsAsFactors = F)
      colnames(insitu.realm.rates.BAMM.strict[[i]])<-c('time.bin','speciation.rate')
    }
    if(length(a1$region[a1$region%in%c(hot.code)])>1){
      insitu.hot.rates.BAMM.strict[[i]]<-aggregate(speciation.rate~time.bin,data=a1[a1$region%in%c(hot.code),],function(x)mean(x))
    }else if (length(a1$region[a1$region%in%c(hot.code)])==1){
      insitu.hot.rates.BAMM.strict[[i]]<-as.data.frame(cbind(a1[a1$region%in%c(hot.code),'time.bin'],a1[a1$region%in%c(hot.code),'speciation.rate']),stringsAsFactors = F)
      colnames(insitu.hot.rates.BAMM.strict[[i]])<-c('time.bin','speciation.rate')
    }else{
      insitu.hot.rates.BAMM.strict[[i]]<-as.data.frame(cbind(1,0),stringsAsFactors = F)
      colnames(insitu.hot.rates.BAMM.strict[[i]])<-c('time.bin','speciation.rate')
    }
    if(length(a1$region[a1$region%in%c(area.codes[grep(realm.code,area.codes)])])>1){
      insitu.realm.rates.BAMM[[i]]<-aggregate(speciation.rate~time.bin,data=a1[a1$region%in%c(area.codes[grep(realm.code,area.codes)]),],function(x)mean(x))
    }else if (length(a1$region[a1$region%in%c(area.codes[grep(realm.code,area.codes)])])==1){
      insitu.realm.rates.BAMM[[i]]<-as.data.frame(cbind(a1[a1$region%in%c(area.codes[grep(realm.code,area.codes)]),'time.bin'],a1[a1$region%in%c(area.codes[grep(realm.code,area.codes)]),'speciation.rate']),stringsAsFactors = F)
      colnames(insitu.realm.rates.BAMM[[i]])<-c('time.bin','speciation.rate')
    }else{
      insitu.realm.rates.BAMM[[i]]<-as.data.frame(cbind(1,0),stringsAsFactors = F)
      colnames(insitu.realm.rates.BAMM[[i]])<-c('time.bin','speciation.rate')
    }
    if(length(a1$region[a1$region%in%c(area.codes[grep(hot.code,area.codes)])])>1){
      insitu.hot.rates.BAMM[[i]]<-aggregate(speciation.rate~time.bin,data=a1[a1$region%in%c(area.codes[grep(hot.code,area.codes)]),],function(x)mean(x))
    }else if (length(a1$region[a1$region%in%c(area.codes[grep(hot.code,area.codes)])])==1){
      insitu.hot.rates.BAMM[[i]]<-as.data.frame(cbind(a1[a1$region%in%c(area.codes[grep(hot.code,area.codes)]),'time.bin'],a1[a1$region%in%c(area.codes[grep(hot.code,area.codes)]),'speciation.rate']),stringsAsFactors = F)
      colnames(insitu.hot.rates.BAMM[[i]])<-c('time.bin','speciation.rate')
    }else{
      insitu.hot.rates.BAMM[[i]]<-as.data.frame(cbind(1,0),stringsAsFactors = F)
      colnames(insitu.hot.rates.BAMM[[i]])<-c('time.bin','speciation.rate')
    }
    
    
  }
  #add 0s to empty bins
  insitu.realm.rates.BAMM<-lapply(insitu.realm.rates.BAMM,function(x) {missing.time.bins<-seq(from=0,to=max.age,by=2)[!c(seq(from=0,to=max.age,by=2)%in%x$time.bin)];missing.df<-as.data.frame(cbind(missing.time.bins,rep(0,length(missing.time.bins))));colnames(missing.df)<-c('time.bin','speciation.rate');new.df<-rbind(x,missing.df);new.df<-new.df[order(new.df$time.bin),];return(rev(new.df$speciation.rate))})
  insitu.hot.rates.BAMM<-lapply(insitu.hot.rates.BAMM,function(x) {missing.time.bins<-seq(from=0,to=max.age,by=2)[!c(seq(from=0,to=max.age,by=2)%in%x$time.bin)];missing.df<-as.data.frame(cbind(missing.time.bins,rep(0,length(missing.time.bins))));colnames(missing.df)<-c('time.bin','speciation.rate');new.df<-rbind(x,missing.df);new.df<-new.df[order(new.df$time.bin),];return(rev(new.df$speciation.rate))})
  #add 0s to empty bins
  insitu.realm.rates.BAMM.strict<-lapply(insitu.realm.rates.BAMM.strict,function(x) {missing.time.bins<-seq(from=0,to=max.age,by=2)[!c(seq(from=0,to=max.age,by=2)%in%x$time.bin)];missing.df<-as.data.frame(cbind(missing.time.bins,rep(0,length(missing.time.bins))));colnames(missing.df)<-c('time.bin','speciation.rate');new.df<-rbind(x,missing.df);new.df<-new.df[order(new.df$time.bin),];return(rev(new.df$speciation.rate))})
  insitu.hot.rates.BAMM.strict<-lapply(insitu.hot.rates.BAMM.strict,function(x) {missing.time.bins<-seq(from=0,to=max.age,by=2)[!c(seq(from=0,to=max.age,by=2)%in%x$time.bin)];missing.df<-as.data.frame(cbind(missing.time.bins,rep(0,length(missing.time.bins))));colnames(missing.df)<-c('time.bin','speciation.rate');new.df<-rbind(x,missing.df);new.df<-new.df[order(new.df$time.bin),];return(rev(new.df$speciation.rate))})
  
  #get the BSM rates CIs to draw on plot
  insitu.realm.rates.BSM.CI<-list()
  insitu.hot.rates.BSM.CI<-list()
  for (i in 1:length(insitu.realm.rates.BSM[[1]])){
    insitu.realm.rates.BSM.CI[[i]]<-sapply(insitu.realm.rates.BSM,"[[",i)
    insitu.hot.rates.BSM.CI[[i]]<-sapply(insitu.hot.rates.BSM,"[[",i)
  }
  insitu.realm.rates.BSM.CI<-lapply(insitu.realm.rates.BSM.CI,function(x)quantile(x,c(0.05,0.95,0.50),na.rm=TRUE))
  insitu.hot.rates.BSM.CI<-lapply(insitu.hot.rates.BSM.CI,function(x)quantile(x,c(0.05,0.95,0.50),na.rm=TRUE))
  #get the BSM dispersal rates CIs to draw on plot
  hot.into.realm.rates.BSM.CI<-list()
  realm.into.hot.rates.BSM.CI<-list()
  for (i in 1:length(hot.into.realm.rates.BSM[[1]])){
    hot.into.realm.rates.BSM.CI[[i]]<-sapply(hot.into.realm.rates.BSM,"[[",i)
    realm.into.hot.rates.BSM.CI[[i]]<-sapply(realm.into.hot.rates.BSM,"[[",i)
  }
  hot.into.realm.rates.BSM.CI<-lapply(hot.into.realm.rates.BSM.CI,function(x)quantile(x,c(0.05,0.95,0.50),na.rm=TRUE))
  realm.into.hot.rates.BSM.CI<-lapply(realm.into.hot.rates.BSM.CI,function(x)quantile(x,c(0.05,0.95,0.50),na.rm=TRUE))
  dispersal.into.realm.rates.BSM.CI<-list()
  dispersal.into.hot.rates.BSM.CI<-list()
  for (i in 1:length(dispersal.into.realm.rates.BSM[[1]])){
    dispersal.into.realm.rates.BSM.CI[[i]]<-sapply(dispersal.into.realm.rates.BSM,"[[",i)
    dispersal.into.hot.rates.BSM.CI[[i]]<-sapply(dispersal.into.hot.rates.BSM,"[[",i)
  }
  dispersal.into.realm.rates.BSM.CI<-lapply(dispersal.into.realm.rates.BSM.CI,function(x)quantile(x,c(0.05,0.95,0.50),na.rm=TRUE))
  dispersal.into.hot.rates.BSM.CI<-lapply(dispersal.into.hot.rates.BSM.CI,function(x)quantile(x,c(0.05,0.95,0.50),na.rm=TRUE))
  #get the BAMM strict rates CIs to draw on plot
  
  insitu.realm.rates.strict.BAMM.CI<-list()
  insitu.hot.rates.strict.BAMM.CI<-list()
  for (i in 1:length(insitu.realm.rates.BAMM.strict[[1]])){
    insitu.realm.rates.strict.BAMM.CI[[i]]<-sapply(insitu.realm.rates.BAMM.strict,"[[",i)
    insitu.hot.rates.strict.BAMM.CI[[i]]<-sapply(insitu.hot.rates.BAMM.strict,"[[",i)
  }
  insitu.realm.rates.strict.BAMM.CI<-lapply(insitu.realm.rates.strict.BAMM.CI,function(x)quantile(x,c(0.05,0.95,0.50),na.rm=TRUE))
  insitu.hot.rates.strict.BAMM.CI<-lapply(insitu.hot.rates.strict.BAMM.CI,function(x)quantile(x,c(0.05,0.95,0.50),na.rm=TRUE))
  #get the BAMM rates CIs to draw on plot
  insitu.realm.rates.BAMM.CI<-list()
  insitu.hot.rates.BAMM.CI<-list()
  for (i in 1:length(insitu.realm.rates.BAMM[[1]])){
    insitu.realm.rates.BAMM.CI[[i]]<-sapply(insitu.realm.rates.BAMM,"[[",i)
    insitu.hot.rates.BAMM.CI[[i]]<-sapply(insitu.hot.rates.BAMM,"[[",i)
  }
  insitu.realm.rates.BAMM.CI<-lapply(insitu.realm.rates.BAMM.CI,function(x)quantile(x,c(0.05,0.95,0.50),na.rm=TRUE))
  insitu.hot.rates.BAMM.CI<-lapply(insitu.hot.rates.BAMM.CI,function(x)quantile(x,c(0.05,0.95,0.50),na.rm=TRUE))
  
  hot.into.realm.time.bins.counts.rev.cumsum.CI<-list()
  realm.into.hot.time.bins.counts.rev.cumsum.CI<-list()
  insitu.hot.time.bins.counts.rev.cumsum.CI<-list()
  insitu.realm.time.bins.counts.rev.cumsum.CI<-list()
  dispersal.into.hot.time.bins.counts.rev.cumsum.CI<-list()
  dispersal.into.realm.time.bins.counts.rev.cumsum.CI<-list()
  dispersal.from.realm.time.bins.counts.rev.cumsum.CI<-list()
  dispersal.from.hot.time.bins.counts.rev.cumsum.CI<-list()
  
  #get cumulative events counts
  for (i in 1:length(dispersal.from.realm.time.bins.counts.rev.cumsum[[1]])){
    dispersal.from.realm.time.bins.counts.rev.cumsum.CI[[i]]<-sapply(dispersal.from.realm.time.bins.counts.rev.cumsum,"[[",i)
    dispersal.from.hot.time.bins.counts.rev.cumsum.CI[[i]]<-sapply(dispersal.from.hot.time.bins.counts.rev.cumsum,"[[",i)
    dispersal.into.realm.time.bins.counts.rev.cumsum.CI[[i]]<-sapply(dispersal.into.realm.time.bins.counts.rev.cumsum,"[[",i)
    dispersal.into.hot.time.bins.counts.rev.cumsum.CI[[i]]<-sapply(dispersal.into.hot.time.bins.counts.rev.cumsum,"[[",i)
    hot.into.realm.time.bins.counts.rev.cumsum.CI[[i]]<-sapply(hot.into.realm.time.bins.counts.rev.cumsum,"[[",i)
    realm.into.hot.time.bins.counts.rev.cumsum.CI[[i]]<-sapply(realm.into.hot.time.bins.counts.rev.cumsum,"[[",i)
    insitu.realm.time.bins.counts.rev.cumsum.CI[[i]]<-sapply(insitu.realm.time.bins.counts.rev.cumsum,"[[",i)
    insitu.hot.time.bins.counts.rev.cumsum.CI[[i]]<-sapply(insitu.hot.time.bins.counts.rev.cumsum,"[[",i)
  }
  
  
  dispersal.from.realm.time.bins.counts.rev.cumsum.CI<-lapply(dispersal.from.realm.time.bins.counts.rev.cumsum.CI,function(x)quantile(x,c(0.05,0.95,0.50)))
  dispersal.from.hot.time.bins.counts.rev.cumsum.CI<-lapply(dispersal.from.hot.time.bins.counts.rev.cumsum.CI,function(x)quantile(x,c(0.05,0.95,0.50)))
  dispersal.into.realm.time.bins.counts.rev.cumsum.CI<-lapply(dispersal.into.realm.time.bins.counts.rev.cumsum.CI,function(x)quantile(x,c(0.05,0.95,0.50)))
  dispersal.into.hot.time.bins.counts.rev.cumsum.CI<-lapply(dispersal.into.hot.time.bins.counts.rev.cumsum.CI,function(x)quantile(x,c(0.05,0.95,0.50)))
  hot.into.realm.time.bins.counts.rev.cumsum.CI<-lapply(hot.into.realm.time.bins.counts.rev.cumsum.CI,function(x)quantile(x,c(0.05,0.95,0.50)))
  realm.into.hot.time.bins.counts.rev.cumsum.CI<-lapply(realm.into.hot.time.bins.counts.rev.cumsum.CI,function(x)quantile(x,c(0.05,0.95,0.50)))
  insitu.realm.time.bins.counts.rev.cumsum.CI<-lapply(insitu.realm.time.bins.counts.rev.cumsum.CI,function(x)quantile(x,c(0.05,0.95,0.50)))
  insitu.hot.time.bins.counts.rev.cumsum.CI<-lapply(insitu.hot.time.bins.counts.rev.cumsum.CI,function(x)quantile(x,c(0.05,0.95,0.50)))
  
  rates.list<-list(insitu.realm.rates.BSM.CI,insitu.hot.rates.BSM.CI,hot.into.realm.rates.BSM.CI,realm.into.hot.rates.BSM.CI,dispersal.into.realm.rates.BSM.CI,dispersal.into.hot.rates.BSM.CI,insitu.realm.rates.strict.BAMM.CI,insitu.hot.rates.strict.BAMM.CI,insitu.realm.rates.BAMM.CI,insitu.hot.rates.BAMM.CI)
  names(rates.list)<-c("insitu.realm.rates.BSM.CI","insitu.hot.rates.BSM.CI","hot.into.realm.rates.BSM.CI","realm.into.hot.rates.BSM.CI","dispersal.into.realm.rates.BSM.CI","dispersal.into.hot.rates.BSM.CI","insitu.realm.rates.strict.BAMM.CI","insitu.hot.rates.strict.BAMM.CI","insitu.realm.rates.BAMM.CI","insitu.hot.rates.BAMM.CI")
  events.list<-list(hot.into.realm.time.bins.counts.rev.cumsum.CI,realm.into.hot.time.bins.counts.rev.cumsum.CI,insitu.hot.time.bins.counts.rev.cumsum.CI,insitu.realm.time.bins.counts.rev.cumsum.CI,dispersal.into.hot.time.bins.counts.rev.cumsum.CI,dispersal.into.realm.time.bins.counts.rev.cumsum.CI,dispersal.from.realm.time.bins.counts.rev.cumsum.CI,dispersal.from.hot.time.bins.counts.rev.cumsum.CI)
  names(events.list)<-c("hot.into.realm.time.bins.counts.rev.cumsum.CI","realm.into.hot.time.bins.counts.rev.cumsum.CI","insitu.hot.time.bins.counts.rev.cumsum.CI","insitu.realm.time.bins.counts.rev.cumsum.CI","dispersal.into.hot.time.bins.counts.rev.cumsum.CI","dispersal.into.realm.time.bins.counts.rev.cumsum.CI","dispersal.from.realm.time.bins.counts.rev.cumsum.CI","dispersal.from.hot.time.bins.counts.rev.cumsum.CI")
  output.list<-list(rates.list,events.list,realm.lineages.counts.BSM,hot.lineages.counts.BSM)
  return(output.list)
}

get_BSM_BAMM_fullrates_events_2my_object_old<-function(results.BSM,name,speciation.rate.treefile){
  if(class(speciation.rate.treefile)=='character'){
    cat('getting speciation rates','\n')
    speciation.rate.tree<-read.tree(speciation.rate.treefile)  
    table.speciation.rate<-prt(speciation.rate.tree,printflag = FALSE)
  }else if (class(speciation.rate.treefile)=='data.frame'){
    table.speciation.rate<-speciation.rate.treefile
  }
  area.codes<-c('O','A','B','C','D','E','F','G','AB','AC','AD','AE','AF','AG','BC','BD','BE','BF','BG','CD','CE','CF','CG','DE','DF','DG','EF','EG','FG')
  #get event names from BSM
  if(name=='afrotrop'){
    realm.code<-'A'
  }else if(name=='austral'){
    realm.code<-'B'
  }else if(name=='indo'){
    realm.code<-'C'
  }else if(name=='nearctic'){
    realm.code<-'D'
  }else if(name=='neotrop'){
    realm.code<-'E'
  }else if(name=='palearctic'){
    realm.code<-'F'
  }
  hot.code<-'G'
  #read clado and ana events table of each BSM replicate
  insitu.realm.rates.BSM<-list()
  insitu.hot.rates.BSM<-list()
  insitu.realm.rates.BAMM.strict<-list()
  insitu.hot.rates.BAMM.strict<-list()
  insitu.realm.rates.BAMM<-list()
  insitu.hot.rates.BAMM<-list()
  dispersal.into.realm.rates.BSM<-list()
  dispersal.into.hot.rates.BSM<-list()
  hot.into.realm.rates.BSM<-list()
  realm.into.hot.rates.BSM<-list()
  
  dispersal.into.hot.time.bins.counts.rev.cumsum<-list()
  dispersal.into.realm.time.bins.counts.rev.cumsum<-list()
  hot.into.realm.time.bins.counts.rev.cumsum<-list()
  realm.into.hot.time.bins.counts.rev.cumsum<-list()
  insitu.hot.time.bins.counts.rev.cumsum<-list()
  insitu.realm.time.bins.counts.rev.cumsum<-list()
  dispersal.from.hot.time.bins.counts.rev.cumsum<-list()
  dispersal.from.realm.time.bins.counts.rev.cumsum<-list()
  realm.into.hot.time.bins.counts.rev.cumsum<-list()
  into.into.realm.time.bins.counts.rev.cumsum<-list()
  realm.lineages.counts.BSM<-list()
  hot.lineages.counts.BSM<-list()
  cat('analysing BSMs','\n')
  for (i in 1:length(results.BSM$RES_clado_events_tables)){
    ###BioGeoBEARS
    #cat(i,'\n')
    a1<-results.BSM$RES_clado_events_tables[[i]]
    a2<-results.BSM$RES_ana_events_tables[[i]]
    #assign events to time bin
    a1$time.bin<-floor(a1$time_bp/2)*2
    a2$time.bin<-floor(a2$time_bp/2)*2
    max.age<-max(a1$time.bin)
    #get time bins of each event type
    #dispersal into hot
    dispersal.into.hot.time.bins<-sort(c(a1[a1$clado_dispersal_to==hot.code,'time.bin'],a2[a2$dispersal_to==hot.code,'time.bin']))
    dispersal.into.hot.time.bins.counts<-sapply(seq(from=0,to=max.age,by=2),function(x)length(dispersal.into.hot.time.bins[dispersal.into.hot.time.bins%in%x]))
    dispersal.into.hot.time.bins.counts.rev<-rev(dispersal.into.hot.time.bins.counts)
    dispersal.into.hot.time.bins.counts.rev.cumsum[[i]]<-cumsum(dispersal.into.hot.time.bins.counts.rev)
    #dispersal into realm
    dispersal.into.realm.time.bins<-sort(c(a1[a1$clado_dispersal_to==realm.code,'time.bin'],a2[a2$dispersal_to==realm.code,'time.bin']))
    dispersal.into.realm.time.bins.counts<-sapply(seq(from=0,to=max.age,by=2),function(x)length(dispersal.into.realm.time.bins[dispersal.into.realm.time.bins%in%x]))
    dispersal.into.realm.time.bins.counts.rev<-rev(dispersal.into.realm.time.bins.counts)
    dispersal.into.realm.time.bins.counts.rev.cumsum[[i]]<-cumsum(dispersal.into.realm.time.bins.counts.rev)
    #hot local extinctions
    hot.local.extinctions.time.bins<-sort(a2[a2$extirpation_from==hot.code,'time.bin'])
    hot.local.extinctions.time.bins.counts<-sapply(seq(from=0,to=max.age,by=2),function(x)length(hot.local.extinctions.time.bins[hot.local.extinctions.time.bins%in%x]))
    hot.local.extinctions.time.bins.counts.rev<-rev(hot.local.extinctions.time.bins.counts)
    hot.local.extinctions.time.bins.counts.rev.cumsum<-cumsum(hot.local.extinctions.time.bins.counts.rev)
    #realm local extinctions
    realm.local.extinctions.time.bins<-sort(a2[a2$extirpation_from==realm.code,'time.bin'])
    realm.local.extinctions.time.bins.counts<-sapply(seq(from=0,to=max.age,by=2),function(x)length(realm.local.extinctions.time.bins[realm.local.extinctions.time.bins%in%x]))
    realm.local.extinctions.time.bins.counts.rev<-rev(realm.local.extinctions.time.bins.counts)
    realm.local.extinctions.time.bins.counts.rev.cumsum<-cumsum(realm.local.extinctions.time.bins.counts.rev)
    #hot into realm dispersal
    hot.into.realm.time.bins<-sort(c(a2[intersect(intersect(grep(hot.code,a2$current_rangetxt),grep(realm.code,a2$current_rangetxt,invert=TRUE)),grep(realm.code,a2$new_rangetxt)),'time.bin'],a1[a1$sampled_states_AT_nodes%in%intersect(grep(hot.code,area.codes),grep(realm.code,area.codes,invert=TRUE))&a1$clado_dispersal_to==realm.code,'time.bin']))
    hot.into.realm.time.bins.counts<-sapply(seq(from=0,to=max.age,by=2),function(x)length(hot.into.realm.time.bins[hot.into.realm.time.bins%in%x]))
    hot.into.realm.time.bins.counts.rev<-rev(hot.into.realm.time.bins.counts)
    hot.into.realm.time.bins.counts.rev.cumsum[[i]]<-cumsum(hot.into.realm.time.bins.counts.rev)
    #realm into hot dispersal
    realm.into.hot.time.bins<-sort(c(a2[intersect(intersect(grep(realm.code,a2$current_rangetxt),grep(hot.code,a2$current_rangetxt,invert=TRUE)),grep(hot.code,a2$new_rangetxt)),'time.bin'],a1[a1$sampled_states_AT_nodes%in%intersect(grep(realm.code,area.codes),grep(hot.code,area.codes,invert=TRUE))&a1$clado_dispersal_to==hot.code,'time.bin']))
    realm.into.hot.time.bins.counts<-sapply(seq(from=0,to=max.age,by=2),function(x)length(realm.into.hot.time.bins[realm.into.hot.time.bins%in%x]))
    realm.into.hot.time.bins.counts.rev<-rev(realm.into.hot.time.bins.counts)
    realm.into.hot.time.bins.counts.rev.cumsum[[i]]<-cumsum(realm.into.hot.time.bins.counts.rev)
    #in situ hot cladogenesis
    insitu.hot.time.bins<-sort(a1[a1$clado_dispersal_to==''&(a1$sampled_states_AT_nodes%in%grep(hot.code,area.codes))&a1$clado_event_txt!='','time.bin'])
    insitu.hot.time.bins.counts<-sapply(seq(from=0,to=max.age,by=2),function(x)length(insitu.hot.time.bins[insitu.hot.time.bins%in%x]))
    insitu.hot.time.bins.counts.rev<-rev(insitu.hot.time.bins.counts)
    insitu.hot.time.bins.counts.rev.cumsum[[i]]<-cumsum(insitu.hot.time.bins.counts.rev)
    #in situ realm cladogenesis
    insitu.realm.time.bins<-sort(a1[a1$clado_dispersal_to==''&(a1$sampled_states_AT_nodes%in%grep(realm.code,area.codes))&a1$clado_event_txt!='','time.bin'])
    insitu.realm.time.bins.counts<-sapply(seq(from=0,to=max.age,by=2),function(x)length(insitu.realm.time.bins[insitu.realm.time.bins%in%x]))
    insitu.realm.time.bins.counts.rev<-rev(insitu.realm.time.bins.counts)
    insitu.realm.time.bins.counts.rev.cumsum[[i]]<-cumsum(insitu.realm.time.bins.counts.rev)
    #dispersal from hot
    dispersal.from.hot.time.bins<-sort(c(a2[a2$extirpation_from=='-'&(a2$current_rangetxt%in%area.codes[grep(hot.code,area.codes)]),'time.bin'],a1[(a1$clado_event_txt!='')&(a1$sampled_states_AT_nodes%in%grep(hot.code,area.codes))&(!(a1$samp_LEFT_dcorner%in%grep(hot.code,area.codes))|!(a1$samp_RIGHT_dcorner%in%grep(hot.code,area.codes))),'time.bin']))
    dispersal.from.hot.time.bins.counts<-sapply(seq(from=0,to=max.age,by=2),function(x)length(dispersal.from.hot.time.bins[dispersal.from.hot.time.bins%in%x]))
    dispersal.from.hot.time.bins.counts.rev<-rev(dispersal.from.hot.time.bins.counts)
    dispersal.from.hot.time.bins.counts.rev.cumsum[[i]]<-cumsum(dispersal.from.hot.time.bins.counts.rev)
    #dispersal from realm
    dispersal.from.realm.time.bins<-sort(c(a2[a2$extirpation_from=='-'&(a2$current_rangetxt%in%area.codes[grep(realm.code,area.codes)]),'time.bin'],a1[(a1$clado_event_txt!='')&(a1$sampled_states_AT_nodes%in%grep(realm.code,area.codes))&(!(a1$samp_LEFT_dcorner%in%grep(realm.code,area.codes))|!(a1$samp_RIGHT_dcorner%in%grep(realm.code,area.codes))),'time.bin']))
    dispersal.from.realm.time.bins.counts<-sapply(seq(from=0,to=max.age,by=2),function(x)length(dispersal.from.realm.time.bins[dispersal.from.realm.time.bins%in%x]))
    dispersal.from.realm.time.bins.counts.rev<-rev(dispersal.from.realm.time.bins.counts)
    dispersal.from.realm.time.bins.counts.rev.cumsum[[i]]<-cumsum(dispersal.from.realm.time.bins.counts.rev)
    #calculate number of lineages in hot and realm at each time bin
    realm.lineages.counts<-insitu.realm.time.bins.counts.rev.cumsum[[i]]+dispersal.into.realm.time.bins.counts.rev.cumsum[[i]]-realm.local.extinctions.time.bins.counts.rev.cumsum
    hot.lineages.counts<-insitu.hot.time.bins.counts.rev.cumsum[[i]]+dispersal.into.hot.time.bins.counts.rev.cumsum[[i]]-hot.local.extinctions.time.bins.counts.rev.cumsum
    realm.lineages.counts.BSM[[i]]<-list()
    hot.lineages.counts.BSM[[i]]<-list()
    realm.lineages.counts.BSM[[i]]<-realm.lineages.counts
    hot.lineages.counts.BSM[[i]]<-hot.lineages.counts
    #calculate in.situ.rates.BSM
    insitu.realm.rates.BSM[[i]]<-unlist(lapply(2:length(realm.lineages.counts),function(x)insitu.realm.time.bins.counts.rev[x]/realm.lineages.counts[x-1]))
    insitu.hot.rates.BSM[[i]]<-unlist(lapply(2:length(hot.lineages.counts),function(x)insitu.hot.time.bins.counts.rev[x]/hot.lineages.counts[x-1]))
    #calculate colonisation rates BSM
    dispersal.into.realm.rates.BSM[[i]]<-unlist(lapply(2:length(realm.lineages.counts),function(x)dispersal.into.realm.time.bins.counts.rev[x]/realm.lineages.counts[x-1]))
    dispersal.into.hot.rates.BSM[[i]]<-unlist(lapply(2:length(hot.lineages.counts),function(x)dispersal.into.hot.time.bins.counts.rev[x]/hot.lineages.counts[x-1]))
    #calculate dispersal rates hot->realm and viceversa BSM
    hot.into.realm.rates.BSM[[i]]<-unlist(lapply(2:length(realm.lineages.counts),function(x)hot.into.realm.time.bins.counts.rev[x]/realm.lineages.counts[x-1]))
    realm.into.hot.rates.BSM[[i]]<-unlist(lapply(2:length(hot.lineages.counts),function(x)realm.into.hot.time.bins.counts.rev[x]/hot.lineages.counts[x-1]))
    
    ###BAMM
    #assign regions to each node
    regions<-vector()
    for(a in 1:nrow(a1)){
      if(a1[a,'node.type']=='root'){
        regions[a]<-area.codes[a1[a,'sampled_states_AT_nodes']]
      }else if((a1[a,'sampled_states_AT_nodes']==a1[a,'sampled_states_AT_brbots'])&(a1[a,'anagenetic_events_txt_below_node']=='none')){
        regions[a]<-area.codes[a1[a,'sampled_states_AT_nodes']]
      }else{
        regions[a]<-NA
      }#this ignores anagenetic changes, only assigns rates to branches that do not change state
      #else{
      #   new.range<-strsplit(a1[a,'anagenetic_events_txt_below_node'],split=',')[[1]][grep('new_rangetxt',strsplit(a1[a,'anagenetic_events_txt_below_node'],split=',')[[1]])[length(grep('new_rangetxt',strsplit(a1[a,'anagenetic_events_txt_below_node'],split=',')[[1]]))]]
      #   regions[a]<-gsub(new.range,pattern='new_rangetxt:',replacement='')
      #  }
    }
    a1$speciation.rate<-table.speciation.rate$edge.length
    a1$region<-regions
    #assign events to time bin
    a1$time.bin<-floor(a1$time_bp/2)*2
    #get mean rates per region (strict = species only can occur in a single region)
    if(length(a1$region[a1$region%in%c(realm.code)])>1){
      insitu.realm.rates.BAMM.strict[[i]]<-aggregate(speciation.rate~time.bin,data=a1[a1$region%in%c(realm.code),],function(x)mean(x))
    }else if (length(a1$region[a1$region%in%c(realm.code)])==1){
      insitu.realm.rates.BAMM.strict[[i]]<-as.data.frame(cbind(a1[a1$region%in%c(realm.code),'time.bin'],a1[a1$region%in%c(realm.code),'speciation.rate']),stringsAsFactors = F)
      colnames(insitu.realm.rates.BAMM.strict[[i]])<-c('time.bin','speciation.rate')
    }else{
      insitu.realm.rates.BAMM.strict[[i]]<-as.data.frame(cbind(1,0),stringsAsFactors = F)
      colnames(insitu.realm.rates.BAMM.strict[[i]])<-c('time.bin','speciation.rate')
    }
    if(length(a1$region[a1$region%in%c(hot.code)])>1){
      insitu.hot.rates.BAMM.strict[[i]]<-aggregate(speciation.rate~time.bin,data=a1[a1$region%in%c(hot.code),],function(x)mean(x))
    }else if (length(a1$region[a1$region%in%c(hot.code)])==1){
      insitu.hot.rates.BAMM.strict[[i]]<-as.data.frame(cbind(a1[a1$region%in%c(hot.code),'time.bin'],a1[a1$region%in%c(hot.code),'speciation.rate']),stringsAsFactors = F)
      colnames(insitu.hot.rates.BAMM.strict[[i]])<-c('time.bin','speciation.rate')
    }else{
      insitu.hot.rates.BAMM.strict[[i]]<-as.data.frame(cbind(1,0),stringsAsFactors = F)
      colnames(insitu.hot.rates.BAMM.strict[[i]])<-c('time.bin','speciation.rate')
    }
    if(length(a1$region[a1$region%in%c(area.codes[grep(realm.code,area.codes)])])>1){
      insitu.realm.rates.BAMM[[i]]<-aggregate(speciation.rate~time.bin,data=a1[a1$region%in%c(area.codes[grep(realm.code,area.codes)]),],function(x)mean(x))
    }else if (length(a1$region[a1$region%in%c(area.codes[grep(realm.code,area.codes)])])==1){
      insitu.realm.rates.BAMM[[i]]<-as.data.frame(cbind(a1[a1$region%in%c(area.codes[grep(realm.code,area.codes)]),'time.bin'],a1[a1$region%in%c(area.codes[grep(realm.code,area.codes)]),'speciation.rate']),stringsAsFactors = F)
      colnames(insitu.realm.rates.BAMM[[i]])<-c('time.bin','speciation.rate')
    }else{
      insitu.realm.rates.BAMM[[i]]<-as.data.frame(cbind(1,0),stringsAsFactors = F)
      colnames(insitu.realm.rates.BAMM[[i]])<-c('time.bin','speciation.rate')
    }
    if(length(a1$region[a1$region%in%c(area.codes[grep(hot.code,area.codes)])])>1){
      insitu.hot.rates.BAMM[[i]]<-aggregate(speciation.rate~time.bin,data=a1[a1$region%in%c(area.codes[grep(hot.code,area.codes)]),],function(x)mean(x))
    }else if (length(a1$region[a1$region%in%c(area.codes[grep(hot.code,area.codes)])])==1){
      insitu.hot.rates.BAMM[[i]]<-as.data.frame(cbind(a1[a1$region%in%c(area.codes[grep(hot.code,area.codes)]),'time.bin'],a1[a1$region%in%c(area.codes[grep(hot.code,area.codes)]),'speciation.rate']),stringsAsFactors = F)
      colnames(insitu.hot.rates.BAMM[[i]])<-c('time.bin','speciation.rate')
    }else{
      insitu.hot.rates.BAMM[[i]]<-as.data.frame(cbind(1,0),stringsAsFactors = F)
      colnames(insitu.hot.rates.BAMM[[i]])<-c('time.bin','speciation.rate')
    }
    
    
  }
  #add 0s to empty bins
  insitu.realm.rates.BAMM<-lapply(insitu.realm.rates.BAMM,function(x) {missing.time.bins<-seq(from=0,to=max.age,by=2)[!c(seq(from=0,to=max.age,by=2)%in%x$time.bin)];missing.df<-as.data.frame(cbind(missing.time.bins,rep(0,length(missing.time.bins))));colnames(missing.df)<-c('time.bin','speciation.rate');new.df<-rbind(x,missing.df);new.df<-new.df[order(new.df$time.bin),];return(rev(new.df$speciation.rate))})
  insitu.hot.rates.BAMM<-lapply(insitu.hot.rates.BAMM,function(x) {missing.time.bins<-seq(from=0,to=max.age,by=2)[!c(seq(from=0,to=max.age,by=2)%in%x$time.bin)];missing.df<-as.data.frame(cbind(missing.time.bins,rep(0,length(missing.time.bins))));colnames(missing.df)<-c('time.bin','speciation.rate');new.df<-rbind(x,missing.df);new.df<-new.df[order(new.df$time.bin),];return(rev(new.df$speciation.rate))})
  #add 0s to empty bins
  insitu.realm.rates.BAMM.strict<-lapply(insitu.realm.rates.BAMM.strict,function(x) {missing.time.bins<-seq(from=0,to=max.age,by=2)[!c(seq(from=0,to=max.age,by=2)%in%x$time.bin)];missing.df<-as.data.frame(cbind(missing.time.bins,rep(0,length(missing.time.bins))));colnames(missing.df)<-c('time.bin','speciation.rate');new.df<-rbind(x,missing.df);new.df<-new.df[order(new.df$time.bin),];return(rev(new.df$speciation.rate))})
  insitu.hot.rates.BAMM.strict<-lapply(insitu.hot.rates.BAMM.strict,function(x) {missing.time.bins<-seq(from=0,to=max.age,by=2)[!c(seq(from=0,to=max.age,by=2)%in%x$time.bin)];missing.df<-as.data.frame(cbind(missing.time.bins,rep(0,length(missing.time.bins))));colnames(missing.df)<-c('time.bin','speciation.rate');new.df<-rbind(x,missing.df);new.df<-new.df[order(new.df$time.bin),];return(rev(new.df$speciation.rate))})
  
  #get the BSM rates CIs to draw on plot
  insitu.realm.rates.BSM.CI<-list()
  insitu.hot.rates.BSM.CI<-list()
  for (i in 1:length(insitu.realm.rates.BSM[[1]])){
    insitu.realm.rates.BSM.CI[[i]]<-sapply(insitu.realm.rates.BSM,"[[",i)
    insitu.hot.rates.BSM.CI[[i]]<-sapply(insitu.hot.rates.BSM,"[[",i)
  }
  #insitu.realm.rates.BSM.CI<-lapply(insitu.realm.rates.BSM.CI,function(x)quantile(x,c(0.05,0.95,0.50),na.rm=TRUE))
  #insitu.hot.rates.BSM.CI<-lapply(insitu.hot.rates.BSM.CI,function(x)quantile(x,c(0.05,0.95,0.50),na.rm=TRUE))
  #get the BSM dispersal rates CIs to draw on plot
  hot.into.realm.rates.BSM.CI<-list()
  realm.into.hot.rates.BSM.CI<-list()
  for (i in 1:length(hot.into.realm.rates.BSM[[1]])){
    hot.into.realm.rates.BSM.CI[[i]]<-sapply(hot.into.realm.rates.BSM,"[[",i)
    realm.into.hot.rates.BSM.CI[[i]]<-sapply(realm.into.hot.rates.BSM,"[[",i)
  }
  #hot.into.realm.rates.BSM.CI<-lapply(hot.into.realm.rates.BSM.CI,function(x)quantile(x,c(0.05,0.95,0.50),na.rm=TRUE))
  #realm.into.hot.rates.BSM.CI<-lapply(realm.into.hot.rates.BSM.CI,function(x)quantile(x,c(0.05,0.95,0.50),na.rm=TRUE))
  dispersal.into.realm.rates.BSM.CI<-list()
  dispersal.into.hot.rates.BSM.CI<-list()
  for (i in 1:length(dispersal.into.realm.rates.BSM[[1]])){
    dispersal.into.realm.rates.BSM.CI[[i]]<-sapply(dispersal.into.realm.rates.BSM,"[[",i)
    dispersal.into.hot.rates.BSM.CI[[i]]<-sapply(dispersal.into.hot.rates.BSM,"[[",i)
  }
  #dispersal.into.realm.rates.BSM.CI<-lapply(dispersal.into.realm.rates.BSM.CI,function(x)quantile(x,c(0.05,0.95,0.50),na.rm=TRUE))
  #dispersal.into.hot.rates.BSM.CI<-lapply(dispersal.into.hot.rates.BSM.CI,function(x)quantile(x,c(0.05,0.95,0.50),na.rm=TRUE))
  #get the BAMM strict rates CIs to draw on plot
  
  insitu.realm.rates.strict.BAMM.CI<-list()
  insitu.hot.rates.strict.BAMM.CI<-list()
  for (i in 1:length(insitu.realm.rates.BAMM.strict[[1]])){
    insitu.realm.rates.strict.BAMM.CI[[i]]<-sapply(insitu.realm.rates.BAMM.strict,"[[",i)
    insitu.hot.rates.strict.BAMM.CI[[i]]<-sapply(insitu.hot.rates.BAMM.strict,"[[",i)
  }
  #insitu.realm.rates.strict.BAMM.CI<-lapply(insitu.realm.rates.strict.BAMM.CI,function(x)quantile(x,c(0.05,0.95,0.50),na.rm=TRUE))
  #insitu.hot.rates.strict.BAMM.CI<-lapply(insitu.hot.rates.strict.BAMM.CI,function(x)quantile(x,c(0.05,0.95,0.50),na.rm=TRUE))
  #get the BAMM rates CIs to draw on plot
  insitu.realm.rates.BAMM.CI<-list()
  insitu.hot.rates.BAMM.CI<-list()
  for (i in 1:length(insitu.realm.rates.BAMM[[1]])){
    insitu.realm.rates.BAMM.CI[[i]]<-sapply(insitu.realm.rates.BAMM,"[[",i)
    insitu.hot.rates.BAMM.CI[[i]]<-sapply(insitu.hot.rates.BAMM,"[[",i)
  }
  #insitu.realm.rates.BAMM.CI<-lapply(insitu.realm.rates.BAMM.CI,function(x)quantile(x,c(0.05,0.95,0.50),na.rm=TRUE))
  #insitu.hot.rates.BAMM.CI<-lapply(insitu.hot.rates.BAMM.CI,function(x)quantile(x,c(0.05,0.95,0.50),na.rm=TRUE))
  
  hot.into.realm.time.bins.counts.rev.cumsum.CI<-list()
  realm.into.hot.time.bins.counts.rev.cumsum.CI<-list()
  insitu.hot.time.bins.counts.rev.cumsum.CI<-list()
  insitu.realm.time.bins.counts.rev.cumsum.CI<-list()
  dispersal.into.hot.time.bins.counts.rev.cumsum.CI<-list()
  dispersal.into.realm.time.bins.counts.rev.cumsum.CI<-list()
  dispersal.from.realm.time.bins.counts.rev.cumsum.CI<-list()
  dispersal.from.hot.time.bins.counts.rev.cumsum.CI<-list()
  
  #get cumulative events counts
  for (i in 1:length(dispersal.from.realm.time.bins.counts.rev.cumsum[[1]])){
    dispersal.from.realm.time.bins.counts.rev.cumsum.CI[[i]]<-sapply(dispersal.from.realm.time.bins.counts.rev.cumsum,"[[",i)
    dispersal.from.hot.time.bins.counts.rev.cumsum.CI[[i]]<-sapply(dispersal.from.hot.time.bins.counts.rev.cumsum,"[[",i)
    dispersal.into.realm.time.bins.counts.rev.cumsum.CI[[i]]<-sapply(dispersal.into.realm.time.bins.counts.rev.cumsum,"[[",i)
    dispersal.into.hot.time.bins.counts.rev.cumsum.CI[[i]]<-sapply(dispersal.into.hot.time.bins.counts.rev.cumsum,"[[",i)
    hot.into.realm.time.bins.counts.rev.cumsum.CI[[i]]<-sapply(hot.into.realm.time.bins.counts.rev.cumsum,"[[",i)
    realm.into.hot.time.bins.counts.rev.cumsum.CI[[i]]<-sapply(realm.into.hot.time.bins.counts.rev.cumsum,"[[",i)
    insitu.realm.time.bins.counts.rev.cumsum.CI[[i]]<-sapply(insitu.realm.time.bins.counts.rev.cumsum,"[[",i)
    insitu.hot.time.bins.counts.rev.cumsum.CI[[i]]<-sapply(insitu.hot.time.bins.counts.rev.cumsum,"[[",i)
  }
  
  
  #dispersal.from.realm.time.bins.counts.rev.cumsum.CI<-lapply(dispersal.from.realm.time.bins.counts.rev.cumsum.CI,function(x)quantile(x,c(0.05,0.95,0.50)))
  #dispersal.from.hot.time.bins.counts.rev.cumsum.CI<-lapply(dispersal.from.hot.time.bins.counts.rev.cumsum.CI,function(x)quantile(x,c(0.05,0.95,0.50)))
  #dispersal.into.realm.time.bins.counts.rev.cumsum.CI<-lapply(dispersal.into.realm.time.bins.counts.rev.cumsum.CI,function(x)quantile(x,c(0.05,0.95,0.50)))
  #dispersal.into.hot.time.bins.counts.rev.cumsum.CI<-lapply(dispersal.into.hot.time.bins.counts.rev.cumsum.CI,function(x)quantile(x,c(0.05,0.95,0.50)))
  #hot.into.realm.time.bins.counts.rev.cumsum.CI<-lapply(hot.into.realm.time.bins.counts.rev.cumsum.CI,function(x)quantile(x,c(0.05,0.95,0.50)))
  #realm.into.hot.time.bins.counts.rev.cumsum.CI<-lapply(realm.into.hot.time.bins.counts.rev.cumsum.CI,function(x)quantile(x,c(0.05,0.95,0.50)))
  #insitu.realm.time.bins.counts.rev.cumsum.CI<-lapply(insitu.realm.time.bins.counts.rev.cumsum.CI,function(x)quantile(x,c(0.05,0.95,0.50)))
  #insitu.hot.time.bins.counts.rev.cumsum.CI<-lapply(insitu.hot.time.bins.counts.rev.cumsum.CI,function(x)quantile(x,c(0.05,0.95,0.50)))
  
  rates.list<-list(insitu.realm.rates.BSM.CI,insitu.hot.rates.BSM.CI,hot.into.realm.rates.BSM.CI,realm.into.hot.rates.BSM.CI,dispersal.into.realm.rates.BSM.CI,dispersal.into.hot.rates.BSM.CI,insitu.realm.rates.strict.BAMM.CI,insitu.hot.rates.strict.BAMM.CI,insitu.realm.rates.BAMM.CI,insitu.hot.rates.BAMM.CI)
  names(rates.list)<-c("insitu.realm.rates.BSM.CI","insitu.hot.rates.BSM.CI","hot.into.realm.rates.BSM.CI","realm.into.hot.rates.BSM.CI","dispersal.into.realm.rates.BSM.CI","dispersal.into.hot.rates.BSM.CI","insitu.realm.rates.strict.BAMM.CI","insitu.hot.rates.strict.BAMM.CI","insitu.realm.rates.BAMM.CI","insitu.hot.rates.BAMM.CI")
  events.list<-list(hot.into.realm.time.bins.counts.rev.cumsum.CI,realm.into.hot.time.bins.counts.rev.cumsum.CI,insitu.hot.time.bins.counts.rev.cumsum.CI,insitu.realm.time.bins.counts.rev.cumsum.CI,dispersal.into.hot.time.bins.counts.rev.cumsum.CI,dispersal.into.realm.time.bins.counts.rev.cumsum.CI,dispersal.from.realm.time.bins.counts.rev.cumsum.CI,dispersal.from.hot.time.bins.counts.rev.cumsum.CI)
  names(events.list)<-c("hot.into.realm.time.bins.counts.rev.cumsum.CI","realm.into.hot.time.bins.counts.rev.cumsum.CI","insitu.hot.time.bins.counts.rev.cumsum.CI","insitu.realm.time.bins.counts.rev.cumsum.CI","dispersal.into.hot.time.bins.counts.rev.cumsum.CI","dispersal.into.realm.time.bins.counts.rev.cumsum.CI","dispersal.from.realm.time.bins.counts.rev.cumsum.CI","dispersal.from.hot.time.bins.counts.rev.cumsum.CI")
  output.list<-list(rates.list,events.list,realm.lineages.counts.BSM,hot.lineages.counts.BSM)
  return(output.list)
}

get_BSM_BAMM_fullrates_events_2my_object<-function(results.BSM,name,speciation.rate.treefile){
  if(class(speciation.rate.treefile)=='character'){
    cat('getting speciation rates','\n')
    speciation.rate.tree<-read.tree(speciation.rate.treefile)  
    table.speciation.rate<-prt(speciation.rate.tree,printflag = FALSE)
  }else if (class(speciation.rate.treefile)=='data.frame'){
    table.speciation.rate<-speciation.rate.treefile
  }
  area.codes<-c('O','A','B','C','D','E','F','G','AB','AC','AD','AE','AF','AG','BC','BD','BE','BF','BG','CD','CE','CF','CG','DE','DF','DG','EF','EG','FG')
  #get event names from BSM
  if(name=='afrotrop'){
    realm.code<-'A'
  }else if(name=='austral'){
    realm.code<-'B'
  }else if(name=='indo'){
    realm.code<-'C'
  }else if(name=='nearctic'){
    realm.code<-'D'
  }else if(name=='neotrop'){
    realm.code<-'E'
  }else if(name=='palearctic'){
    realm.code<-'F'
  }
  hot.code<-'G'
  #read clado and ana events table of each BSM replicate
  insitu.realm.rates.BSM<-list()
  insitu.hot.rates.BSM<-list()
  insitu.realm.rates.BAMM.strict<-list()
  insitu.hot.rates.BAMM.strict<-list()
  insitu.realm.rates.BAMM<-list()
  insitu.hot.rates.BAMM<-list()
  dispersal.into.realm.rates.BSM<-list()
  dispersal.into.hot.rates.BSM<-list()
  dispersal.from.realm.rates.BSM<-list()
  dispersal.from.hot.rates.BSM<-list()
  hot.into.realm.rates.BSM<-list()
  realm.into.hot.rates.BSM<-list()
  
  dispersal.into.hot.time.bins.counts.rev.cumsum<-list()
  dispersal.into.realm.time.bins.counts.rev.cumsum<-list()
  hot.into.realm.time.bins.counts.rev.cumsum<-list()
  realm.into.hot.time.bins.counts.rev.cumsum<-list()
  insitu.hot.time.bins.counts.rev.cumsum<-list()
  insitu.realm.time.bins.counts.rev.cumsum<-list()
  dispersal.from.hot.time.bins.counts.rev.cumsum<-list()
  dispersal.from.realm.time.bins.counts.rev.cumsum<-list()
  realm.into.hot.time.bins.counts.rev.cumsum<-list()
  into.into.realm.time.bins.counts.rev.cumsum<-list()
  realm.lineages.counts.BSM<-list()
  hot.lineages.counts.BSM<-list()
  cat('analysing BSMs','\n')
  for (i in 1:length(results.BSM$RES_clado_events_tables)){
    ###BioGeoBEARS
    #cat(i,'\n')
    a1<-results.BSM$RES_clado_events_tables[[i]]
    a2<-results.BSM$RES_ana_events_tables[[i]]
    #assign events to time bin
    a1$time.bin<-floor(a1$time_bp/2)*2
    a2$time.bin<-floor(a2$time_bp/2)*2
    max.age<-max(a1$time.bin)
    #get time bins of each event type
    #dispersal into hot
    dispersal.into.hot.time.bins<-sort(c(a1[a1$clado_dispersal_to==hot.code,'time.bin'],a2[a2$dispersal_to==hot.code,'time.bin']))
    dispersal.into.hot.time.bins.counts<-sapply(seq(from=0,to=max.age,by=2),function(x)length(dispersal.into.hot.time.bins[dispersal.into.hot.time.bins%in%x]))
    dispersal.into.hot.time.bins.counts.rev<-rev(dispersal.into.hot.time.bins.counts)
    dispersal.into.hot.time.bins.counts.rev.cumsum[[i]]<-cumsum(dispersal.into.hot.time.bins.counts.rev)
    #dispersal into realm
    dispersal.into.realm.time.bins<-sort(c(a1[a1$clado_dispersal_to==realm.code,'time.bin'],a2[a2$dispersal_to==realm.code,'time.bin']))
    dispersal.into.realm.time.bins.counts<-sapply(seq(from=0,to=max.age,by=2),function(x)length(dispersal.into.realm.time.bins[dispersal.into.realm.time.bins%in%x]))
    dispersal.into.realm.time.bins.counts.rev<-rev(dispersal.into.realm.time.bins.counts)
    dispersal.into.realm.time.bins.counts.rev.cumsum[[i]]<-cumsum(dispersal.into.realm.time.bins.counts.rev)
    #hot local extinctions
    hot.local.extinctions.time.bins<-sort(a2[a2$extirpation_from==hot.code,'time.bin'])
    hot.local.extinctions.time.bins.counts<-sapply(seq(from=0,to=max.age,by=2),function(x)length(hot.local.extinctions.time.bins[hot.local.extinctions.time.bins%in%x]))
    hot.local.extinctions.time.bins.counts.rev<-rev(hot.local.extinctions.time.bins.counts)
    hot.local.extinctions.time.bins.counts.rev.cumsum<-cumsum(hot.local.extinctions.time.bins.counts.rev)
    #realm local extinctions
    realm.local.extinctions.time.bins<-sort(a2[a2$extirpation_from==realm.code,'time.bin'])
    realm.local.extinctions.time.bins.counts<-sapply(seq(from=0,to=max.age,by=2),function(x)length(realm.local.extinctions.time.bins[realm.local.extinctions.time.bins%in%x]))
    realm.local.extinctions.time.bins.counts.rev<-rev(realm.local.extinctions.time.bins.counts)
    realm.local.extinctions.time.bins.counts.rev.cumsum<-cumsum(realm.local.extinctions.time.bins.counts.rev)
    #hot into realm dispersal
    hot.into.realm.time.bins<-sort(c(a2[intersect(intersect(grep(hot.code,a2$current_rangetxt),grep(realm.code,a2$current_rangetxt,invert=TRUE)),grep(realm.code,a2$new_rangetxt)),'time.bin'],a1[a1$sampled_states_AT_nodes%in%intersect(grep(hot.code,area.codes),grep(realm.code,area.codes,invert=TRUE))&a1$clado_dispersal_to==realm.code,'time.bin']))
    hot.into.realm.time.bins.counts<-sapply(seq(from=0,to=max.age,by=2),function(x)length(hot.into.realm.time.bins[hot.into.realm.time.bins%in%x]))
    hot.into.realm.time.bins.counts.rev<-rev(hot.into.realm.time.bins.counts)
    hot.into.realm.time.bins.counts.rev.cumsum[[i]]<-cumsum(hot.into.realm.time.bins.counts.rev)
    #realm into hot dispersal
    realm.into.hot.time.bins<-sort(c(a2[intersect(intersect(grep(realm.code,a2$current_rangetxt),grep(hot.code,a2$current_rangetxt,invert=TRUE)),grep(hot.code,a2$new_rangetxt)),'time.bin'],a1[a1$sampled_states_AT_nodes%in%intersect(grep(realm.code,area.codes),grep(hot.code,area.codes,invert=TRUE))&a1$clado_dispersal_to==hot.code,'time.bin']))
    realm.into.hot.time.bins.counts<-sapply(seq(from=0,to=max.age,by=2),function(x)length(realm.into.hot.time.bins[realm.into.hot.time.bins%in%x]))
    realm.into.hot.time.bins.counts.rev<-rev(realm.into.hot.time.bins.counts)
    realm.into.hot.time.bins.counts.rev.cumsum[[i]]<-cumsum(realm.into.hot.time.bins.counts.rev)
    #in situ hot cladogenesis
    insitu.hot.time.bins<-sort(a1[a1$clado_dispersal_to==''&(a1$sampled_states_AT_nodes%in%grep(hot.code,area.codes))&a1$clado_event_txt!='','time.bin'])
    insitu.hot.time.bins.counts<-sapply(seq(from=0,to=max.age,by=2),function(x)length(insitu.hot.time.bins[insitu.hot.time.bins%in%x]))
    insitu.hot.time.bins.counts.rev<-rev(insitu.hot.time.bins.counts)
    insitu.hot.time.bins.counts.rev.cumsum[[i]]<-cumsum(insitu.hot.time.bins.counts.rev)
    #in situ realm cladogenesis
    insitu.realm.time.bins<-sort(a1[a1$clado_dispersal_to==''&(a1$sampled_states_AT_nodes%in%grep(realm.code,area.codes))&a1$clado_event_txt!='','time.bin'])
    insitu.realm.time.bins.counts<-sapply(seq(from=0,to=max.age,by=2),function(x)length(insitu.realm.time.bins[insitu.realm.time.bins%in%x]))
    insitu.realm.time.bins.counts.rev<-rev(insitu.realm.time.bins.counts)
    insitu.realm.time.bins.counts.rev.cumsum[[i]]<-cumsum(insitu.realm.time.bins.counts.rev)
    #dispersal from hot
    dispersal.from.hot.time.bins<-sort(c(a2[a2$extirpation_from=='-'&(a2$current_rangetxt%in%area.codes[grep(hot.code,area.codes)]),'time.bin'],a1[(a1$clado_event_txt!='')&(a1$sampled_states_AT_nodes%in%grep(hot.code,area.codes))&(!(a1$samp_LEFT_dcorner%in%grep(hot.code,area.codes))|!(a1$samp_RIGHT_dcorner%in%grep(hot.code,area.codes))),'time.bin']))
    dispersal.from.hot.time.bins.counts<-sapply(seq(from=0,to=max.age,by=2),function(x)length(dispersal.from.hot.time.bins[dispersal.from.hot.time.bins%in%x]))
    dispersal.from.hot.time.bins.counts.rev<-rev(dispersal.from.hot.time.bins.counts)
    dispersal.from.hot.time.bins.counts.rev.cumsum[[i]]<-cumsum(dispersal.from.hot.time.bins.counts.rev)
    #dispersal from realm
    dispersal.from.realm.time.bins<-sort(c(a2[a2$extirpation_from=='-'&(a2$current_rangetxt%in%area.codes[grep(realm.code,area.codes)]),'time.bin'],a1[(a1$clado_event_txt!='')&(a1$sampled_states_AT_nodes%in%grep(realm.code,area.codes))&(!(a1$samp_LEFT_dcorner%in%grep(realm.code,area.codes))|!(a1$samp_RIGHT_dcorner%in%grep(realm.code,area.codes))),'time.bin']))
    dispersal.from.realm.time.bins.counts<-sapply(seq(from=0,to=max.age,by=2),function(x)length(dispersal.from.realm.time.bins[dispersal.from.realm.time.bins%in%x]))
    dispersal.from.realm.time.bins.counts.rev<-rev(dispersal.from.realm.time.bins.counts)
    dispersal.from.realm.time.bins.counts.rev.cumsum[[i]]<-cumsum(dispersal.from.realm.time.bins.counts.rev)
    #calculate number of lineages in hot and realm at each time bin
    realm.lineages.counts<-insitu.realm.time.bins.counts.rev.cumsum[[i]]+dispersal.into.realm.time.bins.counts.rev.cumsum[[i]]-realm.local.extinctions.time.bins.counts.rev.cumsum
    hot.lineages.counts<-insitu.hot.time.bins.counts.rev.cumsum[[i]]+dispersal.into.hot.time.bins.counts.rev.cumsum[[i]]-hot.local.extinctions.time.bins.counts.rev.cumsum
    realm.lineages.counts.BSM[[i]]<-list()
    hot.lineages.counts.BSM[[i]]<-list()
    realm.lineages.counts.BSM[[i]]<-realm.lineages.counts
    hot.lineages.counts.BSM[[i]]<-hot.lineages.counts
    #calculate in.situ.rates.BSM
    insitu.realm.rates.BSM[[i]]<-unlist(lapply(2:length(realm.lineages.counts),function(x)insitu.realm.time.bins.counts.rev[x]/realm.lineages.counts[x-1]))
    insitu.hot.rates.BSM[[i]]<-unlist(lapply(2:length(hot.lineages.counts),function(x)insitu.hot.time.bins.counts.rev[x]/hot.lineages.counts[x-1]))
    #calculate colonisation rates BSM
    dispersal.into.realm.rates.BSM[[i]]<-unlist(lapply(2:length(realm.lineages.counts),function(x)dispersal.into.realm.time.bins.counts.rev[x]/realm.lineages.counts[x-1]))
    dispersal.into.hot.rates.BSM[[i]]<-unlist(lapply(2:length(hot.lineages.counts),function(x)dispersal.into.hot.time.bins.counts.rev[x]/hot.lineages.counts[x-1]))
    dispersal.from.realm.rates.BSM[[i]]<-unlist(lapply(2:length(realm.lineages.counts),function(x)dispersal.from.realm.time.bins.counts.rev[x]/realm.lineages.counts[x-1]))
    dispersal.from.hot.rates.BSM[[i]]<-unlist(lapply(2:length(hot.lineages.counts),function(x)dispersal.from.hot.time.bins.counts.rev[x]/hot.lineages.counts[x-1]))
    
    #calculate dispersal rates hot->realm and viceversa BSM
    hot.into.realm.rates.BSM[[i]]<-unlist(lapply(2:length(realm.lineages.counts),function(x)hot.into.realm.time.bins.counts.rev[x]/realm.lineages.counts[x-1]))
    realm.into.hot.rates.BSM[[i]]<-unlist(lapply(2:length(hot.lineages.counts),function(x)realm.into.hot.time.bins.counts.rev[x]/hot.lineages.counts[x-1]))
    
    ###BAMM
    #assign regions to each node
    regions<-vector()
    for(a in 1:nrow(a1)){
      if(a1[a,'node.type']=='root'){
        regions[a]<-area.codes[a1[a,'sampled_states_AT_nodes']]
      }else if((a1[a,'sampled_states_AT_nodes']==a1[a,'sampled_states_AT_brbots'])&(a1[a,'anagenetic_events_txt_below_node']=='none')){
        regions[a]<-area.codes[a1[a,'sampled_states_AT_nodes']]
      }else{
        regions[a]<-NA
      }#this ignores anagenetic changes, only assigns rates to branches that do not change state
      #else{
      #   new.range<-strsplit(a1[a,'anagenetic_events_txt_below_node'],split=',')[[1]][grep('new_rangetxt',strsplit(a1[a,'anagenetic_events_txt_below_node'],split=',')[[1]])[length(grep('new_rangetxt',strsplit(a1[a,'anagenetic_events_txt_below_node'],split=',')[[1]]))]]
      #   regions[a]<-gsub(new.range,pattern='new_rangetxt:',replacement='')
      #  }
    }
    a1$speciation.rate<-table.speciation.rate$edge.length
    a1$region<-regions
    #assign events to time bin
    a1$time.bin<-floor(a1$time_bp/2)*2
    #get mean rates per region (strict = species only can occur in a single region)
    if(length(a1$region[a1$region%in%c(realm.code)])>1){
      insitu.realm.rates.BAMM.strict[[i]]<-aggregate(speciation.rate~time.bin,data=a1[a1$region%in%c(realm.code),],function(x)mean(x))
    }else if (length(a1$region[a1$region%in%c(realm.code)])==1){
      insitu.realm.rates.BAMM.strict[[i]]<-as.data.frame(cbind(a1[a1$region%in%c(realm.code),'time.bin'],a1[a1$region%in%c(realm.code),'speciation.rate']),stringsAsFactors = F)
      colnames(insitu.realm.rates.BAMM.strict[[i]])<-c('time.bin','speciation.rate')
    }else{
      insitu.realm.rates.BAMM.strict[[i]]<-as.data.frame(cbind(1,0),stringsAsFactors = F)
      colnames(insitu.realm.rates.BAMM.strict[[i]])<-c('time.bin','speciation.rate')
    }
    if(length(a1$region[a1$region%in%c(hot.code)])>1){
      insitu.hot.rates.BAMM.strict[[i]]<-aggregate(speciation.rate~time.bin,data=a1[a1$region%in%c(hot.code),],function(x)mean(x))
    }else if (length(a1$region[a1$region%in%c(hot.code)])==1){
      insitu.hot.rates.BAMM.strict[[i]]<-as.data.frame(cbind(a1[a1$region%in%c(hot.code),'time.bin'],a1[a1$region%in%c(hot.code),'speciation.rate']),stringsAsFactors = F)
      colnames(insitu.hot.rates.BAMM.strict[[i]])<-c('time.bin','speciation.rate')
    }else{
      insitu.hot.rates.BAMM.strict[[i]]<-as.data.frame(cbind(1,0),stringsAsFactors = F)
      colnames(insitu.hot.rates.BAMM.strict[[i]])<-c('time.bin','speciation.rate')
    }
    if(length(a1$region[a1$region%in%c(area.codes[grep(realm.code,area.codes)])])>1){
      insitu.realm.rates.BAMM[[i]]<-aggregate(speciation.rate~time.bin,data=a1[a1$region%in%c(area.codes[grep(realm.code,area.codes)]),],function(x)mean(x))
    }else if (length(a1$region[a1$region%in%c(area.codes[grep(realm.code,area.codes)])])==1){
      insitu.realm.rates.BAMM[[i]]<-as.data.frame(cbind(a1[a1$region%in%c(area.codes[grep(realm.code,area.codes)]),'time.bin'],a1[a1$region%in%c(area.codes[grep(realm.code,area.codes)]),'speciation.rate']),stringsAsFactors = F)
      colnames(insitu.realm.rates.BAMM[[i]])<-c('time.bin','speciation.rate')
    }else{
      insitu.realm.rates.BAMM[[i]]<-as.data.frame(cbind(1,0),stringsAsFactors = F)
      colnames(insitu.realm.rates.BAMM[[i]])<-c('time.bin','speciation.rate')
    }
    if(length(a1$region[a1$region%in%c(area.codes[grep(hot.code,area.codes)])])>1){
      insitu.hot.rates.BAMM[[i]]<-aggregate(speciation.rate~time.bin,data=a1[a1$region%in%c(area.codes[grep(hot.code,area.codes)]),],function(x)mean(x))
    }else if (length(a1$region[a1$region%in%c(area.codes[grep(hot.code,area.codes)])])==1){
      insitu.hot.rates.BAMM[[i]]<-as.data.frame(cbind(a1[a1$region%in%c(area.codes[grep(hot.code,area.codes)]),'time.bin'],a1[a1$region%in%c(area.codes[grep(hot.code,area.codes)]),'speciation.rate']),stringsAsFactors = F)
      colnames(insitu.hot.rates.BAMM[[i]])<-c('time.bin','speciation.rate')
    }else{
      insitu.hot.rates.BAMM[[i]]<-as.data.frame(cbind(1,0),stringsAsFactors = F)
      colnames(insitu.hot.rates.BAMM[[i]])<-c('time.bin','speciation.rate')
    }
    
    
  }
  #add 0s to empty bins
  insitu.realm.rates.BAMM<-lapply(insitu.realm.rates.BAMM,function(x) {missing.time.bins<-seq(from=0,to=max.age,by=2)[!c(seq(from=0,to=max.age,by=2)%in%x$time.bin)];missing.df<-as.data.frame(cbind(missing.time.bins,rep(0,length(missing.time.bins))));colnames(missing.df)<-c('time.bin','speciation.rate');new.df<-rbind(x,missing.df);new.df<-new.df[order(new.df$time.bin),];return(rev(new.df$speciation.rate))})
  insitu.hot.rates.BAMM<-lapply(insitu.hot.rates.BAMM,function(x) {missing.time.bins<-seq(from=0,to=max.age,by=2)[!c(seq(from=0,to=max.age,by=2)%in%x$time.bin)];missing.df<-as.data.frame(cbind(missing.time.bins,rep(0,length(missing.time.bins))));colnames(missing.df)<-c('time.bin','speciation.rate');new.df<-rbind(x,missing.df);new.df<-new.df[order(new.df$time.bin),];return(rev(new.df$speciation.rate))})
  #add 0s to empty bins
  insitu.realm.rates.BAMM.strict<-lapply(insitu.realm.rates.BAMM.strict,function(x) {missing.time.bins<-seq(from=0,to=max.age,by=2)[!c(seq(from=0,to=max.age,by=2)%in%x$time.bin)];missing.df<-as.data.frame(cbind(missing.time.bins,rep(0,length(missing.time.bins))));colnames(missing.df)<-c('time.bin','speciation.rate');new.df<-rbind(x,missing.df);new.df<-new.df[order(new.df$time.bin),];return(rev(new.df$speciation.rate))})
  insitu.hot.rates.BAMM.strict<-lapply(insitu.hot.rates.BAMM.strict,function(x) {missing.time.bins<-seq(from=0,to=max.age,by=2)[!c(seq(from=0,to=max.age,by=2)%in%x$time.bin)];missing.df<-as.data.frame(cbind(missing.time.bins,rep(0,length(missing.time.bins))));colnames(missing.df)<-c('time.bin','speciation.rate');new.df<-rbind(x,missing.df);new.df<-new.df[order(new.df$time.bin),];return(rev(new.df$speciation.rate))})
  
  #get the BSM rates CIs to draw on plot
  insitu.realm.rates.BSM.CI<-list()
  insitu.hot.rates.BSM.CI<-list()
  for (i in 1:length(insitu.realm.rates.BSM[[1]])){
    insitu.realm.rates.BSM.CI[[i]]<-sapply(insitu.realm.rates.BSM,"[[",i)
    insitu.hot.rates.BSM.CI[[i]]<-sapply(insitu.hot.rates.BSM,"[[",i)
  }
  #insitu.realm.rates.BSM.CI<-lapply(insitu.realm.rates.BSM.CI,function(x)quantile(x,c(0.05,0.95,0.50),na.rm=TRUE))
  #insitu.hot.rates.BSM.CI<-lapply(insitu.hot.rates.BSM.CI,function(x)quantile(x,c(0.05,0.95,0.50),na.rm=TRUE))
  #get the BSM dispersal rates CIs to draw on plot
  hot.into.realm.rates.BSM.CI<-list()
  realm.into.hot.rates.BSM.CI<-list()
  for (i in 1:length(hot.into.realm.rates.BSM[[1]])){
    hot.into.realm.rates.BSM.CI[[i]]<-sapply(hot.into.realm.rates.BSM,"[[",i)
    realm.into.hot.rates.BSM.CI[[i]]<-sapply(realm.into.hot.rates.BSM,"[[",i)
  }
  #hot.into.realm.rates.BSM.CI<-lapply(hot.into.realm.rates.BSM.CI,function(x)quantile(x,c(0.05,0.95,0.50),na.rm=TRUE))
  #realm.into.hot.rates.BSM.CI<-lapply(realm.into.hot.rates.BSM.CI,function(x)quantile(x,c(0.05,0.95,0.50),na.rm=TRUE))
  dispersal.into.realm.rates.BSM.CI<-list()
  dispersal.into.hot.rates.BSM.CI<-list()
  for (i in 1:length(dispersal.into.realm.rates.BSM[[1]])){
    dispersal.into.realm.rates.BSM.CI[[i]]<-sapply(dispersal.into.realm.rates.BSM,"[[",i)
    dispersal.into.hot.rates.BSM.CI[[i]]<-sapply(dispersal.into.hot.rates.BSM,"[[",i)
  }
  dispersal.from.realm.rates.BSM.CI<-list()
  dispersal.from.hot.rates.BSM.CI<-list()
  for (i in 1:length(dispersal.from.realm.rates.BSM[[1]])){
    dispersal.from.realm.rates.BSM.CI[[i]]<-sapply(dispersal.from.realm.rates.BSM,"[[",i)
    dispersal.from.hot.rates.BSM.CI[[i]]<-sapply(dispersal.from.hot.rates.BSM,"[[",i)
  }
  
  #dispersal.into.realm.rates.BSM.CI<-lapply(dispersal.into.realm.rates.BSM.CI,function(x)quantile(x,c(0.05,0.95,0.50),na.rm=TRUE))
  #dispersal.into.hot.rates.BSM.CI<-lapply(dispersal.into.hot.rates.BSM.CI,function(x)quantile(x,c(0.05,0.95,0.50),na.rm=TRUE))
  #get the BAMM strict rates CIs to draw on plot
  
  insitu.realm.rates.strict.BAMM.CI<-list()
  insitu.hot.rates.strict.BAMM.CI<-list()
  for (i in 1:length(insitu.realm.rates.BAMM.strict[[1]])){
    insitu.realm.rates.strict.BAMM.CI[[i]]<-sapply(insitu.realm.rates.BAMM.strict,"[[",i)
    insitu.hot.rates.strict.BAMM.CI[[i]]<-sapply(insitu.hot.rates.BAMM.strict,"[[",i)
  }
  #insitu.realm.rates.strict.BAMM.CI<-lapply(insitu.realm.rates.strict.BAMM.CI,function(x)quantile(x,c(0.05,0.95,0.50),na.rm=TRUE))
  #insitu.hot.rates.strict.BAMM.CI<-lapply(insitu.hot.rates.strict.BAMM.CI,function(x)quantile(x,c(0.05,0.95,0.50),na.rm=TRUE))
  #get the BAMM rates CIs to draw on plot
  insitu.realm.rates.BAMM.CI<-list()
  insitu.hot.rates.BAMM.CI<-list()
  for (i in 1:length(insitu.realm.rates.BAMM[[1]])){
    insitu.realm.rates.BAMM.CI[[i]]<-sapply(insitu.realm.rates.BAMM,"[[",i)
    insitu.hot.rates.BAMM.CI[[i]]<-sapply(insitu.hot.rates.BAMM,"[[",i)
  }
  #insitu.realm.rates.BAMM.CI<-lapply(insitu.realm.rates.BAMM.CI,function(x)quantile(x,c(0.05,0.95,0.50),na.rm=TRUE))
  #insitu.hot.rates.BAMM.CI<-lapply(insitu.hot.rates.BAMM.CI,function(x)quantile(x,c(0.05,0.95,0.50),na.rm=TRUE))
  
  hot.into.realm.time.bins.counts.rev.cumsum.CI<-list()
  realm.into.hot.time.bins.counts.rev.cumsum.CI<-list()
  insitu.hot.time.bins.counts.rev.cumsum.CI<-list()
  insitu.realm.time.bins.counts.rev.cumsum.CI<-list()
  dispersal.into.hot.time.bins.counts.rev.cumsum.CI<-list()
  dispersal.into.realm.time.bins.counts.rev.cumsum.CI<-list()
  dispersal.from.realm.time.bins.counts.rev.cumsum.CI<-list()
  dispersal.from.hot.time.bins.counts.rev.cumsum.CI<-list()
  
  #get cumulative events counts
  for (i in 1:length(dispersal.from.realm.time.bins.counts.rev.cumsum[[1]])){
    dispersal.from.realm.time.bins.counts.rev.cumsum.CI[[i]]<-sapply(dispersal.from.realm.time.bins.counts.rev.cumsum,"[[",i)
    dispersal.from.hot.time.bins.counts.rev.cumsum.CI[[i]]<-sapply(dispersal.from.hot.time.bins.counts.rev.cumsum,"[[",i)
    dispersal.into.realm.time.bins.counts.rev.cumsum.CI[[i]]<-sapply(dispersal.into.realm.time.bins.counts.rev.cumsum,"[[",i)
    dispersal.into.hot.time.bins.counts.rev.cumsum.CI[[i]]<-sapply(dispersal.into.hot.time.bins.counts.rev.cumsum,"[[",i)
    hot.into.realm.time.bins.counts.rev.cumsum.CI[[i]]<-sapply(hot.into.realm.time.bins.counts.rev.cumsum,"[[",i)
    realm.into.hot.time.bins.counts.rev.cumsum.CI[[i]]<-sapply(realm.into.hot.time.bins.counts.rev.cumsum,"[[",i)
    insitu.realm.time.bins.counts.rev.cumsum.CI[[i]]<-sapply(insitu.realm.time.bins.counts.rev.cumsum,"[[",i)
    insitu.hot.time.bins.counts.rev.cumsum.CI[[i]]<-sapply(insitu.hot.time.bins.counts.rev.cumsum,"[[",i)
  }
  
  
  #dispersal.from.realm.time.bins.counts.rev.cumsum.CI<-lapply(dispersal.from.realm.time.bins.counts.rev.cumsum.CI,function(x)quantile(x,c(0.05,0.95,0.50)))
  #dispersal.from.hot.time.bins.counts.rev.cumsum.CI<-lapply(dispersal.from.hot.time.bins.counts.rev.cumsum.CI,function(x)quantile(x,c(0.05,0.95,0.50)))
  #dispersal.into.realm.time.bins.counts.rev.cumsum.CI<-lapply(dispersal.into.realm.time.bins.counts.rev.cumsum.CI,function(x)quantile(x,c(0.05,0.95,0.50)))
  #dispersal.into.hot.time.bins.counts.rev.cumsum.CI<-lapply(dispersal.into.hot.time.bins.counts.rev.cumsum.CI,function(x)quantile(x,c(0.05,0.95,0.50)))
  #hot.into.realm.time.bins.counts.rev.cumsum.CI<-lapply(hot.into.realm.time.bins.counts.rev.cumsum.CI,function(x)quantile(x,c(0.05,0.95,0.50)))
  #realm.into.hot.time.bins.counts.rev.cumsum.CI<-lapply(realm.into.hot.time.bins.counts.rev.cumsum.CI,function(x)quantile(x,c(0.05,0.95,0.50)))
  #insitu.realm.time.bins.counts.rev.cumsum.CI<-lapply(insitu.realm.time.bins.counts.rev.cumsum.CI,function(x)quantile(x,c(0.05,0.95,0.50)))
  #insitu.hot.time.bins.counts.rev.cumsum.CI<-lapply(insitu.hot.time.bins.counts.rev.cumsum.CI,function(x)quantile(x,c(0.05,0.95,0.50)))
  
  rates.list<-list(insitu.realm.rates.BSM.CI,insitu.hot.rates.BSM.CI,hot.into.realm.rates.BSM.CI,realm.into.hot.rates.BSM.CI,dispersal.into.realm.rates.BSM.CI,dispersal.into.hot.rates.BSM.CI,dispersal.from.realm.rates.BSM.CI,dispersal.from.hot.rates.BSM.CI,insitu.realm.rates.strict.BAMM.CI,insitu.hot.rates.strict.BAMM.CI,insitu.realm.rates.BAMM.CI,insitu.hot.rates.BAMM.CI)
  names(rates.list)<-c("insitu.realm.rates.BSM.CI","insitu.hot.rates.BSM.CI","hot.into.realm.rates.BSM.CI","realm.into.hot.rates.BSM.CI","dispersal.into.realm.rates.BSM.CI","dispersal.into.hot.rates.BSM.CI","dispersal.from.realm.rates.BSM.CI","dispersal.from.hot.rates.BSM.CI","insitu.realm.rates.strict.BAMM.CI","insitu.hot.rates.strict.BAMM.CI","insitu.realm.rates.BAMM.CI","insitu.hot.rates.BAMM.CI")
  events.list<-list(hot.into.realm.time.bins.counts.rev.cumsum.CI,realm.into.hot.time.bins.counts.rev.cumsum.CI,insitu.hot.time.bins.counts.rev.cumsum.CI,insitu.realm.time.bins.counts.rev.cumsum.CI,dispersal.into.hot.time.bins.counts.rev.cumsum.CI,dispersal.into.realm.time.bins.counts.rev.cumsum.CI,dispersal.from.realm.time.bins.counts.rev.cumsum.CI,dispersal.from.hot.time.bins.counts.rev.cumsum.CI)
  names(events.list)<-c("hot.into.realm.time.bins.counts.rev.cumsum.CI","realm.into.hot.time.bins.counts.rev.cumsum.CI","insitu.hot.time.bins.counts.rev.cumsum.CI","insitu.realm.time.bins.counts.rev.cumsum.CI","dispersal.into.hot.time.bins.counts.rev.cumsum.CI","dispersal.into.realm.time.bins.counts.rev.cumsum.CI","dispersal.from.realm.time.bins.counts.rev.cumsum.CI","dispersal.from.hot.time.bins.counts.rev.cumsum.CI")
  output.list<-list(rates.list,events.list,realm.lineages.counts.BSM,hot.lineages.counts.BSM)
  return(output.list)
}

plot_rates_events_real<-function(BSMoutput.real.rates,BSMoutput.real.events,age,name,plots){
  #correct age
  age<-length(BSMoutput.real.rates[['insitu.realm.rates.BSM.CI']])-age
  #trim lists to age period
  BSMoutput.real.rates<-lapply(BSMoutput.real.rates,function(x)x[age:length(x)])
  BSMoutput.real.events<-lapply(BSMoutput.real.events,function(x)x[age:length(x)])
  #remove last two million year time bins
  BSMoutput.real.rates<-lapply(BSMoutput.real.rates,function(x) x<-x[-c((length(x)-1),length(x))])
  BSMoutput.real.events<-lapply(BSMoutput.real.events,function(x) x<-x[-c((length(x)-2),(length(x)-1),length(x))])
  #assign variables for easier handling
  hot.into.realm.rates.BSM.CI<-BSMoutput.real.rates[['hot.into.realm.rates.BSM.CI']]
  realm.into.hot.rates.BSM.CI<-BSMoutput.real.rates[['realm.into.hot.rates.BSM.CI']]
  insitu.realm.rates.BAMM.CI<-BSMoutput.real.rates[['insitu.realm.rates.BAMM.CI']]
  insitu.hot.rates.BAMM.CI<-BSMoutput.real.rates[['insitu.hot.rates.BAMM.CI']]
  insitu.realm.rates.BSM.CI<-BSMoutput.real.rates[['insitu.realm.rates.BSM.CI']]
  insitu.hot.rates.BSM.CI<-BSMoutput.real.rates[['insitu.hot.rates.BSM.CI']]
  dispersal.into.realm.rates.BSM.CI<-BSMoutput.real.rates[['dispersal.into.realm.rates.BSM.CI']]
  dispersal.into.hot.rates.BSM.CI<-BSMoutput.real.rates[['dispersal.into.hot.rates.BSM.CI']]
  hot.into.realm.time.bins.counts.rev.cumsum.CI<-BSMoutput.real.events[['hot.into.realm.time.bins.counts.rev.cumsum.CI']]
  realm.into.hot.time.bins.counts.rev.cumsum.CI<-BSMoutput.real.events[['realm.into.hot.time.bins.counts.rev.cumsum.CI']]
  insitu.realm.time.bins.counts.rev.cumsum.CI<-BSMoutput.real.events[['insitu.realm.time.bins.counts.rev.cumsum.CI']]
  insitu.hot.time.bins.counts.rev.cumsum.CI<-BSMoutput.real.events[['insitu.hot.time.bins.counts.rev.cumsum.CI']]
  dispersal.into.realm.time.bins.counts.rev.cumsum.CI<-BSMoutput.real.events[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']]
  dispersal.into.hot.time.bins.counts.rev.cumsum.CI<-BSMoutput.real.events[['dispersal.into.hot.time.bins.counts.rev.cumsum.CI']]
  dispersal.from.realm.time.bins.counts.rev.cumsum.CI<-BSMoutput.real.events[['dispersal.from.realm.time.bins.counts.rev.cumsum.CI']]
  dispersal.from.hot.time.bins.counts.rev.cumsum.CI<-BSMoutput.real.events[['dispersal.from.hot.time.bins.counts.rev.cumsum.CI']]
  
  if(plots=='dispersal.rates'){
    plot(c(1,1),xlim=c(0,length(hot.into.realm.rates.BSM.CI)),ylim=c(0,max(c(unlist(hot.into.realm.rates.BSM.CI),unlist(realm.into.hot.rates.BSM.CI))[is.finite(c(unlist(hot.into.realm.rates.BSM.CI),unlist(realm.into.hot.rates.BSM.CI)))])),type='n',xaxt='n',xlab='age',ylab='dispersal rates',main=paste(name,' dispersal rates',sep=''))
    lines(x=c(1:length(hot.into.realm.rates.BSM.CI)),y=unlist(lapply(hot.into.realm.rates.BSM.CI,function(x)x[3])),col=adjustcolor( "red", alpha.f = 0.9))
    polygon(x=c((1:length(hot.into.realm.rates.BSM.CI)),(rev(1:length(hot.into.realm.rates.BSM.CI)))),y=c((unlist(lapply(hot.into.realm.rates.BSM.CI,function(x)x[1]))),rev(unlist(lapply(hot.into.realm.rates.BSM.CI,function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
    lines(x=c(1:length(realm.into.hot.rates.BSM.CI)),y=unlist(lapply(realm.into.hot.rates.BSM.CI,function(x)x[3])),col=adjustcolor( "blue", alpha.f = 0.9))
    polygon(x=c((1:length(realm.into.hot.rates.BSM.CI)),(rev(1:length(realm.into.hot.rates.BSM.CI)))),y=c((unlist(lapply(realm.into.hot.rates.BSM.CI,function(x)x[1]))),rev(unlist(lapply(realm.into.hot.rates.BSM.CI,function(x)x[2])))),col=adjustcolor( "blue", alpha.f = 0.2))
    axis(1,at=c(0:length(realm.into.hot.rates.BSM.CI)),labels=c((length(realm.into.hot.rates.BSM.CI)+2):2))
    legend('topleft',c('hot.into.realm','realm.into.hot'),col=c('red','blue'),lty=1,cex=.7,bty='n')
  }else if(plots=='cladogenesis.BAMM.rates'){
    #plot BAMM in situ rates
    plot(c(1,1),xlim=c(0,length(insitu.hot.rates.BAMM.CI)),ylim=c(0,max(c(unlist(insitu.hot.rates.BAMM.CI),unlist(insitu.realm.rates.BAMM.CI)))),type='n',xaxt='n',xlab='age',ylab='insitu.rates',main=paste(name,' cladogenesis BAMM rates',sep=''))
    lines(x=c(1:length(insitu.hot.rates.BAMM.CI)),y=unlist(lapply(insitu.hot.rates.BAMM.CI,function(x)x[3])),col=adjustcolor( "red", alpha.f = 0.9))
    polygon(x=c((1:length(insitu.hot.rates.BAMM.CI)),(rev(1:length(insitu.hot.rates.BAMM.CI)))),y=c((unlist(lapply(insitu.hot.rates.BAMM.CI,function(x)x[1]))),rev(unlist(lapply(insitu.hot.rates.BAMM.CI,function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
    lines(x=c(1:length(insitu.realm.rates.BAMM.CI)),y=unlist(lapply(insitu.realm.rates.BAMM.CI,function(x)x[3])),col=adjustcolor( "blue", alpha.f = 0.9))
    polygon(x=c((1:length(insitu.realm.rates.BAMM.CI)),(rev(1:length(insitu.realm.rates.BAMM.CI)))),y=c((unlist(lapply(insitu.realm.rates.BAMM.CI,function(x)x[1]))),rev(unlist(lapply(insitu.realm.rates.BAMM.CI,function(x)x[2])))),col=adjustcolor( "blue", alpha.f = 0.2))
    axis(1,at=c(0:length(insitu.realm.rates.BAMM.CI)),labels=c((length(insitu.realm.rates.BAMM.CI)+2):2))
    legend('topleft',c('insitu.hot','insitu.realm'),col=c('red','blue'),lty=1,cex=.7,bty='n')
  }else if(plots=='cladogenesis.BSM.rates'){
    #plot BSM in situ rates
    plot(c(1,1),xlim=c(0,length(insitu.hot.rates.BSM.CI)),ylim=c(0,max(c(unlist(insitu.hot.rates.BSM.CI),unlist(insitu.realm.rates.BSM.CI))[is.finite(c(unlist(insitu.hot.rates.BSM.CI),unlist(insitu.realm.rates.BSM.CI)))])),type='n',xaxt='n',xlab='age',ylab='insitu.rates',main=paste(name,' cladogenesis BSM rates',sep=''))
    lines(x=c(1:length(insitu.hot.rates.BSM.CI)),y=unlist(lapply(insitu.hot.rates.BSM.CI,function(x)x[3])),col=adjustcolor( "red", alpha.f = 0.9))
    polygon(x=c((1:length(insitu.hot.rates.BSM.CI)),(rev(1:length(insitu.hot.rates.BSM.CI)))),y=c((unlist(lapply(insitu.hot.rates.BSM.CI,function(x)x[1]))),rev(unlist(lapply(insitu.hot.rates.BSM.CI,function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
    lines(x=c(1:length(insitu.realm.rates.BSM.CI)),y=unlist(lapply(insitu.realm.rates.BSM.CI,function(x)x[3])),col=adjustcolor( "blue", alpha.f = 0.9))
    polygon(x=c((1:length(insitu.realm.rates.BSM.CI)),(rev(1:length(insitu.realm.rates.BSM.CI)))),y=c((unlist(lapply(insitu.realm.rates.BSM.CI,function(x)x[1]))),rev(unlist(lapply(insitu.realm.rates.BSM.CI,function(x)x[2])))),col=adjustcolor( "blue", alpha.f = 0.2))
    axis(1,at=c(0:length(insitu.realm.rates.BSM.CI)),labels=c((length(insitu.realm.rates.BSM.CI)+2):2))
    legend('topleft',c('insitu.hot','insitu.realm'),col=c('red','blue'),lty=1,cex=.7,bty='n')
  }else if(plots=='dispersal.events'){
    plot(c(1,1),xlim=c(0,length(hot.into.realm.time.bins.counts.rev.cumsum.CI)),ylim=c(0,max(c(unlist(hot.into.realm.time.bins.counts.rev.cumsum.CI),unlist(realm.into.hot.time.bins.counts.rev.cumsum.CI)))),type='n',xaxt='n',xlab='age',ylab='cumulative # of events',main=paste(name,' dispersal events',sep=''))
    lines(x=c(1:length(hot.into.realm.time.bins.counts.rev.cumsum.CI)),y=unlist(lapply(hot.into.realm.time.bins.counts.rev.cumsum.CI,function(x)x[3])),col=adjustcolor( "red", alpha.f = 0.9))
    polygon(x=c((1:length(hot.into.realm.time.bins.counts.rev.cumsum.CI)),(rev(1:length(hot.into.realm.time.bins.counts.rev.cumsum.CI)))),y=c((unlist(lapply(hot.into.realm.time.bins.counts.rev.cumsum.CI,function(x)x[1]))),rev(unlist(lapply(hot.into.realm.time.bins.counts.rev.cumsum.CI,function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
    lines(x=c(1:length(realm.into.hot.time.bins.counts.rev.cumsum.CI)),y=unlist(lapply(realm.into.hot.time.bins.counts.rev.cumsum.CI,function(x)x[3])),col=adjustcolor( "blue", alpha.f = 0.9))
    polygon(x=c((1:length(realm.into.hot.time.bins.counts.rev.cumsum.CI)),(rev(1:length(realm.into.hot.time.bins.counts.rev.cumsum.CI)))),y=c((unlist(lapply(realm.into.hot.time.bins.counts.rev.cumsum.CI,function(x)x[1]))),rev(unlist(lapply(realm.into.hot.time.bins.counts.rev.cumsum.CI,function(x)x[2])))),col=adjustcolor( "blue", alpha.f = 0.2))
    axis(1,at=c(0:length(BSMoutput.real.rates[['insitu.realm.rates.BSM.CI']])),labels=c((length(BSMoutput.real.rates[['insitu.realm.rates.BSM.CI']])+2):2))
    legend('topleft',c('hot.into.realm','realm.into.hot'),col=c('red','blue'),lty=1,cex=.7,bty='n')
  }else if(plots=='cladogenesis.events'){
    plot(c(1,1),xlim=c(0,length(insitu.hot.time.bins.counts.rev.cumsum.CI)),ylim=c(0,max(c(unlist(insitu.hot.time.bins.counts.rev.cumsum.CI),unlist(insitu.realm.time.bins.counts.rev.cumsum.CI)))),type='n',xaxt='n',xlab='age',ylab='cumulative # of events',main=paste(name,' cladogenetic events',sep=''))
    lines(x=c(1:length(insitu.hot.time.bins.counts.rev.cumsum.CI)),y=unlist(lapply(insitu.hot.time.bins.counts.rev.cumsum.CI,function(x)x[3])),col=adjustcolor( "red", alpha.f = 0.9))
    polygon(x=c((1:length(insitu.hot.time.bins.counts.rev.cumsum.CI)),(rev(1:length(insitu.hot.time.bins.counts.rev.cumsum.CI)))),y=c((unlist(lapply(insitu.hot.time.bins.counts.rev.cumsum.CI,function(x)x[1]))),rev(unlist(lapply(insitu.hot.time.bins.counts.rev.cumsum.CI,function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
    lines(x=c(1:length(insitu.realm.time.bins.counts.rev.cumsum.CI)),y=unlist(lapply(insitu.realm.time.bins.counts.rev.cumsum.CI,function(x)x[3])),col=adjustcolor( "blue", alpha.f = 0.9))
    polygon(x=c((1:length(insitu.realm.time.bins.counts.rev.cumsum.CI)),(rev(1:length(insitu.realm.time.bins.counts.rev.cumsum.CI)))),y=c((unlist(lapply(insitu.realm.time.bins.counts.rev.cumsum.CI,function(x)x[1]))),rev(unlist(lapply(insitu.realm.time.bins.counts.rev.cumsum.CI,function(x)x[2])))),col=adjustcolor( "blue", alpha.f = 0.2))
    axis(1,at=c(0:length(BSMoutput.real.rates[['insitu.realm.rates.BSM.CI']])),labels=c((length(BSMoutput.real.rates[['insitu.realm.rates.BSM.CI']])+2):2))
    legend('topleft',c('insitu.hot','insitu.realm'),col=c('red','blue'),lty=1,cex=.7,bty='n')
  }else if (plots=='dispersal.all.rates'){
    plot(c(1,1),xlim=c(0,length(dispersal.into.realm.rates.BSM.CI)),ylim=c(0,max(c(unlist(dispersal.into.realm.rates.BSM.CI),unlist(dispersal.into.hot.rates.BSM.CI))[is.finite(c(unlist(dispersal.into.realm.rates.BSM.CI),unlist(dispersal.into.hot.rates.BSM.CI)))])),type='n',xaxt='n',xlab='age',ylab='dispersal rates',main=paste(name,' dispersal rates',sep=''))
    lines(x=c(1:length(dispersal.into.realm.rates.BSM.CI)),y=unlist(lapply(dispersal.into.realm.rates.BSM.CI,function(x)x[3])),col=adjustcolor( "red", alpha.f = 0.9))
    polygon(x=c((1:length(dispersal.into.realm.rates.BSM.CI)),(rev(1:length(dispersal.into.realm.rates.BSM.CI)))),y=c((unlist(lapply(dispersal.into.realm.rates.BSM.CI,function(x)x[1]))),rev(unlist(lapply(dispersal.into.realm.rates.BSM.CI,function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
    lines(x=c(1:length(dispersal.into.hot.rates.BSM.CI)),y=unlist(lapply(dispersal.into.hot.rates.BSM.CI,function(x)x[3])),col=adjustcolor( "blue", alpha.f = 0.9))
    polygon(x=c((1:length(dispersal.into.hot.rates.BSM.CI)),(rev(1:length(dispersal.into.hot.rates.BSM.CI)))),y=c((unlist(lapply(dispersal.into.hot.rates.BSM.CI,function(x)x[1]))),rev(unlist(lapply(dispersal.into.hot.rates.BSM.CI,function(x)x[2])))),col=adjustcolor( "blue", alpha.f = 0.2))
    axis(1,at=c(0:length(dispersal.into.hot.rates.BSM.CI)),labels=c((length(dispersal.into.hot.rates.BSM.CI)+2):2))
    legend('topleft',c('dispersal.into.hot','dispersal.into.realm'),col=c('blue','red'),lty=1,cex=.7,bty='n')
  }else if(plots=='dispersal.all.events'){
    plot(c(1,1),xlim=c(0,length(dispersal.into.realm.time.bins.counts.rev.cumsum.CI)),ylim=c(0,max(c(unlist(dispersal.into.realm.time.bins.counts.rev.cumsum.CI),unlist(dispersal.into.hot.time.bins.counts.rev.cumsum.CI))[is.finite(c(unlist(dispersal.into.realm.time.bins.counts.rev.cumsum.CI),unlist(dispersal.into.hot.time.bins.counts.rev.cumsum.CI)))])),type='n',xaxt='n',xlab='age',ylab='cumulative # of events',main=paste(name,' dispersal events',sep=''))
    lines(x=c(1:length(dispersal.into.realm.time.bins.counts.rev.cumsum.CI)),y=unlist(lapply(dispersal.into.realm.time.bins.counts.rev.cumsum.CI,function(x)x[3])),col=adjustcolor( "red", alpha.f = 0.9))
    polygon(x=c((1:length(dispersal.into.realm.time.bins.counts.rev.cumsum.CI)),(rev(1:length(dispersal.into.realm.time.bins.counts.rev.cumsum.CI)))),y=c((unlist(lapply(dispersal.into.realm.time.bins.counts.rev.cumsum.CI,function(x)x[1]))),rev(unlist(lapply(dispersal.into.realm.time.bins.counts.rev.cumsum.CI,function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
    lines(x=c(1:length(dispersal.into.hot.time.bins.counts.rev.cumsum.CI)),y=unlist(lapply(dispersal.into.hot.time.bins.counts.rev.cumsum.CI,function(x)x[3])),col=adjustcolor( "blue", alpha.f = 0.9))
    polygon(x=c((1:length(dispersal.into.hot.time.bins.counts.rev.cumsum.CI)),(rev(1:length(dispersal.into.hot.time.bins.counts.rev.cumsum.CI)))),y=c((unlist(lapply(dispersal.into.hot.time.bins.counts.rev.cumsum.CI,function(x)x[1]))),rev(unlist(lapply(dispersal.into.hot.time.bins.counts.rev.cumsum.CI,function(x)x[2])))),col=adjustcolor( "blue", alpha.f = 0.2))
    axis(1,at=c(0:length(BSMoutput.real.rates[['insitu.realm.rates.BSM.CI']])),labels=c((length(BSMoutput.real.rates[['insitu.realm.rates.BSM.CI']])+2):2))
    legend('topleft',c('dispersal.into.realm','dispersal.into.hot'),col=c('red','blue'),lty=1,cex=.7,bty='n')
  }else if(plots=='dispersal.into.from.events'){
    plot(c(1,1),xlim=c(0,length(dispersal.into.realm.time.bins.counts.rev.cumsum.CI)),ylim=c(0,max(c(unlist(dispersal.into.realm.time.bins.counts.rev.cumsum.CI),unlist(dispersal.into.hot.time.bins.counts.rev.cumsum.CI),unlist(dispersal.from.hot.time.bins.counts.rev.cumsum.CI),unlist(dispersal.from.realm.time.bins.counts.rev.cumsum.CI)))),type='n',xaxt='n',xlab='age',ylab='cumulative # of events',main=paste(name,' dispersal events',sep=''))
    lines(x=c(1:length(dispersal.into.realm.time.bins.counts.rev.cumsum.CI)),y=unlist(lapply(dispersal.into.realm.time.bins.counts.rev.cumsum.CI,function(x)x[3])),col=adjustcolor( "green", alpha.f = 0.9))
    polygon(x=c((1:length(dispersal.into.realm.time.bins.counts.rev.cumsum.CI)),(rev(1:length(dispersal.into.realm.time.bins.counts.rev.cumsum.CI)))),y=c((unlist(lapply(dispersal.into.realm.time.bins.counts.rev.cumsum.CI,function(x)x[1]))),rev(unlist(lapply(dispersal.into.realm.time.bins.counts.rev.cumsum.CI,function(x)x[2])))),col=adjustcolor( "green", alpha.f = 0.2))
    lines(x=c(1:length(dispersal.into.hot.time.bins.counts.rev.cumsum.CI)),y=unlist(lapply(dispersal.into.hot.time.bins.counts.rev.cumsum.CI,function(x)x[3])),col=adjustcolor( "orange", alpha.f = 0.9))
    polygon(x=c((1:length(dispersal.into.hot.time.bins.counts.rev.cumsum.CI)),(rev(1:length(dispersal.into.hot.time.bins.counts.rev.cumsum.CI)))),y=c((unlist(lapply(dispersal.into.hot.time.bins.counts.rev.cumsum.CI,function(x)x[1]))),rev(unlist(lapply(dispersal.into.hot.time.bins.counts.rev.cumsum.CI,function(x)x[2])))),col=adjustcolor( "orange", alpha.f = 0.2))
    lines(x=c(1:length(dispersal.from.realm.time.bins.counts.rev.cumsum.CI)),y=unlist(lapply(dispersal.from.realm.time.bins.counts.rev.cumsum.CI,function(x)x[3])),col=adjustcolor( "blue", alpha.f = 0.9))
    polygon(x=c((1:length(dispersal.from.realm.time.bins.counts.rev.cumsum.CI)),(rev(1:length(dispersal.from.realm.time.bins.counts.rev.cumsum.CI)))),y=c((unlist(lapply(dispersal.from.realm.time.bins.counts.rev.cumsum.CI,function(x)x[1]))),rev(unlist(lapply(dispersal.from.realm.time.bins.counts.rev.cumsum.CI,function(x)x[2])))),col=adjustcolor( "blue", alpha.f = 0.2))
    lines(x=c(1:length(dispersal.from.hot.time.bins.counts.rev.cumsum.CI)),y=unlist(lapply(dispersal.from.hot.time.bins.counts.rev.cumsum.CI,function(x)x[3])),col=adjustcolor( "red", alpha.f = 0.9))
    polygon(x=c((1:length(dispersal.from.hot.time.bins.counts.rev.cumsum.CI)),(rev(1:length(dispersal.from.hot.time.bins.counts.rev.cumsum.CI)))),y=c((unlist(lapply(dispersal.from.hot.time.bins.counts.rev.cumsum.CI,function(x)x[1]))),rev(unlist(lapply(dispersal.from.hot.time.bins.counts.rev.cumsum.CI,function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
    axis(1,at=c(0:length(BSMoutput.real.rates[['insitu.realm.rates.BSM.CI']])),labels=c((length(BSMoutput.real.rates[['insitu.realm.rates.BSM.CI']])+2):2))
    legend('topleft',c('dispersal.from.realm','dispersal.from.hot','dispersal.into.realm','dispersal.into.hot'),col=c('blue','red','green','orange'),lty=1,cex=.7,bty='n')
  }
  
  #plot BSM dispersal rates hot into realm
  
  ##plot BAMM in situ rates
  #par(new=T)
  #plot(c(1,1),xlim=c(0,length(insitu.hot.rates.BAMM.CI)),ylim=c(0,max(c(unlist(insitu.hot.rates.BAMM.CI),unlist(insitu.realm.rates.BAMM.CI)))),type='n',xaxt='n',xlab='age',yaxt='n')
  ##polygon(x=c(c(age:length(insitu.hot.rates.BAMM.CI)),c(length(insitu.hot.rates.BAMM.CI):age)),y=c(unlist(lapply(insitu.hot.rates.BAMM.CI[c(age:length(insitu.hot.rates.BAMM.CI))],function(x)x[1])),rev(unlist(lapply(insitu.hot.rates.BAMM.CI[c(age:length(insitu.hot.rates.BAMM.CI))],function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
  #lines(x=c(1:length(insitu.hot.rates.BAMM.CI)),y=unlist(lapply(insitu.hot.rates.BAMM.CI,function(x)x[3])),col=adjustcolor( "red", alpha.f = 0.9))
  #polygon(x=c((1:length(insitu.hot.rates.BAMM.CI)),(rev(1:length(insitu.hot.rates.BAMM.CI)))),y=c((unlist(lapply(insitu.hot.rates.BAMM.CI,function(x)x[1]))),rev(unlist(lapply(insitu.hot.rates.BAMM.CI,function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
  #lines(x=c(1:length(insitu.realm.rates.BAMM.CI)),y=unlist(lapply(insitu.realm.rates.BAMM.CI,function(x)x[3])),col=adjustcolor( "blue", alpha.f = 0.9))
  #polygon(x=c((1:length(insitu.realm.rates.BAMM.CI)),(rev(1:length(insitu.realm.rates.BAMM.CI)))),y=c((unlist(lapply(insitu.realm.rates.BAMM.CI,function(x)x[1]))),rev(unlist(lapply(insitu.realm.rates.BAMM.CI,function(x)x[2])))),col=adjustcolor( "blue", alpha.f = 0.2))
  ##axis(1,at=c(0:length(BSMoutput.real.rates[['insitu.realm.rates.BSM.CI']])),labels=c((length(BSMoutput.real.rates[['insitu.realm.rates.BSM.CI']])+2):2))
  #axis(4)
  #legend('topright',c('insitu.hot','insitu.realm'),col=c('red','blue'),lty=1,cex=.7,bty='n')
}
plot_rates_events_real_2mybin_difference_both<-function(BSMoutput.real.rates.mammals,BSMoutput.real.events.mammals,BSMoutput.real.rates.birds,BSMoutput.real.events.birds,age,name,plots){
  #correct age
  age.mammals<-(length(BSMoutput.real.rates.mammals[['insitu.realm.rates.BSM.CI']])*2-age)/2
  #trim lists to age period
  BSMoutput.real.rates.mammals<-lapply(BSMoutput.real.rates.mammals,function(x)x[age.mammals:length(x)])
  BSMoutput.real.events.mammals<-lapply(BSMoutput.real.events.mammals,function(x)x[age.mammals:length(x)])
  #delete Inf in rates and events
  BSMoutput.real.rates.mammals<-lapply(BSMoutput.real.rates.mammals,function(x) lapply(x,function(x) {x[!is.finite(x)]<-NA;return(x)}))
  BSMoutput.real.events.mammals<-lapply(BSMoutput.real.events.mammals,function(x) lapply(x,function(x) {x[!is.finite(x)]<-NA;return(x)}))
  #remove the last bin - two million year time bin
  BSMoutput.real.rates.mammals<-lapply(BSMoutput.real.rates.mammals,function(x) x<-x[-c(length(x))])
  BSMoutput.real.events.mammals<-lapply(BSMoutput.real.events.mammals,function(x) x<-x[-c((length(x)-1),length(x))])
  #assign variables for easier handling
  hot.into.realm.rates.BSM.mammals<-BSMoutput.real.rates.mammals[['hot.into.realm.rates.BSM.CI']]
  realm.into.hot.rates.BSM.mammals<-BSMoutput.real.rates.mammals[['realm.into.hot.rates.BSM.CI']]
  insitu.realm.rates.BAMM.mammals<-BSMoutput.real.rates.mammals[['insitu.realm.rates.BAMM.CI']]
  insitu.hot.rates.BAMM.mammals<-BSMoutput.real.rates.mammals[['insitu.hot.rates.BAMM.CI']]
  insitu.realm.rates.BSM.mammals<-BSMoutput.real.rates.mammals[['insitu.realm.rates.BSM.CI']]
  insitu.hot.rates.BSM.mammals<-BSMoutput.real.rates.mammals[['insitu.hot.rates.BSM.CI']]
  dispersal.into.realm.rates.BSM.mammals<-BSMoutput.real.rates.mammals[['dispersal.into.realm.rates.BSM.CI']]
  dispersal.into.hot.rates.BSM.mammals<-BSMoutput.real.rates.mammals[['dispersal.into.hot.rates.BSM.CI']]
  dispersal.from.realm.rates.BSM.mammals<-BSMoutput.real.rates.mammals[['dispersal.from.realm.rates.BSM.CI']]
  dispersal.from.hot.rates.BSM.mammals<-BSMoutput.real.rates.mammals[['dispersal.from.hot.rates.BSM.CI']]
  
  hot.into.realm.time.bins.counts.rev.cumsum.mammals<-BSMoutput.real.events.mammals[['hot.into.realm.time.bins.counts.rev.cumsum.CI']]
  realm.into.hot.time.bins.counts.rev.cumsum.mammals<-BSMoutput.real.events.mammals[['realm.into.hot.time.bins.counts.rev.cumsum.CI']]
  insitu.realm.time.bins.counts.rev.cumsum.mammals<-BSMoutput.real.events.mammals[['insitu.realm.time.bins.counts.rev.cumsum.CI']]
  insitu.hot.time.bins.counts.rev.cumsum.mammals<-BSMoutput.real.events.mammals[['insitu.hot.time.bins.counts.rev.cumsum.CI']]
  dispersal.into.realm.time.bins.counts.rev.cumsum.mammals<-BSMoutput.real.events.mammals[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']]
  dispersal.into.hot.time.bins.counts.rev.cumsum.mammals<-BSMoutput.real.events.mammals[['dispersal.into.hot.time.bins.counts.rev.cumsum.CI']]
  dispersal.from.realm.time.bins.counts.rev.cumsum.mammals<-BSMoutput.real.events.mammals[['dispersal.from.realm.time.bins.counts.rev.cumsum.CI']]
  dispersal.from.hot.time.bins.counts.rev.cumsum.mammals<-BSMoutput.real.events.mammals[['dispersal.from.hot.time.bins.counts.rev.cumsum.CI']]
  
  age.birds<-(length(BSMoutput.real.rates.birds[['insitu.realm.rates.BSM.CI']])*2-age)/2
  #trim lists to age period
  BSMoutput.real.rates.birds<-lapply(BSMoutput.real.rates.birds,function(x)x[age.birds:length(x)])
  BSMoutput.real.events.birds<-lapply(BSMoutput.real.events.birds,function(x)x[age.birds:length(x)])
  #delete Inf in rates and events
  BSMoutput.real.rates.birds<-lapply(BSMoutput.real.rates.birds,function(x) lapply(x,function(x) {x[!is.finite(x)]<-NA;return(x)}))
  BSMoutput.real.events.birds<-lapply(BSMoutput.real.events.birds,function(x) lapply(x,function(x) {x[!is.finite(x)]<-NA;return(x)}))
  #remove the last bin - two million year time bin
  BSMoutput.real.rates.birds<-lapply(BSMoutput.real.rates.birds,function(x) x<-x[-c(length(x))])
  BSMoutput.real.events.birds<-lapply(BSMoutput.real.events.birds,function(x) x<-x[-c((length(x)-1),length(x))])
  #assign variables for easier handling
  hot.into.realm.rates.BSM.birds<-BSMoutput.real.rates.birds[['hot.into.realm.rates.BSM.CI']]
  realm.into.hot.rates.BSM.birds<-BSMoutput.real.rates.birds[['realm.into.hot.rates.BSM.CI']]
  insitu.realm.rates.BAMM.birds<-BSMoutput.real.rates.birds[['insitu.realm.rates.BAMM.CI']]
  insitu.hot.rates.BAMM.birds<-BSMoutput.real.rates.birds[['insitu.hot.rates.BAMM.CI']]
  insitu.realm.rates.BSM.birds<-BSMoutput.real.rates.birds[['insitu.realm.rates.BSM.CI']]
  insitu.hot.rates.BSM.birds<-BSMoutput.real.rates.birds[['insitu.hot.rates.BSM.CI']]
  dispersal.into.realm.rates.BSM.birds<-BSMoutput.real.rates.birds[['dispersal.into.realm.rates.BSM.CI']]
  dispersal.into.hot.rates.BSM.birds<-BSMoutput.real.rates.birds[['dispersal.into.hot.rates.BSM.CI']]
  dispersal.from.realm.rates.BSM.birds<-BSMoutput.real.rates.birds[['dispersal.from.realm.rates.BSM.CI']]
  dispersal.from.hot.rates.BSM.birds<-BSMoutput.real.rates.birds[['dispersal.from.hot.rates.BSM.CI']]
  
  hot.into.realm.time.bins.counts.rev.cumsum.birds<-BSMoutput.real.events.birds[['hot.into.realm.time.bins.counts.rev.cumsum.CI']]
  realm.into.hot.time.bins.counts.rev.cumsum.birds<-BSMoutput.real.events.birds[['realm.into.hot.time.bins.counts.rev.cumsum.CI']]
  insitu.realm.time.bins.counts.rev.cumsum.birds<-BSMoutput.real.events.birds[['insitu.realm.time.bins.counts.rev.cumsum.CI']]
  insitu.hot.time.bins.counts.rev.cumsum.birds<-BSMoutput.real.events.birds[['insitu.hot.time.bins.counts.rev.cumsum.CI']]
  dispersal.into.realm.time.bins.counts.rev.cumsum.birds<-BSMoutput.real.events.birds[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']]
  dispersal.into.hot.time.bins.counts.rev.cumsum.birds<-BSMoutput.real.events.birds[['dispersal.into.hot.time.bins.counts.rev.cumsum.CI']]
  dispersal.from.realm.time.bins.counts.rev.cumsum.birds<-BSMoutput.real.events.birds[['dispersal.from.realm.time.bins.counts.rev.cumsum.CI']]
  dispersal.from.hot.time.bins.counts.rev.cumsum.birds<-BSMoutput.real.events.birds[['dispersal.from.hot.time.bins.counts.rev.cumsum.CI']]
  
  if(plots=='within.realm.dispersal.rates'){
    #difference.into.hotminusrealm<-lapply(c(1:length(dispersal.into.realm.rates.BSM)),function(x)sample(dispersal.into.hot.rates.BSM[[x]],size = 100,replace=TRUE)-sample(dispersal.into.realm.rates.BSM[[x]],size = 100,replace=TRUE))
    difference.hotintorealm.realmintohot.mammals<-lapply(c(1:length(hot.into.realm.rates.BSM.mammals)),function(x)sample(hot.into.realm.rates.BSM.mammals[[x]],size = 50,replace=FALSE)-sample(realm.into.hot.rates.BSM.mammals[[x]],size = 50,replace=FALSE))
    difference.hotintorealm.realmintohot.birds<-lapply(c(1:length(hot.into.realm.rates.BSM.birds)),function(x)sample(hot.into.realm.rates.BSM.birds[[x]],size = 50,replace=FALSE)-sample(realm.into.hot.rates.BSM.birds[[x]],size = 50,replace=FALSE))
    median.mammals<-unlist(lapply(difference.hotintorealm.realmintohot.mammals,function(x)median(x,na.rm=TRUE)))
    lower.mammals<-unlist(lapply(difference.hotintorealm.realmintohot.mammals,function(x)quantile(x,c(0.1),na.rm=TRUE)))
    upper.mammals<-unlist(lapply(difference.hotintorealm.realmintohot.mammals,function(x)quantile(x,c(0.9),na.rm=TRUE)))
    median.birds<-unlist(lapply(difference.hotintorealm.realmintohot.birds,function(x)median(x,na.rm=TRUE)))
    lower.birds<-unlist(lapply(difference.hotintorealm.realmintohot.birds,function(x)quantile(x,c(0.1),na.rm=TRUE)))
    upper.birds<-unlist(lapply(difference.hotintorealm.realmintohot.birds,function(x)quantile(x,c(0.9),na.rm=TRUE)))
    #lower<-unlist(lapply(difference.hotintorealm.realmintohot,function(x)quantile(x,c(0.025),na.rm=TRUE)))
    #upper<-unlist(lapply(difference.hotintorealm.realmintohot,function(x)quantile(x,c(0.975),na.rm=TRUE)))
    plot(x=c(1:length(median.mammals)),y=median.mammals,ylim=c(-max(abs(c(median.mammals,upper.mammals,lower.mammals,median.birds,upper.birds,lower.birds))),max(abs(c(median.mammals,upper.mammals,lower.mammals,median.birds,upper.birds,lower.birds)))),type='n',xaxt='n',yaxt='n',xlab='age',ylab='diff. dispersal hot.to.realm - realm.to.hot',main=paste(name,' diff in dispersal rates',sep=''))
    abline(h=0,lty=2)
    polygon(x=c(c(1:length(median.mammals),length(median.mammals):1)),y=c(lower.mammals,rev(upper.mammals)),col=adjustcolor(col='#009999',alpha.f=0.25))
    lines(x=c(1:length(median.mammals)),y=median.mammals,col='#009999') 
    polygon(x=c(c(1:length(median.birds),length(median.birds):1)),y=c(lower.birds,rev(upper.birds)),col=adjustcolor(col="#FF7400",alpha.f=0.25))
    lines(x=c(1:length(median.birds)),y=median.birds,col="#FF7400") 
    axis(1,at=c(1:length(difference.hotintorealm.realmintohot.mammals)),labels=c(seq(from=length(difference.hotintorealm.realmintohot.mammals)*2,to=2,by=-2)))
    axis(2,las=2)
    #plot(c(1,1),xlim=c(0,length(hot.into.realm.rates.BSM.CI)),ylim=c(0,max(c(unlist(hot.into.realm.rates.BSM.CI),unlist(realm.into.hot.rates.BSM.CI))[is.finite(c(unlist(hot.into.realm.rates.BSM.CI),unlist(realm.into.hot.rates.BSM.CI)))])),type='n',xaxt='n',xlab='age',ylab='dispersal rates',main=paste(name,' dispersal rates',sep=''))
    #lines(x=c(1:length(hot.into.realm.rates.BSM.CI)),y=unlist(lapply(hot.into.realm.rates.BSM.CI,function(x)x[3])),col=adjustcolor( "red", alpha.f = 0.9))
    #polygon(x=c((1:length(hot.into.realm.rates.BSM.CI)),(rev(1:length(hot.into.realm.rates.BSM.CI)))),y=c((unlist(lapply(hot.into.realm.rates.BSM.CI,function(x)x[1]))),rev(unlist(lapply(hot.into.realm.rates.BSM.CI,function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
    #lines(x=c(1:length(realm.into.hot.rates.BSM.CI)),y=unlist(lapply(realm.into.hot.rates.BSM.CI,function(x)x[3])),col=adjustcolor( "blue", alpha.f = 0.9))
    #polygon(x=c((1:length(realm.into.hot.rates.BSM.CI)),(rev(1:length(realm.into.hot.rates.BSM.CI)))),y=c((unlist(lapply(realm.into.hot.rates.BSM.CI,function(x)x[1]))),rev(unlist(lapply(realm.into.hot.rates.BSM.CI,function(x)x[2])))),col=adjustcolor( "blue", alpha.f = 0.2))
    #
    legend('topleft',c('mammal.diff','bird.diff'),col=c('#009999',"#FF7400"),lty=1,cex=.7,bty='n')
  }else if(plots=='within.realm.cladogenesis.BAMM.rates'){
    difference.insituhot.insiturealm.mammals<-lapply(c(1:length(insitu.hot.rates.BAMM.mammals)),function(x)sample(insitu.hot.rates.BAMM.mammals[[x]],size = 100,replace=TRUE)-sample(insitu.realm.rates.BAMM.mammals[[x]],size = 100,replace=TRUE))
    difference.insituhot.insiturealm.birds<-lapply(c(1:length(insitu.hot.rates.BAMM.birds)),function(x)sample(insitu.hot.rates.BAMM.birds[[x]],size = 100,replace=TRUE)-sample(insitu.realm.rates.BAMM.birds[[x]],size = 100,replace=TRUE))
    median.mammals<-unlist(lapply(difference.insituhot.insiturealm.mammals,function(x)median(x,na.rm=TRUE)))
    lower.mammals<-unlist(lapply(difference.insituhot.insiturealm.mammals,function(x)quantile(x,c(0.1),na.rm=TRUE)))
    upper.mammals<-unlist(lapply(difference.insituhot.insiturealm.mammals,function(x)quantile(x,c(0.9),na.rm=TRUE)))
    median.birds<-unlist(lapply(difference.insituhot.insiturealm.birds,function(x)median(x,na.rm=TRUE)))
    lower.birds<-unlist(lapply(difference.insituhot.insiturealm.birds,function(x)quantile(x,c(0.1),na.rm=TRUE)))
    upper.birds<-unlist(lapply(difference.insituhot.insiturealm.birds,function(x)quantile(x,c(0.9),na.rm=TRUE)))
    #lower<-unlist(lapply(difference.insituhot.insiturealm,function(x)quantile(x,c(0.025),na.rm=TRUE)))
    #upper<-unlist(lapply(difference.insituhot.insiturealm,function(x)quantile(x,c(0.975),na.rm=TRUE)))
    plot(x=c(1:length(median.mammals)),y=median.mammals,ylim=c(-max(abs(c(median.mammals,upper.mammals,lower.mammals,median.birds,upper.birds,lower.birds))),max(abs(c(median.mammals,upper.mammals,lower.mammals,median.birds,upper.birds,lower.birds)))),type='n',xaxt='n',yaxt='n',xlab='age',ylab='diff. insitu.hot - insitu.realm BAMM',main=paste(name,' diff in insitu.BAMM.rates',sep=''))
    abline(h=0,lty=2)
    polygon(x=c(c(1:length(median.mammals),length(median.mammals):1)),y=c(lower.mammals,rev(upper.mammals)),col=adjustcolor(col='#009999',alpha.f=0.25))
    lines(x=c(1:length(median.mammals)),y=median.mammals,col='#009999') 
    polygon(x=c(c(1:length(median.birds),length(median.birds):1)),y=c(lower.birds,rev(upper.birds)),col=adjustcolor(col="#FF7400",alpha.f=0.25))
    lines(x=c(1:length(median.birds)),y=median.birds,col="#FF7400") 
    axis(1,at=c(1:length(difference.insituhot.insiturealm.mammals)),labels=c(seq(from=length(difference.insituhot.insiturealm.mammals)*2,to=2,by=-2)))
    axis(2,las=2)
    #plot(c(1,1),xlim=c(0,length(insitu.hot.rates.BAMM.CI)),ylim=c(0,max(c(unlist(insitu.hot.rates.BAMM.CI),unlist(insitu.realm.rates.BAMM.CI))[is.finite(c(unlist(insitu.hot.rates.BAMM.CI),unlist(insitu.realm.rates.BAMM.CI)))])),type='n',xaxt='n',xlab='age',ylab='dispersal rates',main=paste(name,' dispersal rates',sep=''))
    #lines(x=c(1:length(insitu.hot.rates.BAMM.CI)),y=unlist(lapply(insitu.hot.rates.BAMM.CI,function(x)x[3])),col=adjustcolor( "red", alpha.f = 0.9))
    #polygon(x=c((1:length(insitu.hot.rates.BAMM.CI)),(rev(1:length(insitu.hot.rates.BAMM.CI)))),y=c((unlist(lapply(insitu.hot.rates.BAMM.CI,function(x)x[1]))),rev(unlist(lapply(insitu.hot.rates.BAMM.CI,function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
    #lines(x=c(1:length(insitu.realm.rates.BAMM.CI)),y=unlist(lapply(insitu.realm.rates.BAMM.CI,function(x)x[3])),col=adjustcolor( "blue", alpha.f = 0.9))
    #polygon(x=c((1:length(insitu.realm.rates.BAMM.CI)),(rev(1:length(insitu.realm.rates.BAMM.CI)))),y=c((unlist(lapply(insitu.realm.rates.BAMM.CI,function(x)x[1]))),rev(unlist(lapply(insitu.realm.rates.BAMM.CI,function(x)x[2])))),col=adjustcolor( "blue", alpha.f = 0.2))
    #
    legend('topleft',c('mammal.diff','bird.diff'),col=c('#009999',"#FF7400"),lty=1,cex=.7,bty='n')
  }else if(plots=='within.realm.cladogenesis.BSM.rates'){
    difference.insituhot.insiturealm.mammals<-lapply(c(1:length(insitu.hot.rates.BSM.mammals)),function(x)sample(insitu.hot.rates.BSM.mammals[[x]],size = 100,replace=TRUE)-sample(insitu.realm.rates.BSM.mammals[[x]],size = 100,replace=TRUE))
    difference.insituhot.insiturealm.birds<-lapply(c(1:length(insitu.hot.rates.BSM.birds)),function(x)sample(insitu.hot.rates.BSM.birds[[x]],size = 100,replace=TRUE)-sample(insitu.realm.rates.BSM.birds[[x]],size = 100,replace=TRUE))
    median.mammals<-unlist(lapply(difference.insituhot.insiturealm.mammals,function(x)median(x,na.rm=TRUE)))
    lower.mammals<-unlist(lapply(difference.insituhot.insiturealm.mammals,function(x)quantile(x,c(0.1),na.rm=TRUE)))
    upper.mammals<-unlist(lapply(difference.insituhot.insiturealm.mammals,function(x)quantile(x,c(0.9),na.rm=TRUE)))
    median.birds<-unlist(lapply(difference.insituhot.insiturealm.birds,function(x)median(x,na.rm=TRUE)))
    lower.birds<-unlist(lapply(difference.insituhot.insiturealm.birds,function(x)quantile(x,c(0.1),na.rm=TRUE)))
    upper.birds<-unlist(lapply(difference.insituhot.insiturealm.birds,function(x)quantile(x,c(0.9),na.rm=TRUE)))
    #lower<-unlist(lapply(difference.insituhot.insiturealm,function(x)quantile(x,c(0.025),na.rm=TRUE)))
    #upper<-unlist(lapply(difference.insituhot.insiturealm,function(x)quantile(x,c(0.975),na.rm=TRUE)))
    plot(x=c(1:length(median.mammals)),y=median.mammals,ylim=c(-max(abs(c(median.mammals,upper.mammals,lower.mammals,median.birds,upper.birds,lower.birds))),max(abs(c(median.mammals,upper.mammals,lower.mammals,median.birds,upper.birds,lower.birds)))),type='n',xaxt='n',yaxt='n',xlab='age',ylab='diff. insitu.hot - insitu.realm BSM',main=paste(name,' diff in insitu.BSM.rates',sep=''))
    abline(h=0,lty=2)
    polygon(x=c(c(1:length(median.mammals),length(median.mammals):1)),y=c(lower.mammals,rev(upper.mammals)),col=adjustcolor(col='#009999',alpha.f=0.25))
    lines(x=c(1:length(median.mammals)),y=median.mammals,col='#009999') 
    polygon(x=c(c(1:length(median.birds),length(median.birds):1)),y=c(lower.birds,rev(upper.birds)),col=adjustcolor(col="#FF7400",alpha.f=0.25))
    lines(x=c(1:length(median.birds)),y=median.birds,col="#FF7400") 
    axis(1,at=c(1:length(difference.insituhot.insiturealm.mammals)),labels=c(seq(from=length(difference.insituhot.insiturealm.mammals)*2,to=2,by=-2)))
    axis(2,las=2)
    #plot(c(1,1),xlim=c(0,length(insitu.hot.rates.BSM.CI)),ylim=c(0,max(c(unlist(insitu.hot.rates.BSM.CI),unlist(insitu.realm.rates.BSM.CI))[is.finite(c(unlist(insitu.hot.rates.BSM.CI),unlist(insitu.realm.rates.BSM.CI)))])),type='n',xaxt='n',xlab='age',ylab='dispersal rates',main=paste(name,' dispersal rates',sep=''))
    #lines(x=c(1:length(insitu.hot.rates.BSM.CI)),y=unlist(lapply(insitu.hot.rates.BSM.CI,function(x)x[3])),col=adjustcolor( "red", alpha.f = 0.9))
    #polygon(x=c((1:length(insitu.hot.rates.BSM.CI)),(rev(1:length(insitu.hot.rates.BSM.CI)))),y=c((unlist(lapply(insitu.hot.rates.BSM.CI,function(x)x[1]))),rev(unlist(lapply(insitu.hot.rates.BSM.CI,function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
    #lines(x=c(1:length(insitu.realm.rates.BSM.CI)),y=unlist(lapply(insitu.realm.rates.BSM.CI,function(x)x[3])),col=adjustcolor( "blue", alpha.f = 0.9))
    #polygon(x=c((1:length(insitu.realm.rates.BSM.CI)),(rev(1:length(insitu.realm.rates.BSM.CI)))),y=c((unlist(lapply(insitu.realm.rates.BSM.CI,function(x)x[1]))),rev(unlist(lapply(insitu.realm.rates.BSM.CI,function(x)x[2])))),col=adjustcolor( "blue", alpha.f = 0.2))
    #
    legend('topleft',c('mammal.diff','bird.diff'),col=c('#009999',"#FF7400"),lty=1,cex=.7,bty='n')
  }else if(plots=='global.dispersal.rates'){
    #difference.into.hotminusrealm<-lapply(c(1:length(dispersal.into.realm.rates.BSM)),function(x)sample(dispersal.into.hot.rates.BSM[[x]],size = 100,replace=TRUE)-sample(dispersal.into.realm.rates.BSM[[x]],size = 100,replace=TRUE))
    difference.hotrealm.from<-lapply(c(1:length(dispersal.from.realm.rates.BSM)),function(x)sample(dispersal.from.hot.rates.BSM[[x]],size = 100,replace=TRUE)-sample(dispersal.from.realm.rates.BSM[[x]],size = 100,replace=TRUE))
    difference.hotrealm.into<-lapply(c(1:length(dispersal.into.realm.rates.BSM)),function(x)sample(dispersal.into.hot.rates.BSM[[x]],size = 100,replace=TRUE)-sample(dispersal.into.realm.rates.BSM[[x]],size = 100,replace=TRUE))
    median.from<-unlist(lapply(difference.hotrealm.from,function(x)median(x,na.rm=TRUE)))
    lower.from<-unlist(lapply(difference.hotrealm.from,function(x)quantile(x,c(0.1),na.rm=TRUE)))
    upper.from<-unlist(lapply(difference.hotrealm.from,function(x)quantile(x,c(0.9),na.rm=TRUE)))
    median.into<-unlist(lapply(difference.hotrealm.into,function(x)median(x,na.rm=TRUE)))
    lower.into<-unlist(lapply(difference.hotrealm.into,function(x)quantile(x,c(0.1),na.rm=TRUE)))
    upper.into<-unlist(lapply(difference.hotrealm.into,function(x)quantile(x,c(0.9),na.rm=TRUE)))
    
    plot(x=c(1:length(median.from)),y=median.from,ylim=c(-max(abs(c(median.from,lower.from,upper.from,median.into,lower.into,upper.into))),max(abs(c(median.from,lower.from,upper.from,median.into,lower.into,upper.into)))),type='n',xaxt='n',xlab='age',ylab='diff. dispersal hot - realm',main=paste(name,' diff in global dispersal rates',sep=''))
    polygon(x=c(c(1:length(median.from),length(median.from):1)),y=c(lower.from,rev(upper.from)),col=adjustcolor('red',alpha.f = 0.25))
    lines(x=c(1:length(median.from)),y=median.from,col='red') 
    abline(h=0,lty=2)
    polygon(x=c(c(1:length(median.into),length(median.into):1)),y=c(lower.into,rev(upper.into)),col=adjustcolor('blue',alpha.f = 0.25))
    lines(x=c(1:length(median.into)),y=median.into,col='blue') 
    axis(1,at=c(1:length(difference.hotintorealm.realmintohot)),labels=c(seq(from=length(difference.hotintorealm.realmintohot)*2,to=2,by=-2)))
    legend('topleft',c('dispersal.from.diff','dispersal.into.diff'),col=c('red','blue'),lty=1,cex=.7,bty='n')
    #plot(c(1,1),xlim=c(0,length(hot.into.realm.rates.BSM.CI)),ylim=c(0,max(c(unlist(hot.into.realm.rates.BSM.CI),unlist(realm.into.hot.rates.BSM.CI))[is.finite(c(unlist(hot.into.realm.rates.BSM.CI),unlist(realm.into.hot.rates.BSM.CI)))])),type='n',xaxt='n',xlab='age',ylab='dispersal rates',main=paste(name,' dispersal rates',sep=''))
    #lines(x=c(1:length(hot.into.realm.rates.BSM.CI)),y=unlist(lapply(hot.into.realm.rates.BSM.CI,function(x)x[3])),col=adjustcolor( "red", alpha.f = 0.9))
    #polygon(x=c((1:length(hot.into.realm.rates.BSM.CI)),(rev(1:length(hot.into.realm.rates.BSM.CI)))),y=c((unlist(lapply(hot.into.realm.rates.BSM.CI,function(x)x[1]))),rev(unlist(lapply(hot.into.realm.rates.BSM.CI,function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
    #lines(x=c(1:length(realm.into.hot.rates.BSM.CI)),y=unlist(lapply(realm.into.hot.rates.BSM.CI,function(x)x[3])),col=adjustcolor( "blue", alpha.f = 0.9))
    #polygon(x=c((1:length(realm.into.hot.rates.BSM.CI)),(rev(1:length(realm.into.hot.rates.BSM.CI)))),y=c((unlist(lapply(realm.into.hot.rates.BSM.CI,function(x)x[1]))),rev(unlist(lapply(realm.into.hot.rates.BSM.CI,function(x)x[2])))),col=adjustcolor( "blue", alpha.f = 0.2))
    #
  }else if(plots=='cladogenesis.BAMM.rates'){
    #plot BAMM in situ rates
    plot(c(1,1),xlim=c(0,length(insitu.hot.rates.BAMM.CI)),ylim=c(0,max(c(unlist(insitu.hot.rates.BAMM.CI),unlist(insitu.realm.rates.BAMM.CI)))),type='n',xaxt='n',xlab='age',ylab='insitu.rates',main=paste(name,' cladogenesis BAMM rates',sep=''))
    lines(x=c(1:length(insitu.hot.rates.BAMM.CI)),y=unlist(lapply(insitu.hot.rates.BAMM.CI,function(x)x[3])),col=adjustcolor( "red", alpha.f = 0.9))
    polygon(x=c((1:length(insitu.hot.rates.BAMM.CI)),(rev(1:length(insitu.hot.rates.BAMM.CI)))),y=c((unlist(lapply(insitu.hot.rates.BAMM.CI,function(x)x[1]))),rev(unlist(lapply(insitu.hot.rates.BAMM.CI,function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
    lines(x=c(1:length(insitu.realm.rates.BAMM.CI)),y=unlist(lapply(insitu.realm.rates.BAMM.CI,function(x)x[3])),col=adjustcolor( "blue", alpha.f = 0.9))
    polygon(x=c((1:length(insitu.realm.rates.BAMM.CI)),(rev(1:length(insitu.realm.rates.BAMM.CI)))),y=c((unlist(lapply(insitu.realm.rates.BAMM.CI,function(x)x[1]))),rev(unlist(lapply(insitu.realm.rates.BAMM.CI,function(x)x[2])))),col=adjustcolor( "blue", alpha.f = 0.2))
    axis(1,at=c(1:length(insitu.realm.rates.BAMM.CI)),labels=c(seq(from=length(insitu.realm.rates.BAMM.CI)*2,to=2,by=-2)))
    legend('topleft',c('insitu.hot','insitu.realm'),col=c('red','blue'),lty=1,cex=.7,bty='n')
  }else if(plots=='cladogenesis.BSM.rates'){
    #plot BSM in situ rates
    plot(c(1,1),xlim=c(0,length(insitu.hot.rates.BSM.CI)),ylim=c(0,max(c(unlist(insitu.hot.rates.BSM.CI),unlist(insitu.realm.rates.BSM.CI))[is.finite(c(unlist(insitu.hot.rates.BSM.CI),unlist(insitu.realm.rates.BSM.CI)))])),type='n',xaxt='n',xlab='age',ylab='insitu.rates',main=paste(name,' cladogenesis BSM rates',sep=''))
    lines(x=c(1:length(insitu.hot.rates.BSM.CI)),y=unlist(lapply(insitu.hot.rates.BSM.CI,function(x)x[3])),col=adjustcolor( "red", alpha.f = 0.9))
    polygon(x=c((1:length(insitu.hot.rates.BSM.CI)),(rev(1:length(insitu.hot.rates.BSM.CI)))),y=c((unlist(lapply(insitu.hot.rates.BSM.CI,function(x)x[1]))),rev(unlist(lapply(insitu.hot.rates.BSM.CI,function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
    lines(x=c(1:length(insitu.realm.rates.BSM.CI)),y=unlist(lapply(insitu.realm.rates.BSM.CI,function(x)x[3])),col=adjustcolor( "blue", alpha.f = 0.9))
    polygon(x=c((1:length(insitu.realm.rates.BSM.CI)),(rev(1:length(insitu.realm.rates.BSM.CI)))),y=c((unlist(lapply(insitu.realm.rates.BSM.CI,function(x)x[1]))),rev(unlist(lapply(insitu.realm.rates.BSM.CI,function(x)x[2])))),col=adjustcolor( "blue", alpha.f = 0.2))
    axis(1,at=c(1:length(insitu.realm.rates.BSM.CI)),labels=c(seq(from=length(insitu.realm.rates.BSM.CI)*2,to=2,by=-2)))
    legend('topleft',c('insitu.hot','insitu.realm'),col=c('red','blue'),lty=1,cex=.7,bty='n')
  }else if(plots=='dispersal.events'){
    plot(c(1,1),xlim=c(0,length(hot.into.realm.time.bins.counts.rev.cumsum.CI)),ylim=c(0,max(c(unlist(hot.into.realm.time.bins.counts.rev.cumsum.CI),unlist(realm.into.hot.time.bins.counts.rev.cumsum.CI)))),type='n',xaxt='n',xlab='age',ylab='cumulative # of events',main=paste(name,' dispersal events',sep=''))
    lines(x=c(1:length(hot.into.realm.time.bins.counts.rev.cumsum.CI)),y=unlist(lapply(hot.into.realm.time.bins.counts.rev.cumsum.CI,function(x)x[3])),col=adjustcolor( "red", alpha.f = 0.9))
    polygon(x=c((1:length(hot.into.realm.time.bins.counts.rev.cumsum.CI)),(rev(1:length(hot.into.realm.time.bins.counts.rev.cumsum.CI)))),y=c((unlist(lapply(hot.into.realm.time.bins.counts.rev.cumsum.CI,function(x)x[1]))),rev(unlist(lapply(hot.into.realm.time.bins.counts.rev.cumsum.CI,function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
    lines(x=c(1:length(realm.into.hot.time.bins.counts.rev.cumsum.CI)),y=unlist(lapply(realm.into.hot.time.bins.counts.rev.cumsum.CI,function(x)x[3])),col=adjustcolor( "blue", alpha.f = 0.9))
    polygon(x=c((1:length(realm.into.hot.time.bins.counts.rev.cumsum.CI)),(rev(1:length(realm.into.hot.time.bins.counts.rev.cumsum.CI)))),y=c((unlist(lapply(realm.into.hot.time.bins.counts.rev.cumsum.CI,function(x)x[1]))),rev(unlist(lapply(realm.into.hot.time.bins.counts.rev.cumsum.CI,function(x)x[2])))),col=adjustcolor( "blue", alpha.f = 0.2))
    axis(1,at=c(1:length(realm.into.hot.time.bins.counts.rev.cumsum.CI)),labels=c(seq(from=length(realm.into.hot.time.bins.counts.rev.cumsum.CI)*2,to=2,by=-2)))
    
    legend('topleft',c('hot.into.realm','realm.into.hot'),col=c('red','blue'),lty=1,cex=.7,bty='n')
    
  }else if(plots=='cladogenesis.events'){
    plot(c(1,1),xlim=c(0,length(insitu.hot.time.bins.counts.rev.cumsum.CI)),ylim=c(0,max(c(unlist(insitu.hot.time.bins.counts.rev.cumsum.CI),unlist(insitu.realm.time.bins.counts.rev.cumsum.CI)))),type='n',xaxt='n',xlab='age',ylab='cumulative # of events',main=paste(name,' cladogenetic events',sep=''))
    lines(x=c(1:length(insitu.hot.time.bins.counts.rev.cumsum.CI)),y=unlist(lapply(insitu.hot.time.bins.counts.rev.cumsum.CI,function(x)x[3])),col=adjustcolor( "red", alpha.f = 0.9))
    polygon(x=c((1:length(insitu.hot.time.bins.counts.rev.cumsum.CI)),(rev(1:length(insitu.hot.time.bins.counts.rev.cumsum.CI)))),y=c((unlist(lapply(insitu.hot.time.bins.counts.rev.cumsum.CI,function(x)x[1]))),rev(unlist(lapply(insitu.hot.time.bins.counts.rev.cumsum.CI,function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
    lines(x=c(1:length(insitu.realm.time.bins.counts.rev.cumsum.CI)),y=unlist(lapply(insitu.realm.time.bins.counts.rev.cumsum.CI,function(x)x[3])),col=adjustcolor( "blue", alpha.f = 0.9))
    polygon(x=c((1:length(insitu.realm.time.bins.counts.rev.cumsum.CI)),(rev(1:length(insitu.realm.time.bins.counts.rev.cumsum.CI)))),y=c((unlist(lapply(insitu.realm.time.bins.counts.rev.cumsum.CI,function(x)x[1]))),rev(unlist(lapply(insitu.realm.time.bins.counts.rev.cumsum.CI,function(x)x[2])))),col=adjustcolor( "blue", alpha.f = 0.2))
    axis(1,at=c(1:length(insitu.realm.time.bins.counts.rev.cumsum.CI)),labels=c(seq(from=length(insitu.realm.time.bins.counts.rev.cumsum.CI)*2,to=2,by=-2)))
    legend('topleft',c('insitu.hot','insitu.realm'),col=c('red','blue'),lty=1,cex=.7,bty='n')
  }else if (plots=='dispersal.all.rates'){
    plot(c(1,1),xlim=c(0,length(dispersal.into.realm.rates.BSM.CI)),ylim=c(0,max(c(unlist(dispersal.into.realm.rates.BSM.CI),unlist(dispersal.into.hot.rates.BSM.CI))[is.finite(c(unlist(dispersal.into.realm.rates.BSM.CI),unlist(dispersal.into.hot.rates.BSM.CI)))])),type='n',xaxt='n',xlab='age',ylab='dispersal rates',main=paste(name,' dispersal rates',sep=''))
    lines(x=c(1:length(dispersal.into.realm.rates.BSM.CI)),y=unlist(lapply(dispersal.into.realm.rates.BSM.CI,function(x)x[3])),col=adjustcolor( "red", alpha.f = 0.9))
    polygon(x=c((1:length(dispersal.into.realm.rates.BSM.CI)),(rev(1:length(dispersal.into.realm.rates.BSM.CI)))),y=c((unlist(lapply(dispersal.into.realm.rates.BSM.CI,function(x)x[1]))),rev(unlist(lapply(dispersal.into.realm.rates.BSM.CI,function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
    lines(x=c(1:length(dispersal.into.hot.rates.BSM.CI)),y=unlist(lapply(dispersal.into.hot.rates.BSM.CI,function(x)x[3])),col=adjustcolor( "blue", alpha.f = 0.9))
    polygon(x=c((1:length(dispersal.into.hot.rates.BSM.CI)),(rev(1:length(dispersal.into.hot.rates.BSM.CI)))),y=c((unlist(lapply(dispersal.into.hot.rates.BSM.CI,function(x)x[1]))),rev(unlist(lapply(dispersal.into.hot.rates.BSM.CI,function(x)x[2])))),col=adjustcolor( "blue", alpha.f = 0.2))
    axis(1,at=c(1:length(dispersal.into.hot.rates.BSM.CI)),labels=c(seq(from=length(dispersal.into.hot.rates.BSM.CI)*2,to=2,by=-2)))
    
    legend('topleft',c('dispersal.into.hot','dispersal.into.realm'),col=c('blue','red'),lty=1,cex=.7,bty='n')
  }else if(plots=='dispersal.all.events'){
    plot(c(1,1),xlim=c(0,length(dispersal.into.realm.time.bins.counts.rev.cumsum.CI)),ylim=c(0,max(c(unlist(dispersal.into.realm.time.bins.counts.rev.cumsum.CI),unlist(dispersal.into.hot.time.bins.counts.rev.cumsum.CI))[is.finite(c(unlist(dispersal.into.realm.time.bins.counts.rev.cumsum.CI),unlist(dispersal.into.hot.time.bins.counts.rev.cumsum.CI)))])),type='n',xaxt='n',xlab='age',ylab='cumulative # of events',main=paste(name,' dispersal events',sep=''))
    lines(x=c(1:length(dispersal.into.realm.time.bins.counts.rev.cumsum.CI)),y=unlist(lapply(dispersal.into.realm.time.bins.counts.rev.cumsum.CI,function(x)x[3])),col=adjustcolor( "red", alpha.f = 0.9))
    polygon(x=c((1:length(dispersal.into.realm.time.bins.counts.rev.cumsum.CI)),(rev(1:length(dispersal.into.realm.time.bins.counts.rev.cumsum.CI)))),y=c((unlist(lapply(dispersal.into.realm.time.bins.counts.rev.cumsum.CI,function(x)x[1]))),rev(unlist(lapply(dispersal.into.realm.time.bins.counts.rev.cumsum.CI,function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
    lines(x=c(1:length(dispersal.into.hot.time.bins.counts.rev.cumsum.CI)),y=unlist(lapply(dispersal.into.hot.time.bins.counts.rev.cumsum.CI,function(x)x[3])),col=adjustcolor( "blue", alpha.f = 0.9))
    polygon(x=c((1:length(dispersal.into.hot.time.bins.counts.rev.cumsum.CI)),(rev(1:length(dispersal.into.hot.time.bins.counts.rev.cumsum.CI)))),y=c((unlist(lapply(dispersal.into.hot.time.bins.counts.rev.cumsum.CI,function(x)x[1]))),rev(unlist(lapply(dispersal.into.hot.time.bins.counts.rev.cumsum.CI,function(x)x[2])))),col=adjustcolor( "blue", alpha.f = 0.2))
    axis(1,at=c(1:length(dispersal.into.hot.time.bins.counts.rev.cumsum.CI)),labels=c(seq(from=length(dispersal.into.hot.time.bins.counts.rev.cumsum.CI)*2,to=2,by=-2)))
    legend('topleft',c('dispersal.into.realm','dispersal.into.hot'),col=c('red','blue'),lty=1,cex=.7,bty='n')
  }else if(plots=='dispersal.into.from.events'){
    plot(c(1,1),xlim=c(0,length(dispersal.into.realm.time.bins.counts.rev.cumsum.CI)),ylim=c(0,max(c(unlist(dispersal.into.realm.time.bins.counts.rev.cumsum.CI),unlist(dispersal.into.hot.time.bins.counts.rev.cumsum.CI),unlist(dispersal.from.hot.time.bins.counts.rev.cumsum.CI),unlist(dispersal.from.realm.time.bins.counts.rev.cumsum.CI)))),type='n',xaxt='n',xlab='age',ylab='cumulative # of events',main=paste(name,' dispersal events',sep=''))
    lines(x=c(1:length(dispersal.into.realm.time.bins.counts.rev.cumsum.CI)),y=unlist(lapply(dispersal.into.realm.time.bins.counts.rev.cumsum.CI,function(x)x[3])),col=adjustcolor( "green", alpha.f = 0.9))
    polygon(x=c((1:length(dispersal.into.realm.time.bins.counts.rev.cumsum.CI)),(rev(1:length(dispersal.into.realm.time.bins.counts.rev.cumsum.CI)))),y=c((unlist(lapply(dispersal.into.realm.time.bins.counts.rev.cumsum.CI,function(x)x[1]))),rev(unlist(lapply(dispersal.into.realm.time.bins.counts.rev.cumsum.CI,function(x)x[2])))),col=adjustcolor( "green", alpha.f = 0.2))
    lines(x=c(1:length(dispersal.into.hot.time.bins.counts.rev.cumsum.CI)),y=unlist(lapply(dispersal.into.hot.time.bins.counts.rev.cumsum.CI,function(x)x[3])),col=adjustcolor( "orange", alpha.f = 0.9))
    polygon(x=c((1:length(dispersal.into.hot.time.bins.counts.rev.cumsum.CI)),(rev(1:length(dispersal.into.hot.time.bins.counts.rev.cumsum.CI)))),y=c((unlist(lapply(dispersal.into.hot.time.bins.counts.rev.cumsum.CI,function(x)x[1]))),rev(unlist(lapply(dispersal.into.hot.time.bins.counts.rev.cumsum.CI,function(x)x[2])))),col=adjustcolor( "orange", alpha.f = 0.2))
    lines(x=c(1:length(dispersal.from.realm.time.bins.counts.rev.cumsum.CI)),y=unlist(lapply(dispersal.from.realm.time.bins.counts.rev.cumsum.CI,function(x)x[3])),col=adjustcolor( "blue", alpha.f = 0.9))
    polygon(x=c((1:length(dispersal.from.realm.time.bins.counts.rev.cumsum.CI)),(rev(1:length(dispersal.from.realm.time.bins.counts.rev.cumsum.CI)))),y=c((unlist(lapply(dispersal.from.realm.time.bins.counts.rev.cumsum.CI,function(x)x[1]))),rev(unlist(lapply(dispersal.from.realm.time.bins.counts.rev.cumsum.CI,function(x)x[2])))),col=adjustcolor( "blue", alpha.f = 0.2))
    lines(x=c(1:length(dispersal.from.hot.time.bins.counts.rev.cumsum.CI)),y=unlist(lapply(dispersal.from.hot.time.bins.counts.rev.cumsum.CI,function(x)x[3])),col=adjustcolor( "red", alpha.f = 0.9))
    polygon(x=c((1:length(dispersal.from.hot.time.bins.counts.rev.cumsum.CI)),(rev(1:length(dispersal.from.hot.time.bins.counts.rev.cumsum.CI)))),y=c((unlist(lapply(dispersal.from.hot.time.bins.counts.rev.cumsum.CI,function(x)x[1]))),rev(unlist(lapply(dispersal.from.hot.time.bins.counts.rev.cumsum.CI,function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
    axis(1,at=c(1:length(dispersal.from.hot.time.bins.counts.rev.cumsum.CI)),labels=c(seq(from=length(dispersal.from.hot.time.bins.counts.rev.cumsum.CI)*2,to=2,by=-2)))
    legend('topleft',c('dispersal.from.realm','dispersal.from.hot','dispersal.into.realm','dispersal.into.hot'),col=c('blue','red','green','orange'),lty=1,cex=.7,bty='n')
  }
  
  #plot BSM dispersal rates hot into realm
  
  ##plot BAMM in situ rates
  #par(new=T)
  #plot(c(1,1),xlim=c(0,length(insitu.hot.rates.BAMM.CI)),ylim=c(0,max(c(unlist(insitu.hot.rates.BAMM.CI),unlist(insitu.realm.rates.BAMM.CI)))),type='n',xaxt='n',xlab='age',yaxt='n')
  ##polygon(x=c(c(age:length(insitu.hot.rates.BAMM.CI)),c(length(insitu.hot.rates.BAMM.CI):age)),y=c(unlist(lapply(insitu.hot.rates.BAMM.CI[c(age:length(insitu.hot.rates.BAMM.CI))],function(x)x[1])),rev(unlist(lapply(insitu.hot.rates.BAMM.CI[c(age:length(insitu.hot.rates.BAMM.CI))],function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
  #lines(x=c(1:length(insitu.hot.rates.BAMM.CI)),y=unlist(lapply(insitu.hot.rates.BAMM.CI,function(x)x[3])),col=adjustcolor( "red", alpha.f = 0.9))
  #polygon(x=c((1:length(insitu.hot.rates.BAMM.CI)),(rev(1:length(insitu.hot.rates.BAMM.CI)))),y=c((unlist(lapply(insitu.hot.rates.BAMM.CI,function(x)x[1]))),rev(unlist(lapply(insitu.hot.rates.BAMM.CI,function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
  #lines(x=c(1:length(insitu.realm.rates.BAMM.CI)),y=unlist(lapply(insitu.realm.rates.BAMM.CI,function(x)x[3])),col=adjustcolor( "blue", alpha.f = 0.9))
  #polygon(x=c((1:length(insitu.realm.rates.BAMM.CI)),(rev(1:length(insitu.realm.rates.BAMM.CI)))),y=c((unlist(lapply(insitu.realm.rates.BAMM.CI,function(x)x[1]))),rev(unlist(lapply(insitu.realm.rates.BAMM.CI,function(x)x[2])))),col=adjustcolor( "blue", alpha.f = 0.2))
  ##axis(1,at=c(0:length(BSMoutput.real.rates[['insitu.realm.rates.BSM.CI']])),labels=c((length(BSMoutput.real.rates[['insitu.realm.rates.BSM.CI']])+2):2))
  #axis(4)
  #legend('topright',c('insitu.hot','insitu.realm'),col=c('red','blue'),lty=1,cex=.7,bty='n')
}

plot_rates_events_real_2mybin_difference_both_scaled_fixedyaxis_average<-function(BSMoutput.real.rates.mammals,BSMoutput.real.events.mammals,BSMoutput.real.rates.birds,BSMoutput.real.events.birds,age,name,plots,ylim){
  #correct age
  age.mammals<-(length(BSMoutput.real.rates.mammals[['insitu.realm.rates.BSM.CI']])*2-age)/2
  #trim lists to age period
  BSMoutput.real.rates.mammals<-lapply(BSMoutput.real.rates.mammals,function(x)x[age.mammals:length(x)])
  BSMoutput.real.events.mammals<-lapply(BSMoutput.real.events.mammals,function(x)x[age.mammals:length(x)])
  #delete Inf in rates and events
  BSMoutput.real.rates.mammals<-lapply(BSMoutput.real.rates.mammals,function(x) lapply(x,function(x) {x[!is.finite(x)]<-NA;return(x)}))
  BSMoutput.real.events.mammals<-lapply(BSMoutput.real.events.mammals,function(x) lapply(x,function(x) {x[!is.finite(x)]<-NA;return(x)}))
  #remove the last bin - two million year time bin
  BSMoutput.real.rates.mammals<-lapply(BSMoutput.real.rates.mammals,function(x) x<-x[-c(length(x))])
  BSMoutput.real.events.mammals<-lapply(BSMoutput.real.events.mammals,function(x) x<-x[-c((length(x)-1),length(x))])
  #assign variables for easier handling
  hot.into.realm.rates.BSM.mammals<-BSMoutput.real.rates.mammals[['hot.into.realm.rates.BSM.CI']]
  realm.into.hot.rates.BSM.mammals<-BSMoutput.real.rates.mammals[['realm.into.hot.rates.BSM.CI']]
  insitu.realm.rates.BAMM.mammals<-BSMoutput.real.rates.mammals[['insitu.realm.rates.BAMM.CI']]
  insitu.hot.rates.BAMM.mammals<-BSMoutput.real.rates.mammals[['insitu.hot.rates.BAMM.CI']]
  insitu.realm.rates.BSM.mammals<-BSMoutput.real.rates.mammals[['insitu.realm.rates.BSM.CI']]
  insitu.hot.rates.BSM.mammals<-BSMoutput.real.rates.mammals[['insitu.hot.rates.BSM.CI']]
  dispersal.into.realm.rates.BSM.mammals<-BSMoutput.real.rates.mammals[['dispersal.into.realm.rates.BSM.CI']]
  dispersal.into.hot.rates.BSM.mammals<-BSMoutput.real.rates.mammals[['dispersal.into.hot.rates.BSM.CI']]
  dispersal.from.realm.rates.BSM.mammals<-BSMoutput.real.rates.mammals[['dispersal.from.realm.rates.BSM.CI']]
  dispersal.from.hot.rates.BSM.mammals<-BSMoutput.real.rates.mammals[['dispersal.from.hot.rates.BSM.CI']]
  
  hot.into.realm.time.bins.counts.rev.cumsum.mammals<-BSMoutput.real.events.mammals[['hot.into.realm.time.bins.counts.rev.cumsum.CI']]
  realm.into.hot.time.bins.counts.rev.cumsum.mammals<-BSMoutput.real.events.mammals[['realm.into.hot.time.bins.counts.rev.cumsum.CI']]
  insitu.realm.time.bins.counts.rev.cumsum.mammals<-BSMoutput.real.events.mammals[['insitu.realm.time.bins.counts.rev.cumsum.CI']]
  insitu.hot.time.bins.counts.rev.cumsum.mammals<-BSMoutput.real.events.mammals[['insitu.hot.time.bins.counts.rev.cumsum.CI']]
  dispersal.into.realm.time.bins.counts.rev.cumsum.mammals<-BSMoutput.real.events.mammals[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']]
  dispersal.into.hot.time.bins.counts.rev.cumsum.mammals<-BSMoutput.real.events.mammals[['dispersal.into.hot.time.bins.counts.rev.cumsum.CI']]
  dispersal.from.realm.time.bins.counts.rev.cumsum.mammals<-BSMoutput.real.events.mammals[['dispersal.from.realm.time.bins.counts.rev.cumsum.CI']]
  dispersal.from.hot.time.bins.counts.rev.cumsum.mammals<-BSMoutput.real.events.mammals[['dispersal.from.hot.time.bins.counts.rev.cumsum.CI']]
  
  age.birds<-(length(BSMoutput.real.rates.birds[['insitu.realm.rates.BSM.CI']])*2-age)/2
  #trim lists to age period
  BSMoutput.real.rates.birds<-lapply(BSMoutput.real.rates.birds,function(x)x[age.birds:length(x)])
  BSMoutput.real.events.birds<-lapply(BSMoutput.real.events.birds,function(x)x[age.birds:length(x)])
  #delete Inf in rates and events
  BSMoutput.real.rates.birds<-lapply(BSMoutput.real.rates.birds,function(x) lapply(x,function(x) {x[!is.finite(x)]<-NA;return(x)}))
  BSMoutput.real.events.birds<-lapply(BSMoutput.real.events.birds,function(x) lapply(x,function(x) {x[!is.finite(x)]<-NA;return(x)}))
  #remove the last bin - two million year time bin
  BSMoutput.real.rates.birds<-lapply(BSMoutput.real.rates.birds,function(x) x<-x[-c(length(x))])
  BSMoutput.real.events.birds<-lapply(BSMoutput.real.events.birds,function(x) x<-x[-c((length(x)-1),length(x))])
  #assign variables for easier handling
  hot.into.realm.rates.BSM.birds<-BSMoutput.real.rates.birds[['hot.into.realm.rates.BSM.CI']]
  realm.into.hot.rates.BSM.birds<-BSMoutput.real.rates.birds[['realm.into.hot.rates.BSM.CI']]
  insitu.realm.rates.BAMM.birds<-BSMoutput.real.rates.birds[['insitu.realm.rates.BAMM.CI']]
  insitu.hot.rates.BAMM.birds<-BSMoutput.real.rates.birds[['insitu.hot.rates.BAMM.CI']]
  insitu.realm.rates.BSM.birds<-BSMoutput.real.rates.birds[['insitu.realm.rates.BSM.CI']]
  insitu.hot.rates.BSM.birds<-BSMoutput.real.rates.birds[['insitu.hot.rates.BSM.CI']]
  dispersal.into.realm.rates.BSM.birds<-BSMoutput.real.rates.birds[['dispersal.into.realm.rates.BSM.CI']]
  dispersal.into.hot.rates.BSM.birds<-BSMoutput.real.rates.birds[['dispersal.into.hot.rates.BSM.CI']]
  dispersal.from.realm.rates.BSM.birds<-BSMoutput.real.rates.birds[['dispersal.from.realm.rates.BSM.CI']]
  dispersal.from.hot.rates.BSM.birds<-BSMoutput.real.rates.birds[['dispersal.from.hot.rates.BSM.CI']]
  
  hot.into.realm.time.bins.counts.rev.cumsum.birds<-BSMoutput.real.events.birds[['hot.into.realm.time.bins.counts.rev.cumsum.CI']]
  realm.into.hot.time.bins.counts.rev.cumsum.birds<-BSMoutput.real.events.birds[['realm.into.hot.time.bins.counts.rev.cumsum.CI']]
  insitu.realm.time.bins.counts.rev.cumsum.birds<-BSMoutput.real.events.birds[['insitu.realm.time.bins.counts.rev.cumsum.CI']]
  insitu.hot.time.bins.counts.rev.cumsum.birds<-BSMoutput.real.events.birds[['insitu.hot.time.bins.counts.rev.cumsum.CI']]
  dispersal.into.realm.time.bins.counts.rev.cumsum.birds<-BSMoutput.real.events.birds[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']]
  dispersal.into.hot.time.bins.counts.rev.cumsum.birds<-BSMoutput.real.events.birds[['dispersal.into.hot.time.bins.counts.rev.cumsum.CI']]
  dispersal.from.realm.time.bins.counts.rev.cumsum.birds<-BSMoutput.real.events.birds[['dispersal.from.realm.time.bins.counts.rev.cumsum.CI']]
  dispersal.from.hot.time.bins.counts.rev.cumsum.birds<-BSMoutput.real.events.birds[['dispersal.from.hot.time.bins.counts.rev.cumsum.CI']]
  
  if(plots=='within.realm.dispersal.rates'){
    #difference.into.hotminusrealm<-lapply(c(1:length(dispersal.into.realm.rates.BSM)),function(x)sample(dispersal.into.hot.rates.BSM[[x]],size = 100,replace=TRUE)-sample(dispersal.into.realm.rates.BSM[[x]],size = 100,replace=TRUE))
    difference.hotintorealm.realmintohot.mammals<-lapply(c(1:length(hot.into.realm.rates.BSM.mammals)),function(x)sample(hot.into.realm.rates.BSM.mammals[[x]],size = 50,replace=FALSE)-sample(realm.into.hot.rates.BSM.mammals[[x]],size = 50,replace=FALSE))
    difference.hotintorealm.realmintohot.birds<-lapply(c(1:length(hot.into.realm.rates.BSM.birds)),function(x)sample(hot.into.realm.rates.BSM.birds[[x]],size = 50,replace=FALSE)-sample(realm.into.hot.rates.BSM.birds[[x]],size = 50,replace=FALSE))
    #global.mean<-mean(c(c(unlist(difference.hotintorealm.realmintohot.mammals),unlist(difference.hotintorealm.realmintohot.birds))),na.rm=TRUE)
    global.sd<-sd(c(c(unlist(difference.hotintorealm.realmintohot.mammals),unlist(difference.hotintorealm.realmintohot.birds))),na.rm=TRUE)
    #difference.hotintorealm.realmintohot.mammals<-lapply(difference.hotintorealm.realmintohot.mammals,function(x)(x-global.mean)/global.sd)
    #difference.hotintorealm.realmintohot.birds<-lapply(difference.hotintorealm.realmintohot.birds,function(x)(x-global.mean)/global.sd)
    difference.hotintorealm.realmintohot.mammals<-lapply(difference.hotintorealm.realmintohot.mammals,function(x)x/global.sd)
    difference.hotintorealm.realmintohot.birds<-lapply(difference.hotintorealm.realmintohot.birds,function(x)x/global.sd)
    
    median.mammals<-unlist(lapply(difference.hotintorealm.realmintohot.mammals,function(x)median(x,na.rm=TRUE)))
    lower.mammals<-unlist(lapply(difference.hotintorealm.realmintohot.mammals,function(x)quantile(x,c(0.1),na.rm=TRUE)))
    upper.mammals<-unlist(lapply(difference.hotintorealm.realmintohot.mammals,function(x)quantile(x,c(0.9),na.rm=TRUE)))
    median.birds<-unlist(lapply(difference.hotintorealm.realmintohot.birds,function(x)median(x,na.rm=TRUE)))
    lower.birds<-unlist(lapply(difference.hotintorealm.realmintohot.birds,function(x)quantile(x,c(0.1),na.rm=TRUE)))
    upper.birds<-unlist(lapply(difference.hotintorealm.realmintohot.birds,function(x)quantile(x,c(0.9),na.rm=TRUE)))
    list.means<-list(c(mean(lower.mammals,na.rm=T),mean(median.mammals,na.rm=T),mean(upper.mammals,na.rm=T)),c(mean(lower.birds,na.rm=T),mean(median.birds,na.rm=T),mean(upper.birds,na.rm=T)))
    names(list.means)<-c('mammals','birds')
    names(list.means$mammals)<-c('q10','median','q90')
    names(list.means$birds)<-c('q10','median','q90')
    return(list.means)
  }else if(plots=='within.realm.cladogenesis.BAMM.rates'){
    difference.insituhot.insiturealm.mammals<-lapply(c(1:length(insitu.hot.rates.BAMM.mammals)),function(x)sample(insitu.hot.rates.BAMM.mammals[[x]],size = 100,replace=TRUE)-sample(insitu.realm.rates.BAMM.mammals[[x]],size = 100,replace=TRUE))
    difference.insituhot.insiturealm.birds<-lapply(c(1:length(insitu.hot.rates.BAMM.birds)),function(x)sample(insitu.hot.rates.BAMM.birds[[x]],size = 100,replace=TRUE)-sample(insitu.realm.rates.BAMM.birds[[x]],size = 100,replace=TRUE))
    global.sd<-sd(c(c(unlist(difference.insituhot.insiturealm.mammals),unlist(difference.insituhot.insiturealm.birds))),na.rm=TRUE)
    difference.insituhot.insiturealm.mammals<-lapply(difference.insituhot.insiturealm.mammals,function(x)x/global.sd)
    difference.insituhot.insiturealm.birds<-lapply(difference.insituhot.insiturealm.birds,function(x)x/global.sd)
    median.mammals<-unlist(lapply(difference.insituhot.insiturealm.mammals,function(x)median(x,na.rm=TRUE)))
    lower.mammals<-unlist(lapply(difference.insituhot.insiturealm.mammals,function(x)quantile(x,c(0.1),na.rm=TRUE)))
    upper.mammals<-unlist(lapply(difference.insituhot.insiturealm.mammals,function(x)quantile(x,c(0.9),na.rm=TRUE)))
    median.birds<-unlist(lapply(difference.insituhot.insiturealm.birds,function(x)median(x,na.rm=TRUE)))
    lower.birds<-unlist(lapply(difference.insituhot.insiturealm.birds,function(x)quantile(x,c(0.1),na.rm=TRUE)))
    upper.birds<-unlist(lapply(difference.insituhot.insiturealm.birds,function(x)quantile(x,c(0.9),na.rm=TRUE)))
    list.means<-list(c(mean(lower.mammals,na.rm=T),mean(median.mammals,na.rm=T),mean(upper.mammals,na.rm=T)),c(mean(lower.birds,na.rm=T),mean(median.birds,na.rm=T),mean(upper.birds,na.rm=T)))
    names(list.means)<-c('mammals','birds')
    names(list.means$mammals)<-c('q10','median','q90')
    names(list.means$birds)<-c('q10','median','q90')
    return(list.means)
  }else if(plots=='within.realm.cladogenesis.BSM.rates'){
    difference.insituhot.insiturealm.mammals<-lapply(c(1:length(insitu.hot.rates.BSM.mammals)),function(x)sample(insitu.hot.rates.BSM.mammals[[x]],size = 100,replace=TRUE)-sample(insitu.realm.rates.BSM.mammals[[x]],size = 100,replace=TRUE))
    difference.insituhot.insiturealm.birds<-lapply(c(1:length(insitu.hot.rates.BSM.birds)),function(x)sample(insitu.hot.rates.BSM.birds[[x]],size = 100,replace=TRUE)-sample(insitu.realm.rates.BSM.birds[[x]],size = 100,replace=TRUE))
    median.mammals<-unlist(lapply(difference.insituhot.insiturealm.mammals,function(x)median(x,na.rm=TRUE)))
    lower.mammals<-unlist(lapply(difference.insituhot.insiturealm.mammals,function(x)quantile(x,c(0.1),na.rm=TRUE)))
    upper.mammals<-unlist(lapply(difference.insituhot.insiturealm.mammals,function(x)quantile(x,c(0.9),na.rm=TRUE)))
    median.birds<-unlist(lapply(difference.insituhot.insiturealm.birds,function(x)median(x,na.rm=TRUE)))
    lower.birds<-unlist(lapply(difference.insituhot.insiturealm.birds,function(x)quantile(x,c(0.1),na.rm=TRUE)))
    upper.birds<-unlist(lapply(difference.insituhot.insiturealm.birds,function(x)quantile(x,c(0.9),na.rm=TRUE)))
    list.means<-list(c(mean(lower.mammals,na.rm=T),mean(median.mammals,na.rm=T),mean(upper.mammals,na.rm=T)),c(mean(lower.birds,na.rm=T),mean(median.birds,na.rm=T),mean(upper.birds,na.rm=T)))
    names(list.means)<-c('mammals','birds')
    names(list.means$mammals)<-c('q10','median','q90')
    names(list.means$birds)<-c('q10','median','q90')
    return(list.means)
    
  }else if(plots=='global.dispersal.rates'){
    #difference.into.hotminusrealm<-lapply(c(1:length(dispersal.into.realm.rates.BSM)),function(x)sample(dispersal.into.hot.rates.BSM[[x]],size = 100,replace=TRUE)-sample(dispersal.into.realm.rates.BSM[[x]],size = 100,replace=TRUE))
    difference.hotrealm.from<-lapply(c(1:length(dispersal.from.realm.rates.BSM)),function(x)sample(dispersal.from.hot.rates.BSM[[x]],size = 100,replace=TRUE)-sample(dispersal.from.realm.rates.BSM[[x]],size = 100,replace=TRUE))
    difference.hotrealm.into<-lapply(c(1:length(dispersal.into.realm.rates.BSM)),function(x)sample(dispersal.into.hot.rates.BSM[[x]],size = 100,replace=TRUE)-sample(dispersal.into.realm.rates.BSM[[x]],size = 100,replace=TRUE))
    median.from<-unlist(lapply(difference.hotrealm.from,function(x)median(x,na.rm=TRUE)))
    lower.from<-unlist(lapply(difference.hotrealm.from,function(x)quantile(x,c(0.1),na.rm=TRUE)))
    upper.from<-unlist(lapply(difference.hotrealm.from,function(x)quantile(x,c(0.9),na.rm=TRUE)))
    median.into<-unlist(lapply(difference.hotrealm.into,function(x)median(x,na.rm=TRUE)))
    lower.into<-unlist(lapply(difference.hotrealm.into,function(x)quantile(x,c(0.1),na.rm=TRUE)))
    upper.into<-unlist(lapply(difference.hotrealm.into,function(x)quantile(x,c(0.9),na.rm=TRUE)))
    
    plot(x=c(1:length(median.from)),y=median.from,ylim=c(-max(abs(c(median.from,lower.from,upper.from,median.into,lower.into,upper.into))),max(abs(c(median.from,lower.from,upper.from,median.into,lower.into,upper.into)))),type='n',xaxt='n',xlab='age',ylab='diff. dispersal hot - realm',main=paste(name,' diff in global dispersal rates',sep=''))
    polygon(x=c(c(1:length(median.from),length(median.from):1)),y=c(lower.from,rev(upper.from)),col=adjustcolor('red',alpha.f = 0.25))
    lines(x=c(1:length(median.from)),y=median.from,col='red') 
    abline(h=0,lty=2)
    polygon(x=c(c(1:length(median.into),length(median.into):1)),y=c(lower.into,rev(upper.into)),col=adjustcolor('blue',alpha.f = 0.25))
    lines(x=c(1:length(median.into)),y=median.into,col='blue') 
    axis(1,at=c(1:length(difference.hotintorealm.realmintohot)),labels=c(seq(from=length(difference.hotintorealm.realmintohot)*2,to=2,by=-2)))
    legend('topleft',c('dispersal.from.diff','dispersal.into.diff'),col=c('red','blue'),lty=1,cex=.7,bty='n')
    #plot(c(1,1),xlim=c(0,length(hot.into.realm.rates.BSM.CI)),ylim=c(0,max(c(unlist(hot.into.realm.rates.BSM.CI),unlist(realm.into.hot.rates.BSM.CI))[is.finite(c(unlist(hot.into.realm.rates.BSM.CI),unlist(realm.into.hot.rates.BSM.CI)))])),type='n',xaxt='n',xlab='age',ylab='dispersal rates',main=paste(name,' dispersal rates',sep=''))
    #lines(x=c(1:length(hot.into.realm.rates.BSM.CI)),y=unlist(lapply(hot.into.realm.rates.BSM.CI,function(x)x[3])),col=adjustcolor( "red", alpha.f = 0.9))
    #polygon(x=c((1:length(hot.into.realm.rates.BSM.CI)),(rev(1:length(hot.into.realm.rates.BSM.CI)))),y=c((unlist(lapply(hot.into.realm.rates.BSM.CI,function(x)x[1]))),rev(unlist(lapply(hot.into.realm.rates.BSM.CI,function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
    #lines(x=c(1:length(realm.into.hot.rates.BSM.CI)),y=unlist(lapply(realm.into.hot.rates.BSM.CI,function(x)x[3])),col=adjustcolor( "blue", alpha.f = 0.9))
    #polygon(x=c((1:length(realm.into.hot.rates.BSM.CI)),(rev(1:length(realm.into.hot.rates.BSM.CI)))),y=c((unlist(lapply(realm.into.hot.rates.BSM.CI,function(x)x[1]))),rev(unlist(lapply(realm.into.hot.rates.BSM.CI,function(x)x[2])))),col=adjustcolor( "blue", alpha.f = 0.2))
    #
  }else if(plots=='cladogenesis.BAMM.rates'){
    #plot BAMM in situ rates
    plot(c(1,1),xlim=c(0,length(insitu.hot.rates.BAMM.CI)),ylim=c(0,max(c(unlist(insitu.hot.rates.BAMM.CI),unlist(insitu.realm.rates.BAMM.CI)))),type='n',xaxt='n',xlab='age',ylab='insitu.rates',main=paste(name,' cladogenesis BAMM rates',sep=''))
    lines(x=c(1:length(insitu.hot.rates.BAMM.CI)),y=unlist(lapply(insitu.hot.rates.BAMM.CI,function(x)x[3])),col=adjustcolor( "red", alpha.f = 0.9))
    polygon(x=c((1:length(insitu.hot.rates.BAMM.CI)),(rev(1:length(insitu.hot.rates.BAMM.CI)))),y=c((unlist(lapply(insitu.hot.rates.BAMM.CI,function(x)x[1]))),rev(unlist(lapply(insitu.hot.rates.BAMM.CI,function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
    lines(x=c(1:length(insitu.realm.rates.BAMM.CI)),y=unlist(lapply(insitu.realm.rates.BAMM.CI,function(x)x[3])),col=adjustcolor( "blue", alpha.f = 0.9))
    polygon(x=c((1:length(insitu.realm.rates.BAMM.CI)),(rev(1:length(insitu.realm.rates.BAMM.CI)))),y=c((unlist(lapply(insitu.realm.rates.BAMM.CI,function(x)x[1]))),rev(unlist(lapply(insitu.realm.rates.BAMM.CI,function(x)x[2])))),col=adjustcolor( "blue", alpha.f = 0.2))
    axis(1,at=c(1:length(insitu.realm.rates.BAMM.CI)),labels=c(seq(from=length(insitu.realm.rates.BAMM.CI)*2,to=2,by=-2)))
    legend('topleft',c('insitu.hot','insitu.realm'),col=c('red','blue'),lty=1,cex=.7,bty='n')
  }else if(plots=='cladogenesis.BSM.rates'){
    #plot BSM in situ rates
    plot(c(1,1),xlim=c(0,length(insitu.hot.rates.BSM.CI)),ylim=c(0,max(c(unlist(insitu.hot.rates.BSM.CI),unlist(insitu.realm.rates.BSM.CI))[is.finite(c(unlist(insitu.hot.rates.BSM.CI),unlist(insitu.realm.rates.BSM.CI)))])),type='n',xaxt='n',xlab='age',ylab='insitu.rates',main=paste(name,' cladogenesis BSM rates',sep=''))
    lines(x=c(1:length(insitu.hot.rates.BSM.CI)),y=unlist(lapply(insitu.hot.rates.BSM.CI,function(x)x[3])),col=adjustcolor( "red", alpha.f = 0.9))
    polygon(x=c((1:length(insitu.hot.rates.BSM.CI)),(rev(1:length(insitu.hot.rates.BSM.CI)))),y=c((unlist(lapply(insitu.hot.rates.BSM.CI,function(x)x[1]))),rev(unlist(lapply(insitu.hot.rates.BSM.CI,function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
    lines(x=c(1:length(insitu.realm.rates.BSM.CI)),y=unlist(lapply(insitu.realm.rates.BSM.CI,function(x)x[3])),col=adjustcolor( "blue", alpha.f = 0.9))
    polygon(x=c((1:length(insitu.realm.rates.BSM.CI)),(rev(1:length(insitu.realm.rates.BSM.CI)))),y=c((unlist(lapply(insitu.realm.rates.BSM.CI,function(x)x[1]))),rev(unlist(lapply(insitu.realm.rates.BSM.CI,function(x)x[2])))),col=adjustcolor( "blue", alpha.f = 0.2))
    axis(1,at=c(1:length(insitu.realm.rates.BSM.CI)),labels=c(seq(from=length(insitu.realm.rates.BSM.CI)*2,to=2,by=-2)))
    legend('topleft',c('insitu.hot','insitu.realm'),col=c('red','blue'),lty=1,cex=.7,bty='n')
  }else if(plots=='dispersal.events'){
    plot(c(1,1),xlim=c(0,length(hot.into.realm.time.bins.counts.rev.cumsum.CI)),ylim=c(0,max(c(unlist(hot.into.realm.time.bins.counts.rev.cumsum.CI),unlist(realm.into.hot.time.bins.counts.rev.cumsum.CI)))),type='n',xaxt='n',xlab='age',ylab='cumulative # of events',main=paste(name,' dispersal events',sep=''))
    lines(x=c(1:length(hot.into.realm.time.bins.counts.rev.cumsum.CI)),y=unlist(lapply(hot.into.realm.time.bins.counts.rev.cumsum.CI,function(x)x[3])),col=adjustcolor( "red", alpha.f = 0.9))
    polygon(x=c((1:length(hot.into.realm.time.bins.counts.rev.cumsum.CI)),(rev(1:length(hot.into.realm.time.bins.counts.rev.cumsum.CI)))),y=c((unlist(lapply(hot.into.realm.time.bins.counts.rev.cumsum.CI,function(x)x[1]))),rev(unlist(lapply(hot.into.realm.time.bins.counts.rev.cumsum.CI,function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
    lines(x=c(1:length(realm.into.hot.time.bins.counts.rev.cumsum.CI)),y=unlist(lapply(realm.into.hot.time.bins.counts.rev.cumsum.CI,function(x)x[3])),col=adjustcolor( "blue", alpha.f = 0.9))
    polygon(x=c((1:length(realm.into.hot.time.bins.counts.rev.cumsum.CI)),(rev(1:length(realm.into.hot.time.bins.counts.rev.cumsum.CI)))),y=c((unlist(lapply(realm.into.hot.time.bins.counts.rev.cumsum.CI,function(x)x[1]))),rev(unlist(lapply(realm.into.hot.time.bins.counts.rev.cumsum.CI,function(x)x[2])))),col=adjustcolor( "blue", alpha.f = 0.2))
    axis(1,at=c(1:length(realm.into.hot.time.bins.counts.rev.cumsum.CI)),labels=c(seq(from=length(realm.into.hot.time.bins.counts.rev.cumsum.CI)*2,to=2,by=-2)))
    
    legend('topleft',c('hot.into.realm','realm.into.hot'),col=c('red','blue'),lty=1,cex=.7,bty='n')
    
  }else if(plots=='cladogenesis.events'){
    plot(c(1,1),xlim=c(0,length(insitu.hot.time.bins.counts.rev.cumsum.CI)),ylim=c(0,max(c(unlist(insitu.hot.time.bins.counts.rev.cumsum.CI),unlist(insitu.realm.time.bins.counts.rev.cumsum.CI)))),type='n',xaxt='n',xlab='age',ylab='cumulative # of events',main=paste(name,' cladogenetic events',sep=''))
    lines(x=c(1:length(insitu.hot.time.bins.counts.rev.cumsum.CI)),y=unlist(lapply(insitu.hot.time.bins.counts.rev.cumsum.CI,function(x)x[3])),col=adjustcolor( "red", alpha.f = 0.9))
    polygon(x=c((1:length(insitu.hot.time.bins.counts.rev.cumsum.CI)),(rev(1:length(insitu.hot.time.bins.counts.rev.cumsum.CI)))),y=c((unlist(lapply(insitu.hot.time.bins.counts.rev.cumsum.CI,function(x)x[1]))),rev(unlist(lapply(insitu.hot.time.bins.counts.rev.cumsum.CI,function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
    lines(x=c(1:length(insitu.realm.time.bins.counts.rev.cumsum.CI)),y=unlist(lapply(insitu.realm.time.bins.counts.rev.cumsum.CI,function(x)x[3])),col=adjustcolor( "blue", alpha.f = 0.9))
    polygon(x=c((1:length(insitu.realm.time.bins.counts.rev.cumsum.CI)),(rev(1:length(insitu.realm.time.bins.counts.rev.cumsum.CI)))),y=c((unlist(lapply(insitu.realm.time.bins.counts.rev.cumsum.CI,function(x)x[1]))),rev(unlist(lapply(insitu.realm.time.bins.counts.rev.cumsum.CI,function(x)x[2])))),col=adjustcolor( "blue", alpha.f = 0.2))
    axis(1,at=c(1:length(insitu.realm.time.bins.counts.rev.cumsum.CI)),labels=c(seq(from=length(insitu.realm.time.bins.counts.rev.cumsum.CI)*2,to=2,by=-2)))
    legend('topleft',c('insitu.hot','insitu.realm'),col=c('red','blue'),lty=1,cex=.7,bty='n')
  }else if (plots=='dispersal.all.rates'){
    plot(c(1,1),xlim=c(0,length(dispersal.into.realm.rates.BSM.CI)),ylim=c(0,max(c(unlist(dispersal.into.realm.rates.BSM.CI),unlist(dispersal.into.hot.rates.BSM.CI))[is.finite(c(unlist(dispersal.into.realm.rates.BSM.CI),unlist(dispersal.into.hot.rates.BSM.CI)))])),type='n',xaxt='n',xlab='age',ylab='dispersal rates',main=paste(name,' dispersal rates',sep=''))
    lines(x=c(1:length(dispersal.into.realm.rates.BSM.CI)),y=unlist(lapply(dispersal.into.realm.rates.BSM.CI,function(x)x[3])),col=adjustcolor( "red", alpha.f = 0.9))
    polygon(x=c((1:length(dispersal.into.realm.rates.BSM.CI)),(rev(1:length(dispersal.into.realm.rates.BSM.CI)))),y=c((unlist(lapply(dispersal.into.realm.rates.BSM.CI,function(x)x[1]))),rev(unlist(lapply(dispersal.into.realm.rates.BSM.CI,function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
    lines(x=c(1:length(dispersal.into.hot.rates.BSM.CI)),y=unlist(lapply(dispersal.into.hot.rates.BSM.CI,function(x)x[3])),col=adjustcolor( "blue", alpha.f = 0.9))
    polygon(x=c((1:length(dispersal.into.hot.rates.BSM.CI)),(rev(1:length(dispersal.into.hot.rates.BSM.CI)))),y=c((unlist(lapply(dispersal.into.hot.rates.BSM.CI,function(x)x[1]))),rev(unlist(lapply(dispersal.into.hot.rates.BSM.CI,function(x)x[2])))),col=adjustcolor( "blue", alpha.f = 0.2))
    axis(1,at=c(1:length(dispersal.into.hot.rates.BSM.CI)),labels=c(seq(from=length(dispersal.into.hot.rates.BSM.CI)*2,to=2,by=-2)))
    
    legend('topleft',c('dispersal.into.hot','dispersal.into.realm'),col=c('blue','red'),lty=1,cex=.7,bty='n')
  }else if(plots=='dispersal.all.events'){
    plot(c(1,1),xlim=c(0,length(dispersal.into.realm.time.bins.counts.rev.cumsum.CI)),ylim=c(0,max(c(unlist(dispersal.into.realm.time.bins.counts.rev.cumsum.CI),unlist(dispersal.into.hot.time.bins.counts.rev.cumsum.CI))[is.finite(c(unlist(dispersal.into.realm.time.bins.counts.rev.cumsum.CI),unlist(dispersal.into.hot.time.bins.counts.rev.cumsum.CI)))])),type='n',xaxt='n',xlab='age',ylab='cumulative # of events',main=paste(name,' dispersal events',sep=''))
    lines(x=c(1:length(dispersal.into.realm.time.bins.counts.rev.cumsum.CI)),y=unlist(lapply(dispersal.into.realm.time.bins.counts.rev.cumsum.CI,function(x)x[3])),col=adjustcolor( "red", alpha.f = 0.9))
    polygon(x=c((1:length(dispersal.into.realm.time.bins.counts.rev.cumsum.CI)),(rev(1:length(dispersal.into.realm.time.bins.counts.rev.cumsum.CI)))),y=c((unlist(lapply(dispersal.into.realm.time.bins.counts.rev.cumsum.CI,function(x)x[1]))),rev(unlist(lapply(dispersal.into.realm.time.bins.counts.rev.cumsum.CI,function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
    lines(x=c(1:length(dispersal.into.hot.time.bins.counts.rev.cumsum.CI)),y=unlist(lapply(dispersal.into.hot.time.bins.counts.rev.cumsum.CI,function(x)x[3])),col=adjustcolor( "blue", alpha.f = 0.9))
    polygon(x=c((1:length(dispersal.into.hot.time.bins.counts.rev.cumsum.CI)),(rev(1:length(dispersal.into.hot.time.bins.counts.rev.cumsum.CI)))),y=c((unlist(lapply(dispersal.into.hot.time.bins.counts.rev.cumsum.CI,function(x)x[1]))),rev(unlist(lapply(dispersal.into.hot.time.bins.counts.rev.cumsum.CI,function(x)x[2])))),col=adjustcolor( "blue", alpha.f = 0.2))
    axis(1,at=c(1:length(dispersal.into.hot.time.bins.counts.rev.cumsum.CI)),labels=c(seq(from=length(dispersal.into.hot.time.bins.counts.rev.cumsum.CI)*2,to=2,by=-2)))
    legend('topleft',c('dispersal.into.realm','dispersal.into.hot'),col=c('red','blue'),lty=1,cex=.7,bty='n')
  }else if(plots=='dispersal.into.from.events'){
    plot(c(1,1),xlim=c(0,length(dispersal.into.realm.time.bins.counts.rev.cumsum.CI)),ylim=c(0,max(c(unlist(dispersal.into.realm.time.bins.counts.rev.cumsum.CI),unlist(dispersal.into.hot.time.bins.counts.rev.cumsum.CI),unlist(dispersal.from.hot.time.bins.counts.rev.cumsum.CI),unlist(dispersal.from.realm.time.bins.counts.rev.cumsum.CI)))),type='n',xaxt='n',xlab='age',ylab='cumulative # of events',main=paste(name,' dispersal events',sep=''))
    lines(x=c(1:length(dispersal.into.realm.time.bins.counts.rev.cumsum.CI)),y=unlist(lapply(dispersal.into.realm.time.bins.counts.rev.cumsum.CI,function(x)x[3])),col=adjustcolor( "green", alpha.f = 0.9))
    polygon(x=c((1:length(dispersal.into.realm.time.bins.counts.rev.cumsum.CI)),(rev(1:length(dispersal.into.realm.time.bins.counts.rev.cumsum.CI)))),y=c((unlist(lapply(dispersal.into.realm.time.bins.counts.rev.cumsum.CI,function(x)x[1]))),rev(unlist(lapply(dispersal.into.realm.time.bins.counts.rev.cumsum.CI,function(x)x[2])))),col=adjustcolor( "green", alpha.f = 0.2))
    lines(x=c(1:length(dispersal.into.hot.time.bins.counts.rev.cumsum.CI)),y=unlist(lapply(dispersal.into.hot.time.bins.counts.rev.cumsum.CI,function(x)x[3])),col=adjustcolor( "orange", alpha.f = 0.9))
    polygon(x=c((1:length(dispersal.into.hot.time.bins.counts.rev.cumsum.CI)),(rev(1:length(dispersal.into.hot.time.bins.counts.rev.cumsum.CI)))),y=c((unlist(lapply(dispersal.into.hot.time.bins.counts.rev.cumsum.CI,function(x)x[1]))),rev(unlist(lapply(dispersal.into.hot.time.bins.counts.rev.cumsum.CI,function(x)x[2])))),col=adjustcolor( "orange", alpha.f = 0.2))
    lines(x=c(1:length(dispersal.from.realm.time.bins.counts.rev.cumsum.CI)),y=unlist(lapply(dispersal.from.realm.time.bins.counts.rev.cumsum.CI,function(x)x[3])),col=adjustcolor( "blue", alpha.f = 0.9))
    polygon(x=c((1:length(dispersal.from.realm.time.bins.counts.rev.cumsum.CI)),(rev(1:length(dispersal.from.realm.time.bins.counts.rev.cumsum.CI)))),y=c((unlist(lapply(dispersal.from.realm.time.bins.counts.rev.cumsum.CI,function(x)x[1]))),rev(unlist(lapply(dispersal.from.realm.time.bins.counts.rev.cumsum.CI,function(x)x[2])))),col=adjustcolor( "blue", alpha.f = 0.2))
    lines(x=c(1:length(dispersal.from.hot.time.bins.counts.rev.cumsum.CI)),y=unlist(lapply(dispersal.from.hot.time.bins.counts.rev.cumsum.CI,function(x)x[3])),col=adjustcolor( "red", alpha.f = 0.9))
    polygon(x=c((1:length(dispersal.from.hot.time.bins.counts.rev.cumsum.CI)),(rev(1:length(dispersal.from.hot.time.bins.counts.rev.cumsum.CI)))),y=c((unlist(lapply(dispersal.from.hot.time.bins.counts.rev.cumsum.CI,function(x)x[1]))),rev(unlist(lapply(dispersal.from.hot.time.bins.counts.rev.cumsum.CI,function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
    axis(1,at=c(1:length(dispersal.from.hot.time.bins.counts.rev.cumsum.CI)),labels=c(seq(from=length(dispersal.from.hot.time.bins.counts.rev.cumsum.CI)*2,to=2,by=-2)))
    legend('topleft',c('dispersal.from.realm','dispersal.from.hot','dispersal.into.realm','dispersal.into.hot'),col=c('blue','red','green','orange'),lty=1,cex=.7,bty='n')
  }
  
  #plot BSM dispersal rates hot into realm
  
  ##plot BAMM in situ rates
  #par(new=T)
  #plot(c(1,1),xlim=c(0,length(insitu.hot.rates.BAMM.CI)),ylim=c(0,max(c(unlist(insitu.hot.rates.BAMM.CI),unlist(insitu.realm.rates.BAMM.CI)))),type='n',xaxt='n',xlab='age',yaxt='n')
  ##polygon(x=c(c(age:length(insitu.hot.rates.BAMM.CI)),c(length(insitu.hot.rates.BAMM.CI):age)),y=c(unlist(lapply(insitu.hot.rates.BAMM.CI[c(age:length(insitu.hot.rates.BAMM.CI))],function(x)x[1])),rev(unlist(lapply(insitu.hot.rates.BAMM.CI[c(age:length(insitu.hot.rates.BAMM.CI))],function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
  #lines(x=c(1:length(insitu.hot.rates.BAMM.CI)),y=unlist(lapply(insitu.hot.rates.BAMM.CI,function(x)x[3])),col=adjustcolor( "red", alpha.f = 0.9))
  #polygon(x=c((1:length(insitu.hot.rates.BAMM.CI)),(rev(1:length(insitu.hot.rates.BAMM.CI)))),y=c((unlist(lapply(insitu.hot.rates.BAMM.CI,function(x)x[1]))),rev(unlist(lapply(insitu.hot.rates.BAMM.CI,function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
  #lines(x=c(1:length(insitu.realm.rates.BAMM.CI)),y=unlist(lapply(insitu.realm.rates.BAMM.CI,function(x)x[3])),col=adjustcolor( "blue", alpha.f = 0.9))
  #polygon(x=c((1:length(insitu.realm.rates.BAMM.CI)),(rev(1:length(insitu.realm.rates.BAMM.CI)))),y=c((unlist(lapply(insitu.realm.rates.BAMM.CI,function(x)x[1]))),rev(unlist(lapply(insitu.realm.rates.BAMM.CI,function(x)x[2])))),col=adjustcolor( "blue", alpha.f = 0.2))
  ##axis(1,at=c(0:length(BSMoutput.real.rates[['insitu.realm.rates.BSM.CI']])),labels=c((length(BSMoutput.real.rates[['insitu.realm.rates.BSM.CI']])+2):2))
  #axis(4)
  #legend('topright',c('insitu.hot','insitu.realm'),col=c('red','blue'),lty=1,cex=.7,bty='n')
}


plot_rates_events_real_2mybin_difference_both_scaled_fixedyaxis<-function(BSMoutput.real.rates.mammals,BSMoutput.real.events.mammals,BSMoutput.real.rates.birds,BSMoutput.real.events.birds,age,name,plots,ylim){
  #correct age
  age.mammals<-(length(BSMoutput.real.rates.mammals[['insitu.realm.rates.BSM.CI']])*2-age)/2
  #trim lists to age period
  BSMoutput.real.rates.mammals<-lapply(BSMoutput.real.rates.mammals,function(x)x[age.mammals:length(x)])
  BSMoutput.real.events.mammals<-lapply(BSMoutput.real.events.mammals,function(x)x[age.mammals:length(x)])
  #delete Inf in rates and events
  BSMoutput.real.rates.mammals<-lapply(BSMoutput.real.rates.mammals,function(x) lapply(x,function(x) {x[!is.finite(x)]<-NA;return(x)}))
  BSMoutput.real.events.mammals<-lapply(BSMoutput.real.events.mammals,function(x) lapply(x,function(x) {x[!is.finite(x)]<-NA;return(x)}))
  #remove the last bin - two million year time bin
  BSMoutput.real.rates.mammals<-lapply(BSMoutput.real.rates.mammals,function(x) x<-x[-c(length(x))])
  BSMoutput.real.events.mammals<-lapply(BSMoutput.real.events.mammals,function(x) x<-x[-c((length(x)-1),length(x))])
  #assign variables for easier handling
  hot.into.realm.rates.BSM.mammals<-BSMoutput.real.rates.mammals[['hot.into.realm.rates.BSM.CI']]
  realm.into.hot.rates.BSM.mammals<-BSMoutput.real.rates.mammals[['realm.into.hot.rates.BSM.CI']]
  insitu.realm.rates.BAMM.mammals<-BSMoutput.real.rates.mammals[['insitu.realm.rates.BAMM.CI']]
  insitu.hot.rates.BAMM.mammals<-BSMoutput.real.rates.mammals[['insitu.hot.rates.BAMM.CI']]
  insitu.realm.rates.BSM.mammals<-BSMoutput.real.rates.mammals[['insitu.realm.rates.BSM.CI']]
  insitu.hot.rates.BSM.mammals<-BSMoutput.real.rates.mammals[['insitu.hot.rates.BSM.CI']]
  dispersal.into.realm.rates.BSM.mammals<-BSMoutput.real.rates.mammals[['dispersal.into.realm.rates.BSM.CI']]
  dispersal.into.hot.rates.BSM.mammals<-BSMoutput.real.rates.mammals[['dispersal.into.hot.rates.BSM.CI']]
  dispersal.from.realm.rates.BSM.mammals<-BSMoutput.real.rates.mammals[['dispersal.from.realm.rates.BSM.CI']]
  dispersal.from.hot.rates.BSM.mammals<-BSMoutput.real.rates.mammals[['dispersal.from.hot.rates.BSM.CI']]
  
  hot.into.realm.time.bins.counts.rev.cumsum.mammals<-BSMoutput.real.events.mammals[['hot.into.realm.time.bins.counts.rev.cumsum.CI']]
  realm.into.hot.time.bins.counts.rev.cumsum.mammals<-BSMoutput.real.events.mammals[['realm.into.hot.time.bins.counts.rev.cumsum.CI']]
  insitu.realm.time.bins.counts.rev.cumsum.mammals<-BSMoutput.real.events.mammals[['insitu.realm.time.bins.counts.rev.cumsum.CI']]
  insitu.hot.time.bins.counts.rev.cumsum.mammals<-BSMoutput.real.events.mammals[['insitu.hot.time.bins.counts.rev.cumsum.CI']]
  dispersal.into.realm.time.bins.counts.rev.cumsum.mammals<-BSMoutput.real.events.mammals[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']]
  dispersal.into.hot.time.bins.counts.rev.cumsum.mammals<-BSMoutput.real.events.mammals[['dispersal.into.hot.time.bins.counts.rev.cumsum.CI']]
  dispersal.from.realm.time.bins.counts.rev.cumsum.mammals<-BSMoutput.real.events.mammals[['dispersal.from.realm.time.bins.counts.rev.cumsum.CI']]
  dispersal.from.hot.time.bins.counts.rev.cumsum.mammals<-BSMoutput.real.events.mammals[['dispersal.from.hot.time.bins.counts.rev.cumsum.CI']]
  
  age.birds<-(length(BSMoutput.real.rates.birds[['insitu.realm.rates.BSM.CI']])*2-age)/2
  #trim lists to age period
  BSMoutput.real.rates.birds<-lapply(BSMoutput.real.rates.birds,function(x)x[age.birds:length(x)])
  BSMoutput.real.events.birds<-lapply(BSMoutput.real.events.birds,function(x)x[age.birds:length(x)])
  #delete Inf in rates and events
  BSMoutput.real.rates.birds<-lapply(BSMoutput.real.rates.birds,function(x) lapply(x,function(x) {x[!is.finite(x)]<-NA;return(x)}))
  BSMoutput.real.events.birds<-lapply(BSMoutput.real.events.birds,function(x) lapply(x,function(x) {x[!is.finite(x)]<-NA;return(x)}))
  #remove the last bin - two million year time bin
  BSMoutput.real.rates.birds<-lapply(BSMoutput.real.rates.birds,function(x) x<-x[-c(length(x))])
  BSMoutput.real.events.birds<-lapply(BSMoutput.real.events.birds,function(x) x<-x[-c((length(x)-1),length(x))])
  #assign variables for easier handling
  hot.into.realm.rates.BSM.birds<-BSMoutput.real.rates.birds[['hot.into.realm.rates.BSM.CI']]
  realm.into.hot.rates.BSM.birds<-BSMoutput.real.rates.birds[['realm.into.hot.rates.BSM.CI']]
  insitu.realm.rates.BAMM.birds<-BSMoutput.real.rates.birds[['insitu.realm.rates.BAMM.CI']]
  insitu.hot.rates.BAMM.birds<-BSMoutput.real.rates.birds[['insitu.hot.rates.BAMM.CI']]
  insitu.realm.rates.BSM.birds<-BSMoutput.real.rates.birds[['insitu.realm.rates.BSM.CI']]
  insitu.hot.rates.BSM.birds<-BSMoutput.real.rates.birds[['insitu.hot.rates.BSM.CI']]
  dispersal.into.realm.rates.BSM.birds<-BSMoutput.real.rates.birds[['dispersal.into.realm.rates.BSM.CI']]
  dispersal.into.hot.rates.BSM.birds<-BSMoutput.real.rates.birds[['dispersal.into.hot.rates.BSM.CI']]
  dispersal.from.realm.rates.BSM.birds<-BSMoutput.real.rates.birds[['dispersal.from.realm.rates.BSM.CI']]
  dispersal.from.hot.rates.BSM.birds<-BSMoutput.real.rates.birds[['dispersal.from.hot.rates.BSM.CI']]
  
  hot.into.realm.time.bins.counts.rev.cumsum.birds<-BSMoutput.real.events.birds[['hot.into.realm.time.bins.counts.rev.cumsum.CI']]
  realm.into.hot.time.bins.counts.rev.cumsum.birds<-BSMoutput.real.events.birds[['realm.into.hot.time.bins.counts.rev.cumsum.CI']]
  insitu.realm.time.bins.counts.rev.cumsum.birds<-BSMoutput.real.events.birds[['insitu.realm.time.bins.counts.rev.cumsum.CI']]
  insitu.hot.time.bins.counts.rev.cumsum.birds<-BSMoutput.real.events.birds[['insitu.hot.time.bins.counts.rev.cumsum.CI']]
  dispersal.into.realm.time.bins.counts.rev.cumsum.birds<-BSMoutput.real.events.birds[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']]
  dispersal.into.hot.time.bins.counts.rev.cumsum.birds<-BSMoutput.real.events.birds[['dispersal.into.hot.time.bins.counts.rev.cumsum.CI']]
  dispersal.from.realm.time.bins.counts.rev.cumsum.birds<-BSMoutput.real.events.birds[['dispersal.from.realm.time.bins.counts.rev.cumsum.CI']]
  dispersal.from.hot.time.bins.counts.rev.cumsum.birds<-BSMoutput.real.events.birds[['dispersal.from.hot.time.bins.counts.rev.cumsum.CI']]
  
  if(plots=='within.realm.dispersal.rates'){
    #difference.into.hotminusrealm<-lapply(c(1:length(dispersal.into.realm.rates.BSM)),function(x)sample(dispersal.into.hot.rates.BSM[[x]],size = 100,replace=TRUE)-sample(dispersal.into.realm.rates.BSM[[x]],size = 100,replace=TRUE))
    difference.hotintorealm.realmintohot.mammals<-lapply(c(1:length(hot.into.realm.rates.BSM.mammals)),function(x)sample(hot.into.realm.rates.BSM.mammals[[x]],size = 50,replace=FALSE)-sample(realm.into.hot.rates.BSM.mammals[[x]],size = 50,replace=FALSE))
    difference.hotintorealm.realmintohot.birds<-lapply(c(1:length(hot.into.realm.rates.BSM.birds)),function(x)sample(hot.into.realm.rates.BSM.birds[[x]],size = 50,replace=FALSE)-sample(realm.into.hot.rates.BSM.birds[[x]],size = 50,replace=FALSE))
    #global.mean<-mean(c(c(unlist(difference.hotintorealm.realmintohot.mammals),unlist(difference.hotintorealm.realmintohot.birds))),na.rm=TRUE)
    global.sd<-sd(c(c(unlist(difference.hotintorealm.realmintohot.mammals),unlist(difference.hotintorealm.realmintohot.birds))),na.rm=TRUE)
    #difference.hotintorealm.realmintohot.mammals<-lapply(difference.hotintorealm.realmintohot.mammals,function(x)(x-global.mean)/global.sd)
    #difference.hotintorealm.realmintohot.birds<-lapply(difference.hotintorealm.realmintohot.birds,function(x)(x-global.mean)/global.sd)
    difference.hotintorealm.realmintohot.mammals<-lapply(difference.hotintorealm.realmintohot.mammals,function(x)x/global.sd)
    difference.hotintorealm.realmintohot.birds<-lapply(difference.hotintorealm.realmintohot.birds,function(x)x/global.sd)
    
    median.mammals<-unlist(lapply(difference.hotintorealm.realmintohot.mammals,function(x)median(x,na.rm=TRUE)))
    lower.mammals<-unlist(lapply(difference.hotintorealm.realmintohot.mammals,function(x)quantile(x,c(0.1),na.rm=TRUE)))
    upper.mammals<-unlist(lapply(difference.hotintorealm.realmintohot.mammals,function(x)quantile(x,c(0.9),na.rm=TRUE)))
    median.birds<-unlist(lapply(difference.hotintorealm.realmintohot.birds,function(x)median(x,na.rm=TRUE)))
    lower.birds<-unlist(lapply(difference.hotintorealm.realmintohot.birds,function(x)quantile(x,c(0.1),na.rm=TRUE)))
    upper.birds<-unlist(lapply(difference.hotintorealm.realmintohot.birds,function(x)quantile(x,c(0.9),na.rm=TRUE)))
    #lower<-unlist(lapply(difference.hotintorealm.realmintohot,function(x)quantile(x,c(0.025),na.rm=TRUE)))
    #upper<-unlist(lapply(difference.hotintorealm.realmintohot,function(x)quantile(x,c(0.975),na.rm=TRUE)))
    plot(x=c(1:length(median.mammals)),y=median.mammals,ylim=c(-ylim,ylim),type='n',xaxt='n',yaxt='n',xlab='age',ylab='diff. dispersal hot.to.realm - realm.to.hot',main=paste(name,' diff in dispersal rates',sep=''))
    #plot(x=c(1:length(median.mammals)),y=median.mammals,ylim=c(-max(abs(c(median.mammals,upper.mammals,lower.mammals,median.birds,upper.birds,lower.birds))),max(abs(c(median.mammals,upper.mammals,lower.mammals,median.birds,upper.birds,lower.birds)))),type='n',xaxt='n',yaxt='n',xlab='age',ylab='diff. dispersal hot.to.realm - realm.to.hot',main=paste(name,' diff in dispersal rates',sep=''))
    abline(h=(0/global.sd),lty=2)
    #abline(h=(0-global.mean)/global.sd,lty=2)
    polygon(x=c(c(1:length(median.mammals),length(median.mammals):1)),y=c(lower.mammals,rev(upper.mammals)),col=adjustcolor(col='#009999',alpha.f=0.25))
    lines(x=c(1:length(median.mammals)),y=median.mammals,col='#009999') 
    polygon(x=c(c(1:length(median.birds),length(median.birds):1)),y=c(lower.birds,rev(upper.birds)),col=adjustcolor(col="#FF7400",alpha.f=0.25))
    lines(x=c(1:length(median.birds)),y=median.birds,col="#FF7400") 
    axis(1,at=c(1:length(difference.hotintorealm.realmintohot.mammals)),labels=c(seq(from=length(difference.hotintorealm.realmintohot.mammals)*2,to=2,by=-2)))
    axis(2,las=2)
    #plot(c(1,1),xlim=c(0,length(hot.into.realm.rates.BSM.CI)),ylim=c(0,max(c(unlist(hot.into.realm.rates.BSM.CI),unlist(realm.into.hot.rates.BSM.CI))[is.finite(c(unlist(hot.into.realm.rates.BSM.CI),unlist(realm.into.hot.rates.BSM.CI)))])),type='n',xaxt='n',xlab='age',ylab='dispersal rates',main=paste(name,' dispersal rates',sep=''))
    #lines(x=c(1:length(hot.into.realm.rates.BSM.CI)),y=unlist(lapply(hot.into.realm.rates.BSM.CI,function(x)x[3])),col=adjustcolor( "red", alpha.f = 0.9))
    #polygon(x=c((1:length(hot.into.realm.rates.BSM.CI)),(rev(1:length(hot.into.realm.rates.BSM.CI)))),y=c((unlist(lapply(hot.into.realm.rates.BSM.CI,function(x)x[1]))),rev(unlist(lapply(hot.into.realm.rates.BSM.CI,function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
    #lines(x=c(1:length(realm.into.hot.rates.BSM.CI)),y=unlist(lapply(realm.into.hot.rates.BSM.CI,function(x)x[3])),col=adjustcolor( "blue", alpha.f = 0.9))
    #polygon(x=c((1:length(realm.into.hot.rates.BSM.CI)),(rev(1:length(realm.into.hot.rates.BSM.CI)))),y=c((unlist(lapply(realm.into.hot.rates.BSM.CI,function(x)x[1]))),rev(unlist(lapply(realm.into.hot.rates.BSM.CI,function(x)x[2])))),col=adjustcolor( "blue", alpha.f = 0.2))
    #
    legend('topleft',c('mammal.diff','bird.diff'),col=c('#009999',"#FF7400"),lty=1,cex=.7,bty='n')
  }else if(plots=='within.realm.cladogenesis.BAMM.rates'){
    difference.insituhot.insiturealm.mammals<-lapply(c(1:length(insitu.hot.rates.BAMM.mammals)),function(x)sample(insitu.hot.rates.BAMM.mammals[[x]],size = 100,replace=TRUE)-sample(insitu.realm.rates.BAMM.mammals[[x]],size = 100,replace=TRUE))
    difference.insituhot.insiturealm.birds<-lapply(c(1:length(insitu.hot.rates.BAMM.birds)),function(x)sample(insitu.hot.rates.BAMM.birds[[x]],size = 100,replace=TRUE)-sample(insitu.realm.rates.BAMM.birds[[x]],size = 100,replace=TRUE))
    global.sd<-sd(c(c(unlist(difference.insituhot.insiturealm.mammals),unlist(difference.insituhot.insiturealm.birds))),na.rm=TRUE)
    difference.insituhot.insiturealm.mammals<-lapply(difference.insituhot.insiturealm.mammals,function(x)x/global.sd)
    difference.insituhot.insiturealm.birds<-lapply(difference.insituhot.insiturealm.birds,function(x)x/global.sd)
    median.mammals<-unlist(lapply(difference.insituhot.insiturealm.mammals,function(x)median(x,na.rm=TRUE)))
    lower.mammals<-unlist(lapply(difference.insituhot.insiturealm.mammals,function(x)quantile(x,c(0.1),na.rm=TRUE)))
    upper.mammals<-unlist(lapply(difference.insituhot.insiturealm.mammals,function(x)quantile(x,c(0.9),na.rm=TRUE)))
    median.birds<-unlist(lapply(difference.insituhot.insiturealm.birds,function(x)median(x,na.rm=TRUE)))
    lower.birds<-unlist(lapply(difference.insituhot.insiturealm.birds,function(x)quantile(x,c(0.1),na.rm=TRUE)))
    upper.birds<-unlist(lapply(difference.insituhot.insiturealm.birds,function(x)quantile(x,c(0.9),na.rm=TRUE)))
    #lower<-unlist(lapply(difference.insituhot.insiturealm,function(x)quantile(x,c(0.025),na.rm=TRUE)))
    #upper<-unlist(lapply(difference.insituhot.insiturealm,function(x)quantile(x,c(0.975),na.rm=TRUE)))
    plot(x=c(1:length(median.mammals)),y=median.mammals,ylim=c(-ylim,ylim),type='n',xaxt='n',yaxt='n',xlab='age',ylab='diff. insitu.hot - insitu.realm',main=paste(name,' diff in dispersal rates',sep=''))
    #plot(x=c(1:length(median.mammals)),y=median.mammals,ylim=c(-max(abs(c(median.mammals,upper.mammals,lower.mammals,median.birds,upper.birds,lower.birds))),max(abs(c(median.mammals,upper.mammals,lower.mammals,median.birds,upper.birds,lower.birds)))),type='n',xaxt='n',yaxt='n',xlab='age',ylab='diff. dispersal hot.to.realm - realm.to.hot',main=paste(name,' diff in dispersal rates',sep=''))
    abline(h=(0/global.sd),lty=2)
    #abline(h=(0-global.mean)/global.sd,lty=2)
    polygon(x=c(c(1:length(median.mammals),length(median.mammals):1)),y=c(lower.mammals,rev(upper.mammals)),col=adjustcolor(col='#009999',alpha.f=0.25))
    lines(x=c(1:length(median.mammals)),y=median.mammals,col='#009999') 
    polygon(x=c(c(1:length(median.birds),length(median.birds):1)),y=c(lower.birds,rev(upper.birds)),col=adjustcolor(col="#FF7400",alpha.f=0.25))
    lines(x=c(1:length(median.birds)),y=median.birds,col="#FF7400") 
    axis(1,at=c(1:length(difference.insituhot.insiturealm.mammals)),labels=c(seq(from=length(difference.insituhot.insiturealm.mammals)*2,to=2,by=-2)))
    axis(2,las=2)
    
    #plot(c(1,1),xlim=c(0,length(insitu.hot.rates.BAMM.CI)),ylim=c(0,max(c(unlist(insitu.hot.rates.BAMM.CI),unlist(insitu.realm.rates.BAMM.CI))[is.finite(c(unlist(insitu.hot.rates.BAMM.CI),unlist(insitu.realm.rates.BAMM.CI)))])),type='n',xaxt='n',xlab='age',ylab='dispersal rates',main=paste(name,' dispersal rates',sep=''))
    #lines(x=c(1:length(insitu.hot.rates.BAMM.CI)),y=unlist(lapply(insitu.hot.rates.BAMM.CI,function(x)x[3])),col=adjustcolor( "red", alpha.f = 0.9))
    #polygon(x=c((1:length(insitu.hot.rates.BAMM.CI)),(rev(1:length(insitu.hot.rates.BAMM.CI)))),y=c((unlist(lapply(insitu.hot.rates.BAMM.CI,function(x)x[1]))),rev(unlist(lapply(insitu.hot.rates.BAMM.CI,function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
    #lines(x=c(1:length(insitu.realm.rates.BAMM.CI)),y=unlist(lapply(insitu.realm.rates.BAMM.CI,function(x)x[3])),col=adjustcolor( "blue", alpha.f = 0.9))
    #polygon(x=c((1:length(insitu.realm.rates.BAMM.CI)),(rev(1:length(insitu.realm.rates.BAMM.CI)))),y=c((unlist(lapply(insitu.realm.rates.BAMM.CI,function(x)x[1]))),rev(unlist(lapply(insitu.realm.rates.BAMM.CI,function(x)x[2])))),col=adjustcolor( "blue", alpha.f = 0.2))
    #
    legend('topleft',c('mammal.diff','bird.diff'),col=c('#009999',"#FF7400"),lty=1,cex=.7,bty='n')
  }else if(plots=='within.realm.cladogenesis.BSM.rates'){
    difference.insituhot.insiturealm.mammals<-lapply(c(1:length(insitu.hot.rates.BSM.mammals)),function(x)sample(insitu.hot.rates.BSM.mammals[[x]],size = 100,replace=TRUE)-sample(insitu.realm.rates.BSM.mammals[[x]],size = 100,replace=TRUE))
    difference.insituhot.insiturealm.birds<-lapply(c(1:length(insitu.hot.rates.BSM.birds)),function(x)sample(insitu.hot.rates.BSM.birds[[x]],size = 100,replace=TRUE)-sample(insitu.realm.rates.BSM.birds[[x]],size = 100,replace=TRUE))
    median.mammals<-unlist(lapply(difference.insituhot.insiturealm.mammals,function(x)median(x,na.rm=TRUE)))
    lower.mammals<-unlist(lapply(difference.insituhot.insiturealm.mammals,function(x)quantile(x,c(0.1),na.rm=TRUE)))
    upper.mammals<-unlist(lapply(difference.insituhot.insiturealm.mammals,function(x)quantile(x,c(0.9),na.rm=TRUE)))
    median.birds<-unlist(lapply(difference.insituhot.insiturealm.birds,function(x)median(x,na.rm=TRUE)))
    lower.birds<-unlist(lapply(difference.insituhot.insiturealm.birds,function(x)quantile(x,c(0.1),na.rm=TRUE)))
    upper.birds<-unlist(lapply(difference.insituhot.insiturealm.birds,function(x)quantile(x,c(0.9),na.rm=TRUE)))
    #lower<-unlist(lapply(difference.insituhot.insiturealm,function(x)quantile(x,c(0.025),na.rm=TRUE)))
    #upper<-unlist(lapply(difference.insituhot.insiturealm,function(x)quantile(x,c(0.975),na.rm=TRUE)))
    plot(x=c(1:length(median.mammals)),y=median.mammals,ylim=c(-max(abs(c(median.mammals,upper.mammals,lower.mammals,median.birds,upper.birds,lower.birds))),max(abs(c(median.mammals,upper.mammals,lower.mammals,median.birds,upper.birds,lower.birds)))),type='n',xaxt='n',yaxt='n',xlab='age',ylab='diff. insitu.hot - insitu.realm BSM',main=paste(name,' diff in insitu.BSM.rates',sep=''))
    abline(h=0,lty=2)
    polygon(x=c(c(1:length(median.mammals),length(median.mammals):1)),y=c(lower.mammals,rev(upper.mammals)),col=adjustcolor(col='#009999',alpha.f=0.25))
    lines(x=c(1:length(median.mammals)),y=median.mammals,col='#009999') 
    polygon(x=c(c(1:length(median.birds),length(median.birds):1)),y=c(lower.birds,rev(upper.birds)),col=adjustcolor(col="#FF7400",alpha.f=0.25))
    lines(x=c(1:length(median.birds)),y=median.birds,col="#FF7400") 
    axis(1,at=c(1:length(difference.insituhot.insiturealm.mammals)),labels=c(seq(from=length(difference.insituhot.insiturealm.mammals)*2,to=2,by=-2)))
    axis(2,las=2)
    #plot(c(1,1),xlim=c(0,length(insitu.hot.rates.BSM.CI)),ylim=c(0,max(c(unlist(insitu.hot.rates.BSM.CI),unlist(insitu.realm.rates.BSM.CI))[is.finite(c(unlist(insitu.hot.rates.BSM.CI),unlist(insitu.realm.rates.BSM.CI)))])),type='n',xaxt='n',xlab='age',ylab='dispersal rates',main=paste(name,' dispersal rates',sep=''))
    #lines(x=c(1:length(insitu.hot.rates.BSM.CI)),y=unlist(lapply(insitu.hot.rates.BSM.CI,function(x)x[3])),col=adjustcolor( "red", alpha.f = 0.9))
    #polygon(x=c((1:length(insitu.hot.rates.BSM.CI)),(rev(1:length(insitu.hot.rates.BSM.CI)))),y=c((unlist(lapply(insitu.hot.rates.BSM.CI,function(x)x[1]))),rev(unlist(lapply(insitu.hot.rates.BSM.CI,function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
    #lines(x=c(1:length(insitu.realm.rates.BSM.CI)),y=unlist(lapply(insitu.realm.rates.BSM.CI,function(x)x[3])),col=adjustcolor( "blue", alpha.f = 0.9))
    #polygon(x=c((1:length(insitu.realm.rates.BSM.CI)),(rev(1:length(insitu.realm.rates.BSM.CI)))),y=c((unlist(lapply(insitu.realm.rates.BSM.CI,function(x)x[1]))),rev(unlist(lapply(insitu.realm.rates.BSM.CI,function(x)x[2])))),col=adjustcolor( "blue", alpha.f = 0.2))
    #
    legend('topleft',c('mammal.diff','bird.diff'),col=c('#009999',"#FF7400"),lty=1,cex=.7,bty='n')
  }else if(plots=='global.dispersal.rates'){
    #difference.into.hotminusrealm<-lapply(c(1:length(dispersal.into.realm.rates.BSM)),function(x)sample(dispersal.into.hot.rates.BSM[[x]],size = 100,replace=TRUE)-sample(dispersal.into.realm.rates.BSM[[x]],size = 100,replace=TRUE))
    difference.hotrealm.from<-lapply(c(1:length(dispersal.from.realm.rates.BSM)),function(x)sample(dispersal.from.hot.rates.BSM[[x]],size = 100,replace=TRUE)-sample(dispersal.from.realm.rates.BSM[[x]],size = 100,replace=TRUE))
    difference.hotrealm.into<-lapply(c(1:length(dispersal.into.realm.rates.BSM)),function(x)sample(dispersal.into.hot.rates.BSM[[x]],size = 100,replace=TRUE)-sample(dispersal.into.realm.rates.BSM[[x]],size = 100,replace=TRUE))
    median.from<-unlist(lapply(difference.hotrealm.from,function(x)median(x,na.rm=TRUE)))
    lower.from<-unlist(lapply(difference.hotrealm.from,function(x)quantile(x,c(0.1),na.rm=TRUE)))
    upper.from<-unlist(lapply(difference.hotrealm.from,function(x)quantile(x,c(0.9),na.rm=TRUE)))
    median.into<-unlist(lapply(difference.hotrealm.into,function(x)median(x,na.rm=TRUE)))
    lower.into<-unlist(lapply(difference.hotrealm.into,function(x)quantile(x,c(0.1),na.rm=TRUE)))
    upper.into<-unlist(lapply(difference.hotrealm.into,function(x)quantile(x,c(0.9),na.rm=TRUE)))
    
    plot(x=c(1:length(median.from)),y=median.from,ylim=c(-max(abs(c(median.from,lower.from,upper.from,median.into,lower.into,upper.into))),max(abs(c(median.from,lower.from,upper.from,median.into,lower.into,upper.into)))),type='n',xaxt='n',xlab='age',ylab='diff. dispersal hot - realm',main=paste(name,' diff in global dispersal rates',sep=''))
    polygon(x=c(c(1:length(median.from),length(median.from):1)),y=c(lower.from,rev(upper.from)),col=adjustcolor('red',alpha.f = 0.25))
    lines(x=c(1:length(median.from)),y=median.from,col='red') 
    abline(h=0,lty=2)
    polygon(x=c(c(1:length(median.into),length(median.into):1)),y=c(lower.into,rev(upper.into)),col=adjustcolor('blue',alpha.f = 0.25))
    lines(x=c(1:length(median.into)),y=median.into,col='blue') 
    axis(1,at=c(1:length(difference.hotintorealm.realmintohot)),labels=c(seq(from=length(difference.hotintorealm.realmintohot)*2,to=2,by=-2)))
    legend('topleft',c('dispersal.from.diff','dispersal.into.diff'),col=c('red','blue'),lty=1,cex=.7,bty='n')
    #plot(c(1,1),xlim=c(0,length(hot.into.realm.rates.BSM.CI)),ylim=c(0,max(c(unlist(hot.into.realm.rates.BSM.CI),unlist(realm.into.hot.rates.BSM.CI))[is.finite(c(unlist(hot.into.realm.rates.BSM.CI),unlist(realm.into.hot.rates.BSM.CI)))])),type='n',xaxt='n',xlab='age',ylab='dispersal rates',main=paste(name,' dispersal rates',sep=''))
    #lines(x=c(1:length(hot.into.realm.rates.BSM.CI)),y=unlist(lapply(hot.into.realm.rates.BSM.CI,function(x)x[3])),col=adjustcolor( "red", alpha.f = 0.9))
    #polygon(x=c((1:length(hot.into.realm.rates.BSM.CI)),(rev(1:length(hot.into.realm.rates.BSM.CI)))),y=c((unlist(lapply(hot.into.realm.rates.BSM.CI,function(x)x[1]))),rev(unlist(lapply(hot.into.realm.rates.BSM.CI,function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
    #lines(x=c(1:length(realm.into.hot.rates.BSM.CI)),y=unlist(lapply(realm.into.hot.rates.BSM.CI,function(x)x[3])),col=adjustcolor( "blue", alpha.f = 0.9))
    #polygon(x=c((1:length(realm.into.hot.rates.BSM.CI)),(rev(1:length(realm.into.hot.rates.BSM.CI)))),y=c((unlist(lapply(realm.into.hot.rates.BSM.CI,function(x)x[1]))),rev(unlist(lapply(realm.into.hot.rates.BSM.CI,function(x)x[2])))),col=adjustcolor( "blue", alpha.f = 0.2))
    #
  }else if(plots=='cladogenesis.BAMM.rates'){
    #plot BAMM in situ rates
    plot(c(1,1),xlim=c(0,length(insitu.hot.rates.BAMM.CI)),ylim=c(0,max(c(unlist(insitu.hot.rates.BAMM.CI),unlist(insitu.realm.rates.BAMM.CI)))),type='n',xaxt='n',xlab='age',ylab='insitu.rates',main=paste(name,' cladogenesis BAMM rates',sep=''))
    lines(x=c(1:length(insitu.hot.rates.BAMM.CI)),y=unlist(lapply(insitu.hot.rates.BAMM.CI,function(x)x[3])),col=adjustcolor( "red", alpha.f = 0.9))
    polygon(x=c((1:length(insitu.hot.rates.BAMM.CI)),(rev(1:length(insitu.hot.rates.BAMM.CI)))),y=c((unlist(lapply(insitu.hot.rates.BAMM.CI,function(x)x[1]))),rev(unlist(lapply(insitu.hot.rates.BAMM.CI,function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
    lines(x=c(1:length(insitu.realm.rates.BAMM.CI)),y=unlist(lapply(insitu.realm.rates.BAMM.CI,function(x)x[3])),col=adjustcolor( "blue", alpha.f = 0.9))
    polygon(x=c((1:length(insitu.realm.rates.BAMM.CI)),(rev(1:length(insitu.realm.rates.BAMM.CI)))),y=c((unlist(lapply(insitu.realm.rates.BAMM.CI,function(x)x[1]))),rev(unlist(lapply(insitu.realm.rates.BAMM.CI,function(x)x[2])))),col=adjustcolor( "blue", alpha.f = 0.2))
    axis(1,at=c(1:length(insitu.realm.rates.BAMM.CI)),labels=c(seq(from=length(insitu.realm.rates.BAMM.CI)*2,to=2,by=-2)))
    legend('topleft',c('insitu.hot','insitu.realm'),col=c('red','blue'),lty=1,cex=.7,bty='n')
  }else if(plots=='cladogenesis.BSM.rates'){
    #plot BSM in situ rates
    plot(c(1,1),xlim=c(0,length(insitu.hot.rates.BSM.CI)),ylim=c(0,max(c(unlist(insitu.hot.rates.BSM.CI),unlist(insitu.realm.rates.BSM.CI))[is.finite(c(unlist(insitu.hot.rates.BSM.CI),unlist(insitu.realm.rates.BSM.CI)))])),type='n',xaxt='n',xlab='age',ylab='insitu.rates',main=paste(name,' cladogenesis BSM rates',sep=''))
    lines(x=c(1:length(insitu.hot.rates.BSM.CI)),y=unlist(lapply(insitu.hot.rates.BSM.CI,function(x)x[3])),col=adjustcolor( "red", alpha.f = 0.9))
    polygon(x=c((1:length(insitu.hot.rates.BSM.CI)),(rev(1:length(insitu.hot.rates.BSM.CI)))),y=c((unlist(lapply(insitu.hot.rates.BSM.CI,function(x)x[1]))),rev(unlist(lapply(insitu.hot.rates.BSM.CI,function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
    lines(x=c(1:length(insitu.realm.rates.BSM.CI)),y=unlist(lapply(insitu.realm.rates.BSM.CI,function(x)x[3])),col=adjustcolor( "blue", alpha.f = 0.9))
    polygon(x=c((1:length(insitu.realm.rates.BSM.CI)),(rev(1:length(insitu.realm.rates.BSM.CI)))),y=c((unlist(lapply(insitu.realm.rates.BSM.CI,function(x)x[1]))),rev(unlist(lapply(insitu.realm.rates.BSM.CI,function(x)x[2])))),col=adjustcolor( "blue", alpha.f = 0.2))
    axis(1,at=c(1:length(insitu.realm.rates.BSM.CI)),labels=c(seq(from=length(insitu.realm.rates.BSM.CI)*2,to=2,by=-2)))
    legend('topleft',c('insitu.hot','insitu.realm'),col=c('red','blue'),lty=1,cex=.7,bty='n')
  }else if(plots=='dispersal.events'){
    plot(c(1,1),xlim=c(0,length(hot.into.realm.time.bins.counts.rev.cumsum.CI)),ylim=c(0,max(c(unlist(hot.into.realm.time.bins.counts.rev.cumsum.CI),unlist(realm.into.hot.time.bins.counts.rev.cumsum.CI)))),type='n',xaxt='n',xlab='age',ylab='cumulative # of events',main=paste(name,' dispersal events',sep=''))
    lines(x=c(1:length(hot.into.realm.time.bins.counts.rev.cumsum.CI)),y=unlist(lapply(hot.into.realm.time.bins.counts.rev.cumsum.CI,function(x)x[3])),col=adjustcolor( "red", alpha.f = 0.9))
    polygon(x=c((1:length(hot.into.realm.time.bins.counts.rev.cumsum.CI)),(rev(1:length(hot.into.realm.time.bins.counts.rev.cumsum.CI)))),y=c((unlist(lapply(hot.into.realm.time.bins.counts.rev.cumsum.CI,function(x)x[1]))),rev(unlist(lapply(hot.into.realm.time.bins.counts.rev.cumsum.CI,function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
    lines(x=c(1:length(realm.into.hot.time.bins.counts.rev.cumsum.CI)),y=unlist(lapply(realm.into.hot.time.bins.counts.rev.cumsum.CI,function(x)x[3])),col=adjustcolor( "blue", alpha.f = 0.9))
    polygon(x=c((1:length(realm.into.hot.time.bins.counts.rev.cumsum.CI)),(rev(1:length(realm.into.hot.time.bins.counts.rev.cumsum.CI)))),y=c((unlist(lapply(realm.into.hot.time.bins.counts.rev.cumsum.CI,function(x)x[1]))),rev(unlist(lapply(realm.into.hot.time.bins.counts.rev.cumsum.CI,function(x)x[2])))),col=adjustcolor( "blue", alpha.f = 0.2))
    axis(1,at=c(1:length(realm.into.hot.time.bins.counts.rev.cumsum.CI)),labels=c(seq(from=length(realm.into.hot.time.bins.counts.rev.cumsum.CI)*2,to=2,by=-2)))
    
    legend('topleft',c('hot.into.realm','realm.into.hot'),col=c('red','blue'),lty=1,cex=.7,bty='n')
    
  }else if(plots=='cladogenesis.events'){
    plot(c(1,1),xlim=c(0,length(insitu.hot.time.bins.counts.rev.cumsum.CI)),ylim=c(0,max(c(unlist(insitu.hot.time.bins.counts.rev.cumsum.CI),unlist(insitu.realm.time.bins.counts.rev.cumsum.CI)))),type='n',xaxt='n',xlab='age',ylab='cumulative # of events',main=paste(name,' cladogenetic events',sep=''))
    lines(x=c(1:length(insitu.hot.time.bins.counts.rev.cumsum.CI)),y=unlist(lapply(insitu.hot.time.bins.counts.rev.cumsum.CI,function(x)x[3])),col=adjustcolor( "red", alpha.f = 0.9))
    polygon(x=c((1:length(insitu.hot.time.bins.counts.rev.cumsum.CI)),(rev(1:length(insitu.hot.time.bins.counts.rev.cumsum.CI)))),y=c((unlist(lapply(insitu.hot.time.bins.counts.rev.cumsum.CI,function(x)x[1]))),rev(unlist(lapply(insitu.hot.time.bins.counts.rev.cumsum.CI,function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
    lines(x=c(1:length(insitu.realm.time.bins.counts.rev.cumsum.CI)),y=unlist(lapply(insitu.realm.time.bins.counts.rev.cumsum.CI,function(x)x[3])),col=adjustcolor( "blue", alpha.f = 0.9))
    polygon(x=c((1:length(insitu.realm.time.bins.counts.rev.cumsum.CI)),(rev(1:length(insitu.realm.time.bins.counts.rev.cumsum.CI)))),y=c((unlist(lapply(insitu.realm.time.bins.counts.rev.cumsum.CI,function(x)x[1]))),rev(unlist(lapply(insitu.realm.time.bins.counts.rev.cumsum.CI,function(x)x[2])))),col=adjustcolor( "blue", alpha.f = 0.2))
    axis(1,at=c(1:length(insitu.realm.time.bins.counts.rev.cumsum.CI)),labels=c(seq(from=length(insitu.realm.time.bins.counts.rev.cumsum.CI)*2,to=2,by=-2)))
    legend('topleft',c('insitu.hot','insitu.realm'),col=c('red','blue'),lty=1,cex=.7,bty='n')
  }else if (plots=='dispersal.all.rates'){
    plot(c(1,1),xlim=c(0,length(dispersal.into.realm.rates.BSM.CI)),ylim=c(0,max(c(unlist(dispersal.into.realm.rates.BSM.CI),unlist(dispersal.into.hot.rates.BSM.CI))[is.finite(c(unlist(dispersal.into.realm.rates.BSM.CI),unlist(dispersal.into.hot.rates.BSM.CI)))])),type='n',xaxt='n',xlab='age',ylab='dispersal rates',main=paste(name,' dispersal rates',sep=''))
    lines(x=c(1:length(dispersal.into.realm.rates.BSM.CI)),y=unlist(lapply(dispersal.into.realm.rates.BSM.CI,function(x)x[3])),col=adjustcolor( "red", alpha.f = 0.9))
    polygon(x=c((1:length(dispersal.into.realm.rates.BSM.CI)),(rev(1:length(dispersal.into.realm.rates.BSM.CI)))),y=c((unlist(lapply(dispersal.into.realm.rates.BSM.CI,function(x)x[1]))),rev(unlist(lapply(dispersal.into.realm.rates.BSM.CI,function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
    lines(x=c(1:length(dispersal.into.hot.rates.BSM.CI)),y=unlist(lapply(dispersal.into.hot.rates.BSM.CI,function(x)x[3])),col=adjustcolor( "blue", alpha.f = 0.9))
    polygon(x=c((1:length(dispersal.into.hot.rates.BSM.CI)),(rev(1:length(dispersal.into.hot.rates.BSM.CI)))),y=c((unlist(lapply(dispersal.into.hot.rates.BSM.CI,function(x)x[1]))),rev(unlist(lapply(dispersal.into.hot.rates.BSM.CI,function(x)x[2])))),col=adjustcolor( "blue", alpha.f = 0.2))
    axis(1,at=c(1:length(dispersal.into.hot.rates.BSM.CI)),labels=c(seq(from=length(dispersal.into.hot.rates.BSM.CI)*2,to=2,by=-2)))
    
    legend('topleft',c('dispersal.into.hot','dispersal.into.realm'),col=c('blue','red'),lty=1,cex=.7,bty='n')
  }else if(plots=='dispersal.all.events'){
    plot(c(1,1),xlim=c(0,length(dispersal.into.realm.time.bins.counts.rev.cumsum.CI)),ylim=c(0,max(c(unlist(dispersal.into.realm.time.bins.counts.rev.cumsum.CI),unlist(dispersal.into.hot.time.bins.counts.rev.cumsum.CI))[is.finite(c(unlist(dispersal.into.realm.time.bins.counts.rev.cumsum.CI),unlist(dispersal.into.hot.time.bins.counts.rev.cumsum.CI)))])),type='n',xaxt='n',xlab='age',ylab='cumulative # of events',main=paste(name,' dispersal events',sep=''))
    lines(x=c(1:length(dispersal.into.realm.time.bins.counts.rev.cumsum.CI)),y=unlist(lapply(dispersal.into.realm.time.bins.counts.rev.cumsum.CI,function(x)x[3])),col=adjustcolor( "red", alpha.f = 0.9))
    polygon(x=c((1:length(dispersal.into.realm.time.bins.counts.rev.cumsum.CI)),(rev(1:length(dispersal.into.realm.time.bins.counts.rev.cumsum.CI)))),y=c((unlist(lapply(dispersal.into.realm.time.bins.counts.rev.cumsum.CI,function(x)x[1]))),rev(unlist(lapply(dispersal.into.realm.time.bins.counts.rev.cumsum.CI,function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
    lines(x=c(1:length(dispersal.into.hot.time.bins.counts.rev.cumsum.CI)),y=unlist(lapply(dispersal.into.hot.time.bins.counts.rev.cumsum.CI,function(x)x[3])),col=adjustcolor( "blue", alpha.f = 0.9))
    polygon(x=c((1:length(dispersal.into.hot.time.bins.counts.rev.cumsum.CI)),(rev(1:length(dispersal.into.hot.time.bins.counts.rev.cumsum.CI)))),y=c((unlist(lapply(dispersal.into.hot.time.bins.counts.rev.cumsum.CI,function(x)x[1]))),rev(unlist(lapply(dispersal.into.hot.time.bins.counts.rev.cumsum.CI,function(x)x[2])))),col=adjustcolor( "blue", alpha.f = 0.2))
    axis(1,at=c(1:length(dispersal.into.hot.time.bins.counts.rev.cumsum.CI)),labels=c(seq(from=length(dispersal.into.hot.time.bins.counts.rev.cumsum.CI)*2,to=2,by=-2)))
    legend('topleft',c('dispersal.into.realm','dispersal.into.hot'),col=c('red','blue'),lty=1,cex=.7,bty='n')
  }else if(plots=='dispersal.into.from.events'){
    plot(c(1,1),xlim=c(0,length(dispersal.into.realm.time.bins.counts.rev.cumsum.CI)),ylim=c(0,max(c(unlist(dispersal.into.realm.time.bins.counts.rev.cumsum.CI),unlist(dispersal.into.hot.time.bins.counts.rev.cumsum.CI),unlist(dispersal.from.hot.time.bins.counts.rev.cumsum.CI),unlist(dispersal.from.realm.time.bins.counts.rev.cumsum.CI)))),type='n',xaxt='n',xlab='age',ylab='cumulative # of events',main=paste(name,' dispersal events',sep=''))
    lines(x=c(1:length(dispersal.into.realm.time.bins.counts.rev.cumsum.CI)),y=unlist(lapply(dispersal.into.realm.time.bins.counts.rev.cumsum.CI,function(x)x[3])),col=adjustcolor( "green", alpha.f = 0.9))
    polygon(x=c((1:length(dispersal.into.realm.time.bins.counts.rev.cumsum.CI)),(rev(1:length(dispersal.into.realm.time.bins.counts.rev.cumsum.CI)))),y=c((unlist(lapply(dispersal.into.realm.time.bins.counts.rev.cumsum.CI,function(x)x[1]))),rev(unlist(lapply(dispersal.into.realm.time.bins.counts.rev.cumsum.CI,function(x)x[2])))),col=adjustcolor( "green", alpha.f = 0.2))
    lines(x=c(1:length(dispersal.into.hot.time.bins.counts.rev.cumsum.CI)),y=unlist(lapply(dispersal.into.hot.time.bins.counts.rev.cumsum.CI,function(x)x[3])),col=adjustcolor( "orange", alpha.f = 0.9))
    polygon(x=c((1:length(dispersal.into.hot.time.bins.counts.rev.cumsum.CI)),(rev(1:length(dispersal.into.hot.time.bins.counts.rev.cumsum.CI)))),y=c((unlist(lapply(dispersal.into.hot.time.bins.counts.rev.cumsum.CI,function(x)x[1]))),rev(unlist(lapply(dispersal.into.hot.time.bins.counts.rev.cumsum.CI,function(x)x[2])))),col=adjustcolor( "orange", alpha.f = 0.2))
    lines(x=c(1:length(dispersal.from.realm.time.bins.counts.rev.cumsum.CI)),y=unlist(lapply(dispersal.from.realm.time.bins.counts.rev.cumsum.CI,function(x)x[3])),col=adjustcolor( "blue", alpha.f = 0.9))
    polygon(x=c((1:length(dispersal.from.realm.time.bins.counts.rev.cumsum.CI)),(rev(1:length(dispersal.from.realm.time.bins.counts.rev.cumsum.CI)))),y=c((unlist(lapply(dispersal.from.realm.time.bins.counts.rev.cumsum.CI,function(x)x[1]))),rev(unlist(lapply(dispersal.from.realm.time.bins.counts.rev.cumsum.CI,function(x)x[2])))),col=adjustcolor( "blue", alpha.f = 0.2))
    lines(x=c(1:length(dispersal.from.hot.time.bins.counts.rev.cumsum.CI)),y=unlist(lapply(dispersal.from.hot.time.bins.counts.rev.cumsum.CI,function(x)x[3])),col=adjustcolor( "red", alpha.f = 0.9))
    polygon(x=c((1:length(dispersal.from.hot.time.bins.counts.rev.cumsum.CI)),(rev(1:length(dispersal.from.hot.time.bins.counts.rev.cumsum.CI)))),y=c((unlist(lapply(dispersal.from.hot.time.bins.counts.rev.cumsum.CI,function(x)x[1]))),rev(unlist(lapply(dispersal.from.hot.time.bins.counts.rev.cumsum.CI,function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
    axis(1,at=c(1:length(dispersal.from.hot.time.bins.counts.rev.cumsum.CI)),labels=c(seq(from=length(dispersal.from.hot.time.bins.counts.rev.cumsum.CI)*2,to=2,by=-2)))
    legend('topleft',c('dispersal.from.realm','dispersal.from.hot','dispersal.into.realm','dispersal.into.hot'),col=c('blue','red','green','orange'),lty=1,cex=.7,bty='n')
  }
  
  #plot BSM dispersal rates hot into realm
  
  ##plot BAMM in situ rates
  #par(new=T)
  #plot(c(1,1),xlim=c(0,length(insitu.hot.rates.BAMM.CI)),ylim=c(0,max(c(unlist(insitu.hot.rates.BAMM.CI),unlist(insitu.realm.rates.BAMM.CI)))),type='n',xaxt='n',xlab='age',yaxt='n')
  ##polygon(x=c(c(age:length(insitu.hot.rates.BAMM.CI)),c(length(insitu.hot.rates.BAMM.CI):age)),y=c(unlist(lapply(insitu.hot.rates.BAMM.CI[c(age:length(insitu.hot.rates.BAMM.CI))],function(x)x[1])),rev(unlist(lapply(insitu.hot.rates.BAMM.CI[c(age:length(insitu.hot.rates.BAMM.CI))],function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
  #lines(x=c(1:length(insitu.hot.rates.BAMM.CI)),y=unlist(lapply(insitu.hot.rates.BAMM.CI,function(x)x[3])),col=adjustcolor( "red", alpha.f = 0.9))
  #polygon(x=c((1:length(insitu.hot.rates.BAMM.CI)),(rev(1:length(insitu.hot.rates.BAMM.CI)))),y=c((unlist(lapply(insitu.hot.rates.BAMM.CI,function(x)x[1]))),rev(unlist(lapply(insitu.hot.rates.BAMM.CI,function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
  #lines(x=c(1:length(insitu.realm.rates.BAMM.CI)),y=unlist(lapply(insitu.realm.rates.BAMM.CI,function(x)x[3])),col=adjustcolor( "blue", alpha.f = 0.9))
  #polygon(x=c((1:length(insitu.realm.rates.BAMM.CI)),(rev(1:length(insitu.realm.rates.BAMM.CI)))),y=c((unlist(lapply(insitu.realm.rates.BAMM.CI,function(x)x[1]))),rev(unlist(lapply(insitu.realm.rates.BAMM.CI,function(x)x[2])))),col=adjustcolor( "blue", alpha.f = 0.2))
  ##axis(1,at=c(0:length(BSMoutput.real.rates[['insitu.realm.rates.BSM.CI']])),labels=c((length(BSMoutput.real.rates[['insitu.realm.rates.BSM.CI']])+2):2))
  #axis(4)
  #legend('topright',c('insitu.hot','insitu.realm'),col=c('red','blue'),lty=1,cex=.7,bty='n')
}



plot_rates_events_real_2mybin_difference<-function(BSMoutput.real.rates,BSMoutput.real.events,age,name,plots){
  #correct age
  age<-(length(BSMoutput.real.rates[['insitu.realm.rates.BSM.CI']])*2-age)/2
  #trim lists to age period
  BSMoutput.real.rates<-lapply(BSMoutput.real.rates,function(x)x[age:length(x)])
  BSMoutput.real.events<-lapply(BSMoutput.real.events,function(x)x[age:length(x)])
  #delete Inf in rates and events
  BSMoutput.real.rates<-lapply(BSMoutput.real.rates,function(x) lapply(x,function(x) {x[!is.finite(x)]<-NA;return(x)}))
  BSMoutput.real.events<-lapply(BSMoutput.real.events,function(x) lapply(x,function(x) {x[!is.finite(x)]<-NA;return(x)}))
  #remove the last bin - two million year time bin
  BSMoutput.real.rates<-lapply(BSMoutput.real.rates,function(x) x<-x[-c(length(x))])
  BSMoutput.real.events<-lapply(BSMoutput.real.events,function(x) x<-x[-c((length(x)-1),length(x))])
  #assign variables for easier handling
  hot.into.realm.rates.BSM<-BSMoutput.real.rates[['hot.into.realm.rates.BSM.CI']]
  realm.into.hot.rates.BSM<-BSMoutput.real.rates[['realm.into.hot.rates.BSM.CI']]
  insitu.realm.rates.BAMM<-BSMoutput.real.rates[['insitu.realm.rates.BAMM.CI']]
  insitu.hot.rates.BAMM<-BSMoutput.real.rates[['insitu.hot.rates.BAMM.CI']]
  insitu.realm.rates.BSM<-BSMoutput.real.rates[['insitu.realm.rates.BSM.CI']]
  insitu.hot.rates.BSM<-BSMoutput.real.rates[['insitu.hot.rates.BSM.CI']]
  dispersal.into.realm.rates.BSM<-BSMoutput.real.rates[['dispersal.into.realm.rates.BSM.CI']]
  dispersal.into.hot.rates.BSM<-BSMoutput.real.rates[['dispersal.into.hot.rates.BSM.CI']]
  dispersal.from.realm.rates.BSM<-BSMoutput.real.rates[['dispersal.from.realm.rates.BSM.CI']]
  dispersal.from.hot.rates.BSM<-BSMoutput.real.rates[['dispersal.from.hot.rates.BSM.CI']]
  
  hot.into.realm.time.bins.counts.rev.cumsum<-BSMoutput.real.events[['hot.into.realm.time.bins.counts.rev.cumsum.CI']]
  realm.into.hot.time.bins.counts.rev.cumsum<-BSMoutput.real.events[['realm.into.hot.time.bins.counts.rev.cumsum.CI']]
  insitu.realm.time.bins.counts.rev.cumsum<-BSMoutput.real.events[['insitu.realm.time.bins.counts.rev.cumsum.CI']]
  insitu.hot.time.bins.counts.rev.cumsum<-BSMoutput.real.events[['insitu.hot.time.bins.counts.rev.cumsum.CI']]
  dispersal.into.realm.time.bins.counts.rev.cumsum<-BSMoutput.real.events[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']]
  dispersal.into.hot.time.bins.counts.rev.cumsum<-BSMoutput.real.events[['dispersal.into.hot.time.bins.counts.rev.cumsum.CI']]
  dispersal.from.realm.time.bins.counts.rev.cumsum<-BSMoutput.real.events[['dispersal.from.realm.time.bins.counts.rev.cumsum.CI']]
  dispersal.from.hot.time.bins.counts.rev.cumsum<-BSMoutput.real.events[['dispersal.from.hot.time.bins.counts.rev.cumsum.CI']]
  if(plots=='within.realm.dispersal.rates'){
    #difference.into.hotminusrealm<-lapply(c(1:length(dispersal.into.realm.rates.BSM)),function(x)sample(dispersal.into.hot.rates.BSM[[x]],size = 100,replace=TRUE)-sample(dispersal.into.realm.rates.BSM[[x]],size = 100,replace=TRUE))
    difference.hotintorealm.realmintohot<-lapply(c(1:length(hot.into.realm.rates.BSM)),function(x)sample(hot.into.realm.rates.BSM[[x]],size = 100,replace=TRUE)-sample(realm.into.hot.rates.BSM[[x]],size = 100,replace=TRUE))
    median<-unlist(lapply(difference.hotintorealm.realmintohot,function(x)median(x,na.rm=TRUE)))
    lower<-unlist(lapply(difference.hotintorealm.realmintohot,function(x)quantile(x,c(0.1),na.rm=TRUE)))
    upper<-unlist(lapply(difference.hotintorealm.realmintohot,function(x)quantile(x,c(0.9),na.rm=TRUE)))
    #lower<-unlist(lapply(difference.hotintorealm.realmintohot,function(x)quantile(x,c(0.025),na.rm=TRUE)))
    #upper<-unlist(lapply(difference.hotintorealm.realmintohot,function(x)quantile(x,c(0.975),na.rm=TRUE)))
    plot(x=c(1:length(median)),y=median,ylim=c(-max(abs(c(median,upper,lower))),max(abs(c(median,upper,lower)))),type='n',xaxt='n',xlab='age',yaxt='n',ylab='diff. dispersal hot.to.realm - realm.to.hot',main=paste(name,' diff in dispersal rates',sep=''))
    polygon(x=c(c(1:length(median),length(median):1)),y=c(lower,rev(upper)),col=adjustcolor('red',alpha.f = 0.25))
    lines(x=c(1:length(difference.hotintorealm.realmintohot)),y=median,col='red') 
    abline(h=0,lty=2)
    axis(1,at=c(1:length(difference.hotintorealm.realmintohot)),labels=c(seq(from=length(difference.hotintorealm.realmintohot)*2,to=2,by=-2)))
    axis(2,las=2)
    #plot(c(1,1),xlim=c(0,length(hot.into.realm.rates.BSM.CI)),ylim=c(0,max(c(unlist(hot.into.realm.rates.BSM.CI),unlist(realm.into.hot.rates.BSM.CI))[is.finite(c(unlist(hot.into.realm.rates.BSM.CI),unlist(realm.into.hot.rates.BSM.CI)))])),type='n',xaxt='n',xlab='age',ylab='dispersal rates',main=paste(name,' dispersal rates',sep=''))
    #lines(x=c(1:length(hot.into.realm.rates.BSM.CI)),y=unlist(lapply(hot.into.realm.rates.BSM.CI,function(x)x[3])),col=adjustcolor( "red", alpha.f = 0.9))
    #polygon(x=c((1:length(hot.into.realm.rates.BSM.CI)),(rev(1:length(hot.into.realm.rates.BSM.CI)))),y=c((unlist(lapply(hot.into.realm.rates.BSM.CI,function(x)x[1]))),rev(unlist(lapply(hot.into.realm.rates.BSM.CI,function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
    #lines(x=c(1:length(realm.into.hot.rates.BSM.CI)),y=unlist(lapply(realm.into.hot.rates.BSM.CI,function(x)x[3])),col=adjustcolor( "blue", alpha.f = 0.9))
    #polygon(x=c((1:length(realm.into.hot.rates.BSM.CI)),(rev(1:length(realm.into.hot.rates.BSM.CI)))),y=c((unlist(lapply(realm.into.hot.rates.BSM.CI,function(x)x[1]))),rev(unlist(lapply(realm.into.hot.rates.BSM.CI,function(x)x[2])))),col=adjustcolor( "blue", alpha.f = 0.2))
    #
    legend('topleft',c('diff.rates.hot-realm'),col=c('red'),lty=1,cex=.7,bty='n')
  }else if(plots=='global.dispersal.rates'){
    #difference.into.hotminusrealm<-lapply(c(1:length(dispersal.into.realm.rates.BSM)),function(x)sample(dispersal.into.hot.rates.BSM[[x]],size = 100,replace=TRUE)-sample(dispersal.into.realm.rates.BSM[[x]],size = 100,replace=TRUE))
    difference.hotrealm.from<-lapply(c(1:length(dispersal.from.realm.rates.BSM)),function(x)sample(dispersal.from.hot.rates.BSM[[x]],size = 100,replace=TRUE)-sample(dispersal.from.realm.rates.BSM[[x]],size = 100,replace=TRUE))
    difference.hotrealm.into<-lapply(c(1:length(dispersal.into.realm.rates.BSM)),function(x)sample(dispersal.into.hot.rates.BSM[[x]],size = 100,replace=TRUE)-sample(dispersal.into.realm.rates.BSM[[x]],size = 100,replace=TRUE))
    median.from<-unlist(lapply(difference.hotrealm.from,function(x)median(x,na.rm=TRUE)))
    lower.from<-unlist(lapply(difference.hotrealm.from,function(x)quantile(x,c(0.1),na.rm=TRUE)))
    upper.from<-unlist(lapply(difference.hotrealm.from,function(x)quantile(x,c(0.9),na.rm=TRUE)))
    median.into<-unlist(lapply(difference.hotrealm.into,function(x)median(x,na.rm=TRUE)))
    lower.into<-unlist(lapply(difference.hotrealm.into,function(x)quantile(x,c(0.1),na.rm=TRUE)))
    upper.into<-unlist(lapply(difference.hotrealm.into,function(x)quantile(x,c(0.9),na.rm=TRUE)))
    
    plot(x=c(1:length(median.from)),y=median.from,ylim=c(-max(abs(c(median.from,lower.from,upper.from,median.into,lower.into,upper.into))),max(abs(c(median.from,lower.from,upper.from,median.into,lower.into,upper.into)))),type='n',xaxt='n',xlab='age',yaxt='n',ylab='diff. dispersal hot - realm',main=paste(name,' diff in global dispersal rates',sep=''))
    polygon(x=c(c(1:length(median.from),length(median.from):1)),y=c(lower.from,rev(upper.from)),col=adjustcolor('red',alpha.f = 0.25))
    lines(x=c(1:length(median.from)),y=median.from,col='red') 
    abline(h=0,lty=2)
    polygon(x=c(c(1:length(median.into),length(median.into):1)),y=c(lower.into,rev(upper.into)),col=adjustcolor('blue',alpha.f = 0.25))
    lines(x=c(1:length(median.into)),y=median.into,col='blue') 
    axis(1,at=c(1:length(difference.hotintorealm.realmintohot)),labels=c(seq(from=length(difference.hotintorealm.realmintohot)*2,to=2,by=-2)))
    axis(2,las=2)
    legend('topleft',c('dispersal.from.diff','dispersal.into.diff'),col=c('red','blue'),lty=1,cex=.7,bty='n')
    #plot(c(1,1),xlim=c(0,length(hot.into.realm.rates.BSM.CI)),ylim=c(0,max(c(unlist(hot.into.realm.rates.BSM.CI),unlist(realm.into.hot.rates.BSM.CI))[is.finite(c(unlist(hot.into.realm.rates.BSM.CI),unlist(realm.into.hot.rates.BSM.CI)))])),type='n',xaxt='n',xlab='age',ylab='dispersal rates',main=paste(name,' dispersal rates',sep=''))
    #lines(x=c(1:length(hot.into.realm.rates.BSM.CI)),y=unlist(lapply(hot.into.realm.rates.BSM.CI,function(x)x[3])),col=adjustcolor( "red", alpha.f = 0.9))
    #polygon(x=c((1:length(hot.into.realm.rates.BSM.CI)),(rev(1:length(hot.into.realm.rates.BSM.CI)))),y=c((unlist(lapply(hot.into.realm.rates.BSM.CI,function(x)x[1]))),rev(unlist(lapply(hot.into.realm.rates.BSM.CI,function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
    #lines(x=c(1:length(realm.into.hot.rates.BSM.CI)),y=unlist(lapply(realm.into.hot.rates.BSM.CI,function(x)x[3])),col=adjustcolor( "blue", alpha.f = 0.9))
    #polygon(x=c((1:length(realm.into.hot.rates.BSM.CI)),(rev(1:length(realm.into.hot.rates.BSM.CI)))),y=c((unlist(lapply(realm.into.hot.rates.BSM.CI,function(x)x[1]))),rev(unlist(lapply(realm.into.hot.rates.BSM.CI,function(x)x[2])))),col=adjustcolor( "blue", alpha.f = 0.2))
    #
  }else if(plots=='cladogenesis.BAMM.rates'){
    #plot BAMM in situ rates
    plot(c(1,1),xlim=c(0,length(insitu.hot.rates.BAMM.CI)),ylim=c(0,max(c(unlist(insitu.hot.rates.BAMM.CI),unlist(insitu.realm.rates.BAMM.CI)))),type='n',xaxt='n',xlab='age',ylab='insitu.rates',main=paste(name,' cladogenesis BAMM rates',sep=''))
    lines(x=c(1:length(insitu.hot.rates.BAMM.CI)),y=unlist(lapply(insitu.hot.rates.BAMM.CI,function(x)x[3])),col=adjustcolor( "red", alpha.f = 0.9))
    polygon(x=c((1:length(insitu.hot.rates.BAMM.CI)),(rev(1:length(insitu.hot.rates.BAMM.CI)))),y=c((unlist(lapply(insitu.hot.rates.BAMM.CI,function(x)x[1]))),rev(unlist(lapply(insitu.hot.rates.BAMM.CI,function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
    lines(x=c(1:length(insitu.realm.rates.BAMM.CI)),y=unlist(lapply(insitu.realm.rates.BAMM.CI,function(x)x[3])),col=adjustcolor( "blue", alpha.f = 0.9))
    polygon(x=c((1:length(insitu.realm.rates.BAMM.CI)),(rev(1:length(insitu.realm.rates.BAMM.CI)))),y=c((unlist(lapply(insitu.realm.rates.BAMM.CI,function(x)x[1]))),rev(unlist(lapply(insitu.realm.rates.BAMM.CI,function(x)x[2])))),col=adjustcolor( "blue", alpha.f = 0.2))
    axis(1,at=c(1:length(insitu.realm.rates.BAMM.CI)),labels=c(seq(from=length(insitu.realm.rates.BAMM.CI)*2,to=2,by=-2)))
    legend('topleft',c('insitu.hot','insitu.realm'),col=c('red','blue'),lty=1,cex=.7,bty='n')
  }else if(plots=='cladogenesis.BSM.rates'){
    #plot BSM in situ rates
    plot(c(1,1),xlim=c(0,length(insitu.hot.rates.BSM.CI)),ylim=c(0,max(c(unlist(insitu.hot.rates.BSM.CI),unlist(insitu.realm.rates.BSM.CI))[is.finite(c(unlist(insitu.hot.rates.BSM.CI),unlist(insitu.realm.rates.BSM.CI)))])),type='n',xaxt='n',xlab='age',ylab='insitu.rates',main=paste(name,' cladogenesis BSM rates',sep=''))
    lines(x=c(1:length(insitu.hot.rates.BSM.CI)),y=unlist(lapply(insitu.hot.rates.BSM.CI,function(x)x[3])),col=adjustcolor( "red", alpha.f = 0.9))
    polygon(x=c((1:length(insitu.hot.rates.BSM.CI)),(rev(1:length(insitu.hot.rates.BSM.CI)))),y=c((unlist(lapply(insitu.hot.rates.BSM.CI,function(x)x[1]))),rev(unlist(lapply(insitu.hot.rates.BSM.CI,function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
    lines(x=c(1:length(insitu.realm.rates.BSM.CI)),y=unlist(lapply(insitu.realm.rates.BSM.CI,function(x)x[3])),col=adjustcolor( "blue", alpha.f = 0.9))
    polygon(x=c((1:length(insitu.realm.rates.BSM.CI)),(rev(1:length(insitu.realm.rates.BSM.CI)))),y=c((unlist(lapply(insitu.realm.rates.BSM.CI,function(x)x[1]))),rev(unlist(lapply(insitu.realm.rates.BSM.CI,function(x)x[2])))),col=adjustcolor( "blue", alpha.f = 0.2))
    axis(1,at=c(1:length(insitu.realm.rates.BSM.CI)),labels=c(seq(from=length(insitu.realm.rates.BSM.CI)*2,to=2,by=-2)))
    legend('topleft',c('insitu.hot','insitu.realm'),col=c('red','blue'),lty=1,cex=.7,bty='n')
  }else if(plots=='dispersal.events'){
    plot(c(1,1),xlim=c(0,length(hot.into.realm.time.bins.counts.rev.cumsum.CI)),ylim=c(0,max(c(unlist(hot.into.realm.time.bins.counts.rev.cumsum.CI),unlist(realm.into.hot.time.bins.counts.rev.cumsum.CI)))),type='n',xaxt='n',xlab='age',ylab='cumulative # of events',main=paste(name,' dispersal events',sep=''))
    lines(x=c(1:length(hot.into.realm.time.bins.counts.rev.cumsum.CI)),y=unlist(lapply(hot.into.realm.time.bins.counts.rev.cumsum.CI,function(x)x[3])),col=adjustcolor( "red", alpha.f = 0.9))
    polygon(x=c((1:length(hot.into.realm.time.bins.counts.rev.cumsum.CI)),(rev(1:length(hot.into.realm.time.bins.counts.rev.cumsum.CI)))),y=c((unlist(lapply(hot.into.realm.time.bins.counts.rev.cumsum.CI,function(x)x[1]))),rev(unlist(lapply(hot.into.realm.time.bins.counts.rev.cumsum.CI,function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
    lines(x=c(1:length(realm.into.hot.time.bins.counts.rev.cumsum.CI)),y=unlist(lapply(realm.into.hot.time.bins.counts.rev.cumsum.CI,function(x)x[3])),col=adjustcolor( "blue", alpha.f = 0.9))
    polygon(x=c((1:length(realm.into.hot.time.bins.counts.rev.cumsum.CI)),(rev(1:length(realm.into.hot.time.bins.counts.rev.cumsum.CI)))),y=c((unlist(lapply(realm.into.hot.time.bins.counts.rev.cumsum.CI,function(x)x[1]))),rev(unlist(lapply(realm.into.hot.time.bins.counts.rev.cumsum.CI,function(x)x[2])))),col=adjustcolor( "blue", alpha.f = 0.2))
    axis(1,at=c(1:length(realm.into.hot.time.bins.counts.rev.cumsum.CI)),labels=c(seq(from=length(realm.into.hot.time.bins.counts.rev.cumsum.CI)*2,to=2,by=-2)))
    
    legend('topleft',c('hot.into.realm','realm.into.hot'),col=c('red','blue'),lty=1,cex=.7,bty='n')
    
  }else if(plots=='cladogenesis.events'){
    plot(c(1,1),xlim=c(0,length(insitu.hot.time.bins.counts.rev.cumsum.CI)),ylim=c(0,max(c(unlist(insitu.hot.time.bins.counts.rev.cumsum.CI),unlist(insitu.realm.time.bins.counts.rev.cumsum.CI)))),type='n',xaxt='n',xlab='age',ylab='cumulative # of events',main=paste(name,' cladogenetic events',sep=''))
    lines(x=c(1:length(insitu.hot.time.bins.counts.rev.cumsum.CI)),y=unlist(lapply(insitu.hot.time.bins.counts.rev.cumsum.CI,function(x)x[3])),col=adjustcolor( "red", alpha.f = 0.9))
    polygon(x=c((1:length(insitu.hot.time.bins.counts.rev.cumsum.CI)),(rev(1:length(insitu.hot.time.bins.counts.rev.cumsum.CI)))),y=c((unlist(lapply(insitu.hot.time.bins.counts.rev.cumsum.CI,function(x)x[1]))),rev(unlist(lapply(insitu.hot.time.bins.counts.rev.cumsum.CI,function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
    lines(x=c(1:length(insitu.realm.time.bins.counts.rev.cumsum.CI)),y=unlist(lapply(insitu.realm.time.bins.counts.rev.cumsum.CI,function(x)x[3])),col=adjustcolor( "blue", alpha.f = 0.9))
    polygon(x=c((1:length(insitu.realm.time.bins.counts.rev.cumsum.CI)),(rev(1:length(insitu.realm.time.bins.counts.rev.cumsum.CI)))),y=c((unlist(lapply(insitu.realm.time.bins.counts.rev.cumsum.CI,function(x)x[1]))),rev(unlist(lapply(insitu.realm.time.bins.counts.rev.cumsum.CI,function(x)x[2])))),col=adjustcolor( "blue", alpha.f = 0.2))
    axis(1,at=c(1:length(insitu.realm.time.bins.counts.rev.cumsum.CI)),labels=c(seq(from=length(insitu.realm.time.bins.counts.rev.cumsum.CI)*2,to=2,by=-2)))
    legend('topleft',c('insitu.hot','insitu.realm'),col=c('red','blue'),lty=1,cex=.7,bty='n')
  }else if (plots=='dispersal.all.rates'){
    plot(c(1,1),xlim=c(0,length(dispersal.into.realm.rates.BSM.CI)),ylim=c(0,max(c(unlist(dispersal.into.realm.rates.BSM.CI),unlist(dispersal.into.hot.rates.BSM.CI))[is.finite(c(unlist(dispersal.into.realm.rates.BSM.CI),unlist(dispersal.into.hot.rates.BSM.CI)))])),type='n',xaxt='n',xlab='age',ylab='dispersal rates',main=paste(name,' dispersal rates',sep=''))
    lines(x=c(1:length(dispersal.into.realm.rates.BSM.CI)),y=unlist(lapply(dispersal.into.realm.rates.BSM.CI,function(x)x[3])),col=adjustcolor( "red", alpha.f = 0.9))
    polygon(x=c((1:length(dispersal.into.realm.rates.BSM.CI)),(rev(1:length(dispersal.into.realm.rates.BSM.CI)))),y=c((unlist(lapply(dispersal.into.realm.rates.BSM.CI,function(x)x[1]))),rev(unlist(lapply(dispersal.into.realm.rates.BSM.CI,function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
    lines(x=c(1:length(dispersal.into.hot.rates.BSM.CI)),y=unlist(lapply(dispersal.into.hot.rates.BSM.CI,function(x)x[3])),col=adjustcolor( "blue", alpha.f = 0.9))
    polygon(x=c((1:length(dispersal.into.hot.rates.BSM.CI)),(rev(1:length(dispersal.into.hot.rates.BSM.CI)))),y=c((unlist(lapply(dispersal.into.hot.rates.BSM.CI,function(x)x[1]))),rev(unlist(lapply(dispersal.into.hot.rates.BSM.CI,function(x)x[2])))),col=adjustcolor( "blue", alpha.f = 0.2))
    axis(1,at=c(1:length(dispersal.into.hot.rates.BSM.CI)),labels=c(seq(from=length(dispersal.into.hot.rates.BSM.CI)*2,to=2,by=-2)))
    
    legend('topleft',c('dispersal.into.hot','dispersal.into.realm'),col=c('blue','red'),lty=1,cex=.7,bty='n')
  }else if(plots=='dispersal.all.events'){
    plot(c(1,1),xlim=c(0,length(dispersal.into.realm.time.bins.counts.rev.cumsum.CI)),ylim=c(0,max(c(unlist(dispersal.into.realm.time.bins.counts.rev.cumsum.CI),unlist(dispersal.into.hot.time.bins.counts.rev.cumsum.CI))[is.finite(c(unlist(dispersal.into.realm.time.bins.counts.rev.cumsum.CI),unlist(dispersal.into.hot.time.bins.counts.rev.cumsum.CI)))])),type='n',xaxt='n',xlab='age',ylab='cumulative # of events',main=paste(name,' dispersal events',sep=''))
    lines(x=c(1:length(dispersal.into.realm.time.bins.counts.rev.cumsum.CI)),y=unlist(lapply(dispersal.into.realm.time.bins.counts.rev.cumsum.CI,function(x)x[3])),col=adjustcolor( "red", alpha.f = 0.9))
    polygon(x=c((1:length(dispersal.into.realm.time.bins.counts.rev.cumsum.CI)),(rev(1:length(dispersal.into.realm.time.bins.counts.rev.cumsum.CI)))),y=c((unlist(lapply(dispersal.into.realm.time.bins.counts.rev.cumsum.CI,function(x)x[1]))),rev(unlist(lapply(dispersal.into.realm.time.bins.counts.rev.cumsum.CI,function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
    lines(x=c(1:length(dispersal.into.hot.time.bins.counts.rev.cumsum.CI)),y=unlist(lapply(dispersal.into.hot.time.bins.counts.rev.cumsum.CI,function(x)x[3])),col=adjustcolor( "blue", alpha.f = 0.9))
    polygon(x=c((1:length(dispersal.into.hot.time.bins.counts.rev.cumsum.CI)),(rev(1:length(dispersal.into.hot.time.bins.counts.rev.cumsum.CI)))),y=c((unlist(lapply(dispersal.into.hot.time.bins.counts.rev.cumsum.CI,function(x)x[1]))),rev(unlist(lapply(dispersal.into.hot.time.bins.counts.rev.cumsum.CI,function(x)x[2])))),col=adjustcolor( "blue", alpha.f = 0.2))
    axis(1,at=c(1:length(dispersal.into.hot.time.bins.counts.rev.cumsum.CI)),labels=c(seq(from=length(dispersal.into.hot.time.bins.counts.rev.cumsum.CI)*2,to=2,by=-2)))
    legend('topleft',c('dispersal.into.realm','dispersal.into.hot'),col=c('red','blue'),lty=1,cex=.7,bty='n')
  }else if(plots=='dispersal.into.from.events'){
    plot(c(1,1),xlim=c(0,length(dispersal.into.realm.time.bins.counts.rev.cumsum.CI)),ylim=c(0,max(c(unlist(dispersal.into.realm.time.bins.counts.rev.cumsum.CI),unlist(dispersal.into.hot.time.bins.counts.rev.cumsum.CI),unlist(dispersal.from.hot.time.bins.counts.rev.cumsum.CI),unlist(dispersal.from.realm.time.bins.counts.rev.cumsum.CI)))),type='n',xaxt='n',xlab='age',ylab='cumulative # of events',main=paste(name,' dispersal events',sep=''))
    lines(x=c(1:length(dispersal.into.realm.time.bins.counts.rev.cumsum.CI)),y=unlist(lapply(dispersal.into.realm.time.bins.counts.rev.cumsum.CI,function(x)x[3])),col=adjustcolor( "green", alpha.f = 0.9))
    polygon(x=c((1:length(dispersal.into.realm.time.bins.counts.rev.cumsum.CI)),(rev(1:length(dispersal.into.realm.time.bins.counts.rev.cumsum.CI)))),y=c((unlist(lapply(dispersal.into.realm.time.bins.counts.rev.cumsum.CI,function(x)x[1]))),rev(unlist(lapply(dispersal.into.realm.time.bins.counts.rev.cumsum.CI,function(x)x[2])))),col=adjustcolor( "green", alpha.f = 0.2))
    lines(x=c(1:length(dispersal.into.hot.time.bins.counts.rev.cumsum.CI)),y=unlist(lapply(dispersal.into.hot.time.bins.counts.rev.cumsum.CI,function(x)x[3])),col=adjustcolor( "orange", alpha.f = 0.9))
    polygon(x=c((1:length(dispersal.into.hot.time.bins.counts.rev.cumsum.CI)),(rev(1:length(dispersal.into.hot.time.bins.counts.rev.cumsum.CI)))),y=c((unlist(lapply(dispersal.into.hot.time.bins.counts.rev.cumsum.CI,function(x)x[1]))),rev(unlist(lapply(dispersal.into.hot.time.bins.counts.rev.cumsum.CI,function(x)x[2])))),col=adjustcolor( "orange", alpha.f = 0.2))
    lines(x=c(1:length(dispersal.from.realm.time.bins.counts.rev.cumsum.CI)),y=unlist(lapply(dispersal.from.realm.time.bins.counts.rev.cumsum.CI,function(x)x[3])),col=adjustcolor( "blue", alpha.f = 0.9))
    polygon(x=c((1:length(dispersal.from.realm.time.bins.counts.rev.cumsum.CI)),(rev(1:length(dispersal.from.realm.time.bins.counts.rev.cumsum.CI)))),y=c((unlist(lapply(dispersal.from.realm.time.bins.counts.rev.cumsum.CI,function(x)x[1]))),rev(unlist(lapply(dispersal.from.realm.time.bins.counts.rev.cumsum.CI,function(x)x[2])))),col=adjustcolor( "blue", alpha.f = 0.2))
    lines(x=c(1:length(dispersal.from.hot.time.bins.counts.rev.cumsum.CI)),y=unlist(lapply(dispersal.from.hot.time.bins.counts.rev.cumsum.CI,function(x)x[3])),col=adjustcolor( "red", alpha.f = 0.9))
    polygon(x=c((1:length(dispersal.from.hot.time.bins.counts.rev.cumsum.CI)),(rev(1:length(dispersal.from.hot.time.bins.counts.rev.cumsum.CI)))),y=c((unlist(lapply(dispersal.from.hot.time.bins.counts.rev.cumsum.CI,function(x)x[1]))),rev(unlist(lapply(dispersal.from.hot.time.bins.counts.rev.cumsum.CI,function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
    axis(1,at=c(1:length(dispersal.from.hot.time.bins.counts.rev.cumsum.CI)),labels=c(seq(from=length(dispersal.from.hot.time.bins.counts.rev.cumsum.CI)*2,to=2,by=-2)))
    legend('topleft',c('dispersal.from.realm','dispersal.from.hot','dispersal.into.realm','dispersal.into.hot'),col=c('blue','red','green','orange'),lty=1,cex=.7,bty='n')
  }
  
  #plot BSM dispersal rates hot into realm
  
  ##plot BAMM in situ rates
  #par(new=T)
  #plot(c(1,1),xlim=c(0,length(insitu.hot.rates.BAMM.CI)),ylim=c(0,max(c(unlist(insitu.hot.rates.BAMM.CI),unlist(insitu.realm.rates.BAMM.CI)))),type='n',xaxt='n',xlab='age',yaxt='n')
  ##polygon(x=c(c(age:length(insitu.hot.rates.BAMM.CI)),c(length(insitu.hot.rates.BAMM.CI):age)),y=c(unlist(lapply(insitu.hot.rates.BAMM.CI[c(age:length(insitu.hot.rates.BAMM.CI))],function(x)x[1])),rev(unlist(lapply(insitu.hot.rates.BAMM.CI[c(age:length(insitu.hot.rates.BAMM.CI))],function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
  #lines(x=c(1:length(insitu.hot.rates.BAMM.CI)),y=unlist(lapply(insitu.hot.rates.BAMM.CI,function(x)x[3])),col=adjustcolor( "red", alpha.f = 0.9))
  #polygon(x=c((1:length(insitu.hot.rates.BAMM.CI)),(rev(1:length(insitu.hot.rates.BAMM.CI)))),y=c((unlist(lapply(insitu.hot.rates.BAMM.CI,function(x)x[1]))),rev(unlist(lapply(insitu.hot.rates.BAMM.CI,function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
  #lines(x=c(1:length(insitu.realm.rates.BAMM.CI)),y=unlist(lapply(insitu.realm.rates.BAMM.CI,function(x)x[3])),col=adjustcolor( "blue", alpha.f = 0.9))
  #polygon(x=c((1:length(insitu.realm.rates.BAMM.CI)),(rev(1:length(insitu.realm.rates.BAMM.CI)))),y=c((unlist(lapply(insitu.realm.rates.BAMM.CI,function(x)x[1]))),rev(unlist(lapply(insitu.realm.rates.BAMM.CI,function(x)x[2])))),col=adjustcolor( "blue", alpha.f = 0.2))
  ##axis(1,at=c(0:length(BSMoutput.real.rates[['insitu.realm.rates.BSM.CI']])),labels=c((length(BSMoutput.real.rates[['insitu.realm.rates.BSM.CI']])+2):2))
  #axis(4)
  #legend('topright',c('insitu.hot','insitu.realm'),col=c('red','blue'),lty=1,cex=.7,bty='n')
}

plot_rates_events_real_2mybin<-function(BSMoutput.real.rates,BSMoutput.real.events,age,name,plots){
  #correct age
  age<-(length(BSMoutput.real.rates[['insitu.realm.rates.BSM.CI']])*2-age)/2
  #trim lists to age period
  BSMoutput.real.rates<-lapply(BSMoutput.real.rates,function(x)x[age:length(x)])
  BSMoutput.real.events<-lapply(BSMoutput.real.events,function(x)x[age:length(x)])
  #delete Inf in rates and events
  BSMoutput.real.rates<-lapply(BSMoutput.real.rates,function(x) lapply(x,function(x) {finite.hit<-which(!is.finite(x));if(length(finite.hit)==0){return(x)}else{if(finite.hit==2){x<-c(x[1],x[3],x[3])};return(x)}}))
  BSMoutput.real.events<-lapply(BSMoutput.real.events,function(x) lapply(x,function(x) {finite.hit<-which(!is.finite(x));if(length(finite.hit)==0){return(x)}else{if(finite.hit==2){x<-c(x[1],x[3],x[3])};return(x)}}))
  #remove the last bin - two million year time bin
  BSMoutput.real.rates<-lapply(BSMoutput.real.rates,function(x) x<-x[-c(length(x))])
  BSMoutput.real.events<-lapply(BSMoutput.real.events,function(x) x<-x[-c((length(x)-1),length(x))])
  #assign variables for easier handling
  hot.into.realm.rates.BSM.CI<-BSMoutput.real.rates[['hot.into.realm.rates.BSM.CI']]
  realm.into.hot.rates.BSM.CI<-BSMoutput.real.rates[['realm.into.hot.rates.BSM.CI']]
  insitu.realm.rates.BAMM.CI<-BSMoutput.real.rates[['insitu.realm.rates.BAMM.CI']]
  insitu.hot.rates.BAMM.CI<-BSMoutput.real.rates[['insitu.hot.rates.BAMM.CI']]
  insitu.realm.rates.BSM.CI<-BSMoutput.real.rates[['insitu.realm.rates.BSM.CI']]
  insitu.hot.rates.BSM.CI<-BSMoutput.real.rates[['insitu.hot.rates.BSM.CI']]
  dispersal.into.realm.rates.BSM.CI<-BSMoutput.real.rates[['dispersal.into.realm.rates.BSM.CI']]
  dispersal.into.hot.rates.BSM.CI<-BSMoutput.real.rates[['dispersal.into.hot.rates.BSM.CI']]
  hot.into.realm.time.bins.counts.rev.cumsum.CI<-BSMoutput.real.events[['hot.into.realm.time.bins.counts.rev.cumsum.CI']]
  realm.into.hot.time.bins.counts.rev.cumsum.CI<-BSMoutput.real.events[['realm.into.hot.time.bins.counts.rev.cumsum.CI']]
  insitu.realm.time.bins.counts.rev.cumsum.CI<-BSMoutput.real.events[['insitu.realm.time.bins.counts.rev.cumsum.CI']]
  insitu.hot.time.bins.counts.rev.cumsum.CI<-BSMoutput.real.events[['insitu.hot.time.bins.counts.rev.cumsum.CI']]
  dispersal.into.realm.time.bins.counts.rev.cumsum.CI<-BSMoutput.real.events[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']]
  dispersal.into.hot.time.bins.counts.rev.cumsum.CI<-BSMoutput.real.events[['dispersal.into.hot.time.bins.counts.rev.cumsum.CI']]
  dispersal.from.realm.time.bins.counts.rev.cumsum.CI<-BSMoutput.real.events[['dispersal.from.realm.time.bins.counts.rev.cumsum.CI']]
  dispersal.from.hot.time.bins.counts.rev.cumsum.CI<-BSMoutput.real.events[['dispersal.from.hot.time.bins.counts.rev.cumsum.CI']]
  if(plots=='difference.dispersal.rates'){
    difference.rates.hotrealm<-unlist(lapply(c(1:length(hot.into.realm.rates.BSM.CI)),function(x){hot.into.realm.rates.BSM.CI[[x]][3]-realm.into.hot.rates.BSM.CI[[x]][[3]]}))
    plot(x=c(1:length(hot.into.realm.rates.BSM.CI)),y=diff,ylim=c(-max(abs(difference.rates.hotrealm)),max(abs(difference.rates.hotrealm))),type='n',xaxt='n',xlab='age',ylab='diff. dispersal hot.to.realm - realm.to.hot',main=paste(name,' diff in dispersal rates',sep=''))
    #plot(c(1,1),xlim=c(0,length(hot.into.realm.rates.BSM.CI)),ylim=c(0,max(c(unlist(hot.into.realm.rates.BSM.CI),unlist(realm.into.hot.rates.BSM.CI))[is.finite(c(unlist(hot.into.realm.rates.BSM.CI),unlist(realm.into.hot.rates.BSM.CI)))])),type='n',xaxt='n',xlab='age',ylab='dispersal rates',main=paste(name,' dispersal rates',sep=''))
    lines(x=c(1:length(hot.into.realm.rates.BSM.CI)),y=difference.rates.hotrealm)
    abline(h=0,lty=2)
    axis(1,at=c(1:length(realm.into.hot.rates.BSM.CI)),labels=c(seq(from=length(realm.into.hot.rates.BSM.CI)*2,to=2,by=-2)))
    
  }
  if(plots=='dispersal.rates'){
    plot(c(1,1),xlim=c(0,length(hot.into.realm.rates.BSM.CI)),ylim=c(0,max(c(unlist(hot.into.realm.rates.BSM.CI),unlist(realm.into.hot.rates.BSM.CI))[is.finite(c(unlist(hot.into.realm.rates.BSM.CI),unlist(realm.into.hot.rates.BSM.CI)))])),type='n',xaxt='n',xlab='age',ylab='dispersal rates',main=paste(name,' dispersal rates',sep=''))
    lines(x=c(1:length(hot.into.realm.rates.BSM.CI)),y=unlist(lapply(hot.into.realm.rates.BSM.CI,function(x)x[3])),col=adjustcolor( "red", alpha.f = 0.9))
    polygon(x=c((1:length(hot.into.realm.rates.BSM.CI)),(rev(1:length(hot.into.realm.rates.BSM.CI)))),y=c((unlist(lapply(hot.into.realm.rates.BSM.CI,function(x)x[1]))),rev(unlist(lapply(hot.into.realm.rates.BSM.CI,function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
    lines(x=c(1:length(realm.into.hot.rates.BSM.CI)),y=unlist(lapply(realm.into.hot.rates.BSM.CI,function(x)x[3])),col=adjustcolor( "blue", alpha.f = 0.9))
    polygon(x=c((1:length(realm.into.hot.rates.BSM.CI)),(rev(1:length(realm.into.hot.rates.BSM.CI)))),y=c((unlist(lapply(realm.into.hot.rates.BSM.CI,function(x)x[1]))),rev(unlist(lapply(realm.into.hot.rates.BSM.CI,function(x)x[2])))),col=adjustcolor( "blue", alpha.f = 0.2))
    axis(1,at=c(1:length(realm.into.hot.rates.BSM.CI)),labels=c(seq(from=length(realm.into.hot.rates.BSM.CI)*2,to=2,by=-2)))
    legend('topleft',c('hot.into.realm','realm.into.hot'),col=c('red','blue'),lty=1,cex=.7,bty='n')
  }else if(plots=='cladogenesis.BAMM.rates'){
    #plot BAMM in situ rates
    plot(c(1,1),xlim=c(0,length(insitu.hot.rates.BAMM.CI)),ylim=c(0,max(c(unlist(insitu.hot.rates.BAMM.CI),unlist(insitu.realm.rates.BAMM.CI)))),type='n',xaxt='n',xlab='age',ylab='insitu.rates',main=paste(name,' cladogenesis BAMM rates',sep=''))
    lines(x=c(1:length(insitu.hot.rates.BAMM.CI)),y=unlist(lapply(insitu.hot.rates.BAMM.CI,function(x)x[3])),col=adjustcolor( "red", alpha.f = 0.9))
    polygon(x=c((1:length(insitu.hot.rates.BAMM.CI)),(rev(1:length(insitu.hot.rates.BAMM.CI)))),y=c((unlist(lapply(insitu.hot.rates.BAMM.CI,function(x)x[1]))),rev(unlist(lapply(insitu.hot.rates.BAMM.CI,function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
    lines(x=c(1:length(insitu.realm.rates.BAMM.CI)),y=unlist(lapply(insitu.realm.rates.BAMM.CI,function(x)x[3])),col=adjustcolor( "blue", alpha.f = 0.9))
    polygon(x=c((1:length(insitu.realm.rates.BAMM.CI)),(rev(1:length(insitu.realm.rates.BAMM.CI)))),y=c((unlist(lapply(insitu.realm.rates.BAMM.CI,function(x)x[1]))),rev(unlist(lapply(insitu.realm.rates.BAMM.CI,function(x)x[2])))),col=adjustcolor( "blue", alpha.f = 0.2))
    axis(1,at=c(1:length(insitu.realm.rates.BAMM.CI)),labels=c(seq(from=length(insitu.realm.rates.BAMM.CI)*2,to=2,by=-2)))
    legend('topleft',c('insitu.hot','insitu.realm'),col=c('red','blue'),lty=1,cex=.7,bty='n')
  }else if(plots=='cladogenesis.BSM.rates'){
    #plot BSM in situ rates
    plot(c(1,1),xlim=c(0,length(insitu.hot.rates.BSM.CI)),ylim=c(0,max(c(unlist(insitu.hot.rates.BSM.CI),unlist(insitu.realm.rates.BSM.CI))[is.finite(c(unlist(insitu.hot.rates.BSM.CI),unlist(insitu.realm.rates.BSM.CI)))])),type='n',xaxt='n',xlab='age',ylab='insitu.rates',main=paste(name,' cladogenesis BSM rates',sep=''))
    lines(x=c(1:length(insitu.hot.rates.BSM.CI)),y=unlist(lapply(insitu.hot.rates.BSM.CI,function(x)x[3])),col=adjustcolor( "red", alpha.f = 0.9))
    polygon(x=c((1:length(insitu.hot.rates.BSM.CI)),(rev(1:length(insitu.hot.rates.BSM.CI)))),y=c((unlist(lapply(insitu.hot.rates.BSM.CI,function(x)x[1]))),rev(unlist(lapply(insitu.hot.rates.BSM.CI,function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
    lines(x=c(1:length(insitu.realm.rates.BSM.CI)),y=unlist(lapply(insitu.realm.rates.BSM.CI,function(x)x[3])),col=adjustcolor( "blue", alpha.f = 0.9))
    polygon(x=c((1:length(insitu.realm.rates.BSM.CI)),(rev(1:length(insitu.realm.rates.BSM.CI)))),y=c((unlist(lapply(insitu.realm.rates.BSM.CI,function(x)x[1]))),rev(unlist(lapply(insitu.realm.rates.BSM.CI,function(x)x[2])))),col=adjustcolor( "blue", alpha.f = 0.2))
    axis(1,at=c(1:length(insitu.realm.rates.BSM.CI)),labels=c(seq(from=length(insitu.realm.rates.BSM.CI)*2,to=2,by=-2)))
    legend('topleft',c('insitu.hot','insitu.realm'),col=c('red','blue'),lty=1,cex=.7,bty='n')
  }else if(plots=='dispersal.events'){
    plot(c(1,1),xlim=c(0,length(hot.into.realm.time.bins.counts.rev.cumsum.CI)),ylim=c(0,max(c(unlist(hot.into.realm.time.bins.counts.rev.cumsum.CI),unlist(realm.into.hot.time.bins.counts.rev.cumsum.CI)))),type='n',xaxt='n',xlab='age',ylab='cumulative # of events',main=paste(name,' dispersal events',sep=''))
    lines(x=c(1:length(hot.into.realm.time.bins.counts.rev.cumsum.CI)),y=unlist(lapply(hot.into.realm.time.bins.counts.rev.cumsum.CI,function(x)x[3])),col=adjustcolor( "red", alpha.f = 0.9))
    polygon(x=c((1:length(hot.into.realm.time.bins.counts.rev.cumsum.CI)),(rev(1:length(hot.into.realm.time.bins.counts.rev.cumsum.CI)))),y=c((unlist(lapply(hot.into.realm.time.bins.counts.rev.cumsum.CI,function(x)x[1]))),rev(unlist(lapply(hot.into.realm.time.bins.counts.rev.cumsum.CI,function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
    lines(x=c(1:length(realm.into.hot.time.bins.counts.rev.cumsum.CI)),y=unlist(lapply(realm.into.hot.time.bins.counts.rev.cumsum.CI,function(x)x[3])),col=adjustcolor( "blue", alpha.f = 0.9))
    polygon(x=c((1:length(realm.into.hot.time.bins.counts.rev.cumsum.CI)),(rev(1:length(realm.into.hot.time.bins.counts.rev.cumsum.CI)))),y=c((unlist(lapply(realm.into.hot.time.bins.counts.rev.cumsum.CI,function(x)x[1]))),rev(unlist(lapply(realm.into.hot.time.bins.counts.rev.cumsum.CI,function(x)x[2])))),col=adjustcolor( "blue", alpha.f = 0.2))
    axis(1,at=c(1:length(realm.into.hot.time.bins.counts.rev.cumsum.CI)),labels=c(seq(from=length(realm.into.hot.time.bins.counts.rev.cumsum.CI)*2,to=2,by=-2)))
    
    legend('topleft',c('hot.into.realm','realm.into.hot'),col=c('red','blue'),lty=1,cex=.7,bty='n')
    
  }else if(plots=='cladogenesis.events'){
    plot(c(1,1),xlim=c(0,length(insitu.hot.time.bins.counts.rev.cumsum.CI)),ylim=c(0,max(c(unlist(insitu.hot.time.bins.counts.rev.cumsum.CI),unlist(insitu.realm.time.bins.counts.rev.cumsum.CI)))),type='n',xaxt='n',xlab='age',ylab='cumulative # of events',main=paste(name,' cladogenetic events',sep=''))
    lines(x=c(1:length(insitu.hot.time.bins.counts.rev.cumsum.CI)),y=unlist(lapply(insitu.hot.time.bins.counts.rev.cumsum.CI,function(x)x[3])),col=adjustcolor( "red", alpha.f = 0.9))
    polygon(x=c((1:length(insitu.hot.time.bins.counts.rev.cumsum.CI)),(rev(1:length(insitu.hot.time.bins.counts.rev.cumsum.CI)))),y=c((unlist(lapply(insitu.hot.time.bins.counts.rev.cumsum.CI,function(x)x[1]))),rev(unlist(lapply(insitu.hot.time.bins.counts.rev.cumsum.CI,function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
    lines(x=c(1:length(insitu.realm.time.bins.counts.rev.cumsum.CI)),y=unlist(lapply(insitu.realm.time.bins.counts.rev.cumsum.CI,function(x)x[3])),col=adjustcolor( "blue", alpha.f = 0.9))
    polygon(x=c((1:length(insitu.realm.time.bins.counts.rev.cumsum.CI)),(rev(1:length(insitu.realm.time.bins.counts.rev.cumsum.CI)))),y=c((unlist(lapply(insitu.realm.time.bins.counts.rev.cumsum.CI,function(x)x[1]))),rev(unlist(lapply(insitu.realm.time.bins.counts.rev.cumsum.CI,function(x)x[2])))),col=adjustcolor( "blue", alpha.f = 0.2))
    axis(1,at=c(1:length(insitu.realm.time.bins.counts.rev.cumsum.CI)),labels=c(seq(from=length(insitu.realm.time.bins.counts.rev.cumsum.CI)*2,to=2,by=-2)))
    legend('topleft',c('insitu.hot','insitu.realm'),col=c('red','blue'),lty=1,cex=.7,bty='n')
  }else if (plots=='dispersal.all.rates'){
    plot(c(1,1),xlim=c(0,length(dispersal.into.realm.rates.BSM.CI)),ylim=c(0,max(c(unlist(dispersal.into.realm.rates.BSM.CI),unlist(dispersal.into.hot.rates.BSM.CI))[is.finite(c(unlist(dispersal.into.realm.rates.BSM.CI),unlist(dispersal.into.hot.rates.BSM.CI)))])),type='n',xaxt='n',xlab='age',ylab='dispersal rates',main=paste(name,' dispersal rates',sep=''))
    lines(x=c(1:length(dispersal.into.realm.rates.BSM.CI)),y=unlist(lapply(dispersal.into.realm.rates.BSM.CI,function(x)x[3])),col=adjustcolor( "red", alpha.f = 0.9))
    polygon(x=c((1:length(dispersal.into.realm.rates.BSM.CI)),(rev(1:length(dispersal.into.realm.rates.BSM.CI)))),y=c((unlist(lapply(dispersal.into.realm.rates.BSM.CI,function(x)x[1]))),rev(unlist(lapply(dispersal.into.realm.rates.BSM.CI,function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
    lines(x=c(1:length(dispersal.into.hot.rates.BSM.CI)),y=unlist(lapply(dispersal.into.hot.rates.BSM.CI,function(x)x[3])),col=adjustcolor( "blue", alpha.f = 0.9))
    polygon(x=c((1:length(dispersal.into.hot.rates.BSM.CI)),(rev(1:length(dispersal.into.hot.rates.BSM.CI)))),y=c((unlist(lapply(dispersal.into.hot.rates.BSM.CI,function(x)x[1]))),rev(unlist(lapply(dispersal.into.hot.rates.BSM.CI,function(x)x[2])))),col=adjustcolor( "blue", alpha.f = 0.2))
    axis(1,at=c(1:length(dispersal.into.hot.rates.BSM.CI)),labels=c(seq(from=length(dispersal.into.hot.rates.BSM.CI)*2,to=2,by=-2)))
    
    legend('topleft',c('dispersal.into.hot','dispersal.into.realm'),col=c('blue','red'),lty=1,cex=.7,bty='n')
  }else if(plots=='dispersal.all.events'){
    plot(c(1,1),xlim=c(0,length(dispersal.into.realm.time.bins.counts.rev.cumsum.CI)),ylim=c(0,max(c(unlist(dispersal.into.realm.time.bins.counts.rev.cumsum.CI),unlist(dispersal.into.hot.time.bins.counts.rev.cumsum.CI))[is.finite(c(unlist(dispersal.into.realm.time.bins.counts.rev.cumsum.CI),unlist(dispersal.into.hot.time.bins.counts.rev.cumsum.CI)))])),type='n',xaxt='n',xlab='age',ylab='cumulative # of events',main=paste(name,' dispersal events',sep=''))
    lines(x=c(1:length(dispersal.into.realm.time.bins.counts.rev.cumsum.CI)),y=unlist(lapply(dispersal.into.realm.time.bins.counts.rev.cumsum.CI,function(x)x[3])),col=adjustcolor( "red", alpha.f = 0.9))
    polygon(x=c((1:length(dispersal.into.realm.time.bins.counts.rev.cumsum.CI)),(rev(1:length(dispersal.into.realm.time.bins.counts.rev.cumsum.CI)))),y=c((unlist(lapply(dispersal.into.realm.time.bins.counts.rev.cumsum.CI,function(x)x[1]))),rev(unlist(lapply(dispersal.into.realm.time.bins.counts.rev.cumsum.CI,function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
    lines(x=c(1:length(dispersal.into.hot.time.bins.counts.rev.cumsum.CI)),y=unlist(lapply(dispersal.into.hot.time.bins.counts.rev.cumsum.CI,function(x)x[3])),col=adjustcolor( "blue", alpha.f = 0.9))
    polygon(x=c((1:length(dispersal.into.hot.time.bins.counts.rev.cumsum.CI)),(rev(1:length(dispersal.into.hot.time.bins.counts.rev.cumsum.CI)))),y=c((unlist(lapply(dispersal.into.hot.time.bins.counts.rev.cumsum.CI,function(x)x[1]))),rev(unlist(lapply(dispersal.into.hot.time.bins.counts.rev.cumsum.CI,function(x)x[2])))),col=adjustcolor( "blue", alpha.f = 0.2))
    axis(1,at=c(1:length(dispersal.into.hot.time.bins.counts.rev.cumsum.CI)),labels=c(seq(from=length(dispersal.into.hot.time.bins.counts.rev.cumsum.CI)*2,to=2,by=-2)))
    legend('topleft',c('dispersal.into.realm','dispersal.into.hot'),col=c('red','blue'),lty=1,cex=.7,bty='n')
  }else if(plots=='dispersal.into.from.events'){
    plot(c(1,1),xlim=c(0,length(dispersal.into.realm.time.bins.counts.rev.cumsum.CI)),ylim=c(0,max(c(unlist(dispersal.into.realm.time.bins.counts.rev.cumsum.CI),unlist(dispersal.into.hot.time.bins.counts.rev.cumsum.CI),unlist(dispersal.from.hot.time.bins.counts.rev.cumsum.CI),unlist(dispersal.from.realm.time.bins.counts.rev.cumsum.CI)))),type='n',xaxt='n',xlab='age',ylab='cumulative # of events',main=paste(name,' dispersal events',sep=''))
    lines(x=c(1:length(dispersal.into.realm.time.bins.counts.rev.cumsum.CI)),y=unlist(lapply(dispersal.into.realm.time.bins.counts.rev.cumsum.CI,function(x)x[3])),col=adjustcolor( "green", alpha.f = 0.9))
    polygon(x=c((1:length(dispersal.into.realm.time.bins.counts.rev.cumsum.CI)),(rev(1:length(dispersal.into.realm.time.bins.counts.rev.cumsum.CI)))),y=c((unlist(lapply(dispersal.into.realm.time.bins.counts.rev.cumsum.CI,function(x)x[1]))),rev(unlist(lapply(dispersal.into.realm.time.bins.counts.rev.cumsum.CI,function(x)x[2])))),col=adjustcolor( "green", alpha.f = 0.2))
    lines(x=c(1:length(dispersal.into.hot.time.bins.counts.rev.cumsum.CI)),y=unlist(lapply(dispersal.into.hot.time.bins.counts.rev.cumsum.CI,function(x)x[3])),col=adjustcolor( "orange", alpha.f = 0.9))
    polygon(x=c((1:length(dispersal.into.hot.time.bins.counts.rev.cumsum.CI)),(rev(1:length(dispersal.into.hot.time.bins.counts.rev.cumsum.CI)))),y=c((unlist(lapply(dispersal.into.hot.time.bins.counts.rev.cumsum.CI,function(x)x[1]))),rev(unlist(lapply(dispersal.into.hot.time.bins.counts.rev.cumsum.CI,function(x)x[2])))),col=adjustcolor( "orange", alpha.f = 0.2))
    lines(x=c(1:length(dispersal.from.realm.time.bins.counts.rev.cumsum.CI)),y=unlist(lapply(dispersal.from.realm.time.bins.counts.rev.cumsum.CI,function(x)x[3])),col=adjustcolor( "blue", alpha.f = 0.9))
    polygon(x=c((1:length(dispersal.from.realm.time.bins.counts.rev.cumsum.CI)),(rev(1:length(dispersal.from.realm.time.bins.counts.rev.cumsum.CI)))),y=c((unlist(lapply(dispersal.from.realm.time.bins.counts.rev.cumsum.CI,function(x)x[1]))),rev(unlist(lapply(dispersal.from.realm.time.bins.counts.rev.cumsum.CI,function(x)x[2])))),col=adjustcolor( "blue", alpha.f = 0.2))
    lines(x=c(1:length(dispersal.from.hot.time.bins.counts.rev.cumsum.CI)),y=unlist(lapply(dispersal.from.hot.time.bins.counts.rev.cumsum.CI,function(x)x[3])),col=adjustcolor( "red", alpha.f = 0.9))
    polygon(x=c((1:length(dispersal.from.hot.time.bins.counts.rev.cumsum.CI)),(rev(1:length(dispersal.from.hot.time.bins.counts.rev.cumsum.CI)))),y=c((unlist(lapply(dispersal.from.hot.time.bins.counts.rev.cumsum.CI,function(x)x[1]))),rev(unlist(lapply(dispersal.from.hot.time.bins.counts.rev.cumsum.CI,function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
    axis(1,at=c(1:length(dispersal.from.hot.time.bins.counts.rev.cumsum.CI)),labels=c(seq(from=length(dispersal.from.hot.time.bins.counts.rev.cumsum.CI)*2,to=2,by=-2)))
    legend('topleft',c('dispersal.from.realm','dispersal.from.hot','dispersal.into.realm','dispersal.into.hot'),col=c('blue','red','green','orange'),lty=1,cex=.7,bty='n')
  }
  
  #plot BSM dispersal rates hot into realm
  
  ##plot BAMM in situ rates
  #par(new=T)
  #plot(c(1,1),xlim=c(0,length(insitu.hot.rates.BAMM.CI)),ylim=c(0,max(c(unlist(insitu.hot.rates.BAMM.CI),unlist(insitu.realm.rates.BAMM.CI)))),type='n',xaxt='n',xlab='age',yaxt='n')
  ##polygon(x=c(c(age:length(insitu.hot.rates.BAMM.CI)),c(length(insitu.hot.rates.BAMM.CI):age)),y=c(unlist(lapply(insitu.hot.rates.BAMM.CI[c(age:length(insitu.hot.rates.BAMM.CI))],function(x)x[1])),rev(unlist(lapply(insitu.hot.rates.BAMM.CI[c(age:length(insitu.hot.rates.BAMM.CI))],function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
  #lines(x=c(1:length(insitu.hot.rates.BAMM.CI)),y=unlist(lapply(insitu.hot.rates.BAMM.CI,function(x)x[3])),col=adjustcolor( "red", alpha.f = 0.9))
  #polygon(x=c((1:length(insitu.hot.rates.BAMM.CI)),(rev(1:length(insitu.hot.rates.BAMM.CI)))),y=c((unlist(lapply(insitu.hot.rates.BAMM.CI,function(x)x[1]))),rev(unlist(lapply(insitu.hot.rates.BAMM.CI,function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
  #lines(x=c(1:length(insitu.realm.rates.BAMM.CI)),y=unlist(lapply(insitu.realm.rates.BAMM.CI,function(x)x[3])),col=adjustcolor( "blue", alpha.f = 0.9))
  #polygon(x=c((1:length(insitu.realm.rates.BAMM.CI)),(rev(1:length(insitu.realm.rates.BAMM.CI)))),y=c((unlist(lapply(insitu.realm.rates.BAMM.CI,function(x)x[1]))),rev(unlist(lapply(insitu.realm.rates.BAMM.CI,function(x)x[2])))),col=adjustcolor( "blue", alpha.f = 0.2))
  ##axis(1,at=c(0:length(BSMoutput.real.rates[['insitu.realm.rates.BSM.CI']])),labels=c((length(BSMoutput.real.rates[['insitu.realm.rates.BSM.CI']])+2):2))
  #axis(4)
  #legend('topright',c('insitu.hot','insitu.realm'),col=c('red','blue'),lty=1,cex=.7,bty='n')
}

  #

plot_rates_events_realvscontrol<-function(BSMoutput.real.rates,BSMoutput.real.events,BSMoutput.control.rates,BSMoutput.control.events,age,name,plots){
  #correct age
  age<-length(BSMoutput.real.rates[['insitu.realm.rates.BSM.CI']])-age
  #trim lists to age period
  BSMoutput.real.rates<-lapply(BSMoutput.real.rates,function(x)x[age:length(x)])
  BSMoutput.real.events<-lapply(BSMoutput.real.events,function(x)x[age:length(x)])
  BSMoutput.control.rates<-lapply(BSMoutput.control.rates,function(x) lapply(x,function(x) x[age:length(x)]))
  BSMoutput.control.events<-lapply(BSMoutput.control.events,function(x) lapply(x,function(x) x[age:length(x)]))
  #remove last two million year time bins
  BSMoutput.real.rates<-lapply(BSMoutput.real.rates,function(x) x<-x[-c((length(x)-1),length(x))])
  BSMoutput.real.events<-lapply(BSMoutput.real.events,function(x) x<-x[-c((length(x)-2),(length(x)-1),length(x))])
  BSMoutput.control.rates<-lapply(BSMoutput.control.rates,function(x) lapply(x,function(x) x<-x[-c((length(x)-1),length(x))]))
  BSMoutput.control.events<-lapply(BSMoutput.control.events,function(x) lapply(x,function(x) x<-x[-c((length(x)-2),(length(x)-1),length(x))]))
  #assign variables for easier handling
  #real
  hot.into.realm.rates.BSM.CI<-BSMoutput.real.rates[['hot.into.realm.rates.BSM.CI']]
  realm.into.hot.rates.BSM.CI<-BSMoutput.real.rates[['realm.into.hot.rates.BSM.CI']]
  insitu.realm.rates.BAMM.CI<-BSMoutput.real.rates[['insitu.realm.rates.BAMM.CI']]
  insitu.hot.rates.BAMM.CI<-BSMoutput.real.rates[['insitu.hot.rates.BAMM.CI']]
  insitu.realm.rates.BSM.CI<-BSMoutput.real.rates[['insitu.realm.rates.BSM.CI']]
  insitu.hot.rates.BSM.CI<-BSMoutput.real.rates[['insitu.hot.rates.BSM.CI']]
  dispersal.into.realm.rates.BSM.CI<-BSMoutput.real.rates[['dispersal.into.realm.rates.BSM.CI']]
  dispersal.into.hot.rates.BSM.CI<-BSMoutput.real.rates[['dispersal.into.hot.rates.BSM.CI']]
  hot.into.realm.time.bins.counts.rev.cumsum.CI<-BSMoutput.real.events[['hot.into.realm.time.bins.counts.rev.cumsum.CI']]
  realm.into.hot.time.bins.counts.rev.cumsum.CI<-BSMoutput.real.events[['realm.into.hot.time.bins.counts.rev.cumsum.CI']]
  insitu.realm.time.bins.counts.rev.cumsum.CI<-BSMoutput.real.events[['insitu.realm.time.bins.counts.rev.cumsum.CI']]
  insitu.hot.time.bins.counts.rev.cumsum.CI<-BSMoutput.real.events[['insitu.hot.time.bins.counts.rev.cumsum.CI']]
  dispersal.into.realm.time.bins.counts.rev.cumsum.CI<-BSMoutput.real.events[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']]
  dispersal.into.hot.time.bins.counts.rev.cumsum.CI<-BSMoutput.real.events[['dispersal.into.hot.time.bins.counts.rev.cumsum.CI']]
  dispersal.from.realm.time.bins.counts.rev.cumsum.CI<-BSMoutput.real.events[['dispersal.from.realm.time.bins.counts.rev.cumsum.CI']]
  dispersal.from.hot.time.bins.counts.rev.cumsum.CI<-BSMoutput.real.events[['dispersal.from.hot.time.bins.counts.rev.cumsum.CI']]
  #controls
  hot.into.realm.rates.BSM.CI.control<-lapply(BSMoutput.control.rates,function(x)x[['hot.into.realm.rates.BSM.CI']])
  realm.into.hot.rates.BSM.CI.control<-lapply(BSMoutput.control.rates,function(x) x[['realm.into.hot.rates.BSM.CI']])
  insitu.realm.rates.BAMM.CI.control<-lapply(BSMoutput.control.rates,function(x) x[['insitu.realm.rates.BAMM.CI']])
  insitu.hot.rates.BAMM.CI.control<-lapply(BSMoutput.control.rates,function(x) x[['insitu.hot.rates.BAMM.CI']])
  insitu.realm.rates.BSM.CI.control<-lapply(BSMoutput.control.rates,function(x) x[['insitu.realm.rates.BSM.CI']])
  insitu.hot.rates.BSM.CI.control<-lapply(BSMoutput.control.rates,function(x) x[['insitu.hot.rates.BSM.CI']])
  dispersal.into.realm.rates.BSM.CI.control<-lapply(BSMoutput.control.rates,function(x) x[['dispersal.into.realm.rates.BSM.CI']])
  dispersal.into.hot.rates.BSM.CI.control<-lapply(BSMoutput.control.rates,function(x) x[['dispersal.into.hot.rates.BSM.CI']])
  hot.into.realm.time.bins.counts.rev.cumsum.CI.control<-lapply(BSMoutput.control.events,function(x) x[['hot.into.realm.time.bins.counts.rev.cumsum.CI']])
  realm.into.hot.time.bins.counts.rev.cumsum.CI.control<-lapply(BSMoutput.control.events,function(x) x[['realm.into.hot.time.bins.counts.rev.cumsum.CI']])
  insitu.realm.time.bins.counts.rev.cumsum.CI.control<-lapply(BSMoutput.control.events,function(x) x[['insitu.realm.time.bins.counts.rev.cumsum.CI']])
  insitu.hot.time.bins.counts.rev.cumsum.CI.control<-lapply(BSMoutput.control.events,function(x) x[['insitu.hot.time.bins.counts.rev.cumsum.CI']])
  dispersal.into.realm.time.bins.counts.rev.cumsum.CI.control<-lapply(BSMoutput.control.events,function(x) x[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']])
  dispersal.into.hot.time.bins.counts.rev.cumsum.CI.control<-lapply(BSMoutput.control.events,function(x) x[['dispersal.into.hot.time.bins.counts.rev.cumsum.CI']])
  dispersal.from.realm.time.bins.counts.rev.cumsum.CI.control<-lapply(BSMoutput.control.events,function(x) x[['dispersal.from.realm.time.bins.counts.rev.cumsum.CI']])
  dispersal.from.hot.time.bins.counts.rev.cumsum.CI.control<-lapply(BSMoutput.control.events,function(x) x[['dispersal.from.hot.time.bins.counts.rev.cumsum.CI']])
  
  if(plots=='dispersal.rates'){
    plot(c(1,1),xlim=c(0,length(hot.into.realm.rates.BSM.CI)),ylim=c(0,max(c(unlist(hot.into.realm.rates.BSM.CI),unlist(realm.into.hot.rates.BSM.CI))[is.finite(c(unlist(hot.into.realm.rates.BSM.CI),unlist(realm.into.hot.rates.BSM.CI)))])),type='n',xaxt='n',xlab='age',ylab='dispersal rates',main=paste(name,' dispersal rates',sep=''))
    lines(x=c(1:length(hot.into.realm.rates.BSM.CI)),y=unlist(lapply(hot.into.realm.rates.BSM.CI,function(x)x[3])),col=adjustcolor( "red", alpha.f = 0.9))
    #polygon(x=c((1:length(hot.into.realm.rates.BSM.CI)),(rev(1:length(hot.into.realm.rates.BSM.CI)))),y=c((unlist(lapply(hot.into.realm.rates.BSM.CI,function(x)x[1]))),rev(unlist(lapply(hot.into.realm.rates.BSM.CI,function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
    lines(x=c(1:length(realm.into.hot.rates.BSM.CI)),y=unlist(lapply(realm.into.hot.rates.BSM.CI,function(x)x[3])),col=adjustcolor( "blue", alpha.f = 0.9))
    #polygon(x=c((1:length(realm.into.hot.rates.BSM.CI)),(rev(1:length(realm.into.hot.rates.BSM.CI)))),y=c((unlist(lapply(realm.into.hot.rates.BSM.CI,function(x)x[1]))),rev(unlist(lapply(realm.into.hot.rates.BSM.CI,function(x)x[2])))),col=adjustcolor( "blue", alpha.f = 0.2))
    lapply(hot.into.realm.rates.BSM.CI.control,function(x)lines(x=c(1:length(hot.into.realm.rates.BSM.CI)),y=unlist(lapply(x,function(x)x[3])),col=adjustcolor( "orange", alpha.f = 0.9)))
    lapply(realm.into.hot.rates.BSM.CI.control,function(x)lines(x=c(1:length(hot.into.realm.rates.BSM.CI)),y=unlist(lapply(x,function(x)x[3])),col=adjustcolor( "green", alpha.f = 0.9)))
    axis(1,at=c(0:length(realm.into.hot.rates.BSM.CI)),labels=c((length(realm.into.hot.rates.BSM.CI)+2):2))
    legend('topleft',c('hot.into.realm','realm.into.hot','hot.into.realm.control','realm.into.hot.control'),col=c('red','blue','orange','green'),lty=1,cex=.7,bty='n')
  }else if(plots=='cladogenesis.BAMM.rates'){
    #plot BAMM in situ rates
    plot(c(1,1),xlim=c(0,length(insitu.hot.rates.BAMM.CI)),ylim=c(0,max(c(unlist(insitu.hot.rates.BAMM.CI),unlist(insitu.realm.rates.BAMM.CI)))),type='n',xaxt='n',xlab='age',ylab='insitu.rates',main=paste(name,' cladogenesis BAMM rates',sep=''))
    lines(x=c(1:length(insitu.hot.rates.BAMM.CI)),y=unlist(lapply(insitu.hot.rates.BAMM.CI,function(x)x[3])),col=adjustcolor( "red", alpha.f = 0.9))
    #polygon(x=c((1:length(insitu.hot.rates.BAMM.CI)),(rev(1:length(insitu.hot.rates.BAMM.CI)))),y=c((unlist(lapply(insitu.hot.rates.BAMM.CI,function(x)x[1]))),rev(unlist(lapply(insitu.hot.rates.BAMM.CI,function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
    lines(x=c(1:length(insitu.realm.rates.BAMM.CI)),y=unlist(lapply(insitu.realm.rates.BAMM.CI,function(x)x[3])),col=adjustcolor( "blue", alpha.f = 0.9))
    #polygon(x=c((1:length(insitu.realm.rates.BAMM.CI)),(rev(1:length(insitu.realm.rates.BAMM.CI)))),y=c((unlist(lapply(insitu.realm.rates.BAMM.CI,function(x)x[1]))),rev(unlist(lapply(insitu.realm.rates.BAMM.CI,function(x)x[2])))),col=adjustcolor( "blue", alpha.f = 0.2))
    lapply(insitu.hot.rates.BAMM.CI.control,function(x)lines(x=c(1:length(insitu.hot.rates.BAMM.CI)),y=unlist(lapply(x,function(x)x[3])),col=adjustcolor( "orange", alpha.f = 0.9)))
    lapply(insitu.realm.rates.BAMM.CI.control,function(x)lines(x=c(1:length(insitu.hot.rates.BAMM.CI)),y=unlist(lapply(x,function(x)x[3])),col=adjustcolor( "green", alpha.f = 0.9)))
    axis(1,at=c(0:length(insitu.realm.rates.BAMM.CI)),labels=c((length(insitu.realm.rates.BAMM.CI)+2):2))
    legend('topleft',c('insitu.hot','insitu.realm','insitu.hot.control','insitu.realm.control'),col=c('red','blue','orange','green'),lty=1,cex=.7,bty='n')
  }else if(plots=='cladogenesis.BSM.rates'){
    #plot BSM in situ rates
    plot(c(1,1),xlim=c(0,length(insitu.hot.rates.BSM.CI)),ylim=c(0,max(c(unlist(insitu.hot.rates.BSM.CI),unlist(insitu.realm.rates.BSM.CI))[is.finite(c(unlist(insitu.hot.rates.BSM.CI),unlist(insitu.realm.rates.BSM.CI)))])),type='n',xaxt='n',xlab='age',ylab='insitu.rates',main=paste(name,' cladogenesis BSM rates',sep=''))
    lines(x=c(1:length(insitu.hot.rates.BSM.CI)),y=unlist(lapply(insitu.hot.rates.BSM.CI,function(x)x[3])),col=adjustcolor( "red", alpha.f = 0.9))
    #polygon(x=c((1:length(insitu.hot.rates.BSM.CI)),(rev(1:length(insitu.hot.rates.BSM.CI)))),y=c((unlist(lapply(insitu.hot.rates.BSM.CI,function(x)x[1]))),rev(unlist(lapply(insitu.hot.rates.BSM.CI,function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
    lines(x=c(1:length(insitu.realm.rates.BSM.CI)),y=unlist(lapply(insitu.realm.rates.BSM.CI,function(x)x[3])),col=adjustcolor( "blue", alpha.f = 0.9))
    #polygon(x=c((1:length(insitu.realm.rates.BSM.CI)),(rev(1:length(insitu.realm.rates.BSM.CI)))),y=c((unlist(lapply(insitu.realm.rates.BSM.CI,function(x)x[1]))),rev(unlist(lapply(insitu.realm.rates.BSM.CI,function(x)x[2])))),col=adjustcolor( "blue", alpha.f = 0.2))
    lapply(insitu.hot.rates.BSM.CI.control,function(x)lines(x=c(1:length(insitu.realm.rates.BSM.CI)),y=unlist(lapply(x,function(x)x[3])),col=adjustcolor( "orange", alpha.f = 0.9)))
    lapply(insitu.realm.rates.BSM.CI.control,function(x)lines(x=c(1:length(insitu.realm.rates.BSM.CI)),y=unlist(lapply(x,function(x)x[3])),col=adjustcolor( "green", alpha.f = 0.9)))
    axis(1,at=c(0:length(insitu.realm.rates.BSM.CI)),labels=c((length(insitu.realm.rates.BSM.CI)+2):2))
    legend('topleft',c('insitu.hot','insitu.realm','insitu.hot.control','insitu.realm.control'),col=c('red','blue','orange','green'),lty=1,cex=.7,bty='n')
  }else if(plots=='dispersal.events'){
    plot(c(1,1),xlim=c(0,length(hot.into.realm.time.bins.counts.rev.cumsum.CI)),ylim=c(0,max(c(unlist(hot.into.realm.time.bins.counts.rev.cumsum.CI),unlist(realm.into.hot.time.bins.counts.rev.cumsum.CI)))),type='n',xaxt='n',xlab='age',ylab='cumulative # of events',main=paste(name,' dispersal events',sep=''))
    lines(x=c(1:length(hot.into.realm.time.bins.counts.rev.cumsum.CI)),y=unlist(lapply(hot.into.realm.time.bins.counts.rev.cumsum.CI,function(x)x[3])),col=adjustcolor( "red", alpha.f = 0.9))
    #polygon(x=c((1:length(hot.into.realm.time.bins.counts.rev.cumsum.CI)),(rev(1:length(hot.into.realm.time.bins.counts.rev.cumsum.CI)))),y=c((unlist(lapply(hot.into.realm.time.bins.counts.rev.cumsum.CI,function(x)x[1]))),rev(unlist(lapply(hot.into.realm.time.bins.counts.rev.cumsum.CI,function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
    lines(x=c(1:length(realm.into.hot.time.bins.counts.rev.cumsum.CI)),y=unlist(lapply(realm.into.hot.time.bins.counts.rev.cumsum.CI,function(x)x[3])),col=adjustcolor( "blue", alpha.f = 0.9))
    #polygon(x=c((1:length(realm.into.hot.time.bins.counts.rev.cumsum.CI)),(rev(1:length(realm.into.hot.time.bins.counts.rev.cumsum.CI)))),y=c((unlist(lapply(realm.into.hot.time.bins.counts.rev.cumsum.CI,function(x)x[1]))),rev(unlist(lapply(realm.into.hot.time.bins.counts.rev.cumsum.CI,function(x)x[2])))),col=adjustcolor( "blue", alpha.f = 0.2))
    lapply(hot.into.realm.time.bins.counts.rev.cumsum.CI.control,function(x)lines(x=c(1:length(realm.into.hot.time.bins.counts.rev.cumsum.CI)),y=unlist(lapply(x,function(x)x[3])),col=adjustcolor( "orange", alpha.f = 0.9)))
    lapply(realm.into.hot.time.bins.counts.rev.cumsum.CI.control,function(x)lines(x=c(1:length(realm.into.hot.time.bins.counts.rev.cumsum.CI)),y=unlist(lapply(x,function(x)x[3])),col=adjustcolor( "green", alpha.f = 0.9)))
    axis(1,at=c(0:length(BSMoutput.real.rates[['insitu.realm.rates.BSM.CI']])),labels=c((length(BSMoutput.real.rates[['insitu.realm.rates.BSM.CI']])+2):2))
    legend('topleft',c('hot.into.realm','realm.into.hot','hot.into.realm.control','realm.into.hot.control'),col=c('red','blue','orange','green'),lty=1,cex=.7,bty='n')
  }else if(plots=='cladogenesis.events'){
    plot(c(1,1),xlim=c(0,length(insitu.hot.time.bins.counts.rev.cumsum.CI)),ylim=c(0,max(c(unlist(insitu.hot.time.bins.counts.rev.cumsum.CI),unlist(insitu.realm.time.bins.counts.rev.cumsum.CI)))),type='n',xaxt='n',xlab='age',ylab='cumulative # of events',main=paste(name,' cladogenetic events',sep=''))
    lines(x=c(1:length(insitu.hot.time.bins.counts.rev.cumsum.CI)),y=unlist(lapply(insitu.hot.time.bins.counts.rev.cumsum.CI,function(x)x[3])),col=adjustcolor( "red", alpha.f = 0.9))
    #polygon(x=c((1:length(insitu.hot.time.bins.counts.rev.cumsum.CI)),(rev(1:length(insitu.hot.time.bins.counts.rev.cumsum.CI)))),y=c((unlist(lapply(insitu.hot.time.bins.counts.rev.cumsum.CI,function(x)x[1]))),rev(unlist(lapply(insitu.hot.time.bins.counts.rev.cumsum.CI,function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
    lines(x=c(1:length(insitu.realm.time.bins.counts.rev.cumsum.CI)),y=unlist(lapply(insitu.realm.time.bins.counts.rev.cumsum.CI,function(x)x[3])),col=adjustcolor( "blue", alpha.f = 0.9))
    #polygon(x=c((1:length(insitu.realm.time.bins.counts.rev.cumsum.CI)),(rev(1:length(insitu.realm.time.bins.counts.rev.cumsum.CI)))),y=c((unlist(lapply(insitu.realm.time.bins.counts.rev.cumsum.CI,function(x)x[1]))),rev(unlist(lapply(insitu.realm.time.bins.counts.rev.cumsum.CI,function(x)x[2])))),col=adjustcolor( "blue", alpha.f = 0.2))
    lapply(insitu.hot.time.bins.counts.rev.cumsum.CI.control,function(x)lines(x=c(1:length(insitu.hot.time.bins.counts.rev.cumsum.CI)),y=unlist(lapply(x,function(x)x[3])),col=adjustcolor( "orange", alpha.f = 0.9)))
    lapply(insitu.realm.time.bins.counts.rev.cumsum.CI.control,function(x)lines(x=c(1:length(insitu.hot.time.bins.counts.rev.cumsum.CI)),y=unlist(lapply(x,function(x)x[3])),col=adjustcolor( "green", alpha.f = 0.9)))
    axis(1,at=c(0:length(BSMoutput.real.rates[['insitu.realm.rates.BSM.CI']])),labels=c((length(BSMoutput.real.rates[['insitu.realm.rates.BSM.CI']])+2):2))
    legend('topleft',c('insitu.hot','insitu.realm','insitu.hot.control','insitu.realm.control'),col=c('red','blue','orange','green'),lty=1,cex=.7,bty='n')
  }else if (plots=='dispersal.all.rates'){
    plot(c(1,1),xlim=c(0,length(dispersal.into.realm.rates.BSM.CI)),ylim=c(0,max(c(unlist(dispersal.into.realm.rates.BSM.CI),unlist(dispersal.into.hot.rates.BSM.CI))[is.finite(c(unlist(dispersal.into.realm.rates.BSM.CI),unlist(dispersal.into.hot.rates.BSM.CI)))])),type='n',xaxt='n',xlab='age',ylab='dispersal rates',main=paste(name,' dispersal rates',sep=''))
    lines(x=c(1:length(dispersal.into.realm.rates.BSM.CI)),y=unlist(lapply(dispersal.into.realm.rates.BSM.CI,function(x)x[3])),col=adjustcolor( "red", alpha.f = 0.9))
    #polygon(x=c((1:length(dispersal.into.realm.rates.BSM.CI)),(rev(1:length(dispersal.into.realm.rates.BSM.CI)))),y=c((unlist(lapply(dispersal.into.realm.rates.BSM.CI,function(x)x[1]))),rev(unlist(lapply(dispersal.into.realm.rates.BSM.CI,function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
    lines(x=c(1:length(dispersal.into.hot.rates.BSM.CI)),y=unlist(lapply(dispersal.into.hot.rates.BSM.CI,function(x)x[3])),col=adjustcolor( "blue", alpha.f = 0.9))
    #polygon(x=c((1:length(dispersal.into.hot.rates.BSM.CI)),(rev(1:length(dispersal.into.hot.rates.BSM.CI)))),y=c((unlist(lapply(dispersal.into.hot.rates.BSM.CI,function(x)x[1]))),rev(unlist(lapply(dispersal.into.hot.rates.BSM.CI,function(x)x[2])))),col=adjustcolor( "blue", alpha.f = 0.2))
    lapply(dispersal.into.realm.rates.BSM.CI.control,function(x)lines(x=c(1:length(dispersal.into.hot.rates.BSM.CI)),y=unlist(lapply(x,function(x)x[3])),col=adjustcolor( "orange", alpha.f = 0.9)))
    lapply(dispersal.into.hot.rates.BSM.CI.control,function(x)lines(x=c(1:length(dispersal.into.hot.rates.BSM.CI)),y=unlist(lapply(x,function(x)x[3])),col=adjustcolor( "green", alpha.f = 0.9)))
    axis(1,at=c(0:length(dispersal.into.hot.rates.BSM.CI)),labels=c((length(dispersal.into.hot.rates.BSM.CI)+2):2))
    legend('topleft',c('dispersal.into.hot','dispersal.into.realm','dispersal.into.hot.control','dispersal.into.realm.control'),col=c('red','blue','orange','green'),lty=1,cex=.7,bty='n')
  }else if(plots=='dispersal.all.events'){
    plot(c(1,1),xlim=c(0,length(dispersal.into.realm.time.bins.counts.rev.cumsum.CI)),ylim=c(0,max(c(unlist(dispersal.into.realm.time.bins.counts.rev.cumsum.CI),unlist(dispersal.into.hot.time.bins.counts.rev.cumsum.CI))[is.finite(c(unlist(dispersal.into.realm.time.bins.counts.rev.cumsum.CI),unlist(dispersal.into.hot.time.bins.counts.rev.cumsum.CI)))])),type='n',xaxt='n',xlab='age',ylab='cumulative # of events',main=paste(name,' dispersal events',sep=''))
    lines(x=c(1:length(dispersal.into.realm.time.bins.counts.rev.cumsum.CI)),y=unlist(lapply(dispersal.into.realm.time.bins.counts.rev.cumsum.CI,function(x)x[3])),col=adjustcolor( "red", alpha.f = 0.9))
    #polygon(x=c((1:length(dispersal.into.realm.time.bins.counts.rev.cumsum.CI)),(rev(1:length(dispersal.into.realm.time.bins.counts.rev.cumsum.CI)))),y=c((unlist(lapply(dispersal.into.realm.time.bins.counts.rev.cumsum.CI,function(x)x[1]))),rev(unlist(lapply(dispersal.into.realm.time.bins.counts.rev.cumsum.CI,function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
    lines(x=c(1:length(dispersal.into.hot.time.bins.counts.rev.cumsum.CI)),y=unlist(lapply(dispersal.into.hot.time.bins.counts.rev.cumsum.CI,function(x)x[3])),col=adjustcolor( "blue", alpha.f = 0.9))
    #polygon(x=c((1:length(dispersal.into.hot.time.bins.counts.rev.cumsum.CI)),(rev(1:length(dispersal.into.hot.time.bins.counts.rev.cumsum.CI)))),y=c((unlist(lapply(dispersal.into.hot.time.bins.counts.rev.cumsum.CI,function(x)x[1]))),rev(unlist(lapply(dispersal.into.hot.time.bins.counts.rev.cumsum.CI,function(x)x[2])))),col=adjustcolor( "blue", alpha.f = 0.2))
    lapply(dispersal.into.realm.time.bins.counts.rev.cumsum.CI.control,function(x)lines(x=c(1:length(dispersal.into.hot.time.bins.counts.rev.cumsum.CI)),y=unlist(lapply(x,function(x)x[3])),col=adjustcolor( "orange", alpha.f = 0.9)))
    lapply(dispersal.into.hot.time.bins.counts.rev.cumsum.CI.control,function(x)lines(x=c(1:length(dispersal.into.hot.time.bins.counts.rev.cumsum.CI)),y=unlist(lapply(x,function(x)x[3])),col=adjustcolor( "green", alpha.f = 0.9)))
    axis(1,at=c(0:length(BSMoutput.real.rates[['insitu.realm.rates.BSM.CI']])),labels=c((length(BSMoutput.real.rates[['insitu.realm.rates.BSM.CI']])+2):2))
    legend('topleft',c('dispersal.into.realm','dispersal.into.hot','dispersal.into.realm.control','dispersal.into.hot.control'),col=c('red','blue','orange','green'),lty=1,cex=.7,bty='n')
  }else if(plots=='dispersal.into.from.events'){
    plot(c(1,1),xlim=c(0,length(dispersal.into.realm.time.bins.counts.rev.cumsum.CI)),ylim=c(0,max(c(unlist(dispersal.into.realm.time.bins.counts.rev.cumsum.CI),unlist(dispersal.into.hot.time.bins.counts.rev.cumsum.CI),unlist(dispersal.from.hot.time.bins.counts.rev.cumsum.CI),unlist(dispersal.from.realm.time.bins.counts.rev.cumsum.CI)))),type='n',xaxt='n',xlab='age',ylab='cumulative # of events',main=paste(name,' dispersal events',sep=''))
    lines(x=c(1:length(dispersal.into.realm.time.bins.counts.rev.cumsum.CI)),y=unlist(lapply(dispersal.into.realm.time.bins.counts.rev.cumsum.CI,function(x)x[3])),col=adjustcolor( "green", alpha.f = 0.9))
    polygon(x=c((1:length(dispersal.into.realm.time.bins.counts.rev.cumsum.CI)),(rev(1:length(dispersal.into.realm.time.bins.counts.rev.cumsum.CI)))),y=c((unlist(lapply(dispersal.into.realm.time.bins.counts.rev.cumsum.CI,function(x)x[1]))),rev(unlist(lapply(dispersal.into.realm.time.bins.counts.rev.cumsum.CI,function(x)x[2])))),col=adjustcolor( "green", alpha.f = 0.2))
    lines(x=c(1:length(dispersal.from.realm.time.bins.counts.rev.cumsum.CI)),y=unlist(lapply(dispersal.from.realm.time.bins.counts.rev.cumsum.CI,function(x)x[3])),col=adjustcolor( "blue", alpha.f = 0.9))
    polygon(x=c((1:length(dispersal.from.realm.time.bins.counts.rev.cumsum.CI)),(rev(1:length(dispersal.from.realm.time.bins.counts.rev.cumsum.CI)))),y=c((unlist(lapply(dispersal.from.realm.time.bins.counts.rev.cumsum.CI,function(x)x[1]))),rev(unlist(lapply(dispersal.from.realm.time.bins.counts.rev.cumsum.CI,function(x)x[2])))),col=adjustcolor( "blue", alpha.f = 0.2))
    lapply(dispersal.into.realm.time.bins.counts.rev.cumsum.CI.control,function(x)lines(x=c(1:length(dispersal.into.hot.time.bins.counts.rev.cumsum.CI)),y=unlist(lapply(x,function(x)x[3])),col=adjustcolor( "green", alpha.f = 0.9)))
    lapply(dispersal.from.realm.time.bins.counts.rev.cumsum.CI.control,function(x)lines(x=c(1:length(dispersal.from.realm.time.bins.counts.rev.cumsum.CI)),y=unlist(lapply(x,function(x)x[3])),col=adjustcolor( "blue", alpha.f = 0.9)))
    axis(1,at=c(0:length(BSMoutput.real.rates[['insitu.realm.rates.BSM.CI']])),labels=c((length(BSMoutput.real.rates[['insitu.realm.rates.BSM.CI']])+2):2))
    legend('topleft',c('dispersal.from.realm','dispersal.into.realm'),col=c('blue','green'),lty=1,cex=.7,bty='n')
    
    plot(c(1,1),xlim=c(0,length(dispersal.into.realm.time.bins.counts.rev.cumsum.CI)),ylim=c(0,max(c(unlist(dispersal.into.realm.time.bins.counts.rev.cumsum.CI),unlist(dispersal.into.hot.time.bins.counts.rev.cumsum.CI),unlist(dispersal.from.hot.time.bins.counts.rev.cumsum.CI),unlist(dispersal.from.realm.time.bins.counts.rev.cumsum.CI)))),type='n',xaxt='n',xlab='age',ylab='cumulative # of events',main=paste(name,' dispersal events',sep=''))
    lines(x=c(1:length(dispersal.into.hot.time.bins.counts.rev.cumsum.CI)),y=unlist(lapply(dispersal.into.hot.time.bins.counts.rev.cumsum.CI,function(x)x[3])),col=adjustcolor( "orange", alpha.f = 0.9))
    polygon(x=c((1:length(dispersal.into.hot.time.bins.counts.rev.cumsum.CI)),(rev(1:length(dispersal.into.hot.time.bins.counts.rev.cumsum.CI)))),y=c((unlist(lapply(dispersal.into.hot.time.bins.counts.rev.cumsum.CI,function(x)x[1]))),rev(unlist(lapply(dispersal.into.hot.time.bins.counts.rev.cumsum.CI,function(x)x[2])))),col=adjustcolor( "orange", alpha.f = 0.2))
    lines(x=c(1:length(dispersal.from.hot.time.bins.counts.rev.cumsum.CI)),y=unlist(lapply(dispersal.from.hot.time.bins.counts.rev.cumsum.CI,function(x)x[3])),col=adjustcolor( "red", alpha.f = 0.9))
    polygon(x=c((1:length(dispersal.from.hot.time.bins.counts.rev.cumsum.CI)),(rev(1:length(dispersal.from.hot.time.bins.counts.rev.cumsum.CI)))),y=c((unlist(lapply(dispersal.from.hot.time.bins.counts.rev.cumsum.CI,function(x)x[1]))),rev(unlist(lapply(dispersal.from.hot.time.bins.counts.rev.cumsum.CI,function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
    lapply(dispersal.into.hot.time.bins.counts.rev.cumsum.CI.control,function(x)lines(x=c(1:length(dispersal.into.hot.time.bins.counts.rev.cumsum.CI)),y=unlist(lapply(x,function(x)x[3])),col=adjustcolor( "orange", alpha.f = 0.9)))
    lapply(dispersal.from.hot.time.bins.counts.rev.cumsum.CI.control,function(x)lines(x=c(1:length(dispersal.from.hot.time.bins.counts.rev.cumsum.CI)),y=unlist(lapply(x,function(x)x[3])),col=adjustcolor( "red", alpha.f = 0.9)))
    axis(1,at=c(0:length(BSMoutput.real.rates[['insitu.hot.rates.BSM.CI']])),labels=c((length(BSMoutput.real.rates[['insitu.hot.rates.BSM.CI']])+2):2))
    legend('topleft',c('dispersal.from.hot','dispersal.into.hot'),col=c('red','orange'),lty=1,cex=.7,bty='n')
    
  }
  
  
}

plot_rates_events_realvscontrol_difference_2mybin_polygon<-function(BSMoutput.real.rates,BSMoutput.real.events,BSMoutput.control.rates,BSMoutput.control.events,age,name,plots){
  #correct age
  age<-(length(BSMoutput.real.rates[['insitu.realm.rates.BSM.CI']])*2-age)/2
  #trim lists to age period
  BSMoutput.real.rates<-lapply(BSMoutput.real.rates,function(x)x[age:length(x)])
  BSMoutput.real.events<-lapply(BSMoutput.real.events,function(x)x[age:length(x)])
  BSMoutput.control.rates<-lapply(BSMoutput.control.rates,function(x) lapply(x,function(x) x[age:length(x)]))
  BSMoutput.control.events<-lapply(BSMoutput.control.events,function(x) lapply(x,function(x) x[age:length(x)]))
  #delete Inf in rates and events
  BSMoutput.real.rates<-lapply(BSMoutput.real.rates,function(x) lapply(x,function(x) {finite.hit<-which(!is.finite(x));if(length(finite.hit)==0){return(x)}else{if(finite.hit==2){x<-c(x[1],x[3],x[3])};return(x)}}))
  BSMoutput.real.events<-lapply(BSMoutput.real.events,function(x) lapply(x,function(x) {finite.hit<-which(!is.finite(x));if(length(finite.hit)==0){return(x)}else{if(finite.hit==2){x<-c(x[1],x[3],x[3])};return(x)}}))
  BSMoutput.control.rates<-lapply(BSMoutput.control.rates,function(x) lapply(x,function(x) lapply(x,function(x) {finite.hit<-which(!is.finite(x));if(length(finite.hit)==0){return(x)}else if (length(finite.hit)==1){if(finite.hit==2){x<-c(x[1],x[3],x[3]);return(x)}}else if(length(finite.hit)==2){x<-c(x[1],x[1],x[1])}})))                                 
  BSMoutput.control.events<-lapply(BSMoutput.control.events,function(x) lapply(x,function(x) lapply(x,function(x) {finite.hit<-which(!is.finite(x));if(length(finite.hit)==0){return(x)}else if (length(finite.hit)==1){if(finite.hit==2){x<-c(x[1],x[3],x[3]);return(x)}}else if(length(finite.hit)==2){x<-c(x[1],x[1],x[1])}})))                                 
  #remove the last bin - two million year time bin
  BSMoutput.real.rates<-lapply(BSMoutput.real.rates,function(x) x<-x[-c(length(x))])
  BSMoutput.real.events<-lapply(BSMoutput.real.events,function(x) x<-x[-c((length(x)-1),length(x))])
  BSMoutput.control.rates<-lapply(BSMoutput.control.rates,function(x) lapply(x,function(x) x<-x[-c(length(x))]))
  BSMoutput.control.events<-lapply(BSMoutput.control.events,function(x) lapply(x,function(x) x<-x[-c(length(x))]))
  #assign variables for easier handling
  hot.into.realm.rates.BSM.CI<-BSMoutput.real.rates[['hot.into.realm.rates.BSM.CI']]
  realm.into.hot.rates.BSM.CI<-BSMoutput.real.rates[['realm.into.hot.rates.BSM.CI']]
  insitu.realm.rates.BAMM.CI<-BSMoutput.real.rates[['insitu.realm.rates.BAMM.CI']]
  insitu.hot.rates.BAMM.CI<-BSMoutput.real.rates[['insitu.hot.rates.BAMM.CI']]
  insitu.realm.rates.BSM.CI<-BSMoutput.real.rates[['insitu.realm.rates.BSM.CI']]
  insitu.hot.rates.BSM.CI<-BSMoutput.real.rates[['insitu.hot.rates.BSM.CI']]
  dispersal.into.realm.rates.BSM.CI<-BSMoutput.real.rates[['dispersal.into.realm.rates.BSM.CI']]
  dispersal.into.hot.rates.BSM.CI<-BSMoutput.real.rates[['dispersal.into.hot.rates.BSM.CI']]
  hot.into.realm.time.bins.counts.rev.cumsum.CI<-BSMoutput.real.events[['hot.into.realm.time.bins.counts.rev.cumsum.CI']]
  realm.into.hot.time.bins.counts.rev.cumsum.CI<-BSMoutput.real.events[['realm.into.hot.time.bins.counts.rev.cumsum.CI']]
  insitu.realm.time.bins.counts.rev.cumsum.CI<-BSMoutput.real.events[['insitu.realm.time.bins.counts.rev.cumsum.CI']]
  insitu.hot.time.bins.counts.rev.cumsum.CI<-BSMoutput.real.events[['insitu.hot.time.bins.counts.rev.cumsum.CI']]
  dispersal.into.realm.time.bins.counts.rev.cumsum.CI<-BSMoutput.real.events[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']]
  dispersal.into.hot.time.bins.counts.rev.cumsum.CI<-BSMoutput.real.events[['dispersal.into.hot.time.bins.counts.rev.cumsum.CI']]
  dispersal.from.realm.time.bins.counts.rev.cumsum.CI<-BSMoutput.real.events[['dispersal.from.realm.time.bins.counts.rev.cumsum.CI']]
  dispersal.from.hot.time.bins.counts.rev.cumsum.CI<-BSMoutput.real.events[['dispersal.from.hot.time.bins.counts.rev.cumsum.CI']]
  #controls
  hot.into.realm.rates.BSM.CI.control<-lapply(BSMoutput.control.rates,function(x)x[['hot.into.realm.rates.BSM.CI']])
  realm.into.hot.rates.BSM.CI.control<-lapply(BSMoutput.control.rates,function(x) x[['realm.into.hot.rates.BSM.CI']])
  insitu.realm.rates.BAMM.CI.control<-lapply(BSMoutput.control.rates,function(x) x[['insitu.realm.rates.BAMM.CI']])
  insitu.hot.rates.BAMM.CI.control<-lapply(BSMoutput.control.rates,function(x) x[['insitu.hot.rates.BAMM.CI']])
  insitu.realm.rates.BSM.CI.control<-lapply(BSMoutput.control.rates,function(x) x[['insitu.realm.rates.BSM.CI']])
  insitu.hot.rates.BSM.CI.control<-lapply(BSMoutput.control.rates,function(x) x[['insitu.hot.rates.BSM.CI']])
  dispersal.into.realm.rates.BSM.CI.control<-lapply(BSMoutput.control.rates,function(x) x[['dispersal.into.realm.rates.BSM.CI']])
  dispersal.into.hot.rates.BSM.CI.control<-lapply(BSMoutput.control.rates,function(x) x[['dispersal.into.hot.rates.BSM.CI']])
  hot.into.realm.time.bins.counts.rev.cumsum.CI.control<-lapply(BSMoutput.control.events,function(x) x[['hot.into.realm.time.bins.counts.rev.cumsum.CI']])
  realm.into.hot.time.bins.counts.rev.cumsum.CI.control<-lapply(BSMoutput.control.events,function(x) x[['realm.into.hot.time.bins.counts.rev.cumsum.CI']])
  insitu.realm.time.bins.counts.rev.cumsum.CI.control<-lapply(BSMoutput.control.events,function(x) x[['insitu.realm.time.bins.counts.rev.cumsum.CI']])
  insitu.hot.time.bins.counts.rev.cumsum.CI.control<-lapply(BSMoutput.control.events,function(x) x[['insitu.hot.time.bins.counts.rev.cumsum.CI']])
  dispersal.into.realm.time.bins.counts.rev.cumsum.CI.control<-lapply(BSMoutput.control.events,function(x) x[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']])
  dispersal.into.hot.time.bins.counts.rev.cumsum.CI.control<-lapply(BSMoutput.control.events,function(x) x[['dispersal.into.hot.time.bins.counts.rev.cumsum.CI']])
  dispersal.from.realm.time.bins.counts.rev.cumsum.CI.control<-lapply(BSMoutput.control.events,function(x) x[['dispersal.from.realm.time.bins.counts.rev.cumsum.CI']])
  dispersal.from.hot.time.bins.counts.rev.cumsum.CI.control<-lapply(BSMoutput.control.events,function(x) x[['dispersal.from.hot.time.bins.counts.rev.cumsum.CI']])
  #for each time bin, substract the median of the real rate vs all the median of the 50 control polygons
  if(plots=='dispersal.rates'){
    hot.into.realm.diff.real.control<-lapply(c(1:length(hot.into.realm.rates.BSM.CI)),function(x){real<-hot.into.realm.rates.BSM.CI[[x]][3];control<-unlist(lapply(hot.into.realm.rates.BSM.CI.control,function(y)y[[x]][3]));return(real-control)})
    realm.into.hot.diff.real.control<-lapply(c(1:length(realm.into.hot.rates.BSM.CI)),function(x){real<-realm.into.hot.rates.BSM.CI[[x]][3];control<-unlist(lapply(realm.into.hot.rates.BSM.CI.control,function(y)y[[x]][3]));return(real-control)})
    hot.into.realm.diff.real.control<-lapply(hot.into.realm.diff.real.control,function(x)quantile(x,c(0.05,0.95,0.5)))
    realm.into.hot.diff.real.control<-lapply(realm.into.hot.diff.real.control,function(x)quantile(x,c(0.05,0.95,0.5)))
    #plot both rates into one plot
    plot(c(1,1),xlim=c(0,length(realm.into.hot.diff.real.control)+1),ylim=c(min(unlist(c(hot.into.realm.diff.real.control,realm.into.hot.diff.real.control))),max(unlist(c(hot.into.realm.diff.real.control,realm.into.hot.diff.real.control)))),type='n',xlab='age',xaxt='n',ylab='(real-control rate)',main=paste(name,' dispersal rates',sep=''))
    abline(h=0,lty=2)
    polygon(x=c(c(1:length(hot.into.realm.diff.real.control)),rev(c(1:length(hot.into.realm.diff.real.control)))),y=c(unlist(lapply(hot.into.realm.diff.real.control,function(x)x[1])),rev(unlist(lapply(hot.into.realm.diff.real.control,function(x)x[2])))),col=adjustcolor('red',alpha=0.25))
    lines(x=c(1:length(hot.into.realm.diff.real.control)),y=unlist(lapply(hot.into.realm.diff.real.control,function(x)x[3])),col='red')
    polygon(x=c(c(1:length(realm.into.hot.diff.real.control)),rev(c(1:length(realm.into.hot.diff.real.control)))),y=c(unlist(lapply(realm.into.hot.diff.real.control,function(x)x[1])),rev(unlist(lapply(realm.into.hot.diff.real.control,function(x)x[2])))),col=adjustcolor('blue',alpha=0.25))
    lines(x=c(1:length(realm.into.hot.diff.real.control)),y=unlist(lapply(realm.into.hot.diff.real.control,function(x)x[3])),col='blue')
    axis(1,at=c(1:length(realm.into.hot.diff.real.control)),labels=c(seq(from=length(realm.into.hot.diff.real.control)*2,to=2,by=-2)))
    legend('bottomleft',c('hot.into.realm','realm.into.hot'),col=c('red','blue'),lty=1,cex=.7,bty='n')
  }else if(plots=='cladogenesis.BAMM.rates'){
    insitu.hot.diff.real.control<-lapply(c(1:length(insitu.hot.rates.BAMM.CI)),function(x){real<-insitu.hot.rates.BAMM.CI[[x]][3];control<-unlist(lapply(insitu.hot.rates.BAMM.CI.control,function(y)y[[x]][3]));return(real-control)})
    insitu.realm.diff.real.control<-lapply(c(1:length(insitu.realm.rates.BAMM.CI)),function(x){real<-insitu.realm.rates.BAMM.CI[[x]][3];control<-unlist(lapply(insitu.realm.rates.BAMM.CI.control,function(y)y[[x]][3]));return(real-control)})
    insitu.hot.diff.real.control<-lapply(insitu.hot.diff.real.control,function(x)quantile(x,c(0.05,0.95,0.5)))
    insitu.realm.diff.real.control<-lapply(insitu.realm.diff.real.control,function(x)quantile(x,c(0.05,0.95,0.5)))
    #plot both rates into one plot
    plot(c(1,1),xlim=c(0,length(insitu.realm.diff.real.control)+1),ylim=c(min(unlist(c(insitu.hot.diff.real.control,insitu.realm.diff.real.control))),max(unlist(c(insitu.hot.diff.real.control,insitu.realm.diff.real.control)))),type='n',xlab='age',xaxt='n',ylab='(real-control rate)',main=paste(name,' cladogenesis BAMM rates',sep=''))
    abline(h=0,lty=2)
    polygon(x=c(c(1:length(insitu.hot.diff.real.control)),rev(c(1:length(insitu.hot.diff.real.control)))),y=c(unlist(lapply(insitu.hot.diff.real.control,function(x)x[1])),rev(unlist(lapply(insitu.hot.diff.real.control,function(x)x[2])))),col=adjustcolor('red',alpha=0.25))
    lines(x=c(1:length(insitu.hot.diff.real.control)),y=unlist(lapply(insitu.hot.diff.real.control,function(x)x[3])),col='red')
    polygon(x=c(c(1:length(insitu.realm.diff.real.control)),rev(c(1:length(insitu.realm.diff.real.control)))),y=c(unlist(lapply(insitu.realm.diff.real.control,function(x)x[1])),rev(unlist(lapply(insitu.realm.diff.real.control,function(x)x[2])))),col=adjustcolor('blue',alpha=0.25))
    lines(x=c(1:length(insitu.realm.diff.real.control)),y=unlist(lapply(insitu.realm.diff.real.control,function(x)x[3])),col='blue')
    axis(1,at=c(1:length(insitu.realm.diff.real.control)),labels=c(seq(from=length(insitu.realm.diff.real.control)*2,to=2,by=-2)))
    legend('bottomleft',c('insitu.hot','insitu.realm'),col=c('red','blue'),lty=1,cex=.7,bty='n')
  }else if(plots=='cladogenesis.BSM.rates'){
    insitu.hot.diff.real.control<-lapply(c(1:length(insitu.hot.rates.BSM.CI)),function(x){real<-insitu.hot.rates.BSM.CI[[x]][3];control<-unlist(lapply(insitu.hot.rates.BSM.CI.control,function(y)y[[x]][3]));return(real-control)})
    insitu.realm.diff.real.control<-lapply(c(1:length(insitu.realm.rates.BSM.CI)),function(x){real<-insitu.realm.rates.BSM.CI[[x]][3];control<-unlist(lapply(insitu.realm.rates.BSM.CI.control,function(y)y[[x]][3]));return(real-control)})
    insitu.hot.diff.real.control<-lapply(insitu.hot.diff.real.control,function(x)quantile(x,c(0.05,0.95,0.5)))
    insitu.realm.diff.real.control<-lapply(insitu.realm.diff.real.control,function(x)quantile(x,c(0.05,0.95,0.5)))
    #plot both rates into one plot
    plot(c(1,1),xlim=c(0,length(insitu.hot.diff.real.control)+1),ylim=c(min(unlist(insitu.hot.diff.real.control)),max(unlist(insitu.hot.diff.real.control))),type='n',xlab='age',xaxt='n',ylab='diff',main=paste(name,' cladogenesis BSM rates',sep=''))
    abline(h=0,lty=2)
    polygon(x=c(c(1:length(insitu.hot.diff.real.control)),rev(c(1:length(insitu.hot.diff.real.control)))),y=c(unlist(lapply(insitu.hot.diff.real.control,function(x)x[1])),rev(unlist(lapply(insitu.hot.diff.real.control,function(x)x[2])))),col=adjustcolor('red',alpha=0.25))
    lines(x=c(1:length(insitu.hot.diff.real.control)),y=unlist(lapply(insitu.hot.diff.real.control,function(x)x[3])),col='red')
    polygon(x=c(c(1:length(insitu.realm.diff.real.control)),rev(c(1:length(insitu.realm.diff.real.control)))),y=c(unlist(lapply(insitu.realm.diff.real.control,function(x)x[1])),rev(unlist(lapply(insitu.realm.diff.real.control,function(x)x[2])))),col=adjustcolor('blue',alpha=0.25))
    lines(x=c(1:length(insitu.realm.diff.real.control)),y=unlist(lapply(insitu.realm.diff.real.control,function(x)x[3])),col='blue')
    axis(1,at=c(1:length(insitu.realm.diff.real.control)),labels=c(seq(from=length(insitu.realm.diff.real.control)*2,to=2,by=-2)))
    legend('bottomleft',c('insitu.hot','insitu.realm'),col=c('red','blue'),lty=1,cex=.7,bty='n')
  }else if(plots=='dispersal.events'){
    hot.into.realm.diff.real.control<-lapply(c(1:length(hot.into.realm.time.bins.counts.rev.cumsum.CI)),function(x){real<-hot.into.realm.time.bins.counts.rev.cumsum.CI[[x]][3];control<-unlist(lapply(hot.into.realm.time.bins.counts.rev.cumsum.CI.control,function(y)y[[x]][3]));return(real-control)})
    realm.into.hot.diff.real.control<-lapply(c(1:length(realm.into.hot.time.bins.counts.rev.cumsum.CI)),function(x){real<-realm.into.hot.time.bins.counts.rev.cumsum.CI[[x]][3];control<-unlist(lapply(realm.into.hot.time.bins.counts.rev.cumsum.CI.control,function(y)y[[x]][3]));return(real-control)})
    hot.into.realm.diff.real.control<-lapply(hot.into.realm.diff.real.control,function(x)quantile(x,c(0.05,0.95,0.5)))
    realm.into.hot.diff.real.control<-lapply(realm.into.hot.diff.real.control,function(x)quantile(x,c(0.05,0.95,0.5)))
    #plot both rates into one plot
    plot(c(1,1),xlim=c(0,length(realm.into.hot.diff.real.control)+1),ylim=c(min(unlist(c(hot.into.realm.diff.real.control,realm.into.hot.diff.real.control))),max(unlist(c(hot.into.realm.diff.real.control,realm.into.hot.diff.real.control)))),type='n',xlab='age',xaxt='n',ylab='(real-control rate)',main=paste(name,' dispersal events',sep=''))
    abline(h=0,lty=2)
    polygon(x=c(c(1:length(hot.into.realm.diff.real.control)),rev(c(1:length(hot.into.realm.diff.real.control)))),y=c(unlist(lapply(hot.into.realm.diff.real.control,function(x)x[1])),rev(unlist(lapply(hot.into.realm.diff.real.control,function(x)x[2])))),col=adjustcolor('red',alpha=0.25))
    lines(x=c(1:length(hot.into.realm.diff.real.control)),y=unlist(lapply(hot.into.realm.diff.real.control,function(x)x[3])),col='red')
    polygon(x=c(c(1:length(realm.into.hot.diff.real.control)),rev(c(1:length(realm.into.hot.diff.real.control)))),y=c(unlist(lapply(realm.into.hot.diff.real.control,function(x)x[1])),rev(unlist(lapply(realm.into.hot.diff.real.control,function(x)x[2])))),col=adjustcolor('blue',alpha=0.25))
    lines(x=c(1:length(realm.into.hot.diff.real.control)),y=unlist(lapply(realm.into.hot.diff.real.control,function(x)x[3])),col='blue')
    axis(1,at=c(1:length(realm.into.hot.diff.real.control)),labels=c(seq(from=length(realm.into.hot.diff.real.control)*2,to=2,by=-2)))
  }else if(plots=='cladogenesis.events'){
    insitu.hot.diff.real.control<-lapply(c(1:length(insitu.hot.time.bins.counts.rev.cumsum.CI)),function(x){real<-insitu.hot.time.bins.counts.rev.cumsum.CI[[x]][3];control<-unlist(lapply(insitu.hot.time.bins.counts.rev.cumsum.CI.control,function(y)y[[x]][3]));return(real-control)})
    insitu.realm.diff.real.control<-lapply(c(1:length(insitu.realm.time.bins.counts.rev.cumsum.CI)),function(x){real<-insitu.realm.time.bins.counts.rev.cumsum.CI[[x]][3];control<-unlist(lapply(insitu.realm.time.bins.counts.rev.cumsum.CI.control,function(y)y[[x]][3]));return(real-control)})
    insitu.hot.diff.real.control<-lapply(insitu.hot.diff.real.control,function(x)quantile(x,c(0.05,0.95,0.5)))
    insitu.realm.diff.real.control<-lapply(insitu.realm.diff.real.control,function(x)quantile(x,c(0.05,0.95,0.5)))
    #plot both rates into one plot
    plot(c(1,1),xlim=c(0,length(insitu.realm.diff.real.control)+1),ylim=c(min(unlist(c(insitu.hot.diff.real.control,insitu.realm.diff.real.control))),max(unlist(c(insitu.hot.diff.real.control,insitu.realm.diff.real.control)))),type='n',xlab='age',xaxt='n',ylab='(real-control events)',main=paste(name,' cladogenesis.events',sep=''))
    abline(h=0,lty=2)
    polygon(x=c(c(1:length(insitu.hot.diff.real.control)),rev(c(1:length(insitu.hot.diff.real.control)))),y=c(unlist(lapply(insitu.hot.diff.real.control,function(x)x[1])),rev(unlist(lapply(insitu.hot.diff.real.control,function(x)x[2])))),col=adjustcolor('red',alpha=0.25))
    lines(x=c(1:length(insitu.hot.diff.real.control)),y=unlist(lapply(insitu.hot.diff.real.control,function(x)x[3])),col='red')
    polygon(x=c(c(1:length(insitu.realm.diff.real.control)),rev(c(1:length(insitu.realm.diff.real.control)))),y=c(unlist(lapply(insitu.realm.diff.real.control,function(x)x[1])),rev(unlist(lapply(insitu.realm.diff.real.control,function(x)x[2])))),col=adjustcolor('blue',alpha=0.25))
    lines(x=c(1:length(insitu.realm.diff.real.control)),y=unlist(lapply(insitu.realm.diff.real.control,function(x)x[3])),col='blue')
    axis(1,at=c(1:length(insitu.realm.diff.real.control)),labels=c(seq(from=length(insitu.realm.diff.real.control)*2,to=2,by=-2)))
    legend('bottomleft',c('insitu.hot','insitu.realm'),col=c('red','blue'),lty=1,cex=.7,bty='n')
  }else if (plots=='dispersal.all.rates'){
    dispersal.into.realm.diff.real.control<-lapply(c(1:length(dispersal.into.realm.rates.BSM.CI)),function(x){real<-dispersal.into.realm.rates.BSM.CI[[x]][3];control<-unlist(lapply(dispersal.into.realm.rates.BSM.CI.control,function(y)y[[x]][3]));return(real-control)})
    dispersal.into.hot.diff.real.control<-lapply(c(1:length(dispersal.into.hot.rates.BSM.CI)),function(x){real<-dispersal.into.hot.rates.BSM.CI[[x]][3];control<-unlist(lapply(dispersal.into.hot.rates.BSM.CI.control,function(y)y[[x]][3]));return(real-control)})
    dispersal.into.realm.diff.real.control<-lapply(dispersal.into.realm.diff.real.control,function(x)quantile(x,c(0.05,0.95,0.5)))
    dispersal.into.hot.diff.real.control<-lapply(dispersal.into.hot.diff.real.control,function(x)quantile(x,c(0.05,0.95,0.5)))
    #plot both rates into one plot
    plot(c(1,1),xlim=c(0,length(dispersal.into.hot.diff.real.control)+1),ylim=c(min(unlist(c(dispersal.into.realm.diff.real.control,dispersal.into.hot.diff.real.control))),max(unlist(c(dispersal.into.realm.diff.real.control,dispersal.into.hot.diff.real.control)))),type='n',xlab='age',xaxt='n',ylab='(real-control rate)',main=paste(name,' dispersal.all.rates',sep=''))
    abline(h=0,lty=2)
    polygon(x=c(c(1:length(dispersal.into.realm.diff.real.control)),rev(c(1:length(dispersal.into.realm.diff.real.control)))),y=c(unlist(lapply(dispersal.into.realm.diff.real.control,function(x)x[1])),rev(unlist(lapply(dispersal.into.realm.diff.real.control,function(x)x[2])))),col=adjustcolor('red',alpha=0.25))
    lines(x=c(1:length(dispersal.into.realm.diff.real.control)),y=unlist(lapply(dispersal.into.realm.diff.real.control,function(x)x[3])),col='red')
    polygon(x=c(c(1:length(dispersal.into.hot.diff.real.control)),rev(c(1:length(dispersal.into.hot.diff.real.control)))),y=c(unlist(lapply(dispersal.into.hot.diff.real.control,function(x)x[1])),rev(unlist(lapply(dispersal.into.hot.diff.real.control,function(x)x[2])))),col=adjustcolor('blue',alpha=0.25))
    lines(x=c(1:length(dispersal.into.hot.diff.real.control)),y=unlist(lapply(dispersal.into.hot.diff.real.control,function(x)x[3])),col='blue')
    axis(1,at=c(1:length(dispersal.into.realm.diff.real.control)),labels=c(seq(from=length(dispersal.into.realm.diff.real.control)*2,to=2,by=-2)))
    legend('bottomleft',c('dispersal.into.realm','dispersal.into.hot'),col=c('red','blue'),lty=1,cex=.7,bty='n')
  }else if(plots=='dispersal.all.events'){
    dispersal.into.realm.diff.real.control<-lapply(c(1:length(dispersal.into.realm.time.bins.counts.rev.cumsum.CI)),function(x){real<-dispersal.into.realm.time.bins.counts.rev.cumsum.CI[[x]][3];control<-unlist(lapply(dispersal.into.realm.time.bins.counts.rev.cumsum.CI.control,function(y)y[[x]][3]));return(real-control)})
    dispersal.into.hot.diff.real.control<-lapply(c(1:length(dispersal.into.hot.time.bins.counts.rev.cumsum.CI)),function(x){real<-dispersal.into.hot.time.bins.counts.rev.cumsum.CI[[x]][3];control<-unlist(lapply(dispersal.into.hot.time.bins.counts.rev.cumsum.CI.control,function(y)y[[x]][3]));return(real-control)})
    dispersal.into.realm.diff.real.control<-lapply(dispersal.into.realm.diff.real.control,function(x)quantile(x,c(0.05,0.95,0.5)))
    dispersal.into.hot.diff.real.control<-lapply(dispersal.into.hot.diff.real.control,function(x)quantile(x,c(0.05,0.95,0.5)))
    #plot both rates into one plot
    plot(c(1,1),xlim=c(0,length(dispersal.into.hot.diff.real.control)+1),ylim=c(min(unlist(c(dispersal.into.realm.diff.real.control,dispersal.into.hot.diff.real.control))),max(unlist(c(dispersal.into.realm.diff.real.control,dispersal.into.hot.diff.real.control)))),type='n',xlab='age',xaxt='n',ylab='(real-control rate)',main=paste(name,' dispersal.all.events',sep=''))
    abline(h=0,lty=2)
    polygon(x=c(c(1:length(dispersal.into.realm.diff.real.control)),rev(c(1:length(dispersal.into.realm.diff.real.control)))),y=c(unlist(lapply(dispersal.into.realm.diff.real.control,function(x)x[1])),rev(unlist(lapply(dispersal.into.realm.diff.real.control,function(x)x[2])))),col=adjustcolor('red',alpha=0.25))
    lines(x=c(1:length(dispersal.into.realm.diff.real.control)),y=unlist(lapply(dispersal.into.realm.diff.real.control,function(x)x[3])),col='red')
    polygon(x=c(c(1:length(dispersal.into.hot.diff.real.control)),rev(c(1:length(dispersal.into.hot.diff.real.control)))),y=c(unlist(lapply(dispersal.into.hot.diff.real.control,function(x)x[1])),rev(unlist(lapply(dispersal.into.hot.diff.real.control,function(x)x[2])))),col=adjustcolor('blue',alpha=0.25))
    lines(x=c(1:length(dispersal.into.hot.diff.real.control)),y=unlist(lapply(dispersal.into.hot.diff.real.control,function(x)x[3])),col='blue')
    axis(1,at=c(1:length(dispersal.into.realm.diff.real.control)),labels=c(seq(from=length(dispersal.into.realm.diff.real.control)*2,to=2,by=-2)))
    legend('bottomleft',c('dispersal.into.realm','dispersal.into.hot'),col=c('red','blue'),lty=1,cex=.7,bty='n')
  }
  
}
plot_rates_events_realvscontrolOLD<-function(BSMoutput.real.rates,BSMoutput.real.events,BSMoutput.control.rates,BSMoutput.control.events,age,name){
  #correct age
  age<-length(BSMoutput.real.rates[['insitu.realm.rates.BSM.CI']])-age
  #remove last two million year time bins
  BSMoutput.real.rates<-lapply(BSMoutput.real.rates,function(x) x<-x[-c((length(x)-1),length(x))])
  BSMoutput.real.events<-lapply(BSMoutput.real.events,function(x) x<-x[-c((length(x)-2),(length(x)-1),length(x))])
  BSMoutput.control.rates<-lapply(BSMoutput.control.rates,function(x) lapply(x,function(x) x<-x[-c((length(x)-1),length(x))]))
  BSMoutput.control.events<-lapply(BSMoutput.control.events,function(x) lapply(x,function(x) x<-x[-c((length(x)-2),(length(x)-1),length(x))]))
  #plot BSM in situ rates realm
  plot(c(1,1),xlim=c(age,length(BSMoutput.real.rates[['insitu.realm.rates.BSM.CI']])),ylim=c(0,0.4),type='n',xaxt='n',xlab='age',ylab='speciation rates in situ',main=paste(name,' BSM in situ rates',sep=''))
  #polygon(x=c(c(age:length(BSMoutput.real.rates[['insitu.realm.rates.BSM.CI']])),c(length(BSMoutput.real.rates[['insitu.realm.rates.BSM.CI']]):age)),y=c(unlist(lapply(BSMoutput.real.rates[['insitu.realm.rates.BSM.CI']][c(age:length(BSMoutput.real.rates[['insitu.realm.rates.BSM.CI']]))],function(x)x[1])),rev(unlist(lapply(BSMoutput.real.rates[['insitu.realm.rates.BSM.CI']][c(age:length(BSMoutput.real.rates[['insitu.realm.rates.BSM.CI']]))],function(x)x[2])))),col=adjustcolor( "blue", alpha.f = 0.2))
  lines(x=c(age:length(BSMoutput.real.rates[['insitu.realm.rates.BSM.CI']])),y=unlist(lapply(BSMoutput.real.rates[['insitu.realm.rates.BSM.CI']][c(age:length(BSMoutput.real.rates[['insitu.realm.rates.BSM.CI']]))],function(x)x[3])),col=adjustcolor( "blue", alpha.f = 0.9))
  #lapply(BSMoutput.control.rates,function(x)polygon(x=c(c(age:length(x[['insitu.realm.rates.BSM.CI']])),c(length(x[['insitu.realm.rates.BSM.CI']]):age)),y=c(unlist(lapply(x[['insitu.realm.rates.BSM.CI']][c(age:length(x[['insitu.realm.rates.BSM.CI']]))],function(x)x[1])),rev(unlist(lapply(x[['insitu.realm.rates.BSM.CI']][c(age:length(x[['insitu.realm.rates.BSM.CI']]))],function(x)x[2])))),col=adjustcolor( "green", alpha.f = 0.2)))
  lapply(BSMoutput.control.rates,function(x)lines(x=c(age:length(x[['insitu.realm.rates.BSM.CI']])),y=unlist(lapply(x[['insitu.realm.rates.BSM.CI']][c(age:length(x[['insitu.realm.rates.BSM.CI']]))],function(x)x[3])),col=adjustcolor( "green", alpha.f = 0.9)))
  axis(1,at=c(age:length(BSMoutput.real.rates[['insitu.realm.rates.BSM.CI']])),labels=c((length(BSMoutput.real.rates[['insitu.realm.rates.BSM.CI']])-age+2):2))
  legend('topright',c('insitu.nonhot','insitu.nonhot.control'),col=c('blue','green'),lty=1,cex=.7,bty='n')
  #plot BSM in situ rates hot
  plot(c(1,1),xlim=c(age,length(BSMoutput.real.rates[['insitu.realm.rates.BSM.CI']])),ylim=c(0,0.4),type='n',xaxt='n',xlab='age',ylab='speciation rates in situ',main=paste(name,' BSM in situ rates',sep=''))
  #polygon(x=c(c(age:length(BSMoutput.real.rates[['insitu.hot.rates.BSM.CI']])),c(length(BSMoutput.real.rates[['insitu.hot.rates.BSM.CI']]):age)),y=c(unlist(lapply(BSMoutput.real.rates[['insitu.hot.rates.BSM.CI']][c(age:length(BSMoutput.real.rates[['insitu.hot.rates.BSM.CI']]))],function(x)x[1])),rev(unlist(lapply(BSMoutput.real.rates[['insitu.hot.rates.BSM.CI']][c(age:length(BSMoutput.real.rates[['insitu.hot.rates.BSM.CI']]))],function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
  lines(x=c(age:length(BSMoutput.real.rates[['insitu.hot.rates.BSM.CI']])),y=unlist(lapply(BSMoutput.real.rates[['insitu.hot.rates.BSM.CI']][c(age:length(BSMoutput.real.rates[['insitu.hot.rates.BSM.CI']]))],function(x)x[3])),col=adjustcolor( "red", alpha.f = 0.9))
  #lapply(BSMoutput.control.rates,function(x)polygon(x=c(c(age:length(x[['insitu.hot.rates.BSM.CI']])),c(length(x[['insitu.hot.rates.BSM.CI']]):age)),y=c(unlist(lapply(x[['insitu.hot.rates.BSM.CI']][c(age:length(x[['insitu.hot.rates.BSM.CI']]))],function(x)x[1])),rev(unlist(lapply(x[['insitu.hot.rates.BSM.CI']][c(age:length(x[['insitu.hot.rates.BSM.CI']]))],function(x)x[2])))),col=adjustcolor( "orange", alpha.f = 0.2)))
  lapply(BSMoutput.control.rates,function(x)lines(x=c(age:length(x[['insitu.hot.rates.BSM.CI']])),y=unlist(lapply(x[['insitu.hot.rates.BSM.CI']][c(age:length(x[['insitu.hot.rates.BSM.CI']]))],function(x)x[3])),col=adjustcolor( "orange", alpha.f = 0.9)))
  axis(1,at=c(age:length(BSMoutput.real.rates[['insitu.realm.rates.BSM.CI']])),labels=c((length(BSMoutput.real.rates[['insitu.realm.rates.BSM.CI']])-age+2):2))
  legend('topright',c('insitu.hot','insitu.hot.control'),col=c('red','orange'),lty=1,cex=.7,bty='n')
  #legend('topright',c('insitu.hot','insitu.nonhot','insitu.hot.control','insitu.nonhot.control'),col=c('red','blue','orange','green'),lty=1,cex=.7,bty='n')
  
  #plot BSM dispersal rates hot into realm
  plot(c(1,1),xlim=c(age,length(BSMoutput.real.rates[['hot.into.realm.rates.BSM.CI']])),ylim=c(0,0.2),type='n',xaxt='n',xlab='age',ylab='hot-non hot dispersal rates',main=paste(name,' BSM dispersal rates',sep=''))
  #polygon(x=c(c(age:length(BSMoutput.real.rates[['hot.into.realm.rates.BSM.CI']])),c(length(BSMoutput.real.rates[['hot.into.realm.rates.BSM.CI']]):age)),y=c(unlist(lapply(BSMoutput.real.rates[['hot.into.realm.rates.BSM.CI']][c(age:length(BSMoutput.real.rates[['hot.into.realm.rates.BSM.CI']]))],function(x)x[1])),rev(unlist(lapply(BSMoutput.real.rates[['hot.into.realm.rates.BSM.CI']][c(age:length(BSMoutput.real.rates[['hot.into.realm.rates.BSM.CI']]))],function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
  lines(x=c(age:length(BSMoutput.real.rates[['hot.into.realm.rates.BSM.CI']])),y=unlist(lapply(BSMoutput.real.rates[['hot.into.realm.rates.BSM.CI']][c(age:length(BSMoutput.real.rates[['hot.into.realm.rates.BSM.CI']]))],function(x)x[3])),col=adjustcolor( "red", alpha.f = 0.9))
  #lapply(BSMoutput.control.rates,function(x)polygon(x=c(c(age:length(x[['hot.into.realm.rates.BSM.CI']])),c(length(x[['hot.into.realm.rates.BSM.CI']]):age)),y=c(unlist(lapply(x[['hot.into.realm.rates.BSM.CI']][c(age:length(x[['hot.into.realm.rates.BSM.CI']]))],function(x)x[1])),rev(unlist(lapply(x[['hot.into.realm.rates.BSM.CI']][c(age:length(x[['hot.into.realm.rates.BSM.CI']]))],function(x)x[2])))),col=adjustcolor( "orange", alpha.f = 0.2)))
  lapply(BSMoutput.control.rates,function(x)lines(x=c(age:length(x[['hot.into.realm.rates.BSM.CI']])),y=unlist(lapply(x[['hot.into.realm.rates.BSM.CI']][c(age:length(x[['hot.into.realm.rates.BSM.CI']]))],function(x)x[3])),col=adjustcolor( "orange", alpha.f = 0.9)))
  axis(1,at=c(age:length(BSMoutput.real.rates[['hot.into.realm.rates.BSM.CI']])),labels=c((length(BSMoutput.real.rates[['hot.into.realm.rates.BSM.CI']])-age+2):2))
  legend('topright',c('hot.into.realm','hot.into.realm.control'),col=c('red','orange'),lty=1,cex=.7,bty='n')
  #plot BSM dispersal rates realm into hot
  plot(c(1,1),xlim=c(age,length(BSMoutput.real.rates[['hot.into.realm.rates.BSM.CI']])),ylim=c(0,0.2),type='n',xaxt='n',xlab='age',ylab='hot-non hot dispersal rates',main=paste(name,' BSM dispersal rates',sep=''))
  #polygon(x=c(c(age:length(BSMoutput.real.rates[['realm.into.hot.rates.BSM.CI']])),c(length(BSMoutput.real.rates[['realm.into.hot.rates.BSM.CI']]):age)),y=c(unlist(lapply(BSMoutput.real.rates[['realm.into.hot.rates.BSM.CI']][c(age:length(BSMoutput.real.rates[['realm.into.hot.rates.BSM.CI']]))],function(x)x[1])),rev(unlist(lapply(BSMoutput.real.rates[['realm.into.hot.rates.BSM.CI']][c(age:length(BSMoutput.real.rates[['realm.into.hot.rates.BSM.CI']]))],function(x)x[2])))),col=adjustcolor( "blue", alpha.f = 0.2))
  lines(x=c(age:length(BSMoutput.real.rates[['realm.into.hot.rates.BSM.CI']])),y=unlist(lapply(BSMoutput.real.rates[['realm.into.hot.rates.BSM.CI']][c(age:length(BSMoutput.real.rates[['realm.into.hot.rates.BSM.CI']]))],function(x)x[3])),col=adjustcolor( "blue", alpha.f = 0.9))
  #lapply(BSMoutput.control.rates,function(x)polygon(x=c(c(age:length(x[['realm.into.hot.rates.BSM.CI']])),c(length(x[['realm.into.hot.rates.BSM.CI']]):age)),y=c(unlist(lapply(x[['realm.into.hot.rates.BSM.CI']][c(age:length(x[['realm.into.hot.rates.BSM.CI']]))],function(x)x[1])),rev(unlist(lapply(x[['realm.into.hot.rates.BSM.CI']][c(age:length(x[['realm.into.hot.rates.BSM.CI']]))],function(x)x[2])))),col=adjustcolor( "green", alpha.f = 0.2)))
  lapply(BSMoutput.control.rates,function(x)lines(x=c(age:length(x[['realm.into.hot.rates.BSM.CI']])),y=unlist(lapply(x[['realm.into.hot.rates.BSM.CI']][c(age:length(x[['realm.into.hot.rates.BSM.CI']]))],function(x)x[3])),col=adjustcolor( "green", alpha.f = 0.9)))
  axis(1,at=c(age:length(BSMoutput.real.rates[['insitu.realm.rates.BSM.CI']])),labels=c((length(BSMoutput.real.rates[['insitu.realm.rates.BSM.CI']])-age+2):2))
  legend('topright',c('realm.into.hot','realm.into.hot.control'),col=c('blue','green'),lty=1,cex=.7,bty='n')
  #legend('topright',c('hot.into.realm','realm.into.hot','hot.into.realm.control','realm.into.hot.control'),col=c('red','blue','orange','purple'),lty=1,cex=.7,bty='n')
  
  #plot BSM general dispersal rates into realm
  plot(c(1,1),xlim=c(age,length(BSMoutput.real.rates[['dispersal.into.realm.rates.BSM.CI']])),ylim=c(0,0.3),type='n',xaxt='n',xlab='age',ylab='general dispersal rates',main=paste(name,' BSM dispersal rates',sep=''))
  #polygon(x=c(c(age:length(BSMoutput.real.rates[['dispersal.into.realm.rates.BSM.CI']])),c(length(BSMoutput.real.rates[['dispersal.into.realm.rates.BSM.CI']]):age)),y=c(unlist(lapply(BSMoutput.real.rates[['dispersal.into.realm.rates.BSM.CI']][c(age:length(BSMoutput.real.rates[['dispersal.into.realm.rates.BSM.CI']]))],function(x)x[1])),rev(unlist(lapply(BSMoutput.real.rates[['dispersal.into.realm.rates.BSM.CI']][c(age:length(BSMoutput.real.rates[['dispersal.into.realm.rates.BSM.CI']]))],function(x)x[2])))),col=adjustcolor( "blue", alpha.f = 0.2))
  lines(x=c(age:length(BSMoutput.real.rates[['dispersal.into.realm.rates.BSM.CI']])),y=unlist(lapply(BSMoutput.real.rates[['dispersal.into.realm.rates.BSM.CI']][c(age:length(BSMoutput.real.rates[['dispersal.into.realm.rates.BSM.CI']]))],function(x)x[3])),col=adjustcolor( "blue", alpha.f = 0.9))
  #lapply(BSMoutput.control.rates,function(x)polygon(x=c(c(age:length(x[['dispersal.into.realm.rates.BSM.CI']])),c(length(x[['dispersal.into.realm.rates.BSM.CI']]):age)),y=c(unlist(lapply(x[['dispersal.into.realm.rates.BSM.CI']][c(age:length(x[['dispersal.into.realm.rates.BSM.CI']]))],function(x)x[1])),rev(unlist(lapply(x[['dispersal.into.realm.rates.BSM.CI']][c(age:length(x[['dispersal.into.realm.rates.BSM.CI']]))],function(x)x[2])))),col=adjustcolor( "purple", alpha.f = 0.2)))
  lapply(BSMoutput.control.rates,function(x)lines(x=c(age:length(x[['dispersal.into.realm.rates.BSM.CI']])),y=unlist(lapply(x[['dispersal.into.realm.rates.BSM.CI']][c(age:length(x[['dispersal.into.realm.rates.BSM.CI']]))],function(x)x[3])),col=adjustcolor( "green", alpha.f = 0.9)))
  axis(1,at=c(age:length(BSMoutput.real.rates[['dispersal.into.realm.rates.BSM.CI']])),labels=c((length(BSMoutput.real.rates[['dispersal.into.realm.rates.BSM.CI']])-age+2):2))
  legend('topright',c('dispersal.into.realm','dispersal.into.realm.control'),col=c('blue','green'),lty=1,cex=.7,bty='n')
  #plot BSM general dispersal rates into hot
  plot(c(1,1),xlim=c(age,length(BSMoutput.real.rates[['dispersal.into.hot.rates.BSM.CI']])),ylim=c(0,0.3),type='n',xaxt='n',xlab='age',ylab='general dispersal rates',main=paste(name,' BSM dispersal rates',sep=''))
  #polygon(x=c(c(age:length(BSMoutput.real.rates[['dispersal.into.hot.rates.BSM.CI']])),c(length(BSMoutput.real.rates[['dispersal.into.hot.rates.BSM.CI']]):age)),y=c(unlist(lapply(BSMoutput.real.rates[['dispersal.into.hot.rates.BSM.CI']][c(age:length(BSMoutput.real.rates[['dispersal.into.hot.rates.BSM.CI']]))],function(x)x[1])),rev(unlist(lapply(BSMoutput.real.rates[['dispersal.into.hot.rates.BSM.CI']][c(age:length(BSMoutput.real.rates[['dispersal.into.hot.rates.BSM.CI']]))],function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
  lines(x=c(age:length(BSMoutput.real.rates[['dispersal.into.hot.rates.BSM.CI']])),y=unlist(lapply(BSMoutput.real.rates[['dispersal.into.hot.rates.BSM.CI']][c(age:length(BSMoutput.real.rates[['dispersal.into.hot.rates.BSM.CI']]))],function(x)x[3])),col=adjustcolor( "red", alpha.f = 0.9))
  #lapply(BSMoutput.control.rates,function(x)polygon(x=c(c(age:length(x[['dispersal.into.hot.rates.BSM.CI']])),c(length(x[['dispersal.into.hot.rates.BSM.CI']]):age)),y=c(unlist(lapply(x[['dispersal.into.hot.rates.BSM.CI']][c(age:length(x[['dispersal.into.hot.rates.BSM.CI']]))],function(x)x[1])),rev(unlist(lapply(x[['dispersal.into.hot.rates.BSM.CI']][c(age:length(x[['dispersal.into.hot.rates.BSM.CI']]))],function(x)x[2])))),col=adjustcolor( "orange", alpha.f = 0.2)))
  lapply(BSMoutput.control.rates,function(x)lines(x=c(age:length(x[['dispersal.into.hot.rates.BSM.CI']])),y=unlist(lapply(x[['dispersal.into.hot.rates.BSM.CI']][c(age:length(x[['dispersal.into.hot.rates.BSM.CI']]))],function(x)x[3])),col=adjustcolor( "orange", alpha.f = 0.9)))
  axis(1,at=c(age:length(BSMoutput.real.rates[['dispersal.into.hot.rates.BSM.CI']])),labels=c((length(BSMoutput.real.rates[['insitu.realm.rates.BSM.CI']])-age+2):2))
  legend('topright',c('dispersal.into.hot','dispersal.into.hot.control'),col=c('red','orange'),lty=1,cex=.7,bty='n')
  #legend('topright',c('dispersal.into.hot','dispersal.into.realm','dispersal.into.hot.control','dispersal.into.realm.control'),col=c('red','blue','orange','purple'),lty=1,cex=.7,bty='n')
  
  #plot BAMM in situ rates realm
  plot(c(1,1),xlim=c(age,length(BSMoutput.real.rates[['insitu.realm.rates.BAMM.CI']])),ylim=c(0,0.4),type='n',xaxt='n',xlab='age',ylab='speciation rates in situ',main=paste(name,' BAMM in situ rates',sep=''))
  #polygon(x=c(c(age:length(BSMoutput.real.rates[['insitu.realm.rates.BAMM.CI']])),c(length(BSMoutput.real.rates[['insitu.realm.rates.BAMM.CI']]):age)),y=c(unlist(lapply(BSMoutput.real.rates[['insitu.realm.rates.BAMM.CI']][c(age:length(BSMoutput.real.rates[['insitu.realm.rates.BAMM.CI']]))],function(x)x[1])),rev(unlist(lapply(BSMoutput.real.rates[['insitu.realm.rates.BAMM.CI']][c(age:length(BSMoutput.real.rates[['insitu.realm.rates.BAMM.CI']]))],function(x)x[2])))),col=adjustcolor( "blue", alpha.f = 0.2))
  lines(x=c(age:length(BSMoutput.real.rates[['insitu.realm.rates.BAMM.CI']])),y=unlist(lapply(BSMoutput.real.rates[['insitu.realm.rates.BAMM.CI']][c(age:length(BSMoutput.real.rates[['insitu.realm.rates.BAMM.CI']]))],function(x)x[3])),col=adjustcolor( "blue", alpha.f = 0.9))
  #lapply(BSMoutput.control.rates,function(x)polygon(x=c(c(age:length(x[['insitu.realm.rates.BAMM.CI']])),c(length(x[['insitu.realm.rates.BAMM.CI']]):age)),y=c(unlist(lapply(x[['insitu.realm.rates.BAMM.CI']][c(age:length(x[['insitu.realm.rates.BAMM.CI']]))],function(x)x[1])),rev(unlist(lapply(x[['insitu.realm.rates.BAMM.CI']][c(age:length(x[['insitu.realm.rates.BAMM.CI']]))],function(x)x[2])))),col=adjustcolor( "purple", alpha.f = 0.2)))
  lapply(BSMoutput.control.rates,function(x)lines(x=c(age:length(x[['insitu.realm.rates.BAMM.CI']])),y=unlist(lapply(x[['insitu.realm.rates.BAMM.CI']][c(age:length(x[['insitu.realm.rates.BAMM.CI']]))],function(x)x[3])),col=adjustcolor( "green", alpha.f = 0.9)))
  axis(1,at=c(age:length(BSMoutput.real.rates[['insitu.realm.rates.BAMM.CI']])),labels=c((length(BSMoutput.real.rates[['insitu.realm.rates.BAMM.CI']])-age+2):2))
  legend('topright',c('insitu.nonhot','insitu.nonhot.control'),col=c('blue','green'),lty=1,cex=.7,bty='n')
  #plot BAMM in situ rates hot
  plot(c(1,1),xlim=c(age,length(BSMoutput.real.rates[['insitu.realm.rates.BAMM.CI']])),ylim=c(0,0.4),type='n',xaxt='n',xlab='age',ylab='speciation rates in situ',main=paste(name,' BAMM in situ rates',sep=''))
  #polygon(x=c(c(age:length(BSMoutput.real.rates[['insitu.hot.rates.BAMM.CI']])),c(length(BSMoutput.real.rates[['insitu.hot.rates.BAMM.CI']]):age)),y=c(unlist(lapply(BSMoutput.real.rates[['insitu.hot.rates.BAMM.CI']][c(age:length(BSMoutput.real.rates[['insitu.hot.rates.BAMM.CI']]))],function(x)x[1])),rev(unlist(lapply(BSMoutput.real.rates[['insitu.hot.rates.BAMM.CI']][c(age:length(BSMoutput.real.rates[['insitu.hot.rates.BAMM.CI']]))],function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
  lines(x=c(age:length(BSMoutput.real.rates[['insitu.hot.rates.BAMM.CI']])),y=unlist(lapply(BSMoutput.real.rates[['insitu.hot.rates.BAMM.CI']][c(age:length(BSMoutput.real.rates[['insitu.hot.rates.BAMM.CI']]))],function(x)x[3])),col=adjustcolor( "red", alpha.f = 0.9))
  #lapply(BSMoutput.control.rates,function(x)polygon(x=c(c(age:length(x[['insitu.hot.rates.BAMM.CI']])),c(length(x[['insitu.hot.rates.BAMM.CI']]):age)),y=c(unlist(lapply(x[['insitu.hot.rates.BAMM.CI']][c(age:length(x[['insitu.hot.rates.BAMM.CI']]))],function(x)x[1])),rev(unlist(lapply(x[['insitu.hot.rates.BAMM.CI']][c(age:length(x[['insitu.hot.rates.BAMM.CI']]))],function(x)x[2])))),col=adjustcolor( "orange", alpha.f = 0.2)))
  lapply(BSMoutput.control.rates,function(x)lines(x=c(age:length(x[['insitu.hot.rates.BAMM.CI']])),y=unlist(lapply(x[['insitu.hot.rates.BAMM.CI']][c(age:length(x[['insitu.hot.rates.BAMM.CI']]))],function(x)x[3])),col=adjustcolor( "orange", alpha.f = 0.9)))
  axis(1,at=c(age:length(BSMoutput.real.rates[['insitu.hot.rates.BAMM.CI']])),labels=c((length(BSMoutput.real.rates[['insitu.hot.rates.BAMM.CI']])-age+2):2))
  legend('topright',c('insitu.hot','insitu.hot.control'),col=c('red','orange'),lty=1,cex=.7,bty='n')
  #legend('topright',c('insitu.hot','insitu.nonhot','insitu.hot.control','insitu.nonhot.control'),col=c('red','blue','orange','purple'),lty=1,cex=.7,bty='n')
  
  #plot strict.BAMM in situ rates realm
  plot(c(1,1),xlim=c(age,length(BSMoutput.real.rates[['insitu.realm.rates.strict.BAMM.CI']])),ylim=c(0,0.4),type='n',xaxt='n',xlab='age',ylab='speciation rates in situ',main=paste(name,' strict.BAMM in situ rates',sep=''))
  #polygon(x=c(c(age:length(BSMoutput.real.rates[['insitu.realm.rates.strict.BAMM.CI']])),c(length(BSMoutput.real.rates[['insitu.realm.rates.strict.BAMM.CI']]):age)),y=c(unlist(lapply(BSMoutput.real.rates[['insitu.realm.rates.strict.BAMM.CI']][c(age:length(BSMoutput.real.rates[['insitu.realm.rates.strict.BAMM.CI']]))],function(x)x[1])),rev(unlist(lapply(BSMoutput.real.rates[['insitu.realm.rates.strict.BAMM.CI']][c(age:length(BSMoutput.real.rates[['insitu.realm.rates.strict.BAMM.CI']]))],function(x)x[2])))),col=adjustcolor( "blue", alpha.f = 0.2))
  lines(x=c(age:length(BSMoutput.real.rates[['insitu.realm.rates.strict.BAMM.CI']])),y=unlist(lapply(BSMoutput.real.rates[['insitu.realm.rates.strict.BAMM.CI']][c(age:length(BSMoutput.real.rates[['insitu.realm.rates.strict.BAMM.CI']]))],function(x)x[3])),col=adjustcolor( "blue", alpha.f = 0.9))
  #lapply(BSMoutput.control.rates,function(x)polygon(x=c(c(age:length(x[['insitu.realm.rates.strict.BAMM.CI']])),c(length(x[['insitu.realm.rates.strict.BAMM.CI']]):age)),y=c(unlist(lapply(x[['insitu.realm.rates.strict.BAMM.CI']][c(age:length(x[['insitu.realm.rates.strict.BAMM.CI']]))],function(x)x[1])),rev(unlist(lapply(x[['insitu.realm.rates.strict.BAMM.CI']][c(age:length(x[['insitu.realm.rates.strict.BAMM.CI']]))],function(x)x[2])))),col=adjustcolor( "purple", alpha.f = 0.2)))
  lapply(BSMoutput.control.rates,function(x)lines(x=c(age:length(x[['insitu.realm.rates.strict.BAMM.CI']])),y=unlist(lapply(x[['insitu.realm.rates.strict.BAMM.CI']][c(age:length(x[['insitu.realm.rates.strict.BAMM.CI']]))],function(x)x[3])),col=adjustcolor( "green", alpha.f = 0.9)))
  axis(1,at=c(age:length(BSMoutput.real.rates[['insitu.realm.rates.strict.BAMM.CI']])),labels=c((length(BSMoutput.real.rates[['insitu.realm.rates.strict.BAMM.CI']])-age+2):2))
  legend('topright',c('insitu.nonhot','insitu.nonhot.control'),col=c('blue','green'),lty=1,cex=.7,bty='n')
  #plot strict.BAMM in situ rates hot
  plot(c(1,1),xlim=c(age,length(BSMoutput.real.rates[['insitu.realm.rates.strict.BAMM.CI']])),ylim=c(0,0.4),type='n',xaxt='n',xlab='age',ylab='speciation rates in situ',main=paste(name,' strict.BAMM in situ rates',sep=''))
  #polygon(x=c(c(age:length(BSMoutput.real.rates[['insitu.hot.rates.strict.BAMM.CI']])),c(length(BSMoutput.real.rates[['insitu.hot.rates.strict.BAMM.CI']]):age)),y=c(unlist(lapply(BSMoutput.real.rates[['insitu.hot.rates.strict.BAMM.CI']][c(age:length(BSMoutput.real.rates[['insitu.hot.rates.strict.BAMM.CI']]))],function(x)x[1])),rev(unlist(lapply(BSMoutput.real.rates[['insitu.hot.rates.strict.BAMM.CI']][c(age:length(BSMoutput.real.rates[['insitu.hot.rates.strict.BAMM.CI']]))],function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
  lines(x=c(age:length(BSMoutput.real.rates[['insitu.hot.rates.strict.BAMM.CI']])),y=unlist(lapply(BSMoutput.real.rates[['insitu.hot.rates.strict.BAMM.CI']][c(age:length(BSMoutput.real.rates[['insitu.hot.rates.strict.BAMM.CI']]))],function(x)x[3])),col=adjustcolor( "red", alpha.f = 0.9))
  #lapply(BSMoutput.control.rates,function(x)polygon(x=c(c(age:length(x[['insitu.hot.rates.strict.BAMM.CI']])),c(length(x[['insitu.hot.rates.strict.BAMM.CI']]):age)),y=c(unlist(lapply(x[['insitu.hot.rates.strict.BAMM.CI']][c(age:length(x[['insitu.hot.rates.strict.BAMM.CI']]))],function(x)x[1])),rev(unlist(lapply(x[['insitu.hot.rates.strict.BAMM.CI']][c(age:length(x[['insitu.hot.rates.strict.BAMM.CI']]))],function(x)x[2])))),col=adjustcolor( "orange", alpha.f = 0.2)))
  lapply(BSMoutput.control.rates,function(x)lines(x=c(age:length(x[['insitu.hot.rates.strict.BAMM.CI']])),y=unlist(lapply(x[['insitu.hot.rates.strict.BAMM.CI']][c(age:length(x[['insitu.hot.rates.strict.BAMM.CI']]))],function(x)x[3])),col=adjustcolor( "orange", alpha.f = 0.9)))
  axis(1,at=c(age:length(BSMoutput.real.rates[['insitu.hot.rates.strict.BAMM.CI']])),labels=c((length(BSMoutput.real.rates[['insitu.hot.rates.strict.BAMM.CI']])-age+2):2))
  legend('topright',c('insitu.hot','insitu.hot.control'),col=c('red','orange'),lty=1,cex=.7,bty='n')
  #legend('topright',c('insitu.hot','insitu.nonhot','insitu.hot.control','insitu.nonhot.control'),col=c('red','blue','orange','purple'),lty=1,cex=.7,bty='n')
  
  #plot cumulative events
  #plot into and from events into.realm
  plot(c(1,1),xlim=c(age,length(BSMoutput.real.events[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']])),ylim=c(0,max(c(unlist(BSMoutput.real.events[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']]),unlist(BSMoutput.real.events[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']]),unlist(BSMoutput.real.events[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']]),unlist(BSMoutput.real.events[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']])))),type='n',xaxt='n',xlab='age',ylab='cumulative number of dispersal events',main=paste(name,' all events',sep=''))
  #polygon(x=c(c(age:length(BSMoutput.real.events[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']])),c(length(BSMoutput.real.events[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']]):age)),y=c(unlist(lapply(BSMoutput.real.events[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']][c(age:length(BSMoutput.real.events[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']]))],function(x)x[1])),rev(unlist(lapply(BSMoutput.real.events[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']][c(age:length(BSMoutput.real.events[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']]))],function(x)x[2])))),col=adjustcolor( "blue", alpha.f = 0.2))
  lines(x=c(age:length(BSMoutput.real.events[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']])),y=unlist(lapply(BSMoutput.real.events[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']][c(age:length(BSMoutput.real.events[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']]))],function(x)x[3])),col=adjustcolor( "blue", alpha.f = 0.9))
  #lapply(BSMoutput.control.events,function(x)polygon(x=c(c(age:length(x[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']])),c(length(x[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']]):age)),y=c(unlist(lapply(x[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']][c(age:length(x[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']]))],function(x)x[1])),rev(unlist(lapply(x[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']][c(age:length(x[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']]))],function(x)x[2])))),col=adjustcolor( "purple", alpha.f = 0.2)))
  lapply(BSMoutput.control.events,function(x)lines(x=c(age:length(x[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']])),y=unlist(lapply(x[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']][c(age:length(x[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']]))],function(x)x[3])),col=adjustcolor( "green", alpha.f = 0.9)))
  legend('topleft',c('into.non.hot.realm','into.non.hot.realm.control'),col=c('blue','green'),lty=1,cex=.7,bty='n')
  axis(1,at=c(age:length(BSMoutput.real.events[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']])),labels=c((length(BSMoutput.real.events[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']])-age+2):2))
  #plot into and from events into.hot
  plot(c(1,1),xlim=c(age,length(BSMoutput.real.events[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']])),ylim=c(0,max(c(unlist(BSMoutput.real.events[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']]),unlist(BSMoutput.real.events[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']]),unlist(BSMoutput.real.events[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']]),unlist(BSMoutput.real.events[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']])))),type='n',xaxt='n',xlab='age',ylab='cumulative number of dispersal events',main=paste(name,' all events',sep=''))
  #polygon(x=c(c(age:length(BSMoutput.real.events[['dispersal.into.hot.time.bins.counts.rev.cumsum.CI']])),c(length(BSMoutput.real.events[['dispersal.into.hot.time.bins.counts.rev.cumsum.CI']]):age)),y=c(unlist(lapply(BSMoutput.real.events[['dispersal.into.hot.time.bins.counts.rev.cumsum.CI']][c(age:length(BSMoutput.real.events[['dispersal.into.hot.time.bins.counts.rev.cumsum.CI']]))],function(x)x[1])),rev(unlist(lapply(BSMoutput.real.events[['dispersal.into.hot.time.bins.counts.rev.cumsum.CI']][c(age:length(BSMoutput.real.events[['dispersal.into.hot.time.bins.counts.rev.cumsum.CI']]))],function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
  lines(x=c(age:length(BSMoutput.real.events[['dispersal.into.hot.time.bins.counts.rev.cumsum.CI']])),y=unlist(lapply(BSMoutput.real.events[['dispersal.into.hot.time.bins.counts.rev.cumsum.CI']][c(age:length(BSMoutput.real.events[['dispersal.into.hot.time.bins.counts.rev.cumsum.CI']]))],function(x)x[3])),col=adjustcolor( "red", alpha.f = 0.9))
  #lapply(BSMoutput.control.events,function(x)polygon(x=c(c(age:length(x[['dispersal.into.hot.time.bins.counts.rev.cumsum.CI']])),c(length(x[['dispersal.into.hot.time.bins.counts.rev.cumsum.CI']]):age)),y=c(unlist(lapply(x[['dispersal.into.hot.time.bins.counts.rev.cumsum.CI']][c(age:length(x[['dispersal.into.hot.time.bins.counts.rev.cumsum.CI']]))],function(x)x[1])),rev(unlist(lapply(x[['dispersal.into.hot.time.bins.counts.rev.cumsum.CI']][c(age:length(x[['dispersal.into.hot.time.bins.counts.rev.cumsum.CI']]))],function(x)x[2])))),col=adjustcolor( "orange", alpha.f = 0.2)))
  lapply(BSMoutput.control.events,function(x)lines(x=c(age:length(x[['dispersal.into.hot.time.bins.counts.rev.cumsum.CI']])),y=unlist(lapply(x[['dispersal.into.hot.time.bins.counts.rev.cumsum.CI']][c(age:length(x[['dispersal.into.hot.time.bins.counts.rev.cumsum.CI']]))],function(x)x[3])),col=adjustcolor( "orange", alpha.f = 0.9)))
  legend('topleft',c('into.hot','into.hot.control'),col=c('red','orange'),lty=1,cex=.7,bty='n')
  axis(1,at=c(age:length(BSMoutput.real.events[['dispersal.into.hot.time.bins.counts.rev.cumsum.CI']])),labels=c((length(BSMoutput.real.events[['dispersal.into.hot.time.bins.counts.rev.cumsum.CI']])-age+2):2))
  #plot into and from events from.realm
  plot(c(1,1),xlim=c(age,length(BSMoutput.real.events[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']])),ylim=c(0,max(c(unlist(BSMoutput.real.events[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']]),unlist(BSMoutput.real.events[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']]),unlist(BSMoutput.real.events[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']]),unlist(BSMoutput.real.events[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']])))),type='n',xaxt='n',xlab='age',ylab='cumulative number of dispersal events',main=paste(name,' all events',sep=''))
  #polygon(x=c(c(age:length(BSMoutput.real.events[['dispersal.from.realm.time.bins.counts.rev.cumsum.CI']])),c(length(BSMoutput.real.events[['dispersal.from.realm.time.bins.counts.rev.cumsum.CI']]):age)),y=c(unlist(lapply(BSMoutput.real.events[['dispersal.from.realm.time.bins.counts.rev.cumsum.CI']][c(age:length(BSMoutput.real.events[['dispersal.from.realm.time.bins.counts.rev.cumsum.CI']]))],function(x)x[1])),rev(unlist(lapply(BSMoutput.real.events[['dispersal.from.realm.time.bins.counts.rev.cumsum.CI']][c(age:length(BSMoutput.real.events[['dispersal.from.realm.time.bins.counts.rev.cumsum.CI']]))],function(x)x[2])))),col=adjustcolor( "blue", alpha.f = 0.2))
  lines(x=c(age:length(BSMoutput.real.events[['dispersal.from.realm.time.bins.counts.rev.cumsum.CI']])),y=unlist(lapply(BSMoutput.real.events[['dispersal.from.realm.time.bins.counts.rev.cumsum.CI']][c(age:length(BSMoutput.real.events[['dispersal.from.realm.time.bins.counts.rev.cumsum.CI']]))],function(x)x[3])),col=adjustcolor( "blue", alpha.f = 0.9))
  #lapply(BSMoutput.control.events,function(x)polygon(x=c(c(age:length(x[['dispersal.from.realm.time.bins.counts.rev.cumsum.CI']])),c(length(x[['dispersal.from.realm.time.bins.counts.rev.cumsum.CI']]):age)),y=c(unlist(lapply(x[['dispersal.from.realm.time.bins.counts.rev.cumsum.CI']][c(age:length(x[['dispersal.from.realm.time.bins.counts.rev.cumsum.CI']]))],function(x)x[1])),rev(unlist(lapply(x[['dispersal.from.realm.time.bins.counts.rev.cumsum.CI']][c(age:length(x[['dispersal.from.realm.time.bins.counts.rev.cumsum.CI']]))],function(x)x[2])))),col=adjustcolor( "purple", alpha.f = 0.2)))
  lapply(BSMoutput.control.events,function(x)lines(x=c(age:length(x[['dispersal.from.realm.time.bins.counts.rev.cumsum.CI']])),y=unlist(lapply(x[['dispersal.from.realm.time.bins.counts.rev.cumsum.CI']][c(age:length(x[['dispersal.from.realm.time.bins.counts.rev.cumsum.CI']]))],function(x)x[3])),col=adjustcolor( "green", alpha.f = 0.9)))
  legend('topleft',c('from.non.hot.realm','from.non.hot.realm.control'),col=c('blue','green'),lty=1,cex=.7,bty='n')
  axis(1,at=c(age:length(BSMoutput.real.events[['dispersal.from.realm.time.bins.counts.rev.cumsum.CI']])),labels=c((length(BSMoutput.real.events[['dispersal.from.realm.time.bins.counts.rev.cumsum.CI']])-age+2):2))
  #plot into and from events from.hot
  plot(c(1,1),xlim=c(age,length(BSMoutput.real.events[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']])),ylim=c(0,max(c(unlist(BSMoutput.real.events[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']]),unlist(BSMoutput.real.events[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']]),unlist(BSMoutput.real.events[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']]),unlist(BSMoutput.real.events[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']])))),type='n',xaxt='n',xlab='age',ylab='cumulative number of dispersal events',main=paste(name,' all events',sep=''))
  #polygon(x=c(c(age:length(BSMoutput.real.events[['dispersal.from.hot.time.bins.counts.rev.cumsum.CI']])),c(length(BSMoutput.real.events[['dispersal.from.hot.time.bins.counts.rev.cumsum.CI']]):age)),y=c(unlist(lapply(BSMoutput.real.events[['dispersal.from.hot.time.bins.counts.rev.cumsum.CI']][c(age:length(BSMoutput.real.events[['dispersal.from.hot.time.bins.counts.rev.cumsum.CI']]))],function(x)x[1])),rev(unlist(lapply(BSMoutput.real.events[['dispersal.from.hot.time.bins.counts.rev.cumsum.CI']][c(age:length(BSMoutput.real.events[['dispersal.from.hot.time.bins.counts.rev.cumsum.CI']]))],function(x)x[2])))),col=adjustcolor( "blue", alpha.f = 0.2))
  lines(x=c(age:length(BSMoutput.real.events[['dispersal.from.hot.time.bins.counts.rev.cumsum.CI']])),y=unlist(lapply(BSMoutput.real.events[['dispersal.from.hot.time.bins.counts.rev.cumsum.CI']][c(age:length(BSMoutput.real.events[['dispersal.from.hot.time.bins.counts.rev.cumsum.CI']]))],function(x)x[3])),col=adjustcolor( "red", alpha.f = 0.9))
  #lapply(BSMoutput.control.events,function(x)polygon(x=c(c(age:length(x[['dispersal.from.hot.time.bins.counts.rev.cumsum.CI']])),c(length(x[['dispersal.from.hot.time.bins.counts.rev.cumsum.CI']]):age)),y=c(unlist(lapply(x[['dispersal.from.hot.time.bins.counts.rev.cumsum.CI']][c(age:length(x[['dispersal.from.hot.time.bins.counts.rev.cumsum.CI']]))],function(x)x[1])),rev(unlist(lapply(x[['dispersal.from.hot.time.bins.counts.rev.cumsum.CI']][c(age:length(x[['dispersal.from.hot.time.bins.counts.rev.cumsum.CI']]))],function(x)x[2])))),col=adjustcolor( "purple", alpha.f = 0.2)))
  lapply(BSMoutput.control.events,function(x)lines(x=c(age:length(x[['dispersal.from.hot.time.bins.counts.rev.cumsum.CI']])),y=unlist(lapply(x[['dispersal.from.hot.time.bins.counts.rev.cumsum.CI']][c(age:length(x[['dispersal.from.hot.time.bins.counts.rev.cumsum.CI']]))],function(x)x[3])),col=adjustcolor( "orange", alpha.f = 0.9)))
  legend('topleft',c('from.hot','from.hot.control'),col=c('red','orange'),lty=1,cex=.7,bty='n')
  axis(1,at=c(age:length(BSMoutput.real.events[['dispersal.from.hot.time.bins.counts.rev.cumsum.CI']])),labels=c((length(BSMoutput.real.events[['dispersal.from.hot.time.bins.counts.rev.cumsum.CI']])-age+2):2))
  
  #plot hot into realm and viceversa realm.into.hot
  plot(c(1,1),xlim=c(age,length(BSMoutput.real.events[['realm.into.hot.time.bins.counts.rev.cumsum.CI']])),ylim=c(0,max(c(unlist(BSMoutput.control.events[['realm.into.hot.time.bins.counts.rev.cumsum.CI']]),unlist(BSMoutput.control.events[['hot.into.realm.time.bins.counts.rev.cumsum.CI']]),unlist(BSMoutput.real.events[['realm.into.hot.time.bins.counts.rev.cumsum.CI']]),unlist(BSMoutput.real.events[['realm.into.hot.time.bins.counts.rev.cumsum.CI']])))),type='n',xaxt='n',xlab='age',ylab='cumulative number of realm-dispersal events',main=paste(name,' all events',sep=''))
  #polygon(x=c(c(age:length(BSMoutput.real.events[['realm.into.hot.time.bins.counts.rev.cumsum.CI']])),c(length(BSMoutput.real.events[['realm.into.hot.time.bins.counts.rev.cumsum.CI']]):age)),y=c(unlist(lapply(BSMoutput.real.events[['realm.into.hot.time.bins.counts.rev.cumsum.CI']][c(age:length(BSMoutput.real.events[['realm.into.hot.time.bins.counts.rev.cumsum.CI']]))],function(x)x[1])),rev(unlist(lapply(BSMoutput.real.events[['realm.into.hot.time.bins.counts.rev.cumsum.CI']][c(age:length(BSMoutput.real.events[['realm.into.hot.time.bins.counts.rev.cumsum.CI']]))],function(x)x[2])))),col=adjustcolor( "blue", alpha.f = 0.2))
  lines(x=c(age:length(BSMoutput.real.events[['realm.into.hot.time.bins.counts.rev.cumsum.CI']])),y=unlist(lapply(BSMoutput.real.events[['realm.into.hot.time.bins.counts.rev.cumsum.CI']][c(age:length(BSMoutput.real.events[['realm.into.hot.time.bins.counts.rev.cumsum.CI']]))],function(x)x[3])),col=adjustcolor( "blue", alpha.f = 0.9))
  #lapply(BSMoutput.control.events,function(x)polygon(x=c(c(age:length(x[['realm.into.hot.time.bins.counts.rev.cumsum.CI']])),c(length(x[['realm.into.hot.time.bins.counts.rev.cumsum.CI']]):age)),y=c(unlist(lapply(x[['realm.into.hot.time.bins.counts.rev.cumsum.CI']][c(age:length(x[['realm.into.hot.time.bins.counts.rev.cumsum.CI']]))],function(x)x[1])),rev(unlist(lapply(x[['realm.into.hot.time.bins.counts.rev.cumsum.CI']][c(age:length(x[['realm.into.hot.time.bins.counts.rev.cumsum.CI']]))],function(x)x[2])))),col=adjustcolor( "purple", alpha.f = 0.2)))
  lapply(BSMoutput.control.events,function(x)lines(x=c(age:length(x[['realm.into.hot.time.bins.counts.rev.cumsum.CI']])),y=unlist(lapply(x[['realm.into.hot.time.bins.counts.rev.cumsum.CI']][c(age:length(x[['realm.into.hot.time.bins.counts.rev.cumsum.CI']]))],function(x)x[3])),col=adjustcolor( "green", alpha.f = 0.9)))
  legend('topleft',c('realm.into.hot','realm.into.hot.control'),col=c('blue','green'),lty=1,cex=.7,bty='n')
  axis(1,at=c(age:length(BSMoutput.real.events[['realm.into.hot.time.bins.counts.rev.cumsum.CI']])),labels=c((length(BSMoutput.real.events[['realm.into.hot.time.bins.counts.rev.cumsum.CI']])-age+2):2))
  #plot hot into realm and viceversa hot.into.realm
  plot(c(1,1),xlim=c(age,length(BSMoutput.real.events[['realm.into.hot.time.bins.counts.rev.cumsum.CI']])),ylim=c(0,max(c(unlist(BSMoutput.control.events[['realm.into.hot.time.bins.counts.rev.cumsum.CI']]),unlist(BSMoutput.control.events[['hot.into.realm.time.bins.counts.rev.cumsum.CI']]),unlist(BSMoutput.real.events[['realm.into.hot.time.bins.counts.rev.cumsum.CI']]),unlist(BSMoutput.real.events[['realm.into.hot.time.bins.counts.rev.cumsum.CI']])))),type='n',xaxt='n',xlab='age',ylab='cumulative number of realm-dispersal events',main=paste(name,' all events',sep=''))
  #polygon(x=c(c(age:length(BSMoutput.real.events[['hot.into.realm.time.bins.counts.rev.cumsum.CI']])),c(length(BSMoutput.real.events[['hot.into.realm.time.bins.counts.rev.cumsum.CI']]):age)),y=c(unlist(lapply(BSMoutput.real.events[['hot.into.realm.time.bins.counts.rev.cumsum.CI']][c(age:length(BSMoutput.real.events[['hot.into.realm.time.bins.counts.rev.cumsum.CI']]))],function(x)x[1])),rev(unlist(lapply(BSMoutput.real.events[['hot.into.realm.time.bins.counts.rev.cumsum.CI']][c(age:length(BSMoutput.real.events[['hot.into.realm.time.bins.counts.rev.cumsum.CI']]))],function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
  lines(x=c(age:length(BSMoutput.real.events[['hot.into.realm.time.bins.counts.rev.cumsum.CI']])),y=unlist(lapply(BSMoutput.real.events[['hot.into.realm.time.bins.counts.rev.cumsum.CI']][c(age:length(BSMoutput.real.events[['hot.into.realm.time.bins.counts.rev.cumsum.CI']]))],function(x)x[3])),col=adjustcolor( "red", alpha.f = 0.9))
  #lapply(BSMoutput.control.events,function(x)polygon(x=c(c(age:length(x[['hot.into.realm.time.bins.counts.rev.cumsum.CI']])),c(length(x[['hot.into.realm.time.bins.counts.rev.cumsum.CI']]):age)),y=c(unlist(lapply(x[['hot.into.realm.time.bins.counts.rev.cumsum.CI']][c(age:length(x[['hot.into.realm.time.bins.counts.rev.cumsum.CI']]))],function(x)x[1])),rev(unlist(lapply(x[['hot.into.realm.time.bins.counts.rev.cumsum.CI']][c(age:length(x[['hot.into.realm.time.bins.counts.rev.cumsum.CI']]))],function(x)x[2])))),col=adjustcolor( "orange", alpha.f = 0.2)))
  lapply(BSMoutput.control.events,function(x)lines(x=c(age:length(x[['hot.into.realm.time.bins.counts.rev.cumsum.CI']])),y=unlist(lapply(x[['hot.into.realm.time.bins.counts.rev.cumsum.CI']][c(age:length(x[['hot.into.realm.time.bins.counts.rev.cumsum.CI']]))],function(x)x[3])),col=adjustcolor( "orange", alpha.f = 0.9)))
  legend('topleft',c('hot.into.realm','hot.into.realm.control'),col=c('red','orange'),lty=1,cex=.7,bty='n')
  axis(1,at=c(age:length(BSMoutput.real.events[['realm.into.hot.time.bins.counts.rev.cumsum.CI']])),labels=c((length(BSMoutput.real.events[['realm.into.hot.time.bins.counts.rev.cumsum.CI']])-age+2):2))
  #axis(1,at=c(age:length(BSMoutput.real.events[['hot.into.realm.time.bins.counts.rev.cumsum.CI']])),labels=c((length(BSMoutput.real.events[['hot.into.realm.time.bins.counts.rev.cumsum.CI']])-age+2):2))
  
  #plot in situ events realm
  plot(c(1,1),xlim=c(age,length(BSMoutput.real.events[['insitu.realm.time.bins.counts.rev.cumsum.CI']])),ylim=c(0,max(c(unlist(BSMoutput.real.events[['insitu.realm.time.bins.counts.rev.cumsum.CI']]),unlist(BSMoutput.control.events[['insitu.hot.time.bins.counts.rev.cumsum.CI']]),unlist(BSMoutput.control.events[['insitu.realm.time.bins.counts.rev.cumsum.CI']]),unlist(BSMoutput.real.events[['insitu.hot.time.bins.counts.rev.cumsum.CI']])))),type='n',xaxt='n',xlab='age',ylab='cumulative number of insitu cladogenesis events',main=paste(name,' all events',sep=''))
  #polygon(x=c(c(age:length(BSMoutput.real.events[['insitu.realm.time.bins.counts.rev.cumsum.CI']])),c(length(BSMoutput.real.events[['insitu.realm.time.bins.counts.rev.cumsum.CI']]):age)),y=c(unlist(lapply(BSMoutput.real.events[['insitu.realm.time.bins.counts.rev.cumsum.CI']][c(age:length(BSMoutput.real.events[['insitu.realm.time.bins.counts.rev.cumsum.CI']]))],function(x)x[1])),rev(unlist(lapply(BSMoutput.real.events[['insitu.realm.time.bins.counts.rev.cumsum.CI']][c(age:length(BSMoutput.real.events[['insitu.realm.time.bins.counts.rev.cumsum.CI']]))],function(x)x[2])))),col=adjustcolor( "blue", alpha.f = 0.2))
  lines(x=c(age:length(BSMoutput.real.events[['insitu.realm.time.bins.counts.rev.cumsum.CI']])),y=unlist(lapply(BSMoutput.real.events[['insitu.realm.time.bins.counts.rev.cumsum.CI']][c(age:length(BSMoutput.real.events[['insitu.realm.time.bins.counts.rev.cumsum.CI']]))],function(x)x[3])),col=adjustcolor( "blue", alpha.f = 0.9))
  #lapply(BSMoutput.control.events,function(x)polygon(x=c(c(age:length(x[['insitu.realm.time.bins.counts.rev.cumsum.CI']])),c(length(x[['insitu.realm.time.bins.counts.rev.cumsum.CI']]):age)),y=c(unlist(lapply(x[['insitu.realm.time.bins.counts.rev.cumsum.CI']][c(age:length(x[['insitu.realm.time.bins.counts.rev.cumsum.CI']]))],function(x)x[1])),rev(unlist(lapply(x[['insitu.realm.time.bins.counts.rev.cumsum.CI']][c(age:length(x[['insitu.realm.time.bins.counts.rev.cumsum.CI']]))],function(x)x[2])))),col=adjustcolor( "purple", alpha.f = 0.2)))
  lapply(BSMoutput.control.events,function(x)lines(x=c(age:length(x[['insitu.realm.time.bins.counts.rev.cumsum.CI']])),y=unlist(lapply(x[['insitu.realm.time.bins.counts.rev.cumsum.CI']][c(age:length(x[['insitu.realm.time.bins.counts.rev.cumsum.CI']]))],function(x)x[3])),col=adjustcolor( "green", alpha.f = 0.9)))
  legend('topleft',c('insitu.realm','insitu.realm.control'),col=c('blue','green'),lty=1,cex=.7,bty='n')
  axis(1,at=c(age:length(BSMoutput.real.events[['insitu.realm.time.bins.counts.rev.cumsum.CI']])),labels=c((length(BSMoutput.real.events[['insitu.realm.time.bins.counts.rev.cumsum.CI']])-age+2):2))
  #plot in situ events hot
  plot(c(1,1),xlim=c(age,length(BSMoutput.real.events[['insitu.realm.time.bins.counts.rev.cumsum.CI']])),ylim=c(0,max(c(unlist(BSMoutput.real.events[['insitu.realm.time.bins.counts.rev.cumsum.CI']]),unlist(BSMoutput.control.events[['insitu.hot.time.bins.counts.rev.cumsum.CI']]),unlist(BSMoutput.control.events[['insitu.realm.time.bins.counts.rev.cumsum.CI']]),unlist(BSMoutput.real.events[['insitu.hot.time.bins.counts.rev.cumsum.CI']])))),type='n',xaxt='n',xlab='age',ylab='cumulative number of insitu cladogenesis events',main=paste(name,' all events',sep=''))
  #polygon(x=c(c(age:length(BSMoutput.real.events[['insitu.hot.time.bins.counts.rev.cumsum.CI']])),c(length(BSMoutput.real.events[['insitu.hot.time.bins.counts.rev.cumsum.CI']]):age)),y=c(unlist(lapply(BSMoutput.real.events[['insitu.hot.time.bins.counts.rev.cumsum.CI']][c(age:length(BSMoutput.real.events[['insitu.hot.time.bins.counts.rev.cumsum.CI']]))],function(x)x[1])),rev(unlist(lapply(BSMoutput.real.events[['insitu.hot.time.bins.counts.rev.cumsum.CI']][c(age:length(BSMoutput.real.events[['insitu.hot.time.bins.counts.rev.cumsum.CI']]))],function(x)x[2])))),col=adjustcolor( "red", alpha.f = 0.2))
  lines(x=c(age:length(BSMoutput.real.events[['insitu.hot.time.bins.counts.rev.cumsum.CI']])),y=unlist(lapply(BSMoutput.real.events[['insitu.hot.time.bins.counts.rev.cumsum.CI']][c(age:length(BSMoutput.real.events[['insitu.hot.time.bins.counts.rev.cumsum.CI']]))],function(x)x[3])),col=adjustcolor( "red", alpha.f = 0.9))
  #lapply(BSMoutput.control.events,function(x)polygon(x=c(c(age:length(x[['insitu.hot.time.bins.counts.rev.cumsum.CI']])),c(length(x[['insitu.hot.time.bins.counts.rev.cumsum.CI']]):age)),y=c(unlist(lapply(x[['insitu.hot.time.bins.counts.rev.cumsum.CI']][c(age:length(x[['insitu.hot.time.bins.counts.rev.cumsum.CI']]))],function(x)x[1])),rev(unlist(lapply(x[['insitu.hot.time.bins.counts.rev.cumsum.CI']][c(age:length(x[['insitu.hot.time.bins.counts.rev.cumsum.CI']]))],function(x)x[2])))),col=adjustcolor( "orange", alpha.f = 0.2)))
  lapply(BSMoutput.control.events,function(x)lines(x=c(age:length(x[['insitu.hot.time.bins.counts.rev.cumsum.CI']])),y=unlist(lapply(x[['insitu.hot.time.bins.counts.rev.cumsum.CI']][c(age:length(x[['insitu.hot.time.bins.counts.rev.cumsum.CI']]))],function(x)x[3])),col=adjustcolor( "orange", alpha.f = 0.9)))
  #legend('topleft',c('insitu.realm','insitu.hot','insitu.realm.control','insitu.hot.control'),col=c('blue','red','purple','orange'),lty=1,cex=.7,bty='n')
  legend('topleft',c('insitu.hot','insitu.hot.control'),col=c('blue','red','purple','orange'),lty=1,cex=.7,bty='n')
  axis(1,at=c(age:length(BSMoutput.real.events[['insitu.hot.time.bins.counts.rev.cumsum.CI']])),labels=c((length(BSMoutput.real.events[['insitu.hot.time.bins.counts.rev.cumsum.CI']])-age+2):2))
  
  #plot logged cumulative events
  #plot into and from events into.realm
  plot(c(1,1),xlim=c(age,length(BSMoutput.real.events[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']])),ylim=c(0,log10(max(c(unlist(BSMoutput.real.events[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']]),unlist(BSMoutput.real.events[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']]),unlist(BSMoutput.real.events[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']]),unlist(BSMoutput.real.events[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']]))))),yaxt='n',type='n',xaxt='n',xlab='age',ylab='cumulative number of dispersal events',main=paste(name,' all events',sep=''))
  #polygon(x=c(c(age:length(BSMoutput.real.events[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']])),c(length(BSMoutput.real.events[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']]):age)),y=log10(c(unlist(lapply(BSMoutput.real.events[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']][c(age:length(BSMoutput.real.events[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']]))],function(x)x[1])),rev(unlist(lapply(BSMoutput.real.events[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']][c(age:length(BSMoutput.real.events[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']]))],function(x)x[2]))))),col=adjustcolor( "blue", alpha.f = 0.2))
  lines(x=c(age:length(BSMoutput.real.events[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']])),y=log10(unlist(lapply(BSMoutput.real.events[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']][c(age:length(BSMoutput.real.events[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']]))],function(x)x[3]))),col=adjustcolor( "blue", alpha.f = 0.9))
  #lapply(BSMoutput.control.events,function(x)polygon(x=c(c(age:length(x[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']])),c(length(x[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']]):age)),y=log10(c(unlist(lapply(x[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']][c(age:length(x[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']]))],function(x)x[1])),rev(unlist(lapply(x[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']][c(age:length(x[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']]))],function(x)x[2]))))),col=adjustcolor( "purple", alpha.f = 0.2)))
  lapply(BSMoutput.control.events,function(x)lines(x=c(age:length(x[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']])),y=log10(unlist(lapply(x[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']][c(age:length(x[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']]))],function(x)x[3]))),col=adjustcolor( "green", alpha.f = 0.9)))
  legend('topleft',c('into.non.hot.realm','into.non.hot.realm.control'),col=c('blue','green'),lty=1,cex=.7,bty='n')
  axis(1,at=c(age:length(BSMoutput.real.events[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']])),labels=c((length(BSMoutput.real.events[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']])-age+2):2))
  axis(2,at=log10(c(1,5,25,100,500)),labels=c('1','5','25','100','500'),las=2)
  #plot into and from events into.hot
  plot(c(1,1),xlim=c(age,length(BSMoutput.real.events[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']])),ylim=c(0,log10(max(c(unlist(BSMoutput.real.events[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']]),unlist(BSMoutput.real.events[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']]),unlist(BSMoutput.real.events[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']]),unlist(BSMoutput.real.events[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']]))))),yaxt='n',type='n',xaxt='n',xlab='age',ylab='cumulative number of dispersal events',main=paste(name,' all events',sep=''))
  #polygon(x=c(c(age:length(BSMoutput.real.events[['dispersal.into.hot.time.bins.counts.rev.cumsum.CI']])),c(length(BSMoutput.real.events[['dispersal.into.hot.time.bins.counts.rev.cumsum.CI']]):age)),y=log10(c(unlist(lapply(BSMoutput.real.events[['dispersal.into.hot.time.bins.counts.rev.cumsum.CI']][c(age:length(BSMoutput.real.events[['dispersal.into.hot.time.bins.counts.rev.cumsum.CI']]))],function(x)x[1])),rev(unlist(lapply(BSMoutput.real.events[['dispersal.into.hot.time.bins.counts.rev.cumsum.CI']][c(age:length(BSMoutput.real.events[['dispersal.into.hot.time.bins.counts.rev.cumsum.CI']]))],function(x)x[2]))))),col=adjustcolor( "red", alpha.f = 0.2))
  lines(x=c(age:length(BSMoutput.real.events[['dispersal.into.hot.time.bins.counts.rev.cumsum.CI']])),y=log10(unlist(lapply(BSMoutput.real.events[['dispersal.into.hot.time.bins.counts.rev.cumsum.CI']][c(age:length(BSMoutput.real.events[['dispersal.into.hot.time.bins.counts.rev.cumsum.CI']]))],function(x)x[3]))),col=adjustcolor( "red", alpha.f = 0.9))
  #lapply(BSMoutput.control.events,function(x)polygon(x=c(c(age:length(x[['dispersal.into.hot.time.bins.counts.rev.cumsum.CI']])),c(length(x[['dispersal.into.hot.time.bins.counts.rev.cumsum.CI']]):age)),y=log10(c(unlist(lapply(x[['dispersal.into.hot.time.bins.counts.rev.cumsum.CI']][c(age:length(x[['dispersal.into.hot.time.bins.counts.rev.cumsum.CI']]))],function(x)x[1])),rev(unlist(lapply(x[['dispersal.into.hot.time.bins.counts.rev.cumsum.CI']][c(age:length(x[['dispersal.into.hot.time.bins.counts.rev.cumsum.CI']]))],function(x)x[2]))))),col=adjustcolor( "orange", alpha.f = 0.2)))
  lapply(BSMoutput.control.events,function(x)lines(x=c(age:length(x[['dispersal.into.hot.time.bins.counts.rev.cumsum.CI']])),y=log10(unlist(lapply(x[['dispersal.into.hot.time.bins.counts.rev.cumsum.CI']][c(age:length(x[['dispersal.into.hot.time.bins.counts.rev.cumsum.CI']]))],function(x)x[3]))),col=adjustcolor( "orange", alpha.f = 0.9)))
  legend('topleft',c('into.hot','into.hot.control'),col=c('red','orange'),lty=1,cex=.7,bty='n')
  axis(1,at=c(age:length(BSMoutput.real.events[['dispersal.into.hot.time.bins.counts.rev.cumsum.CI']])),labels=c((length(BSMoutput.real.events[['dispersal.into.hot.time.bins.counts.rev.cumsum.CI']])-age+2):2))
  axis(2,at=log10(c(1,5,25,100,500)),labels=c('1','5','25','100','500'),las=2)
  #plot into and from events from.realm
  plot(c(1,1),xlim=c(age,length(BSMoutput.real.events[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']])),ylim=c(0,log10(max(c(unlist(BSMoutput.real.events[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']]),unlist(BSMoutput.real.events[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']]),unlist(BSMoutput.real.events[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']]),unlist(BSMoutput.real.events[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']]))))),yaxt='n',type='n',xaxt='n',xlab='age',ylab='cumulative number of dispersal events',main=paste(name,' all events',sep=''))
  #polygon(x=c(c(age:length(BSMoutput.real.events[['dispersal.from.realm.time.bins.counts.rev.cumsum.CI']])),c(length(BSMoutput.real.events[['dispersal.from.realm.time.bins.counts.rev.cumsum.CI']]):age)),y=log10(c(unlist(lapply(BSMoutput.real.events[['dispersal.from.realm.time.bins.counts.rev.cumsum.CI']][c(age:length(BSMoutput.real.events[['dispersal.from.realm.time.bins.counts.rev.cumsum.CI']]))],function(x)x[1])),rev(unlist(lapply(BSMoutput.real.events[['dispersal.from.realm.time.bins.counts.rev.cumsum.CI']][c(age:length(BSMoutput.real.events[['dispersal.from.realm.time.bins.counts.rev.cumsum.CI']]))],function(x)x[2]))))),col=adjustcolor( "blue", alpha.f = 0.2))
  lines(x=c(age:length(BSMoutput.real.events[['dispersal.from.realm.time.bins.counts.rev.cumsum.CI']])),y=log10(unlist(lapply(BSMoutput.real.events[['dispersal.from.realm.time.bins.counts.rev.cumsum.CI']][c(age:length(BSMoutput.real.events[['dispersal.from.realm.time.bins.counts.rev.cumsum.CI']]))],function(x)x[3]))),col=adjustcolor( "blue", alpha.f = 0.9))
  #lapply(BSMoutput.control.events,function(x)polygon(x=c(c(age:length(x[['dispersal.from.realm.time.bins.counts.rev.cumsum.CI']])),c(length(x[['dispersal.from.realm.time.bins.counts.rev.cumsum.CI']]):age)),y=log10(c(unlist(lapply(x[['dispersal.from.realm.time.bins.counts.rev.cumsum.CI']][c(age:length(x[['dispersal.from.realm.time.bins.counts.rev.cumsum.CI']]))],function(x)x[1])),rev(unlist(lapply(x[['dispersal.from.realm.time.bins.counts.rev.cumsum.CI']][c(age:length(x[['dispersal.from.realm.time.bins.counts.rev.cumsum.CI']]))],function(x)x[2]))))),col=adjustcolor( "purple", alpha.f = 0.2)))
  lapply(BSMoutput.control.events,function(x)lines(x=c(age:length(x[['dispersal.from.realm.time.bins.counts.rev.cumsum.CI']])),y=log10(unlist(lapply(x[['dispersal.from.realm.time.bins.counts.rev.cumsum.CI']][c(age:length(x[['dispersal.from.realm.time.bins.counts.rev.cumsum.CI']]))],function(x)x[3]))),col=adjustcolor( "green", alpha.f = 0.9)))
  legend('topleft',c('from.non.hot.realm','from.non.hot.realm.control'),col=c('blue','green'),lty=1,cex=.7,bty='n')
  axis(1,at=c(age:length(BSMoutput.real.events[['dispersal.from.realm.time.bins.counts.rev.cumsum.CI']])),labels=c((length(BSMoutput.real.events[['dispersal.from.realm.time.bins.counts.rev.cumsum.CI']])-age+2):2))
  axis(2,at=log10(c(1,5,25,100,500)),labels=c('1','5','25','100','500'),las=2)
  
  #plot into and from events from.hot
  plot(c(1,1),xlim=c(age,length(BSMoutput.real.events[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']])),ylim=c(0,log10(max(c(unlist(BSMoutput.real.events[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']]),unlist(BSMoutput.real.events[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']]),unlist(BSMoutput.real.events[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']]),unlist(BSMoutput.real.events[['dispersal.into.realm.time.bins.counts.rev.cumsum.CI']]))))),yaxt='n',type='n',xaxt='n',xlab='age',ylab='cumulative number of dispersal events',main=paste(name,' all events',sep=''))
  #polygon(x=c(c(age:length(BSMoutput.real.events[['dispersal.from.hot.time.bins.counts.rev.cumsum.CI']])),c(length(BSMoutput.real.events[['dispersal.from.hot.time.bins.counts.rev.cumsum.CI']]):age)),y=log10(c(unlist(lapply(BSMoutput.real.events[['dispersal.from.hot.time.bins.counts.rev.cumsum.CI']][c(age:length(BSMoutput.real.events[['dispersal.from.hot.time.bins.counts.rev.cumsum.CI']]))],function(x)x[1])),rev(unlist(lapply(BSMoutput.real.events[['dispersal.from.hot.time.bins.counts.rev.cumsum.CI']][c(age:length(BSMoutput.real.events[['dispersal.from.hot.time.bins.counts.rev.cumsum.CI']]))],function(x)x[2]))))),col=adjustcolor( "blue", alpha.f = 0.2))
  lines(x=c(age:length(BSMoutput.real.events[['dispersal.from.hot.time.bins.counts.rev.cumsum.CI']])),y=log10(unlist(lapply(BSMoutput.real.events[['dispersal.from.hot.time.bins.counts.rev.cumsum.CI']][c(age:length(BSMoutput.real.events[['dispersal.from.hot.time.bins.counts.rev.cumsum.CI']]))],function(x)x[3]))),col=adjustcolor( "red", alpha.f = 0.9))
  #lapply(BSMoutput.control.events,function(x)polygon(x=c(c(age:length(x[['dispersal.from.hot.time.bins.counts.rev.cumsum.CI']])),c(length(x[['dispersal.from.hot.time.bins.counts.rev.cumsum.CI']]):age)),y=log10(c(unlist(lapply(x[['dispersal.from.hot.time.bins.counts.rev.cumsum.CI']][c(age:length(x[['dispersal.from.hot.time.bins.counts.rev.cumsum.CI']]))],function(x)x[1])),rev(unlist(lapply(x[['dispersal.from.hot.time.bins.counts.rev.cumsum.CI']][c(age:length(x[['dispersal.from.hot.time.bins.counts.rev.cumsum.CI']]))],function(x)x[2]))))),col=adjustcolor( "purple", alpha.f = 0.2)))
  lapply(BSMoutput.control.events,function(x)lines(x=c(age:length(x[['dispersal.from.hot.time.bins.counts.rev.cumsum.CI']])),y=log10(unlist(lapply(x[['dispersal.from.hot.time.bins.counts.rev.cumsum.CI']][c(age:length(x[['dispersal.from.hot.time.bins.counts.rev.cumsum.CI']]))],function(x)x[3]))),col=adjustcolor( "orange", alpha.f = 0.9)))
  legend('topleft',c('from.hot','from.hot.control'),col=c('red','orange'),lty=1,cex=.7,bty='n')
  axis(1,at=c(age:length(BSMoutput.real.events[['dispersal.from.hot.time.bins.counts.rev.cumsum.CI']])),labels=c((length(BSMoutput.real.events[['dispersal.from.hot.time.bins.counts.rev.cumsum.CI']])-age+2):2))
  axis(2,at=log10(c(1,5,25,100,500)),labels=c('1','5','25','100','500'),las=2)
  #plot hot into realm and viceversa realm.into.hot
  plot(c(1,1),xlim=c(age,length(BSMoutput.real.events[['realm.into.hot.time.bins.counts.rev.cumsum.CI']])),ylim=c(0,log10(max(c(unlist(BSMoutput.control.events[['realm.into.hot.time.bins.counts.rev.cumsum.CI']]),unlist(BSMoutput.control.events[['hot.into.realm.time.bins.counts.rev.cumsum.CI']]),unlist(BSMoutput.real.events[['realm.into.hot.time.bins.counts.rev.cumsum.CI']]),unlist(BSMoutput.real.events[['realm.into.hot.time.bins.counts.rev.cumsum.CI']]))))),yaxt='n',type='n',xaxt='n',xlab='age',ylab='cumulative number of realm-dispersal events',main=paste(name,' all events',sep=''))
  #polygon(x=c(c(age:length(BSMoutput.real.events[['realm.into.hot.time.bins.counts.rev.cumsum.CI']])),c(length(BSMoutput.real.events[['realm.into.hot.time.bins.counts.rev.cumsum.CI']]):age)),y=log10(c(unlist(lapply(BSMoutput.real.events[['realm.into.hot.time.bins.counts.rev.cumsum.CI']][c(age:length(BSMoutput.real.events[['realm.into.hot.time.bins.counts.rev.cumsum.CI']]))],function(x)x[1])),rev(unlist(lapply(BSMoutput.real.events[['realm.into.hot.time.bins.counts.rev.cumsum.CI']][c(age:length(BSMoutput.real.events[['realm.into.hot.time.bins.counts.rev.cumsum.CI']]))],function(x)x[2]))))),col=adjustcolor( "blue", alpha.f = 0.2))
  lines(x=c(age:length(BSMoutput.real.events[['realm.into.hot.time.bins.counts.rev.cumsum.CI']])),y=log10(unlist(lapply(BSMoutput.real.events[['realm.into.hot.time.bins.counts.rev.cumsum.CI']][c(age:length(BSMoutput.real.events[['realm.into.hot.time.bins.counts.rev.cumsum.CI']]))],function(x)x[3]))),col=adjustcolor( "blue", alpha.f = 0.9))
  #lapply(BSMoutput.control.events,function(x)polygon(x=c(c(age:length(x[['realm.into.hot.time.bins.counts.rev.cumsum.CI']])),c(length(x[['realm.into.hot.time.bins.counts.rev.cumsum.CI']]):age)),y=log10(c(unlist(lapply(x[['realm.into.hot.time.bins.counts.rev.cumsum.CI']][c(age:length(x[['realm.into.hot.time.bins.counts.rev.cumsum.CI']]))],function(x)x[1])),rev(unlist(lapply(x[['realm.into.hot.time.bins.counts.rev.cumsum.CI']][c(age:length(x[['realm.into.hot.time.bins.counts.rev.cumsum.CI']]))],function(x)x[2]))))),col=adjustcolor( "purple", alpha.f = 0.2)))
  lapply(BSMoutput.control.events,function(x)lines(x=c(age:length(x[['realm.into.hot.time.bins.counts.rev.cumsum.CI']])),y=log10(unlist(lapply(x[['realm.into.hot.time.bins.counts.rev.cumsum.CI']][c(age:length(x[['realm.into.hot.time.bins.counts.rev.cumsum.CI']]))],function(x)x[3]))),col=adjustcolor( "green", alpha.f = 0.9)))
  legend('topleft',c('realm.into.hot','realm.into.hot.control'),col=c('blue','green'),lty=1,cex=.7,bty='n')
  axis(1,at=c(age:length(BSMoutput.real.events[['realm.into.hot.time.bins.counts.rev.cumsum.CI']])),labels=c((length(BSMoutput.real.events[['realm.into.hot.time.bins.counts.rev.cumsum.CI']])-age+2):2))
  axis(2,at=log10(c(1,5,25,100,500)),labels=c('1','5','25','100','500'),las=2)
  #plot hot into realm and viceversa hot.into.realm
  plot(c(1,1),xlim=c(age,length(BSMoutput.real.events[['realm.into.hot.time.bins.counts.rev.cumsum.CI']])),ylim=c(0,log10(max(c(unlist(BSMoutput.control.events[['realm.into.hot.time.bins.counts.rev.cumsum.CI']]),unlist(BSMoutput.control.events[['hot.into.realm.time.bins.counts.rev.cumsum.CI']]),unlist(BSMoutput.real.events[['realm.into.hot.time.bins.counts.rev.cumsum.CI']]),unlist(BSMoutput.real.events[['realm.into.hot.time.bins.counts.rev.cumsum.CI']]))))),yaxt='n',type='n',xaxt='n',xlab='age',ylab='cumulative number of realm-dispersal events',main=paste(name,' all events',sep=''))
  #polygon(x=c(c(age:length(BSMoutput.real.events[['hot.into.realm.time.bins.counts.rev.cumsum.CI']])),c(length(BSMoutput.real.events[['hot.into.realm.time.bins.counts.rev.cumsum.CI']]):age)),y=log10(c(unlist(lapply(BSMoutput.real.events[['hot.into.realm.time.bins.counts.rev.cumsum.CI']][c(age:length(BSMoutput.real.events[['hot.into.realm.time.bins.counts.rev.cumsum.CI']]))],function(x)x[1])),rev(unlist(lapply(BSMoutput.real.events[['hot.into.realm.time.bins.counts.rev.cumsum.CI']][c(age:length(BSMoutput.real.events[['hot.into.realm.time.bins.counts.rev.cumsum.CI']]))],function(x)x[2]))))),col=adjustcolor( "red", alpha.f = 0.2))
  lines(x=c(age:length(BSMoutput.real.events[['hot.into.realm.time.bins.counts.rev.cumsum.CI']])),y=log10(unlist(lapply(BSMoutput.real.events[['hot.into.realm.time.bins.counts.rev.cumsum.CI']][c(age:length(BSMoutput.real.events[['hot.into.realm.time.bins.counts.rev.cumsum.CI']]))],function(x)x[3]))),col=adjustcolor( "red", alpha.f = 0.9))
  #lapply(BSMoutput.control.events,function(x)polygon(x=c(c(age:length(x[['hot.into.realm.time.bins.counts.rev.cumsum.CI']])),c(length(x[['hot.into.realm.time.bins.counts.rev.cumsum.CI']]):age)),y=log10(c(unlist(lapply(x[['hot.into.realm.time.bins.counts.rev.cumsum.CI']][c(age:length(x[['hot.into.realm.time.bins.counts.rev.cumsum.CI']]))],function(x)x[1])),rev(unlist(lapply(x[['hot.into.realm.time.bins.counts.rev.cumsum.CI']][c(age:length(x[['hot.into.realm.time.bins.counts.rev.cumsum.CI']]))],function(x)x[2]))))),col=adjustcolor( "orange", alpha.f = 0.2)))
  lapply(BSMoutput.control.events,function(x)lines(x=c(age:length(x[['hot.into.realm.time.bins.counts.rev.cumsum.CI']])),y=log10(unlist(lapply(x[['hot.into.realm.time.bins.counts.rev.cumsum.CI']][c(age:length(x[['hot.into.realm.time.bins.counts.rev.cumsum.CI']]))],function(x)x[3]))),col=adjustcolor( "orange", alpha.f = 0.9)))
  legend('topleft',c('hot.into.realm','hot.into.realm.control'),col=c('red','orange'),lty=1,cex=.7,bty='n')
  axis(1,at=c(age:length(BSMoutput.real.events[['realm.into.hot.time.bins.counts.rev.cumsum.CI']])),labels=c((length(BSMoutput.real.events[['realm.into.hot.time.bins.counts.rev.cumsum.CI']])-age+2):2))
  axis(2,at=log10(c(1,5,25,100,500)),labels=c('1','5','25','100','500'),las=2)
  
  #axis(1,at=c(age:length(BSMoutput.real.events[['hot.into.realm.time.bins.counts.rev.cumsum.CI']])),labels=c((length(BSMoutput.real.events[['hot.into.realm.time.bins.counts.rev.cumsum.CI']])-age+2):2))
  
  #plot in situ events realm
  plot(c(1,1),xlim=c(age,length(BSMoutput.real.events[['insitu.realm.time.bins.counts.rev.cumsum.CI']])),ylim=c(0,log10(max(c(unlist(BSMoutput.real.events[['insitu.realm.time.bins.counts.rev.cumsum.CI']]),unlist(BSMoutput.control.events[['insitu.hot.time.bins.counts.rev.cumsum.CI']]),unlist(BSMoutput.control.events[['insitu.realm.time.bins.counts.rev.cumsum.CI']]),unlist(BSMoutput.real.events[['insitu.hot.time.bins.counts.rev.cumsum.CI']]))))),yaxt='n',type='n',xaxt='n',xlab='age',ylab='cumulative number of insitu cladogenesis events',main=paste(name,' all events',sep=''))
  #polygon(x=c(c(age:length(BSMoutput.real.events[['insitu.realm.time.bins.counts.rev.cumsum.CI']])),c(length(BSMoutput.real.events[['insitu.realm.time.bins.counts.rev.cumsum.CI']]):age)),y=log10(c(unlist(lapply(BSMoutput.real.events[['insitu.realm.time.bins.counts.rev.cumsum.CI']][c(age:length(BSMoutput.real.events[['insitu.realm.time.bins.counts.rev.cumsum.CI']]))],function(x)x[1])),rev(unlist(lapply(BSMoutput.real.events[['insitu.realm.time.bins.counts.rev.cumsum.CI']][c(age:length(BSMoutput.real.events[['insitu.realm.time.bins.counts.rev.cumsum.CI']]))],function(x)x[2]))))),col=adjustcolor( "blue", alpha.f = 0.2))
  lines(x=c(age:length(BSMoutput.real.events[['insitu.realm.time.bins.counts.rev.cumsum.CI']])),y=log10(unlist(lapply(BSMoutput.real.events[['insitu.realm.time.bins.counts.rev.cumsum.CI']][c(age:length(BSMoutput.real.events[['insitu.realm.time.bins.counts.rev.cumsum.CI']]))],function(x)x[3]))),col=adjustcolor( "blue", alpha.f = 0.9))
  #lapply(BSMoutput.control.events,function(x)polygon(x=c(c(age:length(x[['insitu.realm.time.bins.counts.rev.cumsum.CI']])),c(length(x[['insitu.realm.time.bins.counts.rev.cumsum.CI']]):age)),y=log10(c(unlist(lapply(x[['insitu.realm.time.bins.counts.rev.cumsum.CI']][c(age:length(x[['insitu.realm.time.bins.counts.rev.cumsum.CI']]))],function(x)x[1])),rev(unlist(lapply(x[['insitu.realm.time.bins.counts.rev.cumsum.CI']][c(age:length(x[['insitu.realm.time.bins.counts.rev.cumsum.CI']]))],function(x)x[2]))))),col=adjustcolor( "purple", alpha.f = 0.2)))
  lapply(BSMoutput.control.events,function(x)lines(x=c(age:length(x[['insitu.realm.time.bins.counts.rev.cumsum.CI']])),y=log10(unlist(lapply(x[['insitu.realm.time.bins.counts.rev.cumsum.CI']][c(age:length(x[['insitu.realm.time.bins.counts.rev.cumsum.CI']]))],function(x)x[3]))),col=adjustcolor( "green", alpha.f = 0.9)))
  legend('topleft',c('insitu.realm','insitu.realm.control'),col=c('blue','green'),lty=1,cex=.7,bty='n')
  axis(1,at=c(age:length(BSMoutput.real.events[['insitu.realm.time.bins.counts.rev.cumsum.CI']])),labels=c((length(BSMoutput.real.events[['insitu.realm.time.bins.counts.rev.cumsum.CI']])-age+2):2))
  axis(2,at=log10(c(1,5,25,100,500)),labels=c('1','5','25','100','500'),las=2)
  #plot in situ events hot
  plot(c(1,1),xlim=c(age,length(BSMoutput.real.events[['insitu.realm.time.bins.counts.rev.cumsum.CI']])),ylim=c(0,log10(max(c(unlist(BSMoutput.real.events[['insitu.realm.time.bins.counts.rev.cumsum.CI']]),unlist(BSMoutput.control.events[['insitu.hot.time.bins.counts.rev.cumsum.CI']]),unlist(BSMoutput.control.events[['insitu.realm.time.bins.counts.rev.cumsum.CI']]),unlist(BSMoutput.real.events[['insitu.hot.time.bins.counts.rev.cumsum.CI']]))))),yaxt='n',type='n',xaxt='n',xlab='age',ylab='cumulative number of insitu cladogenesis events',main=paste(name,' all events',sep=''))
  #polygon(x=c(c(age:length(BSMoutput.real.events[['insitu.hot.time.bins.counts.rev.cumsum.CI']])),c(length(BSMoutput.real.events[['insitu.hot.time.bins.counts.rev.cumsum.CI']]):age)),y=log10(c(unlist(lapply(BSMoutput.real.events[['insitu.hot.time.bins.counts.rev.cumsum.CI']][c(age:length(BSMoutput.real.events[['insitu.hot.time.bins.counts.rev.cumsum.CI']]))],function(x)x[1])),rev(unlist(lapply(BSMoutput.real.events[['insitu.hot.time.bins.counts.rev.cumsum.CI']][c(age:length(BSMoutput.real.events[['insitu.hot.time.bins.counts.rev.cumsum.CI']]))],function(x)x[2]))))),col=adjustcolor( "red", alpha.f = 0.2))
  lines(x=c(age:length(BSMoutput.real.events[['insitu.hot.time.bins.counts.rev.cumsum.CI']])),y=log10(unlist(lapply(BSMoutput.real.events[['insitu.hot.time.bins.counts.rev.cumsum.CI']][c(age:length(BSMoutput.real.events[['insitu.hot.time.bins.counts.rev.cumsum.CI']]))],function(x)x[3]))),col=adjustcolor( "red", alpha.f = 0.9))
  #lapply(BSMoutput.control.events,function(x)polygon(x=c(c(age:length(x[['insitu.hot.time.bins.counts.rev.cumsum.CI']])),c(length(x[['insitu.hot.time.bins.counts.rev.cumsum.CI']]):age)),y=log10(c(unlist(lapply(x[['insitu.hot.time.bins.counts.rev.cumsum.CI']][c(age:length(x[['insitu.hot.time.bins.counts.rev.cumsum.CI']]))],function(x)x[1])),rev(unlist(lapply(x[['insitu.hot.time.bins.counts.rev.cumsum.CI']][c(age:length(x[['insitu.hot.time.bins.counts.rev.cumsum.CI']]))],function(x)x[2]))))),col=adjustcolor( "orange", alpha.f = 0.2)))
  lapply(BSMoutput.control.events,function(x)lines(x=c(age:length(x[['insitu.hot.time.bins.counts.rev.cumsum.CI']])),y=log10(unlist(lapply(x[['insitu.hot.time.bins.counts.rev.cumsum.CI']][c(age:length(x[['insitu.hot.time.bins.counts.rev.cumsum.CI']]))],function(x)x[3]))),col=adjustcolor( "orange", alpha.f = 0.9)))
  #legend('topleft',c('insitu.realm','insitu.hot','insitu.realm.control','insitu.hot.control'),col=c('blue','red','purple','orange'),lty=1,cex=.7,bty='n')
  legend('topleft',c('insitu.hot','insitu.hot.control'),col=c('red','orange'),lty=1,cex=.7,bty='n')
  axis(1,at=c(age:length(BSMoutput.real.events[['insitu.hot.time.bins.counts.rev.cumsum.CI']])),labels=c((length(BSMoutput.real.events[['insitu.hot.time.bins.counts.rev.cumsum.CI']])-age+2):2))
  axis(2,at=log10(c(1,5,25,100,500)),labels=c('1','5','25','100','500'),las=2)
  
  
}


plot_ages_colonisation_realms<-function(ages.realm.colonisation,ages.hot.colonisation,name){
  plot(c(1,1),xlim=c(0,7),ylim=c(0,max(c(unlist(ages.realm.colonisation),unlist(ages.hot.colonisation)))+5),xlab='',ylab='My',xaxt='n',yaxt='n',type='n',main=paste('colonisation age - ',name,sep=''))
  for(i in 1:length(ages.realm.colonisation)){
    vioplot(ages.hot.colonisation[[i]],col=adjustcolor( "red", alpha.f = 0.2),at=i,add=T,rectCol = adjustcolor( "red", alpha.f = 0.5),lwd=0.001,colMed = adjustcolor( "red", alpha.f = 0.9),pchMed=19)
    vioplot(ages.realm.colonisation[[i]],col=adjustcolor( "blue", alpha.f = 0.2),at=i,add=T,rectCol = adjustcolor( "blue", alpha.f = 0.5),lwd=0.001,colMed = adjustcolor( "blue", alpha.f = 0.9),pchMed=19)
  }
  axis(1,at=c(1:6),labels=c(rep('',6)))
  text(x =c(1:6), par("usr")[3] - 0.2, labels = c('Afrotropical','Australasian','Indo-Malay','Nearctic','Neotropical','Palearctic'), srt = 45, pos = 1, xpd = TRUE)
  axis(2,las=2)
  
}

plot_difference_ages_colonisation_realms<-function(diff.ages.hot.minus.realm.colonisation,name){
  plot(c(1,1),xlim=c(0,7),ylim=c(-82,82),xlab='',ylab='My',xaxt='n',yaxt='n',type='n',main=paste('colonisation age (hotspot-non hotspot) - ',name,sep=''),yaxt='n')
  #plot(c(1,1),xlim=c(0,7),ylim=c(min(unlist(diff.ages.hot.minus.realm.colonisation))-2,max(unlist(diff.ages.hot.minus.realm.colonisation))+2),xlab='',ylab='My',xaxt='n',yaxt='n',type='n',main=paste('colonisation age (hotspot-non hotspot) - ',name,sep=''),yaxt='n')
  for(i in 1:length(diff.ages.hot.minus.realm.colonisation)){
    boxplot(diff.ages.hot.minus.realm.colonisation[[i]],col='white',at=i,add=T,outline=F,yaxt='n')
    points(x=jitter(rep(i,length(diff.ages.hot.minus.realm.colonisation[[i]])),3),y=diff.ages.hot.minus.realm.colonisation[[i]],pch=21,bg='white')
#    vioplot(ages.realm.colonisation[[i]],col=adjustcolor( "blue", alpha.f = 0.2),at=i,add=T,rectCol = adjustcolor( "blue", alpha.f = 0.5),lwd=0.001,colMed = adjustcolor( "blue", alpha.f = 0.9),pchMed=19)
  }
  axis(1,at=c(1:6),labels=c(rep('',6)))
  text(x =c(1:6), par("usr")[3] - 0.2, labels = c('Afrotropical','Australasian','Indo-Malay','Nearctic','Neotropical','Palearctic'), srt = 45, pos = 1, xpd = TRUE)
  axis(2,las=2)
  abline(h=0,lty=2)
  
}




