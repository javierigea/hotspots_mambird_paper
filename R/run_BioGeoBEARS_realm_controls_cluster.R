#run model in hydrogen cluster
.libPaths('/home/ecosys/ji247/R/x86_64-pc-linux-gnu-library/3.3')
setwd('/home/ecosys/ji247/hotspots_vertebrates/')


run_whole_realm_BioGeoBEARS_plusBSM_7areas_inhotoutrealm<-function(path,realm.name,number){
  setwd('/home/ji247/hotspots_vertebrates/')
  source('./R/BioGeoBEARS_run_models.R')
  path<-paste(path,realm.name,'/',number,'/',sep='')
  name<-realm.name
  names.trees<-list.files(path,pattern='.tree$')
  names.geography<-list.files(path,pattern='_geographyfile.txt')
  if(length(names.trees)!=length(names.geography)){
    cat('different lenghts, double check folder','\n')
  }
  list.aicc.files<-list.files(path,pattern='_AICc_rellike_formatted.txt')
  if(length(list.aicc.files)==0){
    for(i in 1:length(names.trees)){
    	run_BioGeoBEARS_models(treefile = names.trees[i],geographyfile =names.geography[i],path = path,name=sub(names.trees[i],pattern=paste('_',name,'.tree',sep=''),replacement=''))
	    setwd('/home/ji247/hotspots_vertebrates/')
    	list.aicc.files<-list.files(path,pattern='_AICc_rellike_formatted.txt') 
  	}
  }
  setwd('/home/ji247/hotspots_vertebrates/')
  

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

#######WRITE THIS IN A FUNCTION
args<-commandArgs(trailingOnly = TRUE)
path<-args[2]
realm.name<-args[3]
start.number<-as.numeric(args[4])
end.number<-as.numeric(args[5])

for (i in c(start.number:end.number)){
  run_whole_realm_BioGeoBEARS_plusBSM_7areas_inhotoutrealm(path=path, realm.name = realm.name,number=i)
}
#run_whole_realm_BioGeoBEARS_plusBSM_7areas_inhotoutrealm(realm.name = 'afro',number='10')
#run_whole_realm_BioGeoBEARS_plusBSM_7areas_inhotoutrealm(realm.name = 'austral',number='10')
#run_whole_realm_BioGeoBEARS_plusBSM_7areas_inhotoutrealm(realm.name = 'indo',number='22')
#run_whole_realm_BioGeoBEARS_plusBSM_7areas_inhotoutrealm(realm.name = 'nearctic',number='13')
#run_whole_realm_BioGeoBEARS_plusBSM_7areas_inhotoutrealm(realm.name = 'neotrop',number='43')
#run_whole_realm_BioGeoBEARS_plusBSM_7areas_inhotoutrealm(realm.name = 'palearctic',number='7')


