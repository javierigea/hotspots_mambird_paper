library(ape)
library(phangorn)
library(picante)
library(phytools)


####create folder structure####
dir.create('./output/')
dir.create('./raw_data/mammals/')
dir.create('./raw_data/birds/')
dir.create('./raw_data/mammals/trees/')
dir.create('./raw_data/birds/trees/')
#/raw_data/mammals/trees/ contains FritzTree.rs200k.100trees.tre (100 resampled trees from Kuhn et al Supp Info at https://doi.org/10.1111/j.2041-210X.2011.00103.x)
#/raw_data/birds/trees/ contains BirdzillaEricsonAllTrees (Ericson backbone; all trees from the Jetz et al 2012 posterior, from https://birdtree.org/downloads/)
dir.create('./output/mammals/')
dir.create('./output/birds/')
dir.create('./output/grids/')
dir.create('./output/grids/tables/')

####1) PHYLOGENETIC DATA PREP####
####prepare MCCtree for mammals####
####Rolland et al Plos Biol: they build a MCC tree with the 100 trees in Fritz file, but they use -keep with the nodes. MCCtree is almost equal to tree #23 in the Fritz 100 file
####system('/Applications2/BEAST/1.8.2/bin/treeannotator ./raw_data/mammals/trees/FritzTree.rs200k.100trees.tre ./output/mammals/trees/mammals_MCC_keep.nexus')
#for mammals
dir.create('./output/mammals/trees/')
system('/Applications2/BEAST/1.8.2/bin/treeannotator -heights ca ./raw_data/mammals/trees/FritzTree.rs200k.100trees.tre ./output/mammals/trees/mammals_MCC.nexus')
#converting MCC from Treeannotator (nexus) to Newick
nexus.tree<-read.nexus('./output/mammals/trees/mammals_MCC.nexus')
write.tree(nexus.tree,'./output/mammals/trees/mammals_MCC.tree')
#calibrate the MCC tree with Meredith dates (as in Rolland Plos Biol)
source('./R/calibrate_mammals_trees.R')
calibrated.tree.MCC<-calibrate_tree_Meredith(treefile='./output/mammals/trees/mammals_MCC.tree')
write.tree(calibrated.tree.MCC,file='./output/mammals/trees/mammals_MCC_calibrated.tree')

####use taxonomy to deal with synonyms for mammals####
#read in trees and use IUCN taxonomy to remove synonyms
mammal.tree<-read.tree('./output/mammals/trees/mammals_MCC_calibrated.tree')
source('./R/taxonomy_IUCN_birdlife.R')
mammal.tree.IUCN<-mammal.tree
mammal.tree.IUCN$tip.label<-iucn.taxonomy.synonyms(species=mammal.tree$tip.label,dictionaryfile = './raw_data/IUCNTaxonomy_Mammalia_30082017_version2017-1.csv')
#remove duplicated tips after correcting synonyms
tip.duplicated.mammal.tree.IUCN<-return.duplicate.tips(mammal.tree.IUCN)
mammal.tree.IUCN<-drop.tip(mammal.tree.IUCN,tip.duplicated.mammal.tree.IUCN)
write.tree(mammal.tree.IUCN,'./output/mammals/trees/mammals_tree_IUCN.tree')

####prepare MCCtree for birds####
dir.create('./output/birds/trees/')
#get 100 random trees from pseudoposterior
trees<-scan('./raw_data/birds/trees/BirdzillaEricsonAllTrees.tre',sep='\n',what='char')
#gets 100 trees evenly distributed from the pseudoposterior
trees.sub<-trees[round(seq(1,length(trees),length.out=100))]
writeLines(trees.sub,'./output/birds/trees/BirdzillaEricsonAllTrees_100.trees')
trees100<-read.tree('./output/birds/trees/BirdzillaEricsonAllTrees_100.trees')
writeNexus(trees100,'./output/birds/trees//BirdzillaEricsonAllTrees_100.nex')
system('/Applications2/BEAST/1.8.2/bin/treeannotator -heights ca ./output/birds/trees/BirdzillaEricsonAllTrees_100.nex ./output/birds/trees/birds_MCC.nexus')
#converting MCC from Treeannotator (nexus) to Newick
nexus.tree<-read.nexus('./output/birds/trees/birds_MCC.nexus')
write.tree(nexus.tree,'./output/birds/trees/birds_MCC.tree')

####use taxonomy to deal with synonyms for birds####
bird.tree<-read.tree('./output/birds/trees/birds_MCC.tree')
source('./R/taxonomy_IUCN_birdlife.R')
bird.tree.IUCN<-bird.tree
bird.tree.IUCN$tip.label<-iucn.taxonomy.synonyms(species=bird.tree$tip.label,dictionaryfile = './raw_data/IUCNTaxonomy_Aves_30082017_version2017-1.csv')
#remove duplicated tips after correcting synonyms
tip.duplicated.bird.tree.IUCN<-return.duplicate.tips(bird.tree.IUCN)
bird.tree.IUCN<-drop.tip(bird.tree.IUCN,tip.duplicated.bird.tree.IUCN)
write.tree(bird.tree.IUCN,'./output/birds/trees/birds_tree_IUCN.tree')

####mammal pseudoposterior time calibration####
#100 trees from the mammal pseudoposterior####
source('./R/calibrate_mammals_trees.R')
calibrate_pseudoposterior_Meredith(treesfile='./raw_data/mammals/trees/FritzTree.rs200k.100trees.tre')
#use IUCN taxonomy on posterior calibrated trees
trees100<-read.nexus('./output/mammals/trees/posterior_calibrated/FritzTree.rs200k.100trees_Meredithdates.trees')
source('./R/taxonomy_IUCN_birdlife.R')
trees100.IUCN<-lapply(trees100,function(x){x$tip.label<-iucn.taxonomy.synonyms(species=x$tip.label,dictionaryfile = './raw_data/IUCNTaxonomy_Mammalia_30082017_version2017-1.csv');tip.duplicated.mammal.tree.IUCN<-return.duplicate.tips(x);x<-drop.tip(x,tip.duplicated.mammal.tree.IUCN);return(x)})
write.nexus(trees100.IUCN,file='./output/mammals/trees/posterior_calibrated/FritzTree.rs200k.100trees_Meredithdates_IUCN.trees')

####bird pseudoposterior####
####use IUCN taxonomy on Jetz pseudoposterior 100 trees
bird100trees<-read.tree('./output/birds/trees/BirdzillaEricsonAllTrees_100.trees')
bird100trees.IUCN<-lapply(bird100trees,function(x){x$tip.label<-iucn.taxonomy.synonyms(species=x$tip.label,dictionaryfile = './raw_data/IUCNTaxonomy_Aves_30082017_version2017-1.csv');tip.duplicated.bird.tree.IUCN<-return.duplicate.tips(x);x<-drop.tip(x,tip.duplicated.bird.tree.IUCN);return(x)})
write.nexus(bird100trees.IUCN,file='./output/birds/trees/BirdzillaEricsonAllTrees_100_IUCN.trees')

#####2) SPATIAL DATA PREP & DELINEATING HOTSPOTS####
#####overlap mammal IUCN ranges layers with grid####
#run './R/overlap_realms_mammals_grid_cluster.R' on hydrogen
#/scripts/conscriptoR /home/ji247/hotspots_vertebrates/overlap_realms_mammals_grid_cluster.R -p32 (#for 32 cores,it takes ~20 hours)
#move the *_realms_species_gridoccurrence_table.txt and *_richness_grid_table.txt to ./output/mammals/tables/
#merge realm based data to world level
source('./R/merging_all_realm_occurrences.R')
merge_realms_speciesdata_into_world(path='./output/mammals/tables/')
#copy grid_*100.rds to /output/grids/
merge_realms_grid_into_world(path='./output/grids/')
#then use IUCN taxonomy on them (this is redundant)
source('./R/taxonomy_IUCN_birdlife.R')
mammals_grid.table<-read.table('./output/mammals/tables/100_all_realms_species_gridoccurrence_table.txt',header=T,sep='\t',stringsAsFactors = F)
mammals_grid.table.IUCN<-iucn.taxonomy.synonyms(species=mammals_grid.table$spp,dictionaryfile = './raw_data/IUCNTaxonomy_Mammalia_30082017_version2017-1.csv')
mammals_grid.table$spp<-mammals_grid.table.IUCN
#deal with duplicates: merge the ranges
duplicated.mammals_grid.table.spp<-mammals_grid.table$spp[duplicated(mammals_grid.table$spp)]
if(length(duplicated.birds_grid.table.spp)>0){
  mammals_grid.table<-merge_duplicates_in_speciesgridoccurrence_table(duplicated.species = duplicated.mammals_grid.table.spp,species.gridoccurrence.table = mammals_grid.table)  
}
write.table(mammals_grid.table,'./output/mammals/tables/100_all_realms_species_gridoccurrence_table.txt',sep='\t',quote=F,row.names=F)

#####overlap bird birdlife range layers with grid####
#BOTW.gbd has to be on './raw_data/birds/BOTW/BOTW.gdb/'; run './R/run_split_birds_gb_into_shp.R
#this will create folders 1,501,1001,1501...18001, each with 500 shapefiles
source('./R/run_split_birds_gdb_into_shp.R')
#zip + copy those folders to cluster
#run './R/overlap_realms_birds_grid_cluster.R' on hydrogen
#/scripts/conscriptoR /home/ji247/hotspots_vertebrates/overlap_realms_birds_grid_cluster.R -p48 (#for 48 cores,it takes 4-5 days)
#and then run commented final parts of overlap_realms_birds_grid_cluster.R
#move the *_realms_species_gridoccurrence_table.txt and *_richness_grid_table.txt to ./output/birds/tables/
#merge realm based data to world level
source('./R/merging_all_realm_occurrences.R')
merge_realms_speciesdata_into_world(path='./output/birds/tables/')
#copy grid_*100.rds to /output/grids/
merge_realms_grid_into_world(path='./output/grids/')
#then use IUCN taxonomy on them (this is redundant)
source('./R/taxonomy_IUCN_birdlife.R')
birds_grid.table<-read.table('./output/birds/tables/100_all_realms_species_gridoccurrence_table.txt',header=T,sep='\t',stringsAsFactors = F)
birds_grid.table.IUCN<-iucn.taxonomy.synonyms(species=birds_grid.table$spp,dictionaryfile = './raw_data/IUCNTaxonomy_Aves_30082017_version2017-1.csv')
birds_grid.table$spp<-birds_grid.table.IUCN
#deal with duplicates: merge the ranges
duplicated.birds_grid.table.spp<-birds_grid.table$spp[duplicated(birds_grid.table$spp)]
if(length(duplicated.birds_grid.table.spp)>0){
  birds_grid.table<-merge_duplicates_in_speciesgridoccurrence_table(duplicated.species = duplicated.birds_grid.table.spp,species.gridoccurrence.table = birds_grid.table)  
}
write.table(birds_grid.table,'./output/birds/tables/100_all_realms_species_gridoccurrence_table.txt',sep='\t',quote=F,row.names=F)

####define WE(weighted endemicity) hotspots for mammals, etc####
#get species names of species occurring in each cell
source('./R/grid_functions.R')
get_speciesnames_grid(species.grid.table = './output/mammals/tables/100_all_realms_species_gridoccurrence_table.txt',cell.grid.table = './output/mammals/tables/100_all_realms_richness_grid_table.txt',path = './output/mammals/tables/',name='mammals_all_realms')
#calculate weighted endemicity
source('./R/grid_functions.R')
generate_grid_weightedendemism(speciesgrid.tablefile = './output/mammals/tables/100_all_realms_species_gridoccurrence_table.txt',cellgrid.tablefile = './output/mammals/tables/100_all_realms_richness_grid_table.txt',path = './output/mammals/tables/',name = 'all_realms')
#define hotspots using weighted endemicity
#this generates a hotspot_grid
define_hotspots(grid.tablefile ='./output/mammals/tables/100_all_realms_realms_richness_wend_grid_table.txt',quantile = 0.8,variable.name = 'number.of.species.wend',path = './output/mammals/tables/',name = 'all_realms')

####define WE(weighted endemicity) hotspots for birds, etc####
#get species names of species occurring in each cell
source('./R/grid_functions.R')
get_speciesnames_grid(species.grid.table = './output/birds/tables/100_all_realms_species_gridoccurrence_table.txt',cell.grid.table = './output/birds/tables/100_all_realms_richness_grid_table.txt',path = './output/birds/tables/',name='birds_all_realms')
#calculate weighted endemicity
source('./R/grid_functions.R')
generate_grid_weightedendemism(speciesgrid.tablefile = './output/birds/tables/100_all_realms_species_gridoccurrence_table.txt',cellgrid.tablefile = './output/birds/tables/100_all_realms_richness_grid_table.txt',path = './output/birds/tables/',name = 'all_realms')
#define hotspots using weighted endemicity
#this generates a hotspot_grid
define_hotspots(grid.tablefile ='./output/birds/tables/100_all_realms_realms_richness_wend_grid_table.txt',quantile = 0.8,variable.name = 'number.of.species.wend',path = './output/birds/tables/',name = 'all_realms')

####3) ENVIRONMENTAL DATA PREP####
####analyse elevation####
source('./R/elevation.R')
#this generates a world dem, only need to run this line here once, it stores it in "./output/rasters/world_dem"
build_world_dem('./raw_data/gtopo30_dems/')
dem<-raster('./output/rasters/world_dem')
get_elevation_variables_grid(gridfile = './output/grids/grid_World_RealmsMerged_100.rds',dem = dem,path = './output/grids/tables/')

####analyse NPP####
library(raster)
source('./R/NPP.R')
NPP<-raster('./raw_data/npp-geotiff/npp_geotiff.tif')
get_variable_grid(gridfile ='./output/grids/grid_World_RealmsMerged_100.rds',raster=NPP,name='NPP',path='./output/grids/tables/')

####analyse number of habitats####
source('./R/nhabitats.R')
Nhabitats.raster<-raster('./raw_data/glccgbe20_tif/gbogegeo20.tif')
get_habitat_count_grid(gridfile = './output/grids/grid_World_RealmsMerged_100.rds',raster = Nhabitats.raster,name = '100_allrealms_nhabitats',path = './output/grids/tables/')

####analyse climate change velocity####
source('./R/climatechange_velocity.R')
#this calculates climate change velocity (LGM to present day) and stores rasters in rasters folder
#for bio1 (mean annual temp),bio12 (mean annual rainfall) and bio7 (mean annual temperature range)
generate_LGMCCV_rasters(path_current_raster = '../raw_data/current_2_5m/bio1.bil',path_LGM_raster = '../raw_data/LGM_2_5m/wc_2_5m_MIROC3.2_21k_bio_1.bil',type = 'temperature',name='bio1_LGMCCV')
generate_LGMCCV_rasters(path_current_raster = '../raw_data/current_2_5m/bio12.bil',path_LGM_raster = '../raw_data/LGM_2_5m/wc_2_5m_MIROC3.2_21k_bio_12.bil',type = 'rainfall',name='bio12_LGMCCV')
generate_LGMCCV_rasters(path_current_raster = '../raw_data/current_2_5m/bio7.bil',path_LGM_raster = '../raw_data/LGM_2_5m/wc_2_5m_MIROC3.2_21k_bio_7.bil',type = 'temperature',name='bio7_LGMCCV')
#getting the data for MAT.CCV
MAT.CCV<-raster('./output/rasters/bio1_LGMCCV')
get_variable_grid(gridfile = './output/grids/grid_World_RealmsMerged_100.rds',raster=MAT.CCV,path='./output/grids/tables/',name='MAT.CCV')
#getting the data for AR.CCV
AR.CCV<-raster('./output/rasters/bio12_LGMCCV')
get_variable_grid(gridfile = './output/grids/grid_World_RealmsMerged_100.rds',raster=AR.CCV,path='./output/grids/tables/',name='AR.CCV')
#getting the data for ATR.CCV
ATR.CCV<-raster('./output/rasters/bio7_LGMCCV')
get_variable_grid(gridfile = './output/grids/grid_World_RealmsMerged_100.rds',raster=ATR.CCV,path='./output/grids/tables/',name='ATR.CCV')

####analyse tectonic movements####
source('./R/tectonic_movements.R')
#this creates a shapefile using the World Grid
create_shp_from_grid(world.gridfile = './output/grids/grid_World_RealmsMerged_100.rds',path = './output/grids/')
#use Gplates to project the shapefile and reconstruct past tectonic movements (see Gplates_README.txt in ./other_scripts/)
#copy the xy reconstructions to './output/grids/reconstructions/'
#process reconstructions + get average standard deviation of distances across time of cell vs its neighbours (takes a couple of hours to run)
#output to table
get_average_stdev_tectonicmovements(reconstruction.path ='./output/grids/reconstructions/',world.gridfile = './output/grids/grid_World_RealmsMerged_100.rds',path = './output/grids/tables/' )

####compare environment in hotspots vs non hotspots for mammals
#for mammals, all environmental variables together + check colinearity
source('./R/linear_models_spatialautocorrelation.R')
source('./R/environmental_variables_sarlm.R')
#NPP
table.NPP.mammals<-prepare_linearmodel_table(vector.of.tables = c('./output/grids/tables/NPP_table.txt','./output/mammals/tables/100_all_realms_number.of.species.wend_0.8.txt'))
#elevation
table.hab.elevation.mammals<-prepare_linearmodel_table(vector.of.tables = c('./output/grids/tables/100_all_realms_elevation_table.txt','./output/grids/tables/100_allrealms_nhabitats_table.txt','./output/mammals/tables/100_all_realms_number.of.species.wend_0.8.txt'))
#past climate
table.CCV.mammals<-prepare_linearmodel_table(vector.of.tables = c('./output/grids/tables/AR.CCV_table.txt','./output/grids/tables/MAT.CCV_table.txt','./output/mammals/tables/100_all_realms_number.of.species.wend_0.8.txt'))
#tectonic movements over 65 Myr
table.tectonic.mammals<-prepare_linearmodel_table(vector.of.tables = c('./output/grids/tables/100_all_realms_tectonicmovement.txt','./output/mammals/tables/100_all_realms_number.of.species.wend_0.8.txt'))
#average distance
table.distance.mammals<-prepare_linearmodel_table(vector.of.tables = c('./output/mammals/tables/averagedistancetoclass_1000neighbours.txt','./output/mammals/tables/100_all_realms_number.of.species.wend_0.8.txt'))
#check correlations among variables
corr_variablesNPP(table.NPP = table.NPP.mammals,table.hab.elevation = table.hab.elevation.mammals,table.CCV = table.CCV.mammals,table.tectonic=table.tectonic.mammals)
#select variables where r < 0.7 (remove collinearity)
#variable.vector<-c('mean.NPP','mean.TRI','n.habitats','mean.AR.CCV','mean.MAT.CCV','tectonic.movement')
source('./R/heatmap_enviromentalvariables.R')

pdf('./output/mammals/plots/mammals_heatmap_allenvironmentscaledNPP_sarlm.pdf',width=10,height=7)
environmentalmodelsNPP_sarlm_variables_heatmap(table.NPP = table.NPP.mammals,table.hab.elevation = table.hab.elevation.mammals,table.CCV = table.CCV.mammals,table.tectonic=table.tectonic.mammals,variable.vector = variable.vector)
dev.off()

####4) DIVERSIFICATION ANALYSES (DR,BAMM)#####
####measure DR in mammal trees####
source('./R/measure_DR.R')
#for mammals
DR.mammals<-measure_DR_tree_table(treefile = "./output/mammals/mammals_tree_IUCN.tree")
#select terrestrial
DR.mammals.marine<-get.marine.species(species = DR.mammals$Species,dictionaryfile = './raw_data/IUCNTaxonomy_Mammalia_30082017_version2017-1.csv')
DR.mammals.marine<-c(DR.mammals.marine,'Platanista_minor')
DR.mammals.terrestrial<-DR.mammals[!(DR.mammals$Species%in%DR.mammals.marine),]
write.table(DR.mammals,'./output/mammals/DR_mammals_tree_IUCN.txt',sep='\t',quote=F,row.names=F)
write.table(DR.mammals.terrestrial,'./output/mammals/DR_mammals_terrestrial_tree_IUCN.txt',sep='\t',quote=F,row.names=F)
#for mammals pseudoposterior
mammals100trees<-read.nexus('./output/mammals/trees/posterior_calibrated/FritzTree.rs200k.100trees_Meredithdates_IUCN.trees')
DR.mammals.100<-lapply(mammals100trees,function(x)measure_DR_tree_table(x))
#select terrestrial
DR.mammals.marine<-get.marine.species(species = DR.mammals.100[[1]]$Species,dictionaryfile = './raw_data/IUCNTaxonomy_Mammalia_30082017_version2017-1.csv')
DR.mammals.marine<-c(DR.mammals.marine,'Platanista_minor')
DR.mammals.terrestrial.100<-lapply(DR.mammals.100,function(x)x[!(x$Species%in%DR.mammals.marine),])
DR.mammals.terrestrial<-DR.mammals[!(DR.mammals$Species%in%DR.mammals.marine),]
#save object
saveRDS(DR.mammals.terrestrial.100,file='./output/mammals/trees/DR.mammals.terrestrial.100.rds')
#get median DR across the pseudoposterior for each species
get_pseudoposterior_median_DRtable(DR.pseudoposterior.file='./output/mammals/trees/DR.mammals.terrestrial.100.rds',path='./output/mammals/tables/',name='mammals')
#compare DR from the MCC tree with median of the pseudoposterior
#this is part of Fig. S2
source('./R/measure_DR.R')
pdf('./output/mammals/plots/mammals_DRMCC_vs_pseudoposterior.pdf')
compare_DR_MCC_vs_medianpseudoposterior(DR.tablefile='./output/mammals/DR_mammals_terrestrial_tree_IUCN.txt',DR.pseudoposterior.tablefile='./all_realms_new_pseudoposteriorDRmedian.txt')
compare_DR_MCC_vs_allpseudoposterior(DR.tablefile = './output/mammals/DR_mammals_terrestrial_tree_IUCN.txt',DR.pseudoposterior.file ='./output/mammals/trees/DR.mammals.terrestrial.100.rds' )
dev.off()

####measure DR in bird trees####
source('./R/measure_DR.R')
#for birds
DR.birds<-measure_DR_tree_table(treefile = "./output/birds/birds_tree_IUCN.tree")
write.table(DR.birds,'./output/birds/DR_birds_tree_IUCN.txt',sep='\t',quote=F,row.names=F)
#for birds pseudoposterior
birds100trees<-read.nexus('./output/birds/trees/BirdzillaEricsonAllTrees_100_IUCN.trees')
counter<-0
DR.birds.100<-lapply(birds100trees,function(x){counter<<-counter+1;cat(counter,'\n');measure_DR_tree_table(x)})
#save object
saveRDS(DR.birds.100,file='./output/birds/trees/DR.birds.100.rds')
#get median DR across the pseudoposterior for each species
get_pseudoposterior_median_DRtable(DR.pseudoposterior.file='./output/birds/trees/DR.birds.100.rds',path='./output/birds/tables/',name='all_realms_birds')
#compare DR from the MCC tree with median of the pseudoposterior
#this is part of Fig. S2
source('./R/measure_DR.R')
pdf('./output/birds/plots/birds_DRMCC_vs_pseudoposterior.pdf')
compare_DR_MCC_vs_pseudoposterior(DR.tablefile='./output/birds/DR_birds_tree_IUCN.txt',DR.pseudoposterior.tablefile='./output/birds/tables/all_realms_birds_pseudoposteriorDRmedian.txt')
compare_DR_MCC_vs_allpseudoposterior(DR.tablefile = './output/birds/DR_birds_tree_IUCN.txt',DR.pseudoposterior.file ='./output/birds/trees/DR.birds.100.rds')
dev.off()

####BAMM analyses####
#prepare BAMM input
source('./R/BAMM_functions.R')
prepare_BAMM_input(treefile='./output/mammals/trees/mammals_tree_IUCN.tree',dictionaryfile = './raw_data/IUCNTaxonomy_Mammalia_30082017_version2017-1.csv',path='./output/mammals/trees/',name='mammals_IUCN_BAMM')
prepare_BAMM_input(treefile='./output/birds/birds_tree_IUCN.tree',dictionaryfile = './raw_data/IUCNTaxonomy_Aves_30082017_version2017-1.csv',path='./output/birds/trees/',name='birds_IUCN_BAMM')  
#for mammals it takes ~40 hours (50 million generations), more than enough for convergence (ESS:1200, could work with much less generations) (it took ~80 hours on 4 cores on node 8 - hydrogen)
#for mammals it takes ~21 hours for 30 million generations; for birds it takes 64 hours for 30 million generations
#check convergence
analyse_BAMM_convergence(mcmcout = './output/mammals/trees/mcmc_out_mammals_30m.txt',burnin=0.25)
analyse_BAMM_convergence(mcmcout = './output/birds/trees/mcmc_out_birds_30m.txt',burnin=0.25)
#get TipRates
get_tipRates_BAMM(treefile='./output/mammals/trees/mammals_tree_IUCN.tree',eventfile='./output/mammals/trees/event_data_mammals_30m.txt',burnin=0.25,path='./output/mammals/tables/',name='mammals_all_realms')
get_tipRates_BAMM(treefile='./output/birds/trees/birds_tree_IUCN.tree',eventfile='./output/birds/trees/event_data_birds_30m.txt',burnin=0.25,path='./output/birds/tables/',name='birds_all_realms')

####quartiles of DR for mammals####
#generate the DR quartile richness grid measures for mammals
source('./R/grid_functions.R')
generate_grid_richness_quartilesDR(dr.table='./output/mammals/DR_mammals_terrestrial_tree_IUCN.txt',species.grid.tablefile = './output/mammals/tables/100_all_realms_species_gridoccurrence_table.txt',cellgrid.tablefile='./output/mammals/tables/100_all_realms_richness_grid_table.txt',path='./output/mammals/tables/',name='all_realms')
#linear regressions of speciation rates (total richness vs richness in quartiles, etc)
#run linear model to predict quartile richness with total richness
build_lm_DRquartiles(dr.quartiles.gridtable='./output/mammals/tables/all_realms_quartilesDR_grid_table.txt',path='./output/mammals/tables/',name='all_realms')
#quartiles of DR for pseudoposterior
#save each table in DRpseudoposterior in a tablefile
DR.mammals.terrestrial.pseudoposterior<-readRDS('./output/mammals/trees/DR.mammals.terrestrial.100.rds')
for(i in 1:100){
  write.table(DR.mammals.terrestrial.pseudoposterior[[i]],file=paste('./output/mammals/DR_mammals_terrestrial_tree_IUCN_',i,'.txt',sep=''),sep='\t',quote=F,row.names=F)
}
source('./R/grid_functions.R')
for (i in 1:100){
  generate_grid_richness_quartilesDR(dr.table=paste('./output/mammals/DR_mammals_terrestrial_tree_IUCN_',i,'.txt',sep=''),species.grid.tablefile = './output/mammals/tables/100_all_realms_species_gridoccurrence_table.txt',cellgrid.tablefile='./output/mammals/tables/100_all_realms_richness_grid_table.txt',path='./output/mammals/tables/',name=paste('all_realms_',i,sep=''))
  #linear regressions of speciation rates (total richness vs richness in quartiles, etc)
  #run linear model to predict quartile richness with total richness
  build_lm_DRquartiles(dr.quartiles.gridtable=paste('./output/mammals/tables/all_realms_',i,'_quartilesDR_grid_table.txt',sep=''),path='./output/mammals/tables/',name=paste('all_realms_',i,sep=''))
  
}
#read all DR.pseudoposterior residual files for mammals
DR.pseudoposterior.residual.tablefiles<-list.files(path = './output/mammals/tables/',pattern = 'all_realms_[0-9]*_quartilesDR_grid_residualstable.txt')
DR.pseudoposterior.residual.tables<-lapply(DR.pseudoposterior.residual.tablefiles,function(x) read.table(file=paste('./output/mammals/tables/',x,sep=''),sep='\t',header=T,stringsAsFactors = F))
DR.real.residual.tables<-read.table('./output/mammals/tables/all_realms_quartilesDR_grid_residualstable.txt',header=T,sep='\t',stringsAsFactors = F)
DR.pseudoposterior.Q1.correlation<-unlist(lapply(DR.pseudoposterior.residual.tables,function(x) cor.test(x$Q1.residuals,DR.real.residual.tables$Q1.residuals)$estimate))
DR.pseudoposterior.Q4.correlation<-unlist(lapply(DR.pseudoposterior.residual.tables,function(x) cor.test(x$Q4.residuals,DR.real.residual.tables$Q4.residuals)$estimate))
#this is part of Fig. S2
pdf('./output/mammals/plots/DR_vs_MCC_residualscorrelation.pdf')
hist(DR.pseudoposterior.Q1.correlation,xaxs='i',yaxs='i',xlim=c(0,1),ylab='',xlab='Pearsons r correlation',main='correlation MCC vs pseudoposterior cells specific residuals',breaks=50)
hist(DR.pseudoposterior.Q4.correlation,xaxs='i',yaxs='i',xlim=c(0,1),ylab='',xlab='Pearsons r correlation',main='correlation MCC vs pseudoposterior cells specific residuals',breaks=50)
dev.off()
#generate the DR (median of the pseudoposterior) quartile richness grid measures for mammals
source('./R/grid_functions.R')
generate_grid_richness_quartilesDR(dr.table='./output/mammals/tables/mammals_pseudoposteriorDRmedian.txt',species.grid.tablefile = './output/mammals/tables/100_all_realms_species_gridoccurrence_table.txt',cellgrid.tablefile='./output/mammals/tables/100_all_realms_richness_grid_table.txt',path='./output/mammals/tables/',name='all_realms_pseudoposterior')
#linear regressions of speciation rates (total richness vs richness in quartiles, etc)
#run linear model to predict quartile richness with total richness
build_lm_DRquartiles(dr.quartiles.gridtable='./output/mammals/tables/all_realms_pseudoposterior_quartilesDR_grid_table.txt',path='./output/mammals/tables/',name='all_realms_pseudoposterior')

####quartiles of DR for birds####
#generate the DR quartile richness grid measures for birds
source('./R/grid_functions.R')
generate_grid_richness_quartilesDR(dr.table='./output/birds/DR_birds_tree_IUCN.txt',species.grid.tablefile = './output/birds/tables/100_all_realms_species_gridoccurrence_table.txt',cellgrid.tablefile='./output/birds/tables/100_all_realms_richness_grid_table.txt',path='./output/birds/tables/',name='all_realms')
#linear regressions of speciation rates (total richness vs richness in quartiles, etc)
#run linear model to predict quartile richness with total richness
build_lm_DRquartiles(dr.quartiles.gridtable='./output/birds/tables/all_realms_quartilesDR_grid_table.txt',path='./output/birds/tables/',name='all_realms')
#####quartiles of DR for pseudoposterior
#save each table in DRpseudoposterior in a tablefile
DR.birds.pseudoposterior<-readRDS('./output/birds/trees/DR.birds.100.rds')
for(i in 1:100){
  write.table(DR.birds.pseudoposterior[[i]],file=paste('./output/birds/DR_birds_tree_IUCN_',i,'.txt',sep=''),sep='\t',quote=F,row.names=F)
}
source('./R/grid_functions.R')
for (i in 1:100){
  generate_grid_richness_quartilesDR(dr.table=paste('./output/birds/DR_birds_tree_IUCN_',i,'.txt',sep=''),species.grid.tablefile = './output/birds/tables/100_all_realms_species_gridoccurrence_table.txt',cellgrid.tablefile='./output/birds/tables/100_all_realms_richness_grid_table.txt',path='./output/birds/tables/',name=paste('all_realms_',i,sep=''))
  #linear regressions of speciation rates (total richness vs richness in quartiles, etc)
  #run linear model to predict quartile richness with total richness
  build_lm_DRquartiles(dr.quartiles.gridtable=paste('./output/birds/tables/all_realms_',i,'_quartilesDR_grid_table.txt',sep=''),path='./output/birds/tables/',name=paste('all_realms_',i,sep=''))
  
}
#read all DR.pseudoposterior residual files
DR.pseudoposterior.residual.tablefiles<-list.files(path = './output/birds/tables/',pattern = 'all_realms_[0-9]*_quartilesDR_grid_residualstable.txt')
DR.pseudoposterior.residual.tables<-lapply(DR.pseudoposterior.residual.tablefiles,function(x) read.table(file=paste('./output/birds/tables/',x,sep=''),sep='\t',header=T,stringsAsFactors = F))
DR.real.residual.tables<-read.table('./output/birds/tables/all_realms_quartilesDR_grid_residualstable.txt',header=T,sep='\t',stringsAsFactors = F)
DR.pseudoposterior.Q1.correlation<-unlist(lapply(DR.pseudoposterior.residual.tables,function(x) cor.test(x$Q1.residuals,DR.real.residual.tables$Q1.residuals)$estimate))
DR.pseudoposterior.Q4.correlation<-unlist(lapply(DR.pseudoposterior.residual.tables,function(x) cor.test(x$Q4.residuals,DR.real.residual.tables$Q4.residuals)$estimate))
#this is part of Fig. S2
pdf('./output/birds/plots/DR_vs_MCC_residualscorrelation.pdf')
hist(DR.pseudoposterior.Q1.correlation,xaxs='i',yaxs='i',xlim=c(0,1),ylab='',xlab='Pearsons r correlation',main='correlation MCC vs pseudoposterior cells specific residuals',breaks=50)
hist(DR.pseudoposterior.Q4.correlation,xaxs='i',yaxs='i',xlim=c(0,1),ylab='',xlab='Pearsons r correlation',main='correlation MCC vs pseudoposterior cells specific residuals',breaks=50)
dev.off()
#generate the DR (median of the pseudoposterior) quartile richness grid measures for birds
source('./R/grid_functions.R')
generate_grid_richness_quartilesDR(dr.table='./output/birds/tables/birds_pseudoposteriorDRmedian.txt',species.grid.tablefile = './output/birds/tables/100_all_realms_species_gridoccurrence_table.txt',cellgrid.tablefile='./output/birds/tables/100_all_realms_richness_grid_table.txt',path='./output/birds/tables/',name='all_realms_pseudoposterior')
#linear regressions of speciation rates (total richness vs richness in quartiles, etc)
#run linear model to predict quartile richness with total richness
build_lm_DRquartiles(dr.quartiles.gridtable='./output/birds/tables/all_realms_pseudoposterior_quartilesDR_grid_table.txt',path='./output/birds/tables/',name='all_realms_pseudoposterior')

####quartiles of BAMM for mammals####
source('./R/grid_functions.R')
generate_grid_richness_quartilesBAMM(BAMM.table='./output/mammals/tables/mammals_all_realms_BAMMTipRates.txt',species.grid.tablefile = './output/mammals/tables/100_all_realms_species_gridoccurrence_table.txt',cellgrid.tablefile='./output/mammals/tables/100_all_realms_richness_grid_table.txt',path='./output/mammals/tables/',name='mammals_all_realms')
#linear regressions of speciation rates (total richness vs richness in quartiles, etc)
#run linear model to predict quartile richness with total richness with lambda.avg
build_lm_BAMMquartiles(BAMM.quartiles.gridtable = './output/mammals/tables/mammals_all_realms_quartilesBAMMTipRates_grid_table.txt',variable.name = 'lambda.avg',path='./output/mammals/tables/',name='mammals_all_realms')
#run linear model to predict quartile richness with total richness with mu.avg
build_lm_BAMMquartiles(BAMM.quartiles.gridtable = './output/mammals/tables/mammals_all_realms_quartilesBAMMTipRates_grid_table.txt',variable.name = 'mu.avg',path='./output/mammals/tables/',name='mammals_all_realms')
#run linear model to predict quartile richness with total richness with netdiv.avg
build_lm_BAMMquartiles(BAMM.quartiles.gridtable = './output/mammals/tables/mammals_all_realms_quartilesBAMMTipRates_grid_table.txt',variable.name = 'netdiv.avg',path='./output/mammals/tables/',name='mammals_all_realms')

####quartiles of BAMM for birds####
source('./R/grid_functions.R')
generate_grid_richness_quartilesBAMM(BAMM.table='./output/birds/tables/birds_all_realms_BAMMTipRates.txt',species.grid.tablefile = './output/birds/tables/100_all_realms_species_gridoccurrence_table.txt',cellgrid.tablefile='./output/birds/tables/100_all_realms_richness_grid_table.txt',path='./output/birds/tables/',name='birds_all_realms')
#linear regressions of speciation rates (total richness vs richness in quartiles, etc)
#run linear model to predict quartile richness with total richness with lambda.avg
build_lm_BAMMquartiles(BAMM.quartiles.gridtable = './output/birds/tables/birds_all_realms_quartilesBAMMTipRates_grid_table.txt',variable.name = 'lambda.avg',path='./output/birds/tables/',name='birds_all_realms')
#run linear model to predict quartile richness with total richness with mu.avg
build_lm_BAMMquartiles(BAMM.quartiles.gridtable = './output/birds/tables/birds_all_realms_quartilesBAMMTipRates_grid_table.txt',variable.name = 'mu.avg',path='./output/birds/tables/',name='birds_all_realms')
#run linear model to predict quartile richness with total richness with netdiv.avg
build_lm_BAMMquartiles(BAMM.quartiles.gridtable = './output/birds/tables/birds_all_realms_quartilesBAMMTipRates_grid_table.txt',variable.name = 'netdiv.avg',path='./output/birds/tables/',name='birds_all_realms')

####correlation between DR and BAMM species-specific estimates####
#this is part of Fig. S3
source('./R/BAMM_functions.R')
pdf('./output/mammals/plots/mammals_BAMM_vs_DR_rates.pdf')
correlation_speciesDR_tipBAMM(DRtable.file ='./output/mammals/DR_mammals_tree_IUCN.txt', BAMMtable.file='./output/mammals/tables/mammals_all_realms_BAMMTipRates.txt')
dev.off()
pdf('./output/birds/plots/birds_BAMM_vs_DR_rates.pdf')
correlation_speciesDR_tipBAMM(DRtable.file ='./output/birds/DR_birds_tree_IUCN.txt', BAMMtable.file='./output/birds/tables/birds_all_realms_BAMMTipRates.txt')
dev.off()

####correlation between DR and BAMM residuals####
#this is part of Fig. S3
source('./R/grid_functions.R')
pdf('./output/mammals/plots/mammals_BAMMnetdivresiduals_vs_DRresiduals.pdf')
correlation_residualsDR_BAMM(DRresidualstable.file='./output/mammals/tables/all_realms_quartilesDR_grid_residualstable.txt',BAMMresidualstable.file='./output/mammals/tables/netdiv.avg_mammals_all_realms_quartilesBAMMTipRates_grid_residualstable.txt')
dev.off()
pdf('./output/birds/plots/birds_BAMMnetdivresiduals_vs_DRresiduals.pdf')
correlation_residualsDR_BAMM(DRresidualstable.file='./output/birds/tables/all_realms_quartilesDR_grid_residualstable.txt',BAMMresidualstable.file='./output/birds/tables/netdiv.avg_birds_all_realms_quartilesBAMMTipRates_grid_residualstable.txt')
dev.off()

####5) HISTORICAL BIOGEOGRAPHY ANALYSES####
####BioGeoBEARS run for mammals with realms####
#build a table for all species with ranges (#cells) in the world and in each realm
#also contains ranges inside predefined hotspot regions (with cells) in world and each realm
#results saved to 100_all_realms_ranges_plus_hotspots.txt
source('./R/BioGeoBEARS_wholetree.R')
build_world_range_hotspots_table(path='./output/mammals/tables/',world.table.file='./output/mammals/tables/100_all_realms_species_gridoccurrence_table.txt',weighted.endemism.file='./output/mammals/tables/100_all_realms_realms_richness_wend_grid_table.txt',quantile.hotspot=0.80,path.grids='./output/grids/')
#prepare BioGeoBEARS input for each realm
afrotrop7areas<-prepare_realm_input_7areas(world.table.file='./output/mammals/tables/100_all_realms_ranges_plus_hotspots.txt',treefile = './output/mammals/trees/mammals_tree_IUCN.tree',name='afrotrop',overlap=0.80)
austral7areas<-prepare_realm_input_7areas(world.table.file='./output/mammals/tables/100_all_realms_ranges_plus_hotspots.txt',treefile = './output/mammals/trees/mammals_tree_IUCN.tree',name='austral',overlap=0.80)
indo7areas<-prepare_realm_input_7areas(world.table.file='./output/mammals/tables/100_all_realms_ranges_plus_hotspots.txt',treefile = './output/mammals/trees/mammals_tree_IUCN.tree',name='indo',overlap=0.80)
nearctic7areas<-prepare_realm_input_7areas(world.table.file='./output/mammals/tables/100_all_realms_ranges_plus_hotspots.txt',treefile = './output/mammals/trees/mammals_tree_IUCN.tree',name='nearctic',overlap=0.80)
neotrop7areas<-prepare_realm_input_7areas(world.table.file='./output/mammals/tables/100_all_realms_ranges_plus_hotspots.txt',treefile = './output/mammals/trees/mammals_tree_IUCN.tree',name='neotrop',overlap=0.80)
palearctic7areas<-prepare_realm_input_7areas(world.table.file='./output/mammals/tables/100_all_realms_ranges_plus_hotspots.txt',treefile = './output/mammals/trees/mammals_tree_IUCN.tree',name='palearctic',overlap=0.80)
#move outputs to another folder
#create folder structure
dir.create('./output/mammals/immigration/')
dir.create('./output/mammals/immigration/whole_realms/')
realms<-c('afrotrop','austral','indo','nearctic','neotrop','palearctic')
#move geographyfile and trees
sapply(realms,function(x) dir.create(paste('./output/mammals/immigration/whole_realms/',x,sep='')))
sapply(realms,function(x) file.rename(from=list.files(pattern=paste(x,'.+','_geographyfile',sep=''))[1],to=paste('./output/mammals/immigration/whole_realms/',x,'/',x,'_7areas_inhotoutrealm_geographyfile.txt',sep='')))
sapply(realms,function(x) file.rename(from=list.files(pattern=paste(x,'*','_7areas_inhotoutrealm.tree',sep=''))[1],to=paste('./output/mammals/immigration/whole_realms/',x,'/',x,'_7areas_inhotoutrealm.tree',sep='')))
#run BioGeoBEARS plus BSM
#for 7 areas inhotoutrealm
setwd('~/Desktop/HOTSPOTS/repositories/hotspots_vertebrates/')
run_whole_realm_BioGeoBEARS_plusBSM_7areas_inhotoutrealm('afrotrop',path='./output/mammals/immigration/whole_realms/afrotrop/')
setwd('~/Desktop/HOTSPOTS/repositories/hotspots_vertebrates/')
run_whole_realm_BioGeoBEARS_plusBSM_7areas_inhotoutrealm('austral',path='./output/mammals/immigration/whole_realms/austral/')
setwd('~/Desktop/HOTSPOTS/repositories/hotspots_vertebrates/')
run_whole_realm_BioGeoBEARS_plusBSM_7areas_inhotoutrealm('indo',path='./output/mammals/immigration/whole_realms/indo/')
setwd('~/Desktop/HOTSPOTS/repositories/hotspots_vertebrates/')
run_whole_realm_BioGeoBEARS_plusBSM_7areas_inhotoutrealm('nearctic',path='./output/mammals/immigration/whole_realms/nearctic/')
setwd('~/Desktop/HOTSPOTS/repositories/hotspots_vertebrates/')
run_whole_realm_BioGeoBEARS_plusBSM_7areas_inhotoutrealm('neotrop',path='./output/mammals/immigration/whole_realms/neotrop/')
setwd('~/Desktop/HOTSPOTS/repositories/hotspots_vertebrates/')
run_whole_realm_BioGeoBEARS_plusBSM_7areas_inhotoutrealm('palearctic',path='./output/mammals/immigration/whole_realms/palearctic/')
setwd('~/Desktop/HOTSPOTS/repositories/hotspots_vertebrates/')

####BioGeoBEARS run for birds with realms####
#build a table for all species with ranges (#cells) in the world and in each realm
#also contains ranges inside predefined hotspot regions (with cells) in world and each realm
#results saved to 100_all_realms_ranges_plus_hotspots.txt
source('./R/BioGeoBEARS_wholetree.R')
build_world_range_hotspots_table(path='./output/birds/tables/',world.table.file='./output/birds/tables/100_all_realms_species_gridoccurrence_table.txt',weighted.endemism.file='./output/birds/tables/100_all_realms_realms_richness_wend_grid_table.txt',quantile.hotspot=0.80,path.grids='./output/grids/')
#prepare BioGeoBEARS input for each realm
afrotrop7areas<-prepare_realm_input_7areas(world.table.file='./output/birds/tables/100_all_realms_ranges_plus_hotspots.txt',treefile = './output/birds/trees/birds_tree_IUCN.tree',name='afrotrop',overlap=0.80)
austral7areas<-prepare_realm_input_7areas(world.table.file='./output/birds/tables/100_all_realms_ranges_plus_hotspots.txt',treefile = './output/birds/trees/birds_tree_IUCN.tree',name='austral',overlap=0.80)
indo7areas<-prepare_realm_input_7areas(world.table.file='./output/birds/tables/100_all_realms_ranges_plus_hotspots.txt',treefile = './output/birds/trees/birds_tree_IUCN.tree',name='indo',overlap=0.80)
nearctic7areas<-prepare_realm_input_7areas(world.table.file='./output/birds/tables/100_all_realms_ranges_plus_hotspots.txt',treefile = './output/birds/trees/birds_tree_IUCN.tree',name='nearctic',overlap=0.80)
neotrop7areas<-prepare_realm_input_7areas(world.table.file='./output/birds/tables/100_all_realms_ranges_plus_hotspots.txt',treefile = './output/birds/trees/birds_tree_IUCN.tree',name='neotrop',overlap=0.80)
palearctic7areas<-prepare_realm_input_7areas(world.table.file='./output/birds/tables/100_all_realms_ranges_plus_hotspots.txt',treefile = './output/birds/trees/birds_tree_IUCN.tree',name='palearctic',overlap=0.80)
#move outputs to another folder
#create folder structure
dir.create('./output/birds/immigration/')
dir.create('./output/birds/immigration/whole_realms/')
realms<-c('afrotrop','austral','indo','nearctic','neotrop','palearctic')
#move geographyfile and trees
sapply(realms,function(x) dir.create(paste('./output/birds/immigration/whole_realms/',x,sep='')))
sapply(realms,function(x) file.rename(from=list.files(pattern=paste(x,'.+','_geographyfile',sep=''))[1],to=paste('./output/birds/immigration/whole_realms/',x,'/',x,'_7areas_inhotoutrealm_geographyfile.txt',sep='')))
sapply(realms,function(x) file.rename(from=list.files(pattern=paste(x,'*','_7areas_inhotoutrealm.tree',sep=''))[1],to=paste('./output/birds/immigration/whole_realms/',x,'/',x,'_7areas_inhotoutrealm.tree',sep='')))
#run BioGeoBEARS plus BSM
#for 7 areas inhotoutrealm
setwd('~/Desktop/HOTSPOTS/repositories/hotspots_vertebrates/')
run_whole_realm_BioGeoBEARS_plusBSM_7areas_inhotoutrealm('afrotrop',path='./output/birds/immigration/whole_realms/afrotrop/')
setwd('~/Desktop/HOTSPOTS/repositories/hotspots_vertebrates/')
run_whole_realm_BioGeoBEARS_plusBSM_7areas_inhotoutrealm('austral',path='./output/birds/immigration/whole_realms/austral/')
setwd('~/Desktop/HOTSPOTS/repositories/hotspots_vertebrates/')
run_whole_realm_BioGeoBEARS_plusBSM_7areas_inhotoutrealm('indo',path='./output/birds/immigration/whole_realms/indo/')
setwd('~/Desktop/HOTSPOTS/repositories/hotspots_vertebrates/')
run_whole_realm_BioGeoBEARS_plusBSM_7areas_inhotoutrealm('nearctic',path='./output/birds/immigration/whole_realms/nearctic/')
setwd('~/Desktop/HOTSPOTS/repositories/hotspots_vertebrates/')
run_whole_realm_BioGeoBEARS_plusBSM_7areas_inhotoutrealm('neotrop',path='./output/birds/immigration/whole_realms/neotrop/')
setwd('~/Desktop/HOTSPOTS/repositories/hotspots_vertebrates/')
run_whole_realm_BioGeoBEARS_plusBSM_7areas_inhotoutrealm('palearctic',path='./output/birds/immigration/whole_realms/palearctic/')
setwd('~/Desktop/HOTSPOTS/repositories/hotspots_vertebrates/')

####simulate control polygons for mammals and run BioGeoBEARS####
source('./R/simulate_realmcontrolpolygons_allrealms.R')
dir.create('./output/mammals/immigration/controls/')
simulate_realms_control_polygons(hotspots.tablefile='./output/mammals/tables/100_all_realms_number.of.species.wend_0.8.txt',grid.file='./output/grids/grid_World_RealmsMerged_100.rds',nreplicates=50,path='./output/mammals/immigration/controls/')
###run the following on the cluster
#first move the /output/mammals/immmigration/controls folder (with points.neotropical.50.RDS, points.nearctic.50.RDS, etc) and /output/grids/*.rds and /output/mammals/trees/mammals_IUCN_tree.tree
#then run /scripts/conscriptoR /home/ji247/hotspots_vertebrates/prepare_BioGeoBEARS_point_input_cluster_allrealms_mammals.R
#then run /scripts/conscriptoR /home/ji247/hotspots_vertebrates/run_BioGeoBEARS_realm_controls_cluster.R 1' /home/ji247/hotspots_vertebrates/output/mammals/immigration/controls/' austral 37 42 (= request 1 cpu, run austral replicates from 37 to 42)
#***the commands are stored in ./R/run_controlpolygons_allrealms_cluster.txt

#this is for local processing of controls (see below for cluster processing, which saves some time)
#when they're over (can take up to 3 days) run rename_BSMRDSs_on_cluster_zip.R on node13 in the cluster (this renames all the BSM files and zips them by realm)
#then copy the zip back to local machine, and unzip in './output/mammals/immigration/controls/afrotrop/' etc (it will create a afrotrop.BSM folder that contains the RDS files)



####simulate control polygons for birds and run BioGeoBEARS####
source('./R/simulate_realmcontrolpolygons_allrealms.R')
dir.create('./output/birds/immigration/controls/')
simulate_realms_control_polygons(hotspots.tablefile='./output/birds/tables/100_all_realms_number.of.species.wend_0.8.txt',grid.file='./output/grids/grid_World_RealmsMerged_100.rds',nreplicates=50,path='./output/birds/immigration/controls/')
###run the following on the cluster
#first move the /output/birds/immmigration/controls folder (with points.neotropical.50.RDS, points.nearctic.50.RDS, etc) and /output/grids/*.rds and /output/birds/trees/birds_IUCN_tree.tree
#then run /scripts/conscriptoR /home/ji247/hotspots_vertebrates/run_BioGeoBEARS_realm_controls_cluster.R 1 '/home/ji247/hotspots_vertebrates/output/birds/immigration/controls/' austral 37 42 (= request 1 cpu, run austral replicates from 37 to 42)
#***the commands are stored in ./R/run_controlpolygons_allrealms_cluster.txt

#this is for local processing of controls (see below for cluster processing, which saves some time)
#when they're over (can take up to 10 days) run rename_BSMRDSs_on_cluster_zip.R on node13 in the cluster (this renames all the BSM files and zips them by realm)
#then copy the zip back to local machine, and unzip in './output/birds/immigration/controls/afrotrop/' etc (it will create a afrotrop.BSM folder that contains the RDS files)

#for cluster processing of controls:
#when mammals and birds runs are done, run ./R/cluster_afrotrop_processcontrolrates, ./R/cluster_austral_processcontrolrates.R, ./R/cluster_indo_processcontrolrates.R, ./R/cluster_nearctic_processcontrolrates.R, ./R/cluster_neotrop_processcontrolrates.R, ,./R/cluster_palearctic_processcontrolrates.R
#this will process the rates of controls and output palearctic_BSMcontrol_output_2my_full.RDS, etc files
#copy these to local to proceed (e.g. to ./output/mammals/immigration/controls/afrotrop/)

####6) MAIN TEXT FIGURES####
####Fig1_DRresiduals####
#this plots Fig. 1 for mammals with tropical realms first, then temperate realms
source('./R/linear_models_spatialautocorrelation.R')
table.dr.residuals.mammals<-prepare_linearmodel_table(vector.of.tables = c('./output/mammals/tables/all_realms_quartilesDR_grid_residualstable.txt','./output/mammals/tables/100_all_realms_number.of.species.wend_0.8.txt'))
realms.DRresiduals.sarlm.predicted.mammals<-predict_sarlm_table(table=table.dr.residuals.mammals,predictor.variable='hotspot',response.variables = c('Q1.residuals','Q2.residuals','Q3.residuals','Q4.residuals'),mode = 'realms')
source('./R/plot_boxplot_predictedtable_map.R')
pdf('./output/mammals/plots/mammals_DRresidualssarlm_diffmeansvioplots_new_troporder_WEhotspots.pdf')
lapply(c('Q1.residuals','Q4.residuals'),function(x)plot_propdifferencesmeans_plus_vioplots_troporder(table.object =realms.DRresiduals.sarlm.predicted.mammals,variable = x))
dev.off()
#this plots Fig. 1 for birds with tropical realms first, then temperate realms
source('./R/linear_models_spatialautocorrelation.R')
table.dr.residuals.birds<-prepare_linearmodel_table(vector.of.tables = c('./output/birds/tables/all_realms_quartilesDR_grid_residualstable.txt','./output/birds/tables/100_all_realms_number.of.species.wend_0.8.txt'))
realms.DRresiduals.sarlm.predicted.birds<-predict_sarlm_table(table=table.dr.residuals.birds,predictor.variable='hotspot',response.variables = c('Q1.residuals','Q2.residuals','Q3.residuals','Q4.residuals'),mode = 'realms')
source('./R/plot_boxplot_predictedtable_map.R')
pdf('./output/birds/plots/birds_DRresidualssarlm_diffmeansvioplots_new_troporder_WEhotspots.pdf')
lapply(c('Q1.residuals','Q4.residuals'),function(x)plot_propdifferencesmeans_plus_vioplots_troporder(table.object =realms.DRresiduals.sarlm.predicted.birds,variable = x))
dev.off()

####Figs2&3_average speciation and dispersal####
#this plots Fig. 2 & Fig 3 (the averages from 26 to 2 Ma) for mammals and birds
source('./R/generate_speciation_dispersal_plots_Fig2and3.R')
generate_dispersal_speciation_averageplots(pathmammals = './output/mammals/immigration/',pathbirds = './output/birds/immigration/',name = 'WEhotspots')
#this plots Fig. 2 & Fig 3 (the time bins from 26 to 2 Ma) for mammals and birds
generate_dispersal_speciation_throughtimeplots(pathmammals = './output/mammals/immigration/',pathbirds = './output/birds/immigration/',name = 'WEhotspots')

####Fig4_environmental differences####
#this plots Fig. 4 for mammals and birds
source('./R/linear_models_spatialautocorrelation.R')
source('./R/environmental_variables_sarlm.R')
###for mammals
#NPP
table.NPP.mammals<-prepare_linearmodel_table(vector.of.tables = c('./output/grids/tables/NPP_table.txt','./output/mammals/tables/100_all_realms_number.of.species.wend_0.8.txt'))
#elevation
table.hab.elevation.mammals<-prepare_linearmodel_table(vector.of.tables = c('./output/grids/tables/100_all_realms_elevation_table.txt','./output/grids/tables/100_allrealms_nhabitats_table.txt','./output/mammals/tables/100_all_realms_number.of.species.wend_0.8.txt'))
#past climate
table.CCV.mammals<-prepare_linearmodel_table(vector.of.tables = c('./output/grids/tables/AR.CCV_table.txt','./output/grids/tables/MAT.CCV_table.txt','./output/mammals/tables/100_all_realms_number.of.species.wend_0.8.txt'))
#tectonic movements over 65 Myr
table.tectonic.mammals<-prepare_linearmodel_table(vector.of.tables = c('./output/grids/tables/100_all_realms_tectonicmovement.txt','./output/mammals/tables/100_all_realms_number.of.species.wend_0.8.txt'))
variable.vector<-c('mean.NPP','mean.TRI','n.habitats','mean.AR.CCV','mean.MAT.CCV','tectonic.movement')
source('./R/heatmap_enviromentalvariables.R')
pdf('./output/mammals/plots/mammals_heatmap_allenvironmentscaledNPP_sarlm_SRhotspots_WEhotspots.pdf',width=10,height=7)
environmentalmodelsNPP_sarlm_variables_heatmap(table.NPP = table.NPP.mammals,table.hab.elevation = table.hab.elevation.mammals,table.CCV = table.CCV.mammals,table.tectonic=table.tectonic.mammals,variable.vector = variable.vector)
dev.off()

###for birds
#NPP
table.NPP.birds<-prepare_linearmodel_table(vector.of.tables = c('./output/grids/tables/NPP_table.txt','./output/birds/tables/100_all_realms_number.of.species.wend_0.8.txt'))
#elevation
table.hab.elevation.birds<-prepare_linearmodel_table(vector.of.tables = c('./output/grids/tables/100_all_realms_elevation_table.txt','./output/grids/tables/100_allrealms_nhabitats_table.txt','./output/birds/tables/100_all_realms_number.of.species.wend_0.8.txt'))
#past climate
table.CCV.birds<-prepare_linearmodel_table(vector.of.tables = c('./output/grids/tables/AR.CCV_table.txt','./output/grids/tables/MAT.CCV_table.txt','./output/birds/tables/100_all_realms_number.of.species.wend_0.8.txt'))
#tectonic movements over 65 Myr
table.tectonic.birds<-prepare_linearmodel_table(vector.of.tables = c('./output/grids/tables/100_all_realms_tectonicmovement.txt','./output/birds/tables/100_all_realms_number.of.species.wend_0.8.txt'))
variable.vector<-c('mean.NPP','mean.TRI','n.habitats','mean.AR.CCV','mean.MAT.CCV','tectonic.movement')
source('./R/heatmap_enviromentalvariables.R')
pdf('./output/birds/plots/birds_heatmap_allenvironmentscaledNPP_sarlm_SRhotspots_WEhotspots.pdf',width=10,height=7)
environmentalmodelsNPP_sarlm_variables_heatmap(table.NPP = table.NPP.birds,table.hab.elevation = table.hab.elevation.birds,table.CCV = table.CCV.birds,table.tectonic=table.tectonic.birds,variable.vector = variable.vector)
dev.off()

####7) ALTERNATIVE HOTSPOTS. SR-BASED####
#define hotspots as the top 20% cells with respect to endemism
source('./R/grid_functions.R')
#for mammals
define_hotspots(grid.tablefile ='./output/mammals/tables/100_all_realms_richness_grid_table.txt',quantile = 0.8,variable.name = 'number.of.species',path = './output/mammals/tables/',name = 'all_realms_SR')
#for birds
define_hotspots(grid.tablefile ='./output/birds/tables/100_all_realms_richness_grid_table.txt',quantile = 0.8,variable.name = 'number.of.species',path = './output/birds/tables/',name = 'all_realms_SR')

####8) ALTERNATIVE HOTSPOTS. NARROW RANGED SPECIES BASED####
#get centres of endemism (cells where narrow ranged - with range <=10 - species occur)
source('./R/grid_functions.R')
#for mammals
define_narrowrangedspecies_hotspots(species.ranges.tablefile = './output/mammals/tables/100_all_realms_species_gridoccurrence_table.txt',threshold=10,path='./output/mammals/tables/')
#for birds
define_narrowrangedspecies_hotspots(species.ranges.tablefile = './output/birds/tables/100_all_realms_species_gridoccurrence_table.txt',threshold=10,path='./output/birds/tables/')



####9) SUPPLEMENTARY ANALYSES####
####SRhotspots BioGeoBEARS run for mammals####

##immigration for mammals
#build a table for all species with ranges (#cells) in the world and in each realm
#also contains ranges inside predefined hotspot regions (with cells) in world and each realm
#results saved to 100_all_realms_ranges_plus_hotspots.txt
source('./R/BioGeoBEARS_wholetree.R')
build_world_range_hotspots_table_predefined_SR(path='./output/mammals/tables/',world.table.file='./output/mammals/tables/100_all_realms_species_gridoccurrence_table.txt',hotspots.file='./output/mammals/tables/100_all_realms_speciesrichness_0.8_hotspots.txt',quantile.hotspot=0.80,path.grids='./output/grids/')
#prepare BioGeoBEARS input for each realm
afrotrop7areas<-prepare_realm_input_7areas_SR(world.table.file='./output/mammals/tables/100_all_realms_ranges_plus_hotspots_SR.txt',treefile = './output/mammals/trees/mammals_tree_IUCN.tree',name='afrotrop',overlap=0.80)
austral7areas<-prepare_realm_input_7areas_SR(world.table.file='./output/mammals/tables/100_all_realms_ranges_plus_hotspots_SR.txt',treefile = './output/mammals/trees/mammals_tree_IUCN.tree',name='austral',overlap=0.80)
indo7areas<-prepare_realm_input_7areas_SR(world.table.file='./output/mammals/tables/100_all_realms_ranges_plus_hotspots_SR.txt',treefile = './output/mammals/trees/mammals_tree_IUCN.tree',name='indo',overlap=0.80)
#error in nearctic
nearctic7areas<-prepare_realm_input_7areas_SR(world.table.file='./output/mammals/tables/100_all_realms_ranges_plus_hotspots_SR.txt',treefile = './output/mammals/trees/mammals_tree_IUCN.tree',name='nearctic',overlap=0.80)
neotrop7areas<-prepare_realm_input_7areas_SR(world.table.file='./output/mammals/tables/100_all_realms_ranges_plus_hotspots_SR.txt',treefile = './output/mammals/trees/mammals_tree_IUCN.tree',name='neotrop',overlap=0.80)
palearctic7areas<-prepare_realm_input_7areas_SR(world.table.file='./output/mammals/tables/100_all_realms_ranges_plus_hotspots_SR.txt',treefile = './output/mammals/trees/mammals_tree_IUCN.tree',name='palearctic',overlap=0.80)

dir.create('./output/mammals/immigration_SR_new/')
dir.create('./output/mammals/immigration_SR_new/whole_realms/')
realms<-c('afrotrop','austral','indo','nearctic','neotrop','palearctic')
#move geographyfile and trees
sapply(realms,function(x) dir.create(paste('./output/mammals/immigration_SR_new/whole_realms/',x,sep='')))
sapply(realms,function(x) file.rename(from=list.files(pattern=paste(x,'.+','_geographyfile',sep=''))[1],to=paste('./output/mammals/immigration_SR_new/whole_realms/',x,'/',x,'_7areas_inhotoutrealm_geographyfile.txt',sep='')))
sapply(realms,function(x) file.rename(from=list.files(pattern=paste(x,'*','_7areas_inhotoutrealm.tree',sep=''))[1],to=paste('./output/mammals/immigration_SR_new/whole_realms/',x,'/',x,'_7areas_inhotoutrealm.tree',sep='')))

#run BioGeoBEARS plus BSM
#for 7 areas inhotoutrealm
setwd('/Users/javier/Documents/Work/HOTSPOTS/repositories/hotspots_vertebrates/')
run_whole_realm_BioGeoBEARS_plusBSM_7areas_inhotoutrealm('afrotrop',path='./output/mammals/immigration_SR_new/whole_realms/afrotrop/')
setwd('/Users/javier/Documents/Work/HOTSPOTS/repositories/hotspots_vertebrates/')
run_whole_realm_BioGeoBEARS_plusBSM_7areas_inhotoutrealm('austral',path='./output/mammals/immigration_SR_new/whole_realms/austral/')
setwd('/Users/javier/Documents/Work/HOTSPOTS/repositories/hotspots_vertebrates/')
run_whole_realm_BioGeoBEARS_plusBSM_7areas_inhotoutrealm('indo',path='./output/mammals/immigration_SR_new/whole_realms/indo/')
setwd('/Users/javier/Documents/Work/HOTSPOTS/repositories/hotspots_vertebrates/')
run_whole_realm_BioGeoBEARS_plusBSM_7areas_inhotoutrealm('nearctic',path='./output/mammals/immigration_SR_new/whole_realms/nearctic/')
setwd('/Users/javier/Documents/Work/HOTSPOTS/repositories/hotspots_vertebrates/')
run_whole_realm_BioGeoBEARS_plusBSM_7areas_inhotoutrealm('neotrop',path='./output/mammals/immigration_SR_new/whole_realms/neotrop/')
setwd('/Users/javier/Documents/Work/HOTSPOTS/repositories/hotspots_vertebrates/')
run_whole_realm_BioGeoBEARS_plusBSM_7areas_inhotoutrealm('palearctic',path='./output/mammals/immigration_SR_new/whole_realms/palearctic/')
setwd('/Users/javier/Documents/Work/HOTSPOTS/repositories/hotspots_vertebrates/')

##immigration for birds
#build a table for all species with ranges (#cells) in the world and in each realm
#also contains ranges inside predefined hotspot regions (with cells) in world and each realm
#results saved to 100_all_realms_ranges_plus_hotspots.txt
source('./R/BioGeoBEARS_wholetree.R')
build_world_range_hotspots_table_predefined_SR(path='./output/birds/tables/',world.table.file='./output/birds/tables/100_all_realms_species_gridoccurrence_table.txt',hotspots.file='./output/birds/tables/100_all_realms_speciesrichness_0.8_hotspots.txt',quantile.hotspot=0.80,path.grids='./output/grids/')
#prepare BioGeoBEARS input for each realm
afrotrop7areas<-prepare_realm_input_7areas_SR(world.table.file='./output/birds/tables/100_all_realms_ranges_plus_hotspots_SR.txt',treefile = './output/birds/trees/birds_tree_IUCN.tree',name='afrotrop',overlap=0.80)
austral7areas<-prepare_realm_input_7areas_SR(world.table.file='./output/birds/tables/100_all_realms_ranges_plus_hotspots_SR.txt',treefile = './output/birds/trees/birds_tree_IUCN.tree',name='austral',overlap=0.80)
indo7areas<-prepare_realm_input_7areas_SR(world.table.file='./output/birds/tables/100_all_realms_ranges_plus_hotspots_SR.txt',treefile = './output/birds/trees/birds_tree_IUCN.tree',name='indo',overlap=0.80)
#error in nearctic
nearctic7areas<-prepare_realm_input_7areas_SR(world.table.file='./output/birds/tables/100_all_realms_ranges_plus_hotspots_SR.txt',treefile = './output/birds/trees/birds_tree_IUCN.tree',name='nearctic',overlap=0.80)
neotrop7areas<-prepare_realm_input_7areas_SR(world.table.file='./output/birds/tables/100_all_realms_ranges_plus_hotspots_SR.txt',treefile = './output/birds/trees/birds_tree_IUCN.tree',name='neotrop',overlap=0.80)
palearctic7areas<-prepare_realm_input_7areas_SR(world.table.file='./output/birds/tables/100_all_realms_ranges_plus_hotspots_SR.txt',treefile = './output/birds/trees/birds_tree_IUCN.tree',name='palearctic',overlap=0.80)

dir.create('./output/birds/immigration_SR_new/')
dir.create('./output/birds/immigration_SR_new/whole_realms/')
realms<-c('afrotrop','austral','indo','nearctic','neotrop','palearctic')
#move geographyfile and trees
sapply(realms,function(x) dir.create(paste('./output/birds/immigration_SR_new/whole_realms/',x,sep='')))
sapply(realms,function(x) file.rename(from=list.files(pattern=paste(x,'.+','_geographyfile',sep=''))[1],to=paste('./output/birds/immigration_SR_new/whole_realms/',x,'/',x,'_7areas_inhotoutrealm_geographyfile.txt',sep='')))
sapply(realms,function(x) file.rename(from=list.files(pattern=paste(x,'*','_7areas_inhotoutrealm.tree',sep=''))[1],to=paste('./output/birds/immigration_SR_new/whole_realms/',x,'/',x,'_7areas_inhotoutrealm.tree',sep='')))

#run BioGeoBEARS plus BSM
#for 7 areas inhotoutrealm
setwd('/Users/javier/Documents/Work/HOTSPOTS/repositories/hotspots_vertebrates/')
run_whole_realm_BioGeoBEARS_plusBSM_7areas_inhotoutrealm('afrotrop',path='./output/birds/immigration_SR_new/whole_realms/afrotrop/')
setwd('/Users/javier/Documents/Work/HOTSPOTS/repositories/hotspots_vertebrates/')
run_whole_realm_BioGeoBEARS_plusBSM_7areas_inhotoutrealm('austral',path='./output/birds/immigration_SR_new/whole_realms/austral/')
setwd('/Users/javier/Documents/Work/HOTSPOTS/repositories/hotspots_vertebrates/')
run_whole_realm_BioGeoBEARS_plusBSM_7areas_inhotoutrealm('indo',path='./output/birds/immigration_SR_new/whole_realms/indo/')
setwd('/Users/javier/Documents/Work/HOTSPOTS/repositories/hotspots_vertebrates/')
run_whole_realm_BioGeoBEARS_plusBSM_7areas_inhotoutrealm('nearctic',path='./output/birds/immigration_SR_new/whole_realms/nearctic/')
setwd('/Users/javier/Documents/Work/HOTSPOTS/repositories/hotspots_vertebrates/')
run_whole_realm_BioGeoBEARS_plusBSM_7areas_inhotoutrealm('neotrop',path='./output/birds/immigration_SR_new/whole_realms/neotrop/')
setwd('/Users/javier/Documents/Work/HOTSPOTS/repositories/hotspots_vertebrates/')
run_whole_realm_BioGeoBEARS_plusBSM_7areas_inhotoutrealm('palearctic',path='./output/birds/immigration_SR_new/whole_realms/palearctic/')
setwd('/Users/javier/Documents/Work/HOTSPOTS/repositories/hotspots_vertebrates/')


####NRhotspots BioGeoBEARS run for mammals####
##immigration for mammals
#build a table for all species with ranges (#cells) in the world and in each realm
#also contains ranges inside predefined hotspot regions (with cells) in world and each realm
#results saved to 100_all_realms_ranges_plus_hotspots.txt
source('./R/BioGeoBEARS_wholetree.R')
build_world_range_hotspots_table_predefined_narrow(path='./output/mammals/tables/',world.table.file='./output/mammals/tables/100_all_realms_species_gridoccurrence_table.txt',hotspots.file='./output/mammals/tables/100_all_realms_narrow.ranged.species_hotspots.txt',quantile.hotspot=0.80,path.grids='./output/grids/')
#prepare BioGeoBEARS input for each realm
afrotrop7areas<-prepare_realm_input_7areas_narrow(world.table.file='./output/mammals/tables/100_all_realms_ranges_plus_hotspots_narrowrangesp.txt',treefile = './output/mammals/trees/mammals_tree_IUCN.tree',name='afrotrop',overlap=0.80)
austral7areas<-prepare_realm_input_7areas_narrow(world.table.file='./output/mammals/tables/100_all_realms_ranges_plus_hotspots_narrowrangesp.txt',treefile = './output/mammals/trees/mammals_tree_IUCN.tree',name='austral',overlap=0.80)
indo7areas<-prepare_realm_input_7areas_narrow(world.table.file='./output/mammals/tables/100_all_realms_ranges_plus_hotspots_narrowrangesp.txt',treefile = './output/mammals/trees/mammals_tree_IUCN.tree',name='indo',overlap=0.80)
nearctic7areas<-prepare_realm_input_7areas_narrow(world.table.file='./output/mammals/tables/100_all_realms_ranges_plus_hotspots_narrowrangesp.txt',treefile = './output/mammals/trees/mammals_tree_IUCN.tree',name='nearctic',overlap=0.80)
neotrop7areas<-prepare_realm_input_7areas_narrow(world.table.file='./output/mammals/tables/100_all_realms_ranges_plus_hotspots_narrowrangesp.txt',treefile = './output/mammals/trees/mammals_tree_IUCN.tree',name='neotrop',overlap=0.80)
palearctic7areas<-prepare_realm_input_7areas_narrow(world.table.file='./output/mammals/tables/100_all_realms_ranges_plus_hotspots_narrowrangesp.txt',treefile = './output/mammals/trees/mammals_tree_IUCN.tree',name='palearctic',overlap=0.80)
dir.create('./output/mammals/immigration_narrow/')
dir.create('./output/mammals/immigration_narrow/whole_realms/')
realms<-c('afrotrop','austral','indo','nearctic','neotrop','palearctic')
#move geographyfile and trees
sapply(realms,function(x) dir.create(paste('./output/mammals/immigration_narrow/whole_realms/',x,sep='')))
sapply(realms,function(x) file.rename(from=list.files(pattern=paste(x,'.+','_geographyfile',sep=''))[1],to=paste('./output/mammals/immigration_narrow/whole_realms/',x,'/',x,'_7areas_inhotoutrealm_geographyfile.txt',sep='')))
sapply(realms,function(x) file.rename(from=list.files(pattern=paste(x,'*','_7areas_inhotoutrealm.tree',sep=''))[1],to=paste('./output/mammals/immigration_narrow/whole_realms/',x,'/',x,'_7areas_inhotoutrealm.tree',sep='')))

#run BioGeoBEARS plus BSM
#for 7 areas inhotoutrealm
setwd('/Users/javier/Documents/Work/HOTSPOTS/repositories/hotspots_vertebrates/')
run_whole_realm_BioGeoBEARS_plusBSM_7areas_inhotoutrealm('afrotrop',path='./output/mammals/immigration_narrow/whole_realms/afrotrop/')
setwd('/Users/javier/Documents/Work/HOTSPOTS/repositories/hotspots_vertebrates/')
run_whole_realm_BioGeoBEARS_plusBSM_7areas_inhotoutrealm('austral',path='./output/mammals/immigration_narrow/whole_realms/austral/')
setwd('/Users/javier/Documents/Work/HOTSPOTS/repositories/hotspots_vertebrates/')
run_whole_realm_BioGeoBEARS_plusBSM_7areas_inhotoutrealm('indo',path='./output/mammals/immigration_narrow/whole_realms/indo/')
setwd('/Users/javier/Documents/Work/HOTSPOTS/repositories/hotspots_vertebrates/')
run_whole_realm_BioGeoBEARS_plusBSM_7areas_inhotoutrealm('nearctic',path='./output/mammals/immigration_narrow/whole_realms/nearctic/')
setwd('/Users/javier/Documents/Work/HOTSPOTS/repositories/hotspots_vertebrates/')
run_whole_realm_BioGeoBEARS_plusBSM_7areas_inhotoutrealm('neotrop',path='./output/mammals/immigration_narrow/whole_realms/neotrop/')
setwd('/Users/javier/Documents/Work/HOTSPOTS/repositories/hotspots_vertebrates/')
run_whole_realm_BioGeoBEARS_plusBSM_7areas_inhotoutrealm('palearctic',path='./output/mammals/immigration_narrow/whole_realms/palearctic/')
setwd('/Users/javier/Documents/Work/HOTSPOTS/repositories/hotspots_vertebrates/')

####NRhotspots BioGeoBEARS run for birds####
##immigration for birds
#build a table for all species with ranges (#cells) in the world and in each realm
#also contains ranges inside predefined hotspot regions (with cells) in world and each realm
#results saved to 100_all_realms_ranges_plus_hotspots.txt
source('./R/BioGeoBEARS_wholetree.R')
build_world_range_hotspots_table_predefined_narrow(path='./output/birds/tables/',world.table.file='./output/birds/tables/100_all_realms_species_gridoccurrence_table.txt',hotspots.file='./output/birds/tables/100_all_realms_narrow.ranged.species_hotspots.txt',quantile.hotspot=0.80,path.grids='./output/grids/')
#prepare BioGeoBEARS input for each realm
afrotrop7areas<-prepare_realm_input_7areas_narrow(world.table.file='./output/birds/tables/100_all_realms_ranges_plus_hotspots_narrowrangesp.txt',treefile = './output/birds/trees/birds_tree_IUCN.tree',name='afrotrop',overlap=0.80)
austral7areas<-prepare_realm_input_7areas_narrow(world.table.file='./output/birds/tables/100_all_realms_ranges_plus_hotspots_narrowrangesp.txt',treefile = './output/birds/trees/birds_tree_IUCN.tree',name='austral',overlap=0.80)
indo7areas<-prepare_realm_input_7areas_narrow(world.table.file='./output/birds/tables/100_all_realms_ranges_plus_hotspots_narrowrangesp.txt',treefile = './output/birds/trees/birds_tree_IUCN.tree',name='indo',overlap=0.80)
nearctic7areas<-prepare_realm_input_7areas_narrow(world.table.file='./output/birds/tables/100_all_realms_ranges_plus_hotspots_narrowrangesp.txt',treefile = './output/birds/trees/birds_tree_IUCN.tree',name='nearctic',overlap=0.80)
neotrop7areas<-prepare_realm_input_7areas_narrow(world.table.file='./output/birds/tables/100_all_realms_ranges_plus_hotspots_narrowrangesp.txt',treefile = './output/birds/trees/birds_tree_IUCN.tree',name='neotrop',overlap=0.80)
palearctic7areas<-prepare_realm_input_7areas_narrow(world.table.file='./output/birds/tables/100_all_realms_ranges_plus_hotspots_narrowrangesp.txt',treefile = './output/birds/trees/birds_tree_IUCN.tree',name='palearctic',overlap=0.80)
dir.create('./output/birds/immigration_narrow/')
dir.create('./output/birds/immigration_narrow/whole_realms/')
realms<-c('afrotrop','austral','indo','nearctic','neotrop','palearctic')
#move geographyfile and trees
sapply(realms,function(x) dir.create(paste('./output/birds/immigration_narrow/whole_realms/',x,sep='')))
sapply(realms,function(x) file.rename(from=list.files(pattern=paste(x,'.+','_geographyfile',sep=''))[1],to=paste('./output/birds/immigration_narrow/whole_realms/',x,'/',x,'_7areas_inhotoutrealm_geographyfile.txt',sep='')))
sapply(realms,function(x) file.rename(from=list.files(pattern=paste(x,'*','_7areas_inhotoutrealm.tree',sep=''))[1],to=paste('./output/birds/immigration_narrow/whole_realms/',x,'/',x,'_7areas_inhotoutrealm.tree',sep='')))

#run BioGeoBEARS plus BSM
#for 7 areas inhotoutrealm
setwd('/Users/javier/Documents/Work/HOTSPOTS/repositories/hotspots_vertebrates/')
run_whole_realm_BioGeoBEARS_plusBSM_7areas_inhotoutrealm('afrotrop',path='./output/birds/immigration_narrow/whole_realms/afrotrop/')
setwd('/Users/javier/Documents/Work/HOTSPOTS/repositories/hotspots_vertebrates/')
run_whole_realm_BioGeoBEARS_plusBSM_7areas_inhotoutrealm('austral',path='./output/birds/immigration_narrow/whole_realms/austral/')
setwd('/Users/javier/Documents/Work/HOTSPOTS/repositories/hotspots_vertebrates/')
run_whole_realm_BioGeoBEARS_plusBSM_7areas_inhotoutrealm('indo',path='./output/birds/immigration_narrow/whole_realms/indo/')
setwd('/Users/javier/Documents/Work/HOTSPOTS/repositories/hotspots_vertebrates/')
run_whole_realm_BioGeoBEARS_plusBSM_7areas_inhotoutrealm('nearctic',path='./output/birds/immigration_narrow/whole_realms/nearctic/')
setwd('/Users/javier/Documents/Work/HOTSPOTS/repositories/hotspots_vertebrates/')
run_whole_realm_BioGeoBEARS_plusBSM_7areas_inhotoutrealm('neotrop',path='./output/birds/immigration_narrow/whole_realms/neotrop/')
setwd('/Users/javier/Documents/Work/HOTSPOTS/repositories/hotspots_vertebrates/')
run_whole_realm_BioGeoBEARS_plusBSM_7areas_inhotoutrealm('palearctic',path='./output/birds/immigration_narrow/whole_realms/palearctic/')
setwd('/Users/javier/Documents/Work/HOTSPOTS/repositories/hotspots_vertebrates/')

####10) SUPPLEMENTARY FIGURES
####FigS1_mapofWEhotspots####
#plot hotspots in world map for mammals
source('./R/grid_plots.R')
pdf('./output/mammals/plots/mammals_hotspotWE0.8_map_realms.pdf')
plot_hotspot_map_realms(worldfile='./output/realm.merge.proj.RDS',gridfile='./output/mammals/grid_World_RealmsMerged_100.rds',tablefile = './output/mammals/tables/100_all_realms_number.of.species.wend_0.8.txt',plot.name = 'hotspots_mammals_WE_0.8_realms')
dev.off()
#plot hotspots in world map for birds
source('./R/grid_plots.R')
pdf('./output/birds/plots/birds_hotspotWE0.8_map_realms.pdf')
plot_hotspot_map_realms(worldfile='./output/realm.merge.proj.RDS',gridfile='./output/birds/grid_World_RealmsMerged_100.rds',tablefile = './output/birds/tables/100_all_realms_number.of.species.wend_0.8.txt',plot.name = 'hotspots_birds_WE_0.8_realms')
dev.off()

####FigS4_colonisationtimes for WEhotspots####
generate_colonisation_times(pathmammals = './output/mammals/immigration/',pathbirds = './output/birds/immigration/',name = 'WEhotspots')

####FigS5 to S9 (control polygons rate comparisons vs real rates)
source('./R/plots_rates_real_vs_controls_mammalsbirds.R')

#outputs will be in './output/mammals/plots'
####FigS10_contiguity analysis####
source('./R/averagedistanceneighboursclass.R')
#for mammals
get_averagedistanceneighboursclass(hotspots.table.file = './output/mammals/tables/100_all_realms_number.of.species.wend_0.8.txt',outputpath='./output/mammals/plots/',name='mammals')
#for birds
get_averagedistanceneighboursclass(hotspots.table.file = './output/birds/tables/100_all_realms_number.of.species.wend_0.8.txt',outputpath='./output/birds/plots/',name='birds')

####FigS11_mapofSRhotspots####
#plot hotspots in world map for mammals
source('./R/grid_plots.R')
pdf('./output/mammals/plots/mammals_hotspotSR0.8_map_realms.pdf')
plot_hotspot_map_realms(worldfile='./output/realm.merge.proj.RDS',gridfile='./output/grids/grid_World_RealmsMerged_100.rds',tablefile = './output/mammals/tables/100_all_realms_SR_number.of.species_0.8.txt',plot.name = 'hotspots_mammals_SR_0.8_realms')
dev.off()
#plot hotspots in world map for birds
source('./R/grid_plots.R')
pdf('./output/birds/plots/birds_hotspotSR0.8_map_realms.pdf')
plot_hotspot_map_realms(worldfile='./output/realm.merge.proj.RDS',gridfile='./output/grids/grid_World_RealmsMerged_100.rds',tablefile = './output/birds/tables/100_all_realms_SR_number.of.species_0.8.txt',plot.name = 'hotspots_birds_SR_0.8_realms')
dev.off()

####FigS12_mapofNRhotspots####
#plot hotspots in world map for mammals
source('./R/grid_plots.R')
pdf('./output/mammals/plots/mammals_hotspotsNR_map_realms.pdf')
plot_hotspot_map_realms(worldfile='./output/realm.merge.proj.RDS',gridfile='./output/grids/grid_World_RealmsMerged_100.rds',tablefile = './output/mammals/tables/100_all_realms_narrow.ranged.species_hotspots.txt',plot.name = 'hotspots_mammals_NR_realms')
dev.off()
#plot hotspots in world map for birds
source('./R/grid_plots.R')
pdf('./output/birds/plots/birds_hotspotsNR_map_realms.pdf')
plot_hotspot_map_realms(worldfile='./output/realm.merge.proj.RDS',gridfile='./output/grids/grid_World_RealmsMerged_100.rds',tablefile = './output/birds/tables/100_all_realms_narrow.ranged.species_hotspots.txt',plot.name = 'hotspots_birds_NR_realms')
dev.off()

#####FigS13_DRresiduals_SRhotspots####
#for mammals
source('./R/linear_models_spatialautocorrelation.R')
table.dr.residuals.mammals<-prepare_linearmodel_table(vector.of.tables = c('./output/mammals/tables/all_realms_quartilesDR_grid_residualstable.txt','./output/mammals/tables/100_all_realms_SR_number.of.species_0.8.txt'))
realms.DRresiduals.sarlm.predicted.mammals<-predict_sarlm_table(table=table.dr.residuals.mammals,predictor.variable='hotspot',response.variables = c('Q1.residuals','Q2.residuals','Q3.residuals','Q4.residuals'),mode = 'realms')
source('./R/plot_boxplot_predictedtable_map.R')
source('./R/plot_boxplot_predictedtable_map.R')
pdf('./output/mammals/plots/mammals_DRresidualssarlm_diffmeansvioplots_new_alphaorder_SRhotspots.pdf')
lapply(c('Q1.residuals','Q4.residuals'),function(x)plot_propdifferencesmeans_plus_vioplots_alphaorder(table.object =realms.DRresiduals.sarlm.predicted.mammals,variable = x))
dev.off()
#for birds
source('./R/linear_models_spatialautocorrelation.R')
table.dr.residuals.birds<-prepare_linearmodel_table(vector.of.tables = c('./output/birds/tables/all_realms_quartilesDR_grid_residualstable.txt','./output/birds/tables/100_all_realms_SR_number.of.species_0.8.txt'))
realms.DRresiduals.sarlm.predicted.birds<-predict_sarlm_table(table=table.dr.residuals.birds,predictor.variable='hotspot',response.variables = c('Q1.residuals','Q2.residuals','Q3.residuals','Q4.residuals'),mode = 'realms')
source('./R/plot_boxplot_predictedtable_map.R')
source('./R/plot_boxplot_predictedtable_map.R')
pdf('./output/birds/plots/birds_DRresidualssarlm_diffmeansvioplots_new_alphaorder_SRhotspots.pdf')
lapply(c('Q1.residuals','Q4.residuals'),function(x)plot_propdifferencesmeans_plus_vioplots_alphaorder(table.object =realms.DRresiduals.sarlm.predicted.birds,variable = x))
dev.off()

####FigS14_SRhotspots_speciationdispersal####
source('./R/generate_speciation_dispersal_plots_Fig2and3.R')
generate_dispersal_speciation_averageplots_NRhot(pathmammals = './output/mammals/immigration_SR_new/',pathbirds = './output/birds/immigration_SR_new/',name = 'SRhotspots')
#this plots Fig. 2 & Fig 3 (the time bins from 26 to 2 Ma) for mammals and birds
generate_dispersal_speciation_averageplots_NRhot(pathmammals = './output/mammals/immigration_SR_new/',pathbirds = './output/birds/immigration_SR_new/',name = 'SRhotspots')


####FigS15_DRresiduals_NRhotspots####
#for mammals
source('./R/linear_models_spatialautocorrelation.R')
table.dr.residuals.mammals<-prepare_linearmodel_table(vector.of.tables = c('./output/mammals/tables/all_realms_quartilesDR_grid_residualstable.txt','./output/mammals/tables/100_all_realms_narrow.ranged.species_hotspots.txt'))
realms.DRresiduals.sarlm.predicted.mammals<-predict_sarlm_table(table=table.dr.residuals.mammals,predictor.variable='hotspot',response.variables = c('Q1.residuals','Q2.residuals','Q3.residuals','Q4.residuals'),mode = 'realms')
source('./R/plot_boxplot_predictedtable_map.R')
pdf('./output/mammals/plots/mammals_DRresidualssarlm_diffmeansvioplots_new_alphaorder_narrowrangedhotspots.pdf')
lapply(c('Q1.residuals','Q4.residuals'),function(x)plot_propdifferencesmeans_plus_vioplots_alphaorder(table.object =realms.DRresiduals.sarlm.predicted.mammals,variable = x))
dev.off()
#for birds
source('./R/linear_models_spatialautocorrelation.R')
table.dr.residuals.birds<-prepare_linearmodel_table(vector.of.tables = c('./output/birds/tables/all_realms_quartilesDR_grid_residualstable.txt','./output/birds/tables/100_all_realms_narrow.ranged.species_hotspots.txt'))
realms.DRresiduals.sarlm.predicted.birds<-predict_sarlm_table(table=table.dr.residuals.birds,predictor.variable='hotspot',response.variables = c('Q1.residuals','Q2.residuals','Q3.residuals','Q4.residuals'),mode = 'realms')
source('./R/plot_boxplot_predictedtable_map.R')
pdf('./output/birds/plots/birds_DRresidualssarlm_diffmeansvioplots_new_alphaorder_narrowrangedhotspots.pdf')
lapply(c('Q1.residuals','Q4.residuals'),function(x)plot_propdifferencesmeans_plus_vioplots_alphaorder(table.object =realms.DRresiduals.sarlm.predicted.birds,variable = x))
dev.off()

####FigS16_NRhotspots_speciationdispersal####
source('./R/generate_speciation_dispersal_plots_Fig2and3.R')
generate_dispersal_speciation_averageplots_NRhot(pathmammals = './output/mammals/immigration_narrow/',pathbirds = './output/birds/immigration_narrow/',name = 'NRhotspots')
#this plots Fig. 2 & Fig 3 (the time bins from 26 to 2 Ma) for mammals and birds
generate_dispersal_speciation_averageplots_NRhot(pathmammals = './output/mammals/immigration_narrow/',pathbirds = './output/birds/immigration_narrow/',name = 'NRhotspots')

####FigS17_map_Myershotspots####
#plot hotspots in world map for mammals
source('./R/grid_plots.R')
pdf('./output/mammals/plots/Myers_hotspots.pdf')
plot_Myershotspots_map_realms(worldfile='./output/realm.merge.proj.RDS',gridfile='./output/grids/grid_World_RealmsMerged_100.rds',tablefile = './output/grids/tables/100_all_realms_Myershotspots_grid_table.txt',plot.name = 'Myershotspots')
dev.off()


####FigS18_map_controlhotspots_Afrotropics####
source('./R/grid_plots.R')
plot_Afrotropics_mammalsWEhotspots_pluscontrolpolygons_map(hotspots.table.file = './output/mammals/tables/100_all_realms_number.of.species.wend_0.8.txt')
