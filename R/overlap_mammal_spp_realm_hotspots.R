library(rgdal)
library(sp)
library(rgeos)
library(maptools)
library(parallel)
library(PBSmapping)
library(ape)
library(geiger)
library(phytools)
library(phangorn)
library(diversitree)

source('./R/process_realm_hotspot.R')
ncores<-4
#realm and hotspots names
afr<-'afrotropical'
afr.hot<-c('mediterranean','guinea','horn','afromontane','coasteastafrica','maputaland','karoo','cape')
ara<-'arabian'
ara.hot<-c('mediterranean','iranoanatolia','himalaya','horn','afromontane')
aus<-'australian'
aus.hot<-c('eastaustralia','swaustralia','newzealand')
nea<-'nearctic'
nea.hot<-c('california','madrean','mesoamerica','caribbean')
neo<-'neotropical'
neo.hot<-c("mesoamerica","caribbean","andes","tumbes","atlanticforest","cerrado","valdivia")
oce<-'oceanina'
oce.hot<-c('pm','emi','newcaledonia','newzealand','wallacea','philippines')
ori<-'oriental'
ori.hot<-c('sundaland','wallacea','westernghats','indoburma','himalaya','mtsswchina','iranoanatolia','philippines')
pal<-'palearctic'
pal.hot<-c('caucasus','iranoanatolia','mediterranean','centralasia')
sin<-'sinojapanese'
sin.hot<-c('japan','mttswchina','himalaya','centralasia')

#read mammal spp data
spp.ranges<-readOGR("../raw_data/GIS/TERRESTRIAL_MAMMALS/", 'TERRESTRIAL_MAMMALS')
#get mammal spp names
spp<-unique(spp.ranges@data$binomial)
#separate by spp name
spp.ranges.separated<-mclapply(mc.cores=ncores,spp, function(x) try(spp.ranges[spp.ranges@data$binomial==x,]))
#merge spp polygons
spp.ranges.separated<-mclapply(mc.cores=ncores,spp.ranges.separated,function (x) try(gUnaryUnion(x)))
#fix polygons
spp.ranges.separated<-mclapply(mc.cores=ncores,spp.ranges.separated, function(x) try(gBuffer(x, width=0)))
#check they're valid
valid<-mclapply(mc.cores=ncores,spp.ranges.separated, function(x) try(gIsValid(x)))
length(valid[valid==F])

#get overlaps of species for each hotspot within each realm and corresponding geosse_input_files
process_realm_hotspot_mammals(afr,afr.hot)
process_realm_hotspot_mammals(ara,ara.hot)
process_realm_hotspot_mammals(aus,aus.hot)
process_realm_hotspot_mammals(nea,nea.hot)
process_realm_hotspot_mammals(neo,neo.hot)
process_realm_hotspot_mammals(oce,oce.hot)
process_realm_hotspot_mammals(ori,ori.hot)
process_realm_hotspot_mammals(pal,pal.hot)
process_realm_hotspot_mammals(sin,sin.hot)
