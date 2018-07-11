#cd to './output/mammals/immigration/' or './output/birds/immigration/'
#BSM.files<-list.files(path='./controls/',pattern='_BSM',recursive=T)
BSM.files<-list.files(path='./',pattern='_BSM',recursive=T)
#BSM.newfiles<-sapply(BSM.files,function(x)sub(x,pattern='_BSM',replacement=paste('_',unlist(strsplit(x,split='/'))[2],'_BSM',sep='')))
BSM.newfiles<-sapply(BSM.files,function(x)sub(x,pattern='_BSM',replacement=paste('_',unlist(strsplit(x,split='/'))[1],'_BSM',sep='')))
BSM.newfiles<-sapply(BSM.newfiles,function(x)sub(x,pattern='.+/',replacement=''))
for(i in 1:length(BSM.files)){
  file.rename(from=paste('./',BSM.files[i],sep=''),to=BSM.newfiles[i])  
  
}

#for(i in 1:length(BSM.files)){
#  file.rename(from=paste('./controls/',BSM.files[i],sep=''),to=BSM.newfiles[i])  
#  
#}
system('zip afrotrop.BSM.zip ./afrotrop*_BSM.RDS')
system('zip austral.BSM.zip ./austral*_BSM.RDS')
system('zip indo.BSM.zip ./indo*_BSM.RDS')
system('zip nearctic.BSM.zip ./nearctic*_BSM.RDS')
system('zip neotrop.BSM.zip ./neotrop*_BSM.RDS')
system('zip palearctic.BSM.zip ./palearctic*_BSM.RDS')


for(i in 1:length(BSM.files)){
  file.rename(from=paste('./',BSM.files[i],sep=''),to=BSM.newfiles[i])  
  
}