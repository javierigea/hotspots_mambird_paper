DR.mammals
DR.mammals.100
identical(as.character(DR.mammals$Species),as.character(DR.mammals.100[[2]]$Species))
species.list<-lapply(DR.mammals.100,function(x)x$Species)
#order DR.mammals with DR.mammals.100[[1]]
DR.mammals<-DR.mammals[match(DR.mammals.100[[1]]$Species,DR.mammals$Species),]
DR.all.df<-as.data.frame(matrix(NA,nrow=nrow(DR.mammals),ncol=length(DR.mammals.100)+2),stringAsFactors=F)
DR.all.df[,1]<-DR.mammals$Species
DR.all.df[,2]<-DR.mammals$DR

colnames(DR.all.df)[c(1,2)]<-c('Species','DR.MCC')
for (i in 1:nrow(DR.all.df)){
  cat(i,'\n')
  DR.all.df[i,c(3:ncol(DR.all.df))]<-unlist(lapply(DR.mammals.100,function(x) x[i,2]))
}
quantile.DR<-apply(DR.all.df,1,function(x)ecdf(x[c(3:length(x))])(x[2]))
hist(quantile.DR)
median.DR.pseudoposterior<-apply(DR.all.df,1,function(x)median(x[c(3:length(x))]))
plot(as.numeric(DR.all.df[,2]),as.numeric(median.DR.pseudoposterior))
cor.test(as.numeric(DR.all.df[,2]),as.numeric(median.DR.pseudoposterior))


identical(as.character(DR.birds$Species),as.character(DR.birds.100[[2]]$Species))
species.list<-lapply(DR.birds.100,function(x)x$Species)
length(species.list)
length(unique(species.list))
#order DR.mammals with DR.mammals.100[[1]]
DR.birds<-DR.birds[match(DR.birds.100[[1]]$Species,DR.birds$Species),]
identical(DR.birds$Species,DR.birds.100[[2]]$Species)
DR.all.df<-as.data.frame(matrix(NA,nrow=nrow(DR.birds),ncol=length(DR.birds.100)+2),stringAsFactors=F)
DR.all.df[,1]<-DR.birds$Species
DR.all.df[,2]<-DR.birds$DR

colnames(DR.all.df)[c(1,2)]<-c('Species','DR.MCC')
for (i in 1:nrow(DR.all.df)){
  cat(i,'\n')
  DR.all.df[i,c(3:ncol(DR.all.df))]<-as.numeric(unlist(lapply(DR.birds.100,function(x) x[i,2])))
}
quantile.DR<-apply(DR.all.df,1,function(x)ecdf(x[c(3:length(x))])(x[2]))
hist(quantile.DR)
classes<-apply(DR.all.df,2,function(x)class(x))
median.DR.pseudoposterior<-apply(DR.all.df,1,function(x) median(as.numeric(x[c(3:length(x))])))

plot(as.numeric(DR.all.df[,2]),as.numeric(median.DR.pseudoposterior))
cor.test(as.numeric(DR.all.df[,2]),as.numeric(median.DR.pseudoposterior))
