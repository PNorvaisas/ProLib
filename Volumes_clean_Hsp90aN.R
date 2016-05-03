#To install missing package just run:
#install.packages('package')
#It is advised to install R in user directory, writing to which does not require sudo privileges
library('ggplot2')
library(ggdendro)
library('gplots')
library('plyr')
library('plot3D')
library('reshape2')
library(tidyr)
library(rgl)
library(cluster)
library(Cairo)

#Read summary of volumes
summr<-read.table('Prepared/Summary.csv',sep=',',quote = '"',header = TRUE)
summr$File<-toupper(summr$File)
summr[,c('File','Chain')]<- transform(summr, File = colsplit(File, "_", names = c('a', 'b')))$File
summr<-summr[,colnames(summ)[c(1,13,2:12)]]

#Read inof from RSCB
pdbinfo<-read.table('Hsp90aN_info.csv',sep=',',quote = '"',header = TRUE)
pdbinfo<-pdbinfo[,c('PDB.ID','Resolution','Structure.MW','Residue.Count')]



#Read summary of PDB info
pdb<-read.table('Summary_sequences.csv',sep=',',quote = '"',header = TRUE)
pdb$PDB.ID<-toupper(pdb$PDB.ID)
pdb<-pdb[,1:7]
pdb[pdb=='']<-NA

#Structures without ligands
subset(pdb,is.na(Heteroatoms))

#Read affinties
#aff<-read.table('CAII_affinities_fixed.csv',sep=',',quote = '"',header = TRUE)
aff<-read.table('Hsp90aN_affinities_fixed.csv',sep=',',quote = '"',header = TRUE)
aff<-aff[!duplicated(aff), ]

afflog<-aff[,1:2]
#afflog[,3:8]<-apply(aff[,3:8],1,log10)
afflog$KiM_avg<-apply(log10(aff[,3:length(colnames(aff))]),1,mean,na.rm=TRUE)
afflog$KiM_sd<-apply(log10(aff[,3:length(colnames(aff))]),1,sd,na.rm=TRUE)
#Remove values for salts and solvents
afflog<-subset(afflog,! HET.ID %in%  c('SO4','SO3','ZN','CL'))

#Show structures with multiple ligands
subset(afflog,PDB.ID %in% afflog[which(duplicated(afflog$PDB.ID)),]$PDB.ID)


#Remove unwanted duplicates for structures with multiple ligands
affclean<-subset(afflog, !(PDB.ID=='2QFO' ))



#Get data together for summary and affinities
summ<-merge(summr,affclean,by.x='File',by.y='PDB.ID',all.x = TRUE)
summ<-merge(summ,pdbinfo,by.x='File',by.y='PDB.ID',all.x = TRUE)
summ$HasLigand<-TRUE
#Mark which structures have ligands
summ[summ$File %in% subset(pdb,is.na(Heteroatoms))$PDB.ID,'HasLigand']<-FALSE

#Read cavity volume data
vols<-read.table('Prepared/Volumes.csv',sep=',',quote = '"',header = TRUE)
vols$File<-toupper(vols$File)
vols[,c('File','Chain')]<- transform(vols, File = colsplit(File, "_", names = c('a', 'b')))$File
vols<-vols[,colnames(vols)[c(1,11,2:10)]]



#Calculate volume ratios
summ$Cav_VdW<-summ$Cavity.volume..A.3/summ$VdW.volume..A.3
summ$Clef_VdW<-summ$Cleft.volume..A.3/summ$VdW.volume..A.3
summ$Empty_VdW<-(summ$Cleft.volume..A.3+summ$Cavity.volume..A.3)/summ$VdW.volume..A.3
summ$Empty<-(summ$Cleft.volume..A.3+summ$Cavity.volume..A.3)


summaf<-subset(summ,!is.na(KiM_avg))
mregression<-lm(KiM_avg~Cavity.volume..A.3+Cleft.volume..A.3+VdW.volume..A.3,data=summaf)
#mregression<-lm(KiM_avg~Cav_VdW+Clef_VdW,data=summaf)


summaf$Ki_predict<-predict(mregression,summaf)
ggplot(summaf, aes(x=KiM_avg,y=Ki_predict))+geom_point()+stat_smooth(method='lm')


smregression<-summary(mregression)

termplot(mregression)
plot(mregression)





#'Melt' data - convert columns to rows with volume type being factor variable
sum_melt<-melt(summ,id.vars=colnames(summ)[c(1:7,14:20)], measure.vars = colnames(summ)[c(8:13,21:24)],
               variable.name = 'Volume_type', value.name='Volume')
sum_melt$Volume_type<-factor(sum_melt$Volume_type,
                             levels = c('Cavity.volume..A.3','Cleft.volume..A.3','Total.surface..A.2',
                                        'Protein.volume..A.3','Solvent.accessible.volume..A.3','VdW.volume..A.3',
                                        'Cav_VdW','Clef_VdW','Empty_VdW','Empty'),
                             labels = c('Cavities','Clefts','Surface','Protein','Solvent accessible','VdW',
                                        'Cavities/VdW','Clefts/VdW','All empty/VdW','All empty'))


voleaff<-ggplot(subset(sum_melt,Volume_type %in% c('Cavities','Clefts','All empty')),
               aes(x=KiM_avg,y=Volume))+
  geom_errorbarh(aes(xmax=KiM_avg+KiM_sd,xmin=KiM_avg-KiM_sd),height=.001,alpha=0.2)+
  geom_point()+facet_grid(. ~Volume_type)+stat_smooth(method = "lm")+
  xlab('Log10 Ki, M')+ylab('Volume, A^3')+ggtitle('Correlation between volumes and ligand affinities')
voleaff
dev.copy2pdf(device=cairo_pdf,file="Figures/Cavities_vs_affinities.pdf",width=8,height=6)


voltaff<-ggplot(subset(sum_melt,Volume_type %in% c('Protein','VdW')),
                aes(x=KiM_avg,y=Volume))+
  geom_errorbarh(aes(xmax=KiM_avg+KiM_sd,xmin=KiM_avg-KiM_sd),height=.001,alpha=0.2)+
  geom_point()+facet_grid(Volume_type~.,scale='free_y')+stat_smooth(method = "lm")+
  xlab('Log10 Ki, M')+ylab('Volume')+ggtitle('Correlation between volumes and ligand affinities')
voltaff
dev.copy2pdf(device=cairo_pdf,file="Figures/Volumes_vs_affinities.pdf",width=8,height=6)

volrat<-subset(sum_melt,Volume_type %in% c('Cavities/VdW','Clefts/VdW','All empty/VdW')& Resolution<2)
volraff<-ggplot( volrat,
                aes(x=KiM_avg,y=Volume))+
  geom_errorbarh(aes(xmax=KiM_avg+KiM_sd,xmin=KiM_avg-KiM_sd),height=.001,alpha=0.2)+
  geom_point()+facet_grid(. ~Volume_type)+stat_smooth(method = "lm")+
  xlab('Log10 Ki, M')+ylab('Volume ratio')+ggtitle('Correlation between volumes and ligand affinities')
volraff
dev.copy2pdf(device=cairo_pdf,file="Figures/Ratios_vs_affinities.pdf",width=8,height=6)




#Analysis of cavities
onlyP<-subset(vols,Fragment=='P')
rownames(onlyP)<-paste(onlyP$File,'-',onlyP$Chain,'_',onlyP$Cavity.number,sep='')
onlyP<-merge(onlyP,pdbinfo,by.x='File',by.y='PDB.ID')
filenames<-as.character(unique(onlyP$File))
onlyP$Index<-paste(onlyP$File,onlyP$Chain,sep='-')


#volcor<-round(cor(onlyP[,c(6,8:10)]),2)
#pc <- princomp(onlyP[,c(6,8:10)], cor=TRUE, scores=TRUE)

#plot(pc,type="lines")
#biplot(pc)



plot3d(onlyP[,c('x','y','z')], col=ifelse(onlyP$Type=='cavity','red','blue'),xlab='X',ylab='Y',zlab='Z')
#rgl.postscript( '3D_cav_clef', fmt = "pdf", drawText = TRUE )



#Clustering part
clusts<-5

#Hierarchical clustering
di <- dist(onlyP[,c(7:9)], method="euclidean")
tree <- hclust(di, method="ward.D")
onlyP$hcluster <- as.factor((cutree(tree, k=clusts)-2) %% clusts +1)
# that modulo business just makes the coming table look nicer

#plot(tree, xlab="Cavities/Clefts")
#rect.hclust(tree, k=clusts, border="red")


plot3d(onlyP[,c('x','y','z')], col=onlyP$hcluster,xlab='X',ylab='Y',zlab='Z')
#rgl.postscript( '3D_clust', fmt = "pdf", drawText = TRUE )

#Volume distribution
voldist<-ggplot(onlyP,aes(x=hcluster,y=Volume,color=Type))+
  geom_point(size=3)+geom_boxplot()+
  labs(color='Volumes')+xlab('Cluster')+ylab(expression(paste("Volume, ",ring(A)^3)))
voldist
dev.copy2pdf(device=cairo_pdf,file=paste('Figures/Cluster_vol_dist_',clusts,'.pdf',sep=''),width=8,height=6)
#geom_point(aes(x=subset(onlyP,File=='1uyl')$hcluster,y=subset(onlyP,File=='1uyl')$Volume,shape="1UYL"),color="red")

clusters_C<-ddply(onlyP, .(hcluster), summarise, x=mean(x), y=mean(y), z=mean(z), Volume=mean(Volume))
write.csv(onlyP,file=paste('Centers/CavCl_centers_',clusts,'.csv',sep=''))
write.csv(clusters_C,file=paste('Centers/Cluster_centers_',clusts,'.csv',sep=''))

#Data prep
file_clusters<-ddply(onlyP, .(File,Index,Chain,hcluster,Resolution), summarise, x=mean(x), y=mean(y), z=mean(z), Volume=sum(Volume))



#Ligand affinities
affinities<-merge(file_clusters,afflog,by.x = 'File',by.y = 'PDB.ID')#,all.x = TRUE,all.y=TRUE
affinities<-subset(affinities,!is.na(KiM_avg))

#subset(affinities,Resolution<=2)
affcor<-ggplot(affinities,aes(y=Volume,x=KiM_avg,color=Resolution))+
  geom_point()+stat_smooth(aes(group = 1),method = "lm")+
  ylab(expression(paste("Volume, ",ring(A)^3)))+xlab('Log Ki, M')+
  facet_grid(.~'Cluster'+hcluster)+labs(color=expression(paste("Resolution, ",ring(A))))
affcor
dev.copy2pdf(device=cairo_pdf,file=paste('Figures/Ki-Vol_corr_',clusts,'.pdf',sep=''),width=8,height=6)


volaff<-dcast(affinities,File+Index+Chain +KiM_avg+KiM_sd ~ hcluster,sum,value.var = 'Volume')
cavregression<-lm(volaff$KiM_avg~volaff$'1'+volaff$'2'+volaff$'3'+volaff$'4'+volaff$'5',data=volaff)
#+volaff$'7' +volaff$'6'
volaff$KiM_predict<-predict(cavregression,volaff)
ggplot(volaff, aes(x=KiM_avg,y=KiM_predict))+
  geom_point()+stat_smooth(method='lm')+
  xlab('Log Ki, M')+ylab('Log Ki_predicted, M')+
  ggtitle('Correlation between Ki and volume predicted Ki')
dev.copy2pdf(device=cairo_pdf,file=paste('Figures/Ki_Ki-pred_',clusts,'.pdf',sep=''),width=8,height=6)




#Dendrogram

dendr    <- dendro_data(tree, type="rectangle") # convert for ggplot
clust.df <- data.frame(label=rownames(onlyP), cluster=onlyP$hcluster)
dendr[["labels"]]   <- merge(dendr[["labels"]],clust.df, by="label")
rect <- aggregate(x~cluster,label(dendr),range)
rect <- data.frame(rect$cluster,rect$x)
ymax <- mean(tree$height[length(tree$height)-((clusts-2):(clusts-1))])

denplot<-ggplot() + 
  geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend)) + 
  geom_text(data=label(dendr), aes(x, y, label=label, hjust=0, color=cluster), 
            size=0.5) +
  geom_rect(data=rect, aes(xmin=X1-.3, xmax=X2+.3, ymin=0, ymax=ymax), 
            color="red", fill=NA,alpha=0.5)+
  geom_hline(yintercept=0.33, color="blue")+labs(color='Cluster') +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  coord_flip() + scale_y_reverse(expand=c(0.2, 0)) + theme_dendro()+
  ggtitle(paste('Hierarchical clustering of cavities and clefts by their coordinates with',clusts,'clusters'))
denplot
dev.copy2pdf(device=cairo_pdf,file=paste('Figures/Cluster_dendrogram_',clusts,'.pdf',sep=''),width=8,height=8)



#For correlation
aligned<-merge(file_clusters[,!names(file_clusters) %in% c('x','y','z')],
               file_clusters[,!names(file_clusters) %in% c('x','y','z')],by = 'hcluster')
aligned<-rename(aligned,c('File.x'='Reference_F','File.y'='Target_F','Volume.x'='Volume_R','Volume.y'='Volume_T',
                          'Resolution.x'='Resolution_R','Resolution.y'='Resolution_T',
                          'Chain.x'='Chain_R','Chain.y'='Chain_T',
                          'Index.x'='Reference','Index.y'='Target'))

filesel<-unique(file_clusters$Index)
#Cluster correlation comparison example
vran<-0:5
onlycor<-subset(aligned,Reference %in% filesel[vran] & Target %in% filesel[vran])
volcor<-ggplot(onlycor,aes(x=Volume_R,y=Volume_T,color=hcluster))+
  geom_point()+stat_smooth(aes(group = 1),method = "lm")+xlim(0,100)+ylim(0,100)+
  facet_grid(Reference~Target)+xlab(expression(paste("Volume, ",ring(A)^3)))+
  ylab(expression(paste("Volume, ",ring(A)^3)))+labs(color='Cluster')
volcor
dev.copy2pdf(device=cairo_pdf,file="Figures/Cluster_comparison_example.pdf",width=8,height=7)




voldist2<-ggplot(file_clusters,aes(x=hcluster,y=Volume))+geom_boxplot()+
  geom_point(data=subset(file_clusters,File %in% subset(summ,HasLigand==FALSE)$File),
             aes(color=File),size=3)+
  labs(color='Structures without ligand')+xlab('Cluster')+
  ylab(expression(paste("Volume, ",ring(A)^3)))
voldist2
dev.copy2pdf(device=cairo_pdf,file="Figures/Cluster_variation.pdf",width=8,height=6)



onlyP$CC<-1

getclusts<-function(onlyP,variable,vols){
  presence<-dcast(onlyP, Index  ~ hcluster, sum, value.var=variable)
  head(presence)
  rownames(presence)<-presence$Index
  presence$Index<-NULL
  presence<-as.matrix(presence)
  mode(presence)<-'numeric'
  if (!vols){
    presence[presence>1]<-1
  }
  #presence[presence<1]<-NA
  presence<-data.frame(presence)
  colnames(presence)<-gsub('X','',colnames(presence))
  
  return(presence)
}

presence<-t(getclusts(onlyP,'CC',FALSE))
clvols<-t(getclusts(onlyP,'Volume',TRUE))





heatmap.2(as.matrix(clvols),Rowv=TRUE,Colv=TRUE,trace='none',col=heat.colors(256),
          xlab='Structure',ylab='Cluster',key=TRUE,scale='none',
          dendrogram='both',na.color="grey",cexRow=1.5,cexCol=0.5,main='Ertmiu klasteriu turiai')
#dev.copy2pdf(device=cairo_pdf,file="Cavity-cluster_volume.pdf",width=16,height=9)

heatmap.2(as.matrix(presence),Rowv=TRUE,Colv=TRUE,trace='none',col=c('white','red'),
          xlab='Structure',ylab='Cluster',key=FALSE,scale='none',
          dendrogram='both',na.color="grey",cexRow=1.5,cexCol=0.5,
          colsep=1:ncol(presence),rowsep=1:nrow(presence),sepcolor="black",sepwidth=c(0.001,0.001),
          main='Cluster presence in structures')
#dev.copy2pdf(device=cairo_pdf,file="Cavity-cluster_presence.pdf",width=16,height=9)


filler <- function(frame,names){
  tot<-length(names)
  a<-data.frame(matrix(NA, nrow=tot, ncol=tot))

  colnames(a)<-names
  rownames(a)<-names
  b<-a
  R<-a
  
  i<-0
  
  align<-data.frame()
  
  for (r in names){
    for (t in names) {

      res<-summary(lm(frame[,r]~frame[,t]))
      b[t,r]<-res$coefficients[1]
      a[t,r]<-res$coefficients[2]
      R[t,r]<-res$r.squared
      i<-i+1
    }
    print(i*100/(tot^2))
  }
  
  return(list(a=a,b=b,R=R))
}





res<-filler(clvols,filesel)
#R<-res$R
#a<-res$a



heatmap.2(as.matrix(res$R),Colv=TRUE,Rowv=TRUE,trace='none',col=heat.colors(64),
          xlab="Structure 1",ylab='Structure 2',key=TRUE,scale='none',
          dendrogram=,na.color="grey",symm = TRUE,cexRow=0.5,cexCol=0.5)
dev.copy2pdf(device=cairo_pdf,file=paste("Figures/Heatmap_R-",clusts,".pdf",sep=''),width=12,height=12)
#CairoSVG(file = paste("Figures/Heatmap_R",clusts,".svg",sep=''), width = 6, height = 6)
