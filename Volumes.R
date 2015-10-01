library('ggplot2')
library(ggdendro)
library('gplots')
library('plyr')
library('plot3D')
library('reshape2')
library(tidyr)
library(rgl)
library(cluster)

vols<-read.table('Volumes.csv',sep=',',quote = '"',header = TRUE)

#vols$Type<-factor(vols$Type)

onlyP<-subset(vols,Fragment=='P')
rownames(onlyP)<-paste(onlyP$File,onlyP$Cavity.number,sep='-')
filenames<-as.character(unique(onlyP$File))


volcor<-round(cor(onlyP[,c(5,7:9)]),2)
pc <- princomp(onlyP[,c(5,7:9)], cor=TRUE, scores=TRUE)

plot(pc,type="lines")

biplot(pc)



plot3d(onlyP[,c('x','y','z')], col=ifelse(onlyP$Type=='cavity','red','blue'),xlab='X',ylab='Y',zlab='Z')
rgl.postscript( '3D_cav_clef', fmt = "pdf", drawText = TRUE )

set.seed(42)
cl <- kmeans(onlyP[,c(7:9)],40)
dst<-dist(onlyP[,c(7:9)])
onlyP$cluster <- as.factor(cl$cluster)



with(onlyP, table(cluster, File))

clusts<-6
di <- dist(onlyP[,c(7:9)], method="euclidean")
tree <- hclust(di, method="ward.D")
onlyP$hcluster <- as.factor((cutree(tree, k=clusts)-2) %% clusts +1)
# that modulo business just makes the coming table look nicer

#plot(tree, xlab="Cavities/Clefts")
#rect.hclust(tree, k=clusts, border="red")

dendr    <- dendro_data(tree, type="rectangle") # convert for ggplot
clust.df <- data.frame(label=rownames(onlyP), cluster=onlyP$hcluster)
dendr[["labels"]]   <- merge(dendr[["labels"]],clust.df, by="label")
rect <- aggregate(x~cluster,label(dendr),range)
rect <- data.frame(rect$cluster,rect$x)
ymax <- mean(tree$height[length(tree$height)-((clusts-2):(clusts-1))])

ggplot() + 
  geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend)) + 
  geom_text(data=label(dendr), aes(x, y, label=label, hjust=0, color=cluster), 
            size=0.8) +
  geom_rect(data=rect, aes(xmin=X1-.3, xmax=X2+.3, ymin=0, ymax=ymax), 
            color="red", fill=NA)+
  geom_hline(yintercept=0.33, color="blue")+
  labs(color='Klasteris') +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  coord_flip() + scale_y_reverse(expand=c(0.2, 0)) + 
  theme_dendro()




plot3d(pc$scores[,1:4], col=onlyP$hcluster)

plot3d(onlyP[,c('x','y','z')], col=onlyP$hcluster,xlab='X',ylab='Y',zlab='Z')
rgl.postscript( '3D_clust', fmt = "pdf", drawText = TRUE )


voldist<-ggplot(onlyP,aes(x=hcluster,y=Volume,fill=Type))+geom_boxplot()+geom_point(data=subset(onlyP,File %in% c('1uyl','2yeg','3t0h')),aes(color=Type,shape=File),size=6)   #geom_point(aes(x=subset(onlyP,File=='1uyl')$hcluster,y=subset(onlyP,File=='1uyl')$Volume,shape="1UYL"),color="red")
voldist
dev.copy2pdf(device=cairo_pdf,file="Cluster_vol_dist.pdf",width=8,height=6)






write.csv(onlyP,file='Centers/CavCl_centers.csv')

file_clusters<-ddply(onlyP, .(File,hcluster), summarise, x=mean(x), y=mean(y), z=mean(z), Volume=sum(Volume))
clusters_C<-ddply(onlyP, .(hcluster), summarise, x=mean(x), y=mean(y), z=mean(z), Volume=mean(Volume))
write.csv(clusters_C,file='Centers/Cluster_centers.csv')

aligned<-merge(file_clusters[,!names(file_clusters) %in% c('x','y','z')],file_clusters[,!names(file_clusters) %in% c('x','y','z')],by = 'hcluster')
aligned<-rename(aligned,c('File.x'='Reference','File.y'='Target','Volume.x'='Volume_R','Volume.y'='Volume_T'))

vran<-10:30
onlycor<-subset(aligned,Reference %in% filenames[vran] & Target %in% filenames[vran])
volcor<-ggplot(onlycor,aes(x=Volume_R,y=Volume_T))+geom_point()+stat_smooth(aes(group = 1),method = "lm")+xlim(0,100)+ylim(0,100)
volcor+facet_grid(Reference~Target)+xlab('Struktūra 1')+ylab('Struktūra 2')



onlyP$CC<-1

getclusts<-function(onlyP,variable,vols){
  presence<-dcast(onlyP, File  ~ hcluster, sum, value.var=variable)
  rownames(presence)<-presence$File
  presence$File<-NULL
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
dev.copy2pdf(device=cairo_pdf,file="Cavity-cluster_volume.pdf",width=16,height=9)

heatmap.2(as.matrix(presence),Rowv=TRUE,Colv=TRUE,trace='none',col=c('white','red'),
          xlab='Structure',ylab='Cluster',key=FALSE,scale='none',
          dendrogram='both',na.color="grey",cexRow=1.5,cexCol=0.5,
          colsep=1:ncol(presence),rowsep=1:nrow(presence),sepcolor="black",sepwidth=c(0.01,0.01),
          main='Ertmiu klusteriu reprezentacija strukturose')
dev.copy2pdf(device=cairo_pdf,file="Cavity-cluster_presence.pdf",width=16,height=9)


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





res<-filler(clvols,filenames)
R<-res$R
a<-res$a



heatmap.2(as.matrix(R),Colv=TRUE,Rowv=TRUE,trace='none',col=heat.colors(256),
          xlab='Reference',ylab='Target',key=TRUE,scale='none',
          dendrogram=,na.color="grey",symm = TRUE)








aff<-read.table('Affinities.csv',sep=',',quote = '"',header = TRUE)


crop<-read.table('Summary.csv',sep=',',quote = '"',header = TRUE)
uncrop<-read.table('Uncropped/Summary.csv',sep=',',quote = '"',header = TRUE)

crop$Crop<-TRUE
uncrop$Crop<-FALSE


allsum<-merge(crop,uncrop, all.x = TRUE, all.y = TRUE)
allsum$Ligand<-FALSE

allsum[allsum$File %in% unique(allsum[allsum$Fragment=='L','File']),]$Ligand<-TRUE

#allsum<-lapply(allsum[colnames(allsum)[c(1:5,13,14)]], as.factor)



allvols<-melt(allsum,id.vars=colnames(allsum)[c(1:5,13:14)], measure.vars = colnames(allsum)[c(6:11)], variable.name = 'Volume_type' , value.name='Volume')
allfrags<-dcast(allvols, File + Ligand.name + Crop + Ligand + Volume_type ~ Fragment, value.var='Volume')


fragcompLPL<-ggplot(subset(allfrags,Volume_type %in% c('Protein.volume..A.3')),aes(x=L,y=PL,color=Crop))+geom_point()+stat_smooth(method = "lm")
fragcompLPL


fragcompPPL<-ggplot(subset(allfrags,Volume_type %in% c('Protein.volume..A.3')),aes(x=P,y=PL,color=Crop))+geom_point()+stat_smooth(method = "lm")
fragcompPPL

allsum<-subset(allsum,Fragment!='L' ) #& File !%in% c('3t10','3t0h')
allvols<-subset(allvols,Fragment!='L' )

crop_p<-ggplot(subset(allvols,Volume_type %in% c("Protein.volume..A.3")),aes(x=Volume, fill=Crop))+geom_histogram(binwidth=100,position='identity',alpha=0.5)
crop_p+facet_grid(.~Fragment)+xlab('Baltymo turis, A^3')
dev.copy2pdf(device=cairo_pdf,file="Protein_volume_bef-af_cut.pdf",width=8,height=6)


Pvols<-ggplot(subset(allvols,Volume_type %in% c("Protein.volume..A.3")),aes(x=Fragment, y=Volume, fill=Crop))+geom_boxplot()+geom_point(data=subset(allvols,File %in% c('1uyl','2yeg','3t0h') & Volume_type %in% c("Protein.volume..A.3")),aes(color=Crop),size=6)   #geom_point(aes(x=subset(onlyP,File=='1uyl')$hcluster,y=subset(onlyP,File=='1uyl')$Volume,shape="1UYL"),color="red")
Pvols+xlab('Baltymo turis, A^3')



crop_c<-ggplot(subset(allvols,Volume_type %in% c("Cleft.volume..A.3","Cavity.volume..A.3")),aes(x=Volume, fill=Crop))+geom_histogram(binwidth=20,position='identity',alpha=0.5)
crop_c+facet_grid(Volume_type~Fragment)+xlab('Turis, A^3')
dev.copy2pdf(device=cairo_pdf,family="DejaVu Sans",file="Cav-Cleft_volume_bef-af_cut.pdf",width=8,height=6)

extrm<-subset(allsum,Fragment=='P' & Crop==TRUE)
head(extrm[with(extrm,order(-Protein.volume..A.3)),])

#pairs(allsum[,c('Fragment','Crop','Ligand','Protein.volume..A.3',"Solvent.accessible.volume..A.3","Cavity.volume..A.3","Cleft.volume..A.3")])


volcomp<-ggplot(allsum,aes(x=Protein.volume..A.3,y=VdW.volume..A.3,color=Ligand))+geom_point()
volcomp+facet_grid(Crop~Fragment)

cavcomp<-ggplot(subset(allsum,Ligand==TRUE),aes(x=Cavity.volume..A.3,y=Cleft.volume..A.3))+geom_point()+stat_smooth(method = "lm")
cavcomp+facet_grid(Crop~Fragment)

vclcomp<-ggplot(allsum,aes(x=Protein.volume..A.3,y=Cleft.volume..A.3,color=Ligand))+geom_point()
vclcomp+facet_grid(Crop~Fragment)

vcacomp<-ggplot(allsum,aes(x=Protein.volume..A.3,y=Cavity.volume..A.3,color=Ligand))+geom_point()
vcacomp+facet_grid(Crop~Fragment)
