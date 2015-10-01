library('ggplot2')
library('gplots')
library('plyr')
library('plot3D')
library('reshape2')
library(tidyr)
library(rgl)



align<-read.table('Alignment.csv',sep=',',quote = '"',header = TRUE)
align<-with(align, align[order(Reference,Target),])
paired<-subset(align,! is.na(R.index) & !is.na(T.index!=''))


pdbs<-as.character(unique(paired$Reference))[1:20] 
als<-ggplot(subset(paired,Reference=='1uyl' & Target=='1uy8' & R.Frag=='P' & T.Frag=='P'),aes(x=R.Volume,y=T.Volume,color=R.Type))+geom_point()+geom_point(aes(color=T.Type),alpha=0.3,size=5)+xlim(0,150)+ylim(0,150)+stat_smooth(aes(group = 1),method = "lm")
als+facet_grid(Reference~Target)






