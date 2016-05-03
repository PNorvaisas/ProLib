library('ggplot2')
#library(ggdendro)
library('gplots')
library('plyr')
#library('plot3D')
library('reshape2')
library(tidyr)
library(rgl)
#library(cluster)

#Skriptas veiks paleidus baltymo Results direktorijoje!

#Read summary of volumes
summ<-read.table('Prepared/Summary.csv',sep=',',quote = '"',header = TRUE)
summ$File<-toupper(summ$File)
#Split column values to extract PDB IDs and Chain IDs
summ[,c('File','Chain')]<- transform(summ, File = colsplit(File, "_", names = c('a', 'b')))$File
summ<-summ[,colnames(summ)[c(1,13,2:12)]]

#Read summary of PDB info
pdb<-read.table('Summary_sequences.csv',sep=',',quote = '"',header = TRUE)
pdb$PDB.ID<-toupper(pdb$PDB.ID)
pdb<-pdb[,1:7]
#Make empty fields as not-available
pdb[pdb=='']<-NA

#Structures without ligands
subset(pdb,is.na(Heteroatoms))



#Read affinties
#Filenames differ for each protein
aff<-read.table('Hsp90aN_affinities_fixed.csv',sep=',',quote = '"',header = TRUE)
afflog<-aff[,1:2]
#Generate log averages and sd
afflog$KiM_avg<-apply(log10(aff[,3:length(colnames(aff))]),1,mean,na.rm=TRUE)
afflog$KiM_sd<-apply(log10(aff[,3:length(colnames(aff))]),1,sd,na.rm=TRUE)
#Remove values for salts and solvents
afflog<-subset(afflog,! HET.ID %in%  c('SO4','SO3','ZN','CL'))

#Show structures with multiple ligands
subset(afflog,PDB.ID %in% afflog[which(duplicated(afflog$PDB.ID)),]$PDB.ID)

#Remove unwanted duplicates for structures with multiple ligands
#Now there will be only one Ki_M affinity value for one structure
#CAII
afflog<-subset(afflog, !(PDB.ID=='1AVN' & HET.ID=='AZI'))
afflog<-subset(afflog, !(PDB.ID=='3CYU' & HET.ID=='0CR'))
afflog<-subset(afflog, !(PDB.ID=='3T83' & HET.ID=='MG5'))
#Hsp90


#CDK


#Get data together for summary and affinities
summ<-merge(summ,afflog,by.x='File',by.y='PDB.ID',all.x = TRUE)
summ$HasLigand<-TRUE
#Mark which structures don't have ligands
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


summaf$Ki_predict<-predict(mregression,summaf)
ggplot(summaf, aes(x=KiM_avg,y=Ki_predict))+geom_point()+stat_smooth(method='lm')


smregression<-summary(mregression)



termplot(mregression)
plot(mregression)


#'Melt' data - convert columns to rows with volume type being factor variable
sum_melt<-melt(summ,id.vars=colnames(summ)[c(1:6,13:17)], measure.vars = colnames(summ)[c(7:12,18:21)], variable.name = 'Volume_type' , value.name='Volume')
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
  geom_point()+facet_grid(. ~Volume_type,scale='free_y')+stat_smooth(method = "lm")+
  xlab('Log10 Ki, M')+ylab('Volume')+ggtitle('Correlation between volumes and ligand affinities')
voltaff
dev.copy2pdf(device=cairo_pdf,file="Figures/Volumes_vs_affinities.pdf",width=8,height=6)

volraff<-ggplot(subset(sum_melt,Volume_type %in% c('Cavities/VdW','Clefts/VdW','All empty/VdW')),
                aes(x=KiM_avg,y=Volume))+
  geom_errorbarh(aes(xmax=KiM_avg+KiM_sd,xmin=KiM_avg-KiM_sd),height=.001,alpha=0.2)+
  geom_point()+facet_grid(. ~Volume_type)+stat_smooth(method = "lm")+
  xlab('Log10 Ki, M')+ylab('Volume ratio')+ggtitle('Correlation between volumes and ligand affinities')
volraff
dev.copy2pdf(device=cairo_pdf,file="Figures/Ratios_vs_affinities.pdf",width=8,height=6)


