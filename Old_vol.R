library("plyr")
library("ggplots")
library("ggplot2")




data<-read.csv('R_P_results.csv')

rownames(data)<-data$File
data$File<-NULL

rownames(data)<-data[,0]
data$File
colanmes(data)
colnames(data)
data[0,]
data[1,]
data <- read.csv('R_P_results.csv', header=TRUE, row.names=1)
data
gsub("X","",colnames(data))
colanmes(data)<-gsub("X","",colnames(data))
colnames(data)<-gsub("X","",colnames(data))
data

ref<-data["1uyl",]
boxplot(data)
b<-boxplot(data)
b<-plot(ref)

ggplot2
ggplot
plot1<-ggplot(data=data)+greom_boxplot()
plot1<-ggplot(data=data)+geom_boxplot()
plot1
plot1<-ggplot(data=data,aes=(x=data[0,],y=rownames(data)))+geom_boxplot()
plot1<-ggplot(data=data,aes=(x=data,y=rownames(data)))+geom_boxplot()
plot1<-ggplot(data=data,aes=(x=data,y=data))+geom_boxplot()
plot1<-ggplot(data=data,aes(x=data,y=data))+geom_boxplot()
plot1
data[0,]
data[1,]
data[0,]
ggplot(data=data,aes(x=data[0,],y=data[,0]))+geom_boxplot()
dataT<-t(data)
dataT
ggplot(data=dataT,aes(x=data[0,],y=data[,0]))+geom_boxplot()
ggplot(data=data,aes(x=data[0,],y=data[,0]))+geom_boxplot()
data
boxplot(data)
ggplot2.boxplot(data)
ggplot(data=data)
ggplot(data=data)+geom_boxplot()
ggplot(data=data, aes(x=File,y="1"))+geom_boxplot()
boxplot(data)
par(new=T)
plot(ref)
plot(data[0,],ref)
plot(data)
par(new=F)
plot(data)
plot(data[0,],data)
plot(data)
data<read.csv('All_alignments.csv')
data<-read.csv('All_alignments.csv')
data
data<-read.csv('All_alignments.csv',header=TRUE)
colnames(data)
data$Type
data<-read.csv('All_alignments.csv',header=TRUE)
colnames(data)
comparison<-data[,Type="R_p"]
comparison<-data[,Type="R_P"]
comparison<-data[which(data$Type=="R_P")]
comparison<-data[which(data$Type=="R_P"),c("R.Number","R.Volume","P.Volume")]

comparison<-data[which(data$Type=="R_P" & data$P.Number!=NA),c("R.Number","R.Volume","P.Volume")]

comparison<-data[which(data$Type=="R_P" & data$P.Number!=NA),c("R.Number","R.Volume","P.Volume")]
comparison<-data[which(data$Type=="R_P" & data$P.Number!="NA"),c("R.Number","R.Volume","P.Volume")]

comparison<-data[which(data$Type=="R_P" & data$P.Number!="NA" & data$R.Number!="NA"),c("R.Number","R.Volume","P.Volume")]

ggplot(data=comparison,aes(x=factor(R.Number),y=P.Volume))+geom_boxplot()+geom_point(aes(x=factor(R.Number),y=R.Volume,shape="1UYL"),color="red")+xlab("Ertmės indeksas kontrolinėje struktūroje 1UYL") + ylab(expression(paste("Tūris, A"^3)))+ggtitle(expression(paste("Ertmių tūrių variacija Hsp90",alpha,"N")))+guides(color=guide_legend('Duomenys'))+stat_summary(fun.data="3",geom="text",colour="red",fun.y=mean)

n_fun<-function(x){return(data.frame(y=median(x),label=paste0("n=",length(x))))}
ggplot(data=comparison,aes(x=factor(R.Number),y=P.Volume))+geom_boxplot()+geom_point(aes(x=factor(R.Number),y=R.Volume,shape="1UYL"),color="red")+xlab("Ertmės indeksas kontrolinėje struktūroje 1UYL") + ylab(expression(paste("Tūris, A"^3)))+ggtitle(expression(paste("Ertmių tūrių variacija Hsp90",alpha,"N")))+guides(color=guide_legend('Duomenys'))+stat_summary(fun.data=n_fun,geom="text",colour="red")
n_fun<-function(x){return(data.frame(y=median(x),label=length(x)))}
ggplot(data=comparison,aes(x=factor(R.Number),y=P.Volume))+geom_boxplot()+geom_point(aes(x=factor(R.Number),y=R.Volume,shape="1UYL"),color="red")+xlab("Ertmės indeksas kontrolinėje struktūroje 1UYL") + ylab(expression(paste("Tūris, A"^3)))+ggtitle(expression(paste("Ertmių tūrių variacija Hsp90",alpha,"N")))+guides(color=guide_legend('Duomenys'))+stat_summary(fun.data=n_fun,geom="text",colour="red")
n_fun<-function(x){return(data.frame(y=median(x)*1.1,label=length(x)))}
ggplot(data=comparison,aes(x=factor(R.Number),y=P.Volume))+geom_boxplot()+geom_point(aes(x=factor(R.Number),y=R.Volume,shape="1UYL"),color="red")+xlab("Ertmės indeksas kontrolinėje struktūroje 1UYL") + ylab(expression(paste("Tūris, A"^3)))+ggtitle(expression(paste("Ertmių tūrių variacija Hsp90",alpha,"N")))+guides(color=guide_legend('Duomenys'))+stat_summary(fun.data=n_fun,geom="text",colour="red")
ggplot(data=comparison,aes(x=factor(R.Number),y=P.Volume))+geom_boxplot()+geom_point(aes(x=factor(R.Number),y=R.Volume,shape="1UYL"),color="red")+xlab("Ertmės indeksas kontrolinėje struktūroje 1UYL") + ylab(expression(paste("Tūris, A"^3)))+ggtitle(expression(paste("Ertmių tūrių variacija Hsp90",alpha,"N")))+guides(color=guide_legend('Duomenys'))+stat_summary(fun.data=n_fun,geom="text",colour="blue")
ggplot(data=comparison,aes(x=factor(R.Number),y=P.Volume))+geom_boxplot()+geom_point(aes(x=factor(R.Number),y=R.Volume,shape="1UYL"),color="red")+xlab("Ertmės indeksas kontrolinėje struktūroje 1UYL") + ylab(expression(paste("Tūris, A"^3)))+ggtitle(expression(paste("Ertmių tūrių variacija Hsp90",alpha,"N")))+guides(color=guide_legend('Duomenys'))+stat_summary(aes(label=..count..)geom="text",colour="blue")
ggplot(data=comparison,aes(x=factor(R.Number),y=P.Volume))+geom_boxplot()+geom_point(aes(x=factor(R.Number),y=R.Volume,shape="1UYL"),color="red")+xlab("Ertmės indeksas kontrolinėje struktūroje 1UYL") + ylab(expression(paste("Tūris, A"^3)))+ggtitle(expression(paste("Ertmių tūrių variacija Hsp90",alpha,"N")))+guides(color=guide_legend('Duomenys'))+stat_summary(aes(label=..count..),geom="text",colour="blue")
ggplot(data=comparison,aes(x=factor(R.Number),y=P.Volume))+geom_boxplot()+geom_point(aes(x=factor(R.Number),y=R.Volume,shape="1UYL"),color="red")+xlab("Ertmės indeksas kontrolinėje struktūroje 1UYL") + ylab(expression(paste("Tūris, A"^3)))+ggtitle(expression(paste("Ertmių tūrių variacija Hsp90",alpha,"N")))+guides(color=guide_legend('Duomenys'))+stat_summary(aes(label=..count..),geom="text",colour="blue",position='identity')
ggplot(data=comparison,aes(x=factor(R.Number),y=P.Volume))+geom_boxplot()+geom_point(aes(x=factor(R.Number),y=R.Volume,shape="1UYL"),color="red")+xlab("Ertmės indeksas kontrolinėje struktūroje 1UYL") + ylab(expression(paste("Tūris, A"^3)))+ggtitle(expression(paste("Ertmių tūrių variacija Hsp90",alpha,"N")))+guides(color=guide_legend('Duomenys'))+stat_summary(fun.data=n_fun,geom="text",colour="blue",size=10)
ggplot(data=comparison,aes(x=factor(R.Number),y=P.Volume))+geom_boxplot()+geom_point(aes(x=factor(R.Number),y=R.Volume,shape="1UYL"),color="red")+xlab("Ertmės indeksas kontrolinėje struktūroje 1UYL") + ylab(expression(paste("Tūris, A"^3)))+ggtitle(expression(paste("Ertmių tūrių variacija Hsp90",alpha,"N")))+guides(color=guide_legend('Duomenys'))+stat_summary(fun.data=n_fun,geom="text",colour="blue",size=2)
ggplot(data=comparison,aes(x=factor(R.Number),y=P.Volume))+geom_boxplot()+geom_point(aes(x=factor(R.Number),y=R.Volume,shape="1UYL"),color="red")+xlab("Ertmės indeksas kontrolinėje struktūroje 1UYL") + ylab(expression(paste("Tūris, A"^3)))+ggtitle(expression(paste("Ertmių tūrių variacija Hsp90",alpha,"N")))+guides(color=guide_legend('Duomenys'))+stat_summary(fun.data=n_fun,geom="text",colour="blue",size=5)
ggplot(data=comparison,aes(x=factor(R.Number),y=P.Volume))+geom_boxplot()+geom_point(aes(x=factor(R.Number),y=R.Volume,shape="1UYL"),color="red")+xlab("Ertmės indeksas kontrolinėje struktūroje 1UYL") + ylab(expression(paste("Tūris, A"^3)))+ggtitle(expression(paste("Ertmių tūrių variacija Hsp90",alpha,"N")))+guides(color=guide_legend('Duomenys'))+stat_summary(fun.data=n_fun,geom="text",colour="blue",size=4)
n_fun<-function(x){return(data.frame(y=median(x)+20,label=length(x)))}
ggplot(data=comparison,aes(x=factor(R.Number),y=P.Volume))+geom_boxplot()+geom_point(aes(x=factor(R.Number),y=R.Volume,shape="1UYL"),color="red")+xlab("Ertmės indeksas kontrolinėje struktūroje 1UYL") + ylab(expression(paste("Tūris, A"^3)))+ggtitle(expression(paste("Ertmių tūrių variacija Hsp90",alpha,"N")))+guides(color=guide_legend('Duomenys'))+stat_summary(fun.data=n_fun,geom="text",colour="blue",size=4)
n_fun<-function(x){return(data.frame(y=-5,label=length(x)))}
ggplot(data=comparison,aes(x=factor(R.Number),y=P.Volume))+geom_boxplot()+geom_point(aes(x=factor(R.Number),y=R.Volume,shape="1UYL"),color="red")+xlab("Ertmės indeksas kontrolinėje struktūroje 1UYL") + ylab(expression(paste("Tūris, A"^3)))+ggtitle(expression(paste("Ertmių tūrių variacija Hsp90",alpha,"N")))+guides(color=guide_legend('Duomenys'))+stat_summary(fun.data=n_fun,geom="text",colour="blue",size=4)
cairo_pdf('1UYL_sulyginimas.pdf',family="DejaVu Sans", width=10)
ggplot(data=comparison,aes(x=factor(R.Number),y=P.Volume))+geom_boxplot()+geom_point(aes(x=factor(R.Number),y=R.Volume,shape="1UYL"),color="red")+xlab("Ertmės indeksas kontrolinėje struktūroje 1UYL") + ylab(expression(paste("Tūris, A"^3)))+ggtitle(expression(paste("Ertmių tūrių variacija Hsp90",alpha,"N")))+guides(color=guide_legend('Duomenys'))+stat_summary(fun.data=n_fun,geom="text",colour="blue",size=4)[A
dev.off()
dev(new)
ggplot(data=comparison,aes(x=factor(R.Number),y=P.Volume))+geom_boxplot()+geom_point(aes(x=factor(R.Number),y=R.Volume,shape="1UYL"),color="red")+xlab("Ertmės indeksas kontrolinėje struktūroje 1UYL") + ylab(expression(paste("Tūris, A"^3)))+ggtitle(expression(paste("Ertmių tūrių variacija Hsp90",alpha,"N")))+guides(color=guide_legend('Duomenys'))+stat_summary(fun.data=n_fun,geom="text",colour="blue",size=4)
ggplot(data=comparison,aes(x=factor(R.Number),y=P.Volume))+geom_boxplot()+geom_point(aes(x=factor(R.Number),y=R.Volume,shape="1UYL"),color="red")+xlab("Ertmės indeksas kontrolinėje struktūroje 1UYL") + ylab(expression(paste("Tūris, A"^3)))+ggtitle(expression(paste("Ertmių tūrių variacija Hsp90",alpha,"N")))+guides(color=guide_legend('Duomenys'))+stat_summary(fun.data=count,geom="text",colour="blue",size=4)
ggplot(data=comparison,aes(x=factor(R.Number),y=P.Volume))+geom_boxplot()+geom_point(aes(x=factor(R.Number),y=R.Volume,shape="1UYL"),color="red")+xlab("Ertmės indeksas kontrolinėje struktūroje 1UYL") + ylab(expression(paste("Tūris, A"^3)))+ggtitle(expression(paste("Ertmių tūrių variacija Hsp90",alpha,"N")))+guides(color=guide_legend('Duomenys'))+stat_summary(fun.data=,geom="text",colour="blue",size=4)
ggplot(data=comparison,aes(x=factor(R.Number),y=P.Volume))+geom_boxplot()+geom_point(aes(x=factor(R.Number),y=R.Volume,shape="1UYL"),color="red")+xlab("Ertmės indeksas kontrolinėje struktūroje 1UYL") + ylab(expression(paste("Tūris, A"^3)))+ggtitle(expression(paste("Ertmių tūrių variacija Hsp90",alpha,"N")))+guides(color=guide_legend(titl='Duomenys'))+stat_summary(fun.data=n_fun,geom="text",colour="blue",size=4)
ggplot(data=comparison,aes(x=factor(R.Number),y=P.Volume))+geom_boxplot()+geom_point(aes(x=factor(R.Number),y=R.Volume,shape="1UYL"),color="red")+xlab("Ertmės indeksas kontrolinėje struktūroje 1UYL") + ylab(expression(paste("Tūris, A"^3)))+ggtitle(expression(paste("Ertmių tūrių variacija Hsp90",alpha,"N")))+guides(color=guide_legend(title='Duomenys'))+stat_summary(fun.data=n_fun,geom="text",colour="blue",size=4)
ggplot(data=comparison,aes(x=factor(R.Number),y=P.Volume,shape="Median"))+geom_boxplot()+geom_point(aes(x=factor(R.Number),y=R.Volume,shape="1UYL"),color="red")+xlab("Ertmės indeksas kontrolinėje struktūroje 1UYL") + ylab(expression(paste("Tūris, A"^3)))+ggtitle(expression(paste("Ertmių tūrių variacija Hsp90",alpha,"N")))+guides(color=guide_legend(title='Duomenys'))+stat_summary(fun.data=n_fun,geom="text",colour="blue",size=4)
ggplot(data=comparison,aes(x=factor(R.Number),y=P.Volume))+geom_boxplot()+geom_point(aes(x=factor(R.Number),y=R.Volume,shape="1UYL"),color="red")+xlab("Ertmės indeksas kontrolinėje struktūroje 1UYL") + ylab(expression(paste("Tūris, A"^3)))+ggtitle(expression(paste("Ertmių tūrių variacija Hsp90",alpha,"N")))+guides(color=guide_legend(title='Duomenys'))+stat_summary(fun.data=n_fun,geom="text",colour="blue",size=4)
ggplot(data=comparison,aes(x=factor(R.Number),y=P.Volume))+geom_boxplot()+geom_point(aes(x=factor(R.Number),y=R.Volume,shape="1UYL"),color="red")+xlab("Ertmės indeksas kontrolinėje struktūroje 1UYL") + ylab(expression(paste("Tūris, A"^3)))+ggtitle(expression(paste("Ertmių tūrių variacija Hsp90",alpha,"N")))+guides(color=guide_legend(title="Duomenys"))+stat_summary(fun.data=n_fun,geom="text",colour="blue",size=4)
ggplot(data=comparison,aes(x=factor(R.Number),y=P.Volume))+geom_boxplot()+geom_point(aes(x=factor(R.Number),y=R.Volume,shape="1UYL"),color="red")+xlab("Ertmės indeksas kontrolinėje struktūroje 1UYL") + ylab(expression(paste("Tūris, A"^3)))+ggtitle(expression(paste("Ertmių tūrių variacija Hsp90",alpha,"N be ligandų")))+guides(color=guide_legend(title="Duomenys"))+stat_summary(fun.data=n_fun,geom="text",colour="blue",size=4)
data
RPL<-data[which(data$Type=="R_PL" & data$R.Volume!="NA", & data$PL.Volume!="NA"),c("R.Number","R.Volume","PL.Volume")]
RPL<-data[which(data$Type=="R_PL" & data$R.Volume!="NA" & data$PL.Volume!="NA"),c("R.Number","R.Volume","PL.Volume")]
RPL
cairo_pdf('1UYL_sulyginimas.pdf',family="DejaVu Sans", width=20)
ggplot(data=comparison,aes(x=factor(R.Number),y=P.Volume))+geom_boxplot()+geom_point(aes(x=factor(R.Number),y=R.Volume,shape="1UYL"),color="red")+xlab("Ertmės indeksas kontrolinėje struktūroje 1UYL") + ylab(expression(paste("Tūris, A"^3)))+ggtitle(expression(paste("Ertmių tūrių variacija Hsp90",alpha,"N be ligandų")))+guides(color=guide_legend(title="Duomenys"))+stat_summary(fun.data=n_fun,geom="text",colour="blue",size=4)
dev.off()
cairo_pdf('1UYL_sulyginimas.pdf',family="DejaVu Sans", width=15)
ggplot(data=comparison,aes(x=factor(R.Number),y=P.Volume))+geom_boxplot()+geom_point(aes(x=factor(R.Number),y=R.Volume,shape="1UYL"),color="red")+xlab("Ertmės indeksas kontrolinėje struktūroje 1UYL") + ylab(expression(paste("Tūris, A"^3)))+ggtitle(expression(paste("Ertmių tūrių variacija Hsp90",alpha,"N be ligandų")))+guides(color=guide_legend(title="Duomenys"))+stat_summary(fun.data=n_fun,geom="text",colour="blue",size=4)
dev.off()
cairo_pdf('1UYL_sulyginimas.pdf',family="DejaVu Sans", width=10)
ggplot(data=comparison,aes(x=factor(R.Number),y=P.Volume))+geom_boxplot()+geom_point(aes(x=factor(R.Number),y=R.Volume,shape="1UYL"),color="red")+xlab("Ertmės indeksas kontrolinėje struktūroje 1UYL") + ylab(expression(paste("Tūris, A"^3)))+ggtitle(expression(paste("Ertmių tūrių variacija Hsp90",alpha,"N be ligandų")))+guides(color=guide_legend(title="Duomenys"))+stat_summary(fun.data=n_fun,geom="text",colour="blue",size=4)
dev.off()
cairo_pdf('R_P_1UYL_sulyginimas.pdf',family="DejaVu Sans", width=10)
ggplot(data=comparison,aes(x=factor(R.Number),y=P.Volume))+geom_boxplot()+geom_point(aes(x=factor(R.Number),y=R.Volume,shape="1UYL"),color="red")+xlab("Ertmės indeksas kontrolinėje struktūroje 1UYL") + ylab(expression(paste("Tūris, A"^3)))+ggtitle(expression(paste("Ertmių tūrių variacija Hsp90",alpha,"N be ligandų")))+guides(color=guide_legend(title="Duomenys"))+stat_summary(fun.data=n_fun,geom="text",colour="blue",size=4)
dev.off()
ggplot(data=comparison,aes(x=factor(R.Number),y=P.Volume))+geom_boxplot()+geom_point(aes(x=factor(R.Number),y=R.Volume,shape="1UYL"),color="red")+xlab("Ertmės indeksas kontrolinėje struktūroje 1UYL") + ylab(expression(paste("Tūris, A"^3)))+ggtitle(expression(paste("Ertmių tūrių variacija Hsp90",alpha,"N be ligandų")))+guides(color=guide_legend(title="Duomenys"))+stat_summary(fun.data=n_fun,geom="text",colour="blue",size=4)+theme(legend.title=element_blank())
cairo_pdf('R_P_1UYL_sulyginimas.pdf',family="DejaVu Sans", width=10)
ggplot(data=comparison,aes(x=factor(R.Number),y=P.Volume))+geom_boxplot()+geom_point(aes(x=factor(R.Number),y=R.Volume,shape="1UYL"),color="red")+xlab("Ertmės indeksas kontrolinėje struktūroje 1UYL") + ylab(expression(paste("Tūris, A"^3)))+ggtitle(expression(paste("Ertmių tūrių variacija Hsp90",alpha,"N be ligandų")))+guides(color=guide_legend(title="Duomenys"))+stat_summary(fun.data=n_fun,geom="text",colour="blue",size=4)+theme(legend.title=element_blank())
dev.off()
ggplot(data=RPL,aes(x=factor(R.Number),y=PL.Volume))+geom_boxplot()+geom_point(aes(x=factor(R.Number),y=R.Volume,shape="1UYL"),color="red")+xlab("Ertmės indeksas kontrolinėje struktūroje 1UYL") + ylab(expression(paste("Tūris, A"^3)))+ggtitle(expression(paste("Ertmių tūrių variacija Hsp90",alpha,"N su ligandais")))+guides(color=guide_legend(title="Duomenys"))+stat_summary(fun.data=n_fun,geom="text",colour="blue",size=4)+theme(legend.title=element_blank())
cairo_pdf('R_PL_1UYL_sulyginimas.pdf',family="DejaVu Sans", width=10)
ggplot(data=RPL,aes(x=factor(R.Number),y=PL.Volume))+geom_boxplot()+geom_point(aes(x=factor(R.Number),y=R.Volume,shape="1UYL"),color="red")+xlab("Ertmės indeksas kontrolinėje struktūroje 1UYL") + ylab(expression(paste("Tūris, A"^3)))+ggtitle(expression(paste("Ertmių tūrių variacija Hsp90",alpha,"N su ligandais")))+guides(color=guide_legend(title="Duomenys"))+stat_summary(fun.data=n_fun,geom="text",colour="blue",size=4)+theme(legend.title=element_blank())
dev.off()
ggplot(data=RPL,aes(x=factor(R.Number),y=PL.Volume))+geom_boxplot()+geom_point(aes(x=factor(R.Number),y=R.Volume,shape="1UYL"),color="red")+xlab("Ertmės indeksas kontrolinėje struktūroje 1UYL") + ylab(expression(paste("Tūris, A"^3)))+ggtitle(expression(paste("Ertmių tūrių variacija Hsp90",alpha,"N su ligandais")))+guides(color=guide_legend(title="Duomenys"))+stat_summary(fun.data=n_fun,geom="text",colour="blue",size=4,shape="N")+theme(legend.title=element_blank())
ggplot(data=RPL,aes(x=factor(R.Number),y=PL.Volume))+geom_boxplot()+geom_point(aes(x=factor(R.Number),y=R.Volume,shape="1UYL"),color="red")+xlab("Ertmės indeksas kontrolinėje struktūroje 1UYL") + ylab(expression(paste("Tūris, A"^3)))+ggtitle(expression(paste("Ertmių tūrių variacija Hsp90",alpha,"N su ligandais")))+guides(color=guide_legend(title="Duomenys"))+stat_summary(fun.data=n_fun,geom="text",colour="blue",size=4)+theme(legend.title=element_blank())
ggplot(data=RPL,aes(x=factor(R.Number),y=PL.Volume))+geom_boxplot()+geom_point(aes(x=factor(R.Number),y=R.Volume,shape="1UYL"),color="red")+xlab("Ertmės indeksas kontrolinėje struktūroje 1UYL") + ylab(expression(paste("Tūris, A"^3)))+ggtitle(expression(paste("Ertmių tūrių variacija Hsp90",alpha,"N su ligandais")))+guides(color=guide_legend(title="Duomenys"))+stat_summary(fun.data=n_fun,geom="text",colour="blue",size=4)+theme(legend.title=element_blank())+annotate("text",x="36",y=-5,label="N=")
ggplot(data=RPL,aes(x=factor(R.Number),y=PL.Volume))+geom_boxplot()+geom_point(aes(x=factor(R.Number),y=R.Volume,shape="1UYL"),color="red")+xlab("Ertmės indeksas kontrolinėje struktūroje 1UYL") + ylab(expression(paste("Tūris, A"^3)))+ggtitle(expression(paste("Ertmių tūrių variacija Hsp90",alpha,"N su ligandais")))+guides(color=guide_legend(title="Duomenys"))+stat_summary(fun.data=n_fun,geom="text",colour="blue",size=4)+theme(legend.title=element_blank())+annotate("text",x="1",y=-10,label="N=")
ggplot(data=RPL,aes(x=factor(R.Number),y=PL.Volume))+geom_boxplot()+geom_point(aes(x=factor(R.Number),y=R.Volume,shape="1UYL"),color="red")+xlab("Ertmės indeksas kontrolinėje struktūroje 1UYL") + ylab(expression(paste("Tūris, A"^3)))+ggtitle(expression(paste("Ertmių tūrių variacija Hsp90",alpha,"N su ligandais")))+guides(color=guide_legend(title="Duomenys"))+stat_summary(fun.data=n_fun,geom="text",colour="blue",size=4)+theme(legend.title=element_blank())+annotate("text",x="1",y=-10,label="N=",color="blue")
ggplot(data=RPL,aes(x=factor(R.Number),y=PL.Volume))+geom_boxplot()+geom_point(aes(x=factor(R.Number),y=R.Volume,shape="1UYL"),color="red")+xlab("Ertmės indeksas kontrolinėje struktūroje 1UYL") + ylab(expression(paste("Tūris, A"^3)))+ggtitle(expression(paste("Ertmių tūrių variacija Hsp90",alpha,"N su ligandais")))+guides(color=guide_legend(title="Duomenys"))+stat_summary(fun.data=n_fun,geom="text",colour="blue",size=4)+theme(legend.title=element_blank())+annotate("text",x="1",y=-10,label="N=",color="blue", size="5")
ggplot(data=RPL,aes(x=factor(R.Number),y=PL.Volume))+geom_boxplot()+geom_point(aes(x=factor(R.Number),y=R.Volume,shape="1UYL"),color="red")+xlab("Ertmės indeksas kontrolinėje struktūroje 1UYL") + ylab(expression(paste("Tūris, A"^3)))+ggtitle(expression(paste("Ertmių tūrių variacija Hsp90",alpha,"N su ligandais")))+guides(color=guide_legend(title="Duomenys"))+stat_summary(fun.data=n_fun,geom="text",colour="blue",size=4)+theme(legend.title=element_blank())+annotate("text",x="1",y=-10,label="N=",color="blue", size=5)
ggplot(data=RPL,aes(x=factor(R.Number),y=PL.Volume))+geom_boxplot()+geom_point(aes(x=factor(R.Number),y=R.Volume,shape="1UYL"),color="red")+xlab("Ertmės indeksas kontrolinėje struktūroje 1UYL") + ylab(expression(paste("Tūris, A"^3)))+ggtitle(expression(paste("Ertmių tūrių variacija Hsp90",alpha,"N su ligandais")))+guides(color=guide_legend(title="Duomenys"))+stat_summary(fun.data=n_fun,geom="text",colour="blue",size=4)+theme(legend.title=element_blank())+annotate("text",x="1",y=-20,label="N=",color="blue", size=4)
ggplot(data=RPL,aes(x=factor(R.Number),y=PL.Volume))+geom_boxplot()+geom_point(aes(x=factor(R.Number),y=R.Volume,shape="1UYL"),color="red")+xlab("Ertmės indeksas kontrolinėje struktūroje 1UYL") + ylab(expression(paste("Tūris, A"^3)))+ggtitle(expression(paste("Ertmių tūrių variacija Hsp90",alpha,"N su ligandais")))+guides(color=guide_legend(title="Duomenys"))+stat_summary(fun.data=n_fun,geom="text",colour="blue",size=4)+theme(legend.title=element_blank())+annotate("text",x="1"+1,y=-20,label="N=",color="blue", size=4)
ggplot(data=RPL,aes(x=factor(R.Number),y=PL.Volume))+geom_boxplot()+geom_point(aes(x=factor(R.Number),y=R.Volume,shape="1UYL"),color="red")+xlab("Ertmės indeksas kontrolinėje struktūroje 1UYL") + ylab(expression(paste("Tūris, A"^3)))+ggtitle(expression(paste("Ertmių tūrių variacija Hsp90",alpha,"N su ligandais")))+guides(color=guide_legend(title="Duomenys"))+stat_summary(fun.data=n_fun,geom="text",colour="blue",size=4)+theme(legend.title=element_blank())+annotate("text",x=37,y=-20,label="N=",color="blue", size=4)
ggplot(data=RPL,aes(x=factor(R.Number),y=PL.Volume))+geom_boxplot()+geom_point(aes(x=factor(R.Number),y=R.Volume,shape="1UYL"),color="red")+xlab("Ertmės indeksas kontrolinėje struktūroje 1UYL") + ylab(expression(paste("Tūris, A"^3)))+ggtitle(expression(paste("Ertmių tūrių variacija Hsp90",alpha,"N su ligandais")))+guides(color=guide_legend(title="Duomenys"))+stat_summary(fun.data=n_fun,geom="text",colour="blue",size=4)+theme(legend.title=element_blank())+annotate("text",x=11,y=-20,label="N=",color="blue", size=4)
ggplot(data=RPL,aes(x=factor(R.Number),y=PL.Volume))+geom_boxplot()+geom_point(aes(x=factor(R.Number),y=R.Volume,shape="1UYL"),color="red")+xlab("Ertmės indeksas kontrolinėje struktūroje 1UYL") + ylab(expression(paste("Tūris, A"^3)))+ggtitle(expression(paste("Ertmių tūrių variacija Hsp90",alpha,"N su ligandais")))+guides(color=guide_legend(title="Duomenys"))+stat_summary(fun.data=n_fun,geom="text",colour="blue",size=4)+theme(legend.title=element_blank())+annotate("text",x=15,y=-20,label="N=",color="blue", size=4)
ggplot(data=RPL,aes(x=factor(R.Number),y=PL.Volume))+geom_boxplot()+geom_point(aes(x=factor(R.Number),y=R.Volume,shape="1UYL"),color="red")+xlab("Ertmės indeksas kontrolinėje struktūroje 1UYL") + ylab(expression(paste("Tūris, A"^3)))+ggtitle(expression(paste("Ertmių tūrių variacija Hsp90",alpha,"N su ligandais")))+guides(color=guide_legend(title="Duomenys"))+stat_summary(fun.data=n_fun,geom="text",colour="blue",size=4)+theme(legend.title=element_blank())+annotate("text",x=118,y=-20,label="N=",color="blue", size=4)
ggplot(data=RPL,aes(x=factor(R.Number),y=PL.Volume))+geom_boxplot()+geom_point(aes(x=factor(R.Number),y=R.Volume,shape="1UYL"),color="red")+xlab("Ertmės indeksas kontrolinėje struktūroje 1UYL") + ylab(expression(paste("Tūris, A"^3)))+ggtitle(expression(paste("Ertmių tūrių variacija Hsp90",alpha,"N su ligandais")))+guides(color=guide_legend(title="Duomenys"))+stat_summary(fun.data=n_fun,geom="text",colour="blue",size=4)+theme(legend.title=element_blank())+annotate("text",x=117,y=-20,label="N=",color="blue", size=4)
ggplot(data=RPL,aes(x=factor(R.Number),y=PL.Volume))+geom_boxplot()+geom_point(aes(x=factor(R.Number),y=R.Volume,shape="1UYL"),color="red")+xlab("Ertmės indeksas kontrolinėje struktūroje 1UYL") + ylab(expression(paste("Tūris, A"^3)))+ggtitle(expression(paste("Ertmių tūrių variacija Hsp90",alpha,"N su ligandais")))+guides(color=guide_legend(title="Duomenys"))+stat_summary(fun.data=n_fun,geom="text",colour="blue",size=4)+theme(legend.title=element_blank())+annotate("text",x=116,y=-20,label="N=",color="blue", size=4)
ggplot(data=RPL,aes(x=factor(R.Number),y=PL.Volume))+geom_boxplot()+geom_point(aes(x=factor(R.Number),y=R.Volume,shape="1UYL"),color="red")+xlab("Ertmės indeksas kontrolinėje struktūroje 1UYL") + ylab(expression(paste("Tūris, A"^3)))+ggtitle(expression(paste("Ertmių tūrių variacija Hsp90",alpha,"N su ligandais")))+guides(color=guide_legend(title="Duomenys"))+stat_summary(fun.data=n_fun,geom="text",colour="blue",size=4)+theme(legend.title=element_blank())+annotate("text",x=17,y=-20,label="N=",color="blue", size=4)
ggplot(data=RPL,aes(x=factor(R.Number),y=PL.Volume))+geom_boxplot()+geom_point(aes(x=factor(R.Number),y=R.Volume,shape="1UYL"),color="red")+xlab("Ertmės indeksas kontrolinėje struktūroje 1UYL") + ylab(expression(paste("Tūris, A"^3)))+ggtitle(expression(paste("Ertmių tūrių variacija Hsp90",alpha,"N su ligandais")))+guides(color=guide_legend(title="Duomenys"))+stat_summary(fun.data=n_fun,geom="text",colour="blue",size=4)+theme(legend.title=element_blank())+annotate("text",x=17,y=-10,label="N=",color="blue", size=4)
ggplot(data=RPL,aes(x=factor(R.Number),y=PL.Volume))+geom_boxplot()+geom_point(aes(x=factor(R.Number),y=R.Volume,shape="1UYL"),color="red")+xlab("Ertmės indeksas kontrolinėje struktūroje 1UYL") + ylab(expression(paste("Tūris, A"^3)))+ggtitle(expression(paste("Ertmių tūrių variacija Hsp90",alpha,"N su ligandais")))+guides(color=guide_legend(title="Duomenys"))+stat_summary(fun.data=n_fun,geom="text",colour="blue",size=4)+theme(legend.title=element_blank())+annotate("text",x=-1,y=-10,label="N=",color="blue", size=4)
ggplot(data=RPL,aes(x=factor(R.Number),y=PL.Volume))+geom_boxplot()+geom_point(aes(x=factor(R.Number),y=R.Volume,shape="1UYL"),color="red")+xlab("Ertmės indeksas kontrolinėje struktūroje 1UYL") + ylab(expression(paste("Tūris, A"^3)))+ggtitle(expression(paste("Ertmių tūrių variacija Hsp90",alpha,"N su ligandais")))+guides(color=guide_legend(title="Duomenys"))+stat_summary(fun.data=n_fun,geom="text",colour="blue",size=4)+theme(legend.title=element_blank())+annotate("text",x=0,y=-10,label="N=",color="blue", size=4)
ggplot(data=RPL,aes(x=factor(R.Number),y=PL.Volume))+geom_boxplot()+geom_point(aes(x=factor(R.Number),y=R.Volume,shape="1UYL"),color="red")+xlab("Ertmės indeksas kontrolinėje struktūroje 1UYL") + ylab(expression(paste("Tūris, A"^3)))+ggtitle(expression(paste("Ertmių tūrių variacija Hsp90",alpha,"N su ligandais")))+guides(color=guide_legend(title="Duomenys"))+stat_summary(fun.data=n_fun,geom="text",colour="blue",size=4)+theme(legend.title=element_blank())+annotate("text",x=18,y=-10,label="N=",color="blue", size=4)
ggplot(data=RPL,aes(x=factor(R.Number),y=PL.Volume))+geom_boxplot()+geom_point(aes(x=factor(R.Number),y=R.Volume,shape="1UYL"),color="red")+xlab("Ertmės indeksas kontrolinėje struktūroje 1UYL") + ylab(expression(paste("Tūris, A"^3)))+ggtitle(expression(paste("Ertmių tūrių variacija Hsp90",alpha,"N su ligandais")))+guides(color=guide_legend(title="Duomenys"))+stat_summary(fun.data=n_fun,geom="text",colour="blue",size=4)+theme(legend.title=element_blank())+annotate("text",x=17,y=-10,label="N=",color="blue", size=4)
ggplot(data=RPL,aes(x=factor(R.Number),y=PL.Volume))+geom_boxplot()+geom_point(aes(x=factor(R.Number),y=R.Volume,shape="1UYL"),color="red")+xlab("Ertmės indeksas kontrolinėje struktūroje 1UYL") + ylab(expression(paste("Tūris, A"^3)))+ggtitle(expression(paste("Ertmių tūrių variacija Hsp90",alpha,"N su ligandais")))+guides(color=guide_legend(title="Duomenys"))+stat_summary(fun.data=n_fun,geom="text",colour="blue",size=4)+theme(legend.title=element_blank())+annotate("text",x=17,y=-5,label="=N",color="blue", size=4)
ggplot(data=RPL,aes(x=factor(R.Number),y=PL.Volume))+geom_boxplot()+geom_point(aes(x=factor(R.Number),y=R.Volume,shape="1UYL"),color="red")+xlab("Ertmės indeksas kontrolinėje struktūroje 1UYL") + ylab(expression(paste("Tūris, A"^3)))+ggtitle(expression(paste("Ertmių tūrių variacija Hsp90",alpha,"N su ligandais")))+guides(color=guide_legend(title="Duomenys"))+stat_summary(fun.data=n_fun,geom="text",colour="blue",size=4)+theme(legend.title=element_blank())+annotate("text",x=17,y=-10,label="N=",color="blue", size=4)
ggplot(data=RPL,aes(x=factor(R.Number),y=PL.Volume))+geom_boxplot()+geom_point(aes(x=factor(R.Number),y=R.Volume,shape="1UYL"),color="red")+xlab("Ertmės indeksas kontrolinėje struktūroje 1UYL") + ylab(expression(paste("Tūris, A"^3)))+ggtitle(expression(paste("Ertmių tūrių variacija Hsp90",alpha,"N su ligandais")))+guides(color=guide_legend(title="Duomenys"))+stat_summary(fun.data=n_fun,geom="text",colour="blue",size=4)+theme(legend.title=element_blank())
p<-ggplot(data=RPL,aes(x=factor(R.Number),y=PL.Volume))+geom_boxplot()+geom_point(aes(x=factor(R.Number),y=R.Volume,shape="1UYL"),color="red")+xlab("Ertmės indeksas kontrolinėje struktūroje 1UYL") + ylab(expression(paste("Tūris, A"^3)))+ggtitle(expression(paste("Ertmių tūrių variacija Hsp90",alpha,"N su ligandais")))+guides(color=guide_legend(title="Duomenys"))+stat_summary(fun.data=n_fun,geom="text",colour="blue",size=4)+theme(legend.title=element_blank())
p
cairo_pdf('R_P_1UYL_sulyginimas.pdf',family="DejaVu Sans", width=10)
ggplot(data=comparison,aes(x=factor(R.Number),y=P.Volume))+geom_boxplot()+geom_point(aes(x=factor(R.Number),y=R.Volume,shape="1UYL"),color="red")+xlab("Ertmės indeksas kontrolinėje struktūroje 1UYL") + ylab(expression(paste("Tūris, A"^3)))+ggtitle(expression(paste("Ertmių tūrių variacija Hsp90",alpha,"N be ligandų")))+guides(color=guide_legend(title="Duomenys"))+stat_summary(fun.data=n_fun,geom="text",colour="blue",size=4)+theme(legend.title=element_blank())
dev.off()
volumes<-read.csv("Summary.csv")
volumes<-read.csv("Summary.csv",header=TRUE)
volumes
colnames(volumes)
volumes$Fragment
volumes$Protein.segment
volumes$Fragment
P_volumes<-volumes[which(volumes$Fragment=="P"),c("File","Fragment","Protein.volume..A.3","Solvent.accessible.volume..A.3","VdW.volume..A.3","Cavity.volume..A.3","Cleft.volume..A.3")]
P_volumes
PL_volumes<-volumes[which(volumes$Fragment=="PL"),c("File","Fragment","Protein.volume..A.3","Solvent.accessible.volume..A.3","VdW.volume..A.3","Cavity.volume..A.3","Cleft.volume..A.3")]
P_vol<-ggpot(P_volumes,aes(P_volumes$Protein.volume..A.3)) +geom_histogram()
P_vol<-ggplot(P_volumes,aes(P_volumes$Protein.volume..A.3)) +geom_histogram()
P_vol
P_vol<-ggplot(P_volumes,aes(P_volumes$Protein.volume..A.3)) +geom_histogram(binwidth=50)
P_vol
P_vol<-ggplot(P_volumes,aes(P_volumes$Protein.volume..A.3)) +geom_histogram(binwidth=100)
P_vol
P_vol+xlab(expression(paste("Baltymo tūris, A",^,"3")))
P_vol+xlab(expression(paste("Baltymo tūris, A",^3)))
P_vol+xlab(expression(paste("Baltymo tūris, A",^3)))
P_vol+xlab( expression(paste("Baltymo tūris, A"^3)) )
P_vol+xlab( expression(paste("Baltymo tūris,"ring(A)^3)) )
P_vol+xlab( expression(paste("Baltymo tūris,", ring(A)^3)) )
P_vol+xlab( expression(paste("Baltymo tūris, ", ring(A)^3)) )
P_vol+geom_histogram(binwidth=100,alpha=I(.2))
P_vol+geom_histogram(binwidth=100,alpha=I(.2),fill=I("blue"))
P_vol<-ggplot(P_volumes,aes(P_volumes$Protein.volume..A.3)) +geom_histogram(binwidth=50,fill=I("blue"))
P_vol
P_vol<-geom_histogram(bandwidth=50, fill=I("green"),alpha=I(0.2))
P_vol
P_vol+geom_histogram(bandwidth=50, fill=I("green"),alpha=I(0.2))
P_vol<-ggplot(P_volumes,aes(P_volumes$Protein.volume..A.3))
P_vol+geom_histogram(bandwidth=50, fill=I("green"),alpha=I(0.2))
P_vol+geom_histogram(bandwidth=50, fill=I("blue"),alpha=I(0.2))
P_vol+geom_histogram(bandwidth=50, fill=I("blue"),alpha=I(1))
P_vol+geom_histogram(bandwidth=50, fill=I("blue"),alpha=I(0.5))
P_vol+geom_histogram(bandwidth=50, fill=I("blue"),alpha=I(0.5),shape="P")
P_vol+geom_histogram(binwidth=100, fill=I("blue"),alpha=I(0.5),shape="P")
P_vol+ggplot(PL_volumes,aes(PL_volumes$Protein.volume..A.3))+geom_histogram(binwidth=100,fill="green")
P_vol+geom_histogram(data=PL_volumes,aes(PL_volumes$Protein.volume..A.3),binwidth=100,fill=I("green"),alpha=I(0.5))
P_vol+geom_histogram(data=PL_volumes,aes(PL_volumes$Protein.volume..A.3),binwidth=100,fill=I("green"),alpha=I(0.5))
P_vol+geom_histogram(data=PL_volumes,aes(PL_volumes$Protein.volume..A.3),binwidth=100,fill="green",alpha=0.5)+geom_histogram(binwidth=100,data=P_volumes,aes(P_volumes$Protein.volume..A.3),fill="blue",alpha=0.5)
P_vol+geom_histogram(data=PL_volumes,aes(PL_volumes$Protein.volume..A.3),binwidth=100,fill=I("green"),alpha=I(0.5))
P_vol+geom_histogram(data=PL_volumes,aes(PL_volumes$Protein.volume..A.3),binwidth=100,fill="green",alpha=0.5)+geom_histogram(binwidth=100,data=P_volumes,aes(P_volumes$Protein.volume..A.3),fill="blue",alpha=0.5)
P_vol+guides(fill=guide_legend("Tūriai"))
P_vol+scale_fill_manual(name="Tūriai",values=c("b"="blue","g"="green"),labels=c("b"="P","g"="PL"))
P_vol
P_vol+geom_histogram(data=PL_volumes,aes(PL_volumes$Protein.volume..A.3),binwidth=100,fill="green",alpha=0.5)+geom_histogram(binwidth=100,data=P_volumes,aes(P_volumes$Protein.volume..A.3),fill="blue",alpha=0.5)
P_vol
P_vol+geom_histogram(data=PL_volumes,aes(PL_volumes$Protein.volume..A.3),binwidth=100,fill="green",alpha=0.5)+geom_histogram(binwidth=100,data=P_volumes,aes(P_volumes$Protein.volume..A.3),fill="blue",alpha=0.5)
P_vol<-geom_histogram(data=PL_volumes,aes(PL_volumes$Protein.volume..A.3),binwidth=100,fill="green",alpha=0.5)+geom_histogram(binwidth=100,data=P_volumes,aes(P_volumes$Protein.volume..A.3),fill="blue",alpha=0.5)
P_vol<-ggplot(),geom_histogram(data=PL_volumes,aes(PL_volumes$Protein.volume..A.3),binwidth=100,fill="green",alpha=0.5)+geom_histogram(binwidth=100,data=P_volumes,aes(P_volumes$Protein.volume..A.3),fill="blue",alpha=0.5)
P_vol<-ggplot()+geom_histogram(data=PL_volumes,aes(PL_volumes$Protein.volume..A.3),binwidth=100,fill="green",alpha=0.5)+geom_histogram(binwidth=100,data=P_volumes,aes(P_volumes$Protein.volume..A.3),fill="blue",alpha=0.5)
P_vol
P_vol+guides(fill=guide_legend("Tūris"))
P_vol
P_vol+scale_fill_manual(name="Tūriai",values=c("b"="blue","g"="green"),labels=c("b"="P","g"="PL"))
P_vol
P_vol+scale_fill_manual(name="Tūriai",values=c("b"="blue","g"="green"),labels=c("b"="P","g"="PL"))
P_vol
P_vol<-ggplot()+geom_histogram(data=PL_volumes,aes(PL_volumes$Protein.volume..A.3),binwidth=100,fill="green",colour="red",alpha=0.5)+geom_histogram(binwidth=100,data=P_volumes,aes(P_volumes$Protein.volume..A.3),fill="blue",colour="red",alpha=0.5)
P_vol
P_vol<-ggplot()+geom_histogram(data=PL_volumes,aes(PL_volumes$Protein.volume..A.3),binwidth=100,fill="green",colour="black",alpha=0.5)+geom_histogram(binwidth=100,data=P_volumes,aes(P_volumes$Protein.volume..A.3),fill="blue",colour="black",alpha=0.5)
P_vol
P_vol<-ggplot()+geom_histogram(data=PL_volumes,aes(PL_volumes$Protein.volume..A.3),binwidth=100,fill="green",colour="black",alpha=0.5,xName="P")+geom_histogram(binwidth=100,data=P_volumes,aes(P_volumes$Protein.volume..A.3),fill="blue",colour="black",alpha=0.5,xName="PL")
P_vol
P_vol<-ggplot()+geom_histogram(data=PL_volumes,aes(PL_volumes$Protein.volume..A.3),binwidth=100,fill="green",colour="black",alpha=0.5,xName="Baltymo tūris")+geom_histogram(binwidth=100,data=P_volumes,aes(P_volumes$Protein.volume..A.3),fill="blue",colour="black",alpha=0.5,xName="Baltymo tūris")
PPL<-volumes[which(volumes$Fragment=="P" | volumes$Fragment=="PL"),]
PPL
P_vol<-ggplot(data=PPL,aes(PPL$Protein.volume..A.3,fill=PPL$Fragment))+geom_histogram(binwidth=100,alpha=0.5)
P_vol
P_vol+theme(legend(title="Tūrio tipas"))
P_vol+theme(guide_legend(title="Tūrio tipas"))
P_volŽtheme(legend.title="Tūrio tipas")
P_vol+theme(legend.title="Tūrio tipas")
P_vol+labs(aesthetic="Tūrio tipas")
P_vol
P_vol+xlab("Baltymo tūris")
P_vol+xlab(expression(paste("Baltymo tūris"^3)))
P_vol+xlab(expression(paste("Baltymo tūris, "ring(A)^3)))
P_vol+xlab(expression(paste("Baltymo tūris,", ring(A)^3)))
P_vol+xlab(expression(paste("Baltymo tūris, ", ring(A)^3)))
P_vol+labs(fill="Tūrio tipas")
P_vol<-P_vol+labs(fill="Tūrio tipas")
P_vol
P_vol<-P_vol+xlab(expression(paste("Baltymo tūris, ", ring(A)^3)))
P_vol
P_vol+ylab("Skaičius")
P_vol2<-ggplot()+geom_bar(data=PL_volumes,aes(PL_volumes$Protein.volume..A.3),binwidth=100,fill="green",colour="black",alpha=0.2,)+geom_bar(binwidth=100,data=P_volumes,aes(P_volumes$Protein.volume..A.3),fill="blue",colour="black",alpha=0.2,xName="Baltymo tūris")
P_vol<-ggplot(data=PPL,aes(PPL$Protein.volume..A.3,fill=PPL$Fragment))+geom_histogram(binwidth=100,alpha=0.2)+xlab(expression(paste("Baltymo tūris, ",ring(A)Š3)))+ylab("Skaičius")+labs(fill="Tūrio tipas")+ggtitle("Baltymo tūrio variacija")
P_vol<-ggplot(data=PPL,aes(PPL$Protein.volume..A.3,fill=PPL$Fragment))+geom_histogram(binwidth=100,alpha=0.2)+xlab(expression(paste("Baltymo tūris, ",ring(A)^3)))+ylab("Skaičius")+labs(fill="Tūrio tipas")+ggtitle("Baltymo tūrio variacija")
P_vol
P_vol<-ggplot(data=PPL,aes(PPL$Protein.volume..A.3,fill=PPL$Fragment))+geom_histogram(binwidth=100,alpha=0.4,colour="black")+xlab(expression(paste("Baltymo tūris, ",ring(A)^3)))+ylab("Skaičius")+labs(fill="Tūrio tipas")+ggtitle("Baltymo tūrio variacija")
P_vol
P_vol<-ggplot(data=PPL,aes(PPL$Protein.volume..A.3,fill=PPL$Fragment))+geom_histogram(binwidth=100,alpha=I(0.2),colour="black")+xlab(expression(paste("Baltymo tūris, ",ring(A)^3)))+ylab("Skaičius")+labs(fill="Tūrio tipas")+ggtitle("Baltymo tūrio variacija")
P_vol2
P_vol<-ggplot(data=PPL,aes(PPL$Protein.volume..A.3,fill=PPL$Fragment))+geom_histogram(binwidth=100,alpha=I(0.2),colour="black")+xlab(expression(paste("Baltymo tūris, ",ring(A)^3)))+ylab("Skaičius")+labs(fill="Tūrio tipas")+ggtitle("Baltymo tūrio variacija")
P_vol
P_vol<-ggplot(data=PPL,aes(PPL$Protein.volume..A.3,fill=PPL$Fragment))+geom_density(alpha=I(0.2))+xlab(expression(paste("Baltymo tūris, ",ring(A)^3)))+ylab("Skaičius")+labs(fill="Tūrio tipas")+ggtitle("Baltymo tūrio variacija")
P_vol
P_vol<-ggplot(data=PPL,aes(PPL$Protein.volume..A.3,fill=PPL$Fragment))+geom_density(alpha=I(0.2),aes(y=..density..))+xlab(expression(paste("Baltymo tūris, ",ring(A)^3)))+ylab("Skaičius")+labs(fill="Tūrio tipas")+ggtitle("Baltymo tūrio variacija")
P_vol
P_vol<-ggplot(data=PPL,aes(PPL$Protein.volume..A.3,fill=PPL$Fragment))+geom_histogram(alpha=I(0.2),aes(y=..density..),binwidth=100)+xlab(expression(paste("Baltymo tūris, ",ring(A)^3)))+ylab("Skaičius")+labs(fill="Tūrio tipas")+ggtitle("Baltymo tūrio variacija")
P_vol
P_vol<-ggplot(data=PPL,aes(PPL$Protein.volume..A.3,fill=PPL$Fragment))+geom_histogram(alpha=0.5,aes(y=..density..),binwidth=100,position='identity')+xlab(expression(paste("Baltymo tūris, ",ring(A)^3)))+ylab("Skaičius")+labs(fill="Tūrio tipas")+ggtitle("Baltymo tūrio variacija")
P_vol
P_vol<-ggplot(data=PPL,aes(PPL$Protein.volume..A.3,fill=PPL$Fragment))+geom_histogram(alpha=0.5,binwidth=100,position='identity')+xlab(expression(paste("Baltymo tūris, ",ring(A)^3)))+ylab("Skaičius")+labs(fill="Tūrio tipas")+ggtitle("Baltymo tūrio variacija")
P_vol
cairo_pdf('Protein_vol.pdf',family="DejaVu Sans", width=10)
dev.off()
P_VdW<-ggplot(data=PPL,aes(PPL$VdW.volume..A.3,fill=PPL$Fragment))+geom_histogram(alpha=0.5,binwidth=100,position='identity')+xlab(expression(paste("Tūris, ",ring(A)^3)))+ylab("Skaičius")+labs(fill="Tūrio tipas")+ggtitle("Baltymo VdW tūrio variacija")
P_VdW
cairo_pdf('Protein_VdW_vol.pdf',family="DejaVu Sans", width=10)
P_VdW
dev.off()
cairo_pdf('Protein_vol.pdf',family="DejaVu Sans", width=10)
P_vol<-ggplot(data=PPL,aes(PPL$Protein.volume..A.3,fill=PPL$Fragment))+geom_histogram(alpha=0.5,binwidth=100,position='identity')+xlab(expression(paste("Tūris, ",ring(A)^3)))+ylab("Skaičius")+labs(fill="Tūrio tipas")+ggtitle("Baltymo tūrio variacija")
P_vol
dev.off()
P_Cav_vol<-ggplot(data=PPL,aes(PPL$Cavity.volume..A.3,fill=PPL$Fragment))+geom_histogram(alpha=0.5,binwidth=100,position='identity')+xlab(expression(paste("Tūris, ",ring(A)^3)))+ylab("Skaičius")+labs(fill="Tūrio tipas")+ggtitle("Cavity tūrio variacija")
P_Cav_vol
P_Cav_vol<-ggplot(data=PPL,aes(PPL$Cavity.volume..A.3,fill=PPL$Fragment))+geom_histogram(alpha=0.5,binwidth=20,position='identity')+xlab(expression(paste("Tūris, ",ring(A)^3)))+ylab("Skaičius")+labs(fill="Tūrio tipas")+ggtitle("Cavity tūrio variacija")
P_Cav_vol
P_Cav_vol<-ggplot(data=PPL,aes(PPL$Cavity.volume..A.3,fill=PPL$Fragment))+geom_histogram(alpha=0.5,binwidth=10,position='identity')+xlab(expression(paste("Tūris, ",ring(A)^3)))+ylab("Skaičius")+labs(fill="Tūrio tipas")+ggtitle("Cavity tūrio variacija")
P_Cav_vol
cairo_pdf('Protein_Cav_vol.pdf',family="DejaVu Sans", width=10)
dev.off
dev.off()
P_Clef_vol<-ggplot(data=PPL,aes(PPL$Cleft.volume..A.3,fill=PPL$Fragment))+geom_histogram(alpha=0.5,binwidth=10,position='identity')+xlab(expression(paste("Tūris, ",ring(A)^3)))+ylab("Skaičius")+labs(fill="Tūrio tipas")+ggtitle("Cleft tūrio variacija")
P_Clef_vol
P_Clef_vol<-ggplot(data=PPL,aes(PPL$Cleft.volume..A.3,fill=PPL$Fragment))+geom_histogram(alpha=0.5,binwidth=200,position='identity')+xlab(expression(paste("Tūris, ",ring(A)^3)))+ylab("Skaičius")+labs(fill="Tūrio tipas")+ggtitle("Cleft tūrio variacija")
P_Clef_vol<-ggplot(data=PPL,aes(PPL$Cleft.volume..A.3,fill=PPL$Fragment))+geom_histogram(alpha=0.5,binwidth=20,position='identity')+xlab(expression(paste("Tūris, ",ring(A)^3)))+ylab("Skaičius")+labs(fill="Tūrio tipas")+ggtitle("Cleft tūrio variacija")
P_Clef_vol
cairo_pdf('Protein_Clef_vol.pdf',family="DejaVu Sans", width=10)
P_Clef_vol
dev.off()
P_Solv_vol<-ggplot(data=PPL,aes(PPL$Solvent.accessible.volume..A.3,fill=PPL$Fragment))+geom_histogram(alpha=0.5,binwidth=50,position='identity')+xlab(expression(paste("Tūris, ",ring(A)^3)))+ylab("Skaičius")+labs(fill="Tūrio tipas")+ggtitle("Tirpikliui prieinamo tūrio variacija")
P_Solv_vol
P_Solv_vol<-ggplot(data=PPL,aes(PPL$Solvent.accessible.volume..A.3,fill=PPL$Fragment))+geom_histogram(alpha=0.5,binwidth=100,position='identity')+xlab(expression(paste("Tūris, ",ring(A)^3)))+ylab("Skaičius")+labs(fill="Tūrio tipas")+ggtitle("Tirpikliui prieinamo tūrio variacija")
P_Solv_vol
cairo_pdf('Protein_Solv_vol.pdf',family="DejaVu Sans", width=10)
P_Solv_vol
dev.off()
savehistory('.Rhistory2')
comparison
data
comparison<-data[which(data$Type=="R_P" & data$R.Number!="NA" $ data$P.Number!="NA"),]
comparison<-data[which(data$Type=="R_P" & data$R.Number!="NA" $ data$P.Number!="NA"),]
comparison<-data[which(data$Type=="R_P" & data$R.Number!="NA" $ data$P.Number!="NA"),]
comparison<-data[which(data$Type=="R_P" & data$R.Number!="NA" & data$P.Number!="NA"),]
comparison
Vol_dist<-ggplot(comparison, aes(x=comparison$R.Number,y=comparison$P.Volume,fill=comparison$P.Type))+geom_boxplot()
Vol_dist
Vol_dist<-ggplot(comparison, aes(x=factor(comparison$R.Number,y=comparison$P.Volume,fill=comparison$P.Type))+geom_boxplot()
Vol_dist<-ggplot(comparison, aes(x=factor(comparison$R.Number),y=comparison$P.Volume,fill=comparison$P.Type))+geom_boxplot()
Vol_dist
Vol_dist<-Vol_dist+geom_point(aes(x=factor(R.Number),y=R.Volume)color="red")
Vol_dist<-Vol_dist+geom_point(aes(x=factor(R.Number),y=R.Volume),color="red")
Vol_dist
Vol_dist<-ggplot(comparison, aes(x=factor(comparison$R.Number),y=comparison$PL.Volume,fill=comparison$P.Type))+geom_boxplot()
Vol_dist
Vol_dist+geom_point(aes(x=factor(R.Number),y=R.Volume,colour=R.Type))
Vol_dist+geom_point(aes(x=factor(R.Number),y=R.Volume,colour=R.Type),size=5)
Vol_dist+geom_point(aes(x=factor(R.Number),y=R.Volume,colour=R.Type),size=3)
Vol_dist+geom_point(aes(x=factor(R.Number),y=R.Volume,colour=R.Type),size=3,shape="diamond")
Vol_dist+geom_point(aes(x=factor(R.Number),y=R.Volume,colour=R.Type),size=3,shape=x)
Vol_dist+geom_point(aes(x=factor(R.Number),y=R.Volume,colour=R.Type),size=3,shape="x")
Vol_dist+geom_point(aes(x=factor(R.Number),y=R.Volume,colour=R.Type),size=4,shape="x")
Vol_dist+geom_point(aes(x=factor(R.Number),y=R.Volume,colour=R.Type),size=4,pch=2)
Vol_dist+geom_point(aes(x=factor(R.Number),y=R.Volume,colour=R.Type),size=4,pch=8)
Vol_dist+geom_point(aes(x=factor(R.Number),y=R.Volume,colour=R.Type),size=4,pch=8,cex=1)
Vol_dist+geom_point(aes(x=factor(R.Number),y=R.Volume,colour=R.Type),size=4,pch=8)
Vol_dist+geom_point(aes(x=factor(R.Number),y=R.Volume,colour=R.Type),size=4,pch=8,col="black")
Vol_dist+geom_point(aes(x=factor(R.Number),y=R.Volume,colour=R.Type),size=4,pch=8)
Vol_dist+geom_point(aes(x=factor(R.Number),y=R.Volume,colour=R.Type),size=4,pch=O,col="black")
Vol_dist+geom_point(aes(x=factor(R.Number),y=R.Volume,colour=R.Type),size=4,pch=o,col="black")
Vol_dist+geom_point(aes(x=factor(R.Number),y=R.Volume,colour=R.Type),size=4,pch="o",col="black")
Vol_dist+geom_point(aes(x=factor(R.Number),y=R.Volume,colour=R.Type),size=4,pch="o")
Vol_dist+geom_point(aes(x=factor(R.Number),y=R.Volume,colour=R.Type),size=5,pch="o")
Vol_dist+geom_point(aes(x=factor(R.Number),y=R.Volume,colour=R.Type),size=5,pch="o")+labs(color="Ertmės tipas")
Vol_dist+geom_point(aes(x=factor(R.Number),y=R.Volume,colour=R.Type),size=5,pch="o")+labs(color="1UYL ertmės tipas ir tūris")
Vol_dist+geom_point(aes(x=factor(R.Number),y=R.Volume,colour=R.Type),size=5,pch="o")+labs(color="1UYL ertmės tipas ir tūris")+labs(fill="Ertmės tipas baltymo struktūroje")
Vol_dist+geom_point(aes(x=factor(R.Number),y=R.Volume,colour=R.Type),size=5,pch="o")+labs(color="1UYL ertmės tipas ir tūris")+labs(fill="Ertmės tipas baltymo struktūroje")+xlab("Ertmės indeksas beligandėje 1UYL struktūroje")
Vol_dist+geom_point(aes(x=factor(R.Number),y=R.Volume,colour=R.Type),size=5,pch="o")+labs(color="1UYL ertmės tipas ir tūris")+labs(fill="Ertmės tipas baltymo struktūroje")+xlab("Ertmės indeksas beligandėje 1UYL struktūroje")+ylab(expression(paste("Tūris, ",ring(A)^3)))
Vol_dist+geom_point(aes(x=factor(R.Number),y=R.Volume,colour=R.Type),size=5,pch="o")+labs(color="1UYL ertmės tipas ir tūris")+labs(fill="Ertmės tipas baltymo struktūroje")+xlab("Ertmės indeksas beligandėje 1UYL struktūroje")+ylab(expression(paste("Tūris, ",ring(A)^3)))+ggtitle("Skirtingų ertmių tipų tūrių variacija beligandėse struktūrose")
V_dist<-Vol_dist+labs(fill="Ertmės tipas baltymo struktūroje")+xlab("Ertmės indeksas beligandėje 1UYL struktūroje")+ylab(expression(paste("Tūris, ",ring(A)^3)))+ggtitle("Skirtingų ertmių tipų tūrių variacija beligandėse struktūrose")
V_dist
Vol_dist+geom_point(aes(x=factor(R.Number),y=R.Volume,colour=R.Type),size=5,pch="o")
Vol_dist+geom_point(aes(x=factor(R.Number),y=R.Volume,colour=R.Type),size=5,pch=23,col="black")
Vol_dist+geom_point(aes(x=factor(R.Number),y=R.Volume,colour=R.Type),size=3,pch=23,col="black")
Vol_dist+geom_point(aes(x=factor(R.Number),y=R.Volume,colour=R.Type),size=3,pch=23,bg="black")
Vol_dist+geom_point(aes(x=factor(R.Number),y=R.Volume,colour=R.Type),size=3,pch=23,fg="black")
Vol_dist+geom_point(aes(x=factor(R.Number),y=R.Volume,colour=R.Type),size=3,pch=23)
Vol_dist+geom_point(aes(x=factor(R.Number),y=R.Volume,colour=R.Type),size=3,pch=23,col="black")
Vol_dist+geom_point(aes(x=factor(R.Number),y=R.Volume,colour=R.Type),size=3,pch=1)
Vol_dist+geom_point(aes(x=factor(R.Number),y=R.Volume,colour=R.Type),size=3,pch=23,fg="black")
Vol_dist+geom_point(aes(x=factor(R.Number),y=R.Volume,fill=R.Type),size=3,pch=23,col="black")
Vol_dist+geom_point(aes(x=factor(R.Number),y=R.Volume,fill=R.Type),size=3,pch=23,col="black")+labs(fill="1UYL")
Vol_dist+geom_point(aes(x=factor(R.Number),y=R.Volume,fill=R.Type),size=3,pch=23,col="black")+labs(fill="1UYL")+labs(color="1UYLL")
Vol_dist+geom_point(aes(x=factor(R.Number),y=R.Volume,color=R.Type),size=3,pch=23,fill="black")+labs(fill="1UYL")+labs(color="1UYLL")
Vol_dist+geom_point(aes(x=factor(R.Number),y=R.Volume,bg=R.Type),size=3,pch=23,col="black")+labs(fill="1UYL")+labs(bg="1UYLL")
Vol_dist+geom_point(aes(x=factor(R.Number),y=R.Volume,bg=R.Type),size=3,pch=23,col="black")+labs(color="1UYL")+labs(bg="1UYLL")
Vol_dist+geom_point(aes(x=factor(R.Number),y=R.Volume,bg=R.Type),size=3,pch=23,col="black")+labs(fill="1UYL")+labs(color="1UYLL")
Vol_dist+geom_point(aes(x=factor(R.Number),y=R.Volume,bg=R.Type),size=3,pch=23)+labs(fill="1UYL")+labs(color="1UYLL")
Vol_dist+geom_point(aes(x=factor(R.Number),y=R.Volume,bg=R.Type),size=3,pch=23)+labs(fill="1UYL")+labs(bg="1UYLL")
Vol_dist+geom_point(aes(x=factor(R.Number),y=R.Volume,bg=R.Type),size=3,pch=23)+labs(color="1UYL")+labs(bg="1UYLL")
Vol_dist+geom_point(aes(x=factor(R.Number),y=R.Volume,color=R.Type),size=3,pch=23)+labs(color="1UYL")+labs(fill="1UYLL")
Vol_dist+geom_point(aes(x=factor(R.Number),y=R.Volume,color=R.Type),size=3,pch=23,bg="grey")+labs(color="1UYL")+labs(fill="1UYLL")
Vol_dist+geom_point(aes(x=factor(R.Number),y=R.Volume,color=R.Type),size=3,pch=23,bg="white")+labs(color="1UYL")+labs(fill="1UYLL")
Vol_dist+geom_point(aes(x=factor(R.Number),y=R.Volume,color=R.Type),size=3,pch=23,bg="white")+labs(color="1UYL ertmės tūris ir tipas")+labs(fill="Ertėms tipas")
Vol_dist+geom_point(aes(x=factor(R.Number),y=R.Volume,color=R.Type),size=3,pch=23,bg="white")+labs(color="1UYL ertmės \ntūris ir tipas")+labs(fill="Ertėms tipas")
Vol_dist+geom_point(aes(x=factor(R.Number),y=R.Volume,color=R.Type),size=3,pch=23,bg="white")+labs(color="1UYL ertmės \ntūris ir tipas")+labs(fill="Ertmės tipas")
Vol_dist<-Vol_dist+geom_point(aes(x=factor(R.Number),y=R.Volume,color=R.Type),size=3,pch=23,bg="white")+labs(color="1UYL ertmės \ntūris ir tipas")+labs(fill="Ertmės tipas")
V_dist<-Vol_dist+labs(fill="Ertmės tipas baltymo struktūroje")+xlab("Ertmės indeksas beligandėje 1UYL struktūroje")+ylab(expression(paste("Tūris, ",ring(A)^3)))+ggt
V_dist<-Vol_dist+labs(fill="Ertmės tipas baltymo struktūroje")+xlab("Ertmės indeksas beligandėje 1UYL struktūroje")+ylab(expression(paste("Tūris, ",ring(A)^3)))+ggtitle("Skirtingų tipų ertmių tūrių variacija be ligandų")
V_dist
V_dist+labs(color="Ertmės tipas")
V_dist+labs(fill="Ertmės tipas")
V_dist<-V_dist+labs(fill="Ertmės tipas")
V_dist
V_dist+stat_summary(fun.data=n_fun,geom="text",colour="blue",size=4)
V_dist+stat_summary(fun.data=n_fun,geom="text",size=4,aes(color=comparison$P.Type))
V_dist+stat_summary(fun.data=n_fun,geom="text",size=4,aes(color=comparison$P.Type,fun.y=c(-5,10)))
V_dist+stat_summary(fun.data=n_fun,geom="text",size=4,aes(color=comparison$P.Type))
n_fun<-n_fun<-function(x){return(data.frame(y=median(x),label=length(x)))}
n_fun<-function(x){return(data.frame(y=median(x),label=length(x)))}
V_dist+stat_summary(fun.data=n_fun,geom="text",size=4,aes(color=comparison$P.Type))
V_dist
cairo_pdf('R_CavCleft_P_sulyginimas.pdf',family="DejaVu Sans", width=10)
V_dist
dev.of()
dev.off()
comparison<-data[which(data$Type=="R_PL" & data$R.Number!="NA" $ data$PL.Number!="NA"),]
comparison<-data[which(data$Type=="R_PL" & data$R.Number!="NA" & data$PL.Number!="NA"),]
Vol_dist<-ggplot()+geom_boxpot(aes(x=factor(R.Number),y=PL.Volume))

Vol_dist<-ggplot(comparison, aes(x=comparison$R.Number,y=comparison$PL.Volume,fill=comparison$PL.Type))+geom_boxplot()
Vol_dist
Vol_dist<-ggplot(comparison, aes(x=factor(comparison$R.Number),y=comparison$PL.Volume,fill=comparison$PL.Type))+geom_boxplot()
Vol_dist
V_dist<-Vol_dist+labs(fill="Ertmės tipas baltymo struktūroje")+xlab("Ertmės indeksas beligandėje 1UYL struktūroje")+ylab(expression(paste("Tūris, ",ring(A)^3)))+ggtitle("Skirtingų tipų ertmių tūrių variacija su ligandais")
Vol_dist
Vol_dist<-Vol_dist+labs(fill="Ertmės tipas baltymo struktūroje")+xlab("Ertmės indeksas beligandėje 1UYL struktūroje")+ylab(expression(paste("Tūris, ",ring(A)^3)))+ggtitle("Skirtingų tipų ertmių tūrių variacija su ligandais")
Vol_dist
Vol_dist<-Vol_dist+labs(fill="Ertmės tipas")+xlab("Ertmės indeksas beligandėje 1UYL struktūroje")+ylab(expression(paste("Tūris, ",ring(A)^3)))+ggtitle("Skirtingų tipų ertmių tūrių variacija su ligandais")
Vol_dist
Vol_dist+geom_point(aes(x=factor(R.Number),y=R.Volume,color=R.Type),size=3,pch=23,bg="white")+labs(color="1UYL ertmės \ntūris ir tipas")+labs(fill="Ertmės tipas")
Vol_dist
Vol_dist+geom_point(aes(x=factor(R.Number),y=R.Volume,color=R.Type),size=3,pch=23,bg="white")+labs(color="1UYL ertmės \ntūris ir tipas")+labs(fill="Ertmės tipas")
Vol_dist<-Vol_dist+geom_point(aes(x=factor(R.Number),y=R.Volume,color=R.Type),size=3,pch=23,bg="white")+labs(color="1UYL ertmės \ntūris ir tipas")+labs(fill="Ertmės tipas")
cairo_pdf('R_CavCleft_PL_sulyginimas.pdf',family="DejaVu Sans", width=10)
Vol_dist<-Vol_dist+geom_point(aes(x=factor(R.Number),y=R.Volume,color=R.Type),size=3,pch=23,bg="white")+labs(color="1UYL ertmės \ntūris ir tipas")+labs(fill="Ertmės tipas")
Vol_dist
dev.off()
savehistory(".Rhistory3")

quit()
