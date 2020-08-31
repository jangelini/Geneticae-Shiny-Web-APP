library(ggplot2)
library(ggthemes)

setwd("G:\\Especializacion en bioinformatica\\Para proyecto final\\Escrito\\ezthesis-master\\Graficos")

# width=590
# height= 474
windowsFonts(Times = windowsFont("Times New Roman"))

# AMMI Biplot
datos<-read.csv("datos_biplotGE.csv", sep=";")

AMMI1<-ggplot(data=datos,aes(x=PC1,y=PC2,group="type"))+
  theme_few() +
  theme(plot.caption = element_text(size=12,hjust=0)) +
  scale_size_manual(values=c(4,4))+
  scale_x_continuous("Componente 1", limits=c(-1, 1), breaks = seq(from = -1, to = 1, by = 0.25))+
  scale_y_continuous("Componente 2", limits=c(-1, 1), breaks = seq(-1, 1, 0.25))+
  geom_hline(yintercept=0)+geom_vline(xintercept=0)

AMMI1+geom_segment(xend=0,yend=0,col=alpha("black",1),linetype = 2, size=0.70,data=subset(datos,type=="env"))+
    geom_text(aes(label=name,size=type),show.legend = FALSE, check_overlap = TRUE, size=3.5,
              hjust=0.25, vjust=-0.25)+
    theme(text=element_text(family="Times", size=12,colour = "black"))


# GGE Biplot
datosGGE<-read.csv("datos_biplotGGE.csv", sep=";")
GGE1<-ggplot(data=datosGGE,aes(x=PC1,y=PC2,group="type"))+
  theme_few() +
  theme(plot.caption = element_text(size=12,hjust=0)) +
  scale_size_manual(values=c(4,4))+
  scale_x_continuous("Componente 1", limits=c(-2, 2), breaks = seq(from = -2, to =21, by = 0.5))+
  scale_y_continuous("Componente 2", limits=c(-2, 2), breaks = seq(-2, 2, 0.5))+
  geom_hline(yintercept=0,linetype = 3)+geom_vline(xintercept=0,linetype = 3)+
  theme(text=element_text(family="Times", size=12,colour = "black"))



# Ranking de genotipos para un ambiente
coordenviroment<-datosGGE[datosGGE$type=="env",]
coordgenotype<-datosGGE[datosGGE$type=="gen",]
coordgenotype<-coordgenotype[,-1]

x1<-NULL
for (i in 1:nrow(coordgenotype))
{
  x <- solve(matrix(c(-coordenviroment[coordenviroment$name=="D",4], coordenviroment[coordenviroment$name=="D",3], coordenviroment[coordenviroment$name=="D",3], coordenviroment[coordenviroment$name=="D",4]), nrow = 2),
             matrix(c(0, coordenviroment[coordenviroment$name=="D",3] * coordgenotype[i, 2] +  coordenviroment[coordenviroment$name=="D",4] * coordgenotype[i, 3]), ncol = 1))

  x1<-rbind(x1,t(x))
}
x1
datosGGE$x1_x<-NA
datosGGE$x1_x[datosGGE$type=="gen"]<-x1[,1]

datosGGE$x1_y<-NA
datosGGE$x1_y[datosGGE$type=="gen"]<-x1[,2]



GGE1+
  geom_abline(slope=-0.75/1.20,intercept=0,size=0.7)+
  geom_abline(slope=1.20/0.75,intercept=0)+
  geom_segment(xend=0,yend=0,data=subset(datosGGE,type=="env"&name=="D"),
               arrow =arrow(ends ="first",length=unit(0.2,"inches") ),size=1) +
  geom_segment(aes(xend=x1_x,yend=x1_y),col="black",data=subset(datosGGE,type=="gen"),linetype=2)+
  geom_text(aes(label=name),show.legend = FALSE,data=subset(datosGGE,type=="gen"),col="black",check_overlap = TRUE, size=3.5,
            hjust=1.3, vjust=0)+
  geom_text(aes(label=name),show.legend = FALSE,data=subset(datosGGE,type=="env"&!name%in%c("D")),col="black",check_overlap = TRUE, size=3.5,
            hjust=-0.5, vjust=0) +
  geom_text(aes(label=name),show.legend = FALSE,data=subset(datosGGE,type=="env"&name%in%c("D")), hjust=-0.5, vjust=0,
            col="black",size=4,fontface="bold")+
  theme(text=element_text(family="Times", size=12,colour = "black"))





# Ranking de genotipos para un ambiente
datosGGE<-read.csv("datos_biplotGGE.csv", sep=";")
coordenviroment<-datosGGE[datosGGE$type=="gen",]
coordgenotype<-datosGGE[datosGGE$type=="env",]
coordgenotype<-coordgenotype[,-1]

x1<-NULL
for (i in 1:nrow(coordgenotype))
{
  x <- solve(matrix(c(-coordenviroment[coordenviroment$name=="2",4], coordenviroment[coordenviroment$name=="2",3], coordenviroment[coordenviroment$name=="2",3], coordenviroment[coordenviroment$name=="2",4]), nrow = 2),
             matrix(c(0, coordenviroment[coordenviroment$name=="2",3] * coordgenotype[i, 2] +  coordenviroment[coordenviroment$name=="2",4] * coordgenotype[i, 3]), ncol = 1))

  x1<-rbind(x1,t(x))
}
x1
datosGGE$x1_x<-NA
datosGGE$x1_x[datosGGE$type=="env"]<-x1[,1]

datosGGE$x1_y<-NA
datosGGE$x1_y[datosGGE$type=="env"]<-x1[,2]



GGE1+
  geom_abline(slope=-1.6/1.5,intercept=0,size=0.7)+
  geom_abline(slope=1.5/1.6,intercept=0)+
  geom_segment(xend=0,yend=0,data=subset(datosGGE,type=="gen"&name=="2"),
               arrow =arrow(ends ="first",length=unit(0.2,"inches") ),size=1) +
  geom_segment(aes(xend=x1_x,yend=x1_y),col="black",data=subset(datosGGE,type=="env"),linetype=2)+
  geom_text(aes(label=name),show.legend = FALSE,data=subset(datosGGE,type=="env"),col="black",check_overlap = TRUE, size=3.5,
            hjust=0, vjust=0)+
  geom_text(aes(label=name),show.legend = FALSE,data=subset(datosGGE,type=="gen"&name%in%c("2")), hjust=-0.5, vjust=0,
          col="black",size=4,fontface="bold")+
  geom_text(aes(label=name),show.legend = FALSE,data=subset(datosGGE,type=="gen"&!name%in%c("2")), hjust=-0.5, vjust=0,
            col="black",size=4)+
  theme(text=element_text(family="Times", size=12,colour = "black"))



# Comparacion de dos genotipos
datosGGE<-read.csv("datos_biplotGGE.csv", sep=";")

GGE1+geom_point(data=subset(datosGGE,type=="gen"&name%in%c("6","8")), size=0.7)+
  geom_segment(x=datosGGE$PC1[datosGGE$name=="6"&datosGGE$type=="gen"],
               xend=datosGGE$PC1[datosGGE$name=="8"&datosGGE$type=="gen"],
               y=datosGGE$PC2[datosGGE$name=="6"&datosGGE$type=="gen"],
               yend=datosGGE$PC2[datosGGE$name=="8"&datosGGE$type=="gen"],col="black",size=0.7)+
  geom_abline(intercept = 0, slope = -(datosGGE$PC1[datosGGE$name=="6"] - datosGGE$PC1[datosGGE$name=="8"])/(datosGGE$PC2[datosGGE$name=="6"] - datosGGE$PC2[datosGGE$name=="8"]), col = "black",size=0.7)+
  geom_text(aes(label=name),show.legend = FALSE,data=subset(datosGGE,type=="env"),col="black",size=4)+
  geom_text(aes(label=name),show.legend = FALSE,data=subset(datosGGE,type=="gen"&!name%in%c("6","8")),col="black",size=4)+
  geom_text(aes(label=name),show.legend = FALSE,data=subset(datosGGE,type=="gen"&name%in%c("6","8")), hjust=-0.5, vjust=0,
            col="black",size=4,fontface="bold")

# Poligono envolvente
datosGGE<-read.csv("datos_biplotGGE.csv", sep=";")
coordgenotype<-datosGGE[datosGGE$type=="gen",-c(1,2)]
coordenviroment<-datosGGE[datosGGE$type=="env",-c(1,2)]
indice = c(grDevices::chull(coordgenotype[, 1], coordgenotype[,2]))
www<-data.frame(coordgenotype[indice,])


indice<-c(indice,indice[1])
######################################
segs<-NULL


limx<-layer_scales(GGE1)$x$limits
limy<-layer_scales(GGE1)$y$limits

i <- 1
while (is.na(indice[i + 1]) == FALSE)
{
  m<-(coordgenotype[indice[i], 2] - coordgenotype[indice[i + 1], 2])/(coordgenotype[indice[i],1]-coordgenotype[indice[i + 1],1])
  mperp<--1/m
  c2<-coordgenotype[indice[i + 1], 2] - m*coordgenotype[indice[i + 1],1]
  xint<--c2/(m-mperp)
  xint<-ifelse(xint<0,min(coordenviroment[, 1],coordgenotype[, 1]), max(coordenviroment[, 1],coordgenotype[, 1]))
  yint<-mperp*xint

  xprop<-ifelse(xint<0,xint/limx[1],xint/limx[2])
  yprop<-ifelse(yint<0,yint/limy[1],yint/limy[2])

  m3<-which(c(xprop,yprop)==max(c(xprop,yprop)))
  m2<-abs(c(xint,yint)[m3])
  if(m3==1&xint<0)sl1<-(c(xint,yint)/m2)*abs(limx[1])
  if(m3==1&xint>0)sl1<-(c(xint,yint)/m2)*limx[2]
  if(m3==2&yint<0)sl1<-(c(xint,yint)/m2)*abs(limy[1])
  if(m3==2&yint>0)sl1<-(c(xint,yint)/m2)*limy[2]

  segs<-rbind(segs,sl1)
  i <- i + 1
}
rownames(segs)<-NULL
colnames(segs)<-NULL
segs<-data.frame(segs)



colnames(www)<-c("X1","X2")

winners<-datosGGE[datosGGE$type=="gen",][indice[-1],]
others<-datosGGE[!rownames(datosGGE)%in%rownames(winners),]

GGE1+geom_polygon(aes(x=X1,y=X2),data=www,fill=NA,col="black",size=0.7)+
  geom_segment(data=segs,aes(x=X1,y=X2),xend=0,yend=0,linetype=2,col="black",size=0.7) +
  geom_segment(xend=0,yend=0,col="black",data=subset(datosGGE,type=="env"))+
  geom_text(aes(label=name,size=type),show.legend = FALSE,data=others)+
  geom_text(aes(label=name,x=PC1,y=PC2),show.legend = FALSE,hjust="outward",vjust="outward",data=winners,col="black",size=4,fontface="bold")


# Mean vs stability
datosGGE<-read.csv("datos_biplotGGE.csv", sep=";")
coordgenotype<-datosGGE[datosGGE$type=="gen",-c(1,2)]
coordenviroment<-datosGGE[datosGGE$type=="env",-c(1,2)]

med1 = mean(coordenviroment[, 1])
med2 = mean(coordenviroment[, 2])


x1<-NULL

for (i in 1:nrow(coordgenotype))
{
  x <- solve(matrix(c(-med2, med1, med1, med2), nrow = 2), matrix(c(0, med2 * coordgenotype[i, 2] + med1 * coordgenotype[i, 1]),ncol = 1))
  x1<-rbind(x1,t(x))

  # segments(coordgenotype[i, dimension1], coordgenotype[i,dimension2], x[1], x[2], lty = "dotted")
}
datosGGE$x1_x<-NA
datosGGE$x1_x[datosGGE$type=="gen"]<-x1[,1]

datosGGE$x1_y<-NA
datosGGE$x1_y[datosGGE$type=="gen"]<-x1[,2]

GGE1+
  geom_abline(intercept=0,slope=med2/med1,col="black",size=0.7)+
  geom_abline(intercept=0,slope=-med1/med2,col="black",size=0.7)+
  geom_segment(aes(xend=x1_x, yend=x1_y),col="black",linetype=2,data=subset(datosGGE,type=="gen"))+
  geom_segment(x=0, y=0,xend=med1,yend=med2,arrow =arrow(length=unit(0.15,"inches") ),size=1,col="black")+
  geom_text(aes(label=name,size=type),show.legend = FALSE)

