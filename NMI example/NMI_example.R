### Empty environment
rm(list=ls());gc()

##################################################
##################################################
##################################################

### Load libraries

require(raster)
require(ade4)
require(rgeos)
require(ecospat)
source("NMI_function.R")

library(ecospat)

### Define workdir
PATH="C:/users/obroenni/Google Drive/ecospat/Mammal NCN-matching/Gitub_files/NMI example"
setwd(PATH)

### import data
# native niche
sp.shp<-shapefile("Alces_alces_iucn_redlist_nov2013.shp")
sp.shp<-gUnaryUnion(sp.shp)

# introductions
intros<-read.csv("Alces_alces_introductions.csv")

# climatic data
clim<-getData('worldclim', var='bio', res=10)[[c(2,4,10,11,16:18)]] #select bio2,4,10:11,16:18
clim.df<-na.exclude(getValues(clim))

### pca
pca<- dudi.pca(clim.df,scannf = FALSE, nf = 2)
bkg.scores<-pca$li
sp.scores<-suprow(pca,extract(clim,sp.shp))$li
intros.scores<-suprow(pca,extract(clim,intros[,2:3]))$li

### background and niche margins
z<-ecospat.grid.clim.dyn(bkg.scores,bkg.scores,sp.scores,R=100)
bkg<-z$Z>0
bkg[bkg==0]<-NA
bkg.pol<-rasterToPolygons(bkg, dissolve=T)

sp<-z$z.uncor>0
sp[sp==0]<-NA
sp.pol<-rasterToPolygons(sp, dissolve=T)

### NMI
bkg.pts<-data.frame(coordinates(bkg))
coordinates(bkg.pts) <- cbind(bkg.pts$x , bkg.pts$y)
bkg.NMI<-NMI(foc.pop = bkg.pts,niche=sp.pol)

intro.pts<-intros.scores
coordinates(intro.pts) <- cbind(intros.scores$Axis1 , intros.scores$Axis2)
intro.NMI<-NMI(foc.pop = intro.pts,niche=sp.pol)

### plot

# plot axes
plot(bkg.scores,type="n")

#plot NMI
in.bkg<-!is.na(values(bkg)) #pixels in bkg
pos<-bkg.NMI$NMI>=0 #pixels with positive values

  #plot innerness
  inner<-bkg  
  inner[in.bkg&pos]<-bkg.NMI$NMI[in.bkg&pos] 
  inner[inner==1]<-NA
  Pal.inner <- colorRampPalette(c('white','blue'))
  plot(inner,col=Pal.inner(100),add=T, legend=FALSE)

  #plot outerness
  outer<-bkg  
  outer[in.bkg&!pos]<-bkg.NMI$NMI[in.bkg&!pos]
  outer[outer==1]<-NA
  Pal.outer <- colorRampPalette(c('darkolivegreen4','yellow','white'))
  plot(outer,col=Pal.outer(100),add=T, legend=FALSE)

# accessible area
plot(bkg.pol,add=T,lty=2) 

# native climatic niche margin
plot(sp.pol,add=T,lwd=2,border="darkgreen") 

#plot intros
col<-intros$intro.success
col[col==0]<-"red" #introduction failures
col[col==1]<-"green" #introduction successes
points(intro.pts,col=col,pch=19)
col[col=="red"]<-"darkred" 
col[col=="green"]<-"darkgreen" 
points(intro.pts,col=col,pch=21)

#boxplot
y<-intro.NMI$NMI
x<-intros$intro.success
b<-boxplot(y~x,ylim=c(-0.5,1), outline=F,horizontal = TRUE, ylab= "establishment sucess", 
           xlab="NMI",border=c("red","darkgreen"),frame=FALSE)



