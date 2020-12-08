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
require(ks)
require(dplyr)
require(sf)
require(nngeo)

#source("NMI_function.R")


#bootstrap.functions
  
# tricky functions for serialisation
k.val<-function(to.pcki,mati) {
  return(apply(to.pcki,1,function(to.pcki) mati$estimate[to.pcki[1],to.pcki[2]]))
}
kd.bin<-function(x,y){
  return(matrix(sapply(x$estimate,function(x)ifelse(x<y,0,1)),nrow=nrow(x$estimate)))
}
  
#core bootstrap function
kd.bootstrap<-function(nrep=100,sp.scores=sp.scores,R = 100, th.quant = 1){
  
  my.df<-dplyr::as_tibble(sp.scores) # convert to tibble to do faster resampling
  bst<-list()
  bst <- vector(mode = "list", length = nrep) %>% 
    lapply(function(x)dplyr::sample_n(my.df,nrow(my.df),replace=T)) # resampled data
  
  fhat.sp<-kde(sp.scores,gridsize = c(R,R))  # kernel on observation
  k.value<-contourLevels(fhat.sp, cont=c(th.quant*100)) # density threshold to include observations
  
  cquant<-contourLines(fhat.sp$eval.points[[1]],fhat.sp$eval.points[[2]],fhat.sp$estimate,level=k.value) #create a linear enveloppe
  lquant<-list()        #geometry handling
  for(k in 1:length(cquant))lquant[[k]]=Polygons(list(Polygon(cquant[[k]][-1])),ID=k)#geometry handling
  sp.pol=SpatialPolygons(lquant)# plottable envelope

  cpoints<-lapply(bst,kde,gridsize = c(R,R),xmin = sapply(fhat.sp$eval.points,min), # compute kernel density on the resampled data
                  xmax = sapply(fhat.sp$eval.points,max)) 
  k.value<-lapply(cpoints,contourLevels,cont=c(th.quant*100)) ## compute density threshold to include observations for each resampled dataset
  bin.bst<-mapply(kd.bin,x=cpoints,y=k.value,SIMPLIFY = F) ## binarized kernel distribution
  sum.bst<-Reduce("+",bin.bst) ## sum of the binarized kernel distribution

  kdraster<-raster(t(sum.bst[,nrow(sum.bst):1]),xmn=fhat.sp$eval.points[[1]][1], xmx=range(fhat.sp$eval.points[[1]])[2], ## rasterization of sum matrix
                   ymn=fhat.sp$eval.points[[2]][1], ymx=range(fhat.sp$eval.points[[2]])[2])
  kdraster<-kdraster/nrep*100 ## scaling in %
  sp.envelope<-st_as_sf(sp.pol) %>% ## convert polygon envelope into line
    st_cast('LINESTRING')
  return(list(kd.raster =kdraster, sp.envelope = sp.envelope, kd.sp=fhat.sp))
}

### Define workdir
#PATH="C:/users/obroenni/Google Drive/ecospat/Mammal NCN-matching/Gitub_files/NMI example"
PATH="F:/UNIL/NMI-margin_bootstrap/NMI example"
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

system.time(
  my.bootstrap<-kd.bootstrap(nrep=100,sp.scores=sp.scores,R = 100, th.quant = 1)
  ) # 50 seconds
#nrep = number of replication of the bootstrap
#sp.scores = pca scores for the species distribution
#R = resolution of the grid for the kernel analysis
#th.quant = percentage of the distribution included within the kernel enveloppe.1 -> include 100% of the observation

#plot margins uncertainty
plot(my.bootstrap$kd.raster/maxValue(my.bootstrap$kd.raster)*100,axes = F) # sum of the stack of the bootstrapped distribution = confidence interval
points(sp.scores,cex=0.1) # distribution of the species occurrences
plot(my.bootstrap$sp.envelope,add=T) # kernel envelop based on the observed distribution
plot(rasterToContour(my.bootstrap$kd.raster,levels=95),col=2,add=T) # kernel envelop based on the 95% confidence interval of the bootstrapped distribution. You can change the level of confidence interval with the "levels" parameter

