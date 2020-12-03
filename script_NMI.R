### Empty environment
rm(list=ls());gc()

################################################################################################
#   libraries and directories                                                                  #
################################################################################################

### Load libraries
require(raster)
require(ade4)
require(rgeos)
require(rgdal)
require(ecospat)
require(ks)
require(R2jags)
require(runjags)
require(ggplot2)
require(bayesplot)
require(boa)
require(loo)

### Define workdir
PATH="C:/Users/Olivier/Google Drive/ecospat/Mammal NCN-matching/Gitub_files"
setwd(PATH)

################################################################################################
#  calculation of NMI           (1.643185 hours mins with CPU 2.6 GHz / RAM 8 Go)              #
################################################################################################

t1<-Sys.time()
source("https://raw.githubusercontent.com/ecospat/NMI/master/functions/NMI_function.R")
source("https://raw.githubusercontent.com/ecospat/NMI/master/functions/mve_function.R")

### choose the parameters of the analysis
grain=10 # resolution of the climatic data. In Broennimann et al. 2020 NCOM: 10,30,60
envelope="kde" # choise of the envelop methode, "kde" or "mve"
level=99 # level of inclusion of rare climatic conditions for kde and mve. In Broennimann et al. 2020 NCOM: 99,95,90

### data preparation for all species

# import shapefile of native distributions
all.shp<-shapefile("data/allsp_simplifpol.shp")

# import introductions
all.intros<-read.delim("data/intros.txt")
sp.list<-unique(all.intros$sp)

# Generate climatic data at the chosen grain
clim<-getData('worldclim', var='bio', res=10)[[c(2,4,10,11,16:18)]] #select bio2,4,10:11,16:18
clim<-aggregate(clim,grain/10)
clim.df<-na.exclude(getValues(clim))

### run PCA on all background points (PCA-env sensu Broennimann et al 2012)
pca<- dudi.pca(clim.df,scannf = FALSE, nf = 2)

### calculate NMI across species
allData<-data.frame()

for (sp in sp.list){
  # select native niche of focal species and extract scores in PCA space
  sp.shp<- all.shp[all.shp$binomial==sub("_"," ",sp),]
  sp.scores<-na.omit(suprow(pca,raster::extract(clim,sp.shp))$li) 
  if(nrow(sp.scores)<10) next
  
  # create scores of introductions in PCA space
  intros<-all.intros[which(all.intros$sp==sub(" ","_",sp)),]
  intros.scores<-suprow(pca,raster::extract(clim,intros[,2:3]))$li
  NArows<-is.na(intros.scores$Axis1)
  intros<-intros[!NArows,]
  intros.scores<-intros.scores[!NArows,]
  
  ### delineate niche margins
  # kde
  if(envelope=="kde"){
    fhat<-kde(sp.scores,compute.cont=TRUE)
    c99<-contourLines(fhat$eval.points[[1]],fhat$eval.points[[2]],fhat$estimate,level=fhat$cont[level])
    l99<-list()        
    for(k in 1:length(c99))l99[[k]]=Polygons(list(Polygon(c99[[k]][-1])),ID=k)
    sp.pol=SpatialPolygons(l99)
    sp.pol<-aggregate(sp.pol)
  }
  # mve
  if(envelope=="mve"){sp.pol <- mve(sp.scores, thresh=level/100, robust=F)}
  
  ### calculate NMI for introductions
  coordinates(intros.scores) <- cbind(intros.scores$Axis1 , intros.scores$Axis2) # DataFrame to SpatialPointsDataFrame
  intro.NMI<-NMI(foc.pop = intros.scores,niche=sp.pol)
  nmi<-intro.NMI$NMI
  allData<-rbind(allData,cbind(intros,nmi))
  cat(sp," ")
}

### Save the data

saveRDS(allData, file=paste0("data/allData_nmi_", grain,"min_",envelope,level))

t2<-Sys.time()
t2-t1 # 1.643185 hours with CPU 2.6 GHz / RAM 8 Go 

################################################################################################
#  Figure 1C -  Illustration of NMI for the Alces alces                                        #
################################################################################################

# select the species
sp<-"Alces alces"

# extract scores of the background in PCA space
bkg.scores<-pca$li

# extract scores of the native niche in PCA space
sp.shp<- all.shp[which(all.shp$binomial==sp),]
sp.scores<-suprow(pca,raster::extract(clim,sp.shp))$li

# extract scores of introductions in PCA space
intros<-all.intros[which(all.intros$sp==sub(" ","_",sp)),c(4:5,14)]
intros.scores<-suprow(pca,extract(clim,intros[,1:2]))$li

intros<-all.intros[which(all.intros$sp==sub(" ","_",sp)),]
intros.scores<-suprow(pca,raster::extract(clim,intros[,2:3]))$li

### delineate background and niche margins

fhat<-kde(bkg.scores,compute.cont=TRUE)
c99<-contourLines(fhat$eval.points[[1]],fhat$eval.points[[2]],fhat$estimate,level=0.0001)
l99<-list()        
for(k in 1:length(c99))l99[[k]]=Polygons(list(Polygon(c99[[k]][-1])),ID=k)
bkg.pol=SpatialPolygons(l99)
bkg.pol<-aggregate(bkg.pol)

fhat<-kde(sp.scores,compute.cont=TRUE)
c99<-contourLines(fhat$eval.points[[1]],fhat$eval.points[[2]],fhat$estimate,level=fhat$cont[length(fhat$cont)])
l99<-list()        
for(k in 1:length(c99))l99[[k]]=Polygons(list(Polygon(c99[[k]][-1])),ID=k)
sp.pol=SpatialPolygons(l99)
sp.pol<-aggregate(sp.pol)

### calculate NMI values

bkg.grid<-raster(nrows=100,ncols=100,extent(bkg.pol))
bkg<- rasterize(bkg.pol,bkg.grid)
bkg.pts<-data.frame(coordinates(bkg))
coordinates(bkg.pts) <- cbind(bkg.pts$x , bkg.pts$y)
bkg.NMI<-NMI(foc.pop = bkg.pts,niche=sp.pol)

coordinates(intros.scores) <- cbind(intros.scores$Axis1 , intros.scores$Axis2)
intro.NMI<-NMI(foc.pop = intros.scores,niche=sp.pol)

### plot Fig 1C
# plot axes
plot(coordinates(bkg.pts),type="n",xlab="PC1",ylab="PC2")

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
col<-intros$success
col[col==0]<-"red" #introduction failures
col[col==1]<-"green" #introduction successes
points(intros.scores,col=col,pch=19)
col[col=="red"]<-"darkred" 
col[col=="green"]<-"darkgreen" 
points(intros.scores,col=col,pch=21)

#boxplot
y<-intro.NMI$NMI
x<-intros$success
b<-boxplot(y~x,ylim=c(-2,1), outline=F,horizontal = TRUE, ylab= "establishment sucess", 
           xlab="NMI",border=c("red","darkgreen"),frame=FALSE,main=sp)

################################################################################################
#  Figure 2B -  boxplots of NMI and establishment outcome                                      #
################################################################################################

y<-allData$nmi
x<-allData$success
b<-boxplot(y~x,ylim=c(-2,1), outline=F,horizontal = TRUE, ylab= "establishment sucess", 
           xlab="NMI",border=c("red","darkgreen"),frame=FALSE)
wilcox.test(y~x)

################################################################################################
#  Bayesian Mixed Model               (15.11152 mins with CPU 2.6 GHz / RAM 8 Go)              #
################################################################################################

t1<-Sys.time()

### Define global settings
CV <- FALSE # no cross validation -> model run with the full dataset
allData<-readRDS(paste0("data/allData_nmi_", grain,"min_",envelope,level))

### Create directories to store JAGS data and JAGS outputs
if(!file.exists("JAGS_outputs")) dir.create("JAGS_outputs")
if(!file.exists("JAGS_data")) dir.create("JAGS_data")

### log transfrom some variables to reduce the skewness of their distribution

allData$Native_range_area <- log(allData$Native_range_area+1)
allData$Introduction_date <- log(allData$Introduction_date+1)
allData$Introduction_effort <- log(allData$Introduction_effort+1)
allData$Weaning_age <- log(allData$Weaning_age+1)
allData$Litter_size<- log(allData$Litter_size+1)

### Scale variables
allData[,c(9:ncol(allData))] <- apply(allData[,c(9:ncol(allData))], 2, scale)

### Define trait names
trait.names <- c("Native range area", "Introduction date", "Introduction effort", "Weaning age",
                 "Litter size", "Litters per year", "CV adult body mass", "CV neonate body mass")

### Define indicator variables for random effects
lev.sp <- levels(factor(allData$sp))
ID1 <- match(allData$sp, lev.sp) # indicator variable indicating to which species the observations belong to
allData$realm[is.na(allData$realm)]<-"NA"
lev.realm <- levels(factor(allData$realm))
ID2 <- match(allData$realm, lev.realm) # indicator variable indicating to which realm the observations belong to
lev.fam <- levels(factor(allData$fam))
ID3 <- match(allData$fam, lev.fam) # indicator variable indicating to which family the observations belong to
uni <- unique(allData[,c("sp", "fam")])
ID4 <- match(uni$fam, lev.fam) # indicator variable indicating to which family each species belongs to (for the nested effect)

### Full data or Sub-sample data for cross-validation
rows <- 1:nrow(allData)
Data <- allData
if(CV){
  samp <- sample(rows, round(0.2*nrow(Data))) # 20% of observations for the testing dataset
  Data[samp, "success"] <- NA # observations replaced by missing values -> treated as parameters by the model
} else {
  samp <- rows 
}
 
### Define indicator variables to compute posterior predictive checks and predicted presence-absence for model evaluation
if(CV) ID.samp <- rows[-samp] else ID.samp <- rows # which observations to use for model fitting
ID.pred <- rows[samp] # which observations to use for model evaluation

### Jags settings
nchains <- 3
nburn <- 5000
nthin <- 20
store <- 1000

###### Run the mixed effect model ######

### List of data to feed JAGS
data.jags <- list(Y=allData$success,
                  suitability=allData$nmi,
                  Nobs=nrow(allData),
                  id1=ID1,
                  Nfam=nlevels(factor(allData$fam)),
                  id2=ID2,
                  Nrealm=nlevels(factor(allData$realm)),
                  id3=ID3,
                  Nsp=nlevels(factor(allData$sp)),
                  id4=ID4,
                  Island=allData$island,
                  Xtrait=allData[,c(9:16)], 
                  Ntrait=ncol(allData[,c(9:16)]),
                  ID.samp=ID.samp, Nsamp=length(ID.samp),
                  ID.pred=ID.pred, Npred=length(ID.pred))

### Save the data
saveRDS(data.jags, file=paste0("JAGS_data/nmi_", grain,"min_",envelope,level,".rData"))

### Variables to monitor
my_var <- c("Y.pred", "alpha1", "alpha2", "alpha3", "mean.alpha", "mean.beta", "beta.trait", "beta.island",
            "sd.fam", "sd.realm", "sd.sp", "sd.trait", "bpvalue", "log.like")

### Run Jags in parallel
jags.model <- jags.parallel(data=data.jags,
                            parameters.to.save=my_var,
                            model.file="functions/Final_model_CV.bug",
                            n.iter=(store*nthin)+nburn, n.burnin=nburn, n.chains=nchains, n.thin=nthin,
                            export_obj_names=c("store", "nthin", "nburn", "nchains", "data.jags"))

### Save model outputs
saveRDS(jags.model, file=paste0("JAGS_outputs/nmi_", grain,"min_",envelope,level,".rData"))

###### Analyse model outputs ######

### Choose the model for which metrics and graphics have to be computed

jags.model <- readRDS(paste0("JAGS_outputs/nmi_", grain,"min_",envelope,level,".rData"))
data.jags <- readRDS(paste0("JAGS_data/nmi_", grain,"min_",envelope,level,".rData"))

#------------ Check model convergence

out <- as.mcmc(jags.model)
hyper.par <- which(colnames(out[[1]]) %in% c("mean.alpha", "mean.beta", "beta.trait", "beta.island", "sd.fam", "sd.realm", "sd.sp"))
sp.par <- grep("alpha1", colnames(out[[1]]))
realm.par <- grep("alpha2", colnames(out[[1]]))
fam.par <- grep("alpha3", colnames(out[[1]]))
diagnostic <- gelman.diag(out[,c(hyper.par, sp.par, realm.par, fam.par)], autoburnin=F, multivariate=F)
any(diagnostic$psrf[,1]>1.1) # All parameters must have values below 1.1

#------------ Evaluate models fit

### Posterior predictive checks (Bayesian p-value) 
comb <- combine.mcmc(out)
mean(comb[,"bpvalue"]) # a value close to 0 or 1 indicates a lack of fit

### compute sensitivity, specificity and TSS
Y.obs <- na.omit(data.jags$Y) # Observed PA data
Y.pred <- comb[,paste0("Y.pred[", 1:length(ID.pred), "]")] # Predicted PA data
Sensitivity <- Specificity <- TSS <- numeric()

# Compute metrics for all samples of the posterior predictive distribution
for(i in 1:nrow(Y.pred)){ 
  tmp <- table(Y.pred[i,], Y.obs)
  if(nrow(tmp)==1){
    if(rownames(tmp)=="1") tmp <- rbind(c(0,0),tmp) else tmp <- rbind(tmp,c(0,0)) 
  }
  if(ncol(tmp)==1){
    if(colnames(tmp)=="1") tmp <- cbind(c(0,0),tmp) else tmp <- cbind(tmp,c(0,0)) 
  }
  Sensitivity[i] <- tmp[2,2]/(tmp[2,2]+tmp[1,2])
  Specificity[i] <- tmp[1,1]/(tmp[2,1]+tmp[1,1])
  TSS[i] <- Sensitivity[i] + Specificity[i] - 1
}
TSS <- median(TSS)
Sensitivity <- median(Sensitivity)
Specificity <- median(Specificity)

#------------ Compute marginal and conditional R-square

### Get data for predictors
Predictors <- with(data.jags, cbind(suitability, Island, Xtrait))

### Get slope coefficients
all.coef <- comb[, c("mean.beta", "beta.island", "beta.trait[1]", "beta.trait[2]", "beta.trait[3]", 
                     "beta.trait[4]", "beta.trait[5]", "beta.trait[6]", "beta.trait[7]", "beta.trait[8]")] 
colnames(all.coef) <- c(foc.resp, "Island", trait.names)

### Get the intercept
intercept <- comb[,"mean.alpha"]

### Predict occurrence probabilities
l1 <- list()
for(i in 1:ncol(all.coef)) l1[[i]] <- outer(all.coef[,i], Predictors[,i], "*")
tmp.pred <- Reduce(`+`, lapply(l1, function(x) replace(x, is.na(x), 0))) # deal with missing trait values
tmp.pred <- tmp.pred*NA^!Reduce(`+`, lapply(l1, function(x) !is.na(x)))
pred.psi <- plogis(intercept + tmp.pred)

### Compute marginal and conditional R²
err <- sweep(pred.psi, 2, data.jags$Y, FUN="-")
varFixed <- apply(pred.psi, 1, var)
varResidual <- apply(err, 1, var) # get the residual variance 
varRandom <- plogis(comb[,"sd.sp"]^2 + comb[,"sd.realm"]^2 + comb[,"sd.fam"]^2)   # get the variance of random effects

marginalR2 <- varFixed/(varFixed + varRandom + varResidual)
conditionalR2 <- (varRandom + varFixed)/(varFixed + varRandom + varResidual)

marginalR2 <- median(marginalR2)
conditionalR2 <- median(conditionalR2)

#------------ Compute posterior probabilities and 95% Highest posterior density intervals of slope coefficients

medians <- apply(all.coef, 2, median) 
posterior.prob <- ifelse(medians<0,
                         apply(all.coef, 2, function(x)length(which(x<0))/length(x)), 
                         apply(all.coef, 2, function(x)length(which(x>0))/length(x)))
HPD95 <- apply(all.coef, 2, function(x)boa.hpd(x, alpha=0.05)) 

#------------ Draw Figure 3 in the main text

posterior.plot <- mcmc_areas(all.coef, prob=0.95)+
  geom_vline(xintercept=0, linetype="dashed")+
  ggtitle("Posterior distributions",  "of fixed effects")+
  theme(axis.text.x=element_text(size=14, color="black", face="plain"), 
        axis.text.y=element_text(size=14, color="black", face="plain"), 
        plot.title=element_text(size=24, color="black", face="plain"), 
        plot.subtitle=element_text(size=18, color="black", face="plain"))
print(posterior.plot)

t2<-Sys.time()
t2-t1 # 15.11152 mins with CPU 2.6 GHz / RAM 8 Go 