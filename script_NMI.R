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
require(mcmcplots)
require(reshape2)

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

saveRDS(allData, file=paste0("data/allData_nmi_", grain,"min_",envelope,level,".rData"))

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
b<-boxplot(y~x,ylim=c(-11,1), outline=T,horizontal = TRUE, ylab= "establishment sucess", 
           xlab="NMI",border=c("red","darkgreen"),frame=FALSE)

# wilcoxon rank sum test U
wilcox.test(y~x)

# number and percent of introductions with NMI>0 
sum(y>0)
sum(y>0)/length(y)

# number and percent of successful introductions with NMI>0 
sum(y[x==1]>0)
sum(y[x==1]>0)/length(y[x==1])

# number and percent of successful introductions with NMI>0 
sum(y[x==0]>0)
sum(y[x==0]>0)/length(y[x==0])

################################################################################################
#  Bayesian Mixed Model               (15.11152 mins with CPU 2.6 GHz / RAM 8 Go)              #
################################################################################################

t1<-Sys.time()

### Define global settings
CV <- FALSE # no cross validation -> model run with the full dataset
allData<-readRDS(paste0("data/allData_nmi_", grain,"min_",envelope,level,".rData"))

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


#------------ Evaluate models fit and explanatory power

### Posterior predictive checks (Bayesian p-value)
comb <- combine.mcmc(out)
mean(comb[,"bpvalue"]) # a value close to 0 or 1 indicates a lack of fit

### Get observed data
Y.obs <- data.jags$Y
Predictors <- with(data.jags, cbind(suitability, Island, Xtrait))

### Get coefficients
all.coef <- comb[, c("mean.beta", "beta.island", "beta.trait[1]", "beta.trait[2]", "beta.trait[3]",
                     "beta.trait[4]", "beta.trait[5]", "beta.trait[6]", "beta.trait[7]", "beta.trait[8]")]
intercept <- comb[, "mean.alpha"]

### Predictions related to fixed effects
l1 <- list()
for(i in 1:ncol(all.coef)) l1[[i]] <- outer(all.coef[,i], Predictors[,i], "*")
tmp.pred <- Reduce(`+`, lapply(l1, function(x) replace(x, is.na(x), 0))) # deal with missing trait values
tmp.pred <- tmp.pred*NA^!Reduce(`+`, lapply(l1, function(x) !is.na(x)))
pred.psi <- intercept + tmp.pred

### Add random effects to predictions
for(i in 1:ncol(pred.psi)){
  pred.psi[,i] <- pred.psi[,i] + comb[, paste0("alpha1[", ID1[i], "]")] + comb[, paste0("alpha2[", ID2[i], "]")] +
    comb[, paste0("alpha3[", ID3[i], "]")]
}

### Inverse logit
pred.psi <- plogis(pred.psi)

### Median of the posterior
pred.psi <- apply(pred.psi, 2, median)

### Compute max sensitivity specificity and TSS
Tresh <- seq((min(pred.psi)+0.05), (max(pred.psi)-0.05), 0.01)
TSS <- sensi <- speci <- numeric()
for(i in 1:length(Tresh)){
  tmp.pres <- ifelse(pred.psi<Tresh[i], 0, 1)
  tmp.tab <- table(tmp.pres, Y.obs)
  sensi[i] <- tmp.tab[2,2]/(tmp.tab[2,2]+tmp.tab[1,2])
  speci[i] <- tmp.tab[1,1]/(tmp.tab[2,1]+tmp.tab[1,1])
  TSS[i] <- sensi[i] + speci[i]  - 1
}
out.tss <- cbind(Tresh, sensi, speci, TSS)
Max.TSS <- out.tss[which.max(out.tss[,4]),]

#------------ Compute posterior probabilities and 95% Highest posterior density intervals of slope coefficients

medians <- apply(all.coef, 2, median) 
posterior.prob <- ifelse(medians<0,
                         apply(all.coef, 2, function(x)length(which(x<0))/length(x)), 
                         apply(all.coef, 2, function(x)length(which(x>0))/length(x)))
HPD95 <- apply(all.coef, 2, function(x)boa.hpd(x, alpha=0.05)) 

#------------ Draw Figure 3

posterior.plot <- mcmc_areas(all.coef, prob=0.95)+
  geom_vline(xintercept=0, linetype="dashed")+
  ggtitle("Posterior distributions",  "of fixed effects")+
  theme(axis.text.x=element_text(size=14, color="black", face="plain"), 
        axis.text.y=element_text(size=14, color="black", face="plain"), 
        plot.title=element_text(size=24, color="black", face="plain"), 
        plot.subtitle=element_text(size=18, color="black", face="plain"))
print(posterior.plot)

#------------ Draw elements of Figure S5 (these elements were then arranged altogether using illustrator)

par(mar=c(5,6,4,2))

############### Panel A - estimates

### Species-wise estimates

sp.estimates <- comb[,grep("alpha1",colnames(comb))]
caterplot(sp.estimates, reorder=F, quantiles=list(outer=c(0.025,0.975), inner=c(0.025,0.975)),
          lwd=c(1.5,1.5), labels=gsub("_"," ",lev.sp), style="plain")
abline(v=0, lty=2)
title(expression(paste("Species coefficients(", alpha[s(fam)], ")")))

### Family-wise estimates

family.estimates <- comb[,grep("alpha3",colnames(comb))]
caterplot(family.estimates, reorder=F, quantiles=list(outer=c(0.025,0.975), inner=c(0.025,0.975)),
          lwd=c(1.5,1.5), labels=gsub("_"," ",lev.fam), style="plain")
abline(v=0, lty=2)
title(expression(paste("Family coefficients(", alpha[fam], ")")))

### Region-wise estimates

region.estimates <- comb[,grep("alpha2",colnames(comb))]
caterplot(region.estimates, reorder=F, quantiles=list(outer=c(0.025,0.975), inner=c(0.025,0.975)),
          lwd=c(1.5,1.5), labels=gsub("_"," ",lev.realm), style="plain")
abline(v=0, lty=2)
title(expression(paste("region coefficients(", alpha[region], ")")))

############### Panel B - 95% Highest Posterior Density Intervals

names <- c("NMI", "Island", trait.names)
caterplot(all.coef, reorder=F, quantiles=list(outer=c(0.025,0.975), inner=c(0.025,0.975), 
                                              lwd=c(1,1)), labels=names, style="plain")
abline(v=0,lty=2)
title("95% Highest Posterior Density Intervals")

############### Panel C - posterior distributions

### Posterior distribution of NMI effect
post.coef <- comb[,"mean.beta"]
densplot(post.coef, main="Posterior distribution of NMI effect", ylab="Density")
abline(v=0, lty=2)
prob <- length(which(post.coef > 0)) / length(post.coef) * 100
text(x=0.6, y=3.0, paste("Posterior probability \n for positive effect = ", round(prob,2), "%"))

### Posterior distribution of success probability
post.inter <- comb[,"mean.alpha"]
densplot(plogis(post.inter), main="Posterior distribution of success probability", ylab="Density")

############### Panel D - posterior relationships with 95% highest posterior density intervals 

### Get original data to plot observations
variables.unscaled <- readRDS(paste0("data/allData_nmi_", grain,"min_",envelope,level,".rData"))
variables.unscaled[,c(9:13)] <- apply(variables.unscaled[,c(9:13)], 2, function(x)log(x+1)) #nat_range"   "intro_date"  "intro_effo"  "weaning_ag" "litter_siz"
variables.unscaled <- variables.unscaled[,c(4,17,9:16)] 

### Scale variables for predictions
variables.scaled <- variables.unscaled
variables.scaled[,2:ncol(variables.scaled)] <- apply(variables.scaled[,2:ncol(variables.scaled)], 2, scale)

### Get coefficients
all.coef <- comb[,grep("beta.trait", colnames(comb))]
all.coef <- cbind(comb[,"mean.beta"], all.coef)

### Get slope coefficients and main intercept
all.coef <- comb[, c("mean.beta", "beta.trait[1]", "beta.trait[2]", "beta.trait[3]", 
                     "beta.trait[4]", "beta.trait[5]", "beta.trait[6]", "beta.trait[7]", "beta.trait[8]")] 

inter <- comb[,"mean.alpha"]

### Perform predictions
df.pred <- NULL
new.names <- c("nmi",trait.names)

for(i in 2:ncol(variables.scaled)){
  
  trait.seq1 <- seq(min(variables.unscaled[,i],na.rm=T), max(variables.unscaled[,i], na.rm=T), length.out=1000)
  trait.seq2 <- seq(min(variables.scaled[,i],na.rm=T), max(variables.scaled[,i], na.rm=T), length.out=1000)
  
  coef <- all.coef[,(i-1)]
  
  post.trait <- NULL
  post.trait <- t(plogis(outer(coef, trait.seq2, "*") + inter))
  
  med.trait <- apply(post.trait, 1, median)
  hpd.trait <- apply(post.trait, 1, quantile, probs=c(.025,.975))
  min.trait <- hpd.trait[1,]
  max.trait <- hpd.trait[2,]
  tmp <- data.frame(cbind(med=med.trait, min=min.trait, max=max.trait, val.trait=trait.seq1))
  tmp1 <- cbind(tmp, Trait=new.names[(i-1)])
  
  if(i==2) df.pred <- tmp1  else df.pred <- rbind(df.pred, tmp1)
  
  print(i)
}

### Re-arrange observed data to feed ggplot
DF1 <- variables.unscaled[,-1]
colnames(DF1) <- new.names
DF1 <- melt(DF1)
DF1 <- cbind(DF1, rep(variables.unscaled$success, 9))
colnames(DF1)[c(1,3)]=c("Trait","success")

Predictions <- ggplot(df.pred)+
  geom_ribbon(aes(ymin=min, ymax=max, x=val.trait), alpha=0.8, linetype=0, fill="lightgrey")+
  geom_line(aes(y=med, x=val.trait))+
  geom_point(data=DF1, aes(y=success, x=value))+
  theme_bw()+
  facet_wrap(~Trait, ncol=3, scales="free_x")+
  theme(strip.background=element_blank(), legend.position="none", strip.text.y=element_blank(), 
        strip.text.x=element_text(size=16), axis.text.x=element_text(size=14, color="black", face="plain"), 
        axis.text.y=element_text(size=14, color="black", face="plain"), 
        plot.title=element_text(size=24, color="black", face="plain"), 
        plot.subtitle=element_text(size=18, color="black", face="plain"), 
        axis.title.x=element_text(size=14), axis.title.y=element_text(size=18), panel.spacing=unit(1.2, "lines"))+
  labs(x="", y="Success probability")+ ggtitle("Posterior relationships", "with 95% highest posterior density intervals")
print(Predictions)

t2<-Sys.time()
t2-t1 # 15.11152 mins with CPU 2.6 GHz / RAM 8 Go 
