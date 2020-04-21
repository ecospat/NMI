### Empty environment
rm(list=ls());gc()

##################################################
##################################################
##################################################

### Load libraries
require(R2jags)
require(runjags)
require(ggplot2)
require(bayesplot)

### Define workdir
PATH="D:/Post-doc Lausanne/Analyses/innerness mammals/"
setwd(PATH)

### Get the data to run JAGS
Data <- get(load("Gitub_files/data_jags.Rdata"))

#-- Define whether the model should be run with NMI or HS as predictor
Data[["focus.var"]] <- scale(Data$NMI)[,1] 
# Data[["focus.var"]] <- scale(Data$HS)[,1]  

### Define the variables to monitor
my_var <- c("alpha1", "alpha2", "alpha3", "mean.alpha", "mean.beta", "beta.trait", "beta.island",
            "sd.fam", "sd.realm", "sd.sp", "sd.trait", "fit", "fit.new", "bpvalue")

### Define MCMC settings
nchains=3
nburn=5000
nthin=20
store=1000

### Run Jags (in paralell)
jags.model <- jags.parallel(data=Data,
                            parameters.to.save=my_var,
                            model.file="Gitub_files/JagsMod.bug",
                            n.iter=(store*nthin)+nburn, n.burnin=nburn, n.chains=nchains, n.thin=nthin,
                            export_obj_names=c("store","nthin","nburn","nchains","Data"))

### Save JAGS outputs
save(jags.model, file="Gitub_files/jags_outputs.Rdata")

##################################################
##################################################
##################################################

out <- get(load("Gitub_files/jags_outputs.Rdata"))
out <- as.mcmc(out)

### Gelman and Rubin convergence diagnostic
conv <- gelman.diag(out, multivariate=F, autoburnin=F)
max(conv$psrf[,1]) # < 1.1 for all parameters

### Posterior predictive checks 
comb <- combine.mcmc(out)
mean(comb[,"bpvalue"]) # 0.21
plot(as.vector(comb[,"fit"]) ~ as.vector(comb[,"fit.new"]))
abline(0,1)

### Postrior distribution of slope coefficients
DF.slopes <- comb[,c("mean.beta", "beta.trait[1]", "beta.trait[2]", "beta.trait[3]", 
                     "beta.trait[4]", "beta.trait[5]", "beta.trait[6]", "beta.trait[7]",
                     "beta.trait[8]", "beta.island")] 

colnames(DF.slopes) <- c("NMI", "Native range area", "Introduction date", "Propagule pressure", "Weaning age",
                         "Litter size", "Litters per year", "CV body mass", "CV neonate body mass", "Islands")

post.plot <- mcmc_areas(DF.slopes, prob=0.95)+
  geom_vline(xintercept=0, linetype="dashed")+
  ggtitle("Posterior distributions", "of fixed effects")+
  theme(axis.text.x = element_text(size=14, color="black", face="plain"),
        axis.text.y = element_text(size=14, color="black", face="plain"),
        plot.title = element_text(size=24, color="black", face="plain"),
        plot.subtitle = element_text(size=18, color="black", face="plain"))

tiff(paste0("Gitub_files/Posterior distributions slope coefficients.tiff"),width=3000,height=3000,res=300)
print(post.plot)
dev.off()
