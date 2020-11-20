### Function to assess niche innerness and outerness based on  distance in two-dimensions ecological space
### Written by Blaise Petitpierre, 15.11.2015 / update 22.03.2018 by Olivier Broennimann
### Requires the packages rgeos and sp 
### foc.pop = SpatialPoints object of the population of interest in the environmental space
### niche = SpatialPolygon object used to delimit the margin of the niche in the environmental space
###
NMI<-function(foc.pop,niche){

  require(sp)
  require (rgeos)
  
  margin<-as(niche, 'SpatialLines') #niche margin
  niche.regpoints <- spsample(niche,n=100000,"regular") #regular points within the niche
  max.inner<-max(gDistance(niche.regpoints,margin, byid=T)) #max distance to margin
  
  inout<-over(foc.pop,niche) # belonging to the niche
  NMI.abs<-gDistance(spgeom1=foc.pop, spgeom2=margin, byid=T)
  NMI.abs[which(is.na(inout))]<- -NMI.abs[which(is.na(inout))]
  NMI<-NMI.abs/max.inner
  
  return (list(NMI = c(NMI),NMI.abs=c(NMI.abs)))
}
