#load pakages
library(raster)
library(marmap)
library(fossil)

#load funcion
refac_depth<-function(sites,bat,km_deg){
  sites_hd<-sites
  sites_hd<-get.depth(bat,sites_hd,locator = F)
  sites_matrix<-list()
  new_depths<-list()
  for (i in 1:dim(sites_hd)[1]){
    sites_matrix[[i]]<-matrix(nrow = 4,ncol=2,data=c(
      sites_hd[i,1] + km_deg,sites_hd[i,1] - km_deg,
      sites_hd[i,1] + km_deg,sites_hd[i,1] - km_deg,
      sites_hd[i,2] - km_deg,sites_hd[i,2] + km_deg,
      sites_hd[i,2] + km_deg,sites_hd[i,2] - km_deg
    ))
    new_depths[[i]]<-get.depth(bat,sites_matrix[[i]],locator=F)
  }
  sel_min<-function(x){x[which.min(x$depth),]}
  x<-do.call("rbind",lapply(new_depths, sel_min))
  # sites_hd[x$depth < 0,] <- x[x$depth < 0,]
  tt<- ifelse( sites_hd$depth > 0 & x$depth < 0, TRUE,FALSE)
  #sites[,c(1,2)] <- ifelse(sites$depth > 0,ifelse(x$depth < 0, x[,c(1,2)],sites[,c(1,2)]),sites[,c(1,2)])
  sites_hd[tt == TRUE,] <- x[tt==TRUE,]
  sites_hd[,c(1,2)]
}

#start the analysis: now running code to create dis matrix of points of map
bat <- getNOAA.bathy(155, 179, -55, -33, res = 1, keep=TRUE) ## Bring in a bathymetric map to allow us to do
 # can plot(bat) to see the bathymetric matrix produce from above
km_deg<-.008983 # This is the relative conversion between 1 degree longitude and km
sites <- read.csv("RiverMouth_10sites_forEachPopulation_20200924.csv")
sites$Loc <- 1:10
sites<-apply(X = sites,FUN=as.numeric,MARGIN = 2) # Make sure the format is correct
sites[,2]<-abs(sites[,2]) # Make sure longitudes are absolute - this is weird, but some things have issues with negative longitudes, so we can just make them non-negative even though its not really accurate. Shouldn't make a huge issue - just double check that your distances seems reasonable when they get calculated.
sites<-as.data.frame(sites)
colnames(sites)<-c("lat","long")

sites1 <- refac_depth(sites[,2:1],bat,km_deg) # Re-factor to within 1km if possible
sites2 <- refac_depth(sites1,bat,2*km_deg) # Re-factor remaining depths to within 2km
sites3 <- refac_depth(sites[,2:1],bat,3*km_deg) # Re-factor remaining depths to within 3km
sites3<-get.depth(bat,sites3,locator=F)

sites<-apply(X = sites3,FUN=as.numeric,MARGIN = 2) # Make sure the format is correct
sites <- as.data.frame(sites)

trans2 <- trans.mat(bat, min.depth=0, max.depth =-200) # Convert bathymetric object to transmition matrix with max depth of 5000m

dist2 <- lc.dist(trans2, sites[,1:2], res = "dist") # calculated least cost distance matrix

dist_mat<-as.matrix(dist2)

# export this matrix to your working directory
write.csv(dist_mat, file = "Coordinate_Euclidean_distance_matrix_10sites_WilliamCodeForLeastCostOceanGoing_20200927.csv")
