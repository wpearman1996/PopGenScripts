library(marmap)
bat <- getNOAA.bathy(169.1, 178.776, -44, -33, res = 1, keep=TRUE) ## Bring in a bathymetric map to allow us to do some least cost path analaysis
Hatfields<-c(-36.569330, 174.708932)
StanmoreBay<-c(-36.597523, 174.738320)
StanmoreBayOld<-c(-36.597523, 174.738321)
BrownsBay<-c(-36.717101, 174.768741)
OpitoBay<-c(-36.689997, 175.816381)
MtMaunganui<-c(-37.613310, 176.160670)
Mahia<-c(-39.074737, 177.933028)
WorserBay<-c(-41.374389, 174.794772)
Kaikoura<-c(-42.425248, 173.717603)
coords<-as.data.frame(rbind(StanmoreBayOld,Kaikoura,Hatfields, BrownsBay,OpitoBay,
                            MtMaunganui,Mahia,StanmoreBay,WorserBay))
colnames(coords)<-c("y","x")
sites<-coords[c(2,1)]
rownames(sites)<-1:9
# Convert to transition object
trans1<-trans.mat(bat)
#Prevent movement in depths > 150m
trans2 <- trans.mat(bat, max.depth =-150)
library(fossil);library(raster)
dist2 <- lc.dist(trans2, sites, res = "dist")
dist_mat<-as.matrix(dist2)
dist_mat[upper.tri(dist_mat,diag = T)] <- NA
nifst<-fst_mat[c(1,3,4,5,6,7,8),c(1,3,4,5,6,7,8)]
ni_dist<-dist_mat[c(1,3,4,5,6,7,8),c(1,3,4,5,6,7,8)]
