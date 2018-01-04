#Plot the climatology of the hail prediction analysis

rm(list=ls())
library(fields)
library(maps)
library(mapdata)
library(oz)
library(RColorBrewer)
library(RcppRoll) #For rolling sums
library(ncdf4)
library(reshape)
library(ggplot2)
library(scales)

#Read the RDS files output from cluster_hail_prediction.R
#The cluster prediction file
#cls.pred <- readRDS("pred_test_all_20clusters.rds")

#cls.pred <- readRDS("pred_test_2008_6clusters.rds")
cls.pred <- readRDS("pred_test_2008_10clusters.rds")
#cls.pred <- readRDS("pred_test_2008_15clusters.rds")
#cls.pred <- readRDS("pred_test_2008_20clusters.rds")
#cls.pred <- readRDS("pred_test_2008_64clusters.rds")
#cls.pred <- readRDS("pred_test_2008_100clusters.rds")
#cls.pred <- readRDS("pred_test_2008_250clusters.rds")
#The hail and lightning probability files
#hlprob <- readRDS("prop_hail_lightning_6clusters.rds")
hlprob <- readRDS("prop_hail_lightning_10clusters.rds")
#hlprob <- readRDS("prop_hail_lightning_15clusters.rds")
#hlprob <- readRDS("prop_hail_lightning_20clusters.rds")
#hlprob <- readRDS("prop_hail_lightning_64clusters.rds")
#hlprob <- readRDS("prop_hail_lightning_100clusters.rds")
#hlprob <- readRDS("prop_hail_lightning_250clusters.rds")

#Read the ERA and the lightning data
#fnm_era <- '/home/jpeter/justinp/rfiles/suncorp/data/eraint/sharpy_profiles/00_18/2008_2014/era_profiles_2008_2014_00_18.nc'

#fnm_gpats <- '/home/jpeter/justinp/rfiles/suncorp/data/gpats/ts_lightning_grid/2008_2015/combined/lightning_eragrid_eratime_Australia_2008_2015_offset.nc'
#Use lightning time corresponding to available ERA-data
fnm_gpats <- '/home/u1066578/data/ts_lightning_grid/2008_2014/analysis/lightning_eragrid_eratime_Australia_2008_2014_offset.nc'

#fnm_lsmask <- '/home/jpeter/justinp/rfiles/suncorp/data/eraint/netcdf_download/lsmask/era_interim_lsmask.nc'
fnm_lsmask <- '/home/u1066578/data/eraint/lsmask/era_interim_lsmask.nc'

#Open the netcdf files
#nc_era   <- nc_open(fnm_era)
#print(paste("SHARPPy profile file has",nc_era$nvars,"variables"))
nc_gpats <- nc_open(fnm_gpats)
print(paste("GPATS file has",nc_gpats$nvars,"variables"))

nc_lsmask <- nc_open(fnm_lsmask)


#Get the dimensions of the ERA file
#We use this for plotting as we need the lat and lon
#nlat_era    <- nc_era$dim$lat$len
#nlon_era    <- nc_era$dim$lon$len
#ntime_era   <- nc_era$dim$time$len

#lat_era <- nc_era$dim$lat$vals
#lon_era <- nc_era$dim$lon$vals
#time_era <- nc_era$dim$time$vals
#rtime_era <- as.POSIXct(time_era*60*60, origin="1900-01-01", tz="GMT")

v2 <- nc_lsmask$var[["lsm"]]
lsmask <- ncvar_get(nc_lsmask,v2)
lat_lsmask <- nc_lsmask$dim$lat$vals
lon_lsmask <- nc_lsmask$dim$lon$vals

#Get the dimensions of the GPATS file
nlat_gpats  <- nc_gpats$dim$lat$len
nlon_gpats  <- nc_gpats$dim$lon$len
ntime_gpats <- nc_gpats$dim$time$len

lat_gpats  <- nc_gpats$dim$lat$vals
lon_gpats  <- nc_gpats$dim$lon$vals
time_gpats <- nc_gpats$dim$time$vals
rtime_gpats<- as.POSIXct(time_gpats*60*60, origin="1900-01-01", tz="GMT")

gap_storm_date <- "2008-11-16 06:00:00 GMT"
gap_storm_doy  <- as.numeric(strftime(gap_storm_date,format="%j"))

#Get the total_lightning variable
#Get the total lightning variable
v1<-nc_gpats$var[["total_lightning"]]
total_lightning <- ncvar_get(nc_gpats,v1)
#total_lightning[,,] <- t(apply(total_lightning[,,],1,rev))

nc_close(nc_gpats)
#nc_close(nc_era)
nc_close(nc_lsmask)

#Plot some information about the hail probability
hlprob$prop.hail.s <- hlprob$prop.hail*10.
hd.l <- melt(hlprob,id.vars="cl")

ggplot(subset(hd.l,variable %in% c("prop.hail.s","prop.lightning")),
       aes(x=cl,y=value,fill=variable))+
       geom_bar(stat="identity",position="dodge")+
       xlab("Cluster")+
       guides(fill=guide_legend(title=NULL))+
       scale_fill_discrete(breaks=c("prop.lightning","prop.hail.s"),
                           labels=c("Lightning","Hail"))+
       theme(axis.title=element_text(size=12),
             axis.text=element_text(size=12),
             legend.text=element_text(size=12))+
       scale_y_continuous(labels=percent)

phd <- ggplot(hd.l,aes(cl,value,fill=variable))+
       geom_bar(stat="identity",position="dodge")+
       xlab("Cluster")+ylab("%")+
        scale_fill_brewer("",
                          labels=c("Hail","No hail",
                                   "Lightning","No lightning","Lightning strikes"),
                                   palette="Set1")




##Dont include the Indonesian locations
#This excludes one lat point at the top of the domain
min_lat <- min(lat_gpats)
#max_lat <- max(lat_gpats)
max_lat <- -10.5
min_lon <- min(lon_gpats)
max_lon <- max(lon_gpats)
#
ex.lat_gpats <- which(lat_gpats >= min_lat & lat_gpats <= max_lat)
ex.lon_gpats <- which(lon_gpats >= min_lon & lon_gpats <= max_lon)
lat_gpats <- lat_gpats[ex.lat_gpats]
lon_gpats <- lon_gpats[ex.lon_gpats]
total_lightning <- total_lightning[ex.lon_gpats,ex.lat_gpats,]
#cls.pred   <- cls.pred[ex.lon_era,ex.lat_era,]
#
ex.lat_lsmask <- which(lat_lsmask >= min_lat & lat_lsmask <= max_lat)
ex.lon_lsmask <- which(lon_lsmask >= min_lon & lon_lsmask <= max_lon)
lat_lsmask <- lat_lsmask[ex.lat_lsmask]
lon_lsmask <- lon_lsmask[ex.lon_lsmask]
lsmask <- lsmask[ex.lon_lsmask,ex.lat_lsmask]
#
##For ease of plotting
lat  <- lat_lsmask
lon  <- lon_lsmask
nlat <- length(lat)
nlon <- length(lon)

ntime_cls <- dim(cls.pred)[3]
ndays_cls <- ntime_cls/4.
nyrs_cls <- ndays_cls/365.25

ndays_gpats <- difftime(rtime_gpats[dim(total_lightning)[3]],
                        rtime_gpats[1])

cluster <- array(NA, dim=c(nlon,nlat,ntime_cls))
#cluster <- cls.pred
#dim(cls.pred) <- dim(cluster)

#Transpose for plotting
cluster[,,] <- t(apply(cls.pred[,,],1,rev))
lsmask[,] <- t(apply(lsmask,1,rev))
#Don't transpose the total lightning! It takes too long
#Only transpose later once calculations have been applied
#total_lightning <- t(apply(total_lightning,1,rev))
#ncls <- length(levels(as.factor(cluster)))

ncls <- length(levels(cls.pred))

#Apply the land-sea mask
#for (i in 1:ntime_cls){
#    cluster[,,i][which(lsmask==0)] <- NA
#}

#We can't apply the mask before transposing as the mask has been transposed.
#for (i in 1:ntime_gpats){
#    total_lightning[,,i][which(lsmask==0)] <- NA
#}


#Plot the cluster distribution across Australia
cat("Evaluating cluster distribution \n")
clus.count <- c()
clus.dist <- vector("list",length=ncls)
for (i in 1:ncls){
    clstemp <- cluster
    clus.count[i] <- length(which(clstemp==i))
    #ifelse(clstemp != i,NA,1)
    clstemp[clstemp != i] <- NA
    clstemp[clstemp == i] <- 1
    #clus.dist[[i]] <- apply(clstemp,c(1,2),sum,na.rm=TRUE)
    clus.dist[[i]] <- apply(clstemp,c(1,2),sum,na.rm=TRUE)/dim(cluster)[3]
    #clus.dist[[i]] <- apply(clstemp,c(1,2),mean,na.rm=TRUE)
}

#Apply the land sea mask to the clus.dist
for (i in 1:ncls){
    clus.dist[[i]][lsmask==0] <- NA
}

#Evaluate the precentage of each cluster
clus.perc <- clus.count*100./sum(clus.count)


ncols=50 #No of colors for plotting
pdf(paste0("Cluster_spatial_distribution_",ncls,"_clusters.pdf"),
    width=10,height=10,paper="special")
#par(mfrow=c(3,ncls/3))
par(mfrow=c(2,2))#20 clusters
#par(mar=c(4,5,3,2,1)+0.1)
par(mar=c(4,4,3,2)+0.1)
for (i in 1:ncls){
    image.plot(lon,rev(lat),clus.dist[[i]],
               main=paste0("Cluster: ",i),
               xlab="",ylab="",
               #xlab="Lon",ylab="Lat",
               #zlim=c(0,0.6),
               col=colorRampPalette(rev(brewer.pal(9,"Spectral")))(ncols),
               horizontal=TRUE)
    text(115,-15,paste(round(clus.perc[i],1),"%"),cex=1.25)
    oz(lwd=2.5,add=T)
}
dev.off()


prob.hail <- array(dim=c(nlon,nlat,ntime_cls))
prob.lightning <- array(dim=c(nlon,nlat,ntime_cls))
lightning.per.day <- array(dim=c(nlon,nlat,ntime_cls))
hail.days.per.year<- array(dim=c(nlon,nlat,ntime_cls))

#Assign a probability based on a cluster
#prob.hail[,,] <- hlprob$prop.hail[cls.pred[,,]]
#prob.lightning[,,] <- hlprob$prop.lightning[cls.pred[,,]]
prob.hail[,,] <- hlprob$prop.hail[cluster[,,]]
prob.lightning[,,] <- hlprob$prop.lightning[cluster[,,]]
lightning.per.day[,,] <- hlprob$prop.lightning.strikes[cluster[,,]]
#hail.days.per.year[,,] <- hlprob$days.hail[cluster[,,]]

#Evaluate the number of hail days per year/cluster=x
#Evaluate the number of time a cluster occurs=n
#the number of hail days per year in each cluster is x*n
#Then sum up each cluster

#image.plot(lon,rev(lat),cluster[,,1],nlevel=ncls)
pdf(paste0("Gap_Storm_",ncls,"_clusters.pdf"),width=10,height=10,paper="special")
image.plot(lon,rev(lat),cluster[,,gap_storm_doy*4],nlevel=ncls,
           main="Cluster map")
map('worldHires',interior=T,lwd=2.5,add=T)
oz(lwd=2.5,add=TRUE)

image.plot(lon,rev(lat),prob.hail[,,gap_storm_doy*4],nlevel=ncls,
           main="Hail probability map")
map('worldHires',interior=T,lwd=2.5,add=T)
oz(lwd=2.5,add=TRUE,coast=FALSE)

image.plot(lon,rev(lat),prob.lightning[,,gap_storm_doy*4],nlevel=ncls,
           main="Lightning probability map")
#map('worldHires',interior=T,lwd=2.5,add=T)
oz(lwd=2.5,add=TRUE)
dev.off()

#Now apply the sum

hp <- apply(prob.hail,c(1,2),sum)/ntime_cls
#hp <- (apply(prob.hail,c(1,2),sum)/1460)*4*365.
lp <- apply(prob.lightning,c(1,2),sum)/ntime_cls
#lpd <- apply(lightning.per.day,c(1,2),sum)*24/1460
#lpd <- apply(lightning.per.day,c(1,2),sum)*(24/4)/1460
#lpd <- (apply(lightning.per.day,c(1,2),sum)/1460)*4
#lpd <- (apply(lightning.per.day,c(1,2),sum)/1460)*4
#lpd <- (apply(lightning.per.day,c(1,2),sum)*24.)/(1460*4)
#lpd <- (apply(lightning.per.day,c(1,2),sum)*24.)/(1460*4)
#lpd0 <- roll_sum(lightning.per.day,4,by=4,align="left")
#lpd <- lpd0[lpd0 > 0]

#Evaluate the rolling sum of each fourth element
#Which corresponds to the sum of 00/06/12/18 predictions
#lpd0 <- apply(lightning.per.day,c(1,2),roll_sum,4,by=4,align="left",fill=NA)
#lpd0 <- apply(lightning.per.day,c(1,2),roll_mean,4,by=4,align="left",fill=NA)
lpd0 <- apply(lightning.per.day,c(1,2),roll_max,4,by=4,align="left",fill=NA)
#We need to extract every fourth element as the roll_sum routine leaves it in
#there
lpd0 <- lpd0[seq(1,dim(lpd0)[1],4),,]

#Now sum and average them over the year
#lpd <- apply(lpd0,c(2,3),sum,na.rm=T)/365.
#lpd <- apply(lpd0,c(2,3),sum,na.rm=T)/ntime_cls
lpd <- apply(lpd0,c(2,3),sum,na.rm=T)/ndays_cls

for (i in 1:ntime_cls){
    lpd[lsmask==0] <- NA
}

#Now the hail days per year
#hdpy <- apply(hail.days.per.year,c(1,2),sum)/ntime_cls
#lp <- lp/sum(lp)
#image.plot(lon,rev(lat),hp,nlevel=10,xlim=c(110,160),ylim=c(-45,-10),
#           col=rev(brewer.pal(10,"Spectral")))
#map('worldHires',interior=T,lwd=2.5,add=T)
#image.plot(lon,rev(lat),lp,nlevel=10,xlim=c(110,160),ylim=c(-45,-10))
#map('worldHires',interior=T,lwd=2.5,add=T)

pdf(paste0("Probability_of_lightning_hail_",ncls,"_clusters.pdf"),
    width=10,height=10,paper="special")
image.plot(lon,rev(lat),(lp*ndays_cls)/nyrs_cls,xlim=c(110,160),ylim=c(-45,-10),
           xlab="Lon",ylab="Lat",
           col=colorRampPalette(rev(brewer.pal(9,"Spectral")))(64),
           legend.lab="Lightning days/year",horizontal=TRUE)
oz(lwd=2,add=T)

image.plot(lon,rev(lat),(hp*ndays_cls)/nyrs_cls,xlim=c(110,160),ylim=c(-45,-10),
           xlab="Lon",ylab="Lat",
           #zlim=c(0,3),
           col=colorRampPalette(rev(brewer.pal(9,"Spectral")))(64),
           legend.lab="Hail days/year",horizontal=TRUE)
oz(lwd=2,add=T)
dev.off()


pdf(paste0("Predicted_lightning_strikes_per_day_",ncls,"_clusters.pdf"),
    width=10,height=10,paper="special")
image.plot(lon,rev(lat),lpd,xlim=c(110,160),ylim=c(-45,-10),
#image.plot(lon,rev(lat),nlevel=10,lpd/4,xlim=c(110,160),ylim=c(-45,-10),
           #zlim=c(0,22),
           #zlim=c(0,7),
           xlab="Lon",ylab="Lat",
           #col=rev(brewer.pal(10,"Spectral")))
           col=colorRampPalette(rev(brewer.pal(9,"Spectral")))(ncols),
           legend.lab="Lightning strikes per day",horizontal=TRUE)
#map('worldHires',interior=T,lwd=2.5,add=T)
oz(lwd=2.5,add=TRUE)
dev.off()

pdf(paste0("Predicted_hail_days_per_year_",ncls,"_clusters.pdf"),
    width=10,height=10,paper="special")
#image.plot(lon,rev(lat),hdpy,nlevel=ncls,xlim=c(110,160),ylim=c(-45,-10),
image.plot(lon,rev(lat),(hp*ndays_cls)/nyrs_cls,nlevel=ncls,xlim=c(110,160),ylim=c(-45,-10),
           xlab="Lon",ylab="Lat",
           #col=rev(brewer.pal(10,"Spectral")))
           col=colorRampPalette(rev(brewer.pal(9,"Spectral")))(64),
           #col=tim.colors(64),
           legend.lab="Hail days per year",horizontal=TRUE)
#map('worldHires',interior=T,lwd=2.5,add=T)
oz(lwd=2.5,add=TRUE)
dev.off()

#Now make a climatology of the GPATS data.
#To compare apples with apples we need to evaluate the maximum
#lightning frequency on any given day
#Finding the max in any 24-hour period
#I think we should use this as to construct the probabilities
#we use the max lightning at any time during the day.
max_lightning <- apply(total_lightning,c(1,2),roll_max,24,by=24,align="left",fill=NA)
#Find the max in window around each of the ERA time
#max_lightning <- apply(total_lightning,c(1,2),roll_max,3,by=6,
#                       align="center",fill=NA,na.rm=TRUE)
#Finding the mean around each of the ERA-time
#max_lightning <- apply(total_lightning,c(1,2),roll_mean,3,by=6,
#                       align="center",fill=NA,na.rm=TRUE)
#max_lightning <- apply(total_lightning,c(1,2),roll_sum,3,by=6,
#                       align="center",fill=NA,na.rm=TRUE)
#max_lightning <- apply(total_lightning,c(1,2),roll_mean,1,by=6,
#                       align="center",fill=NA,na.rm=TRUE)
#Extracting the lightning at each ERA time.
#max_lightning <- total_lightning[,,seq(1,dim(total_lightning)[3],6)]

#The mean lightning counts over the whole period.
mean_lightning <- apply(total_lightning,c(1,2),mean)

#Extract only every appropriate elements
#Extractine when we have found the max. lightning at any point during the day
max_lightning <- max_lightning[seq(1,dim(max_lightning)[1],24),,]
#Extracting when we have founf the max around each ERA time
#max_lightning <- max_lightning[seq(2,dim(max_lightning)[1],6),,]


#Now roll sum by 4 (00/06/12/18) to get daily total.
#We do this if we have extracted the lightning at the ERA times
#We use this if we have evaluated at the ERA time above
#max_lightning <-  apply(max_lightning,c(2,3),roll_sum,4,by=4,
#                        align="left",fill=NA,na.rm=TRUE)
#max_lightning <-  apply(max_lightning,c(1,2),roll_sum,4,by=4,
#                        align="left",fill=NA,na.rm=TRUE)
#max_lightning <- max_lightning[seq(1,dim(max_lightning)[1],4),,]
#max_lightning <- max_lightning[seq(1,dim(max_lightning)[1],24),,]
#Sum and average
mxlpdg <- apply(max_lightning,c(2,3),sum,na.rm=T)/(dim(max_lightning)[1])
#mxlpdg <- apply(max_lightning,c(2,3),sum,na.rm=T)/((dim(max_lightning)[1])*4)

#Transpose the arrays for plotting
mxlpdg <- t(apply(mxlpdg,1,rev))
mnlpdg <- t(apply(mean_lightning,1,rev))

#Apply the land sea mask
mxlpdgm <- mxlpdg
mnlpdgm <- mnlpdg
for (i in 1:length(time)){
    mxlpdgm[lsmask==0] <- NA
    mnlpdgm[lsmask==0] <- NA
}

#total_lightning[,,] <- t(apply(total_lightning[,,],1,rev))

pdf(paste0("GPATS_vs_cluster_lightning_per_day_",ncls,"_clusters.pdf"),
    width=10,height=10,paper="special")
image.plot(lon,rev(lat),mxlpdgm,xlim=c(110,160),ylim=c(-45,-10),
           main="GPATS observations of the average daily maximum flash count",
           zlim=c(0,max(c(lpd,mxlpdgm),na.rm=T)),
           xlab="Lon",ylab="Lat",
           #col=rev(brewer.pal(10,"Spectral")))
           col=colorRampPalette(rev(brewer.pal(9,"Spectral")))(64),
           #col=tim.colors(64),
           legend.lab="Lightning counts per day",horizontal=TRUE)
#map('worldHires',interior=T,lwd=2.5,add=T)
oz(lwd=2.5,add=TRUE)

image.plot(lon,rev(lat),lpd,xlim=c(110,160),ylim=c(-45,-10),
           main="Cluster prediction of the average daily maximum flash count",
           zlim=c(0,max(c(lpd,mxlpdgm),na.rm=T)),
           xlab="Lon",ylab="Lat",
           #col=rev(brewer.pal(10,"Spectral")))
           col=colorRampPalette(rev(brewer.pal(9,"Spectral")))(64),
           #col=tim.colors(64),
           legend.lab="Lightning counts per day",horizontal=TRUE)
oz(lwd=2.5,add=TRUE)

dev.off()

#image.plot(lon_gpats,rev(lat),mnlpdg,nlevel=10,xlim=c(110,160),ylim=c(-45,-10),
#           xlab="Lon",ylab="Lat",
#           #col=rev(brewer.pal(10,"Spectral")))
#           col=colorRampPalette(rev(brewer.pal(9,"Spectral")))(64),
#           #col=tim.colors(64),
#           legend.lab="LIghtning counts/year",horizontal=TRUE)
##map('worldHires',interior=T,lwd=2.5,add=T)
#oz(lwd=2.5,add=TRUE)

#Evaluate the proportion of the days in each grid box that have lightning
sum_lightning <- apply(total_lightning,c(1,2),roll_sum,24,by=24,align="left",fill=NA)
sum_lightning <- sum_lightning[seq(1,dim(sum_lightning)[1],24),,]
#lightning_in_grid <- which(sum_lightning > 0)
#tlgb - total lightning in grid box
tlgb <- apply(sum_lightning[,,],c(2,3),sum,na.rm=T)
tlgb <- t(apply(tlgb,1,rev))
#Now evaluate the proportion that each has a grid box has lightning
lightning_in_grid <- which(sum_lightning[,,] > 0)
plgbs <- sum_lightning[,,]
plgbs[lightning_in_grid] <- 1
plgb <- apply(plgbs,c(2,3),sum,na.rm=T)/dim(plgbs)[1]
plgb <- t(apply(plgb,1,rev))
#Apply the land sea mask
plgbm <- plgb
plgbm[lsmask==0] <- NA

#Best z lim for comparing cluster prediction and GPATS data
#zmax <- c(lp/dim(prob.lightning)[3],plgb)[which.max(c(lp/dim(prob.lightning)[3],plgb))]
#zmax <- c(lp,plgb)[which.max(c(lp,plgb))]
zmax <- c(lp,plgbm)[which.max(c(lp,plgb))]

#pdf(paste0("Cluster_vs_GPATS_proportion_lightning_2008_",ncls,"_cluster.pdf"),
pdf(paste0("Cluster_vs_GPATS_proportion_lightning_",ncls,"_cluster.pdf"),
    width=10,height=10,paper="special")
#par(mfrow=c(2,1))
par(mar=c(3,3,1,1))
#split.screen( rbind(c(0, .8,0,1), c(.8,1,0,1)))
#split.screen( rbind(c(0, .85,0,1), c(.85,1,0,1)))
split.screen( rbind(c(0, .9,0,1), c(.9,1,0,1)))
split.screen(c(2,1), screen=1)-> ind
screen( ind[1])
#image.plot(lon_gpats,rev(lat_gpats),plgb,xlim=c(110,160),ylim=c(-45,-10),
image(lon_gpats,rev(lat_gpats),plgbm,xlim=c(110,160),ylim=c(-45,-10),
           main="GPATS observations",
           xlab="Lon",ylab="Lat",
           zlim=c(0,zmax),
           #col=rev(brewer.pal(10,"Spectral")))
           col=colorRampPalette(rev(brewer.pal(9,"Spectral")))(ncols))
           #col=tim.colors(64),
           #legend.lab="Proportion of days with lightning",horizontal=TRUE)
#map('worldHires',interior=T,lwd=2.5,add=T)
oz(lwd=2.5,add=TRUE)

screen(ind[2])
#image.plot(lon,rev(lat),lp/dim(prob.lightning)[3],xlim=c(110,160),ylim=c(-45,-10),
#image(lon,rev(lat),lp/dim(prob.lightning)[3],xlim=c(110,160),ylim=c(-45,-10),
image(lon,rev(lat),lp,xlim=c(110,160),ylim=c(-45,-10),
           main="Cluster prediction",
           xlab="Lon",ylab="Lat",
           zlim=c(0,zmax),
           #col=rev(brewer.pal(10,"Spectral")))
           col=colorRampPalette(rev(brewer.pal(9,"Spectral")))(ncols))
           #col=tim.colors(64),
           #legend.lab="Proportion of days with lightning",horizontal=TRUE)
#map('worldHires',interior=T,lwd=2.5,add=T)
oz(lwd=2.5,add=TRUE)

screen(2)
#image.plot(zlim=c(0,zmax),legend.only=TRUE,smallplot=c(.85,.9, .1,.9),
#image.plot(zlim=c(0,zmax),legend.only=TRUE,smallplot=c(.2,.5, .2,.8),
image.plot(zlim=c(0,zmax),legend.only=TRUE,smallplot=c(.1,.2, .2,.8),
#image.plot(zlim=c(0,zmax),legend.only=TRUE,smallplot=c(.9,.95, .2,.8),
           col=colorRampPalette(rev(brewer.pal(9,"Spectral")))(ncols))

dev.off()

pdf(paste0("Difference_cluster_GPATS_",ncls,"_cluster.pdf"),
    width=10,height=10,paper="special")
diff <- lp-plgbm
#See http://stackoverflow.com/questions/33750235/plotting-a-raster-with-the-color-ramp-diverging-around-zero
col_levs <- 50
max_abs_val <- max(abs(c(min(diff,na.rm=T),max(diff,na.rm=T))))
col_seq <- seq(-1.*max_abs_val,max_abs_val,length.out=col_levs+1)
image.plot(lon,rev(lat),diff,xlim=c(110,160),ylim=c(-45,-10),
#image.plot(lon,rev(lat),diff/plgbm,xlim=c(110,160),ylim=c(-45,-10),
        main="Cluster prediction minus GPATS",
        xlab="Lon",ylab="Lat",
        #zlim=c(0,zmax),
        #col=colorRampPalette(rev(brewer.pal(9,"Spectral")))(ncols))
        #col=colorRampPalette(rev(brewer.pal(9,"RdBu")))(ncols))
        col=colorRampPalette(rev(brewer.pal(9,"RdBu")))(col_levs),
        breaks=col_seq)
oz(lwd=2.5,add=TRUE)

#Difference as a fraction of the observations
prop <- diff/plgbm
prop[prop > 1] <- NA
max_abs_val <- max(abs(c(min(prop,na.rm=T),max(prop,na.rm=T))))
col_seq <- seq(-1.*max_abs_val,max_abs_val,length.out=col_levs+1)
image.plot(lon,rev(lat),prop,xlim=c(110,160),ylim=c(-45,-10),
#image.plot(lon,rev(lat),diff/plgbm,xlim=c(110,160),ylim=c(-45,-10),
        main="Proportionate difference between cluster and GPATS",
        xlab="Lon",ylab="Lat",
        #zlim=c(-1.,1),
        #col=colorRampPalette(rev(brewer.pal(9,"Spectral")))(ncols))
        #col=colorRampPalette(rev(brewer.pal(9,"RdBu")))(ncols))
        col=colorRampPalette(rev(brewer.pal(9,"RdBu")))(col_levs),
        breaks=col_seq)
oz(lwd=2.5,add=TRUE)
dev.off()


#Now do a seasonal analysis
#Get the cluster data in a daily average
#First we need to extract the GPATS time which correspond to the ERA times
library(lubridate)
rtime <- rtime_gpats[seq(1,ntime_gpats,6)]

yr.2008 <- which(year(rtime) == 2008)
djf <- c(12,1,2)
mam <- c(3,4,5)
jja <- c(6,7,8)
son <- c(9,10,11)

summer <- which(month(rtime) %in% djf)
autumn <- which(month(rtime) %in% mam)
winter <- which(month(rtime) %in% jja)
spring <- which(month(rtime) %in% son)

#First the cluster prediction
lpd0.summer <- apply(lightning.per.day[,,yr.2008[summer]],
                     c(1,2),roll_max,4,by=4,align="left",fill=NA)
lpd0.autumn <- apply(lightning.per.day[,,yr.2008[autumn]],
                     c(1,2),roll_max,4,by=4,align="left",fill=NA)
lpd0.winter <- apply(lightning.per.day[,,yr.2008[winter]],
                     c(1,2),roll_max,4,by=4,align="left",fill=NA)
lpd0.spring <- apply(lightning.per.day[,,yr.2008[spring]],
                     c(1,2),roll_max,4,by=4,align="left",fill=NA)
#We need to extract every fourth element as the roll_sum routine leaves it in
#there
lpd0.summer <- lpd0.summer[seq(1,dim(lpd0)[1],4),,]
lpd0.autumn <- lpd0.autumn[seq(1,dim(lpd0)[1],4),,]
lpd0.winter <- lpd0.winter[seq(1,dim(lpd0)[1],4),,]
lpd0.spring <- lpd0.spring[seq(1,dim(lpd0)[1],4),,]

#Now sum and average them over the year
#lpd <- apply(lpd0,c(2,3),sum,na.rm=T)/365.
#lpd <- apply(lpd0,c(2,3),sum,na.rm=T)/ntime_cls
lpd.summer <- apply(lpd0.summer,c(2,3),sum,na.rm=T)/ndays_cls
lpd.autumn <- apply(lpd0.autumn,c(2,3),sum,na.rm=T)/ndays_cls
lpd.winter <- apply(lpd0.winter,c(2,3),sum,na.rm=T)/ndays_cls
lpd.spring <- apply(lpd0.spring,c(2,3),sum,na.rm=T)/ndays_cls

#for (i in 1:ntime_cls){
    lpd.summer[lsmask==0] <- NA
    lpd.autumn[lsmask==0] <- NA
    lpd.winter[lsmask==0] <- NA
    lpd.spring[lsmask==0] <- NA
#}

pdf(paste0("Cluster_seasonal_cycle",ncls,"_cluster.pdf"),
    width=10,height=10,paper="special")
image.plot(lon,rev(lat),lpd.summer,xlim=c(110,160),ylim=c(-45,-10),
        main="Cluster seasonal cycle",
        xlab="Lon",ylab="Lat",
        col=colorRampPalette(rev(brewer.pal(9,"Spectral")))(ncols))
oz(lwd=2.5,add=TRUE)

image.plot(lon,rev(lat),lpd.autumn,xlim=c(110,160),ylim=c(-45,-10),
        main="Cluster seasonal cycle",
        xlab="Lon",ylab="Lat",
        col=colorRampPalette(rev(brewer.pal(9,"Spectral")))(ncols))
oz(lwd=2.5,add=TRUE)

image.plot(lon,rev(lat),lpd.winter,xlim=c(110,160),ylim=c(-45,-10),
        main="Cluster seasonal cycle",
        xlab="Lon",ylab="Lat",
        col=colorRampPalette(rev(brewer.pal(9,"Spectral")))(ncols))
oz(lwd=2.5,add=TRUE)

image.plot(lon,rev(lat),lpd.spring,xlim=c(110,160),ylim=c(-45,-10),
        main="Cluster seasonal cycle",
        xlab="Lon",ylab="Lat",
        col=colorRampPalette(rev(brewer.pal(9,"Spectral")))(ncols))
oz(lwd=2.5,add=TRUE)


dev.off()



#Now the GPATS observations


##lpd0 <- apply(lightning.per.day,c(1,2),roll_mean,4,by=4,align="left",fill=NA)

#This is the two-hourly averaged lightning at the ERA times.
#Since the GPATS data has been offset so that the lightning count is the 
#number measured during over the hour following the time signature,
#we will include the previous, so we have summed the hour
#preceding and following the ERA time
#mngpats <- apply(total_lightning,c(1,2),roll_mean,2,by=6,align="right",fill=NA)
#mxgpats  <- apply(total_lightning,c(1,2),roll_max,2,by=6,align="right",fill=NA)
#mxgpats <- mxgpats[seq(2,dim(mxgpats)[1],6),,]
mxgpats  <- apply(total_lightning,c(1,2),roll_max,24,by=24,align="left",fill=NA)
mxgpats <- mxgpats[seq(1,dim(mxgpats)[1],24),,]

#mxgpats0.summer <- apply(mxgpats[yr.2008[summer],,],
#                        c(2,3),roll_max,4,by=4,align="left",fill=NA)
#mxgpats0.autumn <- apply(mxgpats[yr.2008[autumn],,],
#                        c(2,3),roll_max,4,by=4,align="left",fill=NA)
#mxgpats0.winter <- apply(mxgpats[yr.2008[winter],,],
#                        c(2,3),roll_max,4,by=4,align="left",fill=NA)
#mxgpats0.spring <- apply(mxgpats[yr.2008[spring],,],
#                        c(2,3),roll_max,4,by=4,align="left",fill=NA)

mxgpats0.summer <- apply(mxgpats[yr.2008[summer],,],
                         c(2,3),sum,na.rm=TRUE)

mxgpats0.autumn <- apply(mxgpats[yr.2008[autumn],,],
                         c(2,3),sum,na.rm=TRUE)
mxgpats0.winter <- apply(mxgpats[yr.2008[winter],,],
                         c(2,3),sum,na.rm=TRUE)
mxgpats0.spring <- apply(mxgpats[yr.2008[spring],,],
                         c(2,3),sum,na.rm=TRUE)

#mxgpats0.summer <- mxgpats0.summer[seq(1,dim(mxgpats0.summer)[1],4),,]
#mxgpats0.autumn <- mxgpats0.autumn[seq(1,dim(mxgpats0.autumn)[1],4),,]
#mxgpats0.winter <- mxgpats0.winter[seq(1,dim(mxgpats0.winter)[1],4),,]
#mxgpats0.spring <- mxgpats0.spring[seq(1,dim(mxgpats0.spring)[1],4),,]

#Now sum and average
#mxgpats.summer <- apply(mxgpats0.summer,c(2,3),sum,na.rm=TRUE)/dim(mxgpats0.summer)
#mxgpats.autumn <- apply(mxgpats0.autumn,c(2,3),sum,na.rm=TRUE)/dim(mxgpats0.autumn)
#mxgpats.winter <- apply(mxgpats0.winter,c(2,3),sum,na.rm=TRUE)/dim(mxgpats0.winter)
#mxgpats.spring <- apply(mxgpats0.spring,c(2,3),sum,na.rm=TRUE)/dim(mxgpats0.spring)

mxgpats.summer <- mxgpats0.summer/length(yr.2008[summer])
mxgpats.autumn <- mxgpats0.autumn/length(yr.2008[autumn])
mxgpats.winter <- mxgpats0.winter/length(yr.2008[winter])
mxgpats.spring <- mxgpats0.spring/length(yr.2008[spring])

#Transpose and apply mask
mxgpats.summer <- t(apply(mxgpats.summer,1,rev))
mxgpats.autumn <- t(apply(mxgpats.autumn,1,rev))
mxgpats.winter <- t(apply(mxgpats.winter,1,rev))
mxgpats.spring <- t(apply(mxgpats.spring,1,rev))
#Apply the land sea mask

mxgpats.summer[lsmask==0] <- NA
mxgpats.autumn[lsmask==0] <- NA
mxgpats.winter[lsmask==0] <- NA
mxgpats.spring[lsmask==0] <- NA

#Now plot all of them together
#par(mar=c(3,3,1,1))
#split.screen( rbind(c(0, .8,0,1), c(.8,1,0,1)))
##split.screen( rbind(c(0, .85,0,1), c(.85,1,0,1)))
#split.screen( rbind(c(0, .9,0,1), c(.9,1,0,1)))
#split.screen(c(4,2), screen=1)-> ind
#screen( ind[1])
#
#image.plot(lon,rev(lat),mxgpats.summer)



#Now apply the mask
#for (i in 1:
