#Script to reproduce the AIRS analysis
#See https://ams.confex.com/ams/27SLS/webprogram/.../extendedAbstract_pdf.pdf
#Justin Peter
#ICACS, USQ, 25 October 2016

rm(list=ls())
gc()

library(ncdf4)
library(fields)
library(maps)
library(mapdata)
library(oz)
library(RColorBrewer)
library(oce)
library(RANN) #For nearest neighbour routine
library(RcppRoll) #For rolling min/max/sum etc.

#----------Read in the data
#ERA, GPATS and BoM hail reports
fnm_era <- '/home/jpeter/justinp/rfiles/suncorp/data/eraint/sharpy_profiles/00_18/2008_2014/era_profiles_2008_2014_00_18.nc'

#GPATS lightning data
fnm_gpats <- '/home/jpeter/justinp/rfiles/suncorp/data/gpats/ts_lightning_grid/2008_2015/combined/lightning_eragrid_eratime_Australia_2008_2015_offset.nc'

#BoM hail reports
fnm_bom <- '/home/jpeter/justinp/rfiles/suncorp/data/bom_hail/reports.txt'


#The land sea mask
fnm_lsmask <- '/home/jpeter/justinp/rfiles/suncorp/data/eraint/netcdf_download/lsmask/era_interim_lsmask.nc'

#Read the hail reports data
hail_data <- read.table(fnm_bom)
hail_data[[1]] <- as.Date(hail_data[[1]])
#hail_data[[1]] <- as.POSIXct(hail_data[[1]],tz="UTC")
names(hail_data) <- c("rdate_hail", "lon","lat")

#Open the netcdf files
nc_era   <- nc_open(fnm_era)
nc_gpats <- nc_open(fnm_gpats)
nc_lsmask <- nc_open(fnm_lsmask)

print(paste("SHARPPy profile file has",nc_era$nvars,"variables"))
print(paste("GPATS lightning file has",nc_gpats$nvars,"variables"))

#Get the dimensions of the ERA file
nlat_era    <- nc_era$dim$lat$len
nlon_era    <- nc_era$dim$lon$len
ntime_era   <- nc_era$dim$time$len

lat_era  <- nc_era$dim$lat$vals
lon_era  <- nc_era$dim$lon$vals
time_era <- nc_era$dim$time$vals
rtime_era   <- as.POSIXct(time_era*60*60, origin="1900-01-01", tz="GMT")#in UTC
rdate_era   <- as.Date(rtime_era,tz="UTC")

#Get the dimensions of the GPATS file
nlat_gpats  <- nc_gpats$dim$lat$len
nlon_gpats  <- nc_gpats$dim$lon$len
ntime_gpats <- nc_gpats$dim$time$len

lat_gpats  <- nc_gpats$dim$lat$vals
lon_gpats  <- nc_gpats$dim$lon$vals
time_gpats <- nc_gpats$dim$time$vals
rtime_gpats<- as.POSIXct(time_gpats*60*60, origin="1900-01-01", tz="GMT")
rdate_gpats <- as.Date(rtime_gpats,tz="UTC")


#Get the variables
#SHIP
v1<-nc_era$var[["ship"]]
#ship.raw <- ncvar_get(nc_era,v1)
ship <- ncvar_get(nc_era,v1)

#S06 and CAPE
v4 <- nc_era$var[["sfc_6km_shear"]]
s06 <- ncvar_get(nc_era,v4) #Is in knots
s06 <- s06*0.514444 #Convert to metres/second
v5 <- nc_era$var[["bplus"]]
cape <- ncvar_get(nc_era,v5)

nc_close(nc_era)

#Get the total lightning variable
v2<-nc_gpats$var[["total_lightning"]]
total_lightning <- ncvar_get(nc_gpats,v2)
nc_close(nc_gpats)

#Land sea mask
v3 <- nc_lsmask$var[["lsm"]]
lsmask <- ncvar_get(nc_lsmask,v3)
lat_lsmask <- nc_lsmask$dim$lat$vals
lon_lsmask <- nc_lsmask$dim$lon$vals
nc_close(nc_lsmask)
#Transpose the land-sea mask
#lsmask <- t(apply(lsmask,1,rev))

#Only use the smallest domain
##Only use the spatial of the smallest file (GPATS)
min_lat <- min(lat_gpats)
max_lat <- max(lat_gpats)
min_lon <- min(lon_gpats)
max_lon <- max(lon_gpats)
#
ex.lat_era <- which(lat_era >= min_lat & lat_era <= max_lat)
ex.lon_era <- which(lon_era >= min_lon & lon_era <= max_lon)
lat_era <- lat_era[ex.lat_era]
lon_era <- lon_era[ex.lon_era]

ex.lat_lsmask <- which(lat_lsmask >= min_lat & lat_lsmask <= max_lat)
ex.lon_lsmask <- which(lon_lsmask >= min_lon & lon_lsmask <= max_lon)
lat_lsmask <- lat_lsmask[ex.lat_lsmask]
lon_lsmask <- lon_lsmask[ex.lon_lsmask]
lsmask <- lsmask[ex.lon_lsmask,ex.lat_lsmask]
#
##For ease of plotting
lat  <- lat_gpats
lon  <- lon_gpats
nlat <- length(lat)
nlon <- length(lon)

#Extract the GPATS domain from the ERA data
ship <- ship[ex.lon_era,ex.lat_era,]
s06 <- s06[ex.lon_era,ex.lat_era,]
cape <- cape[ex.lon_era,ex.lat_era,]

#Only use the union of times
#min.date <- min(rdate_era)
#max.date <- max(rdate_era)
min.date <- "2008-01-01"
min.date.hail <- "2008-01-01"
#max.date <- "2014-12-31"
max.date <- "2012-12-31"
max.date.hail <- "2012-12-31"
nyears <- as.double(difftime(max.date,min.date)/365.25)

ex.rtime_era <- which(rdate_era >= min.date &
                      rdate_era <= max.date)
ex.rtime_gpats <- which(rdate_gpats >= min.date &
                        rdate_gpats <= max.date)
ex.rtime_hail <- which(hail_data$rdate_hail >= min.date.hail &
                       hail_data$rdate_hail <= max.date.hail)

ship <- ship[,,ex.rtime_era]
cape <- cape[,,ex.rtime_era]
s06 <- s06[,,ex.rtime_era]
total_lightning <- total_lightning[,,ex.rtime_gpats]
hail_data <- hail_data[ex.rtime_hail,]


##Apply the appropriate modifications to the SHIP data
##The AIR paper uses s06 >= 10 m/s and SHIP is capped at max of 2000 J/kg
#cape.thresh <- 2000.
#shear.thresh <- 10
cape.thresh <- 2000.
shear.thresh <- 10
s06.lt10 <- which(s06 < shear.thresh)
ship[s06.lt10] <- NA
cape.gt2000 <- which(cape > cape.thresh)
cape.lt2000 <- which(cape <= cape.thresh)
#delta.cape <- array(0,dim=c(nlon,nlat,length(rtime_era)))
#delta.ship <- array(0,dim=c(nlon,nlat,length(rtime_era)))
##Evaluate the proportional decrease in CAPE
##to reduce it to 2000 J/kg.
#delta.cape[cape.gt2000] <- (cape[cape.gt2000]-2000.)/2000.
##delta.cape <- (cape[cape.gt2000]-2000.)/2000.
#delta.cape[cape.lt2000] <- 1.0
##Decrease the SHIP appropriately
##delta.ship[cape.gt2000]<-ship[cape.gt2000]-(delta.cape[cape.gt2000]*ship[cape.gt2000])
#delta.ship[cape.gt2000] <- (delta.cape[cape.gt2000]*ship[cape.gt2000])
#delta.ship[cape.lt2000]<-1.0

#sc.ship is the scaled SHIP. sc is the scaling array
#sc <- array(1,dim=c(nlon,nlat,length(rtime_era)))
sc <- array(1,dim=dim(ship))
#sc.ship <- array(0,dim=c(nlon,nlat,length(rtime_era)))
sc[cape.gt2000] <- cape.thresh/cape[cape.gt2000]
sc.ship <- sc*ship

ship <- sc.ship

#Now extract the grids corresponding to the lat,lon grid

delta_era <- 0.75/2
ex.era_grid <- read.table('/home/jpeter/justinp/rfiles/suncorp/data/bom_hail/lonlat.txt')
names(ex.era_grid) <- c('lon','lat')

ex.lon <- ex.era_grid$lon
ex.lat <- ex.era_grid$lat
ngrids <- length(ex.era_grid$lon)

ex.lon.inds <- vector("numeric",length=ngrids)
ex.lat.inds <- vector("numeric",length=ngrids)
#Extract the indices of the lat and lon
for (i in 1:ngrids){
    ex.lon.inds[i] <- which(lon == ex.lon[i])
    ex.lat.inds[i] <- which(lat == ex.lat[i])
}


#Extract the corresponding hail_data indices for each location.
ex.hail_data.inds <- vector("list",length=ngrids)
for (i in 1:ngrids){
    ex.hail_data.inds[[i]] <- which(hail_data$lon >= ex.lon[i]-delta_era &
                           hail_data$lon <= ex.lon[i]+delta_era &
                           hail_data$lat >= ex.lat[i]-delta_era &
                           hail_data$lat <= ex.lat[i]+delta_era)
}

#Assign the data using the above indices.
ex.hail_data <- vector("list",length=ngrids)
#for (i in 1:length(ngrids)){
for (i in 1:ngrids){
    ex.hail_data[[i]] <- hail_data[ex.hail_data.inds[[i]],]
}


#Extract each grid location and evaluate how many
#hail days occurred.
#ex_grid_hail_data  <- vector("list",length=ngrids)
#ex_grid_ship_data  <- vector("list",length=ngrids)
#ex_grid_gpats_data <- vector("list",length=ngrids)
#The following is wrong. It only contains a subset of the points.
#i.e. length(ex_grid_hail_data[[1]])=
#      length(ex_grid_hail_data[[1]]%in%ex.hail_data.inds[[1]])
#for (i in 1:ngrids){
#    ex_grid_hail_data[[i]] <- which(hail_data$lon <= grid_max_lon[i] & 
#                               hail_data$lon >= grid_min_lon[i] & 
#                               hail_data$lat >= grid_min_lat[i] & 
#                               hail_data$lat <= grid_max_lat[i])
#}

#Find the unique days of hail observations at each location.
days_grid_hail <- vector("list",length=ngrids)
for (i in 1:ngrids){
    #days_grid_hail[[i]]<-length(unique(hail_data$rdate_hail[ex_grid_hail_data[[i]]]))
    #days_grid_hail[[i]]<-unique(hail_data$rdate_hail[ex_grid_hail_data[[i]]])
    days_grid_hail[[i]] <- unique(hail_data$rdate_hail[ex.hail_data.inds[[i]]])
    #The following counts multiple hail reports in a grid box on a day.
    #days_grid_hail[[i]] <- hail_data$rdate_hail[ex.hail_data.inds[[i]]]
}

#Find the intersection of dates in BoM/ERA/GPATS.
days_hail_in_era   <- vector("list",length=ngrids)
days_hail_in_gpats <- vector("list",length=ngrids)
for (i in 1:ngrids){
    days_hail_in_era[[i]]   <- which(rdate_era %in% days_grid_hail[[i]])
    days_hail_in_gpats[[i]] <- which(rdate_gpats %in% days_grid_hail[[i]])
    #days_hail_in_era[[i]]   <- which(days_grid_hail[[i]] %in% rdate_era)
    #days_hail_in_gpats[[i]] <- which(days_grid_hail[[i]] %in% rdate_gpats)
}

#Extract grids from the ship and lightning data
ex.ship_data <- vector("list",length=ngrids)
ex.gpats_data <- vector("list",length=ngrids)
for (i in 1:ngrids){
  #ex.ship_data[[i]] <-ship[ex.lon.inds[i],ex.lat.inds[i],days_hail_in_era[[i]]]
  #ex.gpats_data[[i]]<- total_lightning[ex.lon.inds[i],ex.lat.inds[i],
  #                                    days_hail_in_gpats[[i]]]
  ex.ship_data[[i]] <-ship[ex.lon.inds[i],ex.lat.inds[i],]
  ex.gpats_data[[i]]<- total_lightning[ex.lon.inds[i],ex.lat.inds[i],]
}

#We evaluate the max ship and lightning
#in a day. These have length(4) and length(24) respectively.
max.ship_data <- vector("list",length=ngrids)
max.gpats_data <- vector("list",length=ngrids)
for (i in 1:ngrids){
    cat(i,'\n')
    if(length(ex.gpats_data[[i]] >0)){
    max.ship_data[[i]]<-roll_max(ex.ship_data[[i]],4,by=4,na.rm=T)
    max.gpats_data[[i]]<-roll_max(ex.gpats_data[[i]],24,by=24,na.rm=T)
    max.ship_data[[i]]<-max.ship_data[[i]][seq(1,length(max.ship_data[[i]]),4)]
    max.gpats_data[[i]]<-max.gpats_data[[i]][seq(1,length(max.gpats_data[[i]]),24)]
    }
}
 
#Loop over each of those days and evaluate:
#When SHIP > 0.5 and lightning occurred
ndays.hail <- vector("numeric",length=ngrids)
ndays.int.ship.gpats <- rep(0,length=ngrids)

ship.thresh <- 0.5
gpats.thresh <- 1
for (i in 1:ngrids){
    ndays.int.ship.gpats[[i]] <- length(which(max.ship_data[[i]] > ship.thresh & 
                                   #ex.gpats_data[[i]] > gpats.thresh))
                                   max.gpats_data[[i]] > gpats.thresh))
}

#for (i in 1:ngrids){
#    if(length(ex.gpats_data[[i]]) > 0){
#    for (j in 1:length(max.ship_data[[i]])){
#        if (max.ship_data[[i]][j] > 0.5 & max.gpats_data[[i]][j] > 0){
#            ndays.int.ship.gpats[i] = ndays.int.ship.gpats[i]+1
#        }
#    }
#    }
#    cat(i,'\n')
#}
#
for (i in 1:ngrids){
    #ndays.hail[i] <- length(days_hail_in_era[[i]])/4
    ndays.hail[i] <- length(days_grid_hail[[i]])
}

#for (i in 1:ngrids){
#   max.ship_data[[i]]<-max.ship_data[[i]][which(max.ship_data[[i]] > 0.5)]
#   max.gpats_data[[i]]<-max.gpats_data[[i]][which(max.gpats_data[[i]] > 0)]
#}

#Now evaluate the no. of hail days as a function of SHIP>0.5 and lightning
#model <- lm(ndays.hail ~ ndays.int.ship.gpats-1)
#model <- lm(ndays.hail ~ ndays.int.ship.gpats-1)
#USe a threshold for the number of hail days
hail.thresh0 <- which(ndays.hail >= 0)
hail.thresh1 <- which(ndays.hail >= 1)
hail.threshn <- which(ndays.hail >= nyears)
model0 <- lm(ndays.hail[hail.thresh0] ~ 
             ndays.int.ship.gpats[hail.thresh0]-1)
r20 <- summary(model0)$r.squared
model1 <- lm(ndays.hail[hail.thresh1] ~ 
             ndays.int.ship.gpats[hail.thresh1]-1)
r21 <- summary(model1)$r.squared
modeln <- lm(ndays.hail[hail.threshn] ~ 
            ndays.int.ship.gpats[hail.threshn]-1)
r2n <- summary(modeln)$r.squared

#Create the AIR model
x_airs <- c(1.47,4.62,2.54,16.12,18.60)
y_airs <- c(3.83,6.75,18.07,41.24,48.18)
airs_model <- lm(y_airs ~ x_airs-1)
r2a <- summary(airs_model)$r.squared


#Set the model to be used for the prediction
model <- model1
#model <- airs_model
#hail.thresh <- hail.threshn
hail.thresh <- hail.thresh0


#model <- lm(ndays.hail[hail.gt0] ~ ndays.int.ship.gpats[hail.gt0])
#plot(ndays.int.ship.gpats,ndays.hail)


png("air_regression.png",width=15,height=15,units="cm",res=300)
par(cex=1.2)
par(mar=c(4,4,1,1))
plot(ndays.int.ship.gpats[hail.thresh],ndays.hail[hail.thresh],
     xlab="Days of lightning and SHIP >0.5",
     ylab="Days of reported hail",
     lwd=3,
     xlim=c(0,100),ylim=c(0,50))
#model <- lm(ndays.hail[c(1,3,4,5)] ~ ndays.int.ship.gpats[c(1,3,4,5)]-1)
#plot(ndays.int.ship.gpats[c(1,3,4,5)],ndays.hail[c(1,3,4,5)])
#abline(a=0,b=model$coefficients[["ndays.int.ship.gpats"]])
points(x_airs,y_airs,lwd=3,col="red")
abline(a=0,b=model0$coefficients[["ndays.int.ship.gpats[hail.thresh0]"]],lwd=3)
abline(a=0,b=model1$coefficients[["ndays.int.ship.gpats[hail.thresh1]"]],lwd=3,lty="dashed")
abline(a=0,b=modeln$coefficients[["ndays.int.ship.gpats[hail.threshn]"]],lwd=3,lty="twodash")
abline(a=0,b=airs_model$coefficients[["x_airs"]],lwd=3,col="red")
legend("topright",
       c(paste0("Hail days >= 0; ", format(r20,digits=2)),
         paste0("Hail days >= 1; ",format(r21,digits=2)),
         paste0("At least one hail event per year; ", format(r2n,digits=2)),
         paste0("AIR regression; ", format(r2a,digits=2))),
       lwd=3,
       cex=0.8,
       lty=c("solid","dashed","twodash","solid"),
       col=c("black","black","black","red"))
box()
dev.off()

#Now use only lightning data from 2008-2014
#co.dates <- which(rdate_gpats <= max(rdate_era))
#rdate_gpats <- rdate_gpats[co.dates]
#total_lightning <- total_lightning[,,co.dates]

#Find the max ship and lightning on a day
#max.ship <- apply(ship,c(1,2),roll_max,4,by=4,align="left",fill=NA)
#max.ship <- apply(ship,c(1,2),roll_max,4,by=4,align="left",fill=NA,na.rm=T)
max.ship <- apply(ship,c(1,2),roll_max,4,by=4,na.rm=T)
max.ship <- max.ship[seq(1,dim(max.ship)[1],4),,]
#max.gpats <- apply(total_lightning,c(1,2),roll_max,24,by=24,align="left",fill=NA)
#max.gpats <-apply(total_lightning,c(1,2),roll_max,24,by=24,align="left",fill=NA,na.rm=T)
max.gpats <-apply(total_lightning,c(1,2),roll_max,24,by=24,na.rm=T)
max.gpats <- max.gpats[seq(1,dim(max.gpats)[1],24),,]

#Determine which SHIP obs are greater than 0.5
#and set them to 1. Less than 0.5=NA
sh <- which(max.ship > ship.thresh)
shn <- which(max.ship <= ship.thresh)

max.ship[sh] <- 1
max.ship[shn] <- NA

#Determine which grids have lightning
#and set equal to 1.
gp <- which(max.gpats > gpats.thresh)
gpn <- which(max.gpats <= gpats.thresh)

max.gpats[gp] <- 1
max.gpats[gpn] <- NA

nt <- dim(max.ship)[1] #ntimes
lst <- vector("list",length=nt)
mult <- vector("list",length=nt)
for (i in 1:nt){
    #Multiply the ship and gpats arrays together
    #If the result is != NA then the conditions of SHIP>0.5 and 
    #lightning present have been satisfied.
    lst[[i]] <- list(matrix(max.ship[i,,],nrow=nlon),matrix(max.gpats[i,,],nrow=nlon))
    mult[[i]] <- Reduce("*",lst[[i]])
}

#Now sum up each the grid boxes.
#Can't use Reduce as it cant handle NA
#sm <- Reduce("+",mult)
sm <- apply(simplify2array(mult),c(1,2),sum,na.rm=T)
sm[lsmask==0] <- NA
#sm[lsmask==0] <- 0
#Transpose array for plotting and apply land-sea mask
smt <- t(apply(sm,1,rev))
#smt[lsmask==0] <- NA

image.plot(lon,rev(lat),smt,
           #col=colorRampPalette(rev(brewer.pal(9,"Spectral")))(64),
           col=oce.colorsViridis(20),
           legend.lab="Sum of intersection of SHIP > 0.5 and lightning",
           horizontal=TRUE)
oz(add=TRUE)
box()

##Now apply the regression equation

nhd <- model$coefficients[1]*smt/nyears
#nhd <- (model$coefficients[["ndays.int.ship.gpats[hail.gt0]"]]*smt+


#png("airs_hail_days.png",width=600,height=600)
#
##        model$coefficients[["(Intercept)"]])/(nyears)
#image.plot(lon,rev(lat),nhd,
#           #col=colorRampPalette(rev(brewer.pal(9,"Spectral")))(64),
#           #col=oce.colorsViridis(64),
#           col=oce.colorsViridis(20),
#           legend.lab="Hail days/year",
#           horizontal=TRUE)
#oz(add=TRUE)
#box()
#dev.off()

#Evaluate the standard deviation
#Do this by evaluating number of storms for each year.
ndays_by_year <- c(366,365,365,365,366)
stride <- c(0,cumsum(ndays_by_year))
sm_by_year <- vector("list",length=length(ndays_by_year)) #Sum over all days in year
smt_by_year <- vector("list",length=length(ndays_by_year)) #Transpose for plotting

for (i in 1:(length(ndays_by_year))){
    sm_by_year[[i]] <- apply(simplify2array(mult[(stride[i]+1):(stride[i+1])]),c(1,2),sum,na.rm=T)
    sm_by_year[[i]][lsmask==0] <- NA
}


for (i in 1:length(sm_by_year)){
    smt_by_year[[i]] <- t(apply(sm_by_year[[i]],1,rev))
}

nhd_by_year <- vector("list",length(ndays_by_year))
for (i in 1:length(ndays_by_year)){
    nhd_by_year[[i]] <- model$coefficients[1]*smt_by_year[[i]]
}

#A check. This should be the same as sm, the sum of the intersection
#of ship and lightning over all days. sm2 is summing up the days in each year
#individually and then summing over all the years.
sm2 <- apply(simplify2array(smt_by_year),c(1,2),sum,na.rm=T)
#sm2t <- t(apply(sm2,1,rev))
#Have already taken transpose so need to use transpose of mask
lsmaskt <- t(apply(lsmask,1,rev))
sm2[lsmaskt==0] <- NA

nhd2 <- apply(simplify2array(nhd_by_year),c(1,2),sum,na.rm=T)/nyears
nhd2[lsmaskt==0] <- NA
sd.nhd2 <- apply(simplify2array(nhd_by_year),c(1,2),sd,na.rm=T)
sd.nhd2[lsmaskt==0] <- NA


#pdf("airs_hail_days.pdf",width=10,height=10)
png("airs_hail_days.png",width=10,height=10,units="in",res=300)
#par(mar=c(5.5,3,3,2)+0.1)
par(mar=c(2,3,3,0)+0.1)
par(oma=c(0.1,0.1,0.1,0.1))
#par(mar=c(7.0,3,5,2)+0.5)
par(cex=1.75)
image.plot(lon,rev(lat),nhd,
           #main="Hail days per year",xlab="",ylab="",
           xlab="",ylab="",
           xaxt="n",yaxt="n",
           #cex.main=0.9,
           #col=colorRampPalette(rev(brewer.pal(9,"Spectral")))(64),
           #col=oce.colorsViridis(64),
           #col=oce.colorsViridis(20),
           #col=tim.colors(20),
           col=tim.colors(12),
           #col=tim.colors(180),
           #legend.lab="Hail days/year",
           #horizontal=TRUE, cex.axis=1.25,cex.lab=1.25)
           horizontal=TRUE,
           breaks=seq(0,3,by=0.25))
           #breaks=seq(0,45,by=0.25))
           #legend.args=list(text="Hail days/year",cex=1.25))
           #legend.args=list(text="Hail days/year",cex=1.25))
title("Hail days per year",line=0.5,cex.main=2.0)
text(155,-12.5,"a",cex=3.0)
oz(add=TRUE)
map.cities(minpop=1.e06,cex=2.0)
box(lwd=3.0)
dev.off()

#png("airs_hail_days_stddev.png",width=10,height=10,units="cm",res=72)
#pdf("airs_hail_days_stddev.pdf",width=10,height=10)
png("airs_hail_days_stddev.png",width=10,height=10,units="in",res=300)
#par(mar=c(5.5,3,2,2)+0.1)
#par(mar=c(5.5,3,3,2)+0.1)
par(mar=c(2,0,3,2)+0.1)
par(oma=c(0.1,0.1,0.1,0.1))
par(cex=1.75)
image.plot(lon,rev(lat),sd.nhd2,
           #main="Std. dev hail days per year",xlab="",ylab="",
           xlab="",ylab="",
           xaxt="n",yaxt="n",
           #col=colorRampPalette(rev(brewer.pal(9,"Spectral")))(64),
           #col=oce.colorsViridis(64),
           #col=oce.colorsViridis(20),
           #col=tim.colors(20),
           col=tim.colors(12),
           #legend.lab="Std. dev. hail days/year",
           #horizontal=TRUE,cex.axis=1.25,cex.lab=1.25)
           horizontal=TRUE,
           breaks=seq(0,1.2,by=0.1))
           #legend.args=list(text="Std. dev. hail days/year",cex=1.25))
title("Std. dev. hail days per year",line=0.5,cex.main=2.0)
text(155,-12.5,"b",cex=3.0)
oz(add=TRUE)
map.cities(minpop=1.e06,cex=2.0)
box(lwd=3)
dev.off()



#Need to check std dev by evaluating on ship/gpats int days first.
#sd.sm2 <- apply(simplify2array(smt_by_year),c(1,2),sd,na.rm=T)

