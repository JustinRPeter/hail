#Script to analyse timeseries of the ERA_Interim data in the grid boxes
#around SE Qld.

#Justin R. Peter
#ICACS, Uni. of Southern Qld
#27 July 2016

rm(list=ls())

library(ncdf4)
#library(fields)
library(oce) #Some good color scales
library(RColorBrewer)
library(maps)
library(mapdata)
library(fields)

#fname <- "/home/jpeter/justinp/rfiles/suncorp/data/gpats/ts_lightning_grid/analysis/lightning_eragrid_eratime_Australia_2008_2015_00_18_mro.nc"
fname <- "/home/jpeter/justinp/rfiles/suncorp/data/gpats/ts_lightning_grid/2008_2015/analysis/lightning_eragrid_eratime_Australia_2008_2015_00_18_mro_seq_select.nc"
#fname <- "/home/jpeter/justinp/rfiles/suncorp/data/gpats/ts_lightning_grid/2008_2014/analysis/lightning_eragrid_eratime_Australia_2008_2015_00_18_mro_seq_select.nc"

#library(hdf5)
#library(CircStats) # Could also use circular

#library(maps)
#library(mapdata)
#library(mapproj)


#OPen the netcdf file
nc<-nc_open(fname)

print(paste("The file has",nc$nvars,"variables"))

lat <-
#Get the dimensions
nlat   <- nc$dim$lat$len
nlon   <- nc$dim$lon$len
ntime  <- nc$dim$time$len

lat <- nc$dim$lat$vals
lon <- nc$dim$lon$vals
time <- nc$dim$time$vals

#We need to sort the time variable
time.sort.int <- sort.int(time,index.return=TRUE)

#Get the total lightning variable
v1<-nc$var[["total_lightning"]]

total_lightning <- ncvar_get(nc,v1)

#Now sort the lightning and time
time.sort            <- time[time.sort.int$ix]
total_lightning.sort <- total_lightning[,,time.sort.int$ix]

#We need to transform this matrix
tot_lightning <- array(NA, dim=(c(nlon,nlat,ntime)))

for (i in 1:ntime){
    tot_lightning[,,i] <- t(apply(total_lightning.sort[,,i],1,rev))
}

#tot_lightning <- t(apply(tot_lightning,1,rev))

#Convert time to a date
#Have to convert from hours to seconds (*60*60)
rtime <- as.POSIXct(time*60*60, origin = "1900-01-01", tz = "GMT")  # in UTC
#Exatrct only the "06" data
ex.06.times <- which(substr(rtime,12,13) == "06")


#Average over the grid points in the domain
#to get a time series
tot_lightning.mean <- apply(tot_lightning,3,mean,na.rm=TRUE)

#Average over all times to get a snapshot of lightning
tot_lightning.grid.mean <- apply(tot_lightning,c(1,2),mean,na.rm=TRUE)
#tot_lightning.grid.mean <- t(apply(tot_lightning.grid.mean,1,rev))

#Extract all the ship components corresponding to the 2014-11-27 hail storm
#storm <- which(as.Date(rtime) == "2014-11-27")
#image.plot(lon,rev(lat),ship[,,storm[2]])

#-----Open the insurance data and find all the unique dates

insdir    <- '/home/jpeter/justinp/rfiles/suncorp/data/insurance2/clean/'
insfile   <- 'HOME_all_claims_complete.rds'
#insfile   <- 'HOME_all_claims_complete.csv'
#ins_data <- read.csv(paste0(insdir,insfile),na.strings=c("NA",""))
ins_data <- readRDS(paste0(insdir,insfile))

ins_data <- ins_data[!is.na(ins_data$Postcode),]
ins_data <- ins_data[!is.na(ins_data$Address),]
#ins_data <- ins_data[ins_data$Longitude>=min(lon)&ins_data$Longitude<=max(lon),]
#ins_data <- ins_data[ins_data$Latitude>=min(lat) & ins_data$Latitude<=max(lat),]
ins_data <- ins_data[ins_data$Longitude>=151.875&ins_data$Longitude<=152.265,]
ins_data <- ins_data[ins_data$Latitude>=-28.125 & ins_data$Latitude<=-27.375,]
ins_data <- ins_data[ins_data$Incurred.To.Date > 10,]

ins_data$ins.date <- as.Date(ins_data$Date.Occurred,"%d/%m/%Y")

#Find the unique dates
uniq.ins.dates <- unique(ins_data$ins.date)

#Find the times in the ship data which correspond to damage days
#ex.dates <- which(as.Date(rtime[ex.06.times]) %in% uniq.ins.dates)
ex.dates <- which(as.Date(rtime) %in% uniq.ins.dates)

#Extract pertinent times from the Lightning
#tot_lightning.claim.days <- tot_lightning[,,ex.06.times[ex.dates]]
tot_lightning.claim.days <- tot_lightning[,,ex.dates]

mean.tot_lightning.claim.days <- apply(tot_lightning.claim.days,c(1,2),mean,na.rm=TRUE)

#Only use the "06" data
#ex.06.times <- which(substr(rtime[exdates],12,13) == "06")

#Plot the insurance claims on the extracted region
minx=min(lon)
maxx=max(lon)
miny=min(lat)
maxy=max(lon)

#png("mean_ship_on_claim_days.png",
pdf("mean_lightning_on_claim_days.pdf",width=11,height=10,paper="special")

#Could also use filled.contour
#filled.contour(lon,rev(lat),tot_lightning.grid.mean,plot.axes={axis(1);axis(2);map(xlim=c(minx,maxx),ylim=c(miny,maxy),add=T,lwd=2)},color.palette=oce.colorsCDOM,nlevels=15)


#image.plot(lon,rev(lat),ship.claim.days[,,1])
image.plot(lon,rev(lat),mean.tot_lightning.claim.days[,],
           zlim=c(0,max(mean.tot_lightning.claim.days)),
           xlab="Longitude",ylab="Latitude",
           legend.lab="Flash count/hour",
           #col=oce.colorsViridis(64),
           #horizontal=TRUE,
           col=rev(brewer.pal(64,"Spectral")),
           cex.axis=1.25,cex.lab=1.25)
           #col=topo.colors(64))
           #col=terrain.colors(64))
map('worldHires',interior=T,xlim=c(minx,maxx),ylim=c(miny,maxy),lwd=2.5,add=T)
#map.cities(label=T,minpop=5.e04,col="darkblue",cex=1.35,lwd=2)
points(ins_data$Longitude,ins_data$Latitude)
map.cities(label=T,minpop=5.e04,col="darkblue",cex=1.35,lwd=2)
box()
dev.off()

pdf("mean_bris_ship_tseries_withclaim_days.pdf",width=10,height=10,paper="special")
#Now plot a time series of SHIP and overlay hail claim days.
date.min = min(as.Date(rtime[ex.06.times][ex.dates]))
date.max = max(as.Date(rtime[ex.06.times][ex.dates]))
#See http://stackoverflow.com/questions/4355042/label-x-axis-in-time-series-plot-using-r
x.Date <- as.Date(paste(2009:2014,01, 01, sep = "-"))
#plot(as.Date(rtime[ex.06.times]),ship[3,2,ex.06.times],
#Use the mean of the grid boxes over Brisbane and the grid box north
ship.avg.bris <- apply(ship[3,c(2,3),ex.06.times],2,mean)
#plot(as.Date(rtime[ex.06.times]),ship[3,2,ex.06.times],
plot(as.Date(rtime[ex.06.times]),ship.avg.bris,
     xlim=c(date.min,date.max),
     xaxt="n",type="l",
     xlab="Date",ylab="SHIP")
axis(1,at=x.Date,label=paste(2009:2014),las=0)
#points(as.Date(rtime[ex.06.times][ex.dates]),ship.claim.days[3,3,],col="red",lwd=3)
points(as.Date(rtime[ex.06.times][ex.dates]),ship.claim.days[3,3,],col="red",lwd=3)
dev.off()







