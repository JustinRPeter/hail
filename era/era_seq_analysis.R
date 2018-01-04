#Script to analyse timeseries of the ERA_Interim data in the grid boxes
#around SE Qld.

#Justin R. Peter
#ICACS, Uni. of Southern Qld
#27 July 2016

rm(list=ls())

library(ncdf4)
library(fields)
library(oce) #Some good color scales
library(maps)
library(mapdata)

library(TTR) #For time series analysis

fname <- "/home/jpeter/justinp/rfiles/suncorp/data/eraint/sharpy_profiles/00_18/analysis/era_profiles_1979_2014_00_18_seq_select.nc"

fnamel <- "/home/jpeter/justinp/rfiles/suncorp/data/gpats/ts_lightning_grid/analysis/lightning_eragrid_eratime_Australia_2008_2015_00_18_mro_seq_select.nc"


#OPen the netcdf file
nc<-nc_open(fname)
ncl <- nc_open(fnamel)

print(paste("SHARPPy profile file has",nc$nvars,"variables"))
print(paste("GPATS lightning file has",ncl$nvars,"variables"))

#Get the dimensions
nlat   <- nc$dim$lat$len
nlon   <- nc$dim$lon$len
ntime  <- nc$dim$time$len
ntimel <- ncl$dim$time$len

lat <- nc$dim$lat$vals
lon <- nc$dim$lon$vals
time <- nc$dim$time$vals
timel <- ncl$dim$time$vals

#Convert time to a date
#Have to convert from hours to seconds (*60*60)
rtime  <- as.POSIXct(time*60*60, origin = "1900-01-01", tz = "GMT")  # in UTC
rtimel <- as.POSIXct(timel*60*60, origin = "1900-01-01", tz = "GMT")  # in UTC

#Get the variables
#Get the SHIP variable
v1<-nc$var[["ship"]]
ship.raw <- ncvar_get(nc,v1)
nc_close(nc)

#Get the total lightning variable
v2<-ncl$var[["total_lightning"]]
total_lightning <- ncvar_get(ncl,v2)
nc_close(ncl)

#We need to transform the matrices for processing with R
#We need to transform this matrix
ship <- array(NA, dim=(c(nlon,nlat,ntime)))

for (i in 1:ntime){
    ship[,,i] <- t(apply(ship.raw[,,i],1,rev))
}

#Now sort the lightning and time
#Get the indices of the sorted times
timel.sort.int <- sort.int(timel,index.return=TRUE)
timel.sort            <- timel[timel.sort.int$ix]
total_lightning.sort <- total_lightning[,,timel.sort.int$ix]

#We need to apply an offset the GPATS times.
#This is due to the NCO averaging routine which applies the time at
#the end of the averaging period
timel.sort <- timel.sort-1.5
rtimel <- as.POSIXct(timel.sort*60*60, origin = "1900-01-01", tz = "GMT")  # in UTC
#We need to transform this matrix
tot_lightning <- array(NA, dim=(c(nlon,nlat,ntimel)))

for (i in 1:ntimel){
    tot_lightning[,,i] <- t(apply(total_lightning.sort[,,i],1,rev))
}

#Find the time when the ERA and GPATS correspond
#Find the indice of the ERA time corresponsding to the
#start of the GPAT file
#era.time.start.ind <- which(rtime == rtimel[1])
same.times.era.inds   <- which(rtime %in% rtimel)
same.times.gpats.inds <- which(rtimel %in% rtime)

ship.cotime          <- ship[,,same.times.era.inds]
tot_lightning.cotime <- tot_lightning[,,same.times.gpats.inds]
co.rtime             <- rtime[same.times.era.inds]

#Now just look at positive ship
#ship.gpats <- ship[which(ship.gpats > 1.0)]
#floop<-which(tot_lightning.era[3,3,]>50.0)
#plot(ship.gpats[3,3,floop],tot_lightning.era[3,3,floop])
#fit.ship.gpats <- lm(tot_lightning.era[3,3,floop] ~ ship.gpats[3,3,floop] 
#abline(lm(tot_lightning.era[3,3,floop] ~ ship.gpats[3,3,floop]))

#Average over the grid points in the domain
#ship.mean <- apply(ship,3,mean,na.rm=TRUE)
ship.cotime.mean <- apply(ship.cotime,3,mean,na.rm=TRUE)
tot_lightning.cotime.mean <- apply(tot_lightning.cotime,3,mean,na.rm=TRUE)
#Construct a running mean
run.mean.window <- 7
#run.mean.window <- 1
runmean.ship.cotime.mean <- SMA(ship.cotime.mean,n=run.mean.window)
runmean.tot_lightning.cotime.mean <- SMA(tot_lightning.cotime.mean,n=run.mean.window)
#Convert to a timeseries for analysis
ship.cotime.mean.ts <- ts(ship.cotime.mean,frequency=365.25*4,start=c(2008))
runmean.ship.cotime.mean.ts <- ts(runmean.ship.cotime.mean,frequency=365.25*4,start=c(2008))
tot_lightning.cotime.mean.ts <- ts(tot_lightning.cotime.mean,frequency=365.25*4,start=c(2008))
runmean.tot_lightning.cotime.mean.ts <- ts(runmean.tot_lightning.cotime.mean,frequency=365.25*4,start=c(2008))

#Plot each of the individual time series
plot(ship.cotime.mean.ts,
     xlab="Time",ylab="SHIP")
lines(runmean.ship.cotime.mean.ts,col="red")

plot(tot_lightning.cotime.mean.ts,
     xlab="Time",ylab="Lightning strikes")
lines(runmean.tot_lightning.cotime.mean.ts,col="red")

#Plot them on the same graph
pdf("ship_gpats_time_series_2008_2014.pdf",width=10,height=10,paper="special")
par(mar = c(5,5,2,5))
plot(runmean.ship.cotime.mean.ts,ylab="SHIP",type="l")
par(new = T)
plot(runmean.tot_lightning.cotime.mean.ts,
     col="red",type="l",lwd=3,
     axes=F, xlab=NA, ylab=NA)
axis(side = 4)
mtext(side = 4, line = 3, 'Lightning strikes')
dev.off()

#Extract only the "06" data
#ex.06.times <- which(substr(rtime,12,13) == "06")
ex.06.times <- which(substr(co.rtime,12,13) == "06")

#Only look at when we have a lightning strike
lightning.occurred.all <- which(tot_lightning.cotime.mean >= 10.0)
lightning.occurred.06  <- which(tot_lightning.cotime.mean[ex.06.times] >= 10.0)

#Original data. SHIP vs Lightning stikes
#wehn more than one lightning strike measured.
#Examine all times and 06 times.
plot(tot_lightning.cotime.mean[lightning.occurred.all],
     ship.cotime.mean[lightning.occurred.all])
plot(tot_lightning.cotime.mean[lightning.occurred.06],
     ship.cotime.mean[lightning.occurred.06])
#Same as above but looking at the rolling mean period chosen
pdf("Weekly_averaged_ship_lightning.pdf",width=10,height=10,paper="special")
plot(runmean.tot_lightning.cotime.mean[lightning.occurred.all],
     runmean.ship.cotime.mean[lightning.occurred.all],
     main="Weekly averaged ship vs lightning > 10 strikes/hour",
     xlab="Lightning strikes / hour",
     ylab="SHIP")
#plot(runmean.tot_lightning.cotime.mean[lightning.occurred.06],
#     runmean.ship.cotime.mean[lightning.occurred.06])

#Fit a linear model to the data
gpats.ship.model <- lm(runmean.ship.cotime.mean[lightning.occurred.all] ~ runmean.tot_lightning.cotime.mean[lightning.occurred.all] + 0)
#abline(gpats.ship.model$coefficients[1],gpats.ship.model$coefficients[2])
abline(0,gpats.ship.model$coefficients[1])
box()
dev.off()




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
ins_data <- ins_data[ins_data$Longitude>=152.625 &ins_data$Longitude<=153.375,]
ins_data <- ins_data[ins_data$Latitude>= -28.125 & ins_data$Latitude<= -27.375,]
ins_data <- ins_data[ins_data$Incurred.To.Date > 10,]

ins_data$ins.date <- as.Date(ins_data$Date.Occurred,"%d/%m/%Y")

#Find the unique dates
uniq.ins.dates <- unique(ins_data$ins.date)

#Find the times in the ship data which correspond to damage days
#ex.dates <- which(as.Date(rtime[ex.06.times]) %in% uniq.ins.dates)
#ex.dates <- which(as.Date(co.rtime[ex.06.times]) %in% uniq.ins.dates)
ex.dates <- which(as.Date(co.rtime) %in% uniq.ins.dates)

#Extract pertinent times from the SHIP
#Use the 06 ship sounding but lightning at any time!
#ship.claim.days <- ship[,,ex.06.times[ex.dates]]
#ship.claim.days <- ship.cotime[,,ex.06.times[ex.dates]]
#tot_lightning.claim.days <- tot_lightning.cotime[,,ex.06.times[ex.dates]] 
ship.claim.days <- ship.cotime[,,ex.dates]
tot_lightning.claim.days <- tot_lightning.cotime[,,ex.dates] 

mean.ship.claim.days <- apply(ship.claim.days,c(1,2),mean,na.rm=TRUE)
mean.tot_lightning.claim.days <- apply(tot_lightning.claim.days,c(1,2),mean,na.rm=TRUE)

#Only use the "06" data
#ex.06.times <- which(substr(rtime[exdates],12,13) == "06")

#Plot the insurance claims on the extracted region
minx=min(lon)
maxx=max(lon)
miny=min(lat)
maxy=max(lon)

#png("mean_ship_on_claim_days.png",
pdf("mean_ship_on_claim_days.pdf",width=11,height=10,paper="special")
#image.plot(lon,rev(lat),ship.claim.days[,,1])
image.plot(lon,rev(lat),mean.ship.claim.days[,],
           zlim=c(0,max(mean.ship.claim.days)),
           main="Avg. SHIP on claim days in Brisbane grid",
           xlab="Longitude",ylab="Latitude",
           legend.lab="SHIP",
           col=oce.colorsViridis(64),
           cex.axis=1.25,cex.lab=1.25)
           #col=topo.colors(64))
           #col=terrain.colors(64))
map('worldHires',interior=T,xlim=c(minx,maxx),ylim=c(miny,maxy),lwd=2.5,add=T)
#map.cities(label=T,minpop=5.e04,col="darkblue",cex=1.35,lwd=2)
points(ins_data$Longitude,ins_data$Latitude)
map.cities(label=T,minpop=5.e04,col="darkblue",cex=1.35,lwd=2)
box()
dev.off()

pdf("mean_tot_lightning_on_claim_days.pdf",width=11,height=10,paper="special")
#image.plot(lon,rev(lat),ship.claim.days[,,1])
image.plot(lon,rev(lat),mean.tot_lightning.claim.days[,]/2,
           zlim=c(0,max(mean.tot_lightning.claim.days)/2.),
           main="Avg. lightning strikes on claim days in Brisbane grid",
           xlab="Longitude",ylab="Latitude",
           legend.lab="Lightning strikes / hour",
           col=oce.colorsViridis(64),
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
#date.min = min(as.Date(rtime[ex.06.times][ex.dates]))
#date.max = max(as.Date(rtime[ex.06.times][ex.dates]))
date.min = min(as.Date(co.rtime[ex.dates]))
date.max = max(as.Date(co.rtime[ex.dates]))
#See http://stackoverflow.com/questions/4355042/label-x-axis-in-time-series-plot-using-r
#x.Date <- as.Date(paste(2009:2014,01, 01, sep = "-"))
x.Date <- as.Date(paste(2009:2015,01, 01, sep = "-"))
#plot(as.Date(rtime[ex.06.times]),ship[3,2,ex.06.times],
#Use the mean of the grid boxes over Brisbane and the grid box north
#ship.avg.bris <- apply(ship[3,c(2,3),ex.06.times],2,mean)
#ship.avg.bris <- ship[3,3,ex.06.times]
#ship.avg.bris <- ship.cotime[3,3,ex.06.times]
ship.avg.bris <- ship.cotime[3,3,]
#ship.avg.bris <- ship.cotime[3,3,ex.dates]
#plot(as.Date(rtime[ex.06.times]),ship[3,2,ex.06.times],
#plot(as.Date(rtime[ex.06.times]),ship.avg.bris,
#plot(as.Date(co.rtime[ex.06.times]),ship.avg.bris,
#Plot SHIP at any time on that day
plot(as.Date(co.rtime),ship.avg.bris,
     xlim=c(date.min,date.max),
     main="Time series of SHIP averaged over Brisbane grid",
     xaxt="n",type="l",
     xlab="Date",ylab="SHIP")
axis(1,at=x.Date,label=paste(2009:2015),las=0)
#points(as.Date(rtime[ex.06.times][ex.dates]),ship.claim.days[3,3,],col="red",lwd=3)
#points(as.Date(rtime[ex.06.times][ex.dates]),ship.claim.days[3,3,],col="red",lwd=3)
#points(as.Date(co.rtime[ex.06.times][ex.dates]),ship.claim.days[3,3,],col="red",lwd=3)
points(as.Date(co.rtime[ex.dates]),ship.claim.days[3,3,],col="red",lwd=3)
dev.off()

pdf("mean_bris_lightning_tseries_withclaim_days.pdf",width=10,height=10,paper="special")
#Now plot a time series of SHIP and overlay hail claim days.
#date.min = min(as.Date(rtime[ex.06.times][ex.dates]))
#date.max = max(as.Date(rtime[ex.06.times][ex.dates]))
date.min = min(as.Date(co.rtime[ex.dates]))
date.max = max(as.Date(co.rtime[ex.dates]))
#See http://stackoverflow.com/questions/4355042/label-x-axis-in-time-series-plot-using-r
#x.Date <- as.Date(paste(2009:2014,01, 01, sep = "-"))
x.Date <- as.Date(paste(2009:2015,01, 01, sep = "-"))
#plot(as.Date(rtime[ex.06.times]),ship[3,2,ex.06.times],
#Use the mean of the grid boxes over Brisbane and the grid box north
#ship.avg.bris <- apply(ship[3,c(2,3),ex.06.times],2,mean)
#ship.avg.bris <- ship[3,3,ex.06.times]
#tot_lightning.avg.bris <- tot_lightning.cotime[3,3,ex.06.times]
#tot_lightning.avg.bris <- tot_lightning.cotime[3,3,ex.dates]
tot_lightning.avg.bris <- tot_lightning.cotime[3,3,]
#plot(as.Date(rtime[ex.06.times]),ship[3,2,ex.06.times],
#plot(as.Date(rtime[ex.06.times]),ship.avg.bris,
#plot(as.Date(co.rtime[ex.06.times]),tot_lightning.avg.bris/2.,
#PLot lightning at any time
#plot(as.Date(co.rtime[ex.06.times]),tot_lightning.avg.bris/2.,
plot(as.Date(co.rtime),tot_lightning.avg.bris/2.,
     xlim=c(date.min,date.max),
     main="Time series of lightning strikes averaged over Brisbane grid",
     #log="y",
     #ylim=c(0,20),
     xaxt="n",type="l",
     xlab="Date",ylab="Lightning strikes / hour")
axis(1,at=x.Date,label=paste(2009:2015),las=0)
#points(as.Date(rtime[ex.06.times][ex.dates]),ship.claim.days[3,3,],col="red",lwd=3)
#points(as.Date(rtime[ex.06.times][ex.dates]),ship.claim.days[3,3,],col="red",lwd=3)
#points(as.Date(co.rtime[ex.06.times][ex.dates]),tot_lightning.claim.days[3,3,]/2.,col="red",lwd=3)
points(as.Date(co.rtime[ex.dates]),tot_lightning.claim.days[3,3,]/2.,col="red",lwd=3)
dev.off()

#Plot SHIP versus lightning on claim days
#plot(


pdf("histogram_ship_bne_2009_2014_claims.pdf",width=10,height=10,paper="special")
hist(ship.claim.days,freq=F,col="red",main="Histogram of SHIP over Brisbane grid")
hist(ship.avg.bris,freq=F,add=T,col="blue")
legend("topright",c("2009-2014","claim days"),lty=c(1,1),col=c("blue","red"))
dev.off()







