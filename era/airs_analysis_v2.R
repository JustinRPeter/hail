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
max.date <- "2012-12-31"
nyears <- as.double(difftime(max.date,min.date)/365.25)


ex.rtime_era <- which(rdate_era >= min.date &
                      rdate_era <= max.date)
ex.rtime_gpats <- which(rdate_gpats >= min.date &
                        rdate_gpats <= max.date)
ex.rtime_hail <- which(hail_data$rdate_hail >= min.date &
                       hail_data$rdate_hail <= max.date)

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


##nship <- ship %*% delta.ship
#nship <- ship * delta.ship
#ship <- nship
#ship.lt0 <- which(ship < 0 )
#ship[ship.lt0] <- 0

#Now extract the grid points corresponding to the major cities
#Brisbane, Sydney, Melbourne, Perth, Adelaide
#Next locations are in the above order

grid_lon <- c(153.0251,151.2093,144.9631,115.8605,138.6007)
grid_lat <- c(-27.4698,-33.8688,-37.8136,-31.9505,-34.9285)
#grid_lon <- c(153.0251,144.9631,115.8605,138.6007)
#grid_lat <- c(-27.4698,-37.8136,-31.9505,-34.9285)
ngrids <- length(grid_lon)

#grid_min_lon<-grid_max_lon<-grid_min_lat<-grid_max_lat<-c()

#for (i in 1:ngrids){
#    grid_min_lon[i] <- lon[max(which(lon <= grid_lon[i]))]
#    grid_max_lon[i] <- lon[min(which(lon >= grid_lon[i]))]
#    grid_min_lat[i] <- lat[min(which(lat <= grid_lat[i]))]
#    grid_max_lat[i] <- lat[max(which(lat >= grid_lat[i]))]
#}


#Find the indices of the ERA grid boxes in which the 
#specified grid_lons/lats are located
ll_grid <- expand.grid(ll_lon=lon,ll_lat=lat)
locs <- data.frame(cbind(x=grid_lon,y=grid_lat))
closest <- nn2(ll_grid[,1:2], locs,1)
ex.lon <- ll_grid[closest$nn.idx,]$ll_lon
ex.lat <- ll_grid[closest$nn.idx,]$ll_lat
ex.lon.inds <- vector("numeric",length=ngrids)
ex.lat.inds <- vector("numeric",length=ngrids)
#Extract the indices of the lat and lon
for (i in 1:ngrids){
    ex.lon.inds[i] <- which(lon == ex.lon[i])
    ex.lat.inds[i] <- which(lat == ex.lat[i])
}

delta_era=0.75/2

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
}

#Find the intersection of dates in BoM/ERA/GPATS.
days_hail_in_era   <- vector("list",length=ngrids)
days_hail_in_gpats <- vector("list",length=ngrids)
for (i in 1:ngrids){
    days_hail_in_era[[i]]   <- which(rdate_era %in% days_grid_hail[[i]])
    days_hail_in_gpats[[i]] <- which(rdate_gpats %in% days_grid_hail[[i]])
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
    if(length(ex.gpats_data[[i]] >0)){ #loop bombs without this condition.
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

for (i in 1:ngrids){
    ndays.int.ship.gpats[[i]] <- length(which(max.ship_data[[i]] > 0.5 & 
                                              max.gpats_data[[i]] > 0))
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
    ndays.hail[i] <- length(days_hail_in_era[[i]])/4
    #ndays.hail[i] <- length(days_hail_in_era[[i]])
}

#for (i in 1:ngrids){
#   max.ship_data[[i]]<-max.ship_data[[i]][which(max.ship_data[[i]] > 0.5)]
#   max.gpats_data[[i]]<-max.gpats_data[[i]][which(max.gpats_data[[i]] > 0)]
#}

#Now evaluate the no. of hail days as a function of SHIP>0.5 and lightning
model <- lm(ndays.hail ~ ndays.int.ship.gpats-1)
plot(ndays.int.ship.gpats,ndays.hail)
#model <- lm(ndays.hail[c(1,3,4,5)] ~ ndays.int.ship.gpats[c(1,3,4,5)]-1)
#plot(ndays.int.ship.gpats[c(1,3,4,5)],ndays.hail[c(1,3,4,5)])
abline(a=0,b=model$coefficients[["ndays.int.ship.gpats"]])

#Now use only lightning data from 2008-2014
#co.dates <- which(rdate_gpats <= max(rdate_era))
#rdate_gpats <- rdate_gpats[co.dates]
#total_lightning <- total_lightning[,,co.dates]

#Find the max ship and lightning on a day
max.ship <- apply(ship,c(1,2),roll_max,4,by=4,align="left",fill=NA)
max.ship <- max.ship[seq(1,dim(max.ship)[1],4),,]
max.gpats <- apply(total_lightning,c(1,2),roll_max,24,by=24,align="left",fill=NA)
max.gpats <- max.gpats[seq(1,dim(max.gpats)[1],24),,]

#Determine which SHIP obs are greater than 0.5
#and set them to 1. Less than 0.5=NA
sh <- which(max.ship >= 0.5)
shn <- which(max.ship < 0.5)

max.ship[sh] <- 1
max.ship[shn] <- NA

#Determine which grids have lightning
#and set equal to 1.
gp <- which(max.gpats >= 1)
gpn <- which(max.gpats < 1)

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

#Now apply the regression equation
png("airs_hail_days.png")
#nhd <- model$coefficients[1]*smt/(7*4)
nhd <- model$coefficients[1]*smt/(nyears)
image.plot(lon,rev(lat),nhd,
           #col=colorRampPalette(rev(brewer.pal(9,"Spectral")))(64),
           #col=oce.colorsViridis(64),
           col=oce.colorsViridis(20),
           legend.lab="Hail days/year",
           horizontal=TRUE)
oz(add=TRUE)
box()
dev.off()

#Evaluate the standard deviation
