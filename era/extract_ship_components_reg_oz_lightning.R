#Extract ship components for regression
#so that we can construct a regression for lightning.
#This is similar to "extract_ship_components_reg_data.R"
#However, here we just extract lightning probabilities
#and rather than only using ERA boxes which contain hail reports
#we use a random selection over all land areas.

#Integrate the ship, gpats, and hail reports
#(1) Find the max. ship and lightning on each date and in
#    each ERA grid box where hail was observed (via BoM storm reports)
#(2) Write these times and values to a file and include
#    a flag indicating presence of hail.

#Justin R. Peter
#ICACS, Uni. of Southern Qld
#18 August 2016

rm(list=ls())
gc()


source("/home/u1066578/jpeter/rfiles/lib/thermo/thermo.r")

library(ncdf4)
library(fields)
library(RcppRoll)
library(oce) #Some good color scales
library(RColorBrewer)
library(oz)
library(maps)
library(mapdata)

#library(TTR) #For time series analysis

#The data sets
#ERA, GPATS and BoM hail reports
#fnm_era <- '/home/jpeter/justinp/rfiles/suncorp/data/eraint/sharpy_profiles/00_18/2008_2014/era_profiles_2008_2014_00_18.nc'
fnm_era <- '/home/u1066578/data/era_output/00_18/2008_2014/analysis/era_profiles_2008_2014_00_18.nc'

#Will use the 2-hourly averaged file for now, but will 
#need to use 24 hour eventually
#fnm_gpats <- '/home/jpeter/justinp/rfiles/suncorp/data/gpats/ts_lightning_grid/2008_2014/analysis/lightning_eragrid_eratime_Australia_2008_2014_00_18.nc'
#fnm_gpats <- '/home/jpeter/justinp/rfiles/suncorp/data/gpats/ts_lightning_grid/2008_2015/combined/lightning_eragrid_eratime_Australia_2008_2015_offset.nc'
fnm_gpats <- '/home/u1066578/data/ts_lightning_grid/2008_2015/analysis/lightning_eragrid_eratime_Australia_2008_2015_offset.nc'

#fnm_bom <- '/home/jpeter/justinp/rfiles/suncorp/data/bom_hail/reports.txt'

##Read the hail reports data
#hail_data <- read.table(fnm_bom)
#hail_data[[1]] <- as.Date(hail_data[[1]])
##hail_data[[1]] <- as.POSIXct(hail_data[[1]],tz="UTC")
#names(hail_data) <- c("rdate_hail", "lon","lat")

#Read the mask file
#fnm_lsmask <- '/home/jpeter/justinp/rfiles/suncorp/data/eraint/netcdf_download/lsmask/era_interim_lsmask.nc'
fnm_lsmask <- '/home/u1066578/data/eraint/lsmask/era_interim_lsmask.nc'

#Open the netcdf files
nc_era   <- nc_open(fnm_era)
nc_gpats <- nc_open(fnm_gpats)
nc_lsmask <- nc_open(fnm_lsmask)

print(paste("SHARPPy profile file has",nc_era$nvars,"variables"))
print(paste("GPATS lightning file has",nc_gpats$nvars,"variables"))

#Get the dimensions
nlat_era    <- nc_era$dim$lat$len
nlon_era    <- nc_era$dim$lon$len
ntime_era   <- nc_era$dim$time$len

nlat_gpats  <- nc_gpats$dim$lat$len
nlon_gpats  <- nc_gpats$dim$lon$len
ntime_gpats <- nc_gpats$dim$time$len

lat_era  <- nc_era$dim$lat$vals
lon_era  <- nc_era$dim$lon$vals
time_era <- nc_era$dim$time$vals

lat_gpats  <- nc_gpats$dim$lat$vals
lon_gpats  <- nc_gpats$dim$lon$vals
time_gpats <- nc_gpats$dim$time$vals

#Convert time to a date
#Have to convert from hours to seconds (*60*60)
rtime_era   <- as.POSIXct(time_era*60*60, origin="1900-01-01", tz="GMT")   # in UTC
rtime_gpats <- as.POSIXct(time_gpats*60*60, origin="1900-01-01", tz="GMT") # in UTC
rdate_era   <- as.Date(rtime_era,tz="UTC")
rdate_gpats <- as.Date(rtime_gpats,tz="UTC")

#Get the variables
#Get the SHIP variable
v1<-nc_era$var[["ship"]]
#ship.raw <- ncvar_get(nc_era,v1)
ship <- ncvar_get(nc_era,v1)
#nc_close(nc_era)

#Read the SHIP component variables

v2 <- nc_era$var[["bplus"]]
mucape <- ncvar_get(nc_era,v2)

v3 <- nc_era$var[["lr75"]]
lr75 <- ncvar_get(nc_era,v3)

v4 <- nc_era$var[["h5_temp"]]
h5_temp <- ncvar_get(nc_era,v4)

v5 <- nc_era$var[["sfc_6km_shear"]]
s06 <- ncvar_get(nc_era,v5)

v6 <- nc_era$var[["pres"]]
pres <- ncvar_get(nc_era,v6)

v7 <- nc_era$var[["dwpc"]]
dwpc <- ncvar_get(nc_era,v7)

v10 <- nc_era$var[["cap"]]
cap <- ncvar_get(nc_era,v10)

nc_close(nc_era)


#Evaluate the parcel mixing ratio
pmr <- r_mix_td(pres*100.,dwpc+273.15)

#Get the total lightning variable
v8<-nc_gpats$var[["total_lightning"]]
total_lightning <- ncvar_get(nc_gpats,v8)
nc_close(nc_gpats)

#Get the land-sea mask
v9 <- nc_lsmask$var[["lsm"]]
lsmask <- ncvar_get(nc_lsmask,v9)
lat_lsmask <- nc_lsmask$dim$lat$vals
lon_lsmask <- nc_lsmask$dim$lon$vals


print("Finished reading data")

#Truncate everything to the dimensions of the smallest file (GPATS)
min_lat <- min(lat_gpats)
#max_lat <- max(lat_gpats)
max_lat <- -10.5
min_lon <- min(lon_gpats)
#max_lon <- max(lon_gpats)
max_lon <- 155.0
#
ex.lat_era <- which(lat_era >= min_lat & lat_era <= max_lat)
ex.lon_era <- which(lon_era >= min_lon & lon_era <= max_lon)
lat_era <- lat_era[ex.lat_era]
lon_era <- lon_era[ex.lon_era]

ship    <- ship[ex.lon_era,ex.lat_era,]
mucape  <- mucape[ex.lon_era,ex.lat_era,]
lr75    <- lr75[ex.lon_era,ex.lat_era,]
h5_temp <- h5_temp[ex.lon_era,ex.lat_era,]
s06     <- s06[ex.lon_era,ex.lat_era,]
pres    <- pres[ex.lon_era,ex.lat_era,]
dwpc    <- dwpc[ex.lon_era,ex.lat_era,]
pmr     <- pmr[ex.lon_era,ex.lat_era,]
cap     <- cap[ex.lon_era,ex.lat_era,]
#cls.pred   <- cls.pred[ex.lon_era,ex.lat_era,]
#

ex.lat_gpats <- which(lat_gpats >= min_lat & lat_gpats <= max_lat)
ex.lon_gpats <- which(lon_gpats >= min_lon & lon_gpats <= max_lon)
lat_gpats <- lat_gpats[ex.lat_gpats]
lon_gpats <- lon_gpats[ex.lon_gpats]
total_lightning <- total_lightning[ex.lon_gpats,ex.lat_gpats,]

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

#Also concatenate to the same time
date.min <- as.Date("2008-01-01")
date.max <- as.Date("2014-12-31") 

range.dates_gpats <- which(rdate_gpats >= date.min & rdate_gpats <= date.max)
total_lightning <- total_lightning[,,range.dates_gpats]

#Put the total lightning on the same time grid as the ERA-Interim
#In the next line we have to use dim(total_lightning)
#rather than length(range.dates_gpats) as the
#end points are lost in the previous step.
total_lightning <- apply(total_lightning,c(1,2),roll_max,3,by=6,align="center",na.rm=T)
total_lightning <- total_lightning[seq(1,dim(total_lightning)[1],6),,]



#In this version we randomly draw several grid boxes from Australia
#Only extract those lats and lons which occur over land
#Set the % of samples to use for the training data set
#train.fraction <- 1.0
train.fraction <- 0.9
#train.fraction <- 0.7
#train.fraction <- 0.5
#train.fraction <- 0.1
#train.fraction <- 0.05
ex.all.land <- as.data.frame(which(lsmask==1,arr.ind=T))
len.ex.all.land <- nrow(ex.all.land)
ex.era_grid.inds <- ex.all.land[sample(nrow(ex.all.land),train.fraction*len.ex.all.land),]

ex.era_grid.lons <- lon_era[ex.era_grid.inds$row]
ex.era_grid.lats <- lat_era[ex.era_grid.inds$col]

ex.era_grid <- data.frame(cbind(ex.era_grid.lons,ex.era_grid.lats))
names(ex.era_grid) <- c("lon","lat")

#A plot to see which points are sampled
#png("extracted_samples.png",width=10,height=10,units="in",res=300)
postscript("extracted_samples.eps",onefile=FALSE,
           width=10,height=10,paper="special",horizontal=FALSE)
image(lon,rev(lat),t(apply(lsmask,1,rev)),
      xlab="",ylab="",
      #col=tim.colors(2))
      col=c("turquoise","red4"))
points(ex.era_grid$lon,ex.era_grid$lat,col="gold",lwd=3)
oz(add=TRUE,lwd=3)
#map.cities(minpop=1.e06,cex=2.0,font=2)
map.cities(minpop=1.e06,cex=2.0,lwd=3,label=FALSE)
box(lwd=3.0)
dev.off()


#Also extract just some random times
time.sample <- sort(sample.int(dim(ship)[3],train.fraction*dim(ship)[3]))
#time.sample <- sort(sample.int(dim(ship)[3],0.5*dim(ship)[3]))
ship    <- ship[,,time.sample]
mucape  <- mucape[,,time.sample]
lr75    <- lr75[,,time.sample]
h5_temp <- h5_temp[,,time.sample]
s06     <- s06[,,time.sample]
pres    <- pres[,,time.sample]
dwpc    <- dwpc[,,time.sample]
pmr     <- pmr[,,time.sample]
cap     <- cap[,,time.sample]
total_lightning <- total_lightning[time.sample,,]
total_lightning <- aperm(total_lightning,c(2,3,1))

#lightning.flag <- which(total_lightning > 0)

#Put all the convective variables in a list
scomps <-list(mucape=mucape,lr75=lr75,
            h5_temp=h5_temp,s06=s06,
            #pres=pres,dwpc=dwpc,
            pmr=pmr,
            cap=cap)


#Write to a file

sink(paste0("lightning_predictor_",train.fraction/2.0,".txt"),type='output',append=FALSE)
#sink("hail_predictor_random_oz.txt",type='output',append=FALSE)
#Put the variable names at the header of the file
cat("Lon","Lat","Date","Time",
    "SHIP","MUCAPE","LR75","H5_TEMP","S06","PMR","CAP",
    "GPATS","\n")
    #"Lightning_Occurred", "\n")

for (ill in 1:length(ex.era_grid$lon)){  #For each lat/lon grid box
#for (ill in 2){ 
  #cat("Iteration = ,", ill, "\n")
#for (ill in 1:10){ 

  #ex.era_grid <- data.frame(ex.lon,ex.lat)
  ex.lon <- ex.era_grid$lon[ill]
  ex.lat <- ex.era_grid$lat[ill]

  #cat(ex.lon,ex.lat)
  
  ex.ship.lon <- which(lon_era == ex.lon)
  ex.ship.lat <- which(lat_era == ex.lat)
  ex.ship <- ship[ex.ship.lon,ex.ship.lat,]

  ex.scomps <- vector("list",length=length(scomps))
  for (i in seq_along(scomps)){
      ex.scomps[[i]] <- scomps[[i]][ex.ship.lon,ex.ship.lat,]
  }
  names(ex.scomps) <- paste0("ex.",names(scomps))

  #Extract the lightning data in the grid box from the GPATS data
  ex.gpats.lon       <- which(lon_gpats == ex.lon)
  ex.gpats.lat       <- which(lat_gpats == ex.lat)
  ex.total_lightning <- total_lightning[ex.gpats.lon,ex.gpats.lat,]



  #Now write all the values to a file
  for (i in 1:length(time.sample)){
        cat(ex.lon,ex.lat,
            #rdate_era[time.sample[i]],
            strftime(rtime_era[time.sample[i]],
                     format="%Y-%m-%d %H:%M:%S",
                     tz="UTC"),
            ex.ship[i],
            sapply(ex.scomps,'[[',i),
            ex.total_lightning[i],"\n")
    }
}

  
sink()

#Convert text file to rds for quick reading
raw.data <- read.table(paste0("lightning_predictor_",train.fraction/2.0,".txt"),
                   header=TRUE,
                   na.strings=c("", "NA"),
                   stringsAsFactors=FALSE)
saveRDS(raw.data,paste0("lightning_predictor_",train.fraction/2.0,".rds"))

##Quick plot
#reg_gpats_data <- read.table("lightning_predictor.txt",
#                            header=TRUE,stringsAsFactors=FALSE)
#names(reg_hail_data) <- c("Lon","Lat","Date_ERA","Time_ERA","Date_GPATS","Time_GPATS",
##                          "Max_SHIP","SHIP_Max_GPATS","Max_GPATS","Hail_Occurred")
##attach(reg_hail_data)
#plot(reg_hail_data$Max_SHIP,reg_hail_data$Max_GPATS)
#points(reg_hail_data$SHIP_Max_GPATS,reg_hail_data$Max_GPATS,col="red")


#End of script
