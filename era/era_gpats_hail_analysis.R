#Integrate the ship, gpats, and hail reports
#(1) Find the max. ship and lightning on each date and in
#    each ERA grid box where hail was observed (via BoM storm reports)
#(2) Write these times and values to a file and include
#    a flag indicating of 

#Justin R. Peter
#ICACS, Uni. of Southern Qld
#8 August 2016

rm(list=ls())

library(ncdf4)
library(fields)
library(oce) #Some good color scales
library(RColorBrewer)
library(maps)
library(mapdata)

library(TTR) #For time series analysis

#The data sets
#ERA, GPATS and BoM hail reports
fnm_era <- '/home/jpeter/justinp/rfiles/suncorp/data/eraint/sharpy_profiles/00_18/2008_2014/era_profiles_2008_2014_00_18.nc'

#Will use the 2-hourly averaged file for now, but will 
#need to use 24 hour eventually
#fnm_gpats <- '/home/jpeter/justinp/rfiles/suncorp/data/gpats/ts_lightning_grid/2008_2014/analysis/lightning_eragrid_eratime_Australia_2008_2014_00_18.nc'
fnm_gpats <- '/home/jpeter/justinp/rfiles/suncorp/data/gpats/ts_lightning_grid/2008_2015/combined/lightning_eragrid_eratime_Australia_2008_2015_offset.nc'

fnm_bom <- '/home/jpeter/justinp/rfiles/suncorp/data/bom_hail/reports.txt'

#Read the hail reports data
hail_data <- read.table(fnm_bom)
hail_data[[1]] <- as.Date(hail_data[[1]])
#hail_data[[1]] <- as.POSIXct(hail_data[[1]],tz="UTC")
names(hail_data) <- c("rdate_hail", "lon","lat")

#Open the netcdf files
nc_era   <- nc_open(fnm_era)
nc_gpats <- nc_open(fnm_gpats)


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
nc_close(nc_era)

#Get the total lightning variable
v2<-nc_gpats$var[["total_lightning"]]
total_lightning <- ncvar_get(nc_gpats,v2)
nc_close(nc_gpats)

#We need to transform the matrices for processing with R
#We need to transform this matrix
#ship <- array(NA, dim=(c(nlon_era,nlat_era,ntime_era)))

#for (i in 1:ntime_era){
#    ship[,,i] <- t(apply(ship.raw[,,i],1,rev))
#}

#Use only the hail reports within specified ERA grid box
#Just use Brisbane grid for moment
#ex.lon <- 153.0
#ex.lat <- -27.75
delta_era <- 0.75/2
ex.era_grid <- read.table('/home/jpeter/justinp/rfiles/suncorp/data/bom_hail/lonlat.txt')
names(ex.era_grid) <- c('lon','lat')

#Write to a file
sink("ship_lightning_hail_occurrence.txt",type='output',append=FALSE)
#Put the variable names at the header of the file
cat("Lon","Lat","Date_ERA","Time_ERA","Date_GPATS","Time_GPATS",
    "Max_SHIP","SHIP_Max_GPATS","Max_GPATS","Hail_Occurred", "\n")

for (ill in 1:length(ex.era_grid$lon)){ 
#for (ill in 2){ 
  #cat("Iteration = ,", ill, "\n")
#for (ill in 10){ 

  #ex.era_grid <- data.frame(ex.lon,ex.lat)
  ex.lon <- ex.era_grid$lon[ill]
  ex.lat <- ex.era_grid$lat[ill]

  #cat(ex.lon,ex.lat)
  
  #ex.hail_data.inds <- which(hail_data$lon >= ex.era_grid$ex.lon-delta_era &
  #                           hail_data$lon <= ex.era_grid$ex.lon+delta_era &
  #                           hail_data$lat >= ex.era_grid$ex.lat-delta_era &
  #                           hail_data$lat <= ex.era_grid$ex.lat+delta_era)
  ex.hail_data.inds <- which(hail_data$lon >= ex.lon-delta_era &
                             hail_data$lon <= ex.lon+delta_era &
                             hail_data$lat >= ex.lat-delta_era &
                             hail_data$lat <= ex.lat+delta_era)

  ex.hail_data <- hail_data[ex.hail_data.inds,]

  ex.ship.lon <- which(lon_era == ex.lon)
  ex.ship.lat <- which(lat_era == ex.lat)
  ex.ship <- ship[ex.ship.lon,ex.ship.lat,]
  
  ex.gpats.lon       <- which(lon_gpats == ex.lon)
  ex.gpats.lat       <- which(lat_gpats == ex.lat)
  ex.total_lightning <- total_lightning[ex.gpats.lon,ex.gpats.lat,]
  
  #Now only use times where we are confident
  #of the hail reports
  #date.min <- as.Date("1990-01-01")
  date.min <- as.Date("2008-01-01") #This is when we have lightning data from
  date.max <- as.Date("2012-12-31") #There are missing reports from 2013/2014.
  
  range.dates_era   <- which(rdate_era   >= date.min & rdate_era   <= date.max)
  range.dates_gpats <- which(rdate_gpats >= date.min & rdate_gpats <= date.max)
  range.dates_hail  <- which(ex.hail_data$rdate_hail >= date.min &
                          ex.hail_data$rdate_hail <= date.max)
  #The above coincide will coincide if ERA and GPATS have been
  #averaged to the same grid.
  
  #Now apply the indices to the data
  ex.ship <- ex.ship[range.dates_era]
  ex.total_lightning <- ex.total_lightning[range.dates_gpats]
  ex.hail_data <- ex.hail_data[range.dates_hail,]
  
  #Now find the dates when hail was observed
  hail_dates <- unique(ex.hail_data$rdate_hail)

  #Only proceeed if a hail event was observed in the grid box
  if (length(hail_dates > 0)){

    #Get the days where hail was observed
    ex.dates.hail.ship            <- which(rdate_era[range.dates_era] %in% hail_dates)
    ex.dates.hail.total_lightning <- which(rdate_gpats[range.dates_gpats] %in% hail_dates)
    #The above are same so just do once. They won't be with full GPATS data!
    #ex.dates.hail <- which(rdate_era[range.dates] %in% hail_dates)

    #Get the days when no hail was observed
    ex.dates.no.hail.ship            <- which(!(rdate_era[range.dates_era] %in% hail_dates))
    ex.dates.no.hail.total_lightning <- which(!(rdate_gpats[range.dates_gpats] %in% hail_dates))
    
    #Extract the SHIP and lightning data for each date
    #on which hail was observed in the ERA grid.
    #This requires extracting chunks of length=4 and 24 
    #for the ERA and GPATS data, respectively. These values are
    #set by the granularity of the data sets in a 24-hour period
    #(6-hourly and one-hourly)
    chunk.ntimes.ship  <- 4 #00/06/12/18 UTC
    #chunk.ntimes.gpats <- 4 #No. of lightning obs/day
    chunk.ntimes.gpats <- 24 #No. of lightning obs/day
                            #Will need to change to 24 for full lightning data

    #Each chunk corresponds to a date! Either hail or no hail.
    chunk.start.hail.era  <- seq(1,length(ex.dates.hail.ship),4)
    chunk.end.hail.era    <- seq(chunk.ntimes.ship,length(ex.dates.hail.ship),4)
    chunk.start.no.hail.era  <- seq(1,length(ex.dates.no.hail.ship),4)
    chunk.end.no.hail.era    <- seq(chunk.ntimes.ship,length(ex.dates.no.hail.ship),4)

    chunk.start.hail.gpats <- seq(1,length(ex.dates.hail.total_lightning),24)
    chunk.end.hail.gpats   <- seq(chunk.ntimes.gpats,length(ex.dates.hail.total_lightning),24)
    chunk.start.no.hail.gpats <- seq(1,length(ex.dates.no.hail.total_lightning),24)
    chunk.end.no.hail.gpats   <- seq(chunk.ntimes.gpats,length(ex.dates.no.hail.total_lightning),24)
   
    #The length of the chunks is the number of days that hail was
    #observed in each ERA grid box.
    len.chunks.hail.era   <- length(chunk.start.hail.era)
    len.chunks.hail.gpats <- length(chunk.start.hail.gpats)
    len.chunks.no.hail.era   <- length(chunk.start.no.hail.era)
    len.chunks.no.hail.gpats <- length(chunk.start.no.hail.gpats)
    #These should be the same length
    stopifnot((len.chunks.hail.era == len.chunks.hail.gpats) ||
               (len.chunks.no.hail.era == len.chunks.no.hail.gpats))

    #Create vectors for the max SHIP and Lightning on
    #Hail and No Hail days.
    max.hail.ship            <- vector("numeric",length=len.chunks.hail.era)
    max.hail.ship.time       <- vector("numeric",length=len.chunks.hail.era)
    max.hail.ship.time.ind   <- vector("numeric",length=len.chunks.hail.era)
    max.no.hail.ship            <- vector("numeric",length=len.chunks.no.hail.era)
    max.no.hail.ship.time       <- vector("numeric",length=len.chunks.no.hail.era)
    max.no.hail.ship.time.ind   <- vector("numeric",length=len.chunks.no.hail.era)

    max.hail.total_lightning    <- vector("numeric",length=len.chunks.hail.gpats)
    max.hail.total_lightning.time <- vector("numeric",length=len.chunks.hail.gpats)
    max.hail.total_lightning.time.ind <- vector("numeric",length=len.chunks.hail.gpats)
    max.no.hail.total_lightning    <- vector("numeric",length=len.chunks.no.hail.gpats)
    max.no.hail.total_lightning.time <- vector("numeric",length=len.chunks.no.hail.gpats)
    max.no.hail.total_lightning.time.ind<-vector("numeric",length=len.chunks.no.hail.gpats)

    #The value of the ERA SHIP closest in time to the maximum
    #total no. of lightning strikes (which we use as our proxy of hail).
    ship.at.max.lightning.hail <- vector("numeric",length=len.chunks.hail.era)
    ship.at.max.lightning.no.hail <- vector("numeric",length=len.chunks.no.hail.era)

#---#(1) Evaluate the parameters for all the HAIL days
    #Find the maximum SHIP value in for each hail event date
    for (i in 1:len.chunks.hail.era){
        #Extract the chunk of SHIP data corresponding to each hail event date
        chunk.hail.ship.time.inds <- ex.dates.hail.ship[
                                         chunk.start.hail.era[i]:
                                         chunk.end.hail.era[i]]
        #chunk.hail.ship <- ex.ship[ex.dates.hail.ship[
        #                           chunk.start.hail.era[i]:chunk.end.hail.era[i]]]
        chunk.hail.ship <- ex.ship[chunk.hail.ship.time.inds]
        #max.hail.ship[i] <- max(ex.ship[ex.dates.hail.ship[
        #                                chunk.start.hail.era[i]:chunk.end.hail.era[i]]])
        #max.hail.ship[i] <- max(chunk.hail.ship)
        #max.hail.ship.time[i] <- which(chunk.hail.ship == max(chunk.hail.ship))
        #max.hail.ship.time.ind[i] <- which(chunk.hail.ship == max(chunk.hail.ship))
        #Evaluate the indice of the max. ship
        max.hail.ship.time.ind[i] <- which.max(chunk.hail.ship)
        #max.hail.ship[i] <- max(chunk.hail.ship)
        #Evaluate the SHIP at this indice
        max.hail.ship[i] <- chunk.hail.ship[max.hail.ship.time.ind[i]]
        #Evaluate the time of the maximum SHIP
        max.hail.ship.time[i] <- strftime(rtime_era[range.dates_era[
                                                    chunk.hail.ship.time.inds[
                                                    max.hail.ship.time.ind[i]]]],
                                                    tz="UTC",
                                                    format="%Y-%m-%d %H:%M:%S")
        #print(max.ship[i])
    }
    
    #Find the maximum lightning value for each hail event date
    for (i in 1:len.chunks.hail.gpats){
        #Extract the chunk of GPATS data corresponding to each hail event date
        chunk.hail.gpats.time.inds <- ex.dates.hail.total_lightning[
                                          chunk.start.hail.gpats[i]:
                                          chunk.end.hail.gpats[i]]
        chunk.hail.total_lightning <- ex.total_lightning[chunk.hail.gpats.time.inds]
        #Evaluate the indice of the maximum lightning
        max.hail.total_lightning.time.ind[i] <- which.max(chunk.hail.total_lightning)
        #Evaluate the lightning at this indice
        max.hail.total_lightning[i] <- chunk.hail.total_lightning[
                                           max.hail.total_lightning.time.ind[i]]
        max.hail.total_lightning.time[i] <- strftime(rtime_gpats[range.dates_gpats[
                                                     chunk.hail.gpats.time.inds[
                                                    max.hail.total_lightning.time.ind[i]]]],
                                                    tz="UTC",
                                                    format="%Y-%m-%d %H:%M:%S")
    }

    ##Find the value of the SHIP when the lightning is a maximum
    for (i in 1:len.chunks.hail.gpats){
        chunk.hail.ship.time.inds <- ex.dates.hail.ship[
                                         chunk.start.hail.era[i]:
                                         chunk.end.hail.era[i]]
        chunk.hail.ship <- ex.ship[chunk.hail.ship.time.inds]
    #   #Find the closest ERA time to the max. lightning time
        ship.time.locs.hail <- as.POSIXct(rtime_era[range.dates_era
                                              [chunk.hail.ship.time.inds]],tz="UTC")
        gpats.time.locs.hail <- as.POSIXct(max.hail.total_lightning.time[i],tz="UTC")

        closest.time.ind.hail <- which.min(abs(gpats.time.locs.hail-
                                          ship.time.locs.hail))
        ship.at.max.lightning.hail[i] <- chunk.hail.ship[closest.time.ind.hail]
    }

#---#(2) Evaluate the parameters for all the NO HAIL days
    #Find the maximum SHIP value in for each hail event date
    for (j in 1:len.chunks.no.hail.era){
        #Extract the chunk of SHIP data corresponding to each hail event date
        chunk.no.hail.ship.time.inds <- ex.dates.no.hail.ship[
                                         chunk.start.no.hail.era[j]:
                                         chunk.end.no.hail.era[j]]
        chunk.no.hail.ship <- ex.ship[chunk.no.hail.ship.time.inds]
        #Evaluate the indice of the max. ship
        max.no.hail.ship.time.ind[j] <- which.max(chunk.no.hail.ship)
        #Evaluate the SHIP at this indice
        max.no.hail.ship[j] <- chunk.no.hail.ship[max.no.hail.ship.time.ind[j]]
        #Evaluate the time of the maximum SHIP
        max.no.hail.ship.time[j] <- strftime(rtime_era[range.dates_era[
                                                    chunk.no.hail.ship.time.inds[
                                                    max.no.hail.ship.time.ind[j]]]],
                                                    tz="UTC",
                                                    format="%Y-%m-%d %H:%M:%S")
        #print(max.ship[i])
    }
    
    #Find the maximum lightning value for each hail event date
    for (j in 1:len.chunks.no.hail.gpats){
        #Extract the chunk of GPATS data corresponding to each hail event date
        chunk.no.hail.gpats.time.inds <- ex.dates.no.hail.total_lightning[
                                          chunk.start.no.hail.gpats[j]:
                                          chunk.end.no.hail.gpats[j]]
        chunk.no.hail.total_lightning <- ex.total_lightning[chunk.no.hail.gpats.time.inds]
        #Evaluate the indice of the maximum lightning
        max.no.hail.total_lightning.time.ind[j]<-which.max(chunk.no.hail.total_lightning)
        #Evaluate the lightning at this indice
        max.no.hail.total_lightning[j]<-chunk.no.hail.total_lightning[
                                          max.no.hail.total_lightning.time.ind[j]]
        max.no.hail.total_lightning.time[j] <- strftime(rtime_gpats[range.dates_gpats[
                                                     chunk.no.hail.gpats.time.inds[
                                                    max.no.hail.total_lightning.time.ind[j]]]],
                                                    tz="UTC",
                                                    format="%Y-%m-%d %H:%M:%S")
    }

    ##Find the value of the SHIP when the lightning is a maximum
    for (j in 1:len.chunks.no.hail.gpats){
        chunk.no.hail.ship.time.inds <- ex.dates.no.hail.ship[
                                         chunk.start.no.hail.era[j]:
                                         chunk.end.no.hail.era[j]]
        chunk.no.hail.ship <- ex.ship[chunk.no.hail.ship.time.inds]
    #   #Find the closest ERA time to the max. lightning time
        ship.time.locs.no.hail <- as.POSIXct(rtime_era[range.dates_era
                                              [chunk.no.hail.ship.time.inds]],tz="UTC")
        gpats.time.locs.no.hail <- as.POSIXct(max.no.hail.total_lightning.time[j],tz="UTC")

        closest.time.ind.no.hail <- which.min(abs(gpats.time.locs.no.hail-
                                                  ship.time.locs.no.hail))
        ship.at.max.lightning.no.hail[j] <- chunk.no.hail.ship[closest.time.ind.no.hail]
    }


    #Now write the data to the output file
    #Include a flag to indicate hail/no hail (yes/no)
    #Write all the hail observations
    for (i in 1:len.chunks.hail.era){
        cat(ex.lon,ex.lat,
            max.hail.ship.time[i],
            max.hail.total_lightning.time[i],
            max.hail.ship[i],
            ship.at.max.lightning.hail[i],
            max.hail.total_lightning[i],
            "Yes","\n")
    }
    #Write all the no hail observations
    for (j in 1:len.chunks.no.hail.era){
        cat(ex.lon,ex.lat,
            max.no.hail.ship.time[j],
            max.no.hail.total_lightning.time[j],
            max.no.hail.ship[j],
            ship.at.max.lightning.no.hail[j],
            max.no.hail.total_lightning[j],
            "No","\n")
    }

        #cat(rtime_era[range.dates_era[ex.dates.hail.ship[max.hail.ship.ind[i]]]],
        #    max.hail.ship[i],
        #    max.hail.total_lightning[i],"\n")
    
  }
}
sink()

#Quick plot
reg_hail_data <- read.table("ship_lightning_hail_occurrence.txt",
                            header=TRUE,stringsAsFactors=FALSE)
#names(reg_hail_data) <- c("Lon","Lat","Date_ERA","Time_ERA","Date_GPATS","Time_GPATS",
#                          "Max_SHIP","SHIP_Max_GPATS","Max_GPATS","Hail_Occurred")
#attach(reg_hail_data)
plot(reg_hail_data$Max_SHIP,reg_hail_data$Max_GPATS)
points(reg_hail_data$SHIP_Max_GPATS,reg_hail_data$Max_GPATS,col="red")


#End of script
