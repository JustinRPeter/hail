#Script to analyse timeseries of the ERA_Interim data in the grid boxes
#and the inssurance data.
#Based on era_seq_analysis.R

#Justin R. Peter
#ICACS, Uni. of Southern Qld
#29 June 2017

rm(list=ls())
gc()

library(ncdf4)
library(fields)
library(oce) #Some good color scales
library(maps)
library(mapdata)
library(oz)
library(RColorBrewer)

library(RANN) #For nearest neighbour routine


library(TTR) #For time series analysis

#Specify the location of the data files

fnm_era <- '/home/u1066578/data/era_output/00_18/2008_2014/analysis/era_profiles_2008_2014_00_18.nc'

fnm_gpats <- '/home/u1066578/data/ts_lightning_grid/2008_2014/analysis/lightning_eragrid_eratime_Australia_2008_2014_offset.nc'

fnm_lsmask <- '/home/u1066578/data/eraint/lsmask/era_interim_lsmask.nc'

fnm_ins <- '/home/u1066578/data/insurance2/clean/HOME_all_claims_complete.rds'

fnm_hail <- 'hail_predictor.txt'

#Open and read the ERA data file
#Open the netcdf files
#(1) The ERA-Interim reanalysis data
    nc_era   <- nc_open(fnm_era)
    print(paste("SHARPPy profile file has",nc_era$nvars,"variables"))
    
    #Get the dimensions
    nlat_era    <- nc_era$dim$lat$len
    nlon_era    <- nc_era$dim$lon$len
    ntime_era   <- nc_era$dim$time$len
    
    #Will associate lat, lon with the ERA file
    lat <- nc_era$dim$lat$vals
    lon <- nc_era$dim$lon$vals
    time <- nc_era$dim$time$vals
    rtime  <- as.POSIXct(time*60*60, origin = "1900-01-01", tz = "GMT")  # in UTC
    
    
    #Get the variables
    #Get the SHIP variable
    v1<-nc_era$var[["ship"]]
    #ship.raw <- ncvar_get(nc_era,v1)
    shipl <- ncvar_get(nc_era,v1)
    #nc_close(nc_era)
    
    #Read the SHIP component variables
    
    v2 <- nc_era$var[["bplus"]]
    mucapel <- ncvar_get(nc_era,v2)
    
    v3 <- nc_era$var[["lr75"]]
    lr75l <- ncvar_get(nc_era,v3)
    
    v4 <- nc_era$var[["h5_temp"]]
    h5_templ <- ncvar_get(nc_era,v4)
    
    v5 <- nc_era$var[["sfc_6km_shear"]]
    s06l <- ncvar_get(nc_era,v5)
    
    v6 <- nc_era$var[["pres"]]
    presl <- ncvar_get(nc_era,v6)
    
    v7 <- nc_era$var[["dwpc"]]
    dwpcl <- ncvar_get(nc_era,v7)

    v10 <- nc_era$var[["cap"]]
    capl <- ncvar_get(nc_era,v10)

    nc_close(nc_era)
    
    #Evaluate the parcel mixing ratio
    source("/home/u1066578/jpeter/rfiles/lib/thermo/thermo.r")
    pmrl <- r_mix_td(presl*100.,dwpcl+273.15)*1000.
    
#(2) The GPATS lightning data

    nc_gpats <- nc_open(fnm_gpats)
    print(paste("GPATS lightning file has",nc_gpats$nvars,"variables"))

    #Get the total lightning variable
    v8 <- nc_gpats$var[["total_lightning"]]
    total_lightning <- ncvar_get(nc_gpats,v8)

    #Get the dimensions of the GPATS file
    nlat_gpats  <- nc_gpats$dim$lat$len
    nlon_gpats  <- nc_gpats$dim$lon$len
    ntime_gpats <- nc_gpats$dim$time$len

    lat_gpats  <- nc_gpats$dim$lat$vals
    lon_gpats  <- nc_gpats$dim$lon$vals
    time_gpats <- nc_gpats$dim$time$vals
    rtime_gpats<- as.POSIXct(time_gpats*60*60, origin="1900-01-01", tz="GMT")

    nc_close(nc_gpats)

#(3) The land-sea mask file

    nc_lsmask <- nc_open(fnm_lsmask)
    v9 <- nc_lsmask$var[["lsm"]]
    lsmask <- ncvar_get(nc_lsmask,v9)
    lat_lsmask <- nc_lsmask$dim$lat$vals
    lon_lsmask <- nc_lsmask$dim$lon$vals

    nc_close(nc_lsmask)

#-----Extract only the area around Australia for the analysis
    #Apply the same limits as the GPATS data to ERA data
    #And also apply the land-sea mask to exclude from calculation
    min_lat <- -45.0
    #max_lat <- -9.75
    max_lat <- -10.5
    min_lon <- 109.5
    max_lon <- 160.5
    
    #Extract the correct grid from the ERA data
    ex.lat <- which(lat >= min_lat & lat <= max_lat)
    ex.lon <- which(lon >= min_lon & lon <= max_lon)
    lat <- lat[ex.lat]
    lon <- lon[ex.lon]
    ship     <- shipl[ex.lon,ex.lat,]
    mucape   <- mucapel[ex.lon,ex.lat,]
    lr75     <- lr75l[ex.lon,ex.lat,]
    h5_temp  <- h5_templ[ex.lon,ex.lat,]
    s06      <- s06l[ex.lon,ex.lat,]*0.514444 #Convert to metres/second
    pmr      <- pmrl[ex.lon,ex.lat,]
    
    #Extract the correct grid from the GPATS data
    ex.lat_gpats <- which(lat_gpats >= min_lat & lat_gpats <= max_lat)
    ex.lon_gpats <- which(lon_gpats >= min_lon & lon_gpats <= max_lon)
    lat_gpats <- lat_gpats[ex.lat_gpats]
    lon_gpats <- lon_gpats[ex.lon_gpats]
    total_lightning <- total_lightning[ex.lon_gpats,ex.lat_gpats,]
    
    #Extract the correct grid from the mask file
    ex.lat_lsmask <- which(lat_lsmask >= min_lat & lat_lsmask <= max_lat)
    ex.lon_lsmask <- which(lon_lsmask >= min_lon & lon_lsmask <= max_lon)
    lat_lsmask <- lat_lsmask[ex.lat_lsmask]
    lon_lsmask <- lon_lsmask[ex.lon_lsmask]
    lsmask <- lsmask[ex.lon_lsmask,ex.lat_lsmask]

#(4) The insurance data

    fnm_ins    <- '/home/u1066578/data/insurance2/clean/HOME_all_claims_complete.rds'
#insfile   <- 'HOME_all_claims_complete.csv'
#ins_data <- read.csv(paste0(insdir,insfile),na.strings=c("NA",""))
    ins_data <- readRDS(fnm_ins)

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

#(5) The hail data
    hail_data <- read.table(fnm_hail,
                            header=TRUE,
                            stringsAsFactors=FALSE)
    names(hail_data)[20] <- "HAIL"
    hail_data$HAIL <- ifelse(hail_data$HAIL=="Yes",1,0)
    hail_data$rtime <- as.POSIXct(paste(hail_data$Date_ERA,hail_data$Time_ERA),
                                  format="%Y-%m-%d %H:%M:%S",tz="GMT")

#(6) The cluster prediction for the entire ERA-data set
    #cls.pred <- readRDS("pred_test_all_lightning_cclust10clusters.rds")
    #cls.pred <- readRDS("pred_test_all_lightning_cclust_seg_cape10clusters.rds")
    #cls.pred <- readRDS("pred_test_all_lightning_cclust_seg_cape20clusters.rds")
    #cls.pred <- readRDS("pred_test_all_lightning_cclust_seg_cape50clusters.rds")
    #cls.pred <- readRDS("pred_test_all_lightning_cclust_seg_cape100clusters.rds")

    #Using a subset of variables (MUCAPE, S06,PMR)
    #cls.pred <- readRDS("pred_test_all_lightning_cclust_seg_cape_ex_vars10clusters.rds")
    #cls.pred <- readRDS("pred_test_all_lightning_cclust_seg_cape_ex_vars20clusters.rds")
    #cls.pred <- readRDS("pred_test_all_lightning_cclust_seg_cape_ex_vars50clusters.rds")


    #cls.pred <- readRDS("pred_test_all_20clusters.rds")
    #cls.pred <- readRDS("pred_test_all_ex_vars20clusters.rds")
    cls.pred <- readRDS("pred_test_all_ex_vars_1979_2014_20clusters.rds")


    ncls <- length(levels(cls.pred))
    class(cls.pred) <- "numeric"


#----- Now loop through all of the extracted lat lon pairs
#      and evaluate which cluster it belongs to
    cls.hail <- c()
    for (i in 1:length(hail_data$Lon)){
        ilon <- which(hail_data$Lon[i] == lon)
        ilat <- which(hail_data$Lat[i] == lat)
        itime <- which(hail_data$rtime[i] == rtime)
        cls.hail[i] <- cls.pred[ilon,ilat,itime]
    }
    hail_data$CLS <- cls.hail
    

#----- Find the cluster associated with each entry from
#----- the hail data file.
    hail.inds    <- which(hail_data$HAIL == 1)
    no.hail.inds <- which(hail_data$HAIL == 0)

hail_prob <- c()
for (i in 1:ncls){
    hail_prob[i] <- length(which(hail_data$CLS[hail.inds] == i))/
                    #length(which(hail_data$CLS[no.hail.inds] == i))
                    length(which(hail_data$CLS == i))
}

cls.inds <- c()
clus.dist <- vector("list",length=ncls)
hail.dist <- vector("list",length=ncls)
nyrs <- (dim(cls.pred)[3]/(365.25*4))+1 #The '4' is because we have 4 obs per day in ERA

for (i in 1:ncls){
#for (i in 1){
    clstemp <- cls.pred
    cls.inds <- which(cls.pred == i)
    #clstemp[-cls.inds] <- 0
    #clstemp[cls.inds] <- 1
    clstemp[clstemp != i] <- 0
    clstemp[clstemp == i] <- 1
    #clus.dist[[i]] <- apply(clstemp,c(1,2),sum,na.rm=TRUE)/dim(cls.pred)[3]
    clus.dist[[i]] <- apply(clstemp,c(1,2),sum,na.rm=TRUE)/nyrs
    hail.dist[[i]] <- clus.dist[[i]]*hail_prob[i]
    #hail.dist[[i]] <- clus.dist[[i]]*hail_prob_synth[i]
}

full_clim_hail_pred <- Reduce("+",hail.dist)
for (i in 1:ncls){
    full_clim_hail_pred[lsmask==0] <- NA
}

#Save the RDA file. We will use this to write the data as netCDF file
#save(lat,lon,full_clim_hail_pred,
##save(list=c("lat","lon","full_clim_hail_pred"),
#     file=paste0("National_hail_climatology_",ncls,"_clusters.rda"))
saveRDS(lon,"longitude.rds")
saveRDS(lat,"latitude.rds")
saveRDS(full_clim_hail_pred,"hail_clim_1979_2014_ncls_20_ex_vars.rds")
#saveRDS(full_clim_hail_pred,"hail_clim_ncls_20_ex_vars.rds")

#pdf(paste0("Hail_days_",ncls,"_clusters.pdf"),width=10,height=10,paper="special")
#pdf(paste0("Hail_days_2008_2014_",ncls,"_clusters.pdf"),width=10,height=10,paper="special")
pdf(paste0("Hail_days_1979_2014_",ncls,"_clusters.pdf"),width=10,height=10,paper="special")
image.plot(lon,rev(lat),t(apply(full_clim_hail_pred,1,rev)),
           xlim=c(112,155),ylim=c(-45,-10),
           #zlim=c(0,5),
           zlim=c(0,ceiling(max(full_clim_hail_pred,na.rm=T))),
           #col=colorRampPalette(rev(brewer.pal(9,"Spectral")))(64))
           col=colorRampPalette(rev(brewer.pal(9,"Spectral")))(15))
contour(lon,rev(lat),t(apply(full_clim_hail_pred,1,rev)),
        nlevels=14,labcex=1,lwd=2,add=T)
oz(add=T,lwd=1.5)
dev.off()


#Plot some phase diagrams of the 
#Now extract the grid points corresponding to the major cities
#Brisbane, Sydney, Melbourne, Perth, Adelaide
#Next locations are in the above order

grid_lon <- c(153.0251,151.2093,144.9631,115.8605,138.6007)
grid_lat <- c(-27.4698,-33.8688,-37.8136,-31.9505,-34.9285)
ngrids <- length(grid_lon)
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

#Use time series analysis
#See https://stackoverflow.com/questions/8886677/aggregating-time-series-in-r
library(xts)
lr75.ts <- vector("list",length=ngrids)
pmr.ts <- vector("list",length=ngrids)
cls.ts <- vector("list",length=ngrids)

for (i in 1: ngrids){
     lr75.ts[[i]] <- xts(lr75[ex.lon.inds[i],ex.lat.inds[i],],rtime)
     pmr.ts[[i]]  <- xts(pmr[ex.lon.inds[i],ex.lat.inds[i],],rtime)
     cls.ts[[i]]  <- xts(cls.pred[ex.lon.inds[i],ex.lat.inds[i],],rtime)
}

#Can also evaluate average for each month
#See https://stackoverflow.com/questions/20222124/aggregate-time-series-from-weeks-to-month
mavg.lr75.ts <- vector("list",length=ngrids)
mavg.pmr.ts  <- vector("list",length=ngrids)

for (i in 1:ngrids){
    mavg.lr75.ts[[i]] <- apply.monthly(lr75.ts[[i]],mean)
    mavg.pmr.ts[[i]] <- apply.monthly(pmr.ts[[i]],mean)
}

#Evaluate monthly averages
#See https://stackoverflow.com/questions/23775683/calculate-average-value-over-multiple-years-for-each-hour-and-day
library(lubridate)
mavg.lr75 <- vector("list",length=ngrids)
mavg.pmr  <- vector("list",length=ngrids)
mavg.cls  <- vector("list",length=ngrids)
for (i in 1:ngrids){
    mavg.lr75[[i]] <- aggregate(lr75.ts[[i]],month(index(lr75.ts[[i]])),mean)
    mavg.pmr[[i]] <-  aggregate(pmr.ts[[i]],month(index(pmr.ts[[i]])),mean)
    #Does next line make sense? We averaging a clusetr numbers!
    mavg.cls[[i]] <-  aggregate(cls.ts[[i]],month(index(cls.ts[[i]])),mean)

    #mavg.lr75[[i]] <- aggregate(lr75.ts[[i]],week(index(lr75.ts[[i]])),mean)
    #mavg.pmr[[i]] <-  aggregate(pmr.ts[[i]],week(index(pmr.ts[[i]])),mean)
    ##Does next line make sense?
    #mavg.cls[[i]] <-  aggregate(cls.ts[[i]],week(index(cls.ts[[i]])),mean)

    #mavg.lr75[[i]] <- aggregate(lr75.ts[[i]],yday(index(lr75.ts[[i]])),mean)
    #mavg.pmr[[i]] <-  aggregate(pmr.ts[[i]],yday(index(pmr.ts[[i]])),mean)
    ##Does next line make sense?
    #mavg.cls[[i]] <-  aggregate(cls.ts[[i]],yday(index(cls.ts[[i]])),mean)

}

cols <- palette()
pdf("Phase_diagrams.pdf",width=10,height=10,paper="special")
    plot(mavg.pmr[[1]][1],mavg.lr75[[1]][1],xlim=c(0,15),ylim=c(5.5,7.5),col=cols[1])
    for (i in 1:ngrids){
        #points(mavg.pmr[[i]][1],mavg.lr75[[i]][1],col=cols[i])
        points(mavg.pmr[[i]][1],mavg.lr75[[i]][1],col=cols[i],pch=".")
        #for (j in 2:12){
        for (j in 2:length(mavg.cls[[1]])){
            #points(mavg.pmr[[i]][j],mavg.lr75[[i]][j],col=cols[i])
            points(mavg.pmr[[i]][j],mavg.lr75[[i]][j],col=cols[i],pch=".")
            segments(mavg.pmr[[i]][j-1],mavg.lr75[[i]][j-1],
                     #mavg.pmr[[i]][j],mavg.lr75[[i]][j],col=cols[i],lwd=3)
                     mavg.pmr[[i]][j],mavg.lr75[[i]][j],col=cols[i],lwd=0.5)
        }
    }
dev.off()



#Have commented out below as at 1/12/2017, so that I can run the
#above as a script. Uncomment later?
#Inserted one column of comments (#)


#
##    ndays_clus.dist[[i]][lsmask==0] <- NA
##        clus.dist[[i]][lsmask==0] <- NA
##            full_clim_clus[lsmask==0] <- NA
##                ndays_lightning[lsmask==0] <- NA
##}
#
#
#
#    co.times <- which(rtime %in% hail_data$rtime)
#    hail.times <- hail_data$rtime[which(hail_data$HAIL == 1)]
#    co.hail.times <- which(rtime %in% hail.times)
#    no.hail.times <- hail_data$rtime[which(hail_data$HAIL == 0)]
#    co.no.hail.times <- which(rtime %in% no.hail.times)
#
##Now extract those from the cluster prediction file
#    cls.pred.co.time <- cls.pred[,,co.times]
#    cls.pred.hail.times <- cls.pred[,,co.hail.times]
#    cls.pred.no.hail.times <- cls.pred[,,co.no.hail.times]
#
#
#    
#
#
###We need to transform the matrices for processing with R
###We need to transform this matrix
##ship <- array(NA, dim=(c(nlon,nlat,ntime)))
##
##for (i in 1:ntime){
##    ship[,,i] <- t(apply(ship.raw[,,i],1,rev))
##}
#
###Now sort the lightning and time
###Get the indices of the sorted times
##timel.sort.int <- sort.int(timel,index.return=TRUE)
##timel.sort            <- timel[timel.sort.int$ix]
##total_lightning.sort <- total_lightning[,,timel.sort.int$ix]
##
###We need to apply an offset the GPATS times.
###This is due to the NCO averaging routine which applies the time at
###the end of the averaging period
##timel.sort <- timel.sort-1.5
##rtimel <- as.POSIXct(timel.sort*60*60, origin = "1900-01-01", tz = "GMT")  # in UTC
###We need to transform this matrix
##tot_lightning <- array(NA, dim=(c(nlon,nlat,ntimel)))
##
##for (i in 1:ntimel){
##    tot_lightning[,,i] <- t(apply(total_lightning.sort[,,i],1,rev))
##}
##
##Find the time when the ERA and GPATS correspond
##Find the indice of the ERA time corresponsding to the
##start of the GPAT file
##era.time.start.ind <- which(rtime == rtimel[1])
#same.times.era.inds   <- which(rtime %in% rtimel)
#same.times.gpats.inds <- which(rtimel %in% rtime)
#
#ship.cotime          <- ship[,,same.times.era.inds]
#tot_lightning.cotime <- tot_lightning[,,same.times.gpats.inds]
#co.rtime             <- rtime[same.times.era.inds]
#
##Now just look at positive ship
##ship.gpats <- ship[which(ship.gpats > 1.0)]
##floop<-which(tot_lightning.era[3,3,]>50.0)
##plot(ship.gpats[3,3,floop],tot_lightning.era[3,3,floop])
##fit.ship.gpats <- lm(tot_lightning.era[3,3,floop] ~ ship.gpats[3,3,floop] 
##abline(lm(tot_lightning.era[3,3,floop] ~ ship.gpats[3,3,floop]))
#
##Average over the grid points in the domain
##ship.mean <- apply(ship,3,mean,na.rm=TRUE)
#ship.cotime.mean <- apply(ship.cotime,3,mean,na.rm=TRUE)
#tot_lightning.cotime.mean <- apply(tot_lightning.cotime,3,mean,na.rm=TRUE)
##Construct a running mean
#run.mean.window <- 7
##run.mean.window <- 1
#runmean.ship.cotime.mean <- SMA(ship.cotime.mean,n=run.mean.window)
#runmean.tot_lightning.cotime.mean <- SMA(tot_lightning.cotime.mean,n=run.mean.window)
##Convert to a timeseries for analysis
#ship.cotime.mean.ts <- ts(ship.cotime.mean,frequency=365.25*4,start=c(2008))
#runmean.ship.cotime.mean.ts <- ts(runmean.ship.cotime.mean,frequency=365.25*4,start=c(2008))
#tot_lightning.cotime.mean.ts <- ts(tot_lightning.cotime.mean,frequency=365.25*4,start=c(2008))
#runmean.tot_lightning.cotime.mean.ts <- ts(runmean.tot_lightning.cotime.mean,frequency=365.25*4,start=c(2008))
#
##Plot each of the individual time series
#plot(ship.cotime.mean.ts,
#     xlab="Time",ylab="SHIP")
#lines(runmean.ship.cotime.mean.ts,col="red")
#
#plot(tot_lightning.cotime.mean.ts,
#     xlab="Time",ylab="Lightning strikes")
#lines(runmean.tot_lightning.cotime.mean.ts,col="red")
#
##Plot them on the same graph
#pdf("ship_gpats_time_series_2008_2014.pdf",width=10,height=10,paper="special")
#par(mar = c(5,5,2,5))
#plot(runmean.ship.cotime.mean.ts,ylab="SHIP",type="l")
#par(new = T)
#plot(runmean.tot_lightning.cotime.mean.ts,
#     col="red",type="l",lwd=3,
#     axes=F, xlab=NA, ylab=NA)
#axis(side = 4)
#mtext(side = 4, line = 3, 'Lightning strikes')
#dev.off()
#
##Extract only the "06" data
##ex.06.times <- which(substr(rtime,12,13) == "06")
#ex.06.times <- which(substr(co.rtime,12,13) == "06")
#
##Only look at when we have a lightning strike
#lightning.occurred.all <- which(tot_lightning.cotime.mean >= 10.0)
#lightning.occurred.06  <- which(tot_lightning.cotime.mean[ex.06.times] >= 10.0)
#
##Original data. SHIP vs Lightning stikes
##wehn more than one lightning strike measured.
##Examine all times and 06 times.
#plot(tot_lightning.cotime.mean[lightning.occurred.all],
#     ship.cotime.mean[lightning.occurred.all])
#plot(tot_lightning.cotime.mean[lightning.occurred.06],
#     ship.cotime.mean[lightning.occurred.06])
##Same as above but looking at the rolling mean period chosen
#pdf("Weekly_averaged_ship_lightning.pdf",width=10,height=10,paper="special")
#plot(runmean.tot_lightning.cotime.mean[lightning.occurred.all],
#     runmean.ship.cotime.mean[lightning.occurred.all],
#     main="Weekly averaged ship vs lightning > 10 strikes/hour",
#     xlab="Lightning strikes / hour",
#     ylab="SHIP")
##plot(runmean.tot_lightning.cotime.mean[lightning.occurred.06],
##     runmean.ship.cotime.mean[lightning.occurred.06])
#
##Fit a linear model to the data
#gpats.ship.model <- lm(runmean.ship.cotime.mean[lightning.occurred.all] ~ runmean.tot_lightning.cotime.mean[lightning.occurred.all] + 0)
##abline(gpats.ship.model$coefficients[1],gpats.ship.model$coefficients[2])
#abline(0,gpats.ship.model$coefficients[1])
#box()
#dev.off()
#
#
#
#
##Extract all the ship components corresponding to the 2014-11-27 hail storm
##storm <- which(as.Date(rtime) == "2014-11-27")
##image.plot(lon,rev(lat),ship[,,storm[2]])
#
##-----Open the insurance data and find all the unique dates
#
#insdir    <- '/home/jpeter/justinp/rfiles/suncorp/data/insurance2/clean/'
#insfile   <- 'HOME_all_claims_complete.rds'
##insfile   <- 'HOME_all_claims_complete.csv'
##ins_data <- read.csv(paste0(insdir,insfile),na.strings=c("NA",""))
#ins_data <- readRDS(paste0(insdir,insfile))
#
#ins_data <- ins_data[!is.na(ins_data$Postcode),]
#ins_data <- ins_data[!is.na(ins_data$Address),]
##ins_data <- ins_data[ins_data$Longitude>=min(lon)&ins_data$Longitude<=max(lon),]
##ins_data <- ins_data[ins_data$Latitude>=min(lat) & ins_data$Latitude<=max(lat),]
#ins_data <- ins_data[ins_data$Longitude>=152.625 &ins_data$Longitude<=153.375,]
#ins_data <- ins_data[ins_data$Latitude>= -28.125 & ins_data$Latitude<= -27.375,]
#ins_data <- ins_data[ins_data$Incurred.To.Date > 10,]
#
#ins_data$ins.date <- as.Date(ins_data$Date.Occurred,"%d/%m/%Y")
#
##Find the unique dates
#uniq.ins.dates <- unique(ins_data$ins.date)
#
##Find the times in the ship data which correspond to damage days
##ex.dates <- which(as.Date(rtime[ex.06.times]) %in% uniq.ins.dates)
##ex.dates <- which(as.Date(co.rtime[ex.06.times]) %in% uniq.ins.dates)
#ex.dates <- which(as.Date(co.rtime) %in% uniq.ins.dates)
#
##Extract pertinent times from the SHIP
##Use the 06 ship sounding but lightning at any time!
##ship.claim.days <- ship[,,ex.06.times[ex.dates]]
##ship.claim.days <- ship.cotime[,,ex.06.times[ex.dates]]
##tot_lightning.claim.days <- tot_lightning.cotime[,,ex.06.times[ex.dates]] 
#ship.claim.days <- ship.cotime[,,ex.dates]
#tot_lightning.claim.days <- tot_lightning.cotime[,,ex.dates] 
#
#mean.ship.claim.days <- apply(ship.claim.days,c(1,2),mean,na.rm=TRUE)
#mean.tot_lightning.claim.days <- apply(tot_lightning.claim.days,c(1,2),mean,na.rm=TRUE)
#
##Only use the "06" data
##ex.06.times <- which(substr(rtime[exdates],12,13) == "06")
#
##Plot the insurance claims on the extracted region
#minx=min(lon)
#maxx=max(lon)
#miny=min(lat)
#maxy=max(lon)
#
##png("mean_ship_on_claim_days.png",
#pdf("mean_ship_on_claim_days.pdf",width=11,height=10,paper="special")
##image.plot(lon,rev(lat),ship.claim.days[,,1])
#image.plot(lon,rev(lat),mean.ship.claim.days[,],
#           zlim=c(0,max(mean.ship.claim.days)),
#           main="Avg. SHIP on claim days in Brisbane grid",
#           xlab="Longitude",ylab="Latitude",
#           legend.lab="SHIP",
#           col=oce.colorsViridis(64),
#           cex.axis=1.25,cex.lab=1.25)
#           #col=topo.colors(64))
#           #col=terrain.colors(64))
#map('worldHires',interior=T,xlim=c(minx,maxx),ylim=c(miny,maxy),lwd=2.5,add=T)
##map.cities(label=T,minpop=5.e04,col="darkblue",cex=1.35,lwd=2)
#points(ins_data$Longitude,ins_data$Latitude)
#map.cities(label=T,minpop=5.e04,col="darkblue",cex=1.35,lwd=2)
#box()
#dev.off()
#
#pdf("mean_tot_lightning_on_claim_days.pdf",width=11,height=10,paper="special")
##image.plot(lon,rev(lat),ship.claim.days[,,1])
#image.plot(lon,rev(lat),mean.tot_lightning.claim.days[,]/2,
#           zlim=c(0,max(mean.tot_lightning.claim.days)/2.),
#           main="Avg. lightning strikes on claim days in Brisbane grid",
#           xlab="Longitude",ylab="Latitude",
#           legend.lab="Lightning strikes / hour",
#           col=oce.colorsViridis(64),
#           cex.axis=1.25,cex.lab=1.25)
#           #col=topo.colors(64))
#           #col=terrain.colors(64))
#map('worldHires',interior=T,xlim=c(minx,maxx),ylim=c(miny,maxy),lwd=2.5,add=T)
##map.cities(label=T,minpop=5.e04,col="darkblue",cex=1.35,lwd=2)
#points(ins_data$Longitude,ins_data$Latitude)
#map.cities(label=T,minpop=5.e04,col="darkblue",cex=1.35,lwd=2)
#box()
#dev.off()
#
#
#
#pdf("mean_bris_ship_tseries_withclaim_days.pdf",width=10,height=10,paper="special")
##Now plot a time series of SHIP and overlay hail claim days.
##date.min = min(as.Date(rtime[ex.06.times][ex.dates]))
##date.max = max(as.Date(rtime[ex.06.times][ex.dates]))
#date.min = min(as.Date(co.rtime[ex.dates]))
#date.max = max(as.Date(co.rtime[ex.dates]))
##See http://stackoverflow.com/questions/4355042/label-x-axis-in-time-series-plot-using-r
##x.Date <- as.Date(paste(2009:2014,01, 01, sep = "-"))
#x.Date <- as.Date(paste(2009:2015,01, 01, sep = "-"))
##plot(as.Date(rtime[ex.06.times]),ship[3,2,ex.06.times],
##Use the mean of the grid boxes over Brisbane and the grid box north
##ship.avg.bris <- apply(ship[3,c(2,3),ex.06.times],2,mean)
##ship.avg.bris <- ship[3,3,ex.06.times]
##ship.avg.bris <- ship.cotime[3,3,ex.06.times]
#ship.avg.bris <- ship.cotime[3,3,]
##ship.avg.bris <- ship.cotime[3,3,ex.dates]
##plot(as.Date(rtime[ex.06.times]),ship[3,2,ex.06.times],
##plot(as.Date(rtime[ex.06.times]),ship.avg.bris,
##plot(as.Date(co.rtime[ex.06.times]),ship.avg.bris,
##Plot SHIP at any time on that day
#plot(as.Date(co.rtime),ship.avg.bris,
#     xlim=c(date.min,date.max),
#     main="Time series of SHIP averaged over Brisbane grid",
#     xaxt="n",type="l",
#     xlab="Date",ylab="SHIP")
#axis(1,at=x.Date,label=paste(2009:2015),las=0)
##points(as.Date(rtime[ex.06.times][ex.dates]),ship.claim.days[3,3,],col="red",lwd=3)
##points(as.Date(rtime[ex.06.times][ex.dates]),ship.claim.days[3,3,],col="red",lwd=3)
##points(as.Date(co.rtime[ex.06.times][ex.dates]),ship.claim.days[3,3,],col="red",lwd=3)
#points(as.Date(co.rtime[ex.dates]),ship.claim.days[3,3,],col="red",lwd=3)
#dev.off()
#
#pdf("mean_bris_lightning_tseries_withclaim_days.pdf",width=10,height=10,paper="special")
##Now plot a time series of SHIP and overlay hail claim days.
##date.min = min(as.Date(rtime[ex.06.times][ex.dates]))
##date.max = max(as.Date(rtime[ex.06.times][ex.dates]))
#date.min = min(as.Date(co.rtime[ex.dates]))
#date.max = max(as.Date(co.rtime[ex.dates]))
##See http://stackoverflow.com/questions/4355042/label-x-axis-in-time-series-plot-using-r
##x.Date <- as.Date(paste(2009:2014,01, 01, sep = "-"))
#x.Date <- as.Date(paste(2009:2015,01, 01, sep = "-"))
##plot(as.Date(rtime[ex.06.times]),ship[3,2,ex.06.times],
##Use the mean of the grid boxes over Brisbane and the grid box north
##ship.avg.bris <- apply(ship[3,c(2,3),ex.06.times],2,mean)
##ship.avg.bris <- ship[3,3,ex.06.times]
##tot_lightning.avg.bris <- tot_lightning.cotime[3,3,ex.06.times]
##tot_lightning.avg.bris <- tot_lightning.cotime[3,3,ex.dates]
#tot_lightning.avg.bris <- tot_lightning.cotime[3,3,]
##plot(as.Date(rtime[ex.06.times]),ship[3,2,ex.06.times],
##plot(as.Date(rtime[ex.06.times]),ship.avg.bris,
##plot(as.Date(co.rtime[ex.06.times]),tot_lightning.avg.bris/2.,
##PLot lightning at any time
##plot(as.Date(co.rtime[ex.06.times]),tot_lightning.avg.bris/2.,
#plot(as.Date(co.rtime),tot_lightning.avg.bris/2.,
#     xlim=c(date.min,date.max),
#     main="Time series of lightning strikes averaged over Brisbane grid",
#     #log="y",
#     #ylim=c(0,20),
#     xaxt="n",type="l",
#     xlab="Date",ylab="Lightning strikes / hour")
#axis(1,at=x.Date,label=paste(2009:2015),las=0)
##points(as.Date(rtime[ex.06.times][ex.dates]),ship.claim.days[3,3,],col="red",lwd=3)
##points(as.Date(rtime[ex.06.times][ex.dates]),ship.claim.days[3,3,],col="red",lwd=3)
##points(as.Date(co.rtime[ex.06.times][ex.dates]),tot_lightning.claim.days[3,3,]/2.,col="red",lwd=3)
#points(as.Date(co.rtime[ex.dates]),tot_lightning.claim.days[3,3,]/2.,col="red",lwd=3)
#dev.off()
#
##Plot SHIP versus lightning on claim days
##plot(
#
#
#pdf("histogram_ship_bne_2009_2014_claims.pdf",width=10,height=10,paper="special")
#hist(ship.claim.days,freq=F,col="red",main="Histogram of SHIP over Brisbane grid")
#hist(ship.avg.bris,freq=F,add=T,col="blue")
#legend("topright",c("2009-2014","claim days"),lty=c(1,1),col=c("blue","red"))
#dev.off()
#
#





