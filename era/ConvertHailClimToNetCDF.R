#Script to convert the rda files output by 
#era_insurance_analysis.R to netCDF files

#Justin Peter
#7 November 2017


rm(list=ls())

library(ncdf4)

#Specify the number of clusters
ncls=10

lon <- readRDS("longitude.rds")
lat <- readRDS("latitude.rds")
hail_clim <- readRDS("hail_clim.rds")


#Define dimensions
londim <- ncdim_def("lon","degrees_east",as.double(lon))
latdim <- ncdim_def("lat","degrees_north",as.double(lat))

#Define the variables
fillvalue <- 1e32
dlname <- "Hail days per year"
hail_def <-ncvar_def("hdpy","Days/yr",list(londim,latdim),fillvalue,dlname,prec="double")

#Create the netcdf file and put arrays
ncfname <- "hail_climatology.nc"
ncout <- nc_create(ncfname,list(hail_def),force_v4=T)

#Put the variables
ncvar_put(ncout,hail_def,hail_clim)

#Put some attributes
ncatt_put(ncout,"lon","axis","X")
ncatt_put(ncout,"lat","axis","Y")

ncatt_put(ncout,0,"title","Hail climatology")
ncatt_put(ncout,0,"institution","University of Southern Queensland")
#ncatt_put(ncout,0,"source",datasource$value)
#ncatt_put(ncout,0,"references",references$value)
history <- paste("Justin R. Peter",date(),sep=",")
ncatt_put(ncout,0,"history",history)
#ncatt_put(ncout,0,"Conventions",Conventions$value)


nc_close(ncout)

#nms <- c("lat","lon","hail_clim")
#data <- load("National_hail_climatology_10_clusters.rda")

#names(data) <- nms
