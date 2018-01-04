#Plot the climatology of the prediction of lightning
#based on clustering of convective variables
# ***** Search for this string to change parameters depending
#on the method of clustering cape

rm(list=ls())
gc()
library(fields)
library(oce) #Viridis colors
library(maps)
library(mapdata)
library(oz)
library(RColorBrewer)
library(RcppRoll) #For rolling sums
library(ncdf4)
library(reshape2)
library(ggplot2)
library(scales)
library(Hmisc) #for Ecdf function
library(MASS) #for neg binom glm

#Read the RDS files output from cluster_hail_prediction.R
#The cluster prediction file
#cls.pred <- readRDS("pred_test_all_lightning_6clusters.rds")
#cls.pred <- readRDS("pred_test_all_lightning_20clusters.rds")
#cls.pred <- readRDS("pred_test_all_lightning_kmeans_lloyd20clusters.rds")
#cls.pred <- readRDS("pred_test_all_lightning_cclust10clusters.rds")
#cls.pred <- readRDS("pred_test_all_lightning_cclust_seg_cape10clusters.rds")
cls.pred <- readRDS("pred_test_all_lightning_cclust_seg_cape20clusters.rds")
ncls <- length(levels(cls.pred))


#The hail and lightning probability files
#hlprob <- readRDS("prop_hail_lightning_all_6clusters.rds")
#hlprob <- readRDS("prop_lightning_all_20clusters.rds")
# *****
#hlprob <- readRDS("prop_lightning_all_10clusters.rds")
hlprob <- readRDS("prop_lightning_all_seg_cape10clusters.rds")

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
#Create a time for the cluster data
rtime_cls <- seq(rtime_gpats[1],rtime_gpats[61363],length.out=dim(cls.pred)[3])

#gap_storm_date <- "2008-11-16 06:00:00 GMT"
gap_storm_date <- "2014-11-27 06:00:00 GMT"
gap_storm_doy  <- as.numeric(strftime(gap_storm_date,format="%j"))
ind_gap_storm_date <- which(rtime_cls == as.POSIXlt(gap_storm_date,tz="GMT"))

#Get the total_lightning variable
#Get the total lightning variable
v1<-nc_gpats$var[["total_lightning"]]
total_lightning <- ncvar_get(nc_gpats,v1)
#total_lightning[,,] <- t(apply(total_lightning[,,],1,rev))

nc_close(nc_gpats)
#nc_close(nc_era)
nc_close(nc_lsmask)

#Read the ERA data

#fnm_era <-'/home/u1066578/data/era_output/00_18/2008_2014/analysis/era_profiles_2008_2014_00_18.nc'
fnm_era <- '/home/u1066578/data/era_output/00_18/2008_2014/analysis/era_profiles_2008_2014_00_18.nc'

#Open the netcdf files
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
nc_close(nc_era)

#Evaluate the parcel mixing ratio
source("/home/u1066578/jpeter/rfiles/lib/thermo/thermo.r")
pmrl <- r_mix_td(presl*100.,dwpcl+273.15)*1000.

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
mucape   <- mucapel[ex.lon,ex.lat,]
lr75     <- lr75l[ex.lon,ex.lat,]
h5_temp  <- h5_templ[ex.lon,ex.lat,]
s06      <- s06l[ex.lon,ex.lat,]
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

#Put all of the ERA variables into a data frame
#Put all these variables in a list
#era_data <- list(MUCAPE=as.vector(mucape),LR75=as.vector(lr75),
#                 H5_TEMP=as.vector(h5_temp),S06=as.vector(s06),
#                 PMR=as.vector(pmr))
era_data <- cbind(as.vector(mucape),as.vector(lr75),
                  as.vector(h5_temp),as.vector(s06),
                  as.vector(pmr))

#Next line only needed if era_data is in a list as above.
#era_data <- as.data.frame(era_data)
source("colScale.R")
era_data_n <- colScale(as.matrix(era_data))

mucape_n <- era_data_n[,1]
lr75_n <- era_data_n[,2]
h5_temp_n <- era_data_n[,3]
s06_n <- era_data_n[,4]
pmr_n <- era_data_n[,5]

#Set the dimensions of the vectors so that we can
#use only land points in the regression.
dim(mucape_n) <- dim(lr75_n) <- 
dim(h5_temp_n) <- dim(s06_n) <- 
dim(pmr_n) <- dim(mucape)


#Now extract only the time from the GPATS corresponding
#to the ERA data.
#We will use this to construct a daily proability of lightning
lightning_at_era <- apply(total_lightning,c(1,2),roll_max,3,by=6,
                       #align="center",fill=NA,na.rm=TRUE)
                       align="center",na.rm=TRUE)
lightning_at_era <- lightning_at_era[seq(1,dim(lightning_at_era)[1],6),,]
lightning_at_era <- aperm(lightning_at_era,c(2,3,1))
daily_lightning <- apply(lightning_at_era,c(1,2),roll_sum,4,by=4,
                                align="left",na.rm=TRUE)
daily_lightning <- daily_lightning[seq(1,dim(daily_lightning)[1],4),,]
daily_lightning <- aperm(daily_lightning,c(2,3,1))
daily_lightning[daily_lightning > 0] <- 1
ndays_lightning <- apply(daily_lightning,c(1,2),sum,na.rm=TRUE)

#Also construct a mean lightning
mean_lightning_at_era <- apply(total_lightning,c(1,2),roll_mean,3,by=6,
                             align="center",na.rm=TRUE)
mean_lightning_at_era <- mean_lightning_at_era[seq(1,dim(mean_lightning_at_era)[1],6),,]
mean_lightning_at_era <- aperm(mean_lightning_at_era,c(2,3,1))

#Try the total sum of lightning for the regression
#to get a more representative idea of the daily average.
max_lightning <- apply(total_lightning,c(1,2),roll_sum,6,by=6,
                       #align="center",fill=NA,na.rm=TRUE)
                       align="left",na.rm=TRUE)
max_lightning <- max_lightning[seq(1,dim(max_lightning)[1],6),,]
max_lightning <- aperm(max_lightning,c(2,3,1))


#Because we centred the window for rolling max above,
#which we need to because we want to take the max lightning
#from each hour around the ERA reanalysis
#we lose the last observation (2008-12-31 18:00:00 GMT)
#so we have to remove that from the ERA and cluster obs also
mxt <- dim(max_lightning)[3]

mucape <- mucape[,,-mxt+1]
lr75 <- lr75[,,-mxt+1]
h5_temp <- h5_temp[,,-mxt+1]
s06 <- s06[,,-mxt+1]
pmr <- pmr[,,-mxt+1]

mucape_n <- mucape_n[,,-mxt+1]
lr75_n <- lr75_n[,,-mxt]
h5_temp_n <- h5_temp_n[,,-mxt+1]
s06_n <- s06_n[,,-mxt+1]
pmr_n <- pmr_n[,,-mxt+1]

cls.pred <- cls.pred[,,-mxt+1]

#Now apply the land sea mask to the ERA and GPATS data
#as we only want to include land points in our analysis.
#The regression data only contains land points!

#We need to convert cls.pred to numeric
#as it is a factor
class(cls.pred) <- "numeric"

# Define a mask function
mask <- function(x, y) {
            x[y==0] <- NA 
            x 
        } 

cls.pred <- simplify2array(lapply(1:mxt, function(i) mask(cls.pred[,,i],lsmask)))
mucape_n <- simplify2array(lapply(1:mxt, function(i) mask(mucape_n[,,i],lsmask)))
lr75_n <- simplify2array(lapply(1:mxt, function(i) mask(lr75_n[,,i],lsmask)))
h5_temp_n <- simplify2array(lapply(1:mxt, function(i) mask(h5_temp_n[,,i],lsmask)))
s06_n <- simplify2array(lapply(1:mxt, function(i) mask(s06_n[,,i],lsmask)))
pmr_n <- simplify2array(lapply(1:mxt, function(i) mask(pmr_n[,,i],lsmask)))
max_lightning <- simplify2array(lapply(1:mxt, function(i) mask(max_lightning[,,i],lsmask)))
lightning_at_era <- simplify2array(lapply(1:mxt,
                        function(i) mask(lightning_at_era[,,i],lsmask)))
mean_lightning_at_era <- simplify2array(lapply(1:mxt,
                        function(i) mask(mean_lightning_at_era[,,i],lsmask)))

#Use a flag to use different lightning matrices for the regression
#lightning_matrix <- mean_lightning_at_era
#lightning_matrix <- lightning_at_era
lightning_matrix <- max_lightning

reg_data_n <- data.frame(CLS=as.vector(cls.pred),
                       CAPE=as.vector(mucape_n),
                       LR75=as.vector(lr75_n),
                       H5_T=as.vector(h5_temp_n),
                       S06=as.vector(s06_n),
                       PMR=as.vector(pmr_n),
                       #LGHT=as.vector(max_lightning))
                       #LGHT=as.vector(mean_lightning_era))
                       LGHT=as.vector(lightning_matrix))

#Remove all the NAs from the data set
reg_data_n <- reg_data_n[complete.cases(reg_data_n),]

ncape <- mucape_n*attr(era_data_n,"scaled:scale")[1]+attr(era_data_n,"scaled:center")[1]
nlr75 <- lr75_n*attr(era_data_n,"scaled:scale")[2]+attr(era_data_n,"scaled:center")[2]
nh5_t <- h5_temp_n*attr(era_data_n,"scaled:scale")[3]+attr(era_data_n,"scaled:center")[3]
ns06 <- s06_n*attr(era_data_n,"scaled:scale")[4]+attr(era_data_n,"scaled:center")[4]
npmr <- pmr_n*attr(era_data_n,"scaled:scale")[5]+attr(era_data_n,"scaled:center")[5]

reg_data <- data.frame(CLS=as.vector(cls.pred),
                       CAPE=as.vector(ncape),
                       LR75=as.vector(nlr75),
                       H5_T=as.vector(nh5_t),
                       S06=as.vector(ns06),
                       PMR=as.vector(npmr),
                       #LGHT=as.vector(max_lightning))
                       #LGHT=as.vector(mean_lightning_at_era))
                       LGHT=as.vector(lightning_matrix))
                       #y=as.vector(max_lightning))

#Remove all the NAS from the data set
reg_data <- reg_data[complete.cases(reg_data),]



#For loops are too slow.
#for (i in 1:mxt){
#    cls.pred[,,i][which(lsmask==0)] <- NA
#    mucape_n[,,i][which(lsmask==0)] <- NA
#    lr75_n[,,i][which(lsmask==0)] <- NA
#    h5_temp_n[,,i][which(lsmask==0)] <- NA
#    s06_n[,,i][which(lsmask==0)] <- NA
#    pmr_n[,,i][which(lsmask==0)] <- NA
#    max_lightning[,,i][which(lsmask==0)] <- NA
#}

#Now split into a training and test data set
#library(dplyr)
#train <- sample_frac(reg_data,0.7)
#test <- sample_frac(reg_data,0.3)
#sid<-as.numeric(rownames(train))
#test<-reg_data[-sid,]

#Using the original data
#dt <- sort(sample(nrow(reg_data),nrow(reg_data)*.7))
#train <- reg_data[dt,]
#test <- reg_data[-dt,]

##Using the normalised data
dt <- sort(sample(nrow(reg_data_n),nrow(reg_data_n)*.7))
train <- reg_data_n[dt,]
test <- reg_data_n[-dt,]



#Split the training data set by clusters
train_cls <- vector("list",length=ncls)
test_cls <- vector("list",length=ncls)
for (i in 1:ncls){
    train_inds <- which(train$CLS == i)
    test_inds <- which(test$CLS == i)
    train_cls[[i]] <- train[train_inds,]
    test_cls[[i]] <- test[test_inds,]
}



#---------- Some info about the mean and sd of the training data set
library(data.table)
setDT(train)
mn_lght <- train[,lapply(.SD,mean,na.rm=TRUE),by=CLS,.SDcols="LGHT"]
sd_lght <- train[,lapply(.SD,sd,na.rm=TRUE),by=CLS,.SDcols="LGHT"]

mn_lght <- mn_lght[order(mn_lght[["CLS"]])]
sd_lght <- sd_lght[order(sd_lght[["CLS"]])]

mn_lght_gt0 <- train[LGHT>0,lapply(.SD,mean,na.rm=TRUE),by=CLS,.SDcols="LGHT"]
sd_lght_gt0 <- train[,lapply(.SD,sd,na.rm=TRUE),by=CLS,.SDcols="LGHT"]
mn_lght_gt0 <- mn_lght_gt0[order(mn_lght_gt0[["CLS"]])]
sd_lght_gt0 <- sd_lght_gt0[order(sd_lght_gt0[["CLS"]])]

#mn_lght_cls <- vector("list",length=ncls)
#sd_lght_cls <- vector("list",length=ncls)


#----------
#Add a lightning indicator factor to the data

for (i in 1:ncls){
    ltng <- vector("numeric",length=length(train_cls[[i]]$LGHT))
    linds <- which(train_cls[[i]]$LGHT > 0)
    ltng[linds] <- 1
    ltng[-linds] <- 0
    train_cls[[i]]$Lightning <- ltng
    #train_cls[[i]]$Lightning <- factor(train_cls[[i]]$Lightning)
}

#Evaluate the proability of each cluster producing lightning
prob.lightning <- vector("numeric",length=ncls)
for (i in 1:ncls){
    prob.lightning[i] <- length(which(train_cls[[i]]$Lightning==1))*100./
                           length(train_cls[[i]]$Lightning)
}

#Comment this out as it takes a long time!!!!!
##Plot the pdfs of each of the convective variables
#dt1 <- sort(sample(nrow(reg_data),nrow(reg_data)*.1))
#df <- reg_data[dt1,1:6]
#df$CLS <- factor(df$CLS)
#df$CAPE <- log(df$CAPE)
###test <- reg_data_n[-dt,]
##
###df <- reg_data[,1:6]
#dfg1 <- melt(df,id.vars=c("CLS"),
#               measure.vars=c("CAPE","LR75","H5_T","S06","PMR"))
#
##labels <- c(CAPE="MUCAPE [J/kg]", 
#labels <- c(CAPE="log(MUCAPE)", 
#            LR75="LR75 [C/km]",
#            H5_T= "T_H500 [C]",
#            S06="S06",
#            PMR="PMR [kg/kg]")
#pvars <- ggplot(dfg1,aes(x=value))+
#      xlab("")+
#      #stat_bin(aes(y=..density..,color=cluster),geom="step",size=1.5,position="identity")+
#      geom_freqpoly(aes(y=..density..,color=CLS),size=1.5)+
#      #stat_bin(aes(y=..count..,color=cluster),geom="step",size=1.5)+
#      facet_wrap(~variable,ncol=2,scale="free",
#                 labeller=labeller(variable=labels),
#                 #switch="y")+
#                 strip.position="left")+
#      #scale_colour_manual(values=wes_palette(n=length(ind),name="Zissou"))+
#      #scale_colour_manual(values=rev(oce.colorsViridis(length(ind))))+
#      scale_colour_manual(values=rev(oce.colorsViridis(ncls)))+
#      guides(colour=guide_legend(override.aes=list(size=5)))+
#      theme(axis.text.x=element_text(face="bold",size=12),
#            axis.text.y=element_text(face="bold",size=12),
#            axis.title.x=element_text(face="bold",size=14),
#            strip.text.y=element_text(face="bold",size=12))
#
#
###plot(p3)
###png("hists_vars_by_cls.png",width=xsize,height=ysize,units="in",res=300)
#pdf("hists_vars_by_cls.pdf",width=10,height=10,paper="special")
#print(pvars)
#dev.off()
##




##Now plot the mean lightning per cluster and the probability
##of lightning per cluster
#pdf("Cluster_lightning_probabilities.pdf",width=10,height=10,paper="special")
#par(mfrow=c(2,1))
#par(mar=c(4,4,1,1)+0.1)
#par(cex=1.1)
#barplot(prob.lightning,names=1:ncls,
#        ylab="Percent of cluster with lightning")
##text(1,35,"(a)",cex=1.5,adj=c(1,1))
#text(1,35,"(a)",cex=1.2,adj=c(1,1))
#barplot(mn_lght_gt0$LGHT,names=1:ncls,
#        xlab="Cluster",
#        ylab="Avg. lightning counts per hour")
##text(1,200,"(b)",cex=1.5,adj=c(1,1))
#text(1,200,"(b)",cex=1.2,adj=c(1,1))
#dev.off()
#


#First get rid of all the zeros for the regression
for (i in 1:ncls){
    train_zero_inds <- which(train_cls[[i]]$LGHT == 0)
    test_zero_inds  <- which(test_cls[[i]]$LGHT == 0)
    train_cls[[i]] <- train_cls[[i]][-train_zero_inds,]
    test_cls[[i]] <- test_cls[[i]][-test_zero_inds,]
}


#Construct cumulative histograms

max_lightning_t <- max_lightning
#max_lightning_t[,,] <- t(apply(max_lightning[,,],1,rev))
max_lightning_t[!is.finite(max_lightning_t)] <- NA

#colsh <- tim.colors(ncls) #Colours for the histograms
colsh <- rev(oce.colorsViridis(ncls)) #Colours for the histograms
tot_lght <- vector("list",length=ncls+1) #One extra for all clusters
tot_lght[[1]] <- total_lightning[total_lightning>0]
for (i in 1:ncls){
    #tot_lght[[i+1]] <- max_lightning_t[which(cluster==i)]
    tot_lght[[i+1]] <- max_lightning_t[which(cls.pred==i)]
    tot_lght[[i+1]] <- tot_lght[[i+1]][tot_lght[[i+1]]>0]
}
cum_hists <- vector("list",length=ncls+1) #One extra for all clusters


##Plot the cumulative histogram of the lightning
# *****
# Get the names right here!
#pdf("cum_hist_lightning.pdf",width=10,height=10)
pdf("cum_hist_lightning_seg_cape.pdf",width=10,height=10)
par(cex=1.5,mar=c(5,4,2,2))
#Plot the total lightning
cum_hists[[1]] <- Ecdf(tot_lght[[1]],
                       what="1-F",
                       xlab="Lightning strikes/hour",
                       ylab="Fraction",
                       log="xy",
                       ylim=c(1.e-06,1),
                       #yaxt="n",
                       lwd=5,
                       subtitles=F)
for (i in 1:ncls){
cum_hists[[i+1]] <- Ecdf(tot_lght[[i+1]],
                       what="1-F",
                       lwd=5,
                       subtitles=F,
                       col=colsh[i],add=TRUE)
}
legend("topright",legend=1:ncls,fill=colsh)
dev.off()


#Examine the variance mean ratio of the lightning counts
var_mean_ratio <- c()
for (i in 1:ncls){
    var_mean_ratio[i] <- var(train_cls[[i]]$LGHT)/mean(train_cls[[i]]$LGHT)
    #var_mean_ratio[i] <- sd(train_cls[[i]]$LGHT)/mean(train_cls[[i]]$LGHT)
}



##---------- Step to find best model selection
##Trying Poisson, Quasipoisson and Negative binomial
#
#null.p <- vector("list",length=ncls)
#full.p <- vector("list",length=ncls)
##backward.p <- vector("list",length=ncls)
#forward.p <- vector("list",length=ncls)
##backward.AIC.p <- vector("list",length=ncls)
#forward.AIC.p <- vector("list",length=ncls)
#GPATS_best_glm.p <- vector("list",length=ncls)
#null.qp <- vector("list",length=ncls)
#full.qp <- vector("list",length=ncls)
##backward.qp <- vector("list",length=ncls)
#forward.qp <- vector("list",length=ncls)
#GPATS_best_glm.qp <- vector("list",length=ncls)
#null.nb <- vector("list",length=ncls)
#full.nb <- vector("list",length=ncls)
##backward.nb <- vector("list",length=ncls)
#forward.nb <- vector("list",length=ncls)
#GPATS_best_glm.nb <- vector("list",length=ncls)
#
#
##Use step function to step through
##Poisson
#for (i in 1:ncls){
#    null.p[[i]] <- glm(LGHT~1,family=poisson,data=train_cls[[i]])
#    full.p[[i]] <- glm(LGHT ~ CAPE+LR75+H5_T+S06+PMR,family=poisson,data=train_cls[[i]])
#    forward.p[[i]] <- step(null.p[[i]],data=train_cls[[i]],
#                           scope=list(lower=null.p[[i]], upper=full.p[[i]]),
#                           direction="forward",
#                           trace=TRUE)
#    #GPATS_best_glm.p[[i]] <- backward[[i]]
#    GPATS_best_glm.p[[i]] <- forward.p[[i]]
#}
#
##Quasipoisson
##Can't use AIC criterion for Quasipoisson model fitting
#for (i in 1:ncls){
#    null.qp[[i]] <- glm(log(LGHT)~1,family=quasipoisson,data=train_cls[[i]])
#    full.qp[[i]] <- glm(log(LGHT) ~ CAPE+LR75+H5_T+S06+PMR,family=quasipoisson,data=train_cls[[i]])
#    forward.qp[[i]] <- stats::step(null.qp[[i]],data=train_cls[[i]],
#                           scope=list(lower=null.qp[[i]], upper=full.qp[[i]]),
#                           direction="forward",
#                           trace=TRUE)
#
#    GPATS_best_glm.qp[[i]] <- forward.qp[[i]]
#}
#
##Negative binomial
#for (i in 1:ncls){
##for (i in 10){
#    null.nb[[i]] <- glm.nb(log(LGHT)~1,data=train_cls[[i]])
#    full.nb[[i]] <- glm.nb(log(LGHT) ~ CAPE+LR75+H5_T+S06+PMR,data=train_cls[[i]])
#    forward.nb[[i]] <- step(null.nb[[i]],data=train_cls[[i]],
#                           scope=list(lower=null.nb[[i]], upper=full.nb[[i]]),
#                           direction="forward",
#                           trace=TRUE)
#    GPATS_best_glm.nb[[i]] <- forward.nb[[i]]
#}
#


#fit.p        <- vector("list",length=ncls)
#fit.qp       <- vector("list",length=ncls)
fit.nb       <- vector("list",length=ncls)
#fit.hurdle   <- vector("list",length=ncls)
#fit.zinfl    <- vector("list",length=ncls)
#pred.p <- vector("list",length=ncls)
#pred.qp <- vector("list",length=ncls)
pred.nb <- vector("list",length=ncls)
#pred.hurdle <- vector("list",length=ncls)

#WARNING! Be careful when running the following as ther hurdle model takes a
#while. Do tests just using a single cluster.
for (i in 1:ncls){
#for (i in 10){
    #fit.p[[i]]      <- glm(LGHT~CAPE+LR75+H5_T+S06+PMR,data=train_cls[[i]],family="poisson",trace=TRUE)
    #fit.qp[[i]]     <- glm(LGHT~CAPE+LR75+H5_T+S06+PMR,data=train_cls[[i]],
    #                       family="quasipoisson",trace=TRUE)
    fit.nb[[i]]     <- glm.nb(LGHT~CAPE+LR75+H5_T+S06+PMR,data=train_cls[[i]],trace=TRUE,link=log)
    #fit.hurdle[[i]] <- hurdle(LGHT~CAPE+LR75+H5_T+S06+PMR,data=train_cls[[i]],dist="negbin")
    ##Zero inflated models take too long!
    ##fit.zinfl[[i]]  <- zeroinfl(LGHT~CAPE+LR75+H5_T+S06+PMR,data=train_cls[[i]],dist="negbin")
}

for (i in 1:10){
    ##pred.p[[i]]  <-  predict(fit.p[[i]],test_cls[[i]])
    ##pred.qp[[i]] <- predict(fit.qp[[i]],test_cls[[i]])
    ##pred.nb[[i]] <- predict(fit.nb[[i]],test_cls[[i]])
    #pred.p[[i]]      <- predict(fit.p[[i]],test_cls[[i]],type="response",se.fit=TRUE)
    #pred.qp[[i]]     <- predict(fit.qp[[i]],test_cls[[i]],type="response")
    pred.nb[[i]]     <- predict(fit.nb[[i]],test_cls[[i]],type="response")
    #pred.hurdle[[i]] <- predict(fit.hurdle[[i]],test_cls[[i]],type="response")
}

#Compare the prediction and observations of the mean and sd
##mean.pred <- vector("list",length=ncls)
##mean.obs <- vector("list",length=ncls)
#mean.pred.p  <- vector("numeric",length=ncls)
#sd.pred.p  <- vector("numeric",length=ncls)
mean.pred.nb <- vector("numeric",length=ncls)
sd.pred.nb <- vector("numeric",length=ncls)
mean.obs     <- vector("numeric",length=ncls)
sd.obs     <- vector("numeric",length=ncls)

for (i in 1:10){
    ##mean.pred[[i]] <- mean(pred.nb[[i]])
    ##mean.obs[[i]] <- mean(test_cls[[i]])
    #mean.pred.p[i]  <- mean(pred.p[[i]])
    #sd.pred.p[i]    <- sd(pred.p[[i]])
    mean.pred.nb[i] <- mean(pred.nb[[i]])
    sd.pred.nb[i]   <- sd(pred.nb[[i]])
    mean.obs[i] <- mean(test_cls[[i]]$LGHT)
    sd.obs[i]   <- sd(test_cls[[i]]$LGHT)
}

coefs <- vector("list",length=ncls)
confints <- vector("list",length=ncls)

for (i in 1:ncls){
    confints[[i]] <- confint(fit.nb[[i]])

}

for (i in 1:ncls){
    coefs[[i]] <- fit.nb[[i]]$coefficients
}

mcoefs <- vector("list",length=ncls)
for (i in 1:ncls){
    mcoefs[[i]]<-data.frame(as.data.frame(coefs[[i]]),
                            as.data.frame(confints[[i]]),
                            CLS=i,
                            row.names=c("Intercept",
                                        "CAPE","LR75",
                                        "H5_T","S06",
                                        "PMR"))
    mcoefs[[i]]$coef_id <- rownames(mcoefs[[i]])
}

##Now collapse it all to a data frame
model_coefs <- do.call("rbind",mcoefs)
rownames(model_coefs) <- NULL
names(model_coefs) <- c("Coeff","lower","upper","cls","var")
model_coefs$cls <- as.factor(model_coefs$cls)
model_coefs$var <- as.factor(model_coefs$var)
#library(reshape2)
#model_coefs.melt <- melt(model_coefs,id.vars=c("cls","var"))
#
##Plot some information about the cluster coefficients
#plot(1,fit.nb[[1]]$coefficients[["CAPE"]],
#     xlim=c(1,10),ylim=c(0,2))
#for (i in 2:ncls){
#    points(i,fit.nb[[i]]$coefficients[["CAPE"]])
#}
#

# *****
#Get the names right here
#pdf(paste0("Model_coefficients_",ncls,"_clusters.pdf"),
#    width=10,height=10,paper="special")
pdf(paste0("Model_coefficients_seg_cape_",ncls,"_clusters.pdf"),
    width=10,height=10,paper="special")
#See http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/
#for hints on how to plot the error bars

    coef.plot <- ggplot(model_coefs,aes(x=cls,y=Coeff))+
                     facet_grid(var~.,scales="free")+
                     geom_point()+
                     geom_errorbar(aes(ymin=lower,ymax=upper),width=0.1)+
                     xlab("Cluster")+ylab("Coefficient")
    print(coef.plot)
dev.off()


#Evaluate the proportion of lightning in each cluster using the 
#regression data. This is only land points!
#Whereas using hlprob$prop.lightning has ocean samples too.

land.prop.lightning <- vector("numeric",length=ncls)

for (i in 1:ncls){
    land.prop.lightning[i] <- length(which(reg_data$CLS == i & 
                                           reg_data$LGHT > 0))/
                              length(which(reg_data$CLS == i))
}



##Plot some information about the lightning probability
#hlprob$prop.lightning.s <- hlprob$prop.lightning*100. #As a %
#hlprob$land.prop.lightning <- land.prop.lightning
#hlprob$land.prop.lightning.s <- land.prop.lightning*100.
#hd.l <- melt(hlprob,id.vars="cl")
#
##ggplot(subset(hd.l,variable %in% c("prop.lightning.s","prop.lightning")),
##ggplot(subset(hd.l,variable %in% c("prop.lightning.s")),
#ggplot(subset(hd.l,variable %in% c("prop.lightning","land.prop.lightning")),
#       aes(x=cl,y=value,fill=variable))+
#       #aes(x=cl,y=value,fill="prop.lightning"))+
#       geom_bar(stat="identity",position="dodge")+
#       xlab("Cluster")+ylab("Proportion of cluster that produces lightning")+
#       guides(fill=guide_legend(title=NULL))+
#       #scale_fill_discrete(breaks=c("prop.lightning.s"),
#       #scale_fill_discrete(breaks=c("prop.lightning"),
#       #                    labels=c("Lightning"))+
#       scale_fill_discrete(breaks=c("prop.lightning","land.prop.lightning"),
#                           labels=c("Oceanic and continental","Continental"))+
#       theme(axis.title=element_text(size=12),
#             axis.text=element_text(size=12),
#             legend.text=element_text(size=12))+
#       scale_y_continuous(labels=percent)
#
##phd <- ggplot(hd.l,aes(cl,value,fill=variable))+
##       geom_bar(stat="identity",position="dodge")+
##       xlab("Cluster")+ylab("%")+
##        scale_fill_brewer("",
##                          labels=c("Hail","No hail",
##                                   "Lightning","No lightning","Lightning strikes"),
##                                   palette="Set1")
#



###Dont include the Indonesian locations
##This excludes one lat point at the top of the domain
#min_lat <- min(lat_gpats)
##max_lat <- max(lat_gpats)
#max_lat <- -10.5
#min_lon <- min(lon_gpats)
#max_lon <- max(lon_gpats)
##
#ex.lat_gpats <- which(lat_gpats >= min_lat & lat_gpats <= max_lat)
#ex.lon_gpats <- which(lon_gpats >= min_lon & lon_gpats <= max_lon)
#lat_gpats <- lat_gpats[ex.lat_gpats]
#lon_gpats <- lon_gpats[ex.lon_gpats]
#total_lightning <- total_lightning[ex.lon_gpats,ex.lat_gpats,]
##cls.pred   <- cls.pred[ex.lon_era,ex.lat_era,]
##
#ex.lat_lsmask <- which(lat_lsmask >= min_lat & lat_lsmask <= max_lat)
#ex.lon_lsmask <- which(lon_lsmask >= min_lon & lon_lsmask <= max_lon)
#lat_lsmask <- lat_lsmask[ex.lat_lsmask]
#lon_lsmask <- lon_lsmask[ex.lon_lsmask]
#lsmask <- lsmask[ex.lon_lsmask,ex.lat_lsmask]
##


###For ease of plotting
lat  <- lat_lsmask
lon  <- lon_lsmask
nlat <- length(lat)
nlon <- length(lon)
#
ntime_cls <- dim(cls.pred)[3]
ndays_cls <- ntime_cls/4.
nyrs_cls <- ndays_cls/365.25

ndays_gpats <- difftime(rtime_gpats[dim(total_lightning)[3]],
                        rtime_gpats[1])

#Now read the ERA-data, so that we can perform a GLM regression
#of lightning against each of the ERA variables.

#An array to hold the cluster numbers
cluster <- cls.pred
#An array to hold the predicted number of stirkes
strikes <- vector("list",length=ncls)
for (i in 1:ncls){
    #strikes[[i]] <- array(NA,dim=c(nlon,nlat,ntime_cls))
    strikes[[i]] <- array(0,dim=c(nlon,nlat,ntime_cls))
}

strk_inds <- vector("list",length=ncls)
strk_df <- vector("list",length=ncls)
strk_prd <- vector("list",length=ncls)
sum_strikes <- vector("list",length=ncls)

cat("Calculating predicted number of strikes \n")
for (i in 1:ncls){
    strk_inds[[i]] <- which(cluster == i)
    strk_df[[i]]   <- data.frame(CAPE=mucape_n[strk_inds[[i]]],
                                 LR75=lr75_n[strk_inds[[i]]],
                                 H5_T=h5_temp_n[strk_inds[[i]]],
                                 S06=s06_n[strk_inds[[i]]],
                                 PMR=pmr_n[strk_inds[[i]]])
    #strk_prd[[i]] <- predict(fit.p[[i]],strk_df[[i]],type="response")
    strk_prd[[i]] <- predict(fit.nb[[i]],strk_df[[i]],type="response")
    strikes[[i]][strk_inds[[i]]] <- strk_prd[[i]]*land.prop.lightning[i]
    sum_strikes[[i]] <- apply(strikes[[i]],c(1,2),sum,na.rm=TRUE)
    cat("Finished cluster: ",i,"\n")
}

#This is the full predicted annual climatology
full_clim_pred_strikes <- Reduce("+",sum_strikes)
annual_clim_pred_strikes <- Reduce("+",strikes)


#Now plot the full annual climatology

##Transpose for plotting
#cluster[,,] <- t(apply(cls.pred[,,],1,rev))
##I used cluster_test to check that
##the above was performing the reversal correctly
##i.e., checking that reversing a 3-D array before summing
##was the same as summing and then reversing a 2-D array
##cluster_test <- as.numeric(as.character(cls.pred))
##dim(cluster_test) <- dim(cls.pred)
##cluster_test <- unclass(cls.pred)
#lsmask[,] <- t(apply(lsmask,1,rev))
##Don't transpose the total lightning! It takes too long
##Only transpose later once calculations have been applied
##total_lightning <- t(apply(total_lightning,1,rev))
##ncls <- length(levels(as.factor(cluster)))


#Plot the cluster distribution across Australia
cat("Evaluating cluster distribution \n")
clus.count <- c()
#clus.count_test <- c()
clus.dist <- vector("list",length=ncls)
ndays_clus.dist <- vector("list",length=ncls)
#clus.dist_test <- vector("list",length=ncls)
for (i in 1:ncls){
#for (i in 10){
    clstemp <- cluster
#    clstemp_test <- cluster_test
    clus.count[i] <- length(which(clstemp==i))
#    clus.count_test[i] <- length(which(clstemp_test==i))
    #ifelse(clstemp != i,NA,1)
    #clstemp[clstemp != i] <- NA
    clstemp[clstemp != i] <- 0
    clstemp[clstemp == i] <- 1
    #ndays_clus,dist[[i]] <- apply(clstemp,c(1,2),roll_max,4,by=4,
    #                            align="left",na.rm=TRUE)
    #ndays_clus.dist[[i]] <- (apply(clstemp,c(1,2),sum,na.rm=T))/4.*
    #                        (prob.lightning[i]/100.)
    #clus.dist[[i]] <- ndays_clus.dist[[i]]/dim(cluster)[3]
    ndays_clus.dist[[i]] <- apply(clstemp,c(1,2),roll_max,4,by=4,
                                 align="left",na.rm=TRUE)
    ndays_clus.dist[[i]] <- ndays_clus.dist[[i]][
                             seq(1,dim(ndays_clus.dist[[i]])[1],4),,]
    ndays_clus.dist[[i]] <- aperm(ndays_clus.dist[[i]],c(2,3,1))
    ndays_clus.dist[[i]] <- apply(ndays_clus.dist[[i]],c(1,2),sum,na.rm=T)*
                              prob.lightning[i]/100.

    clus.dist[[i]] <- apply(clstemp,c(1,2),sum,na.rm=TRUE)/dim(cluster)[3]
    #clus.dist[[i]] <- apply(clstemp,c(1,2),mean,na.rm=TRUE)
    cat("Finished cluster: ",i,'\n')
}

#The number of days of lightning predicted by the
#cluster membership
full_clim_clus <- Reduce("+",ndays_clus.dist)

#Apply the land sea mask to the clus.dist
#And convert to number of days of lightning
for (i in 1:ncls){
    ndays_clus.dist[[i]][lsmask==0] <- NA
    clus.dist[[i]][lsmask==0] <- NA
    full_clim_clus[lsmask==0] <- NA
    ndays_lightning[lsmask==0] <- NA
}

#Some statistics to put on the plot
ndays_light_mn <- mean(ndays_lightning/nyrs_cls,na.rm=T)
ndays_light_sd <- sd(ndays_lightning/nyrs_cls,na.rm=T)
full_clim_clus_mn <- mean(full_clim_clus/nyrs_cls,na.rm=T)
full_clim_clus_sd <- sd(full_clim_clus/nyrs_cls,na.rm=T)


#Now plot the number of days of lightning predicted by cluster
#and the number of days of GPATS

# *****
# Get the names right
#pdf(paste0("Ndays_lightning_comparison",ncls,"_clusters.pdf"),
#    width=10,height=15,paper="special")
pdf(paste0("Ndays_lightning_comparison_seg_cape",ncls,"_clusters.pdf"),
    width=10,height=15,paper="special")
#par(mar=c(3,3,1,1),cex=1.25)
#par(cex=1.25)
split.screen( rbind(c(0, .85,0,1), c(.85,1,0,1)))
split.screen(c(2,1), screen=1)-> ind
screen(ind[1])
par(mar=c(0,0,0,0),cex=1.25)
par(oma=c(0,0,0,0))
zmax <- max(c(full_clim_clus/nyrs_cls,ndays_lightning/nyrs_cls),na.rm=T)
image(lon,rev(lat),t(apply(ndays_lightning/nyrs_cls,1,rev)),xlim=c(112,155),ylim=c(-45,-10),
           zlim=c(0,zmax),
           #xlab="Lon",ylab="Lat",
           xlab="",ylab="",
           xaxt="n",yaxt="n",
           col=colorRampPalette(rev(brewer.pal(9,"Spectral")))(64))
text(115,-12,"(a)",cex=2)
#text(135,-12,paste0(ndays_light_mn,%+-%,ndays_light_sd))
#text(135,-12,expression(paste0(round(ndays_light_mn,2),%+-%,round(ndays_light_sd,2))))
oz(lwd=2.5,add=TRUE)

screen(ind[2])
par(mar=c(0,0,0,0),cex=1.25)
par(oma=c(0,0,0,0))
image(lon,rev(lat),t(apply(full_clim_clus/nyrs_cls,1,rev)),xlim=c(112,155),ylim=c(-45,-10),
           #main="GPATS observations of the average hourly flash count",
           zlim=c(0,zmax),
           #xlab="Lon",ylab="Lat",
           xlab="",ylab="",
           xaxt="n",yaxt="n",
           col=colorRampPalette(rev(brewer.pal(9,"Spectral")))(64))
text(115,-12,"(b)",cex=2)
#text(135,-12,expression(full_clim_clus_mn %+-% full_clim_clus_sd))
oz(lwd=2.5,add=TRUE)

screen(2)
par(mar=c(0,0,0,0),cex=1.25)
par(oma=c(0,0,0,0))
image.plot(zlim=c(0,zmax),
           #legend.only=TRUE,smallplot=c(.1,.2, .2,.8),
           legend.only=TRUE,smallplot=c(.1,.3, .2,.8),
           col=colorRampPalette(rev(brewer.pal(9,"Spectral")))(64),
           #legend.lab="Days per year",legend.cex=1.25)
           legend.lab="Days per year")
dev.off()

#Now a plot of the differences to get biases
#bias_lght <- (ndays_lightning-full_clim_clus)/nyrs_cls
bias_lght <- (full_clim_clus-ndays_lightning)/nyrs_cls
col_levs <- 50

#Difference as a fraction of the observations
prop <- bias_lght/(ndays_lightning/nyrs_cls)
prop[prop > 1] <- NA
prop[prop < -1] <- NA

# *****
#Get the names right!
#pdf(paste0("Difference_gpats_cluster_lightning_days_",ncls,"_clusters.pdf"),
#    width=10,height=15,paper="special")
pdf(paste0("Difference_gpats_cluster_lightning_days_seg_cape_",ncls,"_clusters.pdf"),
    width=10,height=15,paper="special")
par(mar=c(3,3,1,1))
par(cex=1.25)
par(mfrow=c(2,1))

max_abs_val <- max(abs(c(min(bias_lght,na.rm=T),max(bias_lght,na.rm=T))))
col_seq <- seq(-1.*max_abs_val,max_abs_val,length.out=col_levs+1)
image.plot(lon,rev(lat),t(apply(bias_lght,1,rev)),
           xlim=c(110,160),ylim=c(-45,-10),
           #main="Observations - model",
           xlab="Lon",ylab="Lat",
           col=colorRampPalette(rev(brewer.pal(9,"RdBu")))(col_levs),
           breaks=col_seq)
           #legend.lab="Lightning counts per day",horizontal=TRUE)
           #legend.lab="Lightning counts per hour",horizontal=TRUE)
text(115,-12,"(a)",cex=2)
oz(lwd=2.5,add=T)

max_abs_val <- max(abs(c(min(prop,na.rm=T),max(prop,na.rm=T))))
col_seq <- seq(-1.*max_abs_val,max_abs_val,length.out=col_levs+1)
image.plot(lon,rev(lat),t(apply(prop,1,rev)),
           xlim=c(110,160),ylim=c(-45,-10),
           #main="Proportionate difference between cluster and GPATS",
           xlab="Lon",ylab="Lat",
           #zlim=c(-1.,1),
           col=colorRampPalette(rev(brewer.pal(9,"RdBu")))(col_levs),
           breaks=col_seq)
text(115,-12,"(b)",cex=2)
oz(lwd=2.5,add=T)
dev.off()

max_abs_val <- max(abs(c(min(prop,na.rm=T),max(prop,na.rm=T))))
col_seq <- seq(-1.*max_abs_val,max_abs_val,length.out=col_levs+1)
pdf(paste0("Proportionate_Difference_gpats_cluster_lightning_days_",ncls,"_clusters.pdf"),
    width=10,height=10,paper="special")
image.plot(lon,rev(lat),t(apply(prop,1,rev)),xlim=c(110,160),ylim=c(-45,-10),
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

#Evaluate the precentage of each cluster
clus.perc <- clus.count*100./sum(clus.count)

#ncols=50 #No of colors for plotting
ncols=15 #No of colors for plotting
#Evaluate the limits for the spatial plots
clus.dist.max <- c()
strike.dist.max <- c()
for (i in 1:ncls){
    clus.dist.max[i] <- max(clus.dist[[i]],na.rm=T)
    strike.dist.max[i] <- max(sum_strikes[[i]],na.rm=T)
}
cls.max <- round(max(clus.dist.max),digits=2)
strike.max <- round(max(strike.dist.max),digits=2)

#The spatial distribution of the clusters with info. on their frequency
# *****
# Get the names right!
#pdf(paste0("Lightning_Cluster_spatial_distribution_",ncls,"_clusters.pdf"),
#    width=10,height=10,paper="special")
pdf(paste0("Lightning_Cluster_spatial_distribution_seg_cape_",ncls,"_clusters.pdf"),
    width=10,height=10,paper="special")
#par(mfrow=c(3,ncls/3))
par(mfrow=c(2,2))#20 clusters
#par(mar=c(4,5,3,2,1)+0.1)
#par(mar=c(4,4,3,2)+0.1)
par(mar=c(4,2,1,1)+0.1)
for (i in 1:ncls){
    #image.plot(lon,rev(lat),clus.dist[[i]],
    image.plot(lon,rev(lat),t(apply(clus.dist[[i]],1,rev)),
               #main=paste0("Cluster: ",i),
               xlab="",ylab="",
               #xlab="Lon",ylab="Lat",
               #zlim=c(0,0.6),
               col=colorRampPalette(rev(brewer.pal(9,"Spectral")))(ncols),
               horizontal=TRUE,zlim=c(0,cls.max))
    text(115,-15,paste0(i),cex=1.5,font=2)
    text(155,-15,paste(round(clus.perc[i],1),"%"),cex=1.25)
    oz(lwd=2.5,add=T)
}
dev.off()

#The spatial distribution of the lightning strikes
pdf(paste0("Lightning_strike_spatial_distribution_",ncls,"_clusters.pdf"),
    width=10,height=10,paper="special")
#par(mfrow=c(3,ncls/3))
par(mfrow=c(2,2))#20 clusters
#par(mar=c(4,5,3,2,1)+0.1)
#par(mar=c(4,4,3,2)+0.1)
par(mar=c(4,2,1,1)+0.1)
for (i in 1:ncls){
    #image.plot(lon,rev(lat),clus.dist[[i]],
    image.plot(lon,rev(lat),t(apply(sum_strikes[[i]],1,rev)),
               #main=paste0("Cluster: ",i),
               xlab="",ylab="",
               #xlab="Lon",ylab="Lat",
               #zlim=c(0,0.6),
               col=colorRampPalette(rev(brewer.pal(9,"Spectral")))(ncols),
               #horizontal=TRUE,zlim=c(0,strike.max))
               horizontal=TRUE)
    text(115,-15,paste0(i),cex=1.5,font=2)
    #text(155,-15,paste(round(clus.perc[i],1),"%"),cex=1.25)
    oz(lwd=2.5,add=T)
}
dev.off()


#Now plot the seasonal climatology of the lightning prediction
lightning.per.day <- array(dim=c(nlon,nlat,ntime_cls))
lightning.per.day[,,] <- prob.lightning[cluster[,,]]/100.
#Get the cluster data in a daily average
#First we need to extract the GPATS time which correspond to the ERA times
library(lubridate)
#rtime <- rtime_gpats[seq(1,ntime_gpats,6)]
#rtime_cls <- rtime_gpats[seq(1,ntime_gpats,6)]
rtime <- rtime_cls[1:mxt]

#lpd0 <- apply(lightning.per.day,c(1,2),roll_mean,4,by=4,align="left",fill=NA)
lpd0 <- apply(lightning.per.day,c(1,2),roll_sum,4,by=4,align="left",fill=NA)

lpd0 <- lpd0[seq(1,dim(lpd0)[1],4),,]
lpd0 <- aperm(lpd0,c(2,3,1))

#Now sum and average them over the year
#lpd <- apply(lpd0,c(2,3),sum,na.rm=T)/365.
#lpd <- apply(lpd0,c(2,3),sum,na.rm=T)/ntime_cls
#lpd <- apply(lpd0,c(2,3),sum,na.rm=T)/ndays_cls


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
#lpd0.summer <- apply(lightning.per.day[,,yr.2008[summer]],
#                     c(1,2),roll_max,4,by=4,align="left",fill=NA)
#lpd0.autumn <- apply(lightning.per.day[,,yr.2008[autumn]],
#                     c(1,2),roll_max,4,by=4,align="left",fill=NA)
#lpd0.winter <- apply(lightning.per.day[,,yr.2008[winter]],
#                     c(1,2),roll_max,4,by=4,align="left",fill=NA)
#lpd0.spring <- apply(lightning.per.day[,,yr.2008[spring]],
#                     c(1,2),roll_max,4,by=4,align="left",fill=NA)

lpd0.summer <- apply(lightning.per.day[,,summer],
                     c(1,2),roll_sum,4,by=4,align="left",fill=NA)
lpd0.autumn <- apply(lightning.per.day[,,autumn],
                     c(1,2),roll_sum,4,by=4,align="left",fill=NA)
lpd0.winter <- apply(lightning.per.day[,,winter],
                     c(1,2),roll_sum,4,by=4,align="left",fill=NA)
lpd0.spring <- apply(lightning.per.day[,,spring],
                     c(1,2),roll_sum,4,by=4,align="left",fill=NA)

#We need to extract every fourth element as the roll_sum routine leaves it in
#there
lpd0.summer <- lpd0.summer[seq(1,dim(lpd0.summer)[1],4),,]
lpd0.autumn <- lpd0.autumn[seq(1,dim(lpd0.autumn)[1],4),,]
lpd0.winter <- lpd0.winter[seq(1,dim(lpd0.winter)[1],4),,]
lpd0.spring <- lpd0.spring[seq(1,dim(lpd0.spring)[1],4),,]

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

# *****
#pdf(paste0("Cluster_seasonal_cycle",ncls,"_cluster.pdf"),
#    width=10,height=10,paper="special")
pdf(paste0("Cluster_seasonal_cycle_seg_cape",ncls,"_cluster.pdf"),
    width=10,height=10,paper="special")
image.plot(lon,rev(lat),t(apply(lpd.summer,1,rev)),xlim=c(110,160),ylim=c(-45,-10),
        main="Cluster seasonal cycle",
        xlab="Lon",ylab="Lat",
        col=colorRampPalette(rev(brewer.pal(9,"Spectral")))(ncols))
oz(lwd=2.5,add=TRUE)

image.plot(lon,rev(lat),t(apply(lpd.autumn,1,rev)),xlim=c(110,160),ylim=c(-45,-10),
        main="Cluster seasonal cycle",
        xlab="Lon",ylab="Lat",
        col=colorRampPalette(rev(brewer.pal(9,"Spectral")))(ncols))
oz(lwd=2.5,add=TRUE)

image.plot(lon,rev(lat),t(apply(lpd.winter,1,rev)),xlim=c(110,160),ylim=c(-45,-10),
        main="Cluster seasonal cycle",
        xlab="Lon",ylab="Lat",
        col=colorRampPalette(rev(brewer.pal(9,"Spectral")))(ncols))
oz(lwd=2.5,add=TRUE)

image.plot(lon,rev(lat),t(apply(lpd.spring,1,rev)),xlim=c(110,160),ylim=c(-45,-10),
        main="Cluster seasonal cycle",
        xlab="Lon",ylab="Lat",
        col=colorRampPalette(rev(brewer.pal(9,"Spectral")))(ncols))
oz(lwd=2.5,add=TRUE)

dev.off()

##########----------##########

#Now do the same for the lightning_strikes

#The observations
spd0.summer <- apply(lightning_matrix[,,summer],
                     c(1,2),roll_sum,4,by=4,align="left",fill=NA)
spd0.autumn <- apply(lightning_matrix[,,autumn],
                     c(1,2),roll_sum,4,by=4,align="left",fill=NA)
spd0.winter <- apply(lightning_matrix[,,winter],
                     c(1,2),roll_sum,4,by=4,align="left",fill=NA)
spd0.spring <- apply(lightning_matrix[,,spring],
                     c(1,2),roll_sum,4,by=4,align="left",fill=NA)

#The prediction
spd1.summer <- apply(annual_clim_pred_strikes[,,summer],
                     c(1,2),roll_sum,4,by=4,align="left",fill=NA)
spd1.autumn <- apply(annual_clim_pred_strikes[,,autumn],
                     c(1,2),roll_sum,4,by=4,align="left",fill=NA)
spd1.winter <- apply(annual_clim_pred_strikes[,,winter],
                     c(1,2),roll_sum,4,by=4,align="left",fill=NA)
spd1.spring <- apply(annual_clim_pred_strikes[,,spring],
                     c(1,2),roll_sum,4,by=4,align="left",fill=NA)


#We need to extract every fourth element as the roll_sum routine leaves it in.
spd0.summer <- spd0.summer[seq(1,dim(spd0.summer)[1],4),,]
spd0.autumn <- spd0.autumn[seq(1,dim(spd0.autumn)[1],4),,]
spd0.winter <- spd0.winter[seq(1,dim(spd0.winter)[1],4),,]
spd0.spring <- spd0.spring[seq(1,dim(spd0.spring)[1],4),,]

spd1.summer <- spd1.summer[seq(1,dim(spd1.summer)[1],4),,]
spd1.autumn <- spd1.autumn[seq(1,dim(spd1.autumn)[1],4),,]
spd1.winter <- spd1.winter[seq(1,dim(spd1.winter)[1],4),,]
spd1.spring <- spd1.spring[seq(1,dim(spd1.spring)[1],4),,]

#Now sum and average them over the year
#lpd <- apply(lpd0,c(2,3),sum,na.rm=T)/365.
#lpd <- apply(lpd0,c(2,3),sum,na.rm=T)/ntime_cls
spd_gpats.summer <- apply(spd0.summer,c(2,3),sum,na.rm=T)/ndays_cls
spd_gpats.autumn <- apply(spd0.autumn,c(2,3),sum,na.rm=T)/ndays_cls
spd_gpats.winter <- apply(spd0.winter,c(2,3),sum,na.rm=T)/ndays_cls
spd_gpats.spring <- apply(spd0.spring,c(2,3),sum,na.rm=T)/ndays_cls

spd_pred.summer <- apply(spd1.summer,c(2,3),sum,na.rm=T)/ndays_cls
spd_pred.autumn <- apply(spd1.autumn,c(2,3),sum,na.rm=T)/ndays_cls
spd_pred.winter <- apply(spd1.winter,c(2,3),sum,na.rm=T)/ndays_cls
spd_pred.spring <- apply(spd1.spring,c(2,3),sum,na.rm=T)/ndays_cls

spd_gpats.summer[lsmask==0] <- NA
spd_gpats.autumn[lsmask==0] <- NA
spd_gpats.winter[lsmask==0] <- NA
spd_gpats.spring[lsmask==0] <- NA

spd_pred.summer[lsmask==0] <- NA
spd_pred.autumn[lsmask==0] <- NA
spd_pred.winter[lsmask==0] <- NA
spd_pred.spring[lsmask==0] <- NA

cat("Plotting seasonal cycle \n")
pdf(paste0("Lightning_seasonal_strike_comparison",ncls,"_clusters.pdf"),
    width=10,height=15,paper="special")
    #par(mar=c(3,3,1,1))
    #par(mar=c(0,0,0,0))
    #par(cex=1.25)
    gpar <- list(mar=c(0,0,0,0),cex=1.25,oma=c(0,0,0,0))
    split.screen( rbind(c(0, .9,0,1), c(.9,1,0,1)))
    split.screen(c(4,2), screen=1)-> ind
    #split.screen(c(4,3), screen=1)-> ind
#Summer
    screen(ind[1])
    par(gpar)
    zmax.summer <- max(c(spd_gpats.summer,spd_pred.summer),na.rm=T)
    image(lon,rev(lat),t(apply(spd_gpats.summer,1,rev)),xlim=c(110,155),ylim=c(-45,-10),
               zlim=c(0,zmax.summer),
               #xlab="Lon",ylab="Lat",
               xlab="",ylab="",
               xaxt="n",yaxt="n",
               col=colorRampPalette(rev(brewer.pal(9,"Spectral")))(64))
    text(115,-12,"DJF",cex=2.0)
    oz(lwd=2.5,add=TRUE)
    
    screen(ind[2])
    par(gpar)
    image(lon,rev(lat),t(apply(spd_pred.summer,1,rev)),xlim=c(110,155),ylim=c(-45,-10),
               #main="GPATS observations of the average hourly flash count",
               zlim=c(0,zmax.summer),
               #xlab="Lon",ylab="Lat",
               xlab="",ylab="",
               xaxt="n",yaxt="n",
               col=colorRampPalette(rev(brewer.pal(9,"Spectral")))(64))
    #text(115,-12,"(b)",cex=2)
    oz(lwd=2.5,add=TRUE)


#Autumn
    screen(ind[3])
    par(gpar)
    zmax.autumn <- max(c(spd_gpats.autumn,spd_pred.autumn),na.rm=T)
    image(lon,rev(lat),t(apply(spd_gpats.autumn,1,rev)),xlim=c(110,155),ylim=c(-45,-10),
               zlim=c(0,zmax.autumn),
               xlab="",ylab="",
               xaxt="n",yaxt="n",
               col=colorRampPalette(rev(brewer.pal(9,"Spectral")))(64))
    text(115,-12,"MAM",cex=2)
    oz(lwd=2.5,add=TRUE)
    
    screen(ind[4])
    par(gpar)
    image(lon,rev(lat),t(apply(spd_pred.autumn,1,rev)),xlim=c(110,155),ylim=c(-45,-10),
               #main="GPATS observations of the average hourly flash count",
               zlim=c(0,zmax.autumn),
               xlab="",ylab="",
               xaxt="n",yaxt="n",
               col=colorRampPalette(rev(brewer.pal(9,"Spectral")))(64))
    #text(115,-12,"(b)",cex=2)
    oz(lwd=2.5,add=TRUE)

#Winter
    screen(ind[5])
    par(gpar)
    zmax.winter <- max(c(spd_gpats.winter,spd_pred.winter),na.rm=T)
    image(lon,rev(lat),t(apply(spd_gpats.winter,1,rev)),xlim=c(110,155),ylim=c(-45,-10),
               zlim=c(0,zmax.winter),
               xlab="",ylab="",
               xaxt="n",yaxt="n",
               col=colorRampPalette(rev(brewer.pal(9,"Spectral")))(64))
    text(115,-12,"JJA",cex=2)
    oz(lwd=2.5,add=TRUE)
    
    screen(ind[6])
    par(gpar)
    image(lon,rev(lat),t(apply(spd_pred.winter,1,rev)),xlim=c(110,155),ylim=c(-45,-10),
               #main="GPATS observations of the average hourly flash count",
               zlim=c(0,zmax.winter),
               xlab="",ylab="",
               xaxt="n",yaxt="n",
               col=colorRampPalette(rev(brewer.pal(9,"Spectral")))(64))
    #text(115,-12,"(b)",cex=2)
    oz(lwd=2.5,add=TRUE)


#Spring
    screen(ind[7])
    par(gpar)
    zmax.spring <- max(c(spd_gpats.spring,spd_pred.spring),na.rm=T)
    image(lon,rev(lat),t(apply(spd_gpats.spring,1,rev)),xlim=c(110,155),ylim=c(-45,-10),
               zlim=c(0,zmax.spring),
               xlab="",ylab="",
               xaxt="n",yaxt="n",
               col=colorRampPalette(rev(brewer.pal(9,"Spectral")))(64))
    text(115,-12,"SON",cex=2)
    oz(lwd=2.5,add=TRUE)
    
    screen(ind[8])
    par(gpar)
    image(lon,rev(lat),t(apply(spd_pred.spring,1,rev)),xlim=c(110,155),ylim=c(-45,-10),
               #main="GPATS observations of the average hourly flash count",
               zlim=c(0,zmax.spring),
               xlab="",ylab="",
               xaxt="n",yaxt="n",
               col=colorRampPalette(rev(brewer.pal(9,"Spectral")))(64))
    #text(115,-12,"(b)",cex=2)
    oz(lwd=2.5,add=TRUE)


    split.screen(c(4,1), screen=2)-> ind
    screen(ind[1])
    par(gpar)
    image.plot(zlim=c(0,zmax.summer),
               #legend.only=TRUE,smallplot=c(.1,.2, .2,.8),
               legend.only=TRUE,smallplot=c(.1,.3, .2,.8),
               col=colorRampPalette(rev(brewer.pal(9,"Spectral")))(64))
               #legend.lab="Days per year",legend.cex=1.25)

    screen(ind[2])
    par(gpar)
    image.plot(zlim=c(0,zmax.autumn),
               #legend.only=TRUE,smallplot=c(.1,.2, .2,.8),
               legend.only=TRUE,smallplot=c(.1,.3, .2,.8),
               col=colorRampPalette(rev(brewer.pal(9,"Spectral")))(64))
               #legend.lab="Days per year",legend.cex=1.25)

    screen(ind[3])
    par(gpar)
    image.plot(zlim=c(0,zmax.winter),
               #legend.only=TRUE,smallplot=c(.1,.2, .2,.8),
               legend.only=TRUE,smallplot=c(.1,.3, .2,.8),
               col=colorRampPalette(rev(brewer.pal(9,"Spectral")))(64))
               #legend.lab="Days per year",legend.cex=1.25)

    screen(ind[4])
    par(gpar)
    image.plot(zlim=c(0,zmax.spring),
               #legend.only=TRUE,smallplot=c(.1,.2, .2,.8),
               legend.only=TRUE,smallplot=c(.1,.3, .2,.8),
               col=colorRampPalette(rev(brewer.pal(9,"Spectral")))(64))
               #legend.lab="Days per year",legend.cex=1.25)

dev.off()

#Plot the difference between the model and the observations

diff.summer <- t(apply(spd_pred.summer,1,rev))-
               t(apply(spd_gpats.summer,1,rev))
diff.autumn <- t(apply(spd_pred.autumn,1,rev))-
               t(apply(spd_gpats.autumn,1,rev))
diff.winter <- t(apply(spd_pred.winter,1,rev))-
               t(apply(spd_gpats.winter,1,rev))
diff.spring <- t(apply(spd_pred.spring,1,rev))-
               t(apply(spd_gpats.spring,1,rev))

prop.diff.summer <- diff.summer/t(apply(spd_gpats.summer,1,rev))
prop.diff.autumn <- diff.autumn/t(apply(spd_gpats.autumn,1,rev))
prop.diff.winter <- diff.winter/t(apply(spd_gpats.winter,1,rev))
prop.diff.spring <- diff.spring/t(apply(spd_gpats.spring,1,rev))

col_levs <- 50
prop.diff.summer[prop.diff.summer > 1] <- NA
prop.diff.autumn[prop.diff.autumn > 1] <- NA
prop.diff.winter[prop.diff.winter > 1] <- NA
prop.diff.spring[prop.diff.spring > 1] <- NA

pdf(paste0("Lightning_seasonal_strike_difference",ncls,"_clusters.pdf"),
    width=10,height=10,paper="special")
    par(mar=c(0,0,0,0)+0.1,bty="n",cex=1.25,oma=c(0,0,0,0.5),mfrow=c(2,2))
    #par(mfrow=c(2,2))
    #par(mar=c(0,0,0,3)+0.1)
    #par(bty="n")
#Summer
    max_abs.summer <- max(abs(c(min(diff.summer,na.rm=T),
                                max(diff.summer,na.rm=T))))
    col_seq <- seq(-1.*max_abs.summer,max_abs.summer,length.out=col_levs+1)
    image.plot(lon,rev(lat),diff.summer,xlim=c(110,155),ylim=c(-45,-10),
               #zlim=c(0,zmax.summer),
               #xlab="Lon",ylab="Lat",
               xlab="",ylab="",
               xaxt="n",yaxt="n",
               col=colorRampPalette(rev(brewer.pal(9,"RdBu")))(col_levs),
               breaks=col_seq,legend.mar=5.5,legend.width=1.5)
    text(115,-12,"DJF",cex=2.0)
    oz(lwd=2.5,add=TRUE)
    
#Autumn
    max_abs.autumn <- max(abs(c(min(diff.autumn,na.rm=T),
                                max(diff.autumn,na.rm=T))))
    col_seq <- seq(-1.*max_abs.autumn,max_abs.autumn,length.out=col_levs+1)
    image.plot(lon,rev(lat),diff.autumn,xlim=c(110,155),ylim=c(-45,-10),
               #zlim=c(0,zmax.summer),
               #xlab="Lon",ylab="Lat",
               xlab="",ylab="",
               xaxt="n",yaxt="n",
               col=colorRampPalette(rev(brewer.pal(9,"RdBu")))(col_levs),
               breaks=col_seq,legend.mar=5.5,legend.width=1.5)
    text(115,-12,"MAM",cex=2.0)
    oz(lwd=2.5,add=TRUE)

#Winter
    max_abs.winter <- max(abs(c(min(diff.winter,na.rm=T),
                                max(diff.winter,na.rm=T))))
    col_seq <- seq(-1.*max_abs.winter,max_abs.winter,length.out=col_levs+1)
    image.plot(lon,rev(lat),diff.winter,xlim=c(110,155),ylim=c(-45,-10),
               #zlim=c(0,zmax.summer),
               #xlab="Lon",ylab="Lat",
               xlab="",ylab="",
               xaxt="n",yaxt="n",
               col=colorRampPalette(rev(brewer.pal(9,"RdBu")))(col_levs),
               breaks=col_seq,legend.mar=5.5,legend.width=1.5)
    text(115,-12,"JJA",cex=2.0)
    oz(lwd=2.5,add=TRUE)

#Spring
    max_abs.spring <- max(abs(c(min(diff.spring,na.rm=T),
                                max(diff.spring,na.rm=T))))
    col_seq <- seq(-1.*max_abs.spring,max_abs.spring,length.out=col_levs+1)
    image.plot(lon,rev(lat),diff.spring,xlim=c(110,155),ylim=c(-45,-10),
               #zlim=c(0,zmax.summer),
               #xlab="Lon",ylab="Lat",
               xlab="",ylab="",
               xaxt="n",yaxt="n",
               col=colorRampPalette(rev(brewer.pal(9,"RdBu")))(col_levs),
               breaks=col_seq,legend.mar=5.5,legend.width=1.5)
    text(115,-12,"SON",cex=2.0)
    oz(lwd=2.5,add=TRUE)

dev.off()
 
    



#Now do the same for the diurnal cycle

era_t1 <- c(0)
era_t2 <- c(6)
era_t3 <- c(12)
era_t4 <- c(18)

t1 <- which(hour(rtime) %in% era_t1)
t2 <- which(hour(rtime) %in% era_t2)
t3 <- which(hour(rtime) %in% era_t3)
t4 <- which(hour(rtime) %in% era_t4)

#The observations
spd0.t1 <- apply(lightning_matrix[,,t1],
#spd0.t1 <- apply(lightning_at_era[,,t1],
                     c(1,2),roll_sum,4,by=4,align="left",fill=NA)
spd0.t2 <- apply(lightning_matrix[,,t2],
#spd0.t2 <- apply(lightning_at_era[,,t2],
                     c(1,2),roll_sum,4,by=4,align="left",fill=NA)
spd0.t3 <- apply(lightning_matrix[,,t3],
#spd0.t3 <- apply(lightning_at_era[,,t3],
                     c(1,2),roll_sum,4,by=4,align="left",fill=NA)
spd0.t4 <- apply(lightning_matrix[,,t4],
#spd0.t4 <- apply(lightning_at_era[,,t4],
                     c(1,2),roll_sum,4,by=4,align="left",fill=NA)

#The prediction
spd1.t1 <- apply(annual_clim_pred_strikes[,,t1],
                     c(1,2),roll_sum,4,by=4,align="left",fill=NA)
spd1.t2 <- apply(annual_clim_pred_strikes[,,t2],
                     c(1,2),roll_sum,4,by=4,align="left",fill=NA)
spd1.t3 <- apply(annual_clim_pred_strikes[,,t3],
                     c(1,2),roll_sum,4,by=4,align="left",fill=NA)
spd1.t4 <- apply(annual_clim_pred_strikes[,,t4],
                     c(1,2),roll_sum,4,by=4,align="left",fill=NA)


#We need to extract every fourth element as the roll_sum routine leaves it in.
spd0.t1 <- spd0.t1[seq(1,dim(spd0.t1)[1],4),,]
spd0.t2 <- spd0.t2[seq(1,dim(spd0.t2)[1],4),,]
spd0.t3 <- spd0.t3[seq(1,dim(spd0.t3)[1],4),,]
spd0.t4 <- spd0.t4[seq(1,dim(spd0.t4)[1],4),,]

spd1.t1 <- spd1.t1[seq(1,dim(spd1.t1)[1],4),,]
spd1.t2 <- spd1.t2[seq(1,dim(spd1.t2)[1],4),,]
spd1.t3 <- spd1.t3[seq(1,dim(spd1.t3)[1],4),,]
spd1.t4 <- spd1.t4[seq(1,dim(spd1.t4)[1],4),,]

#Now sum and average them over the year
#lpd <- apply(lpd0,c(2,3),sum,na.rm=T)/365.
#lpd <- apply(lpd0,c(2,3),sum,na.rm=T)/ntime_cls
spd_gpats.t1 <- apply(spd0.t1,c(2,3),sum,na.rm=T)/ndays_cls
spd_gpats.t2 <- apply(spd0.t2,c(2,3),sum,na.rm=T)/ndays_cls
spd_gpats.t3 <- apply(spd0.t3,c(2,3),sum,na.rm=T)/ndays_cls
spd_gpats.t4 <- apply(spd0.t4,c(2,3),sum,na.rm=T)/ndays_cls

spd_pred.t1 <- apply(spd1.t1,c(2,3),sum,na.rm=T)/ndays_cls
spd_pred.t2 <- apply(spd1.t2,c(2,3),sum,na.rm=T)/ndays_cls
spd_pred.t3 <- apply(spd1.t3,c(2,3),sum,na.rm=T)/ndays_cls
spd_pred.t4 <- apply(spd1.t4,c(2,3),sum,na.rm=T)/ndays_cls

spd_gpats.t1[lsmask==0] <- NA
spd_gpats.t2[lsmask==0] <- NA
spd_gpats.t3[lsmask==0] <- NA
spd_gpats.t4[lsmask==0] <- NA

spd_pred.t1[lsmask==0] <- NA
spd_pred.t2[lsmask==0] <- NA
spd_pred.t3[lsmask==0] <- NA
spd_pred.t4[lsmask==0] <- NA


cat("Plotting diurnal cycle \n")
pdf(paste0("Lightning_diurnal_strike_comparison",ncls,"_clusters.pdf"),
    width=10,height=15,paper="special")
    #par(mar=c(3,3,1,1))
    #par(mar=c(0,0,0,0))
    #par(cex=1.25)
    gpar <- list(mar=c(0,0,0,0),cex=1.25,oma=c(0,0,0,0))
    split.screen( rbind(c(0, .9,0,1), c(.9,1,0,1)))
    split.screen(c(4,2), screen=1)-> ind

#00 UTC
    screen(ind[1])
    par(gpar)
    zmax.t1 <- max(c(spd_gpats.t1,spd_pred.t1),na.rm=T)
    image(lon,rev(lat),t(apply(spd_gpats.t1,1,rev)),xlim=c(110,155),ylim=c(-45,-10),
               zlim=c(0,zmax.t1),
               #xlab="Lon",ylab="Lat",
               xlab="",ylab="",
               xaxt="n",yaxt="n",
               col=colorRampPalette(rev(brewer.pal(9,"Spectral")))(64))
    text(115,-12,"00 UTC",cex=1.5,pos=4)
    oz(lwd=2.5,add=TRUE)
    
    screen(ind[2])
    par(gpar)
    image(lon,rev(lat),t(apply(spd_pred.t1,1,rev)),xlim=c(110,155),ylim=c(-45,-10),
               #main="GPATS observations of the average hourly flash count",
               zlim=c(0,zmax.t1),
               #xlab="Lon",ylab="Lat",
               xlab="",ylab="",
               xaxt="n",yaxt="n",
               col=colorRampPalette(rev(brewer.pal(9,"Spectral")))(64))
    #text(115,-12,"(b)",cex=2)
    oz(lwd=2.5,add=TRUE)

#06 UTC
    screen(ind[3])
    par(gpar)
    zmax.t2 <- max(c(spd_gpats.t2,spd_pred.t2),na.rm=T)
    image(lon,rev(lat),t(apply(spd_gpats.t2,1,rev)),xlim=c(110,155),ylim=c(-45,-10),
               zlim=c(0,zmax.t2),
               xlab="",ylab="",
               xaxt="n",yaxt="n",
               col=colorRampPalette(rev(brewer.pal(9,"Spectral")))(64))
    text(115,-12,"06 UTC",cex=1.5,pos=4)
    oz(lwd=2.5,add=TRUE)
    
    screen(ind[4])
    par(gpar)
    image(lon,rev(lat),t(apply(spd_pred.t2,1,rev)),xlim=c(110,155),ylim=c(-45,-10),
               #main="GPATS observations of the average hourly flash count",
               zlim=c(0,zmax.t2),
               xlab="",ylab="",
               xaxt="n",yaxt="n",
               col=colorRampPalette(rev(brewer.pal(9,"Spectral")))(64))
    #text(115,-12,"(b)",cex=2)
    oz(lwd=2.5,add=TRUE)

#12 UTC
    screen(ind[5])
    par(gpar)
    zmax.t3 <- max(c(spd_gpats.t3,spd_pred.t3),na.rm=T)
    image(lon,rev(lat),t(apply(spd_gpats.t3,1,rev)),xlim=c(110,155),ylim=c(-45,-10),
               zlim=c(0,zmax.t3),
               xlab="",ylab="",
               xaxt="n",yaxt="n",
               col=colorRampPalette(rev(brewer.pal(9,"Spectral")))(64))
    text(115,-12,"12 UTC",cex=1.5,pos=4)
    oz(lwd=2.5,add=TRUE)
    
    screen(ind[6])
    par(gpar)
    image(lon,rev(lat),t(apply(spd_pred.t3,1,rev)),xlim=c(110,155),ylim=c(-45,-10),
               #main="GPATS observations of the average hourly flash count",
               zlim=c(0,zmax.t3),
               xlab="",ylab="",
               xaxt="n",yaxt="n",
               col=colorRampPalette(rev(brewer.pal(9,"Spectral")))(64))
    #text(115,-12,"(b)",cex=2)
    oz(lwd=2.5,add=TRUE)

#18 UTC
    
    screen(ind[7])
    par(gpar)
    zmax.t4 <- max(c(spd_gpats.t4,spd_pred.t4),na.rm=T)
    image(lon,rev(lat),t(apply(spd_gpats.t4,1,rev)),xlim=c(110,155),ylim=c(-45,-10),
               zlim=c(0,zmax.t4),
               xlab="",ylab="",
               xaxt="n",yaxt="n",
               col=colorRampPalette(rev(brewer.pal(9,"Spectral")))(64))
    text(115,-12,"18 UTC",cex=1.5,pos=4)
    oz(lwd=2.5,add=TRUE)
    
    screen(ind[8])
    par(gpar)
    image(lon,rev(lat),t(apply(spd_pred.t4,1,rev)),xlim=c(110,155),ylim=c(-45,-10),
               #main="GPATS observations of the average hourly flash count",
               zlim=c(0,zmax.t4),
               xlab="",ylab="",
               xaxt="n",yaxt="n",
               col=colorRampPalette(rev(brewer.pal(9,"Spectral")))(64))
    #text(115,-12,"(b)",cex=2)
    oz(lwd=2.5,add=TRUE)
    

    split.screen(c(4,1), screen=2)-> ind
    screen(ind[1])
    par(gpar)
    image.plot(zlim=c(0,zmax.t1),
               #legend.only=TRUE,smallplot=c(.1,.2, .2,.8),
               legend.only=TRUE,smallplot=c(.1,.3, .2,.8),
               col=colorRampPalette(rev(brewer.pal(9,"Spectral")))(64))
               #legend.lab="Days per year",legend.cex=1.25)

    screen(ind[2])
    par(gpar)
    image.plot(zlim=c(0,zmax.t2),
               #legend.only=TRUE,smallplot=c(.1,.2, .2,.8),
               legend.only=TRUE,smallplot=c(.1,.3, .2,.8),
               col=colorRampPalette(rev(brewer.pal(9,"Spectral")))(64))
               #legend.lab="Days per year",legend.cex=1.25)

    screen(ind[3])
    par(gpar)
    image.plot(zlim=c(0,zmax.t3),
               #legend.only=TRUE,smallplot=c(.1,.2, .2,.8),
               legend.only=TRUE,smallplot=c(.1,.3, .2,.8),
               col=colorRampPalette(rev(brewer.pal(9,"Spectral")))(64))
               #legend.lab="Days per year",legend.cex=1.25)

    screen(ind[4])
    par(gpar)
    image.plot(zlim=c(0,zmax.t4),
               #legend.only=TRUE,smallplot=c(.1,.2, .2,.8),
               legend.only=TRUE,smallplot=c(.1,.3, .2,.8),
               col=colorRampPalette(rev(brewer.pal(9,"Spectral")))(64))
               #legend.lab="Days per year",legend.cex=1.25)

    #close.screen(all.screens=T)
dev.off()

#Plot the difference between the model and the observations

diff.t1 <- t(apply(spd_pred.t1,1,rev))-
               t(apply(spd_gpats.t1,1,rev))
diff.t2 <- t(apply(spd_pred.t2,1,rev))-
               t(apply(spd_gpats.t2,1,rev))
diff.t3 <- t(apply(spd_pred.t3,1,rev))-
               t(apply(spd_gpats.t3,1,rev))
diff.t4 <- t(apply(spd_pred.t4,1,rev))-
               t(apply(spd_gpats.t4,1,rev))

col_levs <- 50

pdf(paste0("Lightning_diurnal_strike_difference",ncls,"_clusters.pdf"),
    width=10,height=10,paper="special")
    par(mar=c(0,0,0,0)+0.1,bty="n",cex=1.25,oma=c(0,0,0,0.5),mfrow=c(2,2))
#00 UTC
    max_abs.t1 <- max(abs(c(min(diff.t1,na.rm=T),
                                max(diff.t1,na.rm=T))))
    col_seq <- seq(-1.*max_abs.t1,max_abs.t1,length.out=col_levs+1)
    image.plot(lon,rev(lat),diff.t1,xlim=c(110,155),ylim=c(-45,-10),
               #zlim=c(0,zmax.summer),
               #xlab="Lon",ylab="Lat",
               xlab="",ylab="",
               xaxt="n",yaxt="n",
               col=colorRampPalette(rev(brewer.pal(9,"RdBu")))(col_levs),
               breaks=col_seq,legend.mar=5.5,legend.width=1.5)
    text(115,-12,"00 UTC",cex=1.5)
    oz(lwd=2.5,add=TRUE)
    
#06 UTC
    max_abs.t2 <- max(abs(c(min(diff.t2,na.rm=T),
                                max(diff.t2,na.rm=T))))
    col_seq <- seq(-1.*max_abs.t2,max_abs.t2,length.out=col_levs+1)
    image.plot(lon,rev(lat),diff.t2,xlim=c(110,155),ylim=c(-45,-10),
               #zlim=c(0,zmax.summer),
               #xlab="Lon",ylab="Lat",
               xlab="",ylab="",
               xaxt="n",yaxt="n",
               col=colorRampPalette(rev(brewer.pal(9,"RdBu")))(col_levs),
               breaks=col_seq,legend.mar=5.5,legend.width=1.5)
    text(115,-12,"06 UTC",cex=1.5)
    oz(lwd=2.5,add=TRUE)

#12 UTC
    max_abs.t3 <- max(abs(c(min(diff.t3,na.rm=T),
                                max(diff.t3,na.rm=T))))
    col_seq <- seq(-1.*max_abs.t3,max_abs.t3,length.out=col_levs+1)
    image.plot(lon,rev(lat),diff.t3,xlim=c(110,155),ylim=c(-45,-10),
               #zlim=c(0,zmax.summer),
               #xlab="Lon",ylab="Lat",
               xlab="",ylab="",
               xaxt="n",yaxt="n",
               col=colorRampPalette(rev(brewer.pal(9,"RdBu")))(col_levs),
               breaks=col_seq,legend.mar=5.5,legend.width=1.5)
    text(115,-12,"12 UTC",cex=1.5)
    oz(lwd=2.5,add=TRUE)

#18 UTC
    max_abs.t4 <- max(abs(c(min(diff.t4,na.rm=T),
                                max(diff.t4,na.rm=T))))
    col_seq <- seq(-1.*max_abs.t4,max_abs.t4,length.out=col_levs+1)
    image.plot(lon,rev(lat),diff.t4,xlim=c(110,155),ylim=c(-45,-10),
               #zlim=c(0,zmax.summer),
               #xlab="Lon",ylab="Lat",
               xlab="",ylab="",
               xaxt="n",yaxt="n",
               col=colorRampPalette(rev(brewer.pal(9,"RdBu")))(col_levs),
               breaks=col_seq,legend.mar=5.5,legend.width=1.5)
    text(115,-12,"18 UTC",cex=1.5)
    oz(lwd=2.5,add=TRUE)

dev.off()
 
 

#Plot box and whisker plots of the seasonal and diurnal cycles
strikes <- matrix(cbind(as.vector(spd_gpats.summer),
                 as.vector(spd_gpats.autumn),
                 as.vector(spd_gpats.winter),
                 as.vector(spd_gpats.spring),
                 as.vector(spd_pred.summer),
                 as.vector(spd_pred.autumn),
                 as.vector(spd_pred.winter),
                 as.vector(spd_pred.spring)))
obs <- matrix(c(rep("obs",length(spd_gpats.summer)*4),
         rep("model",length(spd_gpats.summer)*4)))
season <- matrix(rep(c(rep("DJF",length(spd_gpats.summer)),
                rep("MAM",length(spd_gpats.autumn)),
                rep("JJA",length(spd_gpats.winter)),
                rep("SON",length(spd_gpats.spring))),2))
seas_df <- data.frame("strikes"=strikes,
                      "obs"=factor(obs,levels=c("obs","model")),
                      "season"=factor(season,levels=c("DJF","MAM","JJA","SON")))
seas_df.m <- melt(seas_df,id.vars=c("obs","season"))


#pdf("boxplot_strikes_season.pdf",width=10,height=10,paper="special")
#p <- ggplot(seas_df.m,aes(x=season,y=value,color=obs))+
#         #geom_boxplot(notch=T)+
#         geom_boxplot()+
#         scale_y_continuous(limits=c(0,40))+
#         scale_color_grey()+
#         labs(x="Season",y="Strikes per year")+
#         theme(text=element_text(size=16))
#print(p)
#dev.off()

strikes.diurnal <- matrix(cbind(as.vector(spd_gpats.t1),
                                as.vector(spd_gpats.t2),
                                as.vector(spd_gpats.t3),
                                as.vector(spd_gpats.t4),
                                as.vector(spd_pred.t1),
                                as.vector(spd_pred.t2),
                                as.vector(spd_pred.t3),
                                as.vector(spd_pred.t4)))

obs <- matrix(c(rep("obs",length(spd_gpats.t1)*4),
         rep("model",length(spd_gpats.t1)*4)))

era_times <- matrix(rep(c(rep("00",length(spd_gpats.t1)),
                    rep("06",length(spd_gpats.t2)),
                    rep("12",length(spd_gpats.t3)),
                    rep("18",length(spd_gpats.t4))),2))

diurnal_df <- data.frame("strikes"=strikes.diurnal,
                         "obs"=factor(obs,levels=c("obs","model")),
                         "Time"=factor(era_times,levels=c("00","06","12","18")))
diurnal_df.m <- melt(diurnal_df,id.vars=c("obs","Time"))

#pdf("boxplot_strikes_diurnal.pdf",width=10,height=10,paper="special")
#p <- ggplot(diurnal_df.m,aes(x=Time,y=value,color=obs))+
#         #geom_boxplot(notch=T)+
#         geom_boxplot()+
#         scale_y_continuous(limits=c(0,40))+
#         scale_color_grey()+
#         labs(x="Time",y="Strikes per year")+
#         theme(text=element_text(size=16))
#print(p)
#dev.off()

pdf("boxplot_strikes_comparison.pdf",width=10,height=10,paper="special")
#par(mfrow=c(2,1))
require(gridExtra)
p1 <- ggplot(seas_df.m,aes(x=season,y=value,color=obs))+
         #geom_boxplot(notch=T)+
         geom_boxplot()+
         scale_y_continuous(limits=c(0,40))+
         scale_color_grey()+
         #labs(x="Season",y="Strikes per year")+
         labs(x="Season",y="Strikes per day")+
         theme(text=element_text(size=16))
p2 <- ggplot(diurnal_df.m,aes(x=Time,y=value,color=obs))+
         #geom_boxplot(notch=T)+
         geom_boxplot()+
         scale_y_continuous(limits=c(0,40))+
         scale_color_grey()+
         #labs(x="Time",y="Strikes per year")+
         labs(x="Time",y="Strikes per day")+
         theme(text=element_text(size=16))
grid.arrange(p1,p2,nrow=2)
dev.off()











#
#
#
#
##The annual climatology of the GPATS observations
##We had to divide by four when we were only looking at the maximum
##for around each ERA time.
##When max_lightning is derived from the SUM, we are regressing
##the sum of lightning that will happen in the next 6 hour interval, given
##the current atmospheric state.
##mxlpd <- apply(max_lightning,c(1,2),sum,na.rm=TRUE)/ndays_cls/4
##mxlpd <- apply(max_lightning,c(1,2),sum,na.rm=TRUE)/ndays_cls/24. #Per hour
#mxlpd <- apply(max_lightning,c(1,2),sum,na.rm=TRUE)/ndays_cls/4. #Per hour
##Apply the mask
#mxlpd[lsmask ==0] <- NA
#
##The annual climatology of the prediction
##mxlpdp <- full_clim/ndays_cls/4.
#mxlpdp <- full_clim/ndays_cls/4. #Per hour
#mxlpdp[lsmask==0] <- NA
#
#
#
##pdf(paste0("GPATS_vs_cluster_lightning_per_day_",ncls,"_clusters.pdf"),
#pdf(paste0("GPATS_vs_cluster_lightning_per_hour_",ncls,"_clusters.pdf"),
#    width=10,height=10,paper="special")
#image.plot(lon,rev(lat),t(apply(mxlpd,1,rev)),xlim=c(110,160),ylim=c(-45,-10),
#           #main="GPATS observations of the average daily maximum flash count",
#           #main="GPATS observations of the average daily flash count",
#           main="GPATS observations of the average hourly flash count",
#           zlim=c(0,max(c(mxlpdp,mxlpd),na.rm=T)),
#           xlab="Lon",ylab="Lat",
#           #col=rev(brewer.pal(10,"Spectral")))
#           #col=colorRampPalette(rev(brewer.pal(9,"Spectral")))(64),
#           #col=tim.colors(64),
#           #legend.lab="Lightning counts per day",horizontal=TRUE)
#           legend.lab="Lightning counts per ERA grid box per hour",horizontal=TRUE)
##map('worldHires',interior=T,lwd=2.5,add=T)
#oz(lwd=2.5,add=TRUE)
#
#image.plot(lon,rev(lat),t(apply(mxlpdp,1,rev)),xlim=c(110,160),ylim=c(-45,-10),
#           main="Cluster prediction of the average daily flash count",
#           zlim=c(0,max(c(mxlpdp,mxlpd),na.rm=T)),
#           xlab="Lon",ylab="Lat",
#           #col=rev(brewer.pal(10,"Spectral")))
#           #col=colorRampPalette(rev(brewer.pal(9,"Spectral")))(64),
#           #col=tim.colors(64),
#           #legend.lab="Lightning counts per day",horizontal=TRUE)
#           legend.lab="Lightning counts per ERA grid box per day",horizontal=TRUE)
#oz(lwd=2.5,add=TRUE)
#
#dev.off()
#
##Now a plot of the differences to get biases
#
#bias_lght <- mxlpd-mxlpdp
#col_levs <- 50
#max_abs_val <- max(abs(c(min(bias_lght,na.rm=T),max(bias_lght,na.rm=T))))
#col_seq <- seq(-1.*max_abs_val,max_abs_val,length.out=col_levs+1)
#
#
#pdf(paste0("Difference_clusterglm_gpats_lightning_counts_",ncls,"_clusters.pdf"),
#    width=10,height=10,paper="special")
#
#image.plot(lon,rev(lat),t(apply(bias_lght,1,rev)),
#           xlim=c(110,160),ylim=c(-45,-10),
#           main="Observations - model",
#           xlab="Lon",ylab="Lat",
#           col=colorRampPalette(rev(brewer.pal(9,"RdBu")))(col_levs),
#           breaks=col_seq,
#           #legend.lab="Lightning counts per day",horizontal=TRUE)
#           legend.lab="Lightning counts per hour",horizontal=TRUE)
#oz(lwd=2.5,add=T)
#
#dev.off()
#
##Difference as a fraction of the observations
#prop <- bias_lght/mxlpd
#prop[prop > 1] <- NA
#prop[prop < -1] <- NA
#max_abs_val <- max(abs(c(min(prop,na.rm=T),max(prop,na.rm=T))))
#col_seq <- seq(-1.*max_abs_val,max_abs_val,length.out=col_levs+1)
#image.plot(lon,rev(lat),t(apply(prop,1,rev)),xlim=c(110,160),ylim=c(-45,-10),
##image.plot(lon,rev(lat),diff/plgbm,xlim=c(110,160),ylim=c(-45,-10),
#        main="Proportionate difference between cluster and GPATS",
#        xlab="Lon",ylab="Lat",
#        #zlim=c(-1.,1),
#        #col=colorRampPalette(rev(brewer.pal(9,"Spectral")))(ncols))
#        #col=colorRampPalette(rev(brewer.pal(9,"RdBu")))(ncols))
#        col=colorRampPalette(rev(brewer.pal(9,"RdBu")))(col_levs),
#        breaks=col_seq)
#oz(lwd=2.5,add=TRUE)
#dev.off()
#
#
#
##Construct cumulative histograms
#
##max_lightning_t <- max_lightning
##max_lightning_t[,,] <- t(apply(max_lightning[,,],1,rev))
##max_lightning_t[!is.finite(max_lightning_t)] <- NA
#max_lightning_hist <- max_lightning
#
##colsh <- tim.colors(ncls) #Colours for the histograms
#colsh <- oce.colorsViridis(ncls) #Colours for the histograms
#tot_lght <- vector("list",length=ncls+1) #One extra for all clusters
#tot_lght[[1]] <- total_lightning[total_lightning>0]
#for (i in 1:ncls){
#    #tot_lght[[i+1]] <- max_lightning_t[which(cluster==i)]
#    #tot_lght[[i+1]] <- tot_lght[[i+1]][tot_lght[[i+1]]>0]
#    tot_lght[[i+1]] <- max_lightning_hist[which(cls.pred==i)]
#    tot_lght[[i+1]] <- tot_lght[[i+1]][tot_lght[[i+1]]>0]
#}
#cum_hists <- vector("list",length=ncls+1) #One extra for all clusters
#
#
##Plot the cumulative histogram of the lightning
#pdf("cum_hist_lightning.pdf",width=10,height=10)
##par(cex=1.5,mar=c(5,4,2,2))
##Plot the total lightning
#cum_hists[[1]] <- Ecdf(tot_lght[[1]],
#                       what="1-F",
#                       xlab="Lightning strikes/hour",
#                       ylab="Fraction",
#                       log="xy",
#                       ylim=c(1.e-06,1),
#                       #yaxt="n",
#                       #yaxt="y",
#                       lwd=5,
#                       subtitles=F)
#for (i in 1:ncls){
#cum_hists[[i+1]] <- Ecdf(tot_lght[[i+1]],
#                       what="1-F",
#                       lwd=5,
#                       subtitles=F,
#                       col=colsh[i],add=TRUE)
#}
#legend("topright",legend=as.character(1:ncls),fill=colsh)
#dev.off()
#
#
#
#
#
#
##prob.hail <- array(dim=c(nlon,nlat,ntime_cls))
#prob.lightning <- array(dim=c(nlon,nlat,ntime_cls))
#lightning.per.day <- array(dim=c(nlon,nlat,ntime_cls))
##hail.days.per.year<- array(dim=c(nlon,nlat,ntime_cls))
#
##Assign a probability based on a cluster
##prob.hail[,,] <- hlprob$prop.hail[cls.pred[,,]]
##prob.lightning[,,] <- hlprob$prop.lightning[cls.pred[,,]]
##prob.hail[,,] <- hlprob$prop.hail[cluster[,,]]
##Previously used all points
##prob.lightning[,,] <- hlprob$prop.lightning[cluster[,,]]
##Just use continental samples
#prob.lightning[,,] <- hlprob$land.prop.lightning[cluster[,,]]
#
#
###lightning.per.day[,,] <- hlprob$prop.lightning.strikes[cluster[,,]]
##lightning.per.day[,,] <- hlprob$lightning.per.cluster[cluster[,,]]
###lightning.per.day[,,] <- hlprob$avg.lightning.per.cluster[cluster[,,]]
###hail.days.per.year[,,] <- hlprob$days.hail[cluster[,,]]
#
##Now use regression relationship to evaluate the lightning per day
#
#
##Evaluate the number of hail days per year/cluster=x
##Evaluate the number of time a cluster occurs=n
##the number of hail days per year in each cluster is x*n
##Then sum up each cluster
#
##Should probably show an example of a specific storm.
##image.plot(lon,rev(lat),cluster[,,1],nlevel=ncls)
#pdf(paste0("Gap_Storm_",ncls,"_clusters.pdf"),width=10,height=10,paper="special")
#image.plot(lon,rev(lat),cluster[,,ind_gap_storm_date],nlevel=ncls,
#           main="Cluster map")
#map('worldHires',interior=T,lwd=2.5,add=T)
#oz(lwd=2.5,add=TRUE)
#
#image.plot(lon,rev(lat),prob.hail[,,ind_gap_storm_date],nlevel=ncls,
#           main="Hail probability map")
#map('worldHires',interior=T,lwd=2.5,add=T)
#oz(lwd=2.5,add=TRUE,coast=FALSE)
#
#image.plot(lon,rev(lat),prob.lightning[,,ind_gap_storm_date],nlevel=ncls,
#           main="Lightning probability map")
##map('worldHires',interior=T,lwd=2.5,add=T)
#oz(lwd=2.5,add=TRUE)
#dev.off()
#
##Now apply the sum
#
##hp <- apply(prob.hail,c(1,2),sum)/ntime_cls
##hp <- (apply(prob.hail,c(1,2),sum)/1460)*4*365.
#lp <- apply(prob.lightning,c(1,2),sum)/ntime_cls
##lpd <- apply(lightning.per.day,c(1,2),sum)*24/1460
##lpd <- apply(lightning.per.day,c(1,2),sum)*(24/4)/1460
##lpd <- (apply(lightning.per.day,c(1,2),sum)/1460)*4
##lpd <- (apply(lightning.per.day,c(1,2),sum)/1460)*4
##lpd <- (apply(lightning.per.day,c(1,2),sum)*24.)/(1460*4)
##lpd <- (apply(lightning.per.day,c(1,2),sum)*24.)/(1460*4)
##lpd0 <- roll_sum(lightning.per.day,4,by=4,align="left")
##lpd <- lpd0[lpd0 > 0]
#
##Evaluate the rolling sum of each fourth element
##Which corresponds to the sum of 00/06/12/18 predictions
##lpd0 <- apply(lightning.per.day,c(1,2),roll_sum,4,by=4,align="left",fill=NA)
#lpd0 <- apply(lightning.per.day,c(1,2),roll_mean,4,by=4,align="left",fill=NA)
##lpd0 <- apply(lightning.per.day,c(1,2),roll_max,4,by=4,align="left",fill=NA)
##We need to extract every fourth element as the roll_sum routine leaves it in
##there
#lpd0 <- lpd0[seq(1,dim(lpd0)[1],4),,]
#
##Now sum and average them over the year
##lpd <- apply(lpd0,c(2,3),sum,na.rm=T)/365.
##lpd <- apply(lpd0,c(2,3),sum,na.rm=T)/ntime_cls
#lpd <- apply(lpd0,c(2,3),sum,na.rm=T)/ndays_cls
##lpd <- apply(lpd0,c(2,3),sum,na.rm=T)*4.0/ndays_cls
#
#for (i in 1:ntime_cls){
#    lpd[lsmask==0] <- NA
#    lp[lsmask==0]<-NA
#}
#
##Now the hail days per year
##hdpy <- apply(hail.days.per.year,c(1,2),sum)/ntime_cls
##lp <- lp/sum(lp)
##image.plot(lon,rev(lat),hp,nlevel=10,xlim=c(110,160),ylim=c(-45,-10),
##           col=rev(brewer.pal(10,"Spectral")))
##map('worldHires',interior=T,lwd=2.5,add=T)
##image.plot(lon,rev(lat),lp,nlevel=10,xlim=c(110,160),ylim=c(-45,-10))
##map('worldHires',interior=T,lwd=2.5,add=T)
#
#pdf(paste0("Probability_of_lightning",ncls,"_clusters.pdf"),
#    width=10,height=10,paper="special")
#image.plot(lon,rev(lat),(lp*ndays_cls)/nyrs_cls,xlim=c(110,160),ylim=c(-45,-10),
#           xlab="Lon",ylab="Lat",
#           col=colorRampPalette(rev(brewer.pal(9,"Spectral")))(64),
#           legend.lab="Lightning days/year",horizontal=TRUE)
#oz(lwd=2,add=T)
#
##image.plot(lon,rev(lat),(hp*ndays_cls)/nyrs_cls,xlim=c(110,160),ylim=c(-45,-10),
##           xlab="Lon",ylab="Lat",
##           #zlim=c(0,3),
##           col=colorRampPalette(rev(brewer.pal(9,"Spectral")))(64),
##           legend.lab="Hail days/year",horizontal=TRUE)
##oz(lwd=2,add=T)
#dev.off()
#
#
#pdf(paste0("Predicted_lightning_strikes_per_day_",ncls,"_clusters.pdf"),
#    width=10,height=10,paper="special")
#image.plot(lon,rev(lat),lpd,xlim=c(110,160),ylim=c(-45,-10),
##image.plot(lon,rev(lat),nlevel=10,lpd/4,xlim=c(110,160),ylim=c(-45,-10),
#           #zlim=c(0,22),
#           #zlim=c(0,7),
#           xlab="Lon",ylab="Lat",
#           #col=rev(brewer.pal(10,"Spectral")))
#           #col=colorRampPalette(rev(brewer.pal(9,"Spectral")))(ncols),
#           col=colorRampPalette(rev(brewer.pal(9,"Spectral")))(64),
#           legend.lab="Lightning strikes per day",horizontal=TRUE)
##map('worldHires',interior=T,lwd=2.5,add=T)
#oz(lwd=2.5,add=TRUE)
#dev.off()
#
##pdf(paste0("Predicted_hail_days_per_year_",ncls,"_clusters.pdf"),
##    width=10,height=10,paper="special")
###image.plot(lon,rev(lat),hdpy,nlevel=ncls,xlim=c(110,160),ylim=c(-45,-10),
##image.plot(lon,rev(lat),(hp*ndays_cls)/nyrs_cls,nlevel=ncls,xlim=c(110,160),ylim=c(-45,-10),
##           xlab="Lon",ylab="Lat",
##           #col=rev(brewer.pal(10,"Spectral")))
##           col=colorRampPalette(rev(brewer.pal(9,"Spectral")))(64),
##           #col=tim.colors(64),
##           legend.lab="Hail days per year",horizontal=TRUE)
###map('worldHires',interior=T,lwd=2.5,add=T)
##oz(lwd=2.5,add=TRUE)
##dev.off()
#
##Now make a climatology of the GPATS data.
##To compare apples with apples we need to evaluate the maximum
##lightning frequency on any given day
##Finding the max in any 24-hour period
##I think we should use this as to construct the probabilities
##we use the max lightning at any time during the day.
##max_lightning <- apply(total_lightning,c(1,2),roll_max,24,by=24,align="left",fill=NA)
#
##Find the max in window around each of the ERA time
##max_lightning <- apply(total_lightning,c(1,2),roll_max,3,by=6,
##                       align="center",fill=NA,na.rm=TRUE)
##Finding the mean around each of the ERA-time
#max_lightning <- apply(total_lightning,c(1,2),roll_mean,3,by=6,
#                       align="center",fill=NA,na.rm=TRUE)
##max_lightning <- apply(total_lightning,c(1,2),roll_sum,3,by=6,
##                       align="center",fill=NA,na.rm=TRUE)
##max_lightning <- apply(total_lightning,c(1,2),roll_mean,1,by=6,
##                       align="center",fill=NA,na.rm=TRUE)
##Extracting the lightning at each ERA time.
##max_lightning <- total_lightning[,,seq(1,dim(total_lightning)[3],6)]
#
##The mean lightning counts over the whole period.
##mean_lightning <- apply(total_lightning,c(1,2),mean)
#
##Extract only every appropriate elements
##Extractine when we have found the max. lightning at any point during the day
##max_lightning <- max_lightning[seq(1,dim(max_lightning)[1],24),,]
#
##Extracting when we have found the max around each ERA time
#max_lightning <- max_lightning[seq(2,dim(max_lightning)[1],6),,]
##max_lightning <- max_lightning[,,seq(2,dim(max_lightning)[1],4)]
#
#
##Now roll sum by 4 (00/06/12/18) to get daily total.
##We do this if we have extracted the lightning at the ERA times
##We use this if we have evaluated at the ERA time above
##max_lightning <-  apply(max_lightning,c(2,3),roll_sum,4,by=4,
##                        align="left",fill=NA,na.rm=TRUE)
##max_lightning <-  apply(max_lightning,c(2,3),roll_max,4,by=4,
##                        align="left",fill=NA,na.rm=TRUE)
#     #!!!!!Why have you done the below?
#     #This is the mean!
#     #You need to check this and fix it!!!! 2 May 2017
#max_lightning <-  apply(max_lightning,c(2,3),roll_mean,4,by=4,
#                        align="left",fill=NA,na.rm=TRUE)
##max_lightning <-  apply(max_lightning,c(1,2),roll_sum,4,by=4,
##                        align="left",fill=NA,na.rm=TRUE)
##max_lightning <-  apply(max_lightning,c(1,2),roll_mean,4,by=4,
##                        align="left",fill=NA,na.rm=TRUE)
#max_lightning <- max_lightning[seq(1,dim(max_lightning)[1],4),,]
##max_lightning <- max_lightning[seq(1,dim(max_lightning)[1],24),,]
##Sum and average
##mxlpdg <- apply(max_lightning,c(2,3),sum,na.rm=T)/(dim(max_lightning)[1])
#mxlpdg <- apply(max_lightning,c(2,3),sum,na.rm=T)/ndays_cls
##mxlpdg <- apply(max_lightning,c(2,3),sum,na.rm=T)/(dim(max_lightning)[3])
##mxlpdg <- apply(max_lightning,c(2,3),sum,na.rm=T)/((dim(max_lightning)[1])*4)
#
##Transpose the arrays for plotting
#mxlpdg <- t(apply(mxlpdg,1,rev))
##mnlpdg <- t(apply(mean_lightning,1,rev))
#
##Apply the land sea mask
#mxlpdgm <- mxlpdg
##mnlpdgm <- mnlpdg
#for (i in 1:length(time)){
#    mxlpdgm[lsmask==0] <- NA
##    mnlpdgm[lsmask==0] <- NA
#}
#
##Evaluate and plot the cumulative histograms of lightning
##for all GPATS observations and for each cluster
#max_lightning <- apply(total_lightning,c(1,2),roll_max,3,by=6,
#                       #align="center",fill=NA,na.rm=TRUE)
#                       align="center",na.rm=TRUE)
#max_lightning <- max_lightning[seq(1,dim(max_lightning)[1],6),,]
#max_lightning <- aperm(max_lightning,c(2,3,1))
##We need to transpose to make sure we have the
##same indices as for the clusters
#max_lightning_t <- max_lightning
#max_lightning_t[,,] <- t(apply(max_lightning[,,],1,rev))
#max_lightning_t[!is.finite(max_lightning_t)] <- NA
##max_lightning_t <- t(apply(max_lightning[,,],c(2,3),rev))
#
#
#
##Construct cumulative histograms
##colsh <- tim.colors(ncls) #Colours for the histograms
#colsh <- oce.colorsViridis(ncls) #Colours for the histograms
#tot_lght <- vector("list",length=ncls+1) #One extra for all clusters
#tot_lght[[1]] <- total_lightning[total_lightning>0]
#for (i in 1:ncls){
#    tot_lght[[i+1]] <- max_lightning_t[which(cluster==i)]
#    tot_lght[[i+1]] <- tot_lght[[i+1]][tot_lght[[i+1]]>0]
#}
#cum_hists <- vector("list",length=ncls+1) #One extra for all clusters
#
#
##Plot the cumulative histogram of the lightning
#pdf("cum_hist_lightning.pdf",width=10,height=10)
##par(cex=1.5,mar=c(5,4,2,2))
##Plot the total lightning
#cum_hists[[1]] <- Ecdf(tot_lght[[1]],
#                       what="1-F",
#                       xlab="Lightning strikes/hour",
#                       ylab="Fraction",
#                       log="xy",
#                       ylim=c(1.e-06,1),
#                       yaxt="n",
#                       lwd=5,
#                       subtitles=F)
#for (i in 1:ncls){
#cum_hists[[i+1]] <- Ecdf(tot_lght[[i+1]],
#                       what="1-F",
#                       lwd=5,
#                       subtitles=F,
#                       col=colsh[i],add=TRUE)
#}
#dev.off()
#
#
#
##Ecdf(total_lightning[total_lightning>0],
##     what="1-F",
##     xlab="Lightning strikes per hour",
##     ylab="Fraction",
##     log="xy",ylim=c(1.e-06,1),yaxt="n",
##     subtitles=F,lwd=3)
##
##mxlt1 <- max_lightning_t[which(cluster==1)]
##mxlt1 <- mxlt1[mxlt1 >0]
##mxlt5 <- max_lightning_t[which(cluster==5)]
##mxlt5 <- mxlt5[mxlt5 >0]
##mxlt10 <- max_lightning_t[which(cluster==10)]
##mxlt10 <- mxlt10[mxlt10 >0]
###Ecdf(max_lightning[which(cluster==1)],
##Ecdf(mxlt1,
##     what="1-F",
##     xlab="Lightning strikes per hour",
##     ylab="Fraction",
##     log="xy",ylim=c(1.e-06,1),yaxt="n",
##     subtitles=F,lwd=3,
##     col="red",add=TRUE)
##     #col="red")
##
##Ecdf(mxlt5,
##     what="1-F",
##     #xlab="Lightning strikes per hour",
##     #ylab="Fraction",
##     #log="xy",ylim=c(1.e-06,1),yaxt="n",
##     subtitles=F,lwd=3,
##     col="blue",add=TRUE)
##
##
###Ecdf(max_lightning[cluster==10],
##Ecdf(mxlt10,
##     what="1-F",
##     #xlab="Lightning strikes per hour",
##     #ylab="Fraction",
##     #log="xy",ylim=c(1.e-06,1),yaxt="n",
##     subtitles=F,lwd=3,
##     col="green",add=TRUE)
##
#
#
##aty <- c(1e-08,1e-07,1e-06,1e-05,1e-04,1e-03,1e-02,1e-01,1e0)
##aty <- c(1e-06,1e-05,1e-04,1e-03,1e-02,1e-01,1e0)
###labelsY <- expression(paste0(rep(10,7),"^",seq(-6,0)))
##labelsY <- parse(text=paste0(rep(10,7),"^",seq(-6,0)))
##axis(2,at=aty,labels=labelsY,las=1)
##grid(col="grey20")
##box()
#
#eccdf <- function(x){return(1.-ecdf(x)(x))}
#lgt0 <- total_lightning[total_lightning>0]
##lecdf <- ecdf(lgt0)
#lecdf <- eccdf(lgt0)
#r <- range(lgt0)
#plot(eccdf)
#
#dev.off()
#
#
#
#
#
##total_lightning[,,] <- t(apply(total_lightning[,,],1,rev))
#
#pdf(paste0("GPATS_vs_cluster_lightning_per_day_",ncls,"_clusters.pdf"),
#    width=10,height=10,paper="special")
#image.plot(lon,rev(lat),mxlpdgm,xlim=c(110,160),ylim=c(-45,-10),
#           main="GPATS observations of the average daily maximum flash count",
#           zlim=c(0,max(c(lpd,mxlpdgm),na.rm=T)),
#           xlab="Lon",ylab="Lat",
#           #col=rev(brewer.pal(10,"Spectral")))
#           #col=colorRampPalette(rev(brewer.pal(9,"Spectral")))(64),
#           #col=tim.colors(64),
#           legend.lab="Lightning counts per day",horizontal=TRUE)
##map('worldHires',interior=T,lwd=2.5,add=T)
#oz(lwd=2.5,add=TRUE)
#
#image.plot(lon,rev(lat),lpd,xlim=c(110,160),ylim=c(-45,-10),
#           main="Cluster prediction of the average daily maximum flash count",
#           zlim=c(0,max(c(lpd,mxlpdgm),na.rm=T)),
#           xlab="Lon",ylab="Lat",
#           #col=rev(brewer.pal(10,"Spectral")))
#           #col=colorRampPalette(rev(brewer.pal(9,"Spectral")))(64),
#           #col=tim.colors(64),
#           legend.lab="Lightning counts per day",horizontal=TRUE)
#oz(lwd=2.5,add=TRUE)
#
#dev.off()
#
#
#difflt <- lpd-mxlpdgm
#col_levs <- 50
#max_abs_val <- max(abs(c(min(difflt,na.rm=T),max(difflt,na.rm=T))))
#col_seq <- seq(-1.*max_abs_val,max_abs_val,length.out=col_levs+1)
#
#pdf(paste0("Difference_cluster_gpats_lightning_counts_",ncls,"_clusters.pdf"),
#    width=10,height=10,paper="special")
#
#image.plot(lon,rev(lat),lpd-mxlpdgm,
#           xlim=c(110,160),ylim=c(-45,-10),
#           xlab="Lon",ylab="Lat",
#           col=colorRampPalette(rev(brewer.pal(9,"RdBu")))(col_levs),
#           breaks=col_seq,
#           legend.lab="Lightning counts per day",horizontal=TRUE)
#oz(lwd=2.5,add=T)
#
#dev.off()
#
##image.plot(lon_gpats,rev(lat),mnlpdg,nlevel=10,xlim=c(110,160),ylim=c(-45,-10),
##           xlab="Lon",ylab="Lat",
##           #col=rev(brewer.pal(10,"Spectral")))
##           col=colorRampPalette(rev(brewer.pal(9,"Spectral")))(64),
##           #col=tim.colors(64),
##           legend.lab="LIghtning counts/year",horizontal=TRUE)
###map('worldHires',interior=T,lwd=2.5,add=T)
##oz(lwd=2.5,add=TRUE)
#
##Evaluate the proportion of the days in each grid box that have lightning
#sum_lightning <- apply(total_lightning,c(1,2),roll_sum,24,by=24,align="left",fill=NA)
#sum_lightning <- sum_lightning[seq(1,dim(sum_lightning)[1],24),,]
##lightning_in_grid <- which(sum_lightning > 0)
##tlgb - total lightning in grid box
#tlgb <- apply(sum_lightning[,,],c(2,3),sum,na.rm=T)
#tlgb <- t(apply(tlgb,1,rev))
##Now evaluate the proportion that each has a grid box has lightning
#lightning_in_grid <- which(sum_lightning[,,] > 0)
#plgbs <- sum_lightning[,,]
#plgbs[lightning_in_grid] <- 1
#plgb <- apply(plgbs,c(2,3),sum,na.rm=T)/dim(plgbs)[1]
#plgb <- t(apply(plgb,1,rev))
##Apply the land sea mask
#plgbm <- plgb
#plgbm[lsmask==0] <- NA
#
##Best z lim for comparing cluster prediction and GPATS data
##zmax <- c(lp/dim(prob.lightning)[3],plgb)[which.max(c(lp/dim(prob.lightning)[3],plgb))]
##zmax <- c(lp,plgb)[which.max(c(lp,plgb))]
#zmax <- c(lp,plgbm)[which.max(c(lp,plgb))]
#
##pdf(paste0("Cluster_vs_GPATS_proportion_lightning_2008_",ncls,"_cluster.pdf"),
#pdf(paste0("Cluster_vs_GPATS_proportion_lightning_",ncls,"_cluster.pdf"),
#    width=10,height=10,paper="special")
##par(mfrow=c(2,1))
#par(mar=c(3,3,1,1))
##split.screen( rbind(c(0, .8,0,1), c(.8,1,0,1)))
##split.screen( rbind(c(0, .85,0,1), c(.85,1,0,1)))
#split.screen( rbind(c(0, .9,0,1), c(.9,1,0,1)))
#split.screen(c(2,1), screen=1)-> ind
#screen( ind[1])
##image.plot(lon_gpats,rev(lat_gpats),plgb,xlim=c(110,160),ylim=c(-45,-10),
#image(lon_gpats,rev(lat_gpats),plgbm,xlim=c(110,160),ylim=c(-45,-10),
#           main="GPATS observations",
#           xlab="Lon",ylab="Lat",
#           zlim=c(0,zmax),
#           #col=rev(brewer.pal(10,"Spectral")))
#           col=colorRampPalette(rev(brewer.pal(9,"Spectral")))(ncols))
#           #col=tim.colors(64),
#           #legend.lab="Proportion of days with lightning",horizontal=TRUE)
##map('worldHires',interior=T,lwd=2.5,add=T)
#oz(lwd=2.5,add=TRUE)
#
#screen(ind[2])
##image.plot(lon,rev(lat),lp/dim(prob.lightning)[3],xlim=c(110,160),ylim=c(-45,-10),
##image(lon,rev(lat),lp/dim(prob.lightning)[3],xlim=c(110,160),ylim=c(-45,-10),
#image(lon,rev(lat),lp,xlim=c(110,160),ylim=c(-45,-10),
#           main="Cluster prediction",
#           xlab="Lon",ylab="Lat",
#           zlim=c(0,zmax),
#           #col=rev(brewer.pal(10,"Spectral")))
#           col=colorRampPalette(rev(brewer.pal(9,"Spectral")))(ncols))
#           #col=tim.colors(64),
#           #legend.lab="Proportion of days with lightning",horizontal=TRUE)
##map('worldHires',interior=T,lwd=2.5,add=T)
#oz(lwd=2.5,add=TRUE)
#
#screen(2)
##image.plot(zlim=c(0,zmax),legend.only=TRUE,smallplot=c(.85,.9, .1,.9),
##image.plot(zlim=c(0,zmax),legend.only=TRUE,smallplot=c(.2,.5, .2,.8),
#image.plot(zlim=c(0,zmax),legend.only=TRUE,smallplot=c(.1,.2, .2,.8),
##image.plot(zlim=c(0,zmax),legend.only=TRUE,smallplot=c(.9,.95, .2,.8),
#           col=colorRampPalette(rev(brewer.pal(9,"Spectral")))(ncols))
#
#dev.off()
#
#pdf(paste0("Difference_cluster_GPATS_",ncls,"_cluster.pdf"),
#    width=10,height=10,paper="special")
#diff <- lp-plgbm
##See http://stackoverflow.com/questions/33750235/plotting-a-raster-with-the-color-ramp-diverging-around-zero
#col_levs <- 50
#max_abs_val <- max(abs(c(min(diff,na.rm=T),max(diff,na.rm=T))))
#col_seq <- seq(-1.*max_abs_val,max_abs_val,length.out=col_levs+1)
#image.plot(lon,rev(lat),diff,xlim=c(110,160),ylim=c(-45,-10),
##image.plot(lon,rev(lat),diff/plgbm,xlim=c(110,160),ylim=c(-45,-10),
#        main="Cluster prediction minus GPATS",
#        xlab="Lon",ylab="Lat",
#        #zlim=c(0,zmax),
#        #col=colorRampPalette(rev(brewer.pal(9,"Spectral")))(ncols))
#        #col=colorRampPalette(rev(brewer.pal(9,"RdBu")))(ncols))
#        col=colorRampPalette(rev(brewer.pal(9,"RdBu")))(col_levs),
#        breaks=col_seq)
#oz(lwd=2.5,add=TRUE)
#
##Difference as a fraction of the observations
#prop <- diff/plgbm
#prop[prop > 1] <- NA
#max_abs_val <- max(abs(c(min(prop,na.rm=T),max(prop,na.rm=T))))
#col_seq <- seq(-1.*max_abs_val,max_abs_val,length.out=col_levs+1)
#image.plot(lon,rev(lat),prop,xlim=c(110,160),ylim=c(-45,-10),
##image.plot(lon,rev(lat),diff/plgbm,xlim=c(110,160),ylim=c(-45,-10),
#        main="Proportionate difference between cluster and GPATS",
#        xlab="Lon",ylab="Lat",
#        #zlim=c(-1.,1),
#        #col=colorRampPalette(rev(brewer.pal(9,"Spectral")))(ncols))
#        #col=colorRampPalette(rev(brewer.pal(9,"RdBu")))(ncols))
#        col=colorRampPalette(rev(brewer.pal(9,"RdBu")))(col_levs),
#        breaks=col_seq)
#oz(lwd=2.5,add=TRUE)
#dev.off()
#
#
##Now do a seasonal analysis
##Get the cluster data in a daily average
##First we need to extract the GPATS time which correspond to the ERA times
#library(lubridate)
#rtime <- rtime_gpats[seq(1,ntime_gpats,6)]
#
#yr.2008 <- which(year(rtime) == 2008)
#djf <- c(12,1,2)
#mam <- c(3,4,5)
#jja <- c(6,7,8)
#son <- c(9,10,11)
#
#summer <- which(month(rtime) %in% djf)
#autumn <- which(month(rtime) %in% mam)
#winter <- which(month(rtime) %in% jja)
#spring <- which(month(rtime) %in% son)
#
##First the cluster prediction
#lpd0.summer <- apply(lightning.per.day[,,yr.2008[summer]],
#                     c(1,2),roll_max,4,by=4,align="left",fill=NA)
#lpd0.autumn <- apply(lightning.per.day[,,yr.2008[autumn]],
#                     c(1,2),roll_max,4,by=4,align="left",fill=NA)
#lpd0.winter <- apply(lightning.per.day[,,yr.2008[winter]],
#                     c(1,2),roll_max,4,by=4,align="left",fill=NA)
#lpd0.spring <- apply(lightning.per.day[,,yr.2008[spring]],
#                     c(1,2),roll_max,4,by=4,align="left",fill=NA)
##We need to extract every fourth element as the roll_sum routine leaves it in
##there
#lpd0.summer <- lpd0.summer[seq(1,dim(lpd0)[1],4),,]
#lpd0.autumn <- lpd0.autumn[seq(1,dim(lpd0)[1],4),,]
#lpd0.winter <- lpd0.winter[seq(1,dim(lpd0)[1],4),,]
#lpd0.spring <- lpd0.spring[seq(1,dim(lpd0)[1],4),,]
#
##Now sum and average them over the year
##lpd <- apply(lpd0,c(2,3),sum,na.rm=T)/365.
##lpd <- apply(lpd0,c(2,3),sum,na.rm=T)/ntime_cls
#lpd.summer <- apply(lpd0.summer,c(2,3),sum,na.rm=T)/ndays_cls
#lpd.autumn <- apply(lpd0.autumn,c(2,3),sum,na.rm=T)/ndays_cls
#lpd.winter <- apply(lpd0.winter,c(2,3),sum,na.rm=T)/ndays_cls
#lpd.spring <- apply(lpd0.spring,c(2,3),sum,na.rm=T)/ndays_cls
#
##for (i in 1:ntime_cls){
#    lpd.summer[lsmask==0] <- NA
#    lpd.autumn[lsmask==0] <- NA
#    lpd.winter[lsmask==0] <- NA
#    lpd.spring[lsmask==0] <- NA
##}
#
#pdf(paste0("Cluster_seasonal_cycle",ncls,"_cluster.pdf"),
#    width=10,height=10,paper="special")
#image.plot(lon,rev(lat),lpd.summer,xlim=c(110,160),ylim=c(-45,-10),
#        main="Cluster seasonal cycle",
#        xlab="Lon",ylab="Lat",
#        col=colorRampPalette(rev(brewer.pal(9,"Spectral")))(ncols))
#oz(lwd=2.5,add=TRUE)
#
#image.plot(lon,rev(lat),lpd.autumn,xlim=c(110,160),ylim=c(-45,-10),
#        main="Cluster seasonal cycle",
#        xlab="Lon",ylab="Lat",
#        col=colorRampPalette(rev(brewer.pal(9,"Spectral")))(ncols))
#oz(lwd=2.5,add=TRUE)
#
#image.plot(lon,rev(lat),lpd.winter,xlim=c(110,160),ylim=c(-45,-10),
#        main="Cluster seasonal cycle",
#        xlab="Lon",ylab="Lat",
#        col=colorRampPalette(rev(brewer.pal(9,"Spectral")))(ncols))
#oz(lwd=2.5,add=TRUE)
#
#image.plot(lon,rev(lat),lpd.spring,xlim=c(110,160),ylim=c(-45,-10),
#        main="Cluster seasonal cycle",
#        xlab="Lon",ylab="Lat",
#        col=colorRampPalette(rev(brewer.pal(9,"Spectral")))(ncols))
#oz(lwd=2.5,add=TRUE)
#
#
#dev.off()
#
#
#
##Now the GPATS observations
#
#
###lpd0 <- apply(lightning.per.day,c(1,2),roll_mean,4,by=4,align="left",fill=NA)
#
##This is the two-hourly averaged lightning at the ERA times.
##Since the GPATS data has been offset so that the lightning count is the 
##number measured during over the hour following the time signature,
##we will include the previous, so we have summed the hour
##preceding and following the ERA time
##mngpats <- apply(total_lightning,c(1,2),roll_mean,2,by=6,align="right",fill=NA)
##mxgpats  <- apply(total_lightning,c(1,2),roll_max,2,by=6,align="right",fill=NA)
##mxgpats <- mxgpats[seq(2,dim(mxgpats)[1],6),,]
#mxgpats  <- apply(total_lightning,c(1,2),roll_max,24,by=24,align="left",fill=NA)
#mxgpats <- mxgpats[seq(1,dim(mxgpats)[1],24),,]
#
##mxgpats0.summer <- apply(mxgpats[yr.2008[summer],,],
##                        c(2,3),roll_max,4,by=4,align="left",fill=NA)
##mxgpats0.autumn <- apply(mxgpats[yr.2008[autumn],,],
##                        c(2,3),roll_max,4,by=4,align="left",fill=NA)
##mxgpats0.winter <- apply(mxgpats[yr.2008[winter],,],
##                        c(2,3),roll_max,4,by=4,align="left",fill=NA)
##mxgpats0.spring <- apply(mxgpats[yr.2008[spring],,],
##                        c(2,3),roll_max,4,by=4,align="left",fill=NA)
#
#mxgpats0.summer <- apply(mxgpats[yr.2008[summer],,],
#                         c(2,3),sum,na.rm=TRUE)
#
#mxgpats0.autumn <- apply(mxgpats[yr.2008[autumn],,],
#                         c(2,3),sum,na.rm=TRUE)
#mxgpats0.winter <- apply(mxgpats[yr.2008[winter],,],
#                         c(2,3),sum,na.rm=TRUE)
#mxgpats0.spring <- apply(mxgpats[yr.2008[spring],,],
#                         c(2,3),sum,na.rm=TRUE)
#
##mxgpats0.summer <- mxgpats0.summer[seq(1,dim(mxgpats0.summer)[1],4),,]
##mxgpats0.autumn <- mxgpats0.autumn[seq(1,dim(mxgpats0.autumn)[1],4),,]
##mxgpats0.winter <- mxgpats0.winter[seq(1,dim(mxgpats0.winter)[1],4),,]
##mxgpats0.spring <- mxgpats0.spring[seq(1,dim(mxgpats0.spring)[1],4),,]
#
##Now sum and average
##mxgpats.summer <- apply(mxgpats0.summer,c(2,3),sum,na.rm=TRUE)/dim(mxgpats0.summer)
##mxgpats.autumn <- apply(mxgpats0.autumn,c(2,3),sum,na.rm=TRUE)/dim(mxgpats0.autumn)
##mxgpats.winter <- apply(mxgpats0.winter,c(2,3),sum,na.rm=TRUE)/dim(mxgpats0.winter)
##mxgpats.spring <- apply(mxgpats0.spring,c(2,3),sum,na.rm=TRUE)/dim(mxgpats0.spring)
#
#mxgpats.summer <- mxgpats0.summer/length(yr.2008[summer])
#mxgpats.autumn <- mxgpats0.autumn/length(yr.2008[autumn])
#mxgpats.winter <- mxgpats0.winter/length(yr.2008[winter])
#mxgpats.spring <- mxgpats0.spring/length(yr.2008[spring])
#
##Transpose and apply mask
#mxgpats.summer <- t(apply(mxgpats.summer,1,rev))
#mxgpats.autumn <- t(apply(mxgpats.autumn,1,rev))
#mxgpats.winter <- t(apply(mxgpats.winter,1,rev))
#mxgpats.spring <- t(apply(mxgpats.spring,1,rev))
##Apply the land sea mask
#
#mxgpats.summer[lsmask==0] <- NA
#mxgpats.autumn[lsmask==0] <- NA
#mxgpats.winter[lsmask==0] <- NA
#mxgpats.spring[lsmask==0] <- NA
#
##Now plot all of them together
##par(mar=c(3,3,1,1))
##split.screen( rbind(c(0, .8,0,1), c(.8,1,0,1)))
###split.screen( rbind(c(0, .85,0,1), c(.85,1,0,1)))
##split.screen( rbind(c(0, .9,0,1), c(.9,1,0,1)))
##split.screen(c(4,2), screen=1)-> ind
##screen( ind[1])
##
##image.plot(lon,rev(lat),mxgpats.summer)
#
#
#
##Now apply the mask
##for (i in 1:
#
#
#
##---------- These are fragments of code used in testing ----------#
##---------- They may come in useful later!              ----------#
#
#
##(1) Regression stuff
###Now perform regression
###library(VGAM)
##
##model3.p <- glm(LGHT~CAPE+LR75+H5_T+S06+PMR,
##                 data=reg_data_10,family=poisson)
##
##
##
##cls10 <- which(reg_data$CLS == 10)
##cls10_10 <- which(reg_data$CLS == 10 & reg_data$LGHT == 10)
##cls10_500 <- which(reg_data$CLS == 10 & reg_data$LGHT == 500)
##cls10_1000 <- which(reg_data$CLS == 10 & reg_data$LGHT == 1000)
##test <- which(reg_data$LGHT[cls10] == 500)
##reg_data_n_10 <- reg_data_n[cls10,]
##reg_data_10 <- reg_data[cls10,]
##
###reg_data2_n <- reg_data_n[cls10,]
###reg_data2 <- reg_data[cls10,]
##
##
##
###Now the regression
###model <- glm(LGHT[gloop]~CAPE[gloop]+LR75[gloop]+H5_T[gloop]+S06[gloop]+PMR[gloop]-1,data=reg_data, family=quasipoisson)
##
##
###model2 <- zeroinfl(LGHT[cls10]~CAPE[cls10]+LR75[cls10]+H5_T[cls10]+S06[cls10]+PMR[cls10]-1,
##                 data=reg_data,dist="negbin")
##
###Try a negative binomial distribution
###model3 <- glm.nb(LGHT[cls10]~CAPE[cls10]+LR75[cls10]+H5_T[cls10]+S06[cls10]+PMR[cls10],
###                 data=reg_data,link=log)
##model3_n <- glm.nb(LGHT~CAPE+LR75+H5_T+S06+PMR,
##                 data=reg_data_n_10,link=log)
##
##model3 <- glm.nb(LGHT~CAPE+LR75+H5_T+S06+PMR,
##                 data=reg_data_10,link=log)
##
###Try Poisson distributions
##model3_n.p <- glm(LGHT~CAPE+LR75+H5_T+S06+PMR,
##                 data=reg_data_n_10,family=poisson)
##
##model3.p <- glm(LGHT~CAPE+LR75+H5_T+S06+PMR,
##                 data=reg_data_10,family=poisson)
##
##
##newdata2_n <- reg_data_n[cls10_10,]
##newdata2 <- reg_data[cls10_10,]
###newdata2_n <- reg_data_n[cls10_500[1:10],]
###newdata2 <- reg_data[cls10_500[1:10],]
###newdata2_n <- reg_data_n[cls10_1000,]
###newdata2 <- reg_data[cls10_1000,]
##
##newdata2_n$plght.nb <- predict(model3_n,newdata=newdata2_n,type="response")
##newdata2$plght.nb <- predict(model3,newdata=newdata2,type="response")
##newdata2_n$plght.p <- predict(model3_n,newdata=newdata2_n,type="response")
##newdata2$plght.p <- predict(model3,newdata=newdata2,type="response")
##
##test_pred <- as.numeric(coefficients(model3))*
##             as.numeric(reg_data[cls10_500[1:10],2:6])
##
##
#
##################### ++++++++++ ####################
#
##(2) Performing the lightning prediction with loop. Too slow!
##strikes <- array(NA,dim=c(nlon,nlat,ntime_cls))
###Do a test using loops
##for (i in 1:nlon){
##    for (j in 1:nlat){
##    #for (k in 1:ntime_cls){
##    for (k in 1:100){
##        plg <- runif(1) #A generated random number
##        ctmp <- cluster[i,j,k]
##        if (is.na(ctmp)){
##            strikes[i,j,k] == NA
##        }
##        else{
##            if (plg >= land.prop.lightning[ctmp]){
##                strikes[i,j,k] = 0
##            }
##            else{
##                tdf <- data.frame(CAPE=mucape_n[i,j,k],
##                                  LR75=lr75_n[i,j,k],
##                                  H5_T=h5_temp_n[i,j,k],
##                                  S06=s06_n[i,j,k],
##                                  PMR=pmr_n[i,j,k])
##                strikes[i,j,k] = predict(fit.p[[ctmp]],tdf)
##            }
##        }
##    }
##    }
##}
##
#
#
##Try a zero inflated Poisson regression
##library(pscl)
#
##-----NOTES!!!!-----
##test.zip <-zeroinfl(LGHT~CAPE+LR75+H5_T+S06+PMR,data=train_cls[[10]],dist="poisson",EM=TRUE)
##The mean of the predicted lightning counts is close to that
##of the observation
##i.e. mean(train_cls[[10]]$LGHT) ~ mean(predict(test.zip,test_cls[[10]]))
#
##Now try a negative binomial. This takes a long while!
##In fact takes too long.
##test.zinb <- zeroinfl(LGHT~CAPE+LR75+H5_T+S06+PMR,data=train_cls[[10]],dist="negbin",EM=TRUE)
##I could use glmmADMB, but I would have to run on my laptop
##See http://glmmadmb.r-forge.r-project.org/glmmADMB.html
##-----NOTES!!!!!-----
#
## The following doesn't work either because of convergence issues!
##Even if we increase maxit.
##i.e. test.zip <-zeroinfl(LGHT~CAPE+LR75+H5_T+S06+PMR,data=train_cls[[5]],dist="poisson",EM=TRUE,maxit=100,trace=TRUE)
#
##fit.zip <- vector("list",length=ncls)
##for (i in 1:ncls){
##    fit.zip[[i]] <- zeroinfl(LGHT~CAPE+LR75+H5_T+S06+PMR,data=train_cls[[i]],dist="poisson",EM=TRUE)
##}
#
#
##Let's use the probability of lightning from the clusters
##and then use zero-truncated regression.
##See https://stats.idre.ucla.edu/r/dae/zero-truncated-negative-binomial/
#
########### ---------- Using bestglm doesn't work ---------- ##########
########### ---------- Need to use step instead ---------- ##########
##Use bestglm package to determine optimal number of parameters
#
###First need to drop the CLS variable
##drops <- "CLS"
##keeps <- names(train)[2:7]
###DF[ , !(names(DF) %in% drops)]
##for (i in 1:ncls){
##    #train_cls[[i]]<-train_cls[[i]][,!(names(train_cls[[i]]) %in% drops)]
##    #train_cls[[i]]<-train_cls[[i]][,2:7]
##    #names(train_cls[[i]]) <- c(names(train_cls[[i]][2:5]),"y")
##    train_cls[[i]] <- subset(train_cls[[i]],select=keeps)
##    names(train_cls[[i]]) <- c(names(train_cls[[i]])[1:5],"y")
##}
##
##GPATS_bestglm <- bestglm(train_cls[[10]], family=poisson, IC="AIC")
###GPATS_bestglm <- bestglm(train_cls[[10]], family=poisson, IC="BIC")
##
##bestglm doesn't seem to be working
##null <- glm(y~1, family=poisson(),data=train_cls[[10]])
#
###Have problems when at cluster 3.
###We may be able to increase the number of iterations
###See https://stats.stackexchange.com/questions/5354/logistic-regression-model-does-not-converge
##GPATS_best_glm.nb <- vector("list",length=ncls)
##
##for (i in 1:ncls){
##    null <- glm.nb(LGHT~1,data=train_cls[[i]])
##    full <- glm.nb(LGHT ~ CAPE+LR75+H5_T+S06+PMR,data=train_cls[[i]])
##    backward <- step(full, data=train_cls[[i]],direction="backward",trace=TRUE)
##    GPATS_best_glm.nb[[i]] <- backward
##}
##    #GPATS_best_glm.nb[[i]] <- glm.nb(LGHT~1,data=train_cls[[
##
##null <- glm(LGHT~1, family=poisson(),data=train_cls[[10]])
##
###full <- glm(y ~ CAPE+LR75+H5_T+S06+PMR,family=poisson(),data=train_cls[[10]])
##full <- glm(LGHT ~ CAPE+LR75+H5_T+S06+PMR,family=poisson(),data=train_cls[[10]])
##
##backward <- step(full, data=train_cls[[10]],direction="backward",trace=FALSE)
#
##Let's try bestglm to see which variables to keep in the regression
##No it doesn't work, even when we get rid of the zeros.
#
####First need to drop the CLS variable
###drops <- "CLS"
##keeps <- names(train)[2:7]
####DF[ , !(names(DF) %in% drops)]
##for (i in 1:ncls){
##    #train_cls[[i]]<-train_cls[[i]][,!(names(train_cls[[i]]) %in% drops)]
##    #train_cls[[i]]<-train_cls[[i]][,2:7]
##    #names(train_cls[[i]]) <- c(names(train_cls[[i]][2:5]),"y")
##    train_cls[[i]] <- subset(train_cls[[i]],select=keeps)
##    names(train_cls[[i]]) <- c(names(train_cls[[i]])[1:5],"y")
##}
##
##GPATS_bestglm <- bestglm(train_cls[[10]], family=poisson, IC="AIC")
###GPATS_bestglm <- bestglm(train_cls[[10]], family=poisson, IC="BIC")
#
##for (i in 1:ncls){
###for (i in 10){
##    null.p[[i]] <- glm(LGHT~1,family=poisson,data=train_cls[[i]])
##    full.p[[i]] <- glm(LGHT ~ CAPE+LR75+H5_T+S06+PMR,family=poisson,data=train_cls[[i]])
##    #null.qp[[i]] <- glm(log(LGHT)~1,family=quasipoisson,data=train_cls[[i]])
##    #full.qp[[i]] <- glm(log(LGHT) ~ CAPE+LR75+H5_T+S06+PMR,family=quasipoisson,data=train_cls[[i]])
##    null.nb[[i]] <- glm.nb(log(LGHT)~1,data=train_cls[[i]])
##    full.nb[[i]] <- glm.nb(log(LGHT) ~ CAPE+LR75+H5_T+S06+PMR,data=train_cls[[i]])
##    ##backward.p[[i]] <- step(full[[i]],data=train_cls[[i]],direction="backward",trace=TRUE)
##    ##backward.AIC.p[[i]] <- stepAIC(full[[i]],data=train_cls[[i]],direction="backward",
##    ##                               trace=TRUE)
##    #forward.p[[i]] <- step(null.p[[i]],data=train_cls[[i]],
##    #                       scope=list(lower=null.p[[i]], upper=full.p[[i]]),
##    #                       direction="forward",
##    #                       trace=TRUE)
##    #forward.qp[[i]] <- step(null.qp[[i]],data=train_cls[[i]],
##    #                       scope=list(lower=null.qp[[i]], upper=full.qp[[i]]),
##    #                       direction="forward",
##    #                       trace=TRUE)
##    #forward.nb[[i]] <- step(null.nb[[i]],data=train_cls[[i]],
##    #                       scope=list(lower=null.nb[[i]], upper=full.nb[[i]]),
##    #                       direction="forward",
##    #                       trace=TRUE)
##
##    #forward.AIC.p[[i]] <- stepAIC(null[[i]],data=train_cls[[i]],
##    #                              scope=list(lower=null[[i]], upper=full.p[[i]]),
##    #                              direction="forward",
##    #                              trace=TRUE)
##    #GPATS_best_glm.nb[[i]] <- backward[[i]]
##    #GPATS_best_glm.p[[i]] <- backward[[i]]
##    #GPATS_best_glm.p[[i]] <- forward.p[[i]]
##    #GPATS_best_glm.qp[[i]] <- forward.qp[[i]]
##    GPATS_best_glm.nb[[i]] <- forward.nb[[i]]
##}
##
#
#
