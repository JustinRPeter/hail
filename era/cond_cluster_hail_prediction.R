#Read in SHIP/GPATS/HAIL occurrence data
#an apply a logistic regression
#See the following for references
#http://www.r-bloggers.com/how-to-perform-a-logistic-regression-in-r/

rm(list=ls())
library(ggplot2)
library(dplyr)
library(reshape)
library(leaps)
library(pscl)
library(ROCR)
library(flexclust) #Used for predicting cluster membership of new data.
library(RColorBrewer)
library(wesanderson)
library(oce)
library(ncdf4)


###############################################################################
#(1) Read the cluster prediction file output from
#cluster_lightning_prediction.R
cls.pred <- readRDS("pred_test_all_lightning_cclust_seg_cape_ex_vars_2008_2014_20clusters.rds")
ncls <- length(levels(cls.pred))
class(cls.pred) <- "numeric"

#Apply the land sea mask to the cluster prediction data


#Read the probability of lightning associated with each cluster
prob.lightning <- readRDS("prop_lightning_all_seg_cape_ex_vars_2008_2014_20clusters.rds")




raw.data <- read.table("hail_predictor.txt",
                   header=TRUE,
                   stringsAsFactors=FALSE)
#Turn the hail occurrence into a factor
#data$Hail_Occurred <- as.factor(data$Hail_Occurred)
#Use a shorter name
names(raw.data)[20] <- "HAIL"
#Replace no/yes with 0/1
raw.data$HAIL <- ifelse(raw.data$HAIL=="Yes",1,0)
#data$loggpats <- log(data$Max_GPATS)
#Transform mixing ratio units to g/kg
raw.data$PMR_SMG <- raw.data$PMR_SMG*1000.


#Extract the values where lightning was observed
#ltng.threshold <- 1000
#ltng.threshold <- 100
#ltng.threshold <- 10
#ltng.threshold <- 1 #Using 1 is indicative that convection occurred.
ltng.threshold <- 0
ltng.obs <- which(raw.data$Max_GPATS >= ltng.threshold)
#ltng.obs <- which(raw.data$HAIL == ltng.threshold)
data <- raw.data[ltng.obs,]

lightning.but.no.hail <- which(data$Max_GPATS > 0 & data$HAIL == 0)
hail.but.no.lightning <- which(data$Max_GPATS == 0 & data$HAIL == 1)

#Remove instances of when hail occurred but there was no lightning
data <- data[-hail.but.no.lightning,]

#Evaluate the length of the data set
data.length <- length(data$Time_GPATS)

#The fraction of data we will use for training
train.fraction <- 0.7

#train.data <- data[sample(nrow(data),train.fraction*data.length),]
train <- sample_frac(data,train.fraction)
test <- subset(data,!rownames(data) %in% rownames(train))


#Perform a cluster analysis of the SHIP components

#Scale the data
#dfn <- scale(data[,7:length(data)])
#dfn <- scale(data[,13:18])
#dfn <- scale(data[,14:19])
#dfn <- scale(data[,13:19])
#We only want to use CAPE, LR75, H5_TEMP, S06 and PMR in cluster calculations
#dfn <- as.data.frame(scale(data[,14:18]))
#names(dfn) <- c(names(data)[14:18])
#Using LR75,H5_TEMP,S06,PMR
dfn <- as.data.frame(scale(data[,15:18]))
names(dfn) <- c(names(data)[15:18])
#pairs(dfn)

#ncmax <- 15
ncmax <- 10
nrep <- 25
#nrep <- 50
wss <- rep(0,ncmax-1)

#ncls=6
#ncls=10
#ncls=15
ncls=20
#ncls=64
#ncls=100
#ncls=250
cat("Calculating clustering with",ncls, "clusters...","\n")
#km.model <- kmeans(dfn[,1:6],centers=5,nstart=nrep,iter.max=250)
#km.model <- kmeans(dfn[,2:6],centers=5,nstart=nrep,iter.max=250)
#km.model <- kmeans(dfn[,1:5],centers=ncls,nstart=nrep,iter.max=250)
#km.model <- kmeans(dfn[,1:5],centers=ncls,nstart=nrep,iter.max=250,algorithm='Lloyd')
km.model <-  kcca(dfn[,1:length(dfn)], k=ncls, family=kccaFamily("kmeans"))
#ncls <- length(km.model$size)
cat('\n',paste0(ncls),'Clusters\n')
cat('Cluster\tSize\tCentres\n')
#ind gives the order of the clusters based on the ordering of CAPE
#ind <- order(km.model$centers[,1])
ind <- order(km.model@centers[,1])
#for (i in 1:ncls) cat(i,km.model$size[ind[i]],km.model$centers[ind[i],],'\n',sep='\t')
for (i in 1:ncls) cat(i,km.model@clusinfo$size[ind[i]],km.model@centers[ind[i],],'\n',sep='\t')
#for (i in 1:3) cat(i,km3$size[i],km3$centers[i,],'\n',sep='\t')

#Plot a barplot of the relative contribution of the clusters
cls.conts <- c()
for (i in 1:ncls){
    #cls.conts[i] <- km.model$size[ind[i]]/sum(km.model$size)
    cls.conts[i] <- km.model@clusinfo$size[ind[i]]/sum(km.model@clusinfo$size)
}
pdf(paste0("cluster_contribution_",ncls,"_clusters.pdf"),width=10,height=10,paper="special")
barplot(cls.conts*100.,names.arg=c(1:ncls),
            xlab="Cluster ordered by CAPE",
            ylab="% contribution")
dev.off()

#Use ggplot2 for plotting

#cluster <- factor(km.model$cluster)
#Order the clusters according to CAPE and name them
#1 through ncls for plotting convenience
#This also means we can reference them naturally
#cluster <- factor(km.model$cluster,
cluster <- factor(km.model@cluster,
                  levels=ind,ordered=TRUE,
                  labels=as.character(1:ncls))

df <- cbind(data[,13:20],cluster)
dfg1 <- melt(df,id.vars=c("SHIP_Max_GPATS","Max_GPATS","cluster"),
               measure.vars=c("MUCAPE_SMG","LR75_SMG","H5_TEMP_SMG","S06_SMG","PMR_SMG"))

labels <- c(SHIP_Max_GPATS="SHIP",
            Max_GPATS="Lightning counts",
            MUCAPE_SMG="MUCAPE [J/kg]", 
            LR75_SMG="LR75 [C/km]",
            H5_TEMP_SMG= "T_H500 [C]",
            S06_SMG="S06",
            PMR_SMG="w [kg/kg]")

#Plot the SHIP as a function of the cluster variables
p1<-ggplot(dfg1,aes(x=SHIP_Max_GPATS,y=value))+
       xlab("SHIP")+ylab("")+
       geom_point(aes(color=cluster))+
       facet_wrap(~variable,ncol=2,scales="free",labeller=labeller(variable=labels[3:7]))+
       scale_colour_manual(values=rev(oce.colorsViridis(length(ind))))+
       guides(colour=guide_legend(override.aes=list(size=5)))
#plot(p1)

##Some plotting parameters
#xsize=ysize=15
##psize=48
#png("SHIP_vs_cls.png",width=xsize,height=ysize,units="cm",res=600)
#print(p1)
#dev.off()

#Plot the Lightning as a function of the cluster variables
p2<-ggplot(dfg1,aes(x=Max_GPATS,y=value))+
       xlab("Lightning counts / hour")+ylab("")+
       geom_point(aes(color=cluster))+
       facet_wrap(~variable,ncol=2,scales="free",labeller=labeller(variable=labels[3:7]))+
       scale_colour_manual(values=rev(oce.colorsViridis(length(ind))))+
       guides(colour=guide_legend(override.aes=list(size=5)))+
       scale_x_log10()
#plot(p2)
#png("GPATS_vs_cls.png",xsize,height=ysize,units="cm",res=600)
#print(p2)
#dev.off()


#Now plot histograms of each of the variables by cluster

#Use the log(Max_GPATS) to differentiate
log.gpats <- log10(data$Max_GPATS)
#temp.minus <- -1*data$H5_TEMP_SMG
dflt <- cbind(data[,c(13,14,15,16,17,18)],log.gpats,cluster)
#names(dflt) <- names(data)[13:19]

llabels <- c(SHIP_Max_GPATS="SHIP",
            log.gpats="log(Lightning counts)",
            MUCAPE_SMG="MUCAPE [J/kg]", 
            LR75_SMG="LR75 [C/km]",
            H5_TEMP_SMG= "T_H500 [C]",
            S06_SMG="S06",
            PMR_SMG="w [kg/kg]")


#dfg2 <- melt(df,id.vars=c("cluster"),
dfg2 <- melt(dflt,id.vars=c("cluster"),
            #measure.vars=c("SHIP_Max_GPATS",
            #               "MUCAPE_SMG","LR75_SMG",
            #               "H5_TEMP_SMG","S06_SMG","PMR_SMG"))
            measure.vars=c("SHIP_Max_GPATS","log.gpats",
                           "MUCAPE_SMG","LR75_SMG",
                           "H5_TEMP_SMG","S06_SMG","PMR_SMG"))
#p3 <- ggplot(dfg,aes(x=SHIP_Max_GPATS,y=value))+
p3 <- ggplot(dfg2,aes(x=value))+
      xlab("")+
      #stat_bin(aes(y=..density..,color=cluster),geom="step",size=1.5,position="identity")+
      geom_freqpoly(aes(y=..density..,color=cluster),size=1.5)+
      #stat_bin(aes(y=..count..,color=cluster),geom="step",size=1.5)+
      facet_wrap(~variable,ncol=2,scale="free",labeller=labeller(variable=llabels))+
      #scale_colour_manual(values=wes_palette(n=length(ind),name="Zissou"))+
      scale_colour_manual(values=rev(oce.colorsViridis(length(ind))))+
      guides(colour=guide_legend(override.aes=list(size=5)))

#plot(p3)
#png("hists_vars_by_cls.png",width=xsize,height=ysize,units="cm",res=600)
#print(p3)
#dev.off()

#Which clusters produce hail/nohail lightning/no lightning?
hail <- which(df$HAIL ==1)
no.hail <- which(df$HAIL == 0)
length(hail)
length(no.hail)
tot.len.hail<-length(hail)+length(no.hail)
no.days.hail <- length(unique(data$Date_GPATS[hail]))
no.days.no.hail <- length(unique(data$Date_ERA[no.hail]))

lightning <- which(df$Max_GPATS > 0)
no.lightning <- which(df$Max_GPATS == 0)
length(lightning)
length(no.lightning)
tot.len.lightning <- length(lightning)+length(no.lightning)

prop.hail <- c()
days.hail <- c()
prop.no.hail <- c()
days.no.hail <- c()
prop.lightning <- c()
prop.lightning.strikes <- c()
days.lightning <- c()
prop.no.lightning <- c()
days.no.lightning <- c()

tot.lightning.strikes <- sum(data$Max_GPATS[lightning])
tot.hours <- tot.len.lightning
lstrks.per.hour <- tot.lightning.strikes/tot.hours

#The number of years
nyrs <- as.numeric(ceiling(difftime(max(data$Date_ERA),
                         min(data$Date_ERA))/365.25))

cat('\n Proportion of each cluster on which hail occurs \n')
for (i in 1:ncls){
    #cat('Cluster ',i,length(which(df$cluster[hail]==i))/length(hail),'\n')
    #prop.hail[i] <- length(which(df$cluster[hail]==i))/length(hail)
    #prop.lightning[i] <- length(which(df$cluster[lightning]==i))/
    #                       length(lightning)
    cat('Cluster ',i,"Hail",
        length(which(df$cluster[hail]==i))/tot.len.hail,'\n')
    cat('Cluster',i,"GPATS",
        length(which(df$cluster[lightning]==i))/tot.len.lightning,'\n')
    #prop.hail[i] <- length(which(df$cluster[hail]==i))/tot.len.hail
    prop.hail[i] <- length(which(df$cluster[hail]==i))/
                    length(which(df$cluster==i))
    #days.hail[i] <- length(which(df$cluster[hail]==i))/nyrs
    #days.hail[i] <- length(unique(data$Date_ERA[which(df$cluster[hail]==i)]))/
                    #nyrs
                    #length(unique(data
    #prop.lightning[i] <- length(which(df$cluster[lightning]==i))/
    #                       tot.len.lightning
    prop.lightning[i] <- length(which(df$cluster[lightning]==i))/
                         length(which(df$cluster==i))
    prop.lightning.strikes[i] <- prop.lightning[i]*(sum(data$Max_GPATS[
                                     lightning[which(df$cluster[lightning]==i)]])/
                                 #tot.lightning.strikes
                                 #tot.len.lightning
                                 #length(lightning))
                                 #length(which(df$cluster==i)))
                                 length(which(df$cluster[lightning]==i)))
                                 #(nyrs*365.25*32)
    #prop.lightning.strikes[i] <- length(which(df$cluster[lightning] == i))
    #days.lightning[i] <- length(which(df$cluster[lightning]==i))/nyrs
}
#plot(df$cluster[hail])
#barplot(prop.hail*100.,names=as.character(1:ncls),
#        xlab="Cluster",ylab="%")


cat('\n Proportion of each cluster on which NO hail/lightning occurs \n')
for (i in 1:ncls){
    #cat('Cluster ',i,length(which(df$cluster[no.hail]==i))/length(no.hail),'\n')
    #prop.no.hail[i] <- length(which(df$cluster[no.hail]==i))/length(no.hail)
    #prop.no.lightning[i] <- length(which(df$cluster[no.lightning]==i))/
    #                        length(no.lightning)
    cat('Cluster ',i,"No Hail",
        length(which(df$cluster[no.hail]==i))/tot.len.hail,'\n')
    cat('Cluster',i,"No GPATS",
        length(which(df$cluster[no.lightning]==i))/tot.len.lightning,'\n')
    
    #prop.no.hail[i] <- length(which(df$cluster[no.hail]==i))/tot.len.hail
    prop.no.hail[i] <- length(which(df$cluster[no.hail]==i))/
                       length(which(df$cluster[no.hail]==i))
    #days.no.hail[i] <- length(unique(data$Date_ERA[which(df$cluster[no.hail]==i)]))/nyrs
    prop.no.lightning[i] <- length(which(df$cluster[no.lightning]==i))/
                             length(which(df$cluster==i))
                            #tot.len.lightning

}
#plot(df$cluster[no.hail])
#barplot(prop.no.hail*100.,names=as.character(1:ncls),
#        xlab="Cluster",ylab="%")

#cl <- as.character(1:ncls)
cl <- factor(1:ncls)
hd <- data.frame(prop.hail,prop.no.hail,
                 #days.hail,days.no.hail,
                 prop.lightning,prop.no.lightning,
                 prop.lightning.strikes,cl)
hd.l <- melt(hd,id.vars="cl")
phd <- ggplot(hd.l,aes(cl,value,fill=variable))+
       geom_bar(stat="identity",position="dodge")+
       xlab("Cluster")+ylab("%")+
       scale_fill_brewer("",
                         labels=c("Hail","No hail","Lightning","No lightning","Lightning strikes"),
                        palette="Set1")
#plot(phd)
#png("hail_lightning_contribution_by_cls.png",width=xsize,height=ysize,units="cm",res=600)
#print(phd)
#dev.off()


#Now use the flexclust package to predict cluster association
#of data

#Read the ERA data for 20008-2014 so that we can 
#compare the probability of lightning

#fnm_era <- '/home/u1066578/data/era_output/00_18/2008_2014/analysis/era_profiles_2008_2014_00_18.nc'
fnm_era <- '/home/u1066578/data/era_output/00_18/1979_2014/analysis/era_profiles_1979_2014_00_18.nc'

#Open the netcdf files
nc_era   <- nc_open(fnm_era)
print(paste("SHARPPy profile file has",nc_era$nvars,"variables"))

#Get the dimensions
nlat_era    <- nc_era$dim$lat$len
nlon_era    <- nc_era$dim$lon$len
ntime_era   <- nc_era$dim$time$len

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

ex.lat <- which(lat >= min_lat & lat <= max_lat)
ex.lon <- which(lon >= min_lon & lon <= max_lon)
lat <- lat[ex.lat]
lon <- lon[ex.lon]
mucape   <- mucapel[ex.lon,ex.lat,]
lr75     <- lr75l[ex.lon,ex.lat,]
h5_temp  <- h5_templ[ex.lon,ex.lat,]
s06      <- s06l[ex.lon,ex.lat,]
pmr      <- pmrl[ex.lon,ex.lat,]



#Apply the land-sea mask
fnm_lsmask <- '/home/u1066578/data/eraint/lsmask/era_interim_lsmask.nc'
nc_lsmask <- nc_open(fnm_lsmask)
v8 <- nc_lsmask$var[["lsm"]]
lsmask <- ncvar_get(nc_lsmask,v8)
lat_lsmask <- nc_lsmask$dim$lat$vals
lon_lsmask <- nc_lsmask$dim$lon$vals
ex.lat_lsmask <- which(lat_lsmask >= min_lat & lat_lsmask <= max_lat)
ex.lon_lsmask <- which(lon_lsmask >= min_lon & lon_lsmask <= max_lon)
lat_lsmask <- lat_lsmask[ex.lat_lsmask]
lon_lsmask <- lon_lsmask[ex.lon_lsmask]

lsmask <- lsmask[ex.lon_lsmask,ex.lat_lsmask]

for (i in 1:length(time)){
    mucape[,,i][which(lsmask==0)] <- NA
    lr75[,,i][which(lsmask==0)] <- NA
    h5_temp[,,i][which(lsmask==0)] <- NA
    s06[,,i][which(lsmask==0)] <- NA
    pmr[,,i][which(lsmask==0)] <- NA
}


nlon <- dim(mucape)[1]
nlat <- dim(mucape)[2]
ntime <- dim(mucape)[3]

#Put all these variables in a list
#era_data <- list(MUCAPE_SMG=as.vector(mucape),LR75_SMG=as.vector(lr75),
#                 H5_TEMP_SMG=as.vector(h5_temp),S06_SMG=as.vector(s06),
#                 PMR_SMG=as.vector(pmr))

era_data <- list(LR75_SMG=as.vector(lr75),
                 H5_TEMP_SMG=as.vector(h5_temp),S06_SMG=as.vector(s06),
                 PMR_SMG=as.vector(pmr))

era_data <- as.data.frame(era_data)



#era_data_n <- as.data.frame(scale(era_data))
source("colScale.R")
cat("Scaling ERA data before cluster prediction: ","\n")
era_data_n <- colScale(as.matrix(era_data))
#era_data_n <- scale(era_data)
era_data_n <- as.data.frame(era_data_n)
names(era_data_n) <- names(era_data)

#era_data <- as.data.frame(as.vector(mucape),as.vector(lr75),
#                  as.vector(h5_temp),as.vector(s06),
#                  as.vector(pmr))


#Clean up some memory
rm(list=c("mucapel","lr75l","h5_templ","s06l","pmrl"))
rm(list="era_data")
gc()


#Convert the kmeans model to kcca data type for processingA

cat("Predicting cluster association for all data: ", "\n")
#cl.model <- as.kcca(km.model,dfn)
#cl.model <- as.kcca(km.model,dfn)
#pred_train <- predict(cl.model)
#pred_test <- predict(cl.model,newdata=era_data[1:184830,]) #1 month
#pred_test <- predict(cl.model,newdata=era_data_n[1:2248765,]) #1 year
#pred_test <- predict(cl.model,newdata=era_data_n[1:8995060,]) #1 year
#pred_test <- predict(km.model,newdata=era_data_n[1:8995060,]) #1 year
#pred_test <- predict(km.model,newdata=era_data_n[1:(nlon*nlat*4*366),]) #1 year
pred_test <- predict(km.model,newdata=era_data_n) #All years

#We need to put the clusters in the same order as for the hail data
pred_test.f <- factor(pred_test,levels=ind,ordered=TRUE,labels=as.character(1:ncls))

#dim(pred_test) <- c(nlon,nlat,366*4)
#dim(pred_test.f) <- c(nlon,nlat,366*4)
dim(pred_test) <- c(nlon,nlat,ntime)
dim(pred_test.f) <- c(nlon,nlat,ntime)

#fptout <- paste0("pred_test_2008_",ncls,"clusters.rds")
#fhpout <- paste0("prop_hail_lightning_",ncls,"clusters.rds")
#fptout <- paste0("pred_test_all_",ncls,"clusters.rds")
#fptout <- paste0("pred_test_all_ex_vars",ncls,"clusters.rds")
#fhpout <- paste0("prop_hail_lightning_all_",ncls,"clusters.rds")
fptout <- paste0("pred_test_all_ex_vars_1979_2014_",ncls,"clusters.rds")
fhpout <- paste0("prop_hail_lightning_all_1979_2014_",ncls,"clusters.rds")
saveRDS(pred_test.f,fptout)
saveRDS(hd,fhpout)

