#Perform cluster analysis
#And predict the cluster membership of the remaining data
# *****
#Show where you need to make changes depending
#on which how you perform cluster based on CAPE
#(1) Using all CAPE values. No transformation
#(2) Use min. CAPE of 100. log transform
#(3) Split CAPE into in equal sized sections

#Remove NA's
#data <- data[complete.cases(data),]
#on the SHIP, SHIP components and GPATS data
#See the following for references
#http://www.r-bloggers.com/how-to-perform-a-logistic-regression-in-r/

rm(list=ls())
gc()
library(ggplot2)
library(dplyr)
library(reshape)
library(leaps)
library(pscl)
library(Hmisc) #for Ecdf function
#library(flexclust) #Used for predicting cluster membership of new data.
library(cclust)    #flexclust also has a cclust function
library(RColorBrewer)
#library(wesanderson)
library(oce)
library(ncdf4)

source("colScale.R")


#raw.data <- read.table("lightning_predictor_0.25.txt",
#                   header=TRUE,
#                   na.strings=c("", "NA"),
#                   stringsAsFactors=FALSE)

#raw.data <- readRDS("lightning_predictor_0.25_ship.rds")
raw.data <- readRDS("lightning_predictor_0.45.rds")

#Add a lightning indicator factor to the data
linds <- which(raw.data$GPATS > 0)
ltng <- vector("numeric",length=length(raw.data$Lat))
ltng[linds] <- 1
ltng[-c(linds)] <- 0
raw.data$Lightning <- ltng

#Convert parcel mixing ratio to g/kg.
raw.data$PMR <- raw.data$PMR*1000.

#Some exploratory analysis
#hist(strptime(raw.data$Time,format="%H:%M:%S"),breaks="hours")


#Extract the values where lightning was observed
#ltng.threshold <- 1000
#ltng.threshold <- 100
#ltng.threshold <- 10
#ltng.threshold <- 1 #Using 1 is indicative that convection occurred.
ltng.threshold <- 0
ltng.obs <- which(raw.data$GPATS >= ltng.threshold)
#ltng.obs <- which(raw.data$HAIL == ltng.threshold)
data <- raw.data[ltng.obs,]

#Remove NA's
data <- data[complete.cases(data),]

#Evaluate the length of the data set
data.length <- length(data$Time)

#The fraction of data we will use for training
train.fraction <- 0.7

#train.data <- data[sample(nrow(data),train.fraction*data.length),]
train <- sample_frac(data,train.fraction)
test <- subset(data,!rownames(data) %in% rownames(train))

#data <- sample_frac(data,0.001)
#data <- sample_frac(data,0.5)
data <- sample_frac(data,1.0)

#  *****
#Use the log of CAPE
data <- data[data$MUCAPE > 0,]
#data <- data[data$MUCAPE > 100,]
#data <- data[data$MUCAPE > 1000,]
#The following to try log cape
#data$MUCAPE <- log(data$MUCAPE)



#Perform a cluster analysis of the SHIP components

#Scale the data
#dfn <- scale(data[,7:length(data)])
#dfn <- scale(data[,13:18])
#dfn <- scale(data[,14:19])
#dfn <- scale(data[,13:19])
#We only want to use CAPE, LR75, H5_TEMP, S06 and PMR in cluster calculations
#This is using all of the predictors
#dfm <- colScale(as.matrix(data[,6:10])) #used for clclust. needs matrix
#Following just using MNUCAPE, S06, PMR
#ex_vars <- c(6,9,10)
#ex_vars <- c(7,9,10) #LR75, S06,PMR
ex_vars <- c(7,8,9,10) #LR75, H5_TEMP, S06,PMR
#ex_vars <- c(7,8,9,10,11) #LR75, H5_TEMP, S06, PMR, CAP
dfm <- colScale(as.matrix(data[,ex_vars])) #used for clclust. needs matrix


##Get the scale and centering vectors
cnt <- attr(dfm,"scaled:center")
scl <- attr(dfm,"scaled:scale")
cat(cnt,scl,"\n")

#dfn <- as.data.frame(scale(data[,6:10]))
dfn <- as.data.frame(dfm)
#names(dfn) <- c(names(data)[6:10])
names(dfn) <- c(names(data)[ex_vars])
#pairs(dfn)

#Used for normal kmeans
#ncmax <- 15
#ncmax <- 10
#nrep <- 25
#nrep <- 50
#wss <- rep(0,ncmax-1)

#Used for cclust
#kmeans specification
#nrep <- 25
#niter <- 500


#Use these for cclust
#ncls=6
#ncls=10
ncls=20
#ncls=50
#niter=500
#ncls=50
#niter=1000
#ncls=100
niter=2000

## *****
##To be used if specifying the cluster centres
#dfn.split <- split(dfn, cut2(dfn$MUCAPE,g=ncls))
##rand.centre <- vector("list",length=dim(dfn)[2])
#rand.centre <- vector("list",length=ncls)
##for (i in 1:nc){
#for (i in 1:ncls){
#    rand.centre[[i]] <-
#    dfn.split[[i]][sample(nrow(dfn.split[[i]]),1),]
#}
#clus.centres <- matrix(unlist(rand.centre),nrow=ncls,byrow=T)


#-----Using kcca
##km.model <- kmeans(dfn[,1:6],centers=5,nstart=nrep,iter.max=250)
##km.model <- kmeans(dfn[,2:6],centers=5,nstart=nrep,iter.max=250)
##km.model <- kmeans(dfn[,1:5],centers=ncls,nstart=nrep,iter.max=250)
##km.model <- kmeans(dfn[,1:5],centers=ncls,nstart=nrep,iter.max=250,algorithm='Lloyd')
#km.model <-  kcca(dfn[,1:5], k=ncls, family=kccaFamily("kmeans"))

#Using cclust from flexclust package
#Tests show that this is same as cclust from cclust
#but want to use cclust so that we can use the predict method it offers.
#km.model <-  cclust(dfm[,1:5], centers=ncls, iter.max=niter,
#                     dist="euclidean",method="kmeans",verbose=TRUE)
#
#
##ncls <- length(km.model$size)
#cat('\n',paste0(ncls),'Clusters\n')
#cat('Cluster\tSize\tCentres\n')
##ind gives the order of the clusters based on the ordering of CAPE
##ind <- order(km.model$centers[,1])
#ind <- order(km.model@centers[,1])
##for (i in 1:ncls) cat(i,km.model$size[ind[i]],km.model$centers[ind[i],],'\n',sep='\t')
#for (i in 1:ncls) cat(i,km.model@clusinfo$size[ind[i]],km.model@centers[ind[i],],'\n',sep='\t')
##for (i in 1:3) cat(i,km3$size[i],km3$centers[i,],'\n',sep='\t')
#
##Plot a barplot of the relative contribution of the clusters
#cls.conts <- c()
#for (i in 1:ncls){
#    #cls.conts[i] <- km.model$size[ind[i]]/sum(km.model$size)
#    cls.conts[i] <- km.model@clusinfo$size[ind[i]]/sum(km.model@clusinfo$size)
#}

#-----Using cclust
# *****
#Change the following if specifying the cluster centres
km <- cclust(dfm,ncls,niter,verbose=T,method="kmeans")
#km <- cclust(dfm,clus.centres,niter,verbose=T,method="kmeans")
ind <- order(km$centers[,1])
cat('Cluster\tSize\tCentres\n')
for (i in 1:ncls) cat(i,km$size[ind[i]],km$centers[ind[i],],'\n',sep='\t')
cls.conts <- c()
for (i in 1:ncls){
    cls.conts[i] <- km$size[ind[i]]/sum(km$size)
}
# *****
#Chnage the following depending on CAPE
#sink(paste0("Cluster_information_",ncls,".txt"),type='output',append=FALSE)
sink(paste0("Cluster_information_seg_cape",ncls,".txt"),type='output',append=FALSE)
cat('Cluster\tSize\tCentres\n')
for (i in 1:ncls) cat(i,km$size[ind[i]],km$centers[ind[i],],'\n',sep='\t')

cat('Normalising parameters\n')
cat('\n')
cat('Scale\tCenter\n')
cat(scl,'\n',sep='\t')
cat(cnt,'\n',sep='\t')

cat('\n')
cat('Normalised cluster centres\n')
ncls_c <- as.data.frame(km$centers[ind,])
write.table(ncls_c)
cat('\n')
cat('Centres in original units\n')
scl_cls_c <- sweep(ncls_c,MARGIN=2,scl,`*`)
cls_c <- sweep(scl_cls_c,MARGIN=2,cnt,'+')
write.table(cls_c)
cat('\n')

sink()



# *****
#pdf(paste0("lightning_cluster_contribution_",ncls,"_clusters.pdf"),width=10,height=10,paper="special")
pdf(paste0("lightning_cluster_contribution_seg_cape",ncls,"_clusters.pdf"),width=10,height=10,paper="special")
par(cex=1.25)
barplot(cls.conts*100.,names.arg=c(1:ncls),
            xlab="Cluster ordered by CAPE",
            #ylab="% contribution",cex.axis=1.25,cex.lab=1.25,cex.names=1.25)
            ylab="Percent contribution")
dev.off()

#Use ggplot2 for plotting

#cluster <- factor(km.model$cluster)
#Order the clusters according to CAPE and name them
#1 through ncls for plotting convenience
#This also means we can reference them naturally
cluster <- factor(km$cluster,
#cluster <- factor(km.model@cluster,
                  levels=ind,ordered=TRUE,
                  labels=as.character(1:ncls))

df <- cbind(data[,5:12],cluster)
dfg1 <- melt(df,id.vars=c("SHIP","GPATS","cluster"),
               measure.vars=c("MUCAPE","LR75","H5_TEMP","S06","PMR"))

labels <- c(SHIP="SHIP",
            GPATS="Lightning counts",
            MUCAPE="MUCAPE [J/kg]", 
            LR75="LR75 [C/km]",
            H5_TEMP= "T_H500 [C]",
            S06="S06",
            PMR="w [kg/kg]")

#Plot the SHIP as a function of the cluster variables
p1<-ggplot(dfg1,aes(x=SHIP,y=value))+
       xlab("SHIP")+ylab("")+
       geom_point(aes(color=cluster))+
       facet_wrap(~variable,ncol=2,
                  scales="free",
                  strip.position = "left",
                  #switch = "y",
                  labeller=labeller(variable=labels[3:7]))+
       scale_colour_manual(values=rev(oce.colorsViridis(length(ind))))+
       guides(colour=guide_legend(override.aes=list(size=5)))+
       theme(axis.text.x=element_text(face="bold",size=12),
             axis.text.y=element_text(face="bold",size=12),
             axis.title.x=element_text(face="bold",size=14),
             strip.text.y=element_text(face="bold",size=14))
       #theme(text=element_text(face="bold",size=12),
       #      strip.text.x=element_text(size=14))
#plot(p1)

##Some plotting parameters
xsize=ysize=10
##psize=48
#png("SHIP_vs_cls.png",width=xsize,height=ysize,units="cm",res=600)
#print(p1)
#dev.off()

#Plot the Lightning as a function of the cluster variables
p2<-ggplot(dfg1,aes(x=GPATS,y=value))+
       xlab("Lightning counts / hour")+ylab("")+
       geom_point(aes(color=cluster))+
       facet_wrap(~variable,ncol=2,scales="free",
                  labeller=labeller(variable=labels[3:7]),
                  switch="y")+
       scale_colour_manual(values=rev(oce.colorsViridis(length(ind))))+
       guides(colour=guide_legend(override.aes=list(size=5)))+
       scale_x_log10()+
       theme(axis.text.x=element_text(face="bold",size=12),
             axis.text.y=element_text(face="bold",size=12),
             axis.title.x=element_text(face="bold",size=14),
             strip.text.y=element_text(face="bold",size=14))
#plot(p2)
#png("GPATS_vs_cls.png",xsize,height=ysize,units="cm",res=600)
#print(p2)
#dev.off()


#Now plot histograms of each of the variables by cluster

#Use the log(Max_GPATS) to differentiate
log.gpats <- log10(data$GPATS)
#temp.minus <- -1*data$H5_TEMP_SMG
dflt <- cbind(data[,5:10],log.gpats,cluster)
#names(dflt) <- names(data)[13:19]

llabels <- c(SHIP="SHIP",
            log.gpats="log(GPATS)",
            MUCAPE="MUCAPE [J/kg]", 
            LR75="LR75 [C/km]",
            H5_TEMP= "T_H500 [C]",
            S06="S06",
            PMR="w [kg/kg]")


#dfg2 <- melt(df,id.vars=c("cluster"),
dfg2 <- melt(dflt,id.vars=c("cluster"),
            #measure.vars=c("SHIP_Max_GPATS",
            #               "MUCAPE_SMG","LR75_SMG",
            #               "H5_TEMP_SMG","S06_SMG","PMR_SMG"))
            measure.vars=c("SHIP","log.gpats",
                           "MUCAPE","LR75",
                           "H5_TEMP","S06","PMR"))
#p3 <- ggplot(dfg,aes(x=SHIP_Max_GPATS,y=value))+
p3 <- ggplot(dfg2,aes(x=value))+
      xlab("")+
      #stat_bin(aes(y=..density..,color=cluster),geom="step",size=1.5,position="identity")+
      geom_freqpoly(aes(y=..density..,color=cluster),size=1.5)+
      #stat_bin(aes(y=..count..,color=cluster),geom="step",size=1.5)+
      facet_wrap(~variable,ncol=2,scale="free",
                 labeller=labeller(variable=llabels),
                 switch="y")+
      #scale_colour_manual(values=wes_palette(n=length(ind),name="Zissou"))+
      scale_colour_manual(values=rev(oce.colorsViridis(length(ind))))+
      guides(colour=guide_legend(override.aes=list(size=5)))+
      theme(axis.text.x=element_text(face="bold",size=12),
            axis.text.y=element_text(face="bold",size=12),
            axis.title.x=element_text(face="bold",size=14),
            strip.text.y=element_text(face="bold",size=12))


#plot(p3)
# *****
#png("hists_vars_by_cls.png",width=xsize,height=ysize,units="in",res=300)
png("hists_vars_by_cls_seg_cape.png",width=xsize,height=ysize,units="in",res=300)
print(p3)
dev.off()

#Which clusters produce lightning/no lightning?
lightning <- which(df$GPATS > 0)
no.lightning <- which(df$GPATS == 0)
length(lightning)
length(no.lightning)
tot.len.lightning <- length(lightning)+length(no.lightning)

prop.lightning <- c()
prop.lightning.strikes <- c()
avg.lightning.per.cluster <- c()
lightning.per.cluster <- c()
#days.lightning <- c()
days.lightning <- vector("list",length=ncls)
days.lightning.per.year <- c()
prop.no.lightning <- c()
days.no.lightning <- c()

tot.lightning.strikes <- sum(data$GPATS[lightning])
tot.hours <- tot.len.lightning
lstrks.per.hour <- tot.lightning.strikes/tot.hours

#The number of years
nyrs <- as.numeric(ceiling(difftime(max(data$Date),
                         min(data$Date))/365.25))

#prop.lightning - the proportion that a cluster produces lightning.
#prop.lightning.strikes - the proportion of total lightning strikes that a
#                         cluster produces.
#avg. lightning per cluster - the average number of lightning strikes a
#                             cluster produces
#lightning.per.cluster - the probability of a cluster producing lightning
#                        multiplied by the average lightning per cluster.
#days.lightning.per.year - the number of days that a cluster is responsible
#                          for lightning at any location in australia
#tot.prob.lightning - the above two multiplied
#i.e. the probability that a cluster produces

cat('\n Proportion of each cluster on which lightning occurs \n')
for (i in 1:ncls){
    cat('Cluster',i,"GPATS",
        #length(which(df$cluster[lightning]==i))/tot.len.lightning,'\n')
        length(which(df$cluster[lightning]==i))/length(which(df$cluster==i)),'\n')
    #prop.hail[i] <- length(which(df$cluster[hail]==i))/tot.len.hail
    #days.hail[i] <- length(which(df$cluster[hail]==i))/nyrs
    #days.hail[i] <- length(unique(data$Date_ERA[which(df$cluster[hail]==i)]))/nyrs
    #prop.lightning[i] <- length(which(df$cluster[lightning]==i))/
    #                       tot.len.lightning
    prop.lightning[i] <- length(which(df$cluster[lightning]==i))/
                         length(which(df$cluster==i))
    #prop.lightning.strikes[i] <- prop.lightning[i]*(sum(data$GPATS[
    prop.lightning.strikes[i] <- sum(data$GPATS[
                                     lightning[which(df$cluster[lightning]==i)]])/
                                 tot.lightning.strikes
                                 #tot.len.lightning
                                 #length(lightning))
                                 #(nyrs*365.25*32)
    avg.lightning.per.cluster[i] <- sum(data$GPATS[lightning[which(df$cluster[lightning]==i)]])/
                                    length(which(df$cluster[lightning]==i))
    #avg.lightning.per.cluster[i] <- sum(data$GPATS[which(df$cluster==i)])/
    #                                length(which(df$cluster==i))
    lightning.per.cluster[i] <- prop.lightning[i]*avg.lightning.per.cluster[i]
    #prop.lightning.strikes[i] <- length(which(df$cluster[lightning] == i))
    #days.lightning[i] <- length(which(df$cluster[lightning]==i))/nyrs
    days.lightning[[i]] <- which(df$cluster[lightning]==i)
    days.lightning.per.year[i] <- length(unique(data$Date[days.lightning[[i]]]))/nyrs
}
#plot(df$cluster[hail])
#barplot(prop.hail*100.,names=as.character(1:ncls),
#        xlab="Cluster",ylab="%")

#Plot some of the lightning characteristics by cluster
# *****
#png("prop_lightning_cluster.png",width=xsize,height=ysize,units="in",res=300)
png("prop_lightning_cluster_seg_cape.png",width=xsize,height=ysize,units="in",res=300)
barplot(prop.lightning*100.,names=as.character(1:ncls),
        xlab="Cluster",ylab="%",
        main="Lightning active percentage by cluster",
        cex.axis=1.5,cex.names=1.5,cex.lab=1.5,cex.main=1.5)
dev.off()
# *****
#png("prop_lightning_strikes.png",width=xsize,height=ysize,units="in",res=300)
png("prop_lightning_strikes_seg_cape.png",width=xsize,height=ysize,units="in",res=300)
barplot(prop.lightning.strikes*100,names=as.character(1:ncls),
        xlab="Cluster",ylab="%",
        main="Proportion of total lightning strikes by cluster",
        cex.axis=1.5,cex.names=1.5,cex.lab=1.5,cex.main=1.5)
dev.off()
png("avg_lightning_strikes_by_cluster.png",width=xsize,height=ysize,units="in",res=300)
barplot(avg.lightning.per.cluster,names=as.character(1:ncls),
        xlab="Cluster",ylab="# of strikes",
        main="Mean lightning strikes per lightning active cluster",
        cex.axis=1.5,cex.names=1.5,cex.lab=1.5,cex.main=1.5)
dev.off()
# *****
#png("avg_lightning_strikes_per_cluster.png",width=xsize,height=ysize,units="in",res=300)
png("avg_lightning_strikes_per_cluster.png",width=xsize,height=ysize,units="in",res=300)
barplot(lightning.per.cluster,names=as.character(1:ncls),
        xlab="Cluster",ylab="# of strikes",
        main="Mean lightning per cluster",
        cex.axis=1.5,cex.names=1.5,cex.lab=1.5,cex.main=1.5)
dev.off()



cat('\n Proportion of each cluster on which NO lightning occurs \n')
for (i in 1:ncls){
    cat('Cluster',i,"No GPATS",
        length(which(df$cluster[no.lightning]==i))/length(which(df$cluster==i)),'\n')
    prop.no.lightning[i] <- length(which(df$cluster[no.lightning]==i))/
                            #tot.len.lightning
                            length(which(df$cluster==i))

}
#plot(df$cluster[no.hail])
#barplot(prop.no.hail*100.,names=as.character(1:ncls),
#        xlab="Cluster",ylab="%")

#cl <- as.character(1:ncls)
cl <- factor(1:ncls)
hd <- data.frame(prop.lightning,prop.no.lightning,
                 prop.lightning.strikes,
                 avg.lightning.per.cluster,
                 lightning.per.cluster,
                 days.lightning.per.year,cl)
hd.l <- melt(hd,id.vars="cl")
#phd <- ggplot(hd.l,aes(cl,value,fill=variable))+
phd <- ggplot(hd.l,aes(cl,value,fill=c("prop.lightning")))+
       geom_bar(stat="identity",position="dodge")+
       xlab("Cluster")+ylab("%")+
       scale_fill_brewer("",
                         labels=c("Lightning","No lightning","Lightning strikes"),
                        palette="Set1")
#plot(phd)
#png("hail_lightning_contribution_by_cls.png",width=xsize,height=ysize,units="cm",res=600)
#print(phd)
#dev.off()



#Now use the flexclust package to predict cluster association
#of data

#Read the ERA data for 20008-2014 so that we can 
#compare the probability of lightning

#fnm_era <- '/home/jpeter/justinp/rfiles/suncorp/data/eraint/sharpy_profiles/00_18/2008_2014/era_profiles_2008_2014_00_18.nc'
#fnm_era <- '/home/u1066578/data/era_output/00_18/2008_2014/analysis/era_profiles_2008_2014_00_18.nc'
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

v9 <- nc_era$var[["cap"]]
capl <- ncvar_get(nc_era,v9)


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
cap      <- capl[ex.lon,ex.lat,]

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

#for (i in 1:length(time)){
#    mucape[,,i][which(lsmask==0)] <- NA
#    lr75[,,i][which(lsmask==0)] <- NA
#    h5_temp[,,i][which(lsmask==0)] <- NA
#    s06[,,i][which(lsmask==0)] <- NA
#    pmr[,,i][which(lsmask==0)] <- NA
#}


nlon <- dim(mucape)[1]
nlat <- dim(mucape)[2]
ntime <- dim(mucape)[3]



#Put all these variables in a list
#era_data <- list(MUCAPE=as.vector(mucape),LR75=as.vector(lr75),
#                 H5_TEMP=as.vector(h5_temp),S06=as.vector(s06),
#                 PMR=as.vector(pmr))
#Only using MUCAPE,S06,PMR
#era_data <- list(MUCAPE=as.vector(mucape),
#                 S06=as.vector(s06),
#                 PMR=as.vector(pmr))

#Only using LR75,S06,PMR
#era_data <- list(LR75=as.vector(lr75),
#                 S06=as.vector(s06),
#                 PMR=as.vector(pmr))

#Only using LR75,H5_TEMP,S06,PMR,CAP
era_data <- list(LR75=as.vector(lr75),
                 h5_TEMP=as.vector(h5_temp),
                 S06=as.vector(s06),
                 PMR=as.vector(pmr))



#Only using LR75,H5_TEMP,S06,PMR,CAP
#era_data <- list(LR75=as.vector(lr75),
#                 h5_TEMP=as.vector(h5_temp),
#                 S06=as.vector(s06),
#                 PMR=as.vector(pmr),
#                 CAP=as.vector(cap))


era_data <- as.data.frame(era_data)

#era_data_n <- as.data.frame(scale(era_data))
source("colScale.R")

#Keep as matrix for ccclust
era_data_n <- colScale(as.matrix(era_data))
#era_data_n <- scale(era_data)
#Use data frame for other prediction methods
#era_data_n <- as.data.frame(era_data_n)
names(era_data_n)

#era_data <- as.data.frame(as.vector(mucape),as.vector(lr75),
#                  as.vector(h5_temp),as.vector(s06),
#                  as.vector(pmr))



#Convert the kmeans model to kcca data type for processing
#cl.model <- as.kcca(km.model,dfn)
#cl.model <- as.kcca(km.model,dfn)
#pred_train <- predict(cl.model)
#pred_test <- predict(cl.model,newdata=era_data[1:184830,]) #1 month
#pred_test <- predict(cl.model,newdata=era_data_n[1:2248765,]) #1 year
#pred_test <- predict(cl.model,newdata=era_data_n[1:8995060,]) #1 year
#pred_test <- predict(km.model,newdata=era_data_n[1:8995060,]) #1 year
#pred_test <- predict(km.model,newdata=era_data_n) #6 years

#Using cclust object
#There is a nan in the matrix. Remove it!
#era_data_n <-era_data_n[!rowSums(!is.finite(era_data_n)),]
era_nan <- which(is.nan(era_data_n))
era_data_n[era_nan] <- mean(c(era_data_n[era_nan-1],era_data_n[era_nan+1]))
pred_test <- predict(km,newdata=era_data_n) #6 years

#We need to put the clusters in the same order as for the hail data
#pred_test.f <- factor(pred_test,levels=ind,ordered=TRUE,labels=as.character(1:ncls))
pred_test.f <- factor(pred_test$cluster,levels=ind,ordered=TRUE,labels=as.character(1:ncls))

#dim(pred_test) <- c(101,61,1460)
#dim(pred_test.f) <- c(101,61,1460)
#saveRDS(pred_test.f,"pred_test_2008_6clusters.rds")
#saveRDS(hd,"prop_hail_lightning_6clusters.rds")


#dim(pred_test) <- c(nlon,nlat,ntime)
dim(pred_test.f) <- c(nlon,nlat,ntime)


#fptout <- paste0("pred_test_all_lightning_kmeans_lloyd",ncls,"clusters.rds")
# *****
#fptout <- paste0("pred_test_all_lightning_cclust",ncls,"clusters.rds")
#fhpout <- paste0("prop_lightning_all_",ncls,"clusters.rds")
#fptout <- paste0("pred_test_all_lightning_cclust_seg_cape",ncls,"clusters.rds")
#fhpout <- paste0("prop_lightning_all_seg_cape",ncls,"clusters.rds")
#The following for only using subset of vars (MUCAPE,S06,PMR)
fptout <- paste0("pred_test_all_lightning_cclust_seg_cape_ex_vars",ncls,"clusters.rds")
fhpout <- paste0("prop_lightning_all_seg_cape_ex_vars",ncls,"clusters.rds")
saveRDS(pred_test.f,fptout)
saveRDS(hd,fhpout)


#len.obs <- 61*101 #One 1.0 observation
#tot.len <- length(era_data_n$MUCAPE_SMG)
#nyrs <- 7
#seq.len <- tot.len/7
##chnk <- seq(1,tot.len,by=seq.len)
#chnk <- seq(0,tot.len,length=14)
#pred_test <- vector("list",length=tot.len)
#for (i in 1:(length(chnk)-1)){
#    cat(i,"\n")
#    pred_test[[i]] <-
#        #predict(cl.model,newdata=era_data_n[chnk[i]:chnk[i+1]-1,])
#        predict(cl.model,newdata=era_data_n[chnk[i]+1:chnk[i+1],])
#}

#pred_test <- unlist(pred_test)
#dim(pred_test) <- c(101,61,10228)

