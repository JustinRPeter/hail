#Read in SHIP/GPATS/HAIL occurrence data
#an apply a logistic regression
#See the following for references
#http://www.r-bloggers.com/how-to-perform-a-logistic-regression-in-r/

rm(list=ls())
library(ggplot2)
library(dplyr)
library(reshape)
#library(leaps)
#library(pscl) #Used for zero-inflated regression
library(Zelig) #Used for rare events logistic regression
library(ROCR)
library(flexclust) #Used for predicting cluster membership of new data.
library(RColorBrewer)
library(wesanderson)
library(oce)
library(ncdf4)
library(fields)
library(oz)


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
#train.fraction <- 1.0
train.fraction <- 0.7
#train.fraction <- 0.5

#####################################CLUSTERING################################
#Scale the data before clustering or regression
#dfn <- as.data.frame(scale(data[,15:18]))
logreg_vars <- c("MUCAPE_SMG","LR75_SMG","H5_TEMP_SMG","S06_SMG","PMR_SMG")
#logreg_vars <- c("LR75_SMG","H5_TEMP_SMG","S06_SMG","PMR_SMG")
#logreg_vars <- c("MUCAPE_SMG","S06_SMG","PMR_SMG")
#logreg_vars <- c("LR75_SMG","S06_SMG","PMR_SMG")
dfn <- as.data.frame(scale(data[,logreg_vars]))
#dfn$HAIL <- data$HAIL
##names(dfn) <- c(names(data)[c(14:18,20)])
#names(dfn) <- c(logreg_vars,"HAIL")


#Perform a k-means clustering of the data.
#We will then conduct a logistic regression on each cluster
#ncls <- 2
#ncls <- 3
#ncls <- 4
#ncls <- 5
ncls <- 10
km.model <- kcca(dfn[,1:length(dfn)], k=ncls, family=kccaFamily("kmeans"))


#Add hail into the normalised data frame for logistic regression
dfn$HAIL <- data$HAIL
names(dfn) <- c(logreg_vars,"HAIL")


cls.data <- vector("list",length=ncls)
for (i in 1:ncls){
    ex.cl <- which(km.model@cluster == i)
    cls.data[[i]] <- dfn[ex.cl,]
}

#Print some information regarding how many of the clusters have
#a hail measurement.
#We only want to use these clusters to construct a regression model
#Furthermore this will get restricted to those clusters which have a hail
#event after we subset a training data set. See below (hl.train.clusters)
hl.ms <- c()
for (i in 1:ncls){
    hl.ms[i] <- length(which(cls.data[[i]]$HAIL == 1))
}
cat("Number of hail events in each cluster: ", hl.ms, "\n")
hl.clusters <- which(hl.ms > 0)
no.hl.clusters <- which(hl.ms == 0)
    
###############################LOGISTIC REGRESSION#############################
#Specify training and test data sets for each cluster
set.seed(200)
train <- vector("list",length=ncls)
test <- vector("list",length=ncls)
#for (i in 1:ncls){
for (i in hl.clusters){
    train[[i]] <- sample_frac(cls.data[[i]],train.fraction)
    test[[i]] <- subset(cls.data[[i]],!rownames(cls.data[[i]]) %in% rownames(train[[i]]))
}

prop.hail.events.cls <- c()
number.hail.events <- c()
for (i in 1:ncls){
    prop.hail.events.cls[i] <- length(which(train[[i]]$HAIL==1))/
                               length(which(train[[i]]$HAIL==0))
    number.hail.events[i] <- length(which(train[[i]]$HAIL==1))
}
#Which of the training data sets in each cluster have a hail observation?
#We can only conduct a regression if there is a hail observation.
#Because we extract train.fraction from each cluster, some of the clusters
#may not have a hail observation to enable the regression.
#hl.train.clusters <- which(prop.hail.events.cls >0)
hl.train.clusters <- which(number.hail.events >= 8)
cat("Clusters with enough hail obs. for regression: ",hl.train.clusters,"\n")

model.null   <- vector("list",length=ncls)
model.full   <- vector("list",length=ncls)
step.model   <- vector("list",length=ncls)
hail.model   <- vector("list",length=ncls)
fitted.model <- vector("list",length=ncls)
pr           <- vector("list",length=ncls)
prf          <- vector("list",length=ncls)
auc          <- vector("list",length=ncls)

#Perform a logistic regression on each cluster of the hail data
#Only perform it on those clusters which have a hail observed
##for (i in 1:ncls){
##for (i in hl.clusters){
#for (i in hl.train.clusters){
#    cat("Performing logistic regression on cluster", i, "\n")
#    model.null[[i]] <- glm(HAIL ~ 1, 
#                     data=train[[i]],
#                     family = binomial(link="logit"),
#                     control = list(maxit=50))
#
#    #model.full[[i]] <- glm(HAIL~MUCAPE_SMG+LR75_SMG+H5_TEMP_SMG+S06_SMG+PMR_SMG,
#    #model.full[[i]] <- glm(HAIL~LR75_SMG+H5_TEMP_SMG+S06_SMG+PMR_SMG,
#    rhs <- paste(logreg_vars,sep="",collapse="+")
#    lhs <- "HAIL ~ "
#    formula <- paste(lhs,rhs,sep="")
#    model.full[[i]] <- glm(formula,
#                      data=train[[i]],
#                      family=binomial(link="logit"),
#                      control = list(maxit=50))
#
#    #This will find the significant paramters in the model
#    #Using all the variables, I found that the most parsomonious model
#    #was one that regressed hail against MUCAPE_SMG, LR75_SMG ang S06_SMG
#    step.model[[i]] <- step(model.null[[i]],scope = list(upper=model.full[[i]]),
#         direction="both",
#         test="Chisq",
#         data=train[[i]])
#
#    hail.model[[i]] <- step.model[[i]]
#}
#

#Fit the models to the test data sets so we can get some error stats
##for (i in 1:ncls){
##for (i in hl.clusters){
#for (i in hl.train.clusters){
#    fitted.model[[i]] <- predict(hail.model[[i]],
#                               newdata=subset(test[[i]],
#                               #select=c(LR75_SMG,H5_TEMP_SMG,S06_SMG,PMR_SMG)),
#                               #select=logreg_vars),
#                               select=tail(names(step.model[[i]]$coefficients),-1)),
#                               type="response")
#}
#

#########################Using Rare Event Logistic Regression##################
#Instead of using normal logistic regression, we will use the Zelig package
#to use the relogit function.
#This is a rare events logistic model
#See the following pages:
#(1) http://docs.zeligproject.org/articles/zelig_relogit.html
#(2) http://zelig.gking.harvard.narkive.com/gm23ySXn/using-zelig-and-amelia-questions-about-odds-ratio-weights-step-method-r-square-adjusted-r-square

x.out <- vector("list",length=ncls)
s.out <- vector("list",length=ncls)
#for (i in hl.clusters){
for (i in hl.train.clusters){
    cat("Performing rare events logistic regression on cluster", i, "\n")

    model.null[[i]] <- glm(HAIL ~ 1, 
                       data=train[[i]],
                       family = binomial(link="logit"),
                       control = list(maxit=50))
    model.null[[i]] <- to_zelig(model.null[[i]])

    rhs <- paste(logreg_vars,sep="",collapse="+")
    lhs <- "HAIL ~ "
    formula <- paste(lhs,rhs,sep="")
    tau <- prop.hail.events.cls[i]
    #tau <- c(1.e-5,prop.hail.events.cls[[i]]+0.001)
    model.full[[i]] <- zelig(formula,
                             data=train[[i]],
                             model="relogit",
                             control=list(maxit=100),
                             tau=118/45553,
                             #tau=c(0.0002,0.003))
                             #tau=NULL,
                             #case.control=c("prior"),
                             case.control=c("weighting"),
                             robust=TRUE,
                             cite=FALSE)

     #This will find the significant paramters in the model
     #Using all the variables, I found that the most parsomonious model
     #was one that regressed hail against MUCAPE_SMG, LR75_SMG ang S06_SMG
     step.model[[i]] <- step(model.null[[i]],scope = list(upper=model.full[[i]]),
          direction="both",
          test="Chisq",
          data=train[[i]])
 
     hail.model[[i]] <- step.model[[i]]

#    hail.model[[i]] <- model.full[[i]]

     #x.out[[i]] <- setx(hail.model[[i]],test[[i]],cond=TRUE)
     x.out[[i]] <- setx(hail.model[[i]],test[[i]])
     s.out[[i]] <- sim(hail.model[[i]],x=x.out[[i]])
}

##for (i in hl.train.clusters){
#for (i in 1){
#    #fitted.model[[i]] <- predict(hail.model[[i]],
#    fitted.model[[i]] <- predict(hail.model[[i]],
#                               newdata=subset(test[[i]],
#                               #select=c(LR75_SMG,H5_TEMP_SMG,S06_SMG,PMR_SMG)),
#                               #select=logreg_vars),
#                               select=tail(names(coef(hail.model[[i]])),-1)),
#                               type="response")
#}

#Can't use predict method.
#See https://sebastiansauer.github.io/convert_logit2prob/
logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}

#Multiply each of the test data frame elements by the 
#hail.model coefficients
for (i in hl.train.clusters){
    pred1 <- test[[i]][,logreg_vars]*coef(hail.model[[i]])[logreg_vars]
    pred2 <- apply(pred1,1,sum)
    pred3 <- pred2+coef(hail.model[[i]])[c("(Intercept)")]
    fitted.model[[i]] <- logit2prob(pred3)
}

#Now construct ROC curves to test model
misClasificError <- c()
#for (i in 1:ncls){
for (i in hl.train.clusters){
    fitted.model[[i]] <- ifelse(fitted.model[[i]] > 0.5,1,0)
    misClasificError[i] <- mean(fitted.model[[i]] != test[[i]]$HAIL)
}


#for (i in hl.clusters){
for (i in hl.train.clusters){
        pr[[i]] <- prediction(fitted.model[[i]], test[[i]]$HAIL)
        prf[[i]] <- performance(pr[[i]], measure = "tpr", x.measure = "fpr")
}

#A variable to determine how many plots
pll <- length(hl.clusters)

#par(mfrow=c(3,2))
#for (i in hl.clusters){
for (i in hl.train.clusters){
        auc[[i]] <- performance(pr[[i]], measure = "auc")
        auc[[i]] <- auc[[i]]@y.values[[1]]
        auc[[i]]
}

misClasificError <- c()
for (i in 1:ncls){
    fitted.model[[i]] <- ifelse(fitted.model[[i]] > 0.5,1,0)
    misClasificError[i] <- mean(fitted.model[[i]] != test[[i]]$HAIL)
}


###############################PREDICTION ON ERA-Interim DATA##################
#Read the ERA data for 20008-2014 so that we can 
#compare the probability of lightning

#fnm_era <- '/home/u1066578/data/era_output/00_18/2008_2014/analysis/era_profiles_2008_2014_00_18.nc'
#fnm_era <- '/home/u1066578/data/era_output/00_18/1979_2014/analysis/era_profiles_1979_2014_00_18.nc'
fnm_era <- '/home/jpeter/Documents/usq/data/era_profiles_2008_2014_00_18.nc'
#fnm_era <- '/home/jpeter/Documents/usq/data/analysis/era_profiles_1979_2014_00_18.nc'

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
#source("/home/u1066578/jpeter/rfiles/lib/thermo/thermo.r")
source("/home/jpeter/Documents/usq/rfiles/lib/thermo/thermo.r")
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
#fnm_lsmask <- '/home/u1066578/data/eraint/lsmask/era_interim_lsmask.nc'
fnm_lsmask <- '/home/jpeter/Documents/usq/data/era_interim_lsmask.nc'
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
era_data <- list(MUCAPE_SMG=as.vector(mucape),LR75_SMG=as.vector(lr75),
                 H5_TEMP_SMG=as.vector(h5_temp),S06_SMG=as.vector(s06),
                 PMR_SMG=as.vector(pmr))

#era_data <- list(LR75_SMG=as.vector(lr75),
#                 H5_TEMP_SMG=as.vector(h5_temp),S06_SMG=as.vector(s06),
#                 PMR_SMG=as.vector(pmr))

era_data <- as.data.frame(era_data)
era_data <- era_data[,logreg_vars]



#era_data_n <- as.data.frame(scale(era_data))
source("colScale.R")
cat("Scaling ERA data before cluster membership prediction: ","\n")
era_data_n <- colScale(as.matrix(era_data))
#era_data_n <- scale(era_data)
era_data_n <- as.data.frame(era_data_n)
names(era_data_n) <- names(era_data)

#era_data <- as.data.frame(as.vector(mucape),as.vector(lr75),
#                  as.vector(h5_temp),as.vector(s06),
#                  as.vector(pmr))


##Clean up some memory
#rm(list=c("mucapel","lr75l","h5_templ","s06l","pmrl"))
#rm(list="era_data")
#gc()


#We have to remove the hail column from dfn before predicting
#on the ERA data set
dfn <- dfn[,logreg_vars]
cat("Predicting cluster membership of ERA data: ", "\n")
cls.pred <- predict(km.model,newdata=era_data_n) #All years

#We need to put the clusters in the same order as for the hail data
#cls.pred.f <- factor(cls.pred,levels=ind,ordered=TRUE,labels=as.character(1:ncls))

#dim(cls.pred) <- c(nlon,nlat,366*4)
#dim(cls.pred.f) <- c(nlon,nlat,366*4)
#dim(cls.pred.f) <- c(nlon,nlat,ntime)

#dim(cls.pred) <- c(nlon,nlat,ntime)

#Place the cluster into the ERA dataframe
era_data_n$CLS <- cls.pred
#Add a column for the hail probability
era_data_n$HAIL_PROB <- NA


#Plot the spatial distribution of each of the clusters

#Plot the cluster distribution across Australia
cat("Evaluating cluster distribution \n")
dim(cls.pred) <- c(nlon,nlat,ntime)
clus.count <- c()
clus.dist <- vector("list",length=ncls)
for (i in 1:ncls){
    clstemp <- cls.pred
    clus.count[i] <- length(which(clstemp==i))
    clstemp[clstemp != i] <- NA
    clstemp[clstemp == i] <- 1
    clus.dist[[i]] <- apply(clstemp,c(1,2),sum,na.rm=TRUE)/ntime
}

#Apply the land sea mask to the clus.dist
for (i in 1:ncls){
    clus.dist[[i]][lsmask==0] <- NA
}

#Evaluate the precentage of each cluster
clus.perc <- clus.count*100./sum(clus.count)


#ncols=50 #No of colors for plotting
ncols=15 #No of colors for plotting
#Evaluate the limits for the spatial plots
clus.dist.max <- c()
for (i in 1:ncls){
    clus.dist.max[i] <- max(clus.dist[[i]],na.rm=T)
}
cls.max <- round(max(clus.dist.max),digits=2)


pdf(paste0("Lightning_Cluster_spatial_distribution_",ncls,"_clusters.pdf"),
    width=10,height=10,paper="special")
#par(mfrow=c(3,ncls/3))
par(mfrow=c(2,2))#20 clusters
#par(mar=c(4,5,3,2,1)+0.1)
#par(mar=c(4,4,3,2)+0.1)
par(mar=c(4,2,1,1)+0.1)
for (i in 1:ncls){
    image.plot(lon,rev(lat),clus.dist[[i]],
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




#Now loop through the clusters in the ERA data and apply the
#logistic regression model to each of the clusters

cat("Predicting hail probability for each cluster: ", "\n")
#First the clusters that produce hail. Use the logistic regression.
#for (i in 1:ncls){
#for (i in hl.clusters){
for (i in hl.train.clusters){
#for (i in 1){
    ex.cl <- which(cls.pred == i)
    era.cls <- era_data_n[ex.cl,]
    #hail.prob.part <- predict(hail.model[[i]],newdata=era.cls,type="response")
    pred1 <- era.cls[,logreg_vars]*coef(hail.model[[i]])[logreg_vars]
    pred2 <- apply(pred1,1,sum)
    pred3 <- pred2+coef(hail.model[[i]])[c("(Intercept)")]
    hail.prob.part <- logit2prob(pred3)

    era_data_n[ex.cl,]$HAIL_PROB <- as.vector(hail.prob.part)
}

if (length(no.hl.clusters) >= 1){
    for (i in no.hl.clusters){
        ex.cl <- which(cls.pred == i)
        era.cls <- era_data_n[ex.cl,]
        #hail.prob.part <- 0
        era_data_n[ex.cl,]$HAIL_PROB <- 0
    }
}

hail.prob <- era_data_n$HAIL_PROB
dim(hail.prob) <- c(nlon,nlat,ntime)

hail_clim <- apply(hail.prob,c(1,2),sum,na.rm=TRUE)/(7*4)
image.plot(t(apply(hail_clim,1,rev)))

#
#fptout <- paste0("Logistic_regression_prediction_2008_2014.rds")
#saveRDS(cls.pred.f,fptout)
#saveRDS(hd,fhpout)




##Perform a cluster analysis of the SHIP components
#
##Scale the data
##dfn <- scale(data[,7:length(data)])
##dfn <- scale(data[,13:18])
##dfn <- scale(data[,14:19])
##dfn <- scale(data[,13:19])
##We only want to use CAPE, LR75, H5_TEMP, S06 and PMR in cluster calculations
##dfn <- as.data.frame(scale(data[,14:18]))
##names(dfn) <- c(names(data)[14:18])
##Using LR75,H5_TEMP,S06,PMR
#dfn <- as.data.frame(scale(data[,15:18]))
#names(dfn) <- c(names(data)[15:18])
##pairs(dfn)
#
##ncmax <- 15
#ncmax <- 10
#nrep <- 25
##nrep <- 50
#wss <- rep(0,ncmax-1)
#
##ncls=6
##ncls=10
##ncls=15
#ncls=20
##ncls=64
##ncls=100
##ncls=250
#cat("Calculating clustering with",ncls, "clusters...","\n")
##km.model <- kmeans(dfn[,1:6],centers=5,nstart=nrep,iter.max=250)
##km.model <- kmeans(dfn[,2:6],centers=5,nstart=nrep,iter.max=250)
##km.model <- kmeans(dfn[,1:5],centers=ncls,nstart=nrep,iter.max=250)
##km.model <- kmeans(dfn[,1:5],centers=ncls,nstart=nrep,iter.max=250,algorithm='Lloyd')
#km.model <-  kcca(dfn[,1:length(dfn)], k=ncls, family=kccaFamily("kmeans"))
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
#pdf(paste0("cluster_contribution_",ncls,"_clusters.pdf"),width=10,height=10,paper="special")
#barplot(cls.conts*100.,names.arg=c(1:ncls),
#            xlab="Cluster ordered by CAPE",
#            ylab="% contribution")
#dev.off()
#
##Use ggplot2 for plotting
#
##cluster <- factor(km.model$cluster)
##Order the clusters according to CAPE and name them
##1 through ncls for plotting convenience
##This also means we can reference them naturally
##cluster <- factor(km.model$cluster,
#cluster <- factor(km.model@cluster,
#                  levels=ind,ordered=TRUE,
#                  labels=as.character(1:ncls))
#
#df <- cbind(data[,13:20],cluster)
#dfg1 <- melt(df,id.vars=c("SHIP_Max_GPATS","Max_GPATS","cluster"),
#               measure.vars=c("MUCAPE_SMG","LR75_SMG","H5_TEMP_SMG","S06_SMG","PMR_SMG"))
#
#labels <- c(SHIP_Max_GPATS="SHIP",
#            Max_GPATS="Lightning counts",
#            MUCAPE_SMG="MUCAPE [J/kg]", 
#            LR75_SMG="LR75 [C/km]",
#            H5_TEMP_SMG= "T_H500 [C]",
#            S06_SMG="S06",
#            PMR_SMG="w [kg/kg]")
#
##Plot the SHIP as a function of the cluster variables
#p1<-ggplot(dfg1,aes(x=SHIP_Max_GPATS,y=value))+
#       xlab("SHIP")+ylab("")+
#       geom_point(aes(color=cluster))+
#       facet_wrap(~variable,ncol=2,scales="free",labeller=labeller(variable=labels[3:7]))+
#       scale_colour_manual(values=rev(oce.colorsViridis(length(ind))))+
#       guides(colour=guide_legend(override.aes=list(size=5)))
##plot(p1)
#
###Some plotting parameters
##xsize=ysize=15
###psize=48
##png("SHIP_vs_cls.png",width=xsize,height=ysize,units="cm",res=600)
##print(p1)
##dev.off()
#
##Plot the Lightning as a function of the cluster variables
#p2<-ggplot(dfg1,aes(x=Max_GPATS,y=value))+
#       xlab("Lightning counts / hour")+ylab("")+
#       geom_point(aes(color=cluster))+
#       facet_wrap(~variable,ncol=2,scales="free",labeller=labeller(variable=labels[3:7]))+
#       scale_colour_manual(values=rev(oce.colorsViridis(length(ind))))+
#       guides(colour=guide_legend(override.aes=list(size=5)))+
#       scale_x_log10()
##plot(p2)
##png("GPATS_vs_cls.png",xsize,height=ysize,units="cm",res=600)
##print(p2)
##dev.off()
#
#
##Now plot histograms of each of the variables by cluster
#
##Use the log(Max_GPATS) to differentiate
#log.gpats <- log10(data$Max_GPATS)
##temp.minus <- -1*data$H5_TEMP_SMG
#dflt <- cbind(data[,c(13,14,15,16,17,18)],log.gpats,cluster)
##names(dflt) <- names(data)[13:19]
#
#llabels <- c(SHIP_Max_GPATS="SHIP",
#            log.gpats="log(Lightning counts)",
#            MUCAPE_SMG="MUCAPE [J/kg]", 
#            LR75_SMG="LR75 [C/km]",
#            H5_TEMP_SMG= "T_H500 [C]",
#            S06_SMG="S06",
#            PMR_SMG="w [kg/kg]")
#
#
##dfg2 <- melt(df,id.vars=c("cluster"),
#dfg2 <- melt(dflt,id.vars=c("cluster"),
#            #measure.vars=c("SHIP_Max_GPATS",
#            #               "MUCAPE_SMG","LR75_SMG",
#            #               "H5_TEMP_SMG","S06_SMG","PMR_SMG"))
#            measure.vars=c("SHIP_Max_GPATS","log.gpats",
#                           "MUCAPE_SMG","LR75_SMG",
#                           "H5_TEMP_SMG","S06_SMG","PMR_SMG"))
##p3 <- ggplot(dfg,aes(x=SHIP_Max_GPATS,y=value))+
#p3 <- ggplot(dfg2,aes(x=value))+
#      xlab("")+
#      #stat_bin(aes(y=..density..,color=cluster),geom="step",size=1.5,position="identity")+
#      geom_freqpoly(aes(y=..density..,color=cluster),size=1.5)+
#      #stat_bin(aes(y=..count..,color=cluster),geom="step",size=1.5)+
#      facet_wrap(~variable,ncol=2,scale="free",labeller=labeller(variable=llabels))+
#      #scale_colour_manual(values=wes_palette(n=length(ind),name="Zissou"))+
#      scale_colour_manual(values=rev(oce.colorsViridis(length(ind))))+
#      guides(colour=guide_legend(override.aes=list(size=5)))
#
##plot(p3)
##png("hists_vars_by_cls.png",width=xsize,height=ysize,units="cm",res=600)
##print(p3)
##dev.off()
#
##Which clusters produce hail/nohail lightning/no lightning?
#hail <- which(df$HAIL ==1)
#no.hail <- which(df$HAIL == 0)
#length(hail)
#length(no.hail)
#tot.len.hail<-length(hail)+length(no.hail)
#no.days.hail <- length(unique(data$Date_GPATS[hail]))
#no.days.no.hail <- length(unique(data$Date_ERA[no.hail]))
#
#lightning <- which(df$Max_GPATS > 0)
#no.lightning <- which(df$Max_GPATS == 0)
#length(lightning)
#length(no.lightning)
#tot.len.lightning <- length(lightning)+length(no.lightning)
#
#prop.hail <- c()
#days.hail <- c()
#prop.no.hail <- c()
#days.no.hail <- c()
#prop.lightning <- c()
#prop.lightning.strikes <- c()
#days.lightning <- c()
#prop.no.lightning <- c()
#days.no.lightning <- c()
#
#tot.lightning.strikes <- sum(data$Max_GPATS[lightning])
#tot.hours <- tot.len.lightning
#lstrks.per.hour <- tot.lightning.strikes/tot.hours
#
##The number of years
#nyrs <- as.numeric(ceiling(difftime(max(data$Date_ERA),
#                         min(data$Date_ERA))/365.25))
#
#cat('\n Proportion of each cluster on which hail occurs \n')
#for (i in 1:ncls){
#    #cat('Cluster ',i,length(which(df$cluster[hail]==i))/length(hail),'\n')
#    #prop.hail[i] <- length(which(df$cluster[hail]==i))/length(hail)
#    #prop.lightning[i] <- length(which(df$cluster[lightning]==i))/
#    #                       length(lightning)
#    cat('Cluster ',i,"Hail",
#        length(which(df$cluster[hail]==i))/tot.len.hail,'\n')
#    cat('Cluster',i,"GPATS",
#        length(which(df$cluster[lightning]==i))/tot.len.lightning,'\n')
#    #prop.hail[i] <- length(which(df$cluster[hail]==i))/tot.len.hail
#    prop.hail[i] <- length(which(df$cluster[hail]==i))/
#                    length(which(df$cluster==i))
#    #days.hail[i] <- length(which(df$cluster[hail]==i))/nyrs
#    #days.hail[i] <- length(unique(data$Date_ERA[which(df$cluster[hail]==i)]))/
#                    #nyrs
#                    #length(unique(data
#    #prop.lightning[i] <- length(which(df$cluster[lightning]==i))/
#    #                       tot.len.lightning
#    prop.lightning[i] <- length(which(df$cluster[lightning]==i))/
#                         length(which(df$cluster==i))
#    prop.lightning.strikes[i] <- prop.lightning[i]*(sum(data$Max_GPATS[
#                                     lightning[which(df$cluster[lightning]==i)]])/
#                                 #tot.lightning.strikes
#                                 #tot.len.lightning
#                                 #length(lightning))
#                                 #length(which(df$cluster==i)))
#                                 length(which(df$cluster[lightning]==i)))
#                                 #(nyrs*365.25*32)
#    #prop.lightning.strikes[i] <- length(which(df$cluster[lightning] == i))
#    #days.lightning[i] <- length(which(df$cluster[lightning]==i))/nyrs
#}
##plot(df$cluster[hail])
##barplot(prop.hail*100.,names=as.character(1:ncls),
##        xlab="Cluster",ylab="%")
#
#
#cat('\n Proportion of each cluster on which NO hail/lightning occurs \n')
#for (i in 1:ncls){
#    #cat('Cluster ',i,length(which(df$cluster[no.hail]==i))/length(no.hail),'\n')
#    #prop.no.hail[i] <- length(which(df$cluster[no.hail]==i))/length(no.hail)
#    #prop.no.lightning[i] <- length(which(df$cluster[no.lightning]==i))/
#    #                        length(no.lightning)
#    cat('Cluster ',i,"No Hail",
#        length(which(df$cluster[no.hail]==i))/tot.len.hail,'\n')
#    cat('Cluster',i,"No GPATS",
#        length(which(df$cluster[no.lightning]==i))/tot.len.lightning,'\n')
#    
#    #prop.no.hail[i] <- length(which(df$cluster[no.hail]==i))/tot.len.hail
#    prop.no.hail[i] <- length(which(df$cluster[no.hail]==i))/
#                       length(which(df$cluster[no.hail]==i))
#    #days.no.hail[i] <- length(unique(data$Date_ERA[which(df$cluster[no.hail]==i)]))/nyrs
#    prop.no.lightning[i] <- length(which(df$cluster[no.lightning]==i))/
#                             length(which(df$cluster==i))
#                            #tot.len.lightning
#
#}
##plot(df$cluster[no.hail])
##barplot(prop.no.hail*100.,names=as.character(1:ncls),
##        xlab="Cluster",ylab="%")
#
##cl <- as.character(1:ncls)
#cl <- factor(1:ncls)
#hd <- data.frame(prop.hail,prop.no.hail,
#                 #days.hail,days.no.hail,
#                 prop.lightning,prop.no.lightning,
#                 prop.lightning.strikes,cl)
#hd.l <- melt(hd,id.vars="cl")
#phd <- ggplot(hd.l,aes(cl,value,fill=variable))+
#       geom_bar(stat="identity",position="dodge")+
#       xlab("Cluster")+ylab("%")+
#       scale_fill_brewer("",
#                         labels=c("Hail","No hail","Lightning","No lightning","Lightning strikes"),
#                        palette="Set1")
##plot(phd)
##png("hail_lightning_contribution_by_cls.png",width=xsize,height=ysize,units="cm",res=600)
##print(phd)
##dev.off()
#
#
##Now use the flexclust package to predict cluster association
##of data
#
#
