#Determine best number of clusters for cluster analysis.
#See http://stackoverflow.com/questions/15376075/cluster-analysis-in-r-determine-the-optimal-number-of-clusters
#for good tutorial on choosing appropriate number of clusters.

#In the following, I try three different formulations:
#(1) Cluster on CAPE
#(2) Cluster on log(CAPE)
#(3) Cluster by segmenting CAPE into equal size
#    partitions in an attempt to fully capture CAPE values.
#The lines labelled ***** are where I need to make changes
#to implement each one!
#Use /\*\*\*\*\* to search for string

rm(list=ls())
gc()
library(ggplot2)
library(dplyr)
#library(reshape)
#library(leaps)
#library(pscl)
#library(ROCR)
library(flexclust) #Used for predicting cluster membership of new data.
library(cclust) #flexclust aslo has a cclust method
#library(RColorBrewer)
#library(wesanderson)
#library(oce)
#library(ncdf4)
library(ff)
library(doParallel)
library(foreach)
library(gtools) #For ordering files on character and numeric 
library(Hmisc) #for cut2

source("colScale.R")
#source("/home/u1066578/jpeter/rfiles/era/colScale.R")


##raw.data <- read.table("lightning_predictor.txt",
#raw.data <- read.table("/home/u1066578/jpeter/rfiles/era/lightning_predictor.txt",
#                   header=TRUE,
#                   na.strings=c("", "NA"),
#                   stringsAsFactors=FALSE)

#raw.data <- readRDS("lightning_predictor.rds")
raw.data <- readRDS("lightning_predictor_0.25.rds")
#raw.data <- readRDS("/home/u1066578/jpeter/rfiles/era/lightning_predictor_0.25.rds")

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

##The fraction of data we will use for training
#train.fraction <- 0.7
#
##train.data <- data[sample(nrow(data),train.fraction*data.length),]
#train <- sample_frac(data,train.fraction)
#test <- subset(data,!rownames(data) %in% rownames(train))

#Use a smaller subset to test parallel computation
#data <- sample_frac(data,0.001)
#data <- sample_frac(data,0.01)
#data <- sample_frac(data,0.05)
#data <- sample_frac(data,0.1)
data <- sample_frac(data,1.0)

#  *****
#Use the log of CAPE
#I used the following for the hail paper
data <- data[data$MUCAPE > 0,]
#I used the following for the lightning paper
#data <- data[data$MUCAPE > 100,]
#The following to try log cape
#data$MUCAPE <- log(data$MUCAPE)

#Use variable to extract the variables
ex_vars <- c(7,8,9,10)

#dfm <- colScale(as.matrix(data[,6:10])) #used for clclust. needs matrix
dfm <- colScale(as.matrix(data[,ex_vars])) #used for clclust. needs matrix
##Get the scale and centering vectors
cnt <- attr(dfm,"scaled:center")
scl <- attr(dfm,"scaled:scale") 
cat(cnt,scl,"\n")
    ##temp <- sweep(clus.centres,MARGIN=2,cnt,"-")
    ##clus.centres.scaled <- sweep(temp,MARGIN=2,scl,"/")
dfn <- as.data.frame(dfm) #used for kmeans. needs data frame
#names(dfn) <- c(names(data)[6:10])
names(dfn) <- c(names(data)[ex_vars])



##Split the data set into n-clusters equally spaced
##along the range of MUCAPE values.
#ncls <- 20
#rng <- range(data$MUCAPE)
##brks <- c(0,exp(seq(log(1),log(rng[2]),length=ncls)))
#brks <- c(0,seq(100,rng[2],length=ncls))
#cuts <- cut(data$MUCAPE,breaks=brks,labels=FALSE,include.lowest=T)
#clus.centres <- matrix(nrow=ncls,ncol=5)
#
#for (i in 1:ncls){
#     cut.select <- which(cuts == i)
#     cat(length(cut.select),"\n")
#     samp <- sample_n(data[cut.select,6:10],1)
#     cat(as.character(samp),"\n")
#     clus.centres[i,]<-c(samp$MUCAPE,samp$LR75,samp$H5_TEMP,samp$S06,samp$PMR)
#}


#Perform a cluster analysis of the SHIP components

#Scale the data
#dfn <- scale(data[,7:length(data)])
#dfn <- scale(data[,13:18])
#dfn <- scale(data[,14:19])
#dfn <- scale(data[,13:19])
#We only want to use CAPE, LR75, H5_TEMP, S06 and PMR in cluster calculations
#dfn <- as.data.frame(scale(data[,6:10]))
#source("colScale.R")
#dfn <- colScale(as.matrix(data[,6:10]))
#names(dfn) <- c(names(data)[6:10])
##Get the scale and centering vectors
#cnt <- attr(dfn,"scaled:center")
#scl <- attr(dfn,"scaled:scale") 
#temp <- sweep(clus.centres,MARGIN=2,cnt,"-")
#clus.centres.scaled <- sweep(temp,MARGIN=2,scl,"/")
#dfn <- as.data.frame(dfn)
##pairs(dfn)

#-----(1) Using the elbow method in the sum of squares
#----- Only use this block of code to find the best number of clusters
#ncmax <- 100
##ncmax <- 20
#ncmax <- 50
##nrep <- 25
#nrep <- 1000
##nrep <- 50
##wss <- rep(0,ncmax-1)
#wss <- rep(NA,ncmax-1)
#wss <- rep(NA,seq(2,ncmax,by=2))
#
###sink('ship.cluster.txt',type='output',split=T)
##
##cat('\nNo.\tWithinSS\tMinSize\n')
##
###for (nc in 2:ncmax) {
###for (nc in seq(2,ncmax,by=10)){
##for (nc in seq(2,ncmax,by=2)){
##for (nc in 30){
#for (nc in 50){

#kmeans specification
nrep <- 25
#nrep <- 1
niter <- 500

#parallel core specs.
ncls <- c(2,3,4,5,6,10,15,20,50,100)
no_cores <- length(ncls)
#cl <- makeCluster(no_cores)
registerDoParallel(cores=no_cores)
#registerDoParallel(cl)
#source("colScale.R")



#sink("within_sum_squares_kmeans_lloyd.txt",append=T)
wss <- rep(NA,max(ncls))
#wss <- vector("list",length=no_cores)
#cluster.model <- vector("list",length=max(ncls))
#cls.model <- vector("list",length=max(ncls))
#for (nc in ncls){
#foreach (nc = c(6,10,15,20,50)) %dopar% {
#cls.model <- foreach (nc = c(6,10,15,20,50),.combine='c') %dopar% {
#cls.model <- foreach (nc = ncls,.combine='c') %dopar% {
#cls.model <- foreach(nc=ncls,.combine=c,.export="cclust",.packages="cclust") %dopar% {
#cls.model <- foreach(nc=ncls,.combine=c) %dopar% {

#If we want to return the results of the foreach function
#we need to place the reuslts in a list
#We will output the results of the cluster analysis individually
#and just return the within sum squares for plotting
comb <- function(...) {
            mapply(cbind, ..., SIMPLIFY=FALSE)
        }

#cls.model <- foreach(nc=ncls,.combine=comb,.multicombine=T,
#                     .init=list(list(), list())) %dopar% {
#cls.model <- foreach(nc=ncls,.combine=c,.multicombine=T,
#                     .init=list(list(),list())) %dopar% {
cls.model <- foreach(nc=ncls,.combine=c) %dopar% {
  
    # *****
    #Only use next 6 lines for option (3)
    #Split the data set into equal size chunks based on CAPE
    #dfn.split <- split(dfn, cut2(dfn$MUCAPE,g=nc))
    #rand.centre <- vector("list",length=dim(dfn)[2])
    #for (i in 1:nc){
    #    rand.centre[[i]] <- dfn.split[[i]][sample(nrow(dfn.split[[i]]),1),]
    #}
    #clus.centres <- matrix(unlist(rand.centre),nrow=nc,byrow=T)
        

    #sink(paste0("within_sum_squares_kmeans_lloyd_",nc,".txt"),append=T)
    #sink(paste0("within_sum_squares_kmeans_lloyd_",nc,".txt"))

    ###Split the data set into n-clusters equally spaced
    ###along the range of MUCAPE values.
    #rng <- range(data$MUCAPE)
    #brks <- c(0,seq(100,rng[2],length=ncls))
    #cuts <- cut(data$MUCAPE,breaks=brks,labels=FALSE,include.lowest=T)
    #clus.centres <- matrix(nrow=nc,ncol=5)
    ##
    #for (i in 1:nc){
    #     cut.select <- which(cuts == i)
    #     cat(length(cut.select),"\n")
    #     samp <- sample_n(data[cut.select,6:10],1)
    #     cat(as.character(samp),"\n")
    #     clus.centres[i,]<-c(samp$MUCAPE,samp$LR75,samp$H5_TEMP,samp$S06,samp$PMR)
    #}
#
    #dfn <- colScale(as.matrix(data[,6:10]))
    #names(dfn) <- c(names(data)[6:10])
    ##Get the scale and centering vectors
    #cnt <- attr(dfn,"scaled:center")
    #scl <- attr(dfn,"scaled:scale") 
    #print(cnt,scl)
    ##temp <- sweep(clus.centres,MARGIN=2,cnt,"-")
    ##clus.centres.scaled <- sweep(temp,MARGIN=2,scl,"/")
    #dfn <- as.data.frame(dfn)

    #Perform the clustering
    #km <- kmeans(dfn[,1:5],centers=clus.centres.scaled,iter.max=500,algorithm="MacQueen")
    #wss <- rep(NA,nc)
    #km <- kmeans(dfn,centers=nc,nstart=nrep,iter.max=niter,algorithm="Lloyd")
    #km <- cclust(dfm,nc,niter,verbose=T,method="kmeans")
    
    # *****
    #The next two lines are the ones to use.
    km <- cclust::cclust(dfm,nc,niter,verbose=TRUE,method="kmeans")
    #km <- cclust::cclust(dfm,clus.centres,niter,verbose=TRUE,method="kmeans")

    #km <- stats::kmeans(dfm,nc,niter,nrep)
    #cluster.model[[nc]] <- km
    #wss[nc] <- sum(km$withinss)
    #kmeans(dfn,centers=nc,nstart=nrep,iter.max=niter,algorithm="Lloyd")
    #cluster.model[[nc]] <- km
    wss[nc]<-sum(km$withinss)
    #wss[[nc]]<-sum(km$withinss)

    #cat("ncls=",nc,"\n",
    #    "within sum squares=", wss[nc],"\n")

    # ***** 
    #Save to the appropriate output files
    #saveRDS(km,paste0("km.model.","kmeans_lloyd_",nc,"_clusters.rds"))
    #saveRDS(km,paste0("km.model.","kmeans_cclust_",nc,"_clusters.rds"))
    #saveRDS(km,paste0("km.model.","kmeans_cclust_log_cape",nc,"_clusters.rds"))
    #saveRDS(km,paste0("km.model.","kmeans_cclust_seg_cape",nc,"_clusters.rds"))
    saveRDS(km,paste0("km.model.","kmeans_cclust_seg_cape_ex_vars",nc,"_clusters.rds"))
    return(wss)
    #return(km)
    #Use next line to output noth in a list but don't really need to do that
    #as we output the results of clustering to rds files
    #list(km,wss)
    #return(res)
    #sink()
}
#sink()
#stopCluster(cl)

# *****
#Get the name right here!
#Output the within sum squares to a file
#sink('Within_sum_squares_output.txt',append=F)
#sink('Within_sum_squares_output_log_cape.txt',append=F)
#sink('Within_sum_squares_output_seg_cape.txt',append=F)
sink('Within_sum_squares_output_seg_cape_ex_vars.txt',append=F)
cls.model.wss <- cls.model[!is.na(cls.model)]
cat('Number of clusters','\n')
cat(ncls,sep='\t','\n','\n')
cat('Within sum squares','\n')
cat(cls.model.wss,sep='\t','\n')
sink()

wss.dat <- data.frame(ncls,cls.model.wss)
ggplot(wss.dat,aes(x=ncls,y=cls.model.wss))+
    geom_point(shape=1)

#Read all the cluster model rds files
#and plot the total amount of variance explained
#This is betwennss/totss
#Need to evaluate the total sum of squares from cclust
dfmsq <- apply(dfm,2,function(x) x^2.)
sumsqdfm <- apply(dfmsq,2,sum)
sumdfm <- apply(dfm,2,sum)
nsamps <- length(dfm[,1])
totss <- sum(sumsqdfm - (1./nsamps)*sumdfm^2.)

# *****
# Get the names right here
#Read all the files
#fls <- list.files(pattern="*kmeans_cclust*")
#fls <- list.files(pattern="*kmeans_cclust_log_cape*")
#fls <- list.files(pattern="*kmeans_cclust_seg_cape*")
fls <- list.files(pattern="*kmeans_cclust_seg_cape_ex_vars*")

fls <- mixedsort(fls)
ndata <- lapply(fls,function(x) readRDS(x))
#The amount of variance is given by
#(totss-tot.withinss=betweenss)/totss
varexp <- vector("numeric",length=no_cores)
for (i in 1:no_cores){
    varexp[i] <- (totss-sum(ndata[[i]]$withinss))*100./totss
}
#totss-sum(ndata[[1]]$withinss)

# *****
# Get the names right here
#Plot the within sum squares and variance explained
#pdf("withinss.pdf",width=10,height=10)
#pdf("withinss_log_cape.pdf",width=10,height=10)
#pdf("withinss_seg_cape.pdf",width=10,height=10)
pdf("withinss_seg_cape_ex_vars.pdf",width=10,height=10)
par(mar = c(5,5,2,5),cex=1.25)
plot(ncls,cls.model.wss,xlim=c(0,20),type="b",lwd=3,
     xlab="No. of clusters",
     ylab="Within cluster sum squares distance",cex=1.25)
par(new=T,cex=1.25)
plot(ncls,varexp,axes=F,xlab=NA,ylab=NA,xlim=c(0,20),
     cex=1.25,type="b",col="blue",lwd=3)
axis(side=4)
mtext(side = 4, line = 3, '% of variance',cex=1.25)
legend("topright",
       legend=c("Within cluster sum of squares","% variance explained"),
       col=c("black","blue"),
       lwd=3)
dev.off()






#To use the NbClust package
#gloop <- sample_n(dfn(,10000)
#res <- NbClust(gloop,distance="euclidean",min.nc=2,max.nc=7,method="complete",index="ch")
#res$All.index
#res$Best.nc
#res$All.CriticalValues
#res$Best.partition


##    #km <- kmeans(dfn,centers=nc,nstart=nrep,iter.max=100)
##    #km <- kmeans(dfn,centers=nc,nstart=nrep,iter.max=250,algorithm="Lloyd")
##     km <- kmeans(dfn[,1:5],centers=nc,nstart=nrep,iter.max=250,algorithm="Lloyd")
#
#     #km.model <-  kcca(dfn[,1:5], k=ncls, family=kccaFamily("kmeans"))
##    #Don't use SHIP in cluster calculations
##    #km <-
##    kmeans(dfn[,2:6],centers=nc,nstart=nrep,iter.max=250,algorithm="Lloyd")
##    #km <-
##    kmeans(dfn[,1:5],centers=nc,nstart=nrep,iter.max=250,algorithm="Lloyd")
##    #km <- kmeans(dfn[,1:5],centers=nc,nstart=nrep,iter.max=500)


#    km.model <- kcca(dfn[,1:5],k=nc,family=kccaFamily("kmeans"))
##    #wss[nc] <- sum(km$withinss)
#    wss[nc] <- info(km.model,"distsum")
#    #Write the cluster output to a file
#    saveRDS(km.model,paste0("km.model.","kcca_",nc,"_clusters.rds"))
#    cat(nc,wss[nc],min(info(km.model,"size")),'\n',sep='\t')
##    #cat(nc,wss[nc],min(km$size),'\n',sep='\t')

#    cat(nc,wss[nc],min(info(km.model,"size")),'\n',sep='\t')
##    plot(km.model)
##    pairs(km.model)
##    image(km.model)
#}
##
#sink()
#Output the sum of squares to a file
#sink("witin_sum_squares_kcca.txt",append=T)
#cat(nc,wss[nc],min(info(km.model,"size")),'\n',sep='\t')
#sink()

#Plot the sum of squares
#pdf('cluster.wss.pdf',width=10,height=10,paper="special")
#plot(wss,main='Cluster SS',xlab='No. Clusters',ylab='Sum of Squares')
#dev.off()

##-----(1) Using the elbow method in the sum of squares
#library(fpc)
#d <- sample_n(dfn,1.e04)
#pamk.best <- pamk(d,krange=seq(2,50,by=4)
#sink("lightning_best_no_clusters_pamk.txt",type='output',append=FALSE)
#cat("number of clusters estimated by optimum average silhouette width:", pamk.best$nc, "\n")
#sink()
#pdf("pamk_clusters.pdf",width=10,height=10,paper="special")
#plot(pam(d, pamk.best$nc))
#dev.off()



