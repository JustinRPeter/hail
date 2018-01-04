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
library(RColorBrewer)
library(wesanderson)
library(oce)


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


#pdf("cdplots.pdf",width=5,height=5)
#   cdplot(as.factor(data$HAIL) ~ data$MESH)
#   abline(v=mean(data$MESH))
#   abline(v=mean(data$MESH)+sd(data$MESH),lty=2)
#   abline(v=mean(data$MESH)-sd(data$MESH),lty=2)
#   
#   cdplot(as.factor(data$HAIL) ~ data$VIL)
#   abline(v=mean(data$VIL))
#   abline(v=mean(data$VIL)+sd(data$VIL),lty=2)
#   abline(v=mean(data$VIL)-sd(data$VIL),lty=2)
#dev.off()

#Evaluate the length of the data set
data.length <- length(data$Time_GPATS)

#The fraction of data we will use for training
train.fraction <- 0.7

#train.data <- data[sample(nrow(data),train.fraction*data.length),]
train <- sample_frac(data,train.fraction)
test <- subset(data,!rownames(data) %in% rownames(train))


#A glm using all variables
glm1 <- glm(HAIL ~ MUCAPE_SMG+LR75_SMG+H5_TEMP_SMG+S06_SMG+PMR_SMG,
            family=binomial(),data=train)

#Use leaps to use less model parameters
leaps=regsubsets(HAIL ~ MUCAPE_SMG+LR75_SMG+H5_TEMP_SMG+S06_SMG+PMR_SMG,data=train)
plot(leaps, scale="adjr2")
plot(leaps, scale="bic")

#Use step function
null <- glm(HAIL ~1, family=binomial(),data=train)
full <- glm1
forward <- step(null, scope=list(lower=null, upper=full), direction="forward")
backward <- step(full, data=train,direction="backward")
#After stepwise analysis
#We get different results for the best glm based on whether we include 
#observations with no lightning
bestglm.all <- glm(formula = HAIL ~ MUCAPE_SMG+S06_SMG+LR75_SMG+H5_TEMP_SMG,
                   family=binomial(),data=train)
bestglm.lightning <- glm(formula = HAIL~MUCAPE_SMG+LR75_SMG+S06_SMG+PMR_SMG,
                         family = binomial(), data=train)


#Perform a cluster analysis of the SHIP components

#Scale the data
#dfn <- scale(data[,7:length(data)])
#dfn <- scale(data[,13:18])
#dfn <- scale(data[,14:19])
#dfn <- scale(data[,13:19])
#We only want to use CAPE, LR75, H5_TEMP, S06 and PMR in cluster calculations
dfn <- as.data.frame(scale(data[,14:18]))
names(dfn) <- c(names(data)[14:18])
#pairs(dfn)

ncmax <- 20
#ncmax <- 10
#nrep <- 25
nrep <- 10
#nrep <- 50
wss <- rep(0,ncmax-1)

#sink('ship.cluster.txt',type='output',split=T)

cat('\nNo.\tWithinSS\tMinSize\n')

for (nc in 2:ncmax) {
    #km <- kmeans(dfn,centers=nc,nstart=nrep,iter.max=100)
    #km <- kmeans(dfn,centers=nc,nstart=nrep,iter.max=250,algorithm="Lloyd")
    #km <- kmeans(dfn[,1:6],centers=nc,nstart=nrep,iter.max=250,algorithm="Lloyd")
    #Don't use SHIP in cluster calculations
    #km <- kmeans(dfn[,2:6],centers=nc,nstart=nrep,iter.max=250,algorithm="Lloyd")
    #km <- kmeans(dfn[,1:5],centers=nc,nstart=nrep,iter.max=250,algorithm="Lloyd")
    km <- kmeans(dfn[,1:5],centers=nc,nstart=nrep,iter.max=500)
    wss[nc] <- sum(km$withinss)
    cat(nc,wss[nc],min(km$size),'\n',sep='\t')
}

png('sng.cluster.wss.png')
#postscript('sng.cluster.wss.ps')
plot(wss,main='Cluster SS',xlab='No. Clusters',ylab='Sum of Squares')
dev.off()

ncls=10
#km.model <- kmeans(dfn[,1:6],centers=5,nstart=nrep,iter.max=250)
#km.model <- kmeans(dfn[,2:6],centers=5,nstart=nrep,iter.max=250)
km.model <- kmeans(dfn[,1:5],centers=ncls,nstart=nrep,iter.max=250)
#ncls <- length(km.model$size)
cat('\n',paste0(ncls),'Clusters\n')
cat('Cluster\tSize\tCentres\n')
#ind gives the order of the clusters based on the ordering of CAPE
ind <- order(km.model$centers[,1])
for (i in 1:ncls) cat(i,km.model$size[ind[i]],km.model$centers[ind[i],],'\n',sep='\t')
#for (i in 1:3) cat(i,km3$size[i],km3$centers[i,],'\n',sep='\t')
png('sng.cluster.var3.png')
#postscript('sng.cluster.var3.ps',paper="special",height=10,width=12,horizontal=F)
par(mfrow=c(3,2))

#The following lines are to plot graphs quickly (no units)
#for (i in 14:18)
#for (i in c(14,15,16,17,18,20))
for (i in c(14,15,16,17,18))
    #plot(data[,13],data[,i],col=km.model$cluster,xlab='SHIP',ylab=names(data[i]),cex.lab=1.25)
    plot(data[,13],data[,i],col=km.model$cluster,xlab='SHIP',ylab=names(data[i]),cex.lab=1.25)

#sink()

#Use ggplot2 for plotting

#cluster <- factor(km.model$cluster)
#Order the clusters according to CAPE and name them
#1 through ncls for plotting convenience
#This also means we can reference them naturally
cluster <- factor(km.model$cluster,
                  levels=ind,ordered=TRUE,
                  #labels=c("1","2","3","4","5"))
                  labels=as.character(1:ncls))
#cluster <- factor(order(km.model$cluster)[ind])
#cluster <- km.model$cluster[order(match(km.model$cluster,ind))]
#cluster <- cluster[ind]
#data <- data[ind,]
df <- cbind(data[,13:20],cluster)
dfg1 <- melt(df,id.vars=c("SHIP_Max_GPATS","Max_GPATS","cluster"),
               measure.vars=c("MUCAPE_SMG","LR75_SMG","H5_TEMP_SMG","S06_SMG","PMR_SMG"))

#levels(dfg1$variable)  <- c(MUCAPE_SMG="MUCAPE [J/kg]", 
            #LR75_SMG=expression(paste0("LR75 [",~degree,"C/km]")),
            #LR75_SMG="LR75 [ *degree C/km]",
labels <- c(SHIP_Max_GPATS="SHIP",
            Max_GPATS="Lightning counts",
            MUCAPE_SMG="MUCAPE [J/kg]", 
            LR75_SMG="LR75 [C/km]",
            H5_TEMP_SMG= "T_H500 [C]",
            S06_SMG="S06",
            PMR_SMG="w [kg/kg]")

#bloop <- expression(paste0("LR75 [",~degree,"C/km]"))
#levels(dfg1$variable) <- list(MUCAPE_SMG="MUCAPE [J/kg]",
#                              "LR75 [ *degreee,
#                              "T_H500 [C]",
#                              "S06",
#                              "w [kg/kg]")

#Plot the SHIP as a function of the cluster variables
p1<-ggplot(dfg1,aes(x=SHIP_Max_GPATS,y=value))+
       xlab("SHIP")+ylab("")+
       geom_point(aes(color=cluster))+
       facet_wrap(~variable,ncol=2,scales="free",labeller=labeller(variable=labels[3:7]))+
       #facet_wrap(~variable,scales="free",labeller=label_parsed)+
       #facet_wrap(~variable,scales="free",labeller=as_labeller(varnames))+
       #facet_grid(. ~ variable,scales="free")+
       #facet_grid(cluster~variable,scales="free")+
       #facet_wrap(~variable,scales="free",switch="y")+
       #theme(strip.background = element_blank())+
       #scale_colour_brewer(palette="Set1")+
       scale_colour_manual(values=rev(oce.colorsViridis(length(ind))))+
       guides(colour=guide_legend(override.aes=list(size=5)))
#plot(p1)

#Some plotting parameters
xsize=ysize=15
#psize=48
png("SHIP_vs_cls.png",width=xsize,height=ysize,units="cm",res=600)
print(p1)
dev.off()

#p1 <- ggplot(data,aes(x=SHIP_Max_GPATS,y=MUCAPE_SMG))+
#          geom_point(aes(color=cluster))+
#          scale_colour_manual(values=rev(oce.colorsViridis(length(ind))))

#Plot the Lightning as a function of the cluster variables
p2<-ggplot(dfg1,aes(x=Max_GPATS,y=value))+
       xlab("Lightning counts / hour")+ylab("")+
       geom_point(aes(color=cluster))+
       facet_wrap(~variable,ncol=2,scales="free",labeller=labeller(variable=labels[3:7]))+
       scale_colour_manual(values=rev(oce.colorsViridis(length(ind))))+
       guides(colour=guide_legend(override.aes=list(size=5)))+
       scale_x_log10()
#plot(p2)
png("GPATS_vs_cls.png",xsize,height=ysize,units="cm",res=600)
print(p2)
dev.off()


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
      #stat_bin(aes(y=..density..,color=cluster),geom="step",size=1.5)+
      #stat_bin(aes(y=..density..,color=cluster),geom="step",size=1.5,position="identity")+
      #geom_bar(aes(color=cluster),size=1.5,position="dodge")+
      geom_freqpoly(aes(y=..density..,color=cluster),size=1.5)+
      #geom_freqpoly(aes(y=..density..,color=cluster),stat="bin",size=1.5)+
      #stat_bin(aes(y=..count..,color=cluster),geom="step",size=1.5)+
      #stat_bin(aes(y=..ndensity..,color=cluster),geom="step",position="jitter",size=1.5)+
      #geom_histogram(aes(y=..density..,color=cluster),stat="identity",position="stack",size=1.5)+
      facet_wrap(~variable,ncol=2,scale="free",labeller=labeller(variable=llabels))+
      #facet_wrap(~variable)+
      #scale_colour_brewer(palette="Set1")+
      #scale_colour_brewer(palette="PuRd")+
      #scale_colour_grey(start=0.5,end=0)+
      #scale_colour_manual(values=topo.colors(5))+
      #scale_colour_manual(values=wes_palette(n=length(ind),name="Zissou"))+
      scale_colour_manual(values=rev(oce.colorsViridis(length(ind))))+
      guides(colour=guide_legend(override.aes=list(size=5)))
      #scale_x_log10()
      #scale_y_log10()
      #coord_trans(y="log10")
#plot(p3)
png("hists_vars_by_cls.png",width=xsize,height=ysize,units="cm",res=600)
print(p3)
dev.off()

#Which clusters produce hail/nohail lightning/no lightning?
hail <- which(df$HAIL ==1)
no.hail <- which(df$HAIL == 0)
length(hail)
length(no.hail)

lightning <- which(df$Max_GPATS > 0)
no.lightning <- which(df$Max_GPATS == 0)
length(lightning)
length(no.lightning)

prop.hail <- c()
prop.no.hail <- c()
prop.lightning <- c()
prop.no.lightning <- c()

cat('\n Proportion of each cluster on which hail occurs \n')
for (i in 1:ncls){
    cat('Cluster ',i,length(which(df$cluster[hail]==i))/length(hail),'\n')
    prop.hail[i] <- length(which(df$cluster[hail]==i))/length(hail)
    prop.lightning[i] <- length(which(df$cluster[lightning]==i))/
                           length(lightning)
}
#plot(df$cluster[hail])
#barplot(prop.hail*100.,names=as.character(1:ncls),
#        xlab="Cluster",ylab="%")


cat('\n Proportion of each cluster on which NO hail/lightning occurs \n')
for (i in 1:ncls){
    cat('Cluster ',i,length(which(df$cluster[no.hail]==i))/length(no.hail),'\n')
    prop.no.hail[i] <- length(which(df$cluster[no.hail]==i))/length(no.hail)
    prop.no.lightning[i] <- length(which(df$cluster[no.lightning]==i))/
                            length(no.lightning)
}
#plot(df$cluster[no.hail])
#barplot(prop.no.hail*100.,names=as.character(1:ncls),
#        xlab="Cluster",ylab="%")

#cl <- as.character(1:ncls)
cl <- factor(1:ncls)
hd <- data.frame(prop.hail,prop.no.hail,prop.lightning,prop.no.lightning,cl)
hd.l <- melt(hd,id.vars="cl")
phd <- ggplot(hd.l,aes(cl,value,fill=variable))+
       geom_bar(stat="identity",position="dodge")+
       xlab("Cluster")+ylab("%")+
       scale_fill_brewer("",
                         labels=c("Hail","No hail","Lightning","No lightning"),
                        palette="Set1")
#plot(phd)
png("hail_lightning_contribution_by_cls.png",width=xsize,height=ysize,units="cm",res=600)
print(phd)
dev.off()



###Develop the logistic model on the training data set
###GLM 1 - Hail Damage
##gpats_glm <- glm(Hail_Occurred ~ Max_GPATS,family=binomial())
###This is probability of hail vs max. ship on the day
##ship_glm_1  <- glm(data$Hail_Occurred ~ data$Max_SHIP,family=binomial())
###This is probability of hail vs. ship at max. gpats
##ship_glm_2  <- glm(data$Hail_Occurred ~ data$SHIP_Max_GPATS,faily=binomial())
#
## Now test predictive ability on the test data set
#glm.ship       <- glm(HAIL ~ Max_SHIP,family=binomial(),data=train)
#glm.shipmg     <- glm(HAIL ~ SHIP_Max_GPATS,family=binomial(),data=train)
#glm.gpats      <- glm(HAIL ~ Max_GPATS, family=binomial(),data=train)
##glm.log.gpats  <- glm(HAIL ~ log(Max_GPATS), family=binomial(),data=train)
##glm.shipmg.log.gpats <- glm(HAIL ~ log(Max_GPATS)+SHIP_Max_GPATS,family=binomial(),data=train)
#glm.shipmg.gpats <- glm(HAIL ~ SHIP_Max_GPATS+Max_GPATS,family=binomial(),data=train)
##glm.shipmg.gpats <- glm(HAIL ~ SHIP_Max_GPATS*Max_GPATS,family=binomial(),data=train)
#
##params <- c("Max_SHIP","SHIP_Max_GPATS","Max_GPATS",
##            #"log(Max_GPATS)","log(Max_GPATS)/SHIP_Max_GPATS")
##            "SHIP_MG_GPATS")
#
#
##Put all of the glms into a list for loop processing
#glms <- list(glm.ship,
#             glm.shipmg,
#             glm.gpats,
#             #glm.log.gpats,
#             #glm.ship.log.gpats,
#             glm.shipmg.gpats)
#names(glms) <- c("SHIP","SHIP_MG","GPATS","SHIP_MG_GPATS")
#
#for (i in 1:length(glms)){
#    summary(glms[[i]])
#    anova(glms[[i]],test="Chisq")
#    confint(glms[[i]])
#    print(exp(coef(glms[[i]])))
#    print(exp(confint(glms[[i]])))
#}
#
##model_vil <- glm(train$HAIL ~ train$VIL,family=binomial(link=logit))
##See http://www.theanalysisfactor.com/r-tutorial-glm1/ for why you have to use
##the lower model description. I have to look into this!
##model_vil <- glm(HAIL ~ VIL,family=binomial, data=test)
#
##Compute some pseudo-R2 measures for the GLM fits
#pr2 <- vector("list",length=length(glms))
#for (i in 1:length(glms)){
#    pr2[[i]]<-pR2(glms[[i]])
#}
#
##Now fit the results to the test case
##Which we nee to evaluate the accuracy of the model
#
#fit.ship.test <- predict(glm.ship,
#                       newdata=subset(test,select=c(7,10)),
#                       type="response")
#
#fit.shipmg.test <- predict(glm.shipmg,
#                              newdata=subset(test,select=c(8,10)),
#                              type="response")
#
#fit.gpats.test <- predict(glm.gpats,
#                      newdata=subset(test,select=c(9,10)),
#                      type="response")
#
#fit.shipmg.gpats.test <- predict(glm.shipmg.gpats,
#                       newdata=subset(test,select=c(8,9,10)),
#                       type="response")
#
##Evaluate the accuracy of each model in determining haail presence
##Requires turning glm model into a binary response
##(>0.5 --> hail and <= 0.5 --> no hail)
#fits.test<-list(fit.ship.test,fit.shipmg.test,fit.gpats.test,fit.shipmg.gpats.test)
#names(fits.test) <- c("SHIP","SHIP_MG","GPATS","SHIP_MG_GPATS")
#fits.test.bi <- vector("list",length=length(fits.test))
#misClasificError      <- vector("numeric",length=length(fits.test))
#accuracy              <- vector("numeric",length=length(fits.test))
#for (i in seq_along(fits.test)){
#    fits.test.bi[[i]] <- ifelse(fits.test[[i]] > 0.5,1,0)
#    misClasificError[i] <- mean(fits.test.bi[[i]] != test$HAIL)
#    #misClasificError <- mean(fitted.results.vil.bi != test$HAIL)
#    #print(paste(names(test)[i],'Accuracy',1-misClasificError[[i]]))
#    accuracy[i] <- 1.-misClasificError[i]
#}
#
##Examine the ROCR curves for the models
##See https://www.r-bloggers.com/how-to-perform-a-logistic-regression-in-r/
#pr  <- vector("list",length=length(fits.test))
#prf <- vector("list",length=length(fits.test))
#for (i in seq_along(fits.test)){
#    pr[[i]]  <- prediction(fits.test[[i]],test$HAIL)
#    prf[[i]] <- performance(pr[[i]],measure="tpr",x.measure="fpr")
#}
##Some performance measures
#auc <- vector("list",length=length(fits.test))
#for (i in seq_along(fits.test)){
#    auc[[i]] <- performance(pr[[i]],measure="auc")
#    auc[[i]] <- auc[[i]]@y.values[[1]]
#}
#
#mycols <- brewer.pal(4,"Set1")
##Plot the ROCR curve
##pdf("rocr_curves.pdf",width=10,height=10,paper="special")
#png("rocr_curves.png")
#plot(prf[[1]],type="n")
#for (i in seq_along(fits.test)){
#    lines(prf[[i]]@x.values[[1]],prf[[i]]@y.values[[1]],
#          lwd=5,col=mycols[i])
#}
#legend("bottomright",legend=paste(names(fits.test),
#                         " (auc = ",
#                         format(auc,digits=2),")",sep=""),
#                         inset=0.1,lwd=5,col=mycols)
#box()
#dev.off()
#
#
#
##Also fit the predicted values across the range of the observed values
##This is to plot the logistic fit
#xgpats <- seq(0,10000,length.out=1000)
##xship <- seq(0,4,length.out=1000)
#xship <- seq(0,10,length.out=1000)
#fit.ship <- predict(glm.ship,
#                       list(Max_SHIP=xship),
#                       type="response")
#
#fit.shipmg <- predict(glm.shipmg,
#                      list(SHIP_Max_GPATS=xship),
#                      type="response")
#
#fit.gpats <- predict(glm.gpats,
#                     list(Max_GPATS=xgpats),
#                     type="response")
#
#fit.shipmg.gpats <- predict(glm.shipmg.gpats,
#                       list(SHIP_Max_GPATS=xship,Max_GPATS=xgpats),
#                       type="response")
#
#
##pdf("logistic_regression_ship_gpats.pdf",
##    width=10,height=10,paper="special")
#png("logistic_regression_ship_gpats%03d.png")
##par(mfrow=c(2,1))
#plot(data$Max_SHIP, data$HAIL,pch=16,col="grey20",
#     xlab="SHIP", ylab="Probability of Hail",
#     xlim=c(0,7))
#lines(xship,fit.ship,lwd=5,col=mycols[1])
#lines(xship,fit.shipmg,lwd=5,col=mycols[2])
#lines(xship,fit.shipmg.gpats,lwd=5,col=mycols[4])
#legend("topright",
#       legend=paste(c("SHIP", "SMG", "SMG+GPATS"),
#                    format(c(accuracy[1],accuracy[2],accuracy[3]),
#                             digits=2),sep=" "),
#       lwd=5,col=mycols[c(1,2,4)],inset=0.1)
#
#plot(data$Max_GPATS,data$HAIL,pch=16,col="grey20",
#     xlab="Lightning strikes / hour",
#     ylab="Probability of hail",
#     xlim=c(0,10000))
#lines(xgpats,fit.gpats,lwd=5,col=mycols[3])
#lines(xgpats,fit.shipmg.gpats,lwd=5,col=mycols[4])
#legend("topright",
#       legend=paste(c("GPATS", "SMG+GPATS"),
#                    format(c(accuracy[3],accuracy[4]),
#                           digits=2),sep=" "),
#       lwd=5,col=mycols[c(3,4)],inset=c(0.05,0.1))
#dev.off()
#
