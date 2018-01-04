#Read in SHIP/GPATS/HAIL occurrence data
#an apply a logistic regression
#See the following for references
#http://www.r-bloggers.com/how-to-perform-a-logistic-regression-in-r/

rm(list=ls())
#library(ggplot2)
library(dplyr)
#library(HSAUR2)
library(pscl)
library(ROCR)
library(RColorBrewer)


data <- read.table("ship_lightning_hail_frequency.txt",
                   header=TRUE,
                   stringsAsFactors=FALSE)
#Turn the hail occurrence into a factor
#data$Hail_Occurred <- as.factor(data$Hail_Occurred)
#Use a shorter name
names(data)[10] <- "HAIL"
#Replace no/yes with 0/1
data$HAIL <- ifelse(data$HAIL=="Yes",1,0)
data$loggpats <- log(data$Max_GPATS)


#Extract the values where lightning was observed
#ltng.threshold <- 1000
#ltng.threshold <- 100
#ltng.threshold <- 10
#ltng.threshold <- 1
ltng.threshold <- 0
ltng.obs <- which(data$Max_GPATS >= ltng.threshold)
data <- data[ltng.obs,]

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


##Develop the logistic model on the training data set
##GLM 1 - Hail Damage
#gpats_glm <- glm(Hail_Occurred ~ Max_GPATS,family=binomial())
##This is probability of hail vs max. ship on the day
#ship_glm_1  <- glm(data$Hail_Occurred ~ data$Max_SHIP,family=binomial())
##This is probability of hail vs. ship at max. gpats
#ship_glm_2  <- glm(data$Hail_Occurred ~ data$SHIP_Max_GPATS,faily=binomial())

# Now test predictive ability on the test data set
glm.ship       <- glm(HAIL ~ Max_SHIP,family=binomial(),data=train)
glm.shipmg     <- glm(HAIL ~ SHIP_Max_GPATS,family=binomial(),data=train)
glm.gpats      <- glm(HAIL ~ Max_GPATS, family=binomial(),data=train)
#glm.log.gpats  <- glm(HAIL ~ log(Max_GPATS), family=binomial(),data=train)
#glm.shipmg.log.gpats <- glm(HAIL ~ log(Max_GPATS)+SHIP_Max_GPATS,family=binomial(),data=train)
glm.shipmg.gpats <- glm(HAIL ~ SHIP_Max_GPATS+Max_GPATS,family=binomial(),data=train)
#glm.shipmg.gpats <- glm(HAIL ~ SHIP_Max_GPATS*Max_GPATS,family=binomial(),data=train)

#params <- c("Max_SHIP","SHIP_Max_GPATS","Max_GPATS",
#            #"log(Max_GPATS)","log(Max_GPATS)/SHIP_Max_GPATS")
#            "SHIP_MG_GPATS")


#Put all of the glms into a list for loop processing
glms <- list(glm.ship,
             glm.shipmg,
             glm.gpats,
             #glm.log.gpats,
             #glm.ship.log.gpats,
             glm.shipmg.gpats)
names(glms) <- c("SHIP","SHIP_MG","GPATS","SHIP_MG_GPATS")

for (i in 1:length(glms)){
    summary(glms[[i]])
    anova(glms[[i]],test="Chisq")
    confint(glms[[i]])
    print(exp(coef(glms[[i]])))
    print(exp(confint(glms[[i]])))
}

#model_vil <- glm(train$HAIL ~ train$VIL,family=binomial(link=logit))
#See http://www.theanalysisfactor.com/r-tutorial-glm1/ for why you have to use
#the lower model description. I have to look into this!
#model_vil <- glm(HAIL ~ VIL,family=binomial, data=test)

#Compute some pseudo-R2 measures for the GLM fits
pr2 <- vector("list",length=length(glms))
for (i in 1:length(glms)){
    pr2[[i]]<-pR2(glms[[i]])
}

#Now fit the results to the test case
#Which we nee to evaluate the accuracy of the model

fit.ship.test <- predict(glm.ship,
                       newdata=subset(test,select=c(7,10)),
                       type="response")

fit.shipmg.test <- predict(glm.shipmg,
                              newdata=subset(test,select=c(8,10)),
                              type="response")

fit.gpats.test <- predict(glm.gpats,
                      newdata=subset(test,select=c(9,10)),
                      type="response")

fit.shipmg.gpats.test <- predict(glm.shipmg.gpats,
                       newdata=subset(test,select=c(8,9,10)),
                       type="response")

#Evaluate the accuracy of each model in determining haail presence
#Requires turning glm model into a binary response
#(>0.5 --> hail and <= 0.5 --> no hail)
fits.test<-list(fit.ship.test,fit.shipmg.test,fit.gpats.test,fit.shipmg.gpats.test)
names(fits.test) <- c("SHIP","SHIP_MG","GPATS","SHIP_MG_GPATS")
fits.test.bi <- vector("list",length=length(fits.test))
misClasificError      <- vector("numeric",length=length(fits.test))
accuracy              <- vector("numeric",length=length(fits.test))
for (i in seq_along(fits.test)){
    fits.test.bi[[i]] <- ifelse(fits.test[[i]] > 0.5,1,0)
    misClasificError[i] <- mean(fits.test.bi[[i]] != test$HAIL)
    #misClasificError <- mean(fitted.results.vil.bi != test$HAIL)
    #print(paste(names(test)[i],'Accuracy',1-misClasificError[[i]]))
    accuracy[i] <- 1.-misClasificError[i]
}

#Examine the ROCR curves for the models
#See https://www.r-bloggers.com/how-to-perform-a-logistic-regression-in-r/
pr  <- vector("list",length=length(fits.test))
prf <- vector("list",length=length(fits.test))
for (i in seq_along(fits.test)){
    pr[[i]]  <- prediction(fits.test[[i]],test$HAIL)
    prf[[i]] <- performance(pr[[i]],measure="tpr",x.measure="fpr")
}
#Some performance measures
auc <- vector("list",length=length(fits.test))
for (i in seq_along(fits.test)){
    auc[[i]] <- performance(pr[[i]],measure="auc")
    auc[[i]] <- auc[[i]]@y.values[[1]]
}

mycols <- brewer.pal(4,"Set1")
#Plot the ROCR curve
#pdf("rocr_curves.pdf",width=10,height=10,paper="special")
png("rocr_curves.png")
plot(prf[[1]],type="n")
for (i in seq_along(fits.test)){
    lines(prf[[i]]@x.values[[1]],prf[[i]]@y.values[[1]],
          lwd=5,col=mycols[i])
}
legend("bottomright",legend=paste(names(fits.test),
                         " (auc = ",
                         format(auc,digits=2),")",sep=""),
                         inset=0.1,lwd=5,col=mycols)
box()
dev.off()



#Also fit the predicted values across the range of the observed values
#This is to plot the logistic fit
xgpats <- seq(0,10000,length.out=1000)
#xship <- seq(0,4,length.out=1000)
xship <- seq(0,10,length.out=1000)
fit.ship <- predict(glm.ship,
                       list(Max_SHIP=xship),
                       type="response")

fit.shipmg <- predict(glm.shipmg,
                      list(SHIP_Max_GPATS=xship),
                      type="response")

fit.gpats <- predict(glm.gpats,
                     list(Max_GPATS=xgpats),
                     type="response")

fit.shipmg.gpats <- predict(glm.shipmg.gpats,
                       list(SHIP_Max_GPATS=xship,Max_GPATS=xgpats),
                       type="response")


#pdf("logistic_regression_ship_gpats.pdf",
#    width=10,height=10,paper="special")
png("logistic_regression_ship_gpats%03d.png")
#par(mfrow=c(2,1))
plot(data$Max_SHIP, data$HAIL,pch=16,col="grey20",
     xlab="SHIP", ylab="Probability of Hail",
     xlim=c(0,7))
lines(xship,fit.ship,lwd=5,col=mycols[1])
lines(xship,fit.shipmg,lwd=5,col=mycols[2])
lines(xship,fit.shipmg.gpats,lwd=5,col=mycols[4])
legend("topright",
       legend=paste(c("SHIP", "SMG", "SMG+GPATS"),
                    format(c(accuracy[1],accuracy[2],accuracy[3]),
                             digits=2),sep=" "),
       lwd=5,col=mycols[c(1,2,4)],inset=0.1)

plot(data$Max_GPATS,data$HAIL,pch=16,col="grey20",
     xlab="Lightning strikes / hour",
     ylab="Probability of hail",
     xlim=c(0,10000))
lines(xgpats,fit.gpats,lwd=5,col=mycols[3])
lines(xgpats,fit.shipmg.gpats,lwd=5,col=mycols[4])
legend("topright",
       legend=paste(c("GPATS", "SMG+GPATS"),
                    format(c(accuracy[3],accuracy[4]),
                           digits=2),sep=" "),
       lwd=5,col=mycols[c(3,4)],inset=c(0.05,0.1))
dev.off()

#plot(test$Max_SHIP,fit.ship.test,
#     xlab="Max. SHIP",
#     ylab="Probability of Hail",
#     ylim=c(0,1),type="n")
#lines(sort(test$Max_SHIP),sort(fit.ship.test),
#      lty=2,lwd=5)
#points(data$Max_SHIP,data$HAIL,col="grey20")
#
#plot(test$SHIP_Max_GPATS,fit.shipmg.test,
#     xlab="Max. SHIP",
#     ylab="Probability of Hail",
#     ylim=c(0,1),type="n")
#lines(sort(test$SHIP_Max_GPATS),sort(fit.shipmg.test),
#      lty=2,lwd=5)
#points(data$Max_SHIP,data$HAIL,col="grey20")
#
#plot(test$Max_GPATS, fit.gpats.test,
#     xlab="Max. GPATS",
#     ylab="Probability of Hail",
#     ylim=c(0,1),type="n")
#lines(sort(test$Max_GPATS),sort(fit.gpats.test),
#      lty=2,lwd=5)
#points(data$Max_GPATS,data$HAIL,col="grey20")
#
#plot(log(test$Max_GPATS), fit.log.gpats.test,
#     xlab="log(Max. GPATS)",
#     ylab="Probability of Hail",
#     ylim=c(0,1),type="n")
#lines(sort(log(test$Max_GPATS)),sort(fit.log.gpats.test),
#      lty=2,lwd=5)
#points(log(data$Max_GPATS),data$HAIL,col="grey20")
#
#plot(test$Max_GPATS, fit.shipmg.gpats.test,
#     xlab="Max. GPATS",
#     ylab="Probability of Hail",
#     ylim=c(0,1),type="n")
#lines(sort(test$Max_GPATS),sort(fit.shipmg.gpats.test),
#      lty=2,lwd=5)
#points(data$Max_GPATS,data$HAIL,col="grey20")
#plot(test$SHIP_Max_GPATS, fit.shipmg.gpats.test,
#     xlab="SHIP at Max. GPATS",
#     ylab="Probability of Hail",
#     ylim=c(0,1),type="n")
#lines(sort(test$SHIP_Max_GPATS),sort(fit.shipmg.gpats.test),
#      lty=2,lwd=5)
#points(data$SHIP_Max_GPATS,data$HAIL,col="grey20")


