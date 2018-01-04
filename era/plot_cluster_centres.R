#Plot the cluster centres
#Justin Peter, ICACS, USQ
#2 June 2017

library(ggplot2)
library(reshape) #reshape has to be loaded after ggplot2

fnm <- "Cluster_information_10.txt"

df <- read.table(fnm,skip=31)
df$cls <- factor(seq(1:10))
df$lcape <- log(df$MUCAPE)

dfm <- melt(df,id.vars="cls")

dfm$variable_f <- factor(dfm$variable,levels=c("lcape","MUCAPE","S06","PMR","H5_TEMP","LR75"))
levels(dfm$variable_f) <-c("log(MUCAPE)","MUCAPE","S06","PMR","H5_TEMP","LR75")


#variable_names <- list("lcape"="log(MUCAPE)",
#                       "LR75"="LR75",
#                       "H5_TEMP"="H5_TEMP",
#                       "S06"="S06",
#                       "PMR"="PMR")
#variable_labeller <- function(variable_f,value){
#    return(variable_names[value])
#}

pdf("Cluster_centres_10.pdf",
    width=10,height=10,paper="special")
p <- ggplot(subset(dfm, variable %in% c("lcape","LR75","H5_TEMP","S06","PMR")),
       aes(x=cls,y=value))+
       geom_point(size=3)+
       geom_line(aes(group=1))+
       #facet_grid(variable_f~.,scales="free",labeller=variable_labeller)+
       facet_grid(variable_f~.,scales="free")+
       xlab("Cluster")+ylab("Centre values")
print(p)
dev.off()


#
#
#coef.plot <- ggplot(df,aes(x=cls,y=Coeff))+
#                    facet_grid(var~.,scales="free")+
#                    geom_point()+
#                    geom_errorbar(aes(ymin=lower,ymax=upper),width=0.1)+
#                    xlab("Cluster")+ylab("Coefficient")
#


