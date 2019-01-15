setwd("/Volumes/jwillou/masterbayes/analyze_WH/")
library(MASS)
library(scales)
library(rstanarm)
library(shinystan)
library(pscl)
library(AER)
library(lme4)
library(languageR)
library(nlme)
library(piecewiseSEM)
library(effects)

#data = read.table("outputbestped0/summary_allstats2.csv", header=TRUE, sep=",")
data = read.table("outputbestped0/summary_allstats2_withgen.csv", sep=",", header=TRUE)
cols = colnames(data)
colors4 = c("olivedrab3", "deepskyblue2") #F, K, 
colors3 = c("navy","dodgerblue3","turquoise3", "grey50")
colors2 = c("olivedrab3", "deepskyblue2", "navy")
colors1 = c("darkorchid3", "firebrick3", "dodgerblue3", "grey30")
#pdata = data.frame(id=data$id, momid = data$momid, dadid=data$dadid, lifespan=data$lifespan, off=data$off, off_survive=data$off_survive)
#write.table(pdata, "purgedata.csv", row.names=FALSE, col.names=TRUE, sep=",")

####gensknown####
# library(kinship2)
# toped2  = data.frame(id=data$id, momid=data$momid, dadid=data$dadid, sex=data$sex)
# rownames(toped2) = seq(1,nrow(toped2),1)
# unassign = as.numeric(as.character(unique(toped2$id[is.na(toped2$sex)])))
# for(u in 1:nrow(toped2)){
#   if(!is.na(match(toped2$momid[u], unassign))){
#     toped2$momid[u] = NA
#     next
#   }
#   if(!is.na(match(toped2$dadid[u], unassign))){
#     toped2$dadid[u] = NA
#   }
# }
# toped2 = toped2[!is.na(toped2$sex),]
# levels(toped2$sex) = c(levels(toped2$sex), "female", "male")
# toped2[toped2$sex=="Female", 4] = "female"
# toped2[toped2$sex=="Male",   4] = "male"  
# 
# #check that either both dad and mom are missing or that both are present
# for(r in 1:nrow(toped2)){
#   if(is.na(toped2$dadid[r]) | is.na(toped2$momid[r]) | toped2$dadid[r]==0 | toped2$momid[r]==0){
#     toped2[r, 2:3] = "us"
#   }
# }
# pedord = pedigree::orderPed(toped2)
# toped2  = toped2[order(pedord),]
# ped = kinship2::pedigree(id=toped2$id, dadid=toped2$dadid, momid=toped2$momid, sex=toped2$sex, missid="us") 
# k = kindepth(id=toped2$id,mom.id=toped2$momid,dad.id=toped2$dadid)
# df = data.frame(id=toped2$id, genback=k)
# data$genback = NULL
# data = merge(data, df, by="id")
# write.table(data, "outputbestped0/summary_allstats2_withgen.csv", sep=",", row.names=FALSE, col.names=TRUE)

par(mfrow=c(1,1))
x = data[data$genback>0,]
hist(x$genback, ylim=c(0,300), xlim=c(0,20), breaks=seq(0,20,1), col=colors1[2], xlab="pedigree depth", main=NA)
mean(x$genback, na.rm=TRUE)

####popchar####
par(mfrow=c(1,2))
x = data[data$genback>cutoff,]
x = x[x$birthyear<=2005,]
x = x[!is.na(x[,1]),]
hist(x$inb, ylim=c(0,250), xlim=c(0,0.5), breaks=seq(0,0.5, 0.01), col=colors2[1], xlab="inbreeding coefficient", main=NA)
hist(x$pkins, ylim=c(0,80), xlim=c(0,0.5), breaks=seq(0,0.5, 0.01), col=colors2[2], xlab="mean kinship between mates", main=NA)

####number of mates####
par(mfrow=c(1,2))
x = data[data$genback>cutoff,]
x = x[x$birthyear<=2005,]
x = x[!is.na(x[,1]),]
hist(x$nmates[x$sex=="Female"], ylim=c(0,60), xlim=c(0,15), breaks=seq(0,15, 1), col=colors1[2], xlab="lifetime number of mates", main=NA)
hist(x$nmates[x$sex=="Male"], ylim=c(0,60), xlim=c(0,15), breaks=seq(0,15, 1), col=colors1[3], xlab="lifetime number of mates", main=NA)
mean(x$nmates[x$sex=="Female"], na.rm=TRUE)
mean(x$nmates[x$sex=="Male"], na.rm=TRUE)

####var in reprosuccess####
par(mfrow=c(1,3))
x = data[data$genback>cutoff,]
x = x[x$birthyear<=2005,]
x = x[x$off>0,]
x = x[!is.na(x[,1]),]
hist(x$off, xlab="number of offspring", ylab="frequency", main=NA, breaks=seq(1, 17, 1), xlim=c(0,20), ylim=c(0,150), col=colors1[1])
hist(x$off[x$sex=="Female"], xlab="number of offspring", ylab="frequency", main=NA, breaks=seq(1, 17, 1), ylim=c(0,80), xlim=c(0,20), col=colors1[2])
hist(x$off[x$sex=="Male"], xlab="number of offspring", ylab="frequency", main=NA, breaks=seq(1, 17, 1), ylim=c(0,80), xlim=c(0,20), col=colors1[3])
