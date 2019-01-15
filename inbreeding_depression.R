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

cutoff = 3 #number generations above parents in data (kinships are from parents back, so add 1 to make total)
data = read.table("outputbestped0/summary_allstats2_withgen.csv", sep=",", header=TRUE)
kinM = read.table(paste("outputbestped0/kinship_moms_c", cutoff, ".csv", sep=""), header=TRUE, sep=",")
kinD = read.table(paste("outputbestped0/kinship_dads_c", cutoff, ".csv", sep=""), header=TRUE, sep=",")
kins = rbind(kinM, kinD)
kdata = merge(data, kins, by="id")
cols = colnames(data)
colors4 = c("olivedrab3", "deepskyblue2") #F, K, 
colors3 = c("navy","dodgerblue3","turquoise3", "grey50")
colors2 = c("olivedrab3", "deepskyblue2", "navy")
colors1 = c("darkorchid3", "firebrick3", "dodgerblue3", "grey30")
#pdata = data.frame(id=data$id, momid = data$momid, dadid=data$dadid, lifespan=data$lifespan, off=data$off, off_survive=data$off_survive)
#write.table(pdata, "purgedata.csv", row.names=FALSE, col.names=TRUE, sep=",")

# p = data.frame(id=data$id, parent_1=data$momid, parent_2=data$dadid, fitness=data$off, generation=data$birthyear-min(data$birthyear, na.rm=T))
# pedord = pedigree::orderPed(p)
# p  = p[order(pedord),]
# p$parent_1[is.na(p$parent_1)] = 0
# p$parent_2[is.na(p$parent_2)] = 0
# write.table(p, "~/Desktop/krat1.csv", row.names=F, col.names=T, sep=",")

####F/K to noffspring####
tdata = kdata[kdata$genback>cutoff,]
tdata = tdata[tdata$birthyear<=2005,]
tdata = tdata[!is.na(tdata[,1]),]
tdata$lifespan = tdata$deathyear-tdata$birthyear

par(mfrow=c(1,1))

#inb
df = data.frame(x = tdata$inb,y = tdata$off, year=tdata$birthyear, sex=tdata$sex, k=tdata$pkins)
df = df[complete.cases(df),]
df = df[df$x<0.5,]
x = df$x
y = df$y
s = df$sex
k = df$k
year = df$year
modelout = glmer(y~x+(1|year), family="poisson")
summary(modelout)
sem.model.fits(modelout)

pdf("~/Desktop/F_off.pdf", width=5, height=5, useDingbats=FALSE, onefile=TRUE)
plot(-100,-100, xlab="inbreeding coefficient", ylab="number of offspring", ylim=c(0, 15), xlim=c(0,0.4), main="off")
plotLMER.fnc(modelout, linecolor=colors4[1], lwd=4, pred="x", n=100, addToExistingPlot=TRUE, fun=exp)
points(x=x, y=y, pch=21, col=alpha("grey30", 0.5), bg=alpha("grey30", 0.5), cex=0.75)
dev.off()

#how many more
cof = mean(df$x) *2
mean(df$y[df$x>cof])/mean(df$y)

#pkins
tdata = kdata
tdata = tdata[tdata$birthyear<=2005,]
tdata = tdata[!is.na(tdata[,1]),]
tdata$lifespan = tdata$deathyear-tdata$birthyear
df = data.frame(x = tdata$pkins,y = tdata$off, year=tdata$birthyear, sex=tdata$sex, i=tdata$inb)
df = df[complete.cases(df),]
df = df[df$x<0.5,]
x = df$x
y = df$y
s = df$sex
i = df$i
year = df$year
modelout = glmer(y~x+(1|year), family="poisson")
summary(modelout)
sem.model.fits(modelout)

pdf("~/Desktop/K_off.pdf", width=5, height=5, useDingbats=FALSE, onefile=TRUE)
plot(-100,-100, xlab="mean mate kinship", ylab="number of offspring", ylim=c(0, 15), xlim=c(0,0.4), main="off")
plotLMER.fnc(modelout, linecolor=colors4[2], lwd=4, pred="x", n=100, addToExistingPlot=TRUE, fun=exp)
points(x=x, y=y, pch=21, col=alpha("grey30", 0.5), bg=alpha("grey30", 0.5), cex=0.75)
dev.off()

#how many more
#cutoff = mean(df$x) *2
#1 - mean(df$y[df$x>cutoff])/mean(df$y)

####F/K to noffspring - males only####
tdata = kdata[kdata$genback>cutoff,]
tdata = tdata[tdata$birthyear<=2005,]
tdata = tdata[tdata$sex=="Male",]
tdata = tdata[!is.na(tdata[,1]),]
tdata$lifespan = tdata$deathyear-tdata$birthyear

par(mfrow=c(1,1))

#inb
df = data.frame(x = tdata$inb,y = tdata$off, year=tdata$birthyear)
df = df[complete.cases(df),]
df = df[df$x<0.5,]
x = df$x
y = df$y
year = df$year
modelout = glmer(y~x+(1|year), family="poisson")
summary(modelout)
sem.model.fits(modelout)

pdf("~/Desktop/F_off_males.pdf", width=5, height=5, useDingbats=FALSE, onefile=TRUE)
plot(-100,-100, xlab="inbreeding coefficient", ylab="number of offspring", ylim=c(0, 15), xlim=c(0,0.4), main="off_males")
plotLMER.fnc(modelout, linecolor=colors4[1], lwd=4, pred="x", n=100, addToExistingPlot=TRUE, fun=exp)
points(x=x, y=y, pch=21, col=alpha("grey30", 0.5), bg=alpha("grey30", 0.5), cex=0.75)
dev.off()

#how many more
#cutoff = mean(df$x) *2
#mean(df$y[df$x>cutoff])/mean(df$y)

#pkins
tdata = kdata
tdata = tdata[tdata$birthyear<=2005,]
tdata = tdata[tdata$sex=="Male",]
tdata = tdata[!is.na(tdata[,1]),]
tdata$lifespan = tdata$deathyear-tdata$birthyear
df = data.frame(x = tdata$pkins,y = tdata$off, year=tdata$birthyear)
df = df[complete.cases(df),]
df = df[df$x<0.5,]
x = df$x
y = df$y
year = df$year
modelout = glmer(y~x+(1|year), family="poisson")
summary(modelout)
sem.model.fits(modelout)

pdf("~/Desktop/K_off_males.pdf", width=5, height=5, useDingbats=FALSE, onefile=TRUE)
plot(-100,-100, xlab="mean mate kinship", ylab="number of offspring", ylim=c(0, 15), xlim=c(0,0.4), main="off_males")
plotLMER.fnc(modelout, linecolor=colors4[2], lwd=4, pred="x", n=100, addToExistingPlot=TRUE, fun=exp)
points(x=x, y=y, pch=21, col=alpha("grey30", 0.5), bg=alpha("grey30", 0.5), cex=0.75)
dev.off()

#how many more
#cutoff = mean(df$x) *2
#1 - mean(df$y[df$x>cutoff])/mean(df$y)

####F/K to noffspring - females only####
tdata = kdata[kdata$genback>cutoff,]
tdata = tdata[tdata$birthyear<=2005,]
tdata = tdata[tdata$sex=="Female",]
tdata = tdata[!is.na(tdata[,1]),]
tdata$lifespan = tdata$deathyear-tdata$birthyear

par(mfrow=c(1,1))

#inb
df = data.frame(x = tdata$inb,y = tdata$off, year=tdata$birthyear)
df = df[complete.cases(df),]
df = df[df$x<0.5,]
x = df$x
y = df$y
year = df$year
modelout = glmer(y~x+(1|year), family="poisson")
summary(modelout)
sem.model.fits(modelout)

pdf("~/Desktop/F_off_females.pdf", width=5, height=5, useDingbats=FALSE, onefile=TRUE)
plot(-100,-100, xlab="inbreeding coefficient", ylab="number of offspring", ylim=c(0, 15), xlim=c(0,0.4), main="off_females")
plotLMER.fnc(modelout, linecolor=colors4[1], lwd=4, pred="x", n=100, addToExistingPlot=TRUE, fun=exp)
points(x=x, y=y, pch=21, col=alpha("grey30", 0.5), bg=alpha("grey30", 0.5), cex=0.75)
dev.off()

#how many more
#cutoff = mean(df$x) *2
#mean(df$y[df$x>cutoff])/mean(df$y)

#pkins
tdata = kdata
tdata = tdata[tdata$birthyear<=2005,]
tdata = tdata[tdata$sex=="Female",]
tdata = tdata[!is.na(tdata[,1]),]
tdata$lifespan = tdata$deathyear-tdata$birthyear
df = data.frame(x = tdata$pkins,y = tdata$off, year=tdata$birthyear)
df = df[complete.cases(df),]
df = df[df$x<0.5,]
x = df$x
y = df$y
year = df$year
modelout = glmer(y~x+(1|year), family="poisson")
summary(modelout)
sem.model.fits(modelout)

pdf("~/Desktop/K_off_females.pdf", width=5, height=5, useDingbats=FALSE, onefile=TRUE)
plot(-100,-100, xlab="mean mate kinship", ylab="number of offspring", ylim=c(0, 15), xlim=c(0,0.4), main="off_females")
plotLMER.fnc(modelout, linecolor=colors4[2], lwd=4, pred="x", n=100, addToExistingPlot=TRUE, fun=exp)
points(x=x, y=y, pch=21, col=alpha("grey30", 0.5), bg=alpha("grey30", 0.5), cex=0.75)
dev.off()

#how many more
#cutoff = mean(df$x) *2
#1 - mean(df$y[df$x>cutoff])/mean(df$y)

####F/K to noffspring survive to age 1####
tdata = kdata[kdata$genback>cutoff,]
tdata = tdata[tdata$birthyear<=2005,]
tdata = tdata[!is.na(tdata[,1]),]
tdata$lifespan = tdata$deathyear-tdata$birthyear

par(mfrow=c(1,1))

#inb
df = data.frame(x = tdata$inb,y = tdata$off_survive, year=tdata$birthyear)
df = df[complete.cases(df),]
df = df[df$x<0.5,]
x = df$x
y = df$y
year = df$year
modelout = glmer(y~x+(1|year), family="poisson")
summary(modelout)
sem.model.fits(modelout)

pdf("~/Desktop/F_offsurv.pdf", width=5, height=5, useDingbats=FALSE, onefile=TRUE)
plot(-100,-100, xlab="inbreeding coefficient", ylab="number of offspring", ylim=c(0, 15), xlim=c(0,0.4), main="off_survive")
plotLMER.fnc(modelout, linecolor=colors4[1], lwd=4, pred="x", n=100, addToExistingPlot=TRUE, fun=exp)
points(x=x, y=jitter(y, 0.25), pch=21, col=alpha("grey30", 0.5), bg=alpha("grey30", 0.5), cex=0.75)
dev.off()

#how many more
cof = mean(df$x) *2
1-mean(df$y[df$x>cof])/mean(df$y)

#pkins
tdata = kdata
tdata = tdata[tdata$birthyear<=2005,]
tdata = tdata[!is.na(tdata[,1]),]
tdata$lifespan = tdata$deathyear-tdata$birthyear
df = data.frame(x = tdata$pkins,y = tdata$off_survive, year=tdata$birthyear)
df = df[complete.cases(df),]
df = df[df$x<0.5,]
x = df$x
y = df$y
year = df$year
modelout = glmer(y~x+(1|year), family="poisson")
summary(modelout)
sem.model.fits(modelout)

pdf("~/Desktop/K_offsurv.pdf", width=5, height=5, useDingbats=FALSE, onefile=TRUE)
plot(-100,-100, xlab="mean mate kinship", ylab="number of offspring", ylim=c(0, 15), xlim=c(0,0.4), main="off_survive")
plotLMER.fnc(modelout, linecolor=colors4[2], lwd=4, pred="x", n=100, addToExistingPlot=TRUE, fun=exp)
points(x=x, y=y, pch=21, col=alpha("grey30", 0.5), bg=alpha("grey30", 0.5), cex=0.75)
dev.off()

#how many more
#cutoff = mean(df$x) *2
#1 - mean(df$y[df$x>cutoff])/mean(df$y)

####F/K to noffspring survive to age 1 - males only####
tdata = kdata[kdata$genback>cutoff,]
tdata = tdata[tdata$birthyear<=2005,]
tdata = tdata[tdata$sex=="Male",]
tdata = tdata[!is.na(tdata[,1]),]
tdata$lifespan = tdata$deathyear-tdata$birthyear

par(mfrow=c(1,1))

#inb
df = data.frame(x = tdata$inb,y = tdata$off_survive, year=tdata$birthyear)
df = df[complete.cases(df),]
df = df[df$x<0.5,]
x = df$x
y = df$y
year = df$year
modelout = glmer(y~x+(1|year), family="poisson")
summary(modelout)
sem.model.fits(modelout)

pdf("~/Desktop/F_offsurv_males.pdf", width=5, height=5, useDingbats=FALSE, onefile=TRUE)
plot(-100,-100, xlab="inbreeding coefficient", ylab="number of offspring", ylim=c(0, 15), xlim=c(0,0.4), main="off_survive_males")
plotLMER.fnc(modelout, linecolor=colors4[1], lwd=4, pred="x", n=100, addToExistingPlot=TRUE, fun=exp)
points(x=x, y=jitter(y, 0.25), pch=21, col=alpha("grey30", 0.5), bg=alpha("grey30", 0.5), cex=0.75)
dev.off()

#how many more
#cutoff = mean(df$x) *2
#mean(df$y[df$x>cutoff])/mean(df$y)

#pkins
tdata = kdata
tdata = tdata[tdata$birthyear<=2005,]
tdata = tdata[tdata$sex=="Male",]
tdata = tdata[!is.na(tdata[,1]),]
tdata$lifespan = tdata$deathyear-tdata$birthyear
df = data.frame(x = tdata$pkins,y = tdata$off_survive, year=tdata$birthyear)
df = df[complete.cases(df),]
df = df[df$x<0.5,]
x = df$x
y = df$y
year = df$year
modelout = glmer(y~x+(1|year), family="poisson")
summary(modelout)
sem.model.fits(modelout)

pdf("~/Desktop/K_offsurv_males.pdf", width=5, height=5, useDingbats=FALSE, onefile=TRUE)
plot(-100,-100, xlab="mean mate kinship", ylab="number of offspring", ylim=c(0, 15), xlim=c(0,0.4), main="off_survive_males")
plotLMER.fnc(modelout, linecolor=colors4[2], lwd=4, pred="x", n=100, addToExistingPlot=TRUE, fun=exp)
points(x=x, y=y, pch=21, col=alpha("grey30", 0.5), bg=alpha("grey30", 0.5), cex=0.75)
dev.off()

#how many more
#cutoff = mean(df$x) *2
#1 - mean(df$y[df$x>cutoff])/mean(df$y)

####F/K to noffspring survive to age 1 - females only####
tdata = kdata[kdata$genback>cutoff,]
tdata = tdata[tdata$birthyear<=2005,]
tdata = tdata[tdata$sex=="Female",]
tdata = tdata[!is.na(tdata[,1]),]
tdata$lifespan = tdata$deathyear-tdata$birthyear

par(mfrow=c(1,1))

#inb
df = data.frame(x = tdata$inb,y = tdata$off_survive, year=tdata$birthyear)
df = df[complete.cases(df),]
df = df[df$x<0.5,]
x = df$x
y = df$y
year = df$year
modelout = glmer(y~x+(1|year), family="poisson")
summary(modelout)
sem.model.fits(modelout)

pdf("~/Desktop/F_offsurv_females.pdf", width=5, height=5, useDingbats=FALSE, onefile=TRUE)
plot(-100,-100, xlab="inbreeding coefficient", ylab="number of offspring", ylim=c(0, 15), xlim=c(0,0.4), main="off_survive_females")
plotLMER.fnc(modelout, linecolor=colors4[1], lwd=4, pred="x", n=100, addToExistingPlot=TRUE, fun=exp)
points(x=x, y=jitter(y, 0.25), pch=21, col=alpha("grey30", 0.5), bg=alpha("grey30", 0.5), cex=0.75)
dev.off()

#how many more
#cutoff = mean(df$x) *2
#mean(df$y[df$x>cutoff])/mean(df$y)

#pkins
tdata = kdata
tdata = tdata[tdata$birthyear<=2005,]
tdata = tdata[tdata$sex=="Female",]
tdata = tdata[!is.na(tdata[,1]),]
tdata$lifespan = tdata$deathyear-tdata$birthyear
df = data.frame(x = tdata$pkins,y = tdata$off_survive, year=tdata$birthyear)
df = df[complete.cases(df),]
df = df[df$x<0.5,]
x = df$x
y = df$y
year = df$year
modelout = glmer(y~x+(1|year), family="poisson")
summary(modelout)
sem.model.fits(modelout)

pdf("~/Desktop/K_offsurv_females.pdf", width=5, height=5, useDingbats=FALSE, onefile=TRUE)
plot(-100,-100, xlab="mean mate kinship", ylab="number of offspring", ylim=c(0, 15), xlim=c(0,0.4), main="off_survive_females")
plotLMER.fnc(modelout, linecolor=colors4[2], lwd=4, pred="x", n=100, addToExistingPlot=TRUE, fun=exp)
points(x=x, y=y, pch=21, col=alpha("grey30", 0.5), bg=alpha("grey30", 0.5), cex=0.75)
dev.off()

#how many more
#cutoff = mean(df$x) *2
#1 - mean(df$y[df$x>cutoff])/mean(df$y)

####F to lifespan####
tdata = kdata[kdata$genback>cutoff,]
tdata = tdata[tdata$birthyear<=2005,]
tdata = tdata[!is.na(tdata[,1]),]
tdata$lifespan = tdata$deathyear-tdata$birthyear
par(mfrow=c(1,1))

#inb
df = data.frame(x = tdata$inb,y = tdata$lifespan, year=tdata$birthyear)
df = df[complete.cases(df),]
df = df[df$x<0.5,]
x = df$x #tdata$inb
y = df$y 
year = df$year
modelout = glmer(y~x+(1|year), family="poisson")
summary(modelout)
sem.model.fits(modelout)

pdf("~/Desktop/F_life.pdf", width=5, height=5, useDingbats=FALSE, onefile=TRUE)
plot(-100,-100, xlab="inbreeding coefficient", ylab="lifespan", ylim=c(0.75, 6), xlim=c(0,0.4))
plotLMER.fnc(modelout, linecolor=colors4[1], lwd=4, pred="x", n=100, addToExistingPlot=TRUE, fun=exp)
points(x=x, y=jitter(y, 0.25), pch=21, col=alpha("grey30", 0.5), bg=alpha("grey30", 0.5), cex=0.75)
dev.off()

#how much shorter
#cutoff = mean(df$x) *2
#mean(df$y[df$x>cutoff])-mean(df$y)

####F to lifespan - males only####
tdata = kdata[kdata$genback>cutoff,]
tdata = tdata[tdata$birthyear<=2005,]
tdata = tdata[tdata$sex=="Male",]
tdata = tdata[!is.na(tdata[,1]),]
tdata$lifespan = tdata$deathyear-tdata$birthyear
par(mfrow=c(1,1))

#inb
df = data.frame(x = tdata$inb,y = tdata$lifespan, year=tdata$birthyear)
df = df[complete.cases(df),]
df = df[df$x<0.5,]
x = df$x #tdata$inb
y = df$y 
year = df$year
modelout = glmer(y~x+(1|year), family="poisson")
summary(modelout)
sem.model.fits(modelout)

pdf("~/Desktop/F_life_males.pdf", width=5, height=5, useDingbats=FALSE, onefile=TRUE)
plot(-100,-100, xlab="inbreeding coefficient", ylab="lifespan", ylim=c(0.75, 6), xlim=c(0,0.4), main="males")
plotLMER.fnc(modelout, linecolor=colors4[1], lwd=4, pred="x", n=100, addToExistingPlot=TRUE, fun=exp)
points(x=x, y=jitter(y, 0.25), pch=21, col=alpha("grey30", 0.5), bg=alpha("grey30", 0.5), cex=0.75)
dev.off()

#how much shorter
#cutoff = mean(df$x) *2
#mean(df$y[df$x>cutoff])-mean(df$y)

####F to lifespan - females only####
tdata = kdata[kdata$genback>cutoff,]
tdata = tdata[tdata$birthyear<=2005,]
tdata = tdata[tdata$sex=="Female",]
tdata = tdata[!is.na(tdata[,1]),]
tdata$lifespan = tdata$deathyear-tdata$birthyear
par(mfrow=c(1,1))

#inb
df = data.frame(x = tdata$inb,y = tdata$lifespan, year=tdata$birthyear)
df = df[complete.cases(df),]
df = df[df$x<0.5,]
x = df$x #tdata$inb
y = df$y 
year = df$year
modelout = glmer(y~x+(1|year), family="poisson")
summary(modelout)
sem.model.fits(modelout)

pdf("~/Desktop/F_life_females.pdf", width=5, height=5, useDingbats=FALSE, onefile=TRUE)
plot(-100,-100, xlab="inbreeding coefficient", ylab="lifespan", ylim=c(0.75, 6), xlim=c(0,0.4), main="females")
plotLMER.fnc(modelout, linecolor=colors4[1], lwd=4, pred="x", n=100, addToExistingPlot=TRUE, fun=exp)
points(x=x, y=jitter(y, 0.25), pch=21, col=alpha("grey30", 0.5), bg=alpha("grey30", 0.5), cex=0.75)
dev.off()

#how much shorter
#cutoff = mean(df$x) *2
#mean(df$y[df$x>cutoff])-mean(df$y)

####sibling pairs offspring####
tdata = data[data$birthyear<=2005,]
tdata = tdata[!is.na(tdata[,1]),]
tdata$lifespan = tdata$deathyear-tdata$birthyear
tdata = tdata[tdata$lifespan>0,]
par(mfrow=c(1,1))

#sibs
parentpairs = NULL
for(r in 1:nrow(tdata)){
  parentpairs = c(parentpairs, paste(tdata$momid[r], tdata$dadid[r], sep="_"))
}
tdata$parentpairs = parentpairs
parentpairs = unique(parentpairs)
ppairs = data.frame(parentpairs=parentpairs, upairs = seq(1, length(parentpairs), 1))
tdata = merge(x=tdata, y=ppairs, by="parentpairs")
temp = table(tdata$upairs)
temp = temp[temp>1]
temp = names(temp)
temp = as.numeric(temp)
tpdata = tdata[tdata$upairs %in% temp,]

df = data.frame(inb = rep(NA, 1000), pkinsdiff = rep(NA, 1000), offdiff = rep(NA, 1000), year = rep(NA, 1000), lifespan = rep(NA, 1000))
count = 0
for(r in 1:length(temp)){
  t = tpdata[tpdata$upairs==temp[r],]
  t = t[!is.na(t$pkins),,drop=FALSE]
  if(nrow(t)<2){next}
  #1-2
  count = count + 1
  df$inb[count]       = t$inb[1]
  df$year[count]      = t$year[1]
  if(t$pkins[1] > t$pkins[2]){
    df$pkinsdiff[count] = t$pkins[1] - t$pkins[2]
    df$offdiff[count]   = t$off[1]/t$off[2]
    df$lifespan[count]  = t$lifespan[1] - t$lifespan[2]
  }else if(t$pkins[1] < t$pkins[2]){
    df$pkinsdiff[count] = t$pkins[2] - t$pkins[1]
    df$offdiff[count]   = t$off[2]/t$off[1]
    df$lifespan[count]  = t$lifespan[2] - t$lifespan[1]
  }
  
  #2-3
  if(nrow(t)==3){
    count = count + 1
    df$inb[count]       = t$inb[1]
    df$year[count]      = t$year[1]
    if(t$pkins[2] > t$pkins[3]){
      df$pkinsdiff[count] = t$pkins[2] - t$pkins[3]
      df$offdiff[count]   = t$off[2]/t$off[3]
      df$lifespan[count]  = t$lifespan[2] - t$lifespan[3]
    }else if(t$pkins[2] < t$pkins[3]){
      df$pkinsdiff[count] = t$pkins[3] - t$pkins[2]
      df$offdiff[count]   = t$off[3]/t$off[2]
      df$lifespan[count]  = t$lifespan[3] - t$lifespan[2]
    }
  }
  
  #1-3
  if(nrow(t)==3){
    count = count + 1
    df$inb[count]       = t$inb[1]
    df$year[count]      = t$year[1]
    if(t$pkins[1] > t$pkins[3]){
      df$pkinsdiff[count] = t$pkins[1] - t$pkins[3]
      df$offdiff[count]   = t$off[1]/t$off[3]
      df$lifespan[count]  = t$lifespan[1] - t$lifespan[3]
    }else if(t$pkins[1] < t$pkins[3]){
      df$pkinsdiff[count] = t$pkins[3] - t$pkins[1]
      df$offdiff[count]   = t$off[3]/t$off[1]
      df$lifespan[count]  = t$lifespan[3] - t$lifespan[1]
    }
  }
}

df = df[1:count,]
df = df[complete.cases(df),]
df = df[df$offdiff<=2,]
x = df$pkinsdiff
y = df$offdiff
f = df$inb
year = df$year
modelout = lme(y~x,random=~1|year)
summary(modelout)
sem.model.fits(modelout)

pdf("~/Desktop/sibs_off.pdf", width=5, height=5, useDingbats=FALSE, onefile=TRUE)
plot(-100,-100, xlab="diff. mate kinship", ylab="relative reproductive success", ylim=c(-0.5, 2.5), xlim=c(0,0.3), main="off")
abline(h=1, lwd=1, col="grey50", lty=2)
segments(x0=min(x), x1=max(x), y0=summary(modelout)$coefficients$fixed[1], y1=(summary(modelout)$coefficients$fixed[2]*max(x))+summary(modelout)$coefficients$fixed[1], col=colors4[2], lwd=3)
points(x=x, y=y, pch=21, col=alpha("grey30", 0.5), bg=alpha("grey30", 0.5), cex=0.75)
dev.off()


####sibling pairs - offspring survive####
tdata = data[data$birthyear<=2005,]
tdata = tdata[!is.na(tdata[,1]),]
tdata$lifespan = tdata$deathyear-tdata$birthyear
tdata = tdata[tdata$lifespan>0,]
par(mfrow=c(1,1))

#sibs
parentpairs = NULL
for(r in 1:nrow(tdata)){
  parentpairs = c(parentpairs, paste(tdata$momid[r], tdata$dadid[r], sep="_"))
}
tdata$parentpairs = parentpairs
parentpairs = unique(parentpairs)
ppairs = data.frame(parentpairs=parentpairs, upairs = seq(1, length(parentpairs), 1))
tdata = merge(x=tdata, y=ppairs, by="parentpairs")
temp = table(tdata$upairs)
temp = temp[temp>1]
temp = names(temp)
temp = as.numeric(temp)
tpdata = tdata[tdata$upairs %in% temp,]

df = data.frame(inb = rep(NA, 1000), pkinsdiff = rep(NA, 1000), offdiff = rep(NA, 1000), year = rep(NA, 1000), lifespan = rep(NA, 1000))
count = 0
for(r in 1:length(temp)){
  t = tpdata[tpdata$upairs==temp[r],]
  t = t[!is.na(t$pkins),,drop=FALSE]
  if(nrow(t)<2){next}
  #1-2
  count = count + 1
  df$inb[count]       = t$inb[1]
  df$year[count]      = t$year[1]
  if(t$pkins[1] > t$pkins[2]){
    df$pkinsdiff[count] = t$pkins[1] - t$pkins[2]
    df$offdiff[count]   = t$off_survive[1]/t$off_survive[2]
    df$lifespan[count]  = t$lifespan[1] - t$lifespan[2]
  }else if(t$pkins[1] < t$pkins[2]){
    df$pkinsdiff[count] = t$pkins[2] - t$pkins[1]
    df$offdiff[count]   = t$off_survive[2]/t$off_survive[1]
    df$lifespan[count]  = t$lifespan[2] - t$lifespan[1]
  }
  
  #2-3
  if(nrow(t)==3){
    count = count + 1
    df$inb[count]       = t$inb[1]
    df$year[count]      = t$year[1]
    if(t$pkins[2] > t$pkins[3]){
      df$pkinsdiff[count] = t$pkins[2] - t$pkins[3]
      df$offdiff[count]   = t$off_survive[2]/t$off_survive[3]
      df$lifespan[count]  = t$lifespan[2] - t$lifespan[3]
    }else if(t$pkins[2] < t$pkins[3]){
      df$pkinsdiff[count] = t$pkins[3] - t$pkins[2]
      df$offdiff[count]   = t$off_survive[3]/t$off_survive[2]
      df$lifespan[count]  = t$lifespan[3] - t$lifespan[2]
    }
  }
  
  #1-3
  if(nrow(t)==3){
    count = count + 1
    df$inb[count]       = t$inb[1]
    df$year[count]      = t$year[1]
    if(t$pkins[1] > t$pkins[3]){
      df$pkinsdiff[count] = t$pkins[1] - t$pkins[3]
      df$offdiff[count]   = t$off_survive[1]/t$off_survive[3]
      df$lifespan[count]  = t$lifespan[1] - t$lifespan[3]
    }else if(t$pkins[1] < t$pkins[3]){
      df$pkinsdiff[count] = t$pkins[3] - t$pkins[1]
      df$offdiff[count]   = t$off_survive[3]/t$off_survive[1]
      df$lifespan[count]  = t$lifespan[3] - t$lifespan[1]
    }
  }
}

df = df[1:count,]
df = df[complete.cases(df),]
df = df[df$offdiff<=2,]
x = df$pkinsdiff
y = df$offdiff
f = df$inb
year = df$year
modelout = lme(y~x,random=~1|year)
summary(modelout)
sem.model.fits(modelout)

pdf("~/Desktop/sibs_offsurv.pdf", width=5, height=5, useDingbats=FALSE, onefile=TRUE)
plot(-100,-100, xlab="diff. mate kinship", ylab="relative reproductive success", ylim=c(-0.5, 2.5), xlim=c(0,0.3), main="off_survive")
abline(h=1, lwd=1, col="grey50", lty=2)
segments(x0=min(x), x1=max(x), y0=summary(modelout)$coefficients$fixed[1], y1=(summary(modelout)$coefficients$fixed[2]*max(x))+summary(modelout)$coefficients$fixed[1], col=colors4[2], lwd=3)
points(x=x, y=y, pch=21, col=alpha("grey30", 0.5), bg=alpha("grey30", 0.5), cex=0.75)
dev.off()

####lethal equivalents####
tdata = data[data$genback>cutoff,]
tdata = tdata[tdata$birthyear<=2006,]
tdata = tdata[!is.na(tdata[,1]),]
tdata$lifespan = tdata$deathyear-tdata$birthyear

par(mfrow=c(1,1))
years = seq(1989,2005,1)
LEoff = data.frame(year=years, slope=rep(NA, length(years)), p=rep(NA, length(years)))
LEoff_survive = data.frame(year=years, slope=rep(NA, length(years)), p=rep(NA, length(years)))
for(y in 1:length(years)){
  temp = tdata[tdata$birthyear<=years[y] & tdata$deathyear>=years[y],,drop=FALSE]
  if(nrow(temp)>10){
    #test for all 0s
    x = temp[temp$inb>0,,drop=FALSE]
    if(nrow(x)>0){
      #all offspring
      plot(temp$inb,log(temp$off), xlim=c(0,0.5), main=paste(years[y], "offspring"))
      LEoff$slope[y] = summary(lm(log(temp$off+1)~temp$inb))$coefficients[2,1]
      LEoff$slope[y] = summary(lm(log(temp$off+1)~temp$inb))$coefficients[2,4]
      
      #survivng offspring
      plot(temp$inb,log(temp$off_survive), xlim=c(0,0.5), main=paste(years[y], "offspring surviving"))
      LEoff_survive$slope[y] = summary(lm(log(temp$off_survive+1)~temp$inb))$coefficients[2,1]
      LEoff_survive$slope[y] = summary(lm(log(temp$off_survive+1)~temp$inb))$coefficients[2,4]
    }
  }
}

mean(LEoff$slope, na.rm=TRUE)
sd(LEoff$slope, na.rm=TRUE)/sqrt(length(LEoff$slope[!is.na(LEoff$slope)]))
mean(LEoff_survive$slope, na.rm=TRUE)
sd(LEoff_survive$slope, na.rm=TRUE)/sqrt(length(LEoff_survive$slope[!is.na(LEoff_survive$slope)]))
