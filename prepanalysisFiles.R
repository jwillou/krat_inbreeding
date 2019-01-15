setwd("/Volumes/jwillou/masterbayes/analyze_WH/")
load(file = "KRATnewdisp.Rdata")

library(coda)
library(MasterBayes)
library(kinship2)
library(pedigree)
library(pedantics)

####background and setup####
phenotypes  = read.csv("input/KRATP.csv", header=TRUE, strip.white=TRUE, as.is=TRUE)
ids = indvs = unique(phenotypes$id)
genotypes   = read.csv("input/KRATG.csv", header=TRUE, strip.white=TRUE, as.is=TRUE)
bdyear      = read.table("outputbestped0/bdyear.csv", header=TRUE, sep=",")
remove(model2, modelFULLnewprelim, modelYRL)
allpeds = modelFULLnewfinal$P

OUTPUT = data.frame(id=as.numeric(as.character(unique(c(ids, allpeds[2:ncol(allpeds)])))))
nopheno = data = NULL
for(i in 1:nrow(OUTPUT)){
  if(as.character(OUTPUT$id[i])==0 | is.na(as.character(OUTPUT$id[i]))){
    next
  }else{
    temp = phenotypes[phenotypes$id==as.character(OUTPUT$id[i]),,drop=FALSE]
    temp$id = as.numeric(temp$id)
    if(nrow(temp)==0){
      temp = data.frame(row=0, id=as.character(OUTPUT$id[i]), offspring=NA, sex=NA, status=NA, terr=NA, lat=NA, long=NA, disp=NA, Lddisp=NA, year=NA)
      nopheno = rbind(nopheno, temp)
      next
    }
    data = rbind(data, temp[temp$year==min(temp$year),])
  }
}

bestped = modeP(modelFULLnewfinal$P, threshold=0, marginal=TRUE, USasNA=TRUE)
pedout = cbind(bestped$P, bestped$prob, bestped$prob.male)
colnames(pedout) = c("id", "momid", "dadid", "prob", "prob.male")
hist(bestped$prob, ylim=c(0,800), breaks=seq(0,1,0.01), col="red", main=NA, xlab="probability of parental assignment", ylab="frequency")
#median(bestped$prob)
#mt = table(bestped$prob)
#names(mt[mt==max(mt)])

#set up output
ALLINDS  = data.frame(id = data$id)
ALLINDS  = merge(ALLINDS, pedout, by="id", all.x=TRUE)

#prepare pedigree
tallpeds = as.data.frame(bestped$P)
colnames(tallpeds) = c("id", "momid", "dadid")
tallpeds$id = as.numeric(as.character(tallpeds$id))
tallpeds$momid = as.numeric(as.character(tallpeds$momid))
tallpeds$dadid = as.numeric(as.character(tallpeds$dadid))

#merge with pheno data, reset columns to be numeric
tallpeds = merge(x=tallpeds, y=data, by="id",all.x=TRUE, all.y=TRUE)
tallpeds$id = as.numeric(as.character(tallpeds$id))
tallpeds$momid = as.numeric(as.character(tallpeds$momid))
tallpeds$dadid = as.numeric(as.character(tallpeds$dadid))  

#estimate pedigree variables  
toped  = data.frame(id=tallpeds$id, momid=tallpeds$momid, dadid=tallpeds$dadid)

#check that either both dad and mom are missing or that both are present
for(r in 1:nrow(toped)){
  if(is.na(toped$dadid[r]) | is.na(toped$momid[r]) | toped$dadid[r]==0 | toped$momid[r]==0){
    toped[r, 2:3] = NA
  }
}
pedord = pedigree::orderPed(toped)
toped  = toped[order(pedord),]
#pedem  = pedigreemm::pedigree(dam=toped[,2], sire=toped[,3], label=toped[,1] )  

#set up pedigree for kinship2 package
toped2  = data.frame(id=tallpeds$id, momid=tallpeds$momid, dadid=tallpeds$dadid, sex=tallpeds$sex)
rownames(toped2) = seq(1,nrow(toped2),1)
unassign = as.numeric(as.character(unique(toped2$id[is.na(toped2$sex)])))
for(u in 1:nrow(toped2)){
  if(!is.na(match(toped2$momid[u], unassign))){
    toped2$momid[u] = NA
    next
  }
  if(!is.na(match(toped2$dadid[u], unassign))){
    toped2$dadid[u] = NA
  }
}
toped2 = toped2[!is.na(toped2$sex),]
levels(toped2$sex) = c(levels(toped2$sex), "female", "male")
toped2[toped2$sex=="Female", 4] = "female"
toped2[toped2$sex=="Male",   4] = "male"  

#check that either both dad and mom are missing or that both are present
for(r in 1:nrow(toped2)){
  if(is.na(toped2$dadid[r]) | is.na(toped2$momid[r]) | toped2$dadid[r]==0 | toped2$momid[r]==0){
    toped2[r, 2:3] = "us"
  }
}
pedord = pedigree::orderPed(toped2)
toped2  = toped2[order(pedord),]
ped = kinship2::pedigree(id=toped2$id, dadid=toped2$dadid, momid=toped2$momid, sex=toped2$sex, missid="us")  
#write.table(ped, "outputbestped0/allkinships.csv", sep=",", row.names=FALSE, col.names=TRUE)
#pdf("~/Desktop/plot_p2000.pdf", width=10, height=10, useDingbats=FALSE, onefile=TRUE)
#setup pedigree
topedtemp = merge(x=toped2, y=data, by="id", all.x=TRUE, all.y=FALSE)
topedtemp$momid[topedtemp$momid=="us"] = NA
topedtemp$dadid[topedtemp$dadid=="us"] = NA
levels(topedtemp$sex.y) = c(levels(topedtemp$sex.y), 0, 1)
topedtemp$sex.y[topedtemp$sex.y=="Female"] = 0
topedtemp$sex.y[topedtemp$sex.y=="Male"] = 1
topedtemp$sex.y = as.numeric(as.character(topedtemp$sex.y))
pedord    = pedigree::orderPed(topedtemp)
topedtemp = topedtemp[order(pedord),]
topedped  = data.frame(id=topedtemp$id, dam=topedtemp$momid, sire=topedtemp$dadid)

#calculate ped stats and save F
st = pedigreeStats(Ped=topedped, cohorts=topedtemp$year, graphicalReport=FALSE)
newF = cbind(st$analyzedPedigree, st$inbreedingCoefficients)
newF$dam = newF$sire = NULL
colnames(newF) = c("id", "newF")
ALLINDS = merge(x=ALLINDS, y=newF, by="id", all.x=TRUE)

#write summary
#write.table(t(names(pedStatSummary(st))), "pedstatssum_pedantics.csv", sep=",", col.names=FALSE, row.names=FALSE, append=TRUE)
#write.table(t(pedStatSummary(st)), "pedstatssum_pedantics.csv", sep=",", col.names=FALSE, row.names=FALSE, append=TRUE)

#drawPedigree(Ped=topedped, cohorts=topedtemp$year, sex=topedtemp$sex.y)
#dev.off()

#analyze pedigree
#estimate inbreeding 
inbped = data.frame(id=toped2$id, dadid=toped2$dadid, momid=toped2$momid)
inb = calcInbreeding(inbped)
inb = data.frame(id = toped2$id, inb = inb)
colnames(inb) = c("id", "inb")
ALLINDS = merge(x=ALLINDS, y=inb, by="id", all.x=TRUE)

#count number of offspring
off = data.frame(id = c(names(table(inbped$momid)), names(table(inbped$dadid))), off = c(as.numeric(table(inbped$momid)), as.numeric(table(inbped$dadid))))
colnames(off) = c("id","off")
ALLINDS = merge(x=ALLINDS, y=off, by="id", all.x=TRUE)

#set up pedigree for kinship2 package withonly survivors to maturity
toped3  = data.frame(id=tallpeds$id, momid=tallpeds$momid, dadid=tallpeds$dadid, sex=tallpeds$sex)
bdyear$deathage = bdyear$deathyear - bdyear$birthyear
tokeep = unique(bdyear$id[bdyear$deathage>=1])
toped3 = toped3[toped3$id %in% tokeep,]
rownames(toped3) = seq(1,nrow(toped3),1)
unassign = as.numeric(as.character(unique(toped3$id[is.na(toped3$sex)])))
for(u in 1:nrow(toped3)){
  if(!is.na(match(toped3$momid[u], unassign))){
    toped3$momid[u] = NA
    next
  }
  if(!is.na(match(toped3$dadid[u], unassign))){
    toped3$dadid[u] = NA
  }
}
toped3 = toped3[!is.na(toped3$sex),]
levels(toped3$sex) = c(levels(toped3$sex), "female", "male")
toped3[toped3$sex=="Female", 4] = "female"
toped3[toped3$sex=="Male",   4] = "male"  

#check that either both dad and mom are missing or that both are present
for(r in 1:nrow(toped3)){
  if(is.na(toped3$dadid[r]) | is.na(toped3$momid[r]) | toped3$dadid[r]==0 | toped3$momid[r]==0){
    toped3[r, 2:3] = "us"
  }
}
pedord = pedigree::orderPed(toped3)
toped3  = toped3[order(pedord),]

#estimate inbreeding 
inbped = data.frame(id=toped3$id, dadid=toped3$dadid, momid=toped3$momid)
off = data.frame(id = c(names(table(inbped$momid)), names(table(inbped$dadid))), off = c(as.numeric(table(inbped$momid)), as.numeric(table(inbped$dadid))))
colnames(off) = c("id","off_survive")
ALLINDS = merge(x=ALLINDS, y=off, by="id", all.x=TRUE)

#estimate kinship between pairs, average for individuals
xx=unique(as.numeric(sapply(toped2$momid, as.character)))
momsactual = xx[!is.na(xx)]
xy=unique(as.numeric(sapply(toped2$dadid, as.character)))
dadsactual = xy[!is.na(xy)]

#calc all pairwise kinship pairs
kins = kinship(ped)
kinshold = kins

#write out kinships to ref later if needed
writekins = cbind(rownames(kins), kins)
colnames(writekins) = c("id", colnames(kins))
#write.table(writekins, "allkins.csv", col.names=TRUE, row.names=FALSE, sep=",")

#find kinship estimates for parental pairs
pkins = data.frame(id = c(momsactual, dadsactual), pkins = rep(NA, length(c(momsactual, dadsactual))), umates = rep(NA, length(c(momsactual, dadsactual))))

#moms
for(m in 1:length(momsactual)){
  ptemp = toped2[toped2$momid==momsactual[m],,drop=FALSE]
  ptemp = ptemp[complete.cases(ptemp),,drop=FALSE]
  #find kinship estimates
  if(nrow(ptemp)>0){
    x=kins[rownames(kins)==ptemp$momid[1],]
    x=x[names(x) %in% ptemp$dadid]
    pkins$pkins[pkins$id==momsactual[m]] = mean(x, na.rm=TRUE)
    pkins$umates[pkins$id==momsactual[m]] = length(unique(ptemp$dadid))
  }
}

#dads
for(d in 1:length(dadsactual)){
  ptemp = toped2[toped2$dadid==dadsactual[d],,drop=FALSE]
  ptemp = ptemp[complete.cases(ptemp),,drop=FALSE]
  #find kinship estimates
  if(nrow(ptemp)>0){
    x=kins[rownames(kins)==ptemp$dadid[1],]
    x=x[names(x) %in% ptemp$momid]
    pkins$pkins[pkins$id==dadsactual[d]] = mean(x, na.rm=TRUE)
    pkins$umates[pkins$id==dadsactual[d]] = length(unique(ptemp$momid))
  }
}
colnames(pkins) = c("id", "pkinships", "nmates")

#kinships
kins   = data.frame(id=pkins$id, pkins=pkins$pkinships)
ALLINDS = merge(x=ALLINDS, y=kins, by="id", all.x=TRUE)


#number of mates
nmates = data.frame(id=pkins$id, nmates=pkins$nmates)
ALLINDS = merge(x=ALLINDS, y=nmates, by="id", all.x=TRUE)

tdata = phenotypes
newdata = data.frame(id=unique(tdata$id), birthyear = rep(0, length(unique(tdata$id))), deathyear = rep(0, length(unique(tdata$id))), maxyear = rep(0, length(unique(tdata$id))), minyear = rep(0, length(unique(tdata$id))))
for(i in 1:length(tdata$id)){
  temp = tdata[tdata$id==as.character(tdata$id[i]),,drop=FALSE]
  if(temp$offspring[temp$year==min(temp$year,na.rm=TRUE)]==0){
    newdata$birthyear[newdata$id==as.character(tdata$id[i])] = min(temp$year,na.rm=TRUE) - 1
  }else{
    newdata$birthyear[newdata$id==as.character(tdata$id[i])] = min(temp$year,na.rm=TRUE)
  }
  if((temp$status[temp$year==max(temp$year,na.rm=TRUE)]=="P" | is.na(temp$status[temp$year==max(temp$year,na.rm=TRUE)])) & temp$offspring[temp$year==min(temp$year,na.rm=TRUE)]!=1){
    newdata$deathyear[newdata$id==as.character(tdata$id[i])] = max(temp$year,na.rm=TRUE) - 1
  }else{
    newdata$deathyear[newdata$id==as.character(tdata$id[i])] = max(temp$year,na.rm=TRUE)
  }
}
remove(i,tdata,temp)
ALLINDS = merge(x=ALLINDS, y=newdata, by="id", all.x=TRUE)

#calc offspring sex ratio
ALLINDS$nFemaleoff = rep(NA, nrow(ALLINDS))
ALLINDS$nMaleoff = rep(NA, nrow(ALLINDS))

allmoms = unique(ALLINDS$momid)
allmoms = allmoms[!is.na(allmoms)]
alldads = unique(ALLINDS$dadid)
alldads = alldads[!is.na(alldads)]

for(m in 1:length(allmoms)){
  temp = ALLINDS[ALLINDS$momid==as.character(allmoms[m]),,drop=FALSE]
  temp = temp[!is.na(temp$id),,drop=FALSE]
  ALLINDS$nFemaleoff[ALLINDS$id==as.character(allmoms[m])] = nrow(temp[temp$sex=="Female",,drop=FALSE])
  ALLINDS$nMaleoff[ALLINDS$id==as.character(allmoms[m])]   = nrow(temp[temp$sex=="Male",,drop=FALSE])
}

for(d in 1:length(alldads)){
  temp = ALLINDS[ALLINDS$dadid==as.character(alldads[d]),,drop=FALSE]
  temp = temp[!is.na(temp$id),,drop=FALSE]
  ALLINDS$nFemaleoff[ALLINDS$id==as.character(alldads[d])] = nrow(temp[temp$sex=="Female",,drop=FALSE])
  ALLINDS$nMaleoff[ALLINDS$id==as.character(alldads[d])]   = nrow(temp[temp$sex=="Male",,drop=FALSE])
}
ALLINDS = merge(x=ALLINDS, y=data, by="id", all.x=TRUE)

#indv kinships
ikins = apply(kinshold, 1, mean, na.rm=TRUE)
ikins = cbind(rownames(kinshold), ikins)
colnames(ikins) = c("id", "mkinship_all")
ALLINDS = merge(x=ALLINDS, y=ikins, by="id", all.x=TRUE)
#write.table(ALLINDS, "outputbestped0/summary_allstats2.csv", sep=",", col.names=TRUE, row.names=FALSE, append=FALSE)

####yearly kinships within pop and groups####
#kinshold
data = ALLINDS
years = seq(range(data$birthyear)[1], range(data$birthyear)[2], 1)
iout = data.frame(id=data$id)
for(y in 1:length(years)){
  temp = data[data$birthyear<=years[y] & data$deathyear>=years[y],]
  tkins = kinshold[rownames(kinshold) %in% temp$id, colnames(kinshold) %in% temp$id]
  tkins = apply(tkins, 1, mean, na.rm=TRUE)
  tkins = data.frame(id=names(tkins), k=as.numeric(as.character(tkins)))
  iout = merge(x=iout, y=tkins, by="id", all.x=TRUE)
  colnames(iout) = c(colnames(iout)[1:(ncol(iout)-1)], paste("K", years[y], sep=""))
}
#write.table(iout, "totalpopkinship_yearly.csv", row.names=FALSE, col.names=TRUE, sep=",")

#groups
plotlocs = data
plotlocs$adjlat  = plotlocs$lat - min(plotlocs$lat, na.rm=TRUE)
plotlocs$adjlong = plotlocs$long - min(plotlocs$long, na.rm=TRUE)
data = plotlocs
years = seq(range(data$birthyear)[1], range(data$birthyear)[2], 1)
g1kins = g2kins = g3kins = g4kins = data.frame(id=data$id)
for(y in 1:length(years)){
  temp = data[data$birthyear<=years[y] & data$deathyear>=years[y],]
  
  g1     = temp[temp$adjlong<425,,drop=FALSE]
  tkins  = kinshold[rownames(kinshold) %in% g1$id, colnames(kinshold) %in% g1$id,drop=FALSE]
  tkins  = apply(tkins, 1, mean, na.rm=TRUE)
  tkins  = data.frame(id=names(tkins), k=as.numeric(as.character(tkins)))
  if(nrow(tkins)>0){
    g1kins = merge(x=g1kins, y=tkins, by="id", all.x=TRUE) 
  }else{
    tkins = data.frame(id = g1$id, K=rep(NA, nrow(g1)))
    g1kins = merge(x=g1kins, y=tkins, by="id", all.x=TRUE) 
  }
  colnames(g1kins) = c(colnames(g1kins)[1:(ncol(g1kins)-1)], paste("K", years[y], sep=""))
  
  g2  = temp[temp$adjlong<635 & temp$adjlong>425,,drop=FALSE]
  tkins  = kinshold[rownames(kinshold) %in% g2$id, colnames(kinshold) %in% g2$id,drop=FALSE]
  tkins  = apply(tkins, 1, mean, na.rm=TRUE)
  tkins  = data.frame(id=names(tkins), k=as.numeric(as.character(tkins)))
  if(nrow(tkins)>0){
    g2kins = merge(x=g2kins, y=tkins, by="id", all.x=TRUE) 
  }else{
    tkins = data.frame(id = g2$id, K=rep(NA, nrow(g2)))
    g2kins = merge(x=g2kins, y=tkins, by="id", all.x=TRUE)  
  }
  colnames(g2kins) = c(colnames(g2kins)[1:(ncol(g2kins)-1)], paste("K", years[y], sep=""))
  
  g3  = temp[temp$adjlong>635 & temp$adjlat<700,,drop=FALSE]
  tkins  = kinshold[rownames(kinshold) %in% g3$id, colnames(kinshold) %in% g3$id,drop=FALSE]
  tkins  = apply(tkins, 1, mean, na.rm=TRUE)
  tkins  = data.frame(id=names(tkins), k=as.numeric(as.character(tkins)))
  if(nrow(tkins)>0){
    g3kins = merge(x=g3kins, y=tkins, by="id", all.x=TRUE) 
  }else{
    tkins = data.frame(id = g3$id, K=rep(NA, nrow(g3)))
    g3kins = merge(x=g3kins, y=tkins, by="id", all.x=TRUE)  
  }
  colnames(g3kins) = c(colnames(g3kins)[1:(ncol(g3kins)-1)], paste("K", years[y], sep=""))
  
  g4  = temp[temp$adjlat>700,,drop=FALSE]
  tkins  = kinshold[rownames(kinshold) %in% g4$id, colnames(kinshold) %in% g4$id,drop=FALSE]
  tkins  = apply(tkins, 1, mean, na.rm=TRUE)
  tkins  = data.frame(id=names(tkins), k=as.numeric(as.character(tkins)))
  if(nrow(tkins)>0){
    g4kins = merge(x=g4kins, y=tkins, by="id", all.x=TRUE) 
  }else{
    tkins = data.frame(id = g4$id, K=rep(NA, nrow(g4)))
    g4kins = merge(x=g4kins, y=tkins, by="id", all.x=TRUE)  
  }
  colnames(g4kins) = c(colnames(g4kins)[1:(ncol(g4kins)-1)], paste("K", years[y], sep=""))
}
g1kins = g1kins[rowSums(is.na(g1kins))!=(ncol(g1kins)-1), ]
g2kins = g2kins[rowSums(is.na(g2kins))!=(ncol(g2kins)-1), ]
g3kins = g3kins[rowSums(is.na(g3kins))!=(ncol(g3kins)-1), ]
g4kins = g4kins[rowSums(is.na(g4kins))!=(ncol(g4kins)-1), ]

#write.table(g1kins, "g1kinship.csv", row.names=FALSE, col.names=TRUE, sep=",")
#write.table(g2kins, "g2kinship.csv", row.names=FALSE, col.names=TRUE, sep=",")
#write.table(g3kins, "g3kinship.csv", row.names=FALSE, col.names=TRUE, sep=",")
#write.table(g4kins, "g4kinship.csv", row.names=FALSE, col.names=TRUE, sep=",")

####yearly summary####
data = ALLINDS
#iterate over years, isolating all pairings to determine F/noffspring for each breeding pair
years = seq(range(data$birthyear)[1], range(data$birthyear)[2], 1)
dataout = NULL
for(y in 1:length(years)){
  temp = data[data$birthyear==years[y],,drop=FALSE]
  moms  = unique(temp$momid)
  dads  = unique(temp$dadid)
  moms = moms[!is.na(moms)]
  if(length(moms[!is.na(moms)])>0){
    for(p in 1:length(moms)){
      pars = temp[temp$momid==as.numeric(as.character(moms[p])),,drop=FALSE]
      pars = pars[complete.cases(pars$id),,drop=FALSE]
      numberoff = nrow(pars)
      finb      = inb$inb[inb$id==as.character(moms[p])]
      writeout = c(years[y], as.character(moms[p]), finb, numberoff)
      dataout = rbind(dataout, writeout)
    }
  }
  if(length(dads[!is.na(dads)])>0){
    for(p in 1:length(dads)){
      pars = temp[temp$dadid==as.numeric(as.character(dads[p])),,drop=FALSE]
      pars = pars[complete.cases(pars$id),,drop=FALSE]
      numberoff = nrow(pars)
      finb      = inb$inb[inb$id==as.character(dads[p])]
      writeout = c(years[y], as.character(dads[p]), finb, numberoff)
      dataout = rbind(dataout, writeout)
    }
  }
}
rownames(dataout) = seq(1,nrow(dataout),1)
colnames(dataout) = c("year","id", "inb", "numberoffspring")
#write.table(dataout, "outputbestped0/inb_offspring.csv", row.names=FALSE, col.names=TRUE, sep=",")

#iterate over years, isolating all pairings to determine kinship/noffspring for each breeding pair
years = seq(range(data$birthyear)[1], range(data$birthyear)[2], 1)
dataout = NULL
for(y in 1:length(years)){
  temp = data$id[data$birthyear==years[y]]
  tpars = toped2[toped2$id %in% temp,,drop=FALSE]
  moms = as.numeric(names(table(tpars$momid)))
  moms = moms[!is.na(moms)]
  dads = as.numeric(names(table(tpars$dadid)))
  dads = dads[!is.na(dads)]
  if(length(moms)>0){
    for(m in 1:length(moms)){
      mtemp = tpars[tpars$momid==moms[m],,drop=FALSE]
      if(nrow(mtemp)>0){
        dads  = unique(mtemp$dadid)
        for(d in 1:length(dads)){
          mdkinship = kinshold[rownames(kinshold)==moms[m], colnames(kinshold)==as.numeric(dads[d])]
          if(mdkinship==0){
            next
          }
          numberoff = nrow(mtemp[mtemp$dadid==as.numeric(dads[d]),,drop=FALSE])
          writeout = c(years[y], moms[m], mdkinship, numberoff)
          dataout = rbind(dataout, writeout)
        }
      }
    }
  }
  if(length(dads)>0){
    for(m in 1:length(dads)){
      mtemp = tpars[tpars$dadid==dads[m],,drop=FALSE]
      if(nrow(mtemp)>0){
        moms  = unique(mtemp$momid)
        for(d in 1:length(moms)){
          mdkinship = kinshold[rownames(kinshold)==dads[m], colnames(kinshold)==as.numeric(moms[d])]
          if(mdkinship==0){
            next
          }
          numberoff = nrow(mtemp[mtemp$momid==as.numeric(moms[d]),,drop=FALSE])
          writeout = c(years[y], dads[m], mdkinship, numberoff)
          dataout = rbind(dataout, writeout)
        }
      }
    }
  }
}
rownames(dataout) = seq(1,nrow(dataout),1)
colnames(dataout) = c("year", "id", "matekinship", "numberoffspring")
write.table(dataout, "outputbestped0/kinships_offspring.csv", row.names=FALSE, col.names=TRUE, sep=",")

#iterate over years, estimate mean F and mean population kinship
years = seq(range(data$birthyear)[1], range(data$birthyear)[2], 1)
yearout = NULL
for(y in 1:length(years)){
 temp = data[data$birthyear<=years[y] & data$deathyear>=years[y],]
 tinb = inb[inb$id %in% temp$id,,drop=FALSE]
 tinb = tinb[tinb$inb>0,,drop=FALSE]
 tinb = mean(tinb$inb, na.rm=TRUE)
 tkins = kinshold[rownames(kinshold) %in% temp$id[temp$pkins>0], colnames(kinshold) %in% temp$id[temp$pkins>0]]
 tkins = apply(tkins, 1, mean, na.rm=TRUE)
 tkins = mean(tkins, na.rm=TRUE)
 tsex  = data[data$id %in% temp$id,]
 yearout = rbind(yearout, c(years[y], tinb, tkins, nrow(temp), nrow(tsex[tsex$sex=="Female",,drop=FALSE]), nrow(tsex[tsex$sex=="Male",,drop=FALSE])))
}
colnames(yearout) = c("year", "inbreeding", "kinship", "N", "females", "males")
#write.table(yearout, "~/Desktop/testout.csv", row.names=FALSE, col.names=TRUE, sep=",")

#define clusters and estimate linear distance between mounds
plotlocs = data
plotlocs$adjlat  = plotlocs$lat - min(plotlocs$lat, na.rm=TRUE)
plotlocs$adjlong = plotlocs$long - min(plotlocs$long, na.rm=TRUE)
data = plotlocs
years = sort(unique(c(plotlocs$birthyear, plotlocs$deathyear)))
colors4 = c("darkorchid3", "dodgerblue3", "goldenrod3", "firebrick3")
for(y in 1:length(years)){
  temp = plotlocs[plotlocs$birthyear<=years[y] & plotlocs$deathyear>=years[y],,drop=FALSE]
  temp  = temp[!is.na(temp$id),]
  if(nrow(temp)==0){next}
  png(file = paste("~/Desktop/maps/year", years[y], ".png", sep=""), bg = "white",width = 480, height = 480)
  plot(-100,-100,xlab="meters", ylab="meters", main=years[y], xlim=c(0,1200), ylim=c(0,1200))
  abline(h=425, lty=3, col="grey30")
  abline(h=635, lty=3, col="grey30")
  abline(v=700, lty=3, col="grey30")
  g1  = temp[temp$adjlong<425,,drop=FALSE]
  g2  = temp[temp$adjlong<635 & temp$adjlong>425,,drop=FALSE]
  g3  = temp[temp$adjlong>635 & temp$adjlat<700,,drop=FALSE]
  g4  = temp[temp$adjlat>700,,drop=FALSE]
  points(y=g1$adjlong, x=g1$adjlat, cex=0.6, col=colors4[1], pch=19)
  points(y=g2$adjlong, x=g2$adjlat, cex=0.6, col=colors4[2], pch=19)
  points(y=g3$adjlong, x=g3$adjlat, cex=0.6, col=colors4[3], pch=19)
  points(y=g4$adjlong, x=g4$adjlat, cex=0.6, col=colors4[4], pch=19)
  dev.off()
}
#compute average distance between adults for each year
euc.dist = function(x1, y1, x2, y2){sqrt(((x1 - x2) ^ 2) + ((y1 - y2) ^ 2))} 
output = NULL
for(y in 1:length(years)){
  temp  = plotlocs[plotlocs$birthyear<years[y] & plotlocs$deathyear>=years[y],,drop=FALSE] #adults
  temp2 = plotlocs[plotlocs$birthyear<=years[y] & plotlocs$deathyear>=years[y],,drop=FALSE] #all individuals
  temp  = temp[!is.na(temp$id),]
  temp2 = temp2[!is.na(temp2$id),]
  if(nrow(temp)==0){next}
  distanceAdults1 = distanceAdults2 = distanceAdults3 = distanceAdults4 = NULL
  distanceAll1    = distanceAll2    = distanceAll3    = distanceAll4    = NULL
  for(i in 1:(nrow(temp)-1)){
    for(i2 in (i+1):nrow(temp)){
      g1  = temp[temp$adjlong<425,,drop=FALSE]
      g12 = temp2[temp2$adjlong<425,,drop=FALSE]
      distanceAdults1 = c(distanceAdults1, euc.dist(g1$adjlat[i],   g1$adjlong[i],   g1$adjlat[i2],  g1$adjlong[i2]))
      distanceAll1    = c(distanceAll1,    euc.dist(g12$adjlat[i],  g12$adjlong[i],  g12$adjlat[i2], g12$adjlong[i2]))
      
      g2  = temp[temp$adjlong<635 & temp$adjlong>425,,drop=FALSE]
      g22 = temp2[temp2$adjlong<635 & temp2$adjlong>425,,drop=FALSE]
      distanceAdults2 = c(distanceAdults2, euc.dist(g2$adjlat[i],   g2$adjlong[i],   g2$adjlat[i2],  g2$adjlong[i2]))
      distanceAll2    = c(distanceAll2,    euc.dist(g22$adjlat[i],  g22$adjlong[i],  g22$adjlat[i2], g22$adjlong[i2]))
      
      g3  = temp[temp$adjlong>635 & temp$adjlat<700,,drop=FALSE]
      g32 = temp2[temp2$adjlong>635 & temp2$adjlat<700,,drop=FALSE]
      distanceAdults3 = c(distanceAdults3, euc.dist(g3$adjlat[i],   g3$adjlong[i],   g3$adjlat[i2],  g3$adjlong[i2]))
      distanceAll3    = c(distanceAll3,    euc.dist(g32$adjlat[i],  g32$adjlong[i],  g32$adjlat[i2], g32$adjlong[i2]))
      
      g4  = temp[temp$adjlat>700,,drop=FALSE]
      g42 = temp2[temp2$adjlat>700,,drop=FALSE]
      distanceAdults4 = c(distanceAdults4, euc.dist(g4$adjlat[i],   g4$adjlong[i],   g4$adjlat[i2],  g4$adjlong[i2]))
      distanceAll4    = c(distanceAll4,    euc.dist(g42$adjlat[i],  g42$adjlong[i],  g42$adjlat[i2], g42$adjlong[i2]))
    }
  }
  output = rbind(output, c(years[y], nrow(g1), nrow(g12), nrow(g2), nrow(g22), nrow(g3), nrow(g32), nrow(g4), nrow(g42), 
                           mean(distanceAdults1, na.rm=TRUE), mean(distanceAll1, na.rm=TRUE),
                           mean(distanceAdults2, na.rm=TRUE), mean(distanceAll2, na.rm=TRUE),
                           mean(distanceAdults3, na.rm=TRUE), mean(distanceAll3, na.rm=TRUE),
                           mean(distanceAdults4, na.rm=TRUE), mean(distanceAll4, na.rm=TRUE)))
}
colnames(output) = c("year", "Nadults1", "Nallindv1", "Nadults2", "Nallindv2", "Nadults3", "Nallindv3", "Nadults4", "Nallindv4", 
                     "avgdistAdults1", "avgdistAll1", "avgdistAdults2", "avgdistAll2", "avgdistAdults3", "avgdistAll3", "avgdistAdults4", "avgdistAll4")
yearout = merge(x=yearout, y=output, by="year", all.x=TRUE)
#write.table(output, "distance_density.csv", col.names=TRUE, row.names=FALSE, sep=",")

####mean F in each geographic group####
fdata = data
fdata = fdata[fdata$inb>0,]
years = sort(unique(c(fdata$birthyear, fdata$deathyear)))
g1Fadult = g2Fadult = g3Fadult = g4Fadult = NULL
g1F = g2F = g3F = g4F = NULL
output = NULL
for(y in 1:length(years)){
  temp  = fdata[fdata$birthyear<years[y] & fdata$deathyear>=years[y],,drop=FALSE] #adults
  temp2 = fdata[fdata$birthyear<=years[y] & fdata$deathyear>=years[y],,drop=FALSE] #all individuals
  temp  = temp[!is.na(temp$id),]
  temp2 = temp2[!is.na(temp2$id),]
  if(nrow(temp2)==0){
    g1Fadult = c(g1Fadult, NA)
    g1F      = c(g1F, NA)
    g2Fadult = c(g2Fadult, NA)
    g2F      = c(g2F, NA)
    g3Fadult = c(g3Fadult, NA)
    g3F      = c(g3F, NA)
    g4Fadult = c(g4Fadult, NA)
    g4F      = c(g4F, NA)
    next
  }else if(nrow(temp)==0){
    g1  = temp[temp$adjlong<425,,drop=FALSE]
    g1Fadult = c(g1Fadult, mean(g1$inb, na.rm=TRUE))
    g2  = temp[temp$adjlong<635 & temp$adjlong>425,,drop=FALSE]
    g2Fadult = c(g2Fadult, mean(g2$inb, na.rm=TRUE))
    g3  = temp[temp$adjlong>635 & temp$adjlat<700,,drop=FALSE]
    g3Fadult = c(g3Fadult, mean(g3$inb, na.rm=TRUE))
    g4  = temp[temp$adjlat>700,,drop=FALSE]
    g4Fadult = c(g4Fadult, mean(g4$inb, na.rm=TRUE))
    
    g1F      = c(g1F, NA)
    g2F      = c(g2F, NA)
    g3F      = c(g3F, NA)
    g4F      = c(g4F, NA)
    next

  }else{
    g1  = temp[temp$adjlong<425,,drop=FALSE]
    g12 = temp2[temp2$adjlong<425,,drop=FALSE]
    g1Fadult = c(g1Fadult, mean(g1$inb, na.rm=TRUE))
    g1F      = c(g1F, mean(g12$inb, na.rm=TRUE))
    
    g2  = temp[temp$adjlong<635 & temp$adjlong>425,,drop=FALSE]
    g22 = temp2[temp2$adjlong<635 & temp2$adjlong>425,,drop=FALSE]
    g2Fadult = c(g2Fadult, mean(g2$inb, na.rm=TRUE))
    g2F      = c(g2F, mean(g22$inb, na.rm=TRUE))
    
    g3  = temp[temp$adjlong>635 & temp$adjlat<700,,drop=FALSE]
    g32 = temp2[temp2$adjlong>635 & temp2$adjlat<700,,drop=FALSE]
    g3Fadult = c(g3Fadult, mean(g3$inb, na.rm=TRUE))
    g3F      = c(g3F, mean(g32$inb, na.rm=TRUE))
    
    g4  = temp[temp$adjlat>700,,drop=FALSE]
    g42 = temp2[temp2$adjlat>700,,drop=FALSE]
    g4Fadult = c(g4Fadult, mean(g4$inb, na.rm=TRUE))
    g4F      = c(g4F, mean(g42$inb, na.rm=TRUE))
  }
}
output = cbind(years, g1Fadult, g1F, g2Fadult, g2F, g3Fadult, g3F, g4Fadult, g4F)
colnames(output) = c("year", "g1Fadult", "g1F", "g2Fadult", "g2F", "g3Fadult", "g3F", "g4Fadult", "g4F")
yearout = merge(x=yearout, y=output, by="year", all.x=TRUE)
#write.table(output, "outputbestped0/F_groups.csv", col.names=TRUE, row.names=FALSE, sep=",")

####mean K and offspring in each geographic group####
fdata = data
fdata = fdata[fdata$pkins>0,]
years = sort(unique(c(fdata$birthyear, fdata$deathyear)))
g1Fadult = g2Fadult = g3Fadult = g4Fadult = NULL
g1F = g2F = g3F = g4F = NULL
output = NULL
for(y in 1:length(years)){
  temp  = fdata[fdata$birthyear<years[y] & fdata$deathyear>=years[y],,drop=FALSE] #adults
  temp2 = fdata[fdata$birthyear<=years[y] & fdata$deathyear>=years[y],,drop=FALSE] #all individuals
  temp  = temp[!is.na(temp$id),]
  temp2 = temp2[!is.na(temp2$id),]
  #if(nrow(temp)==0){next}
  
  g1  = temp[temp$adjlong<425,,drop=FALSE]
  if(nrow(g1)>1){
    g1kins = kinshold[rownames(kinshold) %in% g1$id, colnames(kinshold) %in% g1$id]
    g1kins = apply(g1kins, 1, mean, na.rn=TRUE)
    g1kins = mean(g1kins, na.rm=TRUE)
  }else if(nrow(g1)==1){
    g1kins = 0.5
  }else{
    g1kins = NA
  }
  g1Fadult = c(g1Fadult, g1kins)
  #
  g12 = temp2[temp2$adjlong<425,,drop=FALSE]
  if(nrow(g12)>1){
    g1kinsall = kinshold[rownames(kinshold) %in% g12$id, colnames(kinshold) %in% g12$id]
    g1kinsall = apply(g1kinsall, 1, mean, na.rn=TRUE)
    g1kinsall = mean(g1kinsall, na.rm=TRUE)
  }else if(nrow(g12)==1){
    g1kinsall = 0.5
  }else{
    g1kinsall = NA
  }
  g1F      = c(g1F,g1kinsall)
  
  g2  = temp[temp$adjlong<635 & temp$adjlong>425,,drop=FALSE]
  if(nrow(g2)>1){
    g2kins = kinshold[rownames(kinshold) %in% g2$id, colnames(kinshold) %in% g2$id]
    g2kins = apply(g2kins, 1, mean, na.rn=TRUE)
    g2kins = mean(g2kins, na.rm=TRUE)
  }else if(nrow(g2)==1){
    g2kins = 0.5
  }else{
    g2kins = NA
  }
  g2Fadult = c(g2Fadult, g2kins)
  #
  g22 = temp2[temp2$adjlong<635 & temp2$adjlong>425,,drop=FALSE]
  if(nrow(g22)>1){
    g2kinsall = kinshold[rownames(kinshold) %in% g22$id, colnames(kinshold) %in% g22$id]
    g2kinsall = apply(g2kinsall, 1, mean, na.rn=TRUE)
    g2kinsall = mean(g2kinsall, na.rm=TRUE)
  }else if(nrow(g22)==1){
    g2kinsall = 0.5
  }else{
    g2kinsall = NA
  }
  g2F = c(g2F,g2kinsall)
  
  g3  = temp[temp$adjlong>635 & temp$adjlat<700,,drop=FALSE]
  if(nrow(g3)>1){
    g3kins = kinshold[rownames(kinshold) %in% g3$id, colnames(kinshold) %in% g3$id]
    g3kins = apply(g3kins, 1, mean, na.rn=TRUE)
    g3kins = mean(g3kins, na.rm=TRUE)
  }else if(nrow(g3)==1){
    g3kins = 0.5
  }else{
    g3kins = NA
  }
  g3Fadult = c(g3Fadult, g3kins)
  #
  g32 = temp2[temp2$adjlong>635 & temp2$adjlat<700,,drop=FALSE]
  if(nrow(g32)>1){
    g3kinsall = kinshold[rownames(kinshold) %in% g32$id, colnames(kinshold) %in% g32$id]
    g3kinsall = apply(g3kinsall, 1, mean, na.rn=TRUE)
    g3kinsall = mean(g3kinsall, na.rm=TRUE)
  }else if(nrow(g32)==1){
    g3kinsall = 0.5
  }else{
    g3kinsall = NA
  }
  g3F = c(g3F,g3kinsall)
  
  g4  = temp[temp$adjlat>700,,drop=FALSE]
  if(nrow(g4)>1){
    g4kins = kinshold[rownames(kinshold) %in% g4$id, colnames(kinshold) %in% g4$id]
    g4kins = apply(g4kins, 1, mean, na.rn=TRUE)
    g4kins = mean(g4kins, na.rm=TRUE)
  }else if(nrow(g4)==1){
    g4kins = 0.5
  }else{
    g4kins = NA
  }
  g4Fadult = c(g4Fadult, g4kins)
  #
  g42 = temp2[temp2$adjlat>700,,drop=FALSE]
  if(nrow(g42)>1){
    g4kinsall = kinshold[rownames(kinshold) %in% g42$id, colnames(kinshold) %in% g42$id]
    g4kinsall = apply(g4kinsall, 1, mean, na.rn=TRUE)
    g4kinsall = mean(g4kinsall, na.rm=TRUE)
  }else if(nrow(g42)==1){
    g4kinsall = 0.5
  }else{
    g4kinsall = NA
  }
  g4F = c(g4F,g4kinsall)
  
}
output = cbind(years, g1Fadult, g1F, g2Fadult, g2F, g3Fadult, g3F, g4Fadult, g4F)
colnames(output) = c("year", "g1Kadult", "g1K", "g2Kadult", "g2K", "g3Kadult", "g3K", "g4Kadult", "g4K")
yearout = merge(x=yearout, y=output, by="year", all.x=TRUE)
write.table(yearout, "outputbestped0/yearly_summary.csv", col.names=TRUE, row.names=FALSE, sep=",")

