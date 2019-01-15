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

kins = read.table("outputbestped0/kinships_offspring.csv", sep=",", header=TRUE)
kins = kins[kins$id %in% kdata$id,]

df = data.frame(x = kins$matekinship,y = kins$numberoffspring, year=kins$year, id = kins$id)
df = df[complete.cases(df),]
x = df$x
y = df$y
year = df$year
id = df$id
modelout = glmer(y~x+(1|year)+(1|id), family="poisson")
summary(modelout)
sem.model.fits(modelout)
