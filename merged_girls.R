
getwd()
setwd("U:/My Documents on U/Jay_MUAC_2017/muac_2017/merged_girls")
library(gamlss) # package for fitting GAMLSS (now depends on gamlss.dist and gamlss.data
library(gamlss.dist) # package for all gamlss.family distributions
library(gamlss.data) #  package for all example data used in GAMLSS
library(gamlss.add) # experimental package for new additive terms.
library(gamlss.nl) # package for fitting nonlinear models
library(lattice)
library(car)
library(foreign)
library(epiDisplay)
 

## Modified functions
source("mypercentilespred.R")
source("mycentiles.R")  
source("mycentilespred.R")

# Now we load the hes2hes3nhanes1 data, convert muac from mm to cm and age from years to months.
# We then keep only for children between 5 years and 25 years of age.

nhanes<-read.dta("hes2hes3nhanes1.dta")[,c(1,4,5)]
summ(nhanes)
nhanes$muac<-nhanes$muac/10 # convert to cm
nhanes$agem<-nhanes$agey*12 # convert age to months
## Keep data for ages >=5 to 25 years
# nhanes=nhanes[nhanes$agem>=60 & nhanes$agem<=300,]
summ(nhanes)

# Create seperate datasets for girls only
nhanes$sex<-factor(nhanes$sex,labels=c("Female","Male"))
newf2<-nhanes[nhanes$sex=="Female",]
summ(newf2)

# Now fit BCCG, BCPE and BCT GAMLSS regression models and compare model fit
# Once we select the best model, we use it to compute the Muac Z-scores
# Then we remove individuals considered to be outliers that present with Z scores 
# greater than or less than -4 SD, then merge this data with the WHO simulated data

f1 <- gamlss(muac~pb(agem),sigma.fo=~pb(agem), nu.fo=~pb(agem),tau.fo=~pb(agem),data=newf2, family=BCT)

# f2 <- lms(muac,agem ,data=newf2, method.pb="GAIC", k=log(length(newf2$muac)))
# f2$family
# f3 <- gamlss(muac~pb(agem),sigma.fo=~pb(agem),nu.fo=~pb(agem),tau.fo=~pb(agem), data=newf2, family=BCPE)
# f4 <- lms(muac, agem ,data=newf2, method.pb="GAIC", k=log(length(newf2$muac)), families="BCCG")

# Selection between models
## Residual deviances and AIC
# The GAIC function uses default penalty k = 2, giving the usual Akaike information criterion
# (AIC). Hence the usual AIC [equivalent to GAIC(k = 2)] selects model f1= BCT  as the best model
# (since it has the smallest value of AIC). 

# deviance(f1); deviance(f2);deviance(f3); deviance(f4)
# AIC(f1,f2,f3,f4)
# AIC(f1,f2,f3,f4, k=log(length(newf2$muac)))# Bayesian Information Critetion
# centiles(f4,newf2$agem, main="BCCG model",cent=c(3,15,50,85,97),lwd=3)

## Now we use Model f1 (BCT)  and calculate Z scores from both and eliminate subjects with Z scores +/- 4 SD
gyf <-newf2$muac; length(gyf)
gxf <-newf2$agem; length(gxf)
gzscores <- centiles.pred(f1, xname="agem",xvalues=gxf,yval=gyf, type="z-scores")

head(newf2)

newf2<-cbind(newf2,gzscores)
head(newf2)
max(newf2$gzscores)
min(newf2$gzscores)
hist(gzscores, xlab = "Z scores", main="Nhanes Girls Z-scores")

## Removing newf2 with Z scores  +/- 4SD 

nrow(newf2) # number of observations before dropping
newf<-newf2[gzscores <= 4.0 & gzscores >= -4.0,]
nrow(newf) # number of observations after dropping 

rm(newf2)

hist(newf$gzscores,xlim=c(min(newf$gzscores),max(newf$gzscores+2)), 
     main="Histogram Z scores after dropping 4 outliers")
max(newf$gzscores) # maximum is Z scores = 3
min(newf$gzscores) # minimum is Z scores = -4

# Prepare to merge this data with the WHO simulated data and re-fit the models
who<-read.dta("sim12000.dta") 
summ(who)
who <- who[,c(1, 2,3,7)]
summ(who)

who$agem<-who$age/30 # convert age from days to months
who$agey<-who$age/365 # convert age from days to years
summ(who)

# Create a girls dataset from the WHO
tab1(who$sex,graph=FALSE)
whog<-who[who$sex=="Female",]
summ(whog)

whog$gzscores<-whog$z
summ(whog) 
whog<-whog[,c(2,4,5,7)]
summ(whog)
 
# Merging the girls datasets}
summ(newf)

whog$sex<-1
summ(whog)
summ(newf)
newf <- newf[,c(3,1,4,5)]
summ(newf)
summ(whog)

# Now, we are ready to combine these two girls datasets 
merged.girls=rbind(whog,newf)
summ(merged.girls)

# Re-fit regression models to the merged data and assess goodness of fit and convergence

# h1 <- gamlss(muac~pb(agem),sigma.fo=~pb(agem), nu.fo=~pb(agem),tau.fo=~pb(agem),data=merged.girls, family=BCT)
# h2 <- lms(muac,agem ,data=merged.girls, method.pb="GAIC", k=log(length(merged.girls$muac)))
# h2$family
# h3 <- gamlss(muac~pb(agem),sigma.fo=~pb(agem),nu.fo=~pb(agem),tau.fo=~pb(agem), data=merged.girls, family=BCPE)

h4 <- lms(muac, agem ,data=merged.girls, method.pb="GAIC", k=log(length(merged.girls$muac)), families="BCCG")

# h5 <- gamlss(muac~pb(agem), sigma.formula=~pb(agem), family=BCCG, data=merged.girls) 


## Model selection
# deviance(h1); deviance(h2);deviance(h3); deviance(h4)
# AIC(h4,h5)
# AIC(h4,h5,k=log(length(merged.girls$muac)))


# We now generate Tables with LMS values using BCCG model (h2) and BCPE Model (h3) Using BCPEModel

summ(merged.girls$agem,graph=FALSE)
summ(merged.girls$muac,graph=FALSE)
gyf <-merged.girls$muac; length(gyf)
 
####################
gyf <-merged.girls$muac; length(gyf)
gxf <-merged.girls$agem; length(gxf)

zscores_BCCG1_girls<-mycent.pred(h4, xname="agem",xvalues=24:300,plot=TRUE,
                    legend=FALSE,type="standard-centiles",ylim=c(10,80),xlim=c(24,300),
                    dev=c(-4,-3,-2,-1,0,1,2,3,4),  xlab="Age (Months)", ylab="MUAC (cm)", 
                    main="Predicted Z-values for girls using BCPE model")
legend(315,50, c("-4","-3","-2","-1","0", "1","2","3","4"), fill=c(2:10), horiz=FALSE)

write.csv(zscores_BCCG1_girls, file = "imputed_12000_zscores_girls_BCCG1.csv")

##
################
## Percentiles
############

BCCG1_percentiles <- myper.pred(h4, xname="agem", plot=TRUE,xvalues=24:300,legend=FALSE,
                          xlab="Age (Months)",ylab="MUAC (cm)",cent=c(3,15,50,85,97),
                          main="Predicted girls percentiles using BCCG Model")
legend(315,35, c("3","15","50","85","97"), fill=(2:12), horiz=FALSE)
 

write.csv(BCCG1_percentiles, file = "imputed_12000_BCCG1_girls_percentiles.csv")
#############################
# Diagnostics
###########################
#########################
graphics.off()
plot(muac~agem, data=merged.girls, ylab="MUAC (cm)", xlab="Age (months)", 
     main="Girls BCCG regression line", cex=1.5, cex.lab=1.5, cex.axis = 1.5,cex.main=1.8)
lines(fitted(h4)~merged.girls$agem, col="red",lwd=3)

plot(h4) # BCCGo model
wp(h4, ylim.all = 4,xlim.all = 5,cex.lab=2.5,cex=1.5,cex.axis = 1.5,cex.main=1.8) # BCCG Model

fittedPlot(h4, x=merged.girls$agem) # BCCG

##################################
imputed3 <- read.csv("imputed_36000.csv", header = TRUE,as.is = TRUE)
imputed4 <- read.csv("Imputed_5to6y_only_12000.csv", header = TRUE,as.is = TRUE)
#############################
## imputed 3: 36000 data
##############################

summ(imputed3)
imputed3$agem<-imputed3$age/30 # convert age from days to months
imputed3$agey<-imputed3$age/365 # convert age from days to years
summ(imputed3)

# Create a girls dataset from the imputed1
tab1(imputed3$sex,graph=FALSE)
imputed3b<-imputed3[imputed3$sex=="Female",]
summ(imputed3b)

imputed3b$gzscores<-imputed3b$z
summ(imputed3b) 
imputed3b<-imputed3b[,c(2,4,5,7)]
summ(imputed3b)

# Merging the girls datasets}
summ(newf)

imputed3b$sex<-1
summ(imputed3b)
summ(newf)

# Now, we are ready to combine these two girls datasets 
rm(merged.girls)
merged.girls=rbind(imputed3b,newf)
summ(merged.girls)

# Re-fit regression models to the merged data and assess goodness of fit and convergence

# h1 <- gamlss(muac~pb(agem),sigma.fo=~pb(agem), nu.fo=~pb(agem),tau.fo=~pb(agem),data=merged.girls, family=BCT)
# h2 <- lms(muac,agem ,data=merged.girls, method.pb="GAIC", k=log(length(merged.girls$muac)))
# h2$family
# chose BCTo

# h3 <- gamlss(muac~pb(agem),sigma.fo=~pb(agem),nu.fo=~pb(agem),tau.fo=~pb(agem), data=merged.girls, family=BCPE)
h4 <- lms(muac, agem ,data=merged.girls, method.pb="GAIC",
          k=log(length(merged.girls$muac)), families="BCCG")

# We now generate Tables with LMS values using BCCG model (h2) and BCPE Model (h3) Using BCPEModel

summ(merged.girls$agem,graph=FALSE)
summ(merged.girls$muac,graph=FALSE)
gyf <-merged.girls$muac; length(gyf)

par(mar=c(5.1, 4.1, 4.1, 6.1), xpd=TRUE)

####
gyf <-merged.girls$muac; length(gyf)
gxf <-merged.girls$agem; length(gxf)

zscores_BCCG2<-mycent.pred(h4, xname="agem",xvalues=24:300,plot=TRUE,
                               legend=FALSE,type="standard-centiles",ylim=c(10,80),xlim=c(24,300),
                               dev=c(-4,-3,-2,-1,0,1,2,3,4),  xlab="Age (Months)", ylab="MUAC (cm)", 
                               main="Predicted Z-values for girls using BCCG model")
legend(315,70, c("-4","-3","-2","-1","0", "1","2","3","4"), fill=c(2:10), horiz=FALSE)

write.csv(zscores_BCCG2, file = "imputed_36000_zscores_girls_BCCG2.csv")

#############################
BCCG2percentiles <- myper.pred(h4, xname="agem", plot=TRUE,xvalues=24:300,legend=FALSE,
                              xlab="Age (Months)",ylab="MUAC (cm)",cent=c(3,15,50,85,97),
                              main="Predicted girls percentiles using BCCG Model")
legend(315,35, c("3","15","50","85","97"), fill=(2:12), horiz=FALSE)


write.csv(BCCG2percentiles, file = "imputed_36000_BCCG2_percentiles.csv")

###########################
# Diagnostics
#########################
graphics.off()
plot(muac~agem, data=merged.girls, ylab="MUAC (cm)", xlab="Age (months)", 
     main="Girls BCCG regression line", cex=1.5, cex.lab=1.5, cex.axis = 1.5,cex.main=1.8)
lines(fitted(h4)~merged.girls$agem, col="red",lwd=3)

plot(h4) # BCCGo model
wp(h4, ylim.all = 4,xlim.all = 5,cex.lab=2.5,cex=1.5,cex.axis = 1.5,cex.main=1.8) # BCCG Model

fittedPlot(h4, x=merged.girls$agem) # BCCG


###########################################
### imputed4: only 5-6 years
#########################################
rm(imputed3)
rm(imputed3b)
rm(h4)
rm(merged.girls)
ls()



imputed4 <- read.csv("Imputed_5to6y_only_12000.csv", header = TRUE,as.is = TRUE)

summ(imputed4)
imputed4$agem<-imputed4$age/30 # convert age from days to months
imputed4$agey<-imputed4$age/365 # convert age from days to years
summ(imputed4)

# Create a girls dataset from the imputed1
tab1(imputed4$sex,graph=FALSE)
imputed4b<-imputed4[imputed4$sex=="Female",]
summ(imputed4b)

imputed4b$gzscores<-imputed4b$z
summ(imputed4b) 
imputed4b<-imputed4b[,c(2,4,5,7)]
summ(imputed4b)

# Merging the girls datasets}
summ(newf)
imputed4b$sex<-1
summ(imputed4b)

# Now, we are ready to combine these two girls datasets 
merged.girls=rbind(imputed4b,newf)
summ(merged.girls)

# Re-fit regression models to the merged data and assess goodness of fit and convergence

# h1 <- gamlss(muac~pb(agem),sigma.fo=~pb(agem), nu.fo=~pb(agem),tau.fo=~pb(agem),data=merged.girls, family=BCT)

#h2 <- lms(muac,agem ,data=merged.girls, method.pb="GAIC", k=log(length(merged.girls$muac)))
#h2$family
# chose BCTo

# h3 <- gamlss(muac~pb(agem),sigma.fo=~pb(agem),nu.fo=~pb(agem),tau.fo=~pb(agem), data=merged.girls, family=BCPE)
h4 <- lms(muac, agem ,data=merged.girls, method.pb="GAIC",
          k=log(length(merged.girls$muac)), families="BCCG")

#h5 <- gamlss(muac~pb(agem), sigma.formula=~pb(agem), family=BCCG, data=merged.girls) 
#plot(h)

# Choose model
# deviance(h4); deviance(h5) 
# AIC(h4,h5)
# AIC(h4,h5,k=log(length(merged.girls$muac)))
# We now generate Tables with LMS values using BCCG model (h2) and BCPE Model (h3) Using BCPEModel

summ(merged.girls$agem,graph=FALSE)
summ(merged.girls$muac,graph=FALSE)
gyf <-merged.girls$muac; length(gyf)

par(mar=c(5.1, 4.1, 4.1, 6.1), xpd=TRUE)
###
gyf <-merged.girls$muac; length(gyf)
gxf <-merged.girls$agem; length(gxf)

zscores_girls_BCCG3<- mycent.pred(h4, xname="agem",xvalues=24:300,plot=TRUE,
                                legend=FALSE,type="standard-centiles",ylim=c(10,80),xlim=c(24,300),
                                dev=c(-4,-3,-2,-1,0,1,2,3,4),  xlab="Age (Months)", ylab="MUAC (cm)", 
                                main="Predicted Z-values for girls using BCCG model")
legend(315,70, c("-4","-3","-2","-1","0", "1","2","3","4"), fill=c(2:10), horiz=FALSE)

write.csv(zscores_girls_BCCG3, file = "Imputed_5to6y_only_12000_BCCG3_Zscores.csv")

#################
# Percentiles
#############################
BCCG3_percentiles <- myper.pred(h4, xname="agem", plot=TRUE,xvalues=60:300,legend=FALSE,
                                xlab="Age (Months)",ylab="MUAC (cm)",cent=c(3,15,50,85,97),
                                main="Predicted girls percentiles using BCCG Model")
legend(315,35, c("3","15","50","85","97"), fill=(2:12), horiz=FALSE)


write.csv(BCCG3_percentiles, file = "Imputed_5to6y_only_12000_perce_BCCG3.csv")

#############################
## Diagnostics
############################
graphics.off()
plot(muac~agem, data=merged.girls, ylab="MUAC (cm)", xlab="Age (months)", 
     main="Girls BCCG regression line", cex=1.5, cex.lab=1.5, cex.axis = 1.5,cex.main=1.8)
lines(fitted(h4)~merged.girls$agem, col="red",lwd=3)

plot(h4) # BCCGo model
wp(h4, ylim.all = 4,xlim.all = 5,cex.lab=2.5,cex=1.5,cex.axis = 1.5,cex.main=1.8) # BCCG Model

fittedPlot(h4, x=merged.girls$agem) # BCCG
