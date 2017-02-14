
getwd()
setwd("U:/My Documents on U/Jay_MUAC_2017/muac_2017/merged_boys/Imputed")
library(gamlss) # package for fitting GAMLSS (now depends on gamlss.dist and gamlss.data
library(gamlss.dist) # package for all gamlss.family distributions
library(gamlss.data) #  package for all example data used in GAMLSS
library(gamlss.add) # experimental package for new additive terms.
library(gamlss.nl) # package for fitting nonlinear models
library(lattice)
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

# Create seperate datasets for boys only
nhanes$sex<-factor(nhanes$sex,labels=c("Female","Male"))
newf2 <-nhanes[nhanes$sex=="Male",]
summ(newf2 )
tab1(newf2 $sex, graph=FALSE)

# Now fit BCCG, BCPE and BCT GAMLSS regression models and compare model fit
# Once we select the best model, we use it to compute the Muac Z-scores
# Then we remove individuals considered to be outliers that present with Z scores 
# greater than or less than -4 SD, then merge this data with the imputed1 simulated data

# f1 <- gamlss(muac~pb(agem),sigma.fo=~pb(agem), nu.fo=~pb(agem),tau.fo=~pb(agem),data=newf2 , family=BCT)
# f2 <- lms(muac,agem ,data=newf2 , method.pb="GAIC", k=log(length(newf2 $muac)))
# f2$family
f3 <- gamlss(muac~pb(agem),sigma.fo=~pb(agem),nu.fo=~pb(agem),tau.fo=~pb(agem), data=newf2 , family=BCPE)
# f4 <- lms(muac, agem ,data=newf2 , method.pb="GAIC", k=log(length(newf2 $muac)), families="BCCG")

# Selection between models
## Residual deviances and AIC
# The GAIC function uses default penalty k = 2, giving the usual Akaike information criterion
# (AIC). Hence the usual AIC [equivalent to GAIC(k = 2)] selects model f1= BCT  as the best model
# (since it has the smallest value of AIC). 

# deviance(f1); deviance(f2);deviance(f3); deviance(f4)
# AIC(f1,f2,f3,f4)
# AIC(f1,f2,f3,f4, k=log(length(newf2 $muac)))# Bayesian Information Critetion
# centiles(f4,newf2 $agem, main="BCCG model",cent=c(3,15,50,85,97),lwd=3)

## Now we use Model f3 (BCPE)  and calculate Z scores from both and eliminate subjects with Z scores +/- 4 SD
gyf <-newf2$muac; length(gyf)
gxf <-newf2$agem; length(gxf)
gzscores <- centiles.pred(f3, xname="agem",xvalues=gxf,yval=gyf, type="z-scores")

head(newf2 )

newf2 <-cbind(newf2 ,gzscores)
head(newf2 )
max(newf2 $gzscores)
min(newf2 $gzscores)
hist(gzscores, xlab = "Z scores", main="NHanes Boys Z-scores for 2-26 years",cex.main=1.5,cex.lab=1.5)

## Removing newf2  with Z scores  +/- 4SD 
# There is none to remove sisnce they are all within the range
summ(newf2)



# Prepare to merge this data with the imputed1 simulated data and re-fit the models
dir(pattern = ".csv")
# imputed1 <- read.csv("imputed_12000.csv",header = TRUE,as.is = TRUE)
# imputed2 <- read.csv("imputed_24000.csv", header = TRUE,as.is = TRUE)
imputed3 <- read.csv("imputed_36000.csv", header = TRUE,as.is = TRUE)
imputed4 <- read.csv("Imputed_5to6y_only_12000.csv", header = TRUE,as.is = TRUE)

# summ(imputed1)

imputed1$agem<-imputed1$age/30 # convert age from days to months
imputed1$agey<-imputed1$age/365 # convert age from days to years
summ(imputed1)

# Create a boys dataset from the imputed1
tab1(imputed1$sex,graph=FALSE)
imputed1b<-imputed1[imputed1$sex=="Male",]
summ(imputed1b)

imputed1b$gzscores<-imputed1b$z
summ(imputed1b) 
imputed1b<-imputed1b[,c(2,4,5,7)]
summ(imputed1b)
 
# Merging the boys datasets}
summ(newf2)

imputed1b$sex<-2
summ(imputed1b)
summ(newf2)
newf2 <- newf2[,c(3,1,4,5)]
summ(newf2)
summ(imputed1b)

# Now, we are ready to combine these two boys datasets 
merged.boys=rbind(imputed1b,newf2)
summ(merged.boys)

# Re-fit regression models to the merged data and assess goodness of fit and convergence

# h1 <- gamlss(muac~pb(agem),sigma.fo=~pb(agem), nu.fo=~pb(agem),tau.fo=~pb(agem),data=merged.boys, family=BCT)
h2 <- lms(muac,agem ,data=merged.boys, method.pb="GAIC", k=log(length(merged.boys$muac)))
h2$family
# chose BCTo
# h3 <- gamlss(muac~pb(agem),sigma.fo=~pb(agem),nu.fo=~pb(agem),tau.fo=~pb(agem), data=merged.boys, family=BCPE)
h4 <- lms(muac, agem ,data=merged.boys, method.pb="GAIC", k=log(length(merged.boys$muac)), families="BCCG")


## Model selection
#deviance(h1); deviance(h2);deviance(h3); deviance(h4)
#AIC(h1,h2,h3,h4)
# AIC(h1,h2,h3, h4,k=log(length(merged.boys$muac)))


# We now generate Tables with LMS values using BCCG model (h2) and BCPE Model (h3) Using BCPEModel

summ(merged.boys$agem,graph=FALSE)
summ(merged.boys$muac,graph=FALSE)
gyf <-merged.boys$muac; length(gyf)
 
par(mar=c(5.1, 4.1, 4.1, 6.1), xpd=TRUE)

zscores_boys_BCCG<-centiles.pred(h4, xname="agem",xvalues=24:300,plot=TRUE,
                      legend=FALSE,type="standard-centiles",xlim=c(24,300),
                      dev=c(-4,-3,-2,-1,0,1,2,3,4),
                      xlab="Age (Months)",
                      ylab="MUAC (cm)", 
                      main="Predicted Z-values for boys using BCCG model")
legend(315,50, c("-4","-3","-2","-1","0", "1","2","3","4"), fill=c(2:10), horiz=FALSE)

write.csv(zscores_boys_BCCG,file = "imputed_12000_zscores_boys_BCCG.csv")

gyf <-merged.boys$muac; length(gyf)
gxf <-merged.boys$agem; length(gxf)

zscores_boys_BCTo<-mycent.pred(h2, xname="agem",xvalues=24:300,plot=TRUE,
                    legend=FALSE,type="standard-centiles",ylim=c(10,80),xlim=c(24,300),
                    dev=c(-4,-3,-2,-1,0,1,2,3,4),  xlab="Age (Months)", ylab="MUAC (cm)", 
                    main="Predicted Z-values for boys using BCTo model")
legend(315,70, c("-4","-3","-2","-1","0", "1","2","3","4"), fill=c(2:10), horiz=FALSE)

write.csv(zscores_boys_BCTo, file = "imputed_12000_zscores_boys_BCTo.csv")

################# 

BCCG_percentiles <- centiles.pred(h4, xname="agem", plot=TRUE,xvalues=24:300,legend=FALSE,
                              xlab="Age (Months)",ylab="MUAC (cm)",cent=c(3,15,50,85,97),
                              main="Predicted boys percentiles using BCCG Model")
legend(315,35, c("3","15","50","85","97"), fill=(2:12), horiz=FALSE)
 
write.csv(BCCG_percentiles, file = "imputed_12000_BCCG_percentiles.csv")
#############################
BCTo_percentiles <- myper.pred(h2, xname="agem", plot=TRUE,xvalues=24:300,legend=FALSE,
                               xlab="Age (Months)",ylab="MUAC (cm)",cent=c(3,15,50,85,97),
                               main="Predicted boys percentiles using BCTo Model")
legend(315,35, c("3","15","50","85","97"), fill=(2:12), horiz=FALSE)


write.csv(BCTo_percentiles, file = "imputed_12000_BCPEo_percentiles.csv")

#############################
## Imputed 2: 24000 data
#############################
rm(imputed1)
rm(imputed1b)
ls()

summ(imputed2)
imputed2$agem<-imputed2$age/30 # convert age from days to months
imputed2$agey<-imputed2$age/365 # convert age from days to years
summ(imputed2)

# Create a boys dataset from the imputed1
tab1(imputed2$sex,graph=FALSE)
imputed2b<-imputed2[imputed2$sex=="Male",]
summ(imputed2b)

imputed2b$gzscores<-imputed2b$z
summ(imputed2b) 
imputed2b<-imputed2b[,c(2,4,5,7)]
summ(imputed2b)

# Merging the boys datasets}
summ(newf2)

imputed2b$sex<-2
summ(imputed2b)
summ(newf2)
summ(newf2)
summ(imputed2b)

# Now, we are ready to combine these two boys datasets 
merged.boys=rbind(imputed2b,newf2)
summ(merged.boys)

# Re-fit regression models to the merged data and assess goodness of fit and convergence

# h1 <- gamlss(muac~pb(agem),sigma.fo=~pb(agem), nu.fo=~pb(agem),tau.fo=~pb(agem),data=merged.boys, family=BCT)
h2 <- lms(muac,agem ,data=merged.boys, method.pb="GAIC", k=log(length(merged.boys$muac)))
h2$family
# chose BCTo
# h3 <- gamlss(muac~pb(agem),sigma.fo=~pb(agem),nu.fo=~pb(agem),tau.fo=~pb(agem), data=merged.boys, family=BCPE)
h4 <- lms(muac, agem ,data=merged.boys, method.pb="GAIC", k=log(length(merged.boys$muac)), families="BCCG")

# We now generate Tables with LMS values using BCCG model (h2) and BCPE Model (h3) Using BCPEModel

summ(merged.boys$agem,graph=FALSE)
summ(merged.boys$muac,graph=FALSE)
gyf <-merged.boys$muac; length(gyf)

par(mar=c(5.1, 4.1, 4.1, 6.1), xpd=TRUE)

zscores_boys_BCCG<-centiles.pred(h4, xname="agem",xvalues=24:300,plot=TRUE,
                                 legend=FALSE,type="standard-centiles",xlim=c(24,300),
                                 dev=c(-4,-3,-2,-1,0,1,2,3,4),
                                 xlab="Age (Months)",
                                 ylab="MUAC (cm)", 
                                 main="Predicted Z-values for boys using BCCG model")
legend(315,50, c("-4","-3","-2","-1","0", "1","2","3","4"), fill=c(2:10), horiz=FALSE)

write.csv(zscores_boys_BCCG,file = "imputed_24000_zscores_boys_BCCG.csv")

gyf <-merged.boys$muac; length(gyf)
gxf <-merged.boys$agem; length(gxf)

zscores_boys_BCTo<-mycent.pred(h2, xname="agem",xvalues=24:300,plot=TRUE,
                               legend=FALSE,type="standard-centiles",ylim=c(10,80),xlim=c(24,300),
                               dev=c(-4,-3,-2,-1,0,1,2,3,4),  xlab="Age (Months)", ylab="MUAC (cm)", 
                               main="Predicted Z-values for boys using BCTo model")
legend(315,70, c("-4","-3","-2","-1","0", "1","2","3","4"), fill=c(2:10), horiz=FALSE)

write.csv(zscores_boys_BCTo, file = "imputed_24000_zscores_boys_BCTo.csv")

################# 

BCCG_percentiles <- centiles.pred(h4, xname="agem", plot=TRUE,xvalues=24:300,legend=FALSE,
                                  xlab="Age (Months)",ylab="MUAC (cm)",cent=c(3,15,50,85,97),
                                  main="Predicted boys percentiles using BCCG Model")
legend(315,35, c("3","15","50","85","97"), fill=(2:12), horiz=FALSE)

write.csv(BCCG_percentiles, file = "imputed_24000_BCCG_percentiles.csv")
#############################
BCTo_percentiles <- myper.pred(h2, xname="agem", plot=TRUE,xvalues=24:300,legend=FALSE,
                               xlab="Age (Months)",ylab="MUAC (cm)",cent=c(3,15,50,85,97),
                               main="Predicted boys percentiles using BCTo Model")
legend(315,35, c("3","15","50","85","97"), fill=(2:12), horiz=FALSE)


write.csv(BCTo_percentiles, file = "imputed_24000_BCPEo_percentiles.csv")



#############################
## imputed 3: 36000 data
##############################
rm(imputed2)
rm(imputed2b)
rm(merged.boys)
ls()

summ(imputed3)
imputed3$agem<-imputed3$age/30 # convert age from days to months
imputed3$agey<-imputed3$age/365 # convert age from days to years
summ(imputed3)

# Create a boys dataset from the imputed1
tab1(imputed3$sex,graph=FALSE)
imputed3b<-imputed3[imputed3$sex=="Male",]
summ(imputed3b)

imputed3b$gzscores<-imputed3b$z
summ(imputed3b) 
imputed3b<-imputed3b[,c(2,4,5,7)]
summ(imputed3b)

# Merging the boys datasets}
summ(newf2)

imputed3b$sex<-2
summ(imputed3b)
summ(newf2)
# newf2 <- newf2[,c(3,1,4,5)]

# Now, we are ready to combine these two boys datasets 
merged.boys=rbind(imputed3b,newf2)
summ(merged.boys)

# Re-fit regression models to the merged data and assess goodness of fit and convergence

# h1 <- gamlss(muac~pb(agem),sigma.fo=~pb(agem), nu.fo=~pb(agem),tau.fo=~pb(agem),data=merged.boys, family=BCT)
# h2 <- lms(muac,agem ,data=merged.boys, method.pb="GAIC", k=log(length(merged.boys$muac)))
# h2$family
# chose BCTo

# h3 <- gamlss(muac~pb(agem),sigma.fo=~pb(agem),nu.fo=~pb(agem),tau.fo=~pb(agem), data=merged.boys, family=BCPE)
h4 <- lms(muac, agem ,data=merged.boys, method.pb="GAIC", k=log(length(merged.boys$muac)), families="BCCG")

# We now generate Tables with LMS values using BCCG model (h2) and BCPE Model (h3) Using BCPEModel

summ(merged.boys$agem,graph=FALSE)
summ(merged.boys$muac,graph=FALSE)
gyf <-merged.boys$muac; length(gyf)

par(mar=c(5.1, 4.1, 4.1, 6.1), xpd=TRUE)

####
gyf <-merged.boys$muac; length(gyf)
gxf <-merged.boys$agem; length(gxf)

zscores_boys_BCTo<-mycent.pred(h4, xname="agem",xvalues=24:300,plot=TRUE,
                               legend=FALSE,type="standard-centiles",ylim=c(10,80),xlim=c(24,300),
                               dev=c(-4,-3,-2,-1,0,1,2,3,4),  xlab="Age (Months)", ylab="MUAC (cm)", 
                               main="Predicted Z-values for boys using BCCG model")
legend(315,70, c("-4","-3","-2","-1","0", "1","2","3","4"), fill=c(2:10), horiz=FALSE)

write.csv(zscores_boys_BCTo, file = "imputed_36000_zscores_boys_BCCG.csv")

#############################
BCCGpercentiles <- myper.pred(h4, xname="agem", plot=TRUE,xvalues=24:300,legend=FALSE,
                               xlab="Age (Months)",ylab="MUAC (cm)",cent=c(3,15,50,85,97),
                               main="Predicted boys percentiles using BCCG Model")
legend(315,35, c("3","15","50","85","97"), fill=(2:12), horiz=FALSE)


write.csv(BCCGpercentiles, file = "imputed_36000_BCCGo_percentiles.csv")

###########################
# Diagnostics
#########################
graphics.off()
plot(muac~agem, data=merged.boys, ylab="MUAC (cm)", xlab="Age (months)", 
     main="Boys BCCG regression line", cex=1.5, cex.lab=1.5, cex.axis = 1.5,cex.main=1.8)
lines(fitted(h4)~merged.boys$agem, col="red",lwd=3)

plot(h4) # BCCGo model
wp(h4, ylim.all = 4,xlim.all = 5,cex.lab=2.5,cex=1.5,cex.axis = 1.5,cex.main=1.8) # BCCG Model

fittedPlot(h4, x=merged.boys$agem) # BCCG


###########################################
### imputed4: only 5-6 years
#########################################
rm(imputed3)
rm(imputed3b)
rm(merged.boys)
ls()

imputed4 <- read.csv("Imputed_5to6y_only_12000.csv", header = TRUE,as.is = TRUE)

summ(imputed4)
imputed4$agem<-imputed4$age/30 # convert age from days to months
imputed4$agey<-imputed4$age/365 # convert age from days to years
summ(imputed4)

# Create a boys dataset from the imputed1
tab1(imputed4$sex,graph=FALSE)
imputed4b<-imputed4[imputed4$sex=="Male",]
summ(imputed4b)

imputed4b$gzscores<-imputed4b$z
summ(imputed4b) 
imputed4b<-imputed4b[,c(2,4,5,7)]
summ(imputed4b)

# Merging the boys datasets}
summ(newf2)
imputed4b$sex<-2
summ(imputed4b)

# Now, we are ready to combine these two boys datasets 
merged.boys=rbind(imputed4b,newf2)
summ(merged.boys)

# Re-fit regression models to the merged data and assess goodness of fit and convergence

# h1 <- gamlss(muac~pb(agem),sigma.fo=~pb(agem), nu.fo=~pb(agem),tau.fo=~pb(agem),data=merged.boys, family=BCT)

#h2 <- lms(muac,agem ,data=merged.boys, method.pb="GAIC", k=log(length(merged.boys$muac)))
#h2$family
# chose BCTo

# h3 <- gamlss(muac~pb(agem),sigma.fo=~pb(agem),nu.fo=~pb(agem),tau.fo=~pb(agem), data=merged.boys, family=BCPE)
h4 <- lms(muac, agem ,data=merged.boys, method.pb="GAIC", k=log(length(merged.boys$muac)), families="BCCG")

#h5 <- gamlss(muac~pb(agem), sigma.formula=~pb(agem), family=BCCG, data=merged.boys) 
#plot(h)

# Choose model
# deviance(h4); deviance(h5) 
# AIC(h4,h5)
# AIC(h4,h5,k=log(length(merged.boys$muac)))
# We now generate Tables with LMS values using BCCG model (h2) and BCPE Model (h3) Using BCPEModel

summ(merged.boys$agem,graph=FALSE)
summ(merged.boys$muac,graph=FALSE)
gyf <-merged.boys$muac; length(gyf)

par(mar=c(5.1, 4.1, 4.1, 6.1), xpd=TRUE)
###
gyf <-merged.boys$muac; length(gyf)
gxf <-merged.boys$agem; length(gxf)

zscores_boys_BCCG<- mycent.pred(h4, xname="agem",xvalues=24:300,plot=TRUE,
                               legend=FALSE,type="standard-centiles",ylim=c(10,80),xlim=c(24,300),
                               dev=c(-4,-3,-2,-1,0,1,2,3,4),  xlab="Age (Months)", ylab="MUAC (cm)", 
                               main="Predicted Z-values for boys using BCCG model")
legend(315,70, c("-4","-3","-2","-1","0", "1","2","3","4"), fill=c(2:10), horiz=FALSE)

write.csv(zscores_boys_BCCG, file = "Imputed_5to6y_only_12000_BCCG_Zscores.csv")

#################
# Percentiles
#############################
BCCG1_percentiles <- myper.pred(h4, xname="agem", plot=TRUE,xvalues=60:300,legend=FALSE,
                               xlab="Age (Months)",ylab="MUAC (cm)",cent=c(3,15,50,85,97),
                               main="Predicted boys percentiles using BCCG Model")
legend(315,35, c("3","15","50","85","97"), fill=(2:12), horiz=FALSE)


write.csv(BCCG1_percentiles, file = "Imputed_5to6y_only_12000_perce_BCCG.csv")
 
#############################
## Diagnostics
############################

plot(muac~agem, data=merged.boys, ylab="MUAC (cm)", xlab="Age (months)", 
     main="Boys BCCG regression line", cex=1.5, cex.lab=1.5, cex.axis = 1.5,cex.main=1.8)
lines(fitted(h4)~merged.boys$agem, col="red",lwd=3)

plot(h4) # BCCGo model
wp(h4, ylim.all = 4,xlim.all = 5,cex.lab=2.5,cex=1.5,cex.axis = 1.5,cex.main=1.8) # BCCG Model

 fittedPlot(h4, x=merged.boys$agem) # BCCG


 
 
 