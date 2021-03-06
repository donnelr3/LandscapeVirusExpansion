# "The role of pathogen mediated insect superabundance in the east-African emergence of a plant virus" Fig. 3 Fig. 1 Table 1
#############################################################################################################################
#### COMMENTS                                                                   R. Donnelly 2021
# As decribed in main text and summarised in results Table 1, there are two steps tothe analysis
# a) use LRT to motivate the quadratic regression model (quadratic wrt transect distance)
# b) analyse the data using Bayesian regression based upon the quadratic model (this calls "waveProfileMEmodel.stan")
# From b) we extract credible intervals for quadratic curvature and for turning point
# This process is applied to the following 3 datasets:
################### 1 # Landscape SURVEY fitting CENTRAL TRANSECT ###########################################################
################### 2 # Landscape SURVEY fitting  EASTERN TRANSECT ##########################################################
################### 3 # Landscape EXPERIMENT fitting  #######################################################################
#############################################################################################################################


rm(list=setdiff(ls(), "filepath"))
######################################################################
require(lme4)
library('nlme')  
library('MASS')  
library(pracma) 
library(rstan) 
library(lmtest)
library(bridgesampling)
library(data.table)
options(mc.cores = parallel::detectCores()) 
Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7 -mtune=corei7')
################################################################################################################



################################################################################################################
################### 1 # Landscape survey fitting CENTRAL TRANSECT ##############################################
################################################################################################################
thefile=paste("DigitiseFiles_Survey/SurveyCentralOrientS.csv",sep="")  # data from Legg et al. 1998 digitised in reference folder 
ugData<-read.table(file=thefile, header=TRUE, sep=",")
print(head(ugData))
datasize<-dim(ugData)

# data wrangling
ugIncOnlypre=ugData[ugData$Type=='inc',]
ugAdultOnlypre=ugData[ugData$Type=='adult',]
# aggregating data over 1 yr as described in Supp Inf 4
ugAvrAdValueC=colMeans(rbind(ugAdultOnlypre[ugAdultOnlypre$Year==1,1],ugAdultOnlypre[ugAdultOnlypre$Year==2,1],ugAdultOnlypre[ugAdultOnlypre$Year==3,1],ugAdultOnlypre[ugAdultOnlypre$Year==4,1],ugAdultOnlypre[ugAdultOnlypre$Year==5,1]))
ugAvrIncValueC=colMeans(rbind(ugIncOnlypre[ugIncOnlypre$Year==1,1],ugIncOnlypre[ugIncOnlypre$Year==2,1],ugIncOnlypre[ugIncOnlypre$Year==3,1],ugIncOnlypre[ugIncOnlypre$Year==4,1],ugIncOnlypre[ugIncOnlypre$Year==5,1]))
# wave profile
avrAbunByIncC=ugAvrAdValueC/(ugAvrIncValueC+1)
temp=ugAdultOnlypre[ugAdultOnlypre$Year==1,c(4,7,8,9,10,11)]
row.names(temp)=seq(1,Dim(temp)[1])
ugAdultAvrOnly=data.frame(temp,avrAbunByIncC)
print(head(ugAdultAvrOnly))

# Truncate the data as described in Supp Inf 4 (effect of urban connurbations, Legg '98)
ugAdultAvrOnlyTrunc=ugAdultAvrOnly[ugAdultAvrOnly$Field>2,]
ugAdultAvrOnly$Field=as.factor(ugAdultAvrOnly$Field)
Year=rep("one",Dim(ugAdultAvrOnlyTrunc)[1])
ugAdultAvrOnlyTrunc=data.frame(ugAdultAvrOnlyTrunc,Year)

# likelihood ratio tests motivate quadratic model
regrCubicA <- lm(avrAbunByIncC ~ Fieldkm+Field2km+Field3km, data = ugAdultAvrOnlyTrunc)
regrQuadrA <- lm(avrAbunByIncC ~ Fieldkm+Field2km, data = ugAdultAvrOnlyTrunc)
regrLinA <- lm(avrAbunByIncC ~ Fieldkm, data = ugAdultAvrOnlyTrunc)
linQuadA=lrtest(regrLinA,regrQuadrA)
quadCubA=lrtest(regrQuadrA,regrCubicA)
part1ofA=rbind(linQuadA,quadCubA)
rownames(part1ofA)=c(" ","Table 1 A col.2(i)","","Table 1 A col.2(ii)")
# creating data structure for stan model
adjFields=as.integer(ugAdultAvrOnlyTrunc$Field)-2
XX=cbind(rep(1,8),ugAdultAvrOnlyTrunc$Fieldkm,ugAdultAvrOnlyTrunc$Field2km)
DNF_S=length(unique(ugAdultAvrOnlyTrunc$Field));
DL_S=length(unique(ugAdultAvrOnlyTrunc$Year));
dat2= list(D_Nfields=DNF_S,
           D_levels=DL_S,
           D_K=ncol(XX),
           D_X=t(as.matrix(unique(XX))),
           D_y=matrix(ugAdultAvrOnlyTrunc$avrAbunByIncC, byrow = TRUE, nrow = DL_S, ncol = DNF_S),
           D_subj = adjFields,
           D_intcpt=1   # defaulting to 1, not relevant for the single yr survey data
)

### Setting stan fitting parameters
numWarm=1000
numIter=4000
adaptVal=1-(10^-2)
numChains=4
treeDepth=15

### Running the stan model
fitA = stan(file = "waveProfileMEmodel.stan",
            data = dat2,
            iter = numIter,
            warmup = numWarm,
            chains=numChains,
            control = list(max_treedepth = treeDepth,adapt_delta=adaptVal) )
stanOutA=summary(fitA)
stanDetailA=stanOutA$summary
### summarising stan results
resultsTableA=rbind(stanDetailA[rownames(stanDetailA)=="beta1[1]",][c(1,4,6,8)],
               stanDetailA[rownames(stanDetailA)=="beta1[2]",][c(1,4,6,8)],
               stanDetailA[rownames(stanDetailA)=="beta1[3]",][c(1,4,6,8)],
               stanDetailA[rownames(stanDetailA)=="compoundPam[1]",][c(1,4,6,8)])
row.names(resultsTableA)=c("beta[1]","beta[2]","beta[3]        (Table 1B col.3(i))","compoundPam[1] (Table 1B col.3(ii))")
resultsTableA
confint(regrQuadrA)
################################################################################################################









################################################################################################################
################### 2 # Landscape survey fitting  EASTERN TRANSECT #############################################
################################################################################################################
thefile=paste("DigitiseFiles_Survey/SurveyEasternOrientS.csv",sep="") # data from Legg et al. 1998 digitised in reference folder 
ugData<-read.table(file=thefile, header=TRUE, sep=",")
print(head(ugData))
datasize<-dim(ugData)

# data wrangling
ugIncOnlypre=ugData[ugData$Type=='inc',]
ugAdultOnlypre=ugData[ugData$Type=='adult',]

# aggregating data over 1 yr as described in Supp Inf 4
ugAvrAdValueE=colMeans(rbind(ugAdultOnlypre[ugAdultOnlypre$Year==1,1],ugAdultOnlypre[ugAdultOnlypre$Year==2,1],ugAdultOnlypre[ugAdultOnlypre$Year==3,1],ugAdultOnlypre[ugAdultOnlypre$Year==4,1]))
ugAvrIncValueE=colMeans(rbind(ugIncOnlypre[ugIncOnlypre$Year==1,1],ugIncOnlypre[ugIncOnlypre$Year==2,1],ugIncOnlypre[ugIncOnlypre$Year==3,1],ugIncOnlypre[ugIncOnlypre$Year==4,1]))

# wave profile
avrAbunByIncE=ugAvrAdValueE/(ugAvrIncValueE+1)
temp=ugAdultOnlypre[ugAdultOnlypre$Year==1,c(4,7,8,9,10,11)]
row.names(temp)=seq(1,Dim(temp)[1])
ugAdultAvrOnly=data.frame(temp,avrAbunByIncE)
print(head(ugAdultAvrOnly))

# Truncate the data as described in Supp Inf 4 (effect of urban connurbations, Legg '98)
ugAdultAvrOnlyTrunc=ugAdultAvrOnly[ugAdultAvrOnly$Field>2,]
ugAdultAvrOnly$Field=as.factor(ugAdultAvrOnly$Field)
Year=rep("one",Dim(ugAdultAvrOnlyTrunc)[1])
ugAdultAvrOnlyTrunc=data.frame(ugAdultAvrOnlyTrunc,Year)

# likelihood ratio tests motivate quadratic model
regrCubicB <- lm(avrAbunByIncE ~ Fieldkm+Field2km+Field3km, data = ugAdultAvrOnlyTrunc)
regrQuadrB <- lm(avrAbunByIncE ~ Fieldkm+Field2km, data = ugAdultAvrOnlyTrunc)
regrLinB <- lm(avrAbunByIncE ~ Fieldkm, data = ugAdultAvrOnlyTrunc)
linQuadB=lrtest(regrLinB,regrQuadrB)
quadCubB=lrtest(regrQuadrB,regrCubicB)
part1ofB=rbind(linQuadB,quadCubB)
rownames(part1ofB)=c(" ","Table 1 A col.3(i)","","Table 1 A col.3(ii)")

# creating data structure for stan model
adjFields=as.integer(ugAdultAvrOnlyTrunc$Field)-2    # reorder fields from 2:10 to 1:8
XX=cbind(rep(1,8),ugAdultAvrOnlyTrunc$Fieldkm,ugAdultAvrOnlyTrunc$Field2km)
DNF_S=length(unique(ugAdultAvrOnlyTrunc$Field));
DL_S=length(unique(ugAdultAvrOnlyTrunc$Year));
dat2= list(D_Nfields=DNF_S,
           D_levels=DL_S,
           D_K=ncol(XX),
           D_X=t(as.matrix(unique(XX))),
           D_y=matrix(ugAdultAvrOnlyTrunc$avrAbunByIncE, byrow = TRUE, nrow = DL_S, ncol = DNF_S),
           D_subj = adjFields,
           D_intcpt=1   # defaulting to 1, not relevant for the single yr survey data
)
### Run the stan model 
fitB = stan(file = "waveProfileMEmodel.stan",
            data = dat2,
            iter = numIter,
            warmup = numWarm,
            chains=numChains,
            control = list(max_treedepth = treeDepth,adapt_delta=adaptVal) )
stanOutB=summary(fitB)
stanDetailB=stanOutB$summary
### summarising stan results
resultsTableB=rbind(stanDetailB[rownames(stanDetailB)=="beta1[1]",][c(1,4,6,8)],
               stanDetailB[rownames(stanDetailB)=="beta1[2]",][c(1,4,6,8)],
               stanDetailB[rownames(stanDetailB)=="beta1[3]",][c(1,4,6,8)],
               stanDetailB[rownames(stanDetailB)=="compoundPam[1]",][c(1,4,6,8)])
row.names(resultsTableB)=c("beta[1]","beta[2]","beta[3]        (Table 1B col.4(i))","compoundPam[1] (Table 1B col.4(ii))")

resultsTableB
confint(regrQuadrB)
################################################################################################################














################################################################################################################
################### 3 # Landscape experiment fitting  ##########################################################
################################################################################################################

thefile=paste("DigitiseFiles_Experiment/ExperimentOrientS.csv",sep="") 
expData<-read.table(file=thefile, header=TRUE, sep=",")

print(head(expData))
datasize<-dim(expData)

# data wrangling
expDataProfile=expData
expIncOnly=expData[expData$Type=='inc',]
expAdultOnly=expData[expData$Type=='adult',]
expNymphOnly=expData[expData$Type=='nymph',]

# wave-profile
expadultAbunByInc=expAdultOnly$Value/(expIncOnly$Value+1)
expnymAbunByInc=expNymphOnly$Value/(expIncOnly$Value+1)
tempy=expAdultOnly[,c(3,4,7,8,9)]
row.names(tempy)=seq(1,Dim(tempy)[1])
expAdultOnlyDF=data.frame(tempy,expadultAbunByInc)
expNymphOnlyDF=data.frame(tempy,expnymAbunByInc)

# likelihood ratio tests motivate quadratic model
regLinME <- lmer(expadultAbunByInc ~ Year:(Fieldkm)+(1|Field), data = expAdultOnlyDF, REML=FALSE)
regQuadrME <- lmer(expadultAbunByInc ~ Year:(Fieldkm+Field2km)+(1|Field), data = expAdultOnlyDF, REML=FALSE)
regCubicME <- lmer(expadultAbunByInc ~ Year:(Fieldkm+Field2km+Field3km)+(1|Field), data = expAdultOnlyDF, REML=FALSE)
linQuadC=lrtest(regLinME,regQuadrME) # lrt not appropriate for comparing nested FE models in ME setting
quadCubC=lrtest(regQuadrME,regCubicME) # lrt not appropriate for comparing nested FE models in ME setting
part1ofC=rbind(linQuadC,quadCubC)
rownames(part1ofC)=c(" ","Table 1 A col.1(i)","","Table 1 A col.1(ii)")

# CAPTION REPORTING
regQuadrME_SI <- lmer(expadultAbunByInc ~ Year+Year:(Fieldkm+Field2km)+(1|Field), data = expAdultOnlyDF,REML=FALSE)
regQuadrME_CI <- lmer(expadultAbunByInc ~ Year:(Fieldkm+Field2km)+(1|Field), data = expAdultOnlyDF,REML=FALSE)
lrtest(regQuadrME_CI,regQuadrME_SI) # lrt not appropriate for comparing nested FE models in ME setting

# ME fits for comparison with Bayesian
regQuadrME_SI_REML <- lmer(expadultAbunByInc ~ Year+Year:(Fieldkm+Field2km)+(1|Field), data = expAdultOnlyDF)
regQuadrME_CI_REML <- lmer(expadultAbunByInc ~ Year:(Fieldkm+Field2km)+(1|Field), data = expAdultOnlyDF)

# creating data structure for stan model
numFields=size(unique(expAdultOnlyDF$Field))[2]
XX=cbind(rep(1,numFields),expAdultOnlyDF$Fieldkm,expAdultOnlyDF$Field2km)
DNF_E=length(unique(expAdultOnlyDF$Field));
DL_E=length(unique(expAdultOnlyDF$Year));
dat2= list(D_Nfields=DNF_E,
           D_levels=DL_E,
           D_K=ncol(XX),
           D_X=t(as.matrix(unique(XX))),
           D_y=matrix(expAdultOnlyDF$expadultAbunByInc, byrow = TRUE, nrow = DL_E, ncol = DNF_E),
           D_subj = as.integer(expAdultOnlyDF$Field),
           D_intcpt=1    #shared intercept over years has value 1; value 2 has separate intercepts
)
### Run the stan model
fitC = stan(file = "waveProfileMEmodel.stan",
            data = dat2,
            iter = numIter,
            warmup = numWarm,
            chains=numChains,
            control = list(max_treedepth = treeDepth,adapt_delta=0.99) )
### summarising stan results
stanOutC=summary(fitC)
stanDetailC=stanOutC$summary
resultsTableC=rbind(stanDetailC[rownames(stanDetailC)=="beta1[1]",][c(1,4,6,8)],
                    stanDetailC[rownames(stanDetailC)=="beta1[2]",][c(1,4,6,8)],
                    stanDetailC[rownames(stanDetailC)=="beta2[1]",][c(1,4,6,8)],
                    stanDetailC[rownames(stanDetailC)=="beta1[3]",][c(1,4,6,8)],
                    stanDetailC[rownames(stanDetailC)=="beta2[2]",][c(1,4,6,8)],
                    stanDetailC[rownames(stanDetailC)=="compoundPam[1]",][c(1,4,6,8)],
                    stanDetailC[rownames(stanDetailC)=="compoundPam[2]",][c(1,4,6,8)],
                    stanDetailC[rownames(stanDetailC)=="sigma_eps",][c(1,4,6,8)],
                    stanDetailC[rownames(stanDetailC)=="sigma_id[1]",][c(1,4,6,8)])
row.names(resultsTableC)=c("beta1[1]","beta1[2]","beta2[1]","beta1[3]       (Table 1B col.1(i))","beta2[2]       (Table 1B col.2(i))","compoundPam[1] (Table 1B col.1(ii))","compoundPam[2] (Table 1B col.2(ii))","sigma_eps","sigma_id")

sum_regQuadrME_CI_REML=summary(regQuadrME_CI_REML)
sum_regQuadrME_CI_REML$coefficients[,1]
################################################################################################################


#RESULTS - as presented in Table 1 main text (SEE ROWNAMES IN PRINTED TABLES BELOW)
part1ofA
part1ofB
part1ofC

resultsTableA
resultsTableB
resultsTableC

summary(fitA, probs = c(0.04, 0.05, 0.95, 0.96))$summary

pname <- "sigma_id[1]"
muc <- rstan::extract(fitC, pars=pname,  permuted=FALSE, inc_warmup=FALSE)
mdf <- melt(muc)
ggplot(mdf,aes(x=iterations,y=value,color=chains)) + geom_line() + ylab(mdf$parameters[1])



save.image()