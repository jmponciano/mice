# Comparisons for Mike
source("CompareTools.R")

intranasal.dat <- data.frame(read.csv("~/Documents/UFL/RESEARCH/JasonBlackburn/MikeNorris/Datasets/IntranasalChallenge.csv", header=TRUE))

subcutaneous.dat <- data.frame(read.csv("~/Documents/UFL/RESEARCH/JasonBlackburn/MikeNorris/Datasets/SubcutaneousChallenge.csv", header=TRUE))



# Testing
guess <- my.guess <- c(0.5, 0.6)

intranasalW <- intranasal.dat[intranasal.dat$Strain=="W",]
intranasalM <- intranasal.dat[intranasal.dat$Strain=="M",]
intranasalC <- intranasal.dat[intranasal.dat$Strain=="C",]

subcutanW <- subcutaneous.dat[subcutaneous.dat$Strain=="W",]
subcutanM <- subcutaneous.dat[subcutaneous.dat$Strain=="M",]
subcutanC <- subcutaneous.dat[subcutaneous.dat$Strain=="C",]

# Intranasal and subcutaneous W vs M, all doses: 
intranasalWVsM <- joint.vs.sep.fit(dataframe1=intranasalW, dataframe2=intranasalM)

subcutanWVsM <- joint.vs.sep.fit(dataframe1=subcutanW, dataframe2=subcutanM)

######## Intranasal and subcutaneous W vs M, dose 1000000:
intranasalW1e6 <- intranasalW[intranasalW$Dose==1000000,] 
intranasalM1e6 <- intranasalM[intranasalM$Dose==1000000,] 
intranasalWVsM.1e6 <- joint.vs.sep.fit(dataframe1=intranasalW1e6, dataframe2=intranasalM1e6)

subcutanW1e6 <- subcutanW[subcutanW$Dose==1000000,] 
subcutanM1e6 <- subcutanM[subcutanM$Dose==1000000,] 
subcutanWVsM.1e6 <- joint.vs.sep.fit(dataframe1=subcutanW1e6, dataframe2=subcutanM1e6)

######## Intranasal and subcutaneous W vs M, dose 10000:
intranasalW1e4 <- intranasalW[intranasalW$Dose==10000,] 
intranasalM1e4 <- intranasalM[intranasalM$Dose==10000,] 
intranasalWVsM.1e4 <- joint.vs.sep.fit(dataframe1=intranasalW1e4, dataframe2=intranasalM1e4)

subcutanW1e4 <- subcutanW[subcutanW$Dose==10000,] 
subcutanM1e4 <- subcutanM[subcutanM$Dose==10000,] 
subcutanWVsM.1e4 <- joint.vs.sep.fit(dataframe1=subcutanW1e4, dataframe2=subcutanM1e4)

######## Intranasal and subcutaneous W vs M, dose 1000:
intranasalW1e3 <- intranasalW[intranasalW$Dose==1000,] 
intranasalM1e3 <- intranasalM[intranasalM$Dose==1000,] 
intranasalWVsM.1e3 <- joint.vs.sep.fit(dataframe1=intranasalW1e3, dataframe2=intranasalM1e3)

subcutanW1e3 <- subcutanW[subcutanW$Dose==1000,] 
subcutanM1e3 <- subcutanM[subcutanM$Dose==1000,] 
subcutanWVsM.1e3 <- joint.vs.sep.fit(dataframe1=subcutanW1e3, dataframe2=subcutanM1e3)

######## Intranasal and subcutaneous W vs M, dose 10:
intranasalW10 <- intranasalW[intranasalW$Dose==10,] 
intranasalM10 <- intranasalM[intranasalM$Dose==10,] 
intranasalWVsM.10 <- joint.vs.sep.fit(dataframe1=intranasalW10, dataframe2=intranasalM10)

subcutanW10 <- subcutanW[subcutanW$Dose==10,] 
subcutanM10 <- subcutanM[subcutanM$Dose==10,] 
subcutanWVsM.10 <- joint.vs.sep.fit(dataframe1=subcutanW10, dataframe2=subcutanM10)

############### Subcutaneous 1000 -replicates analyses ###############

subcutanreps1K.dat <- data.frame(read.csv("~/Documents/UFL/RESEARCH/JasonBlackburn/MikeNorris/Datasets/OneK-CFU-datareps.csv", header=TRUE))

subcutan1KW <- subcutanreps1K.dat[subcutanreps1K.dat$Strain=="W",]
subcutan1KM <- subcutanreps1K.dat[subcutanreps1K.dat$Strain=="M",]

subcutanWVsM1K.reps <- joint.vs.sep.fit(dataframe1=subcutan1KW, dataframe2=subcutan1KM,Nreps=TRUE)

Wmles <- subcutanWVsM1K.reps$data1CIs.mat
Mmles <- subcutanWVsM1K.reps$data2CIs.mat

write.csv(Wmles, file="Wmles.csv")
write.csv(Mmles, file="Mmles.csv")

Wpreds <- Expected.pred(mlesmat=Wmles, no=5,max.days=16,dose=1000)
Mpreds <- Expected.pred(mlesmat=Mmles, no=5,max.days=16,dose=1000)

no<- 5
tvec <- Wpreds[,1]
Wpred.hat <- Wpreds[,3]
Mpred.hat <- Mpreds[,3]
par(oma=c(1,1,1,1), mar=c(5,5,3,3))
plot(tvec, Wpred.hat, type="l", col="blue", lwd=2, ylim=c(0,no), ylab="Number of alive mice", xlab="Days", cex.lab=1.5, main="Predicted deaths for the 1000 CFU data replicates")
points(tvec, Wpreds[,2], type="l", col="blue",lty=2, lwd=1)
points(tvec, Wpreds[,4], type="l", col="blue",lty=2, lwd=1)
points(tvec, Mpred.hat, type="l", col="red", lwd=2)
points(tvec, Mpreds[,2], type="l", col="red",lty=2, lwd=1)
points(tvec, Mpreds[,4], type="l", col="red",lty=2, lwd=1)
legend("bottomleft", c("Wild Type", "Mutant"), lty=c(1,1), lwd=c(2,2), col=c("blue","red"), cex=1.25, bty="n")


####### And plot simulations side by side these estimations with the two sets of parameters:

####  Simulating the exact time of death
Wsim <- deaths.sims(mles=Wmles[,2], no=5,dose=1000)
Msim <- deaths.sims(mles=Mmles[,2], no=5, dose=1000)

####  Simulating the observations every day, at the same time
Wsim.samps <- Sampling.DP(AllDeaths=Wsim$Nalive, cumtimes=Wsim$cumtimes, len=16,tau=1)
Msim.samps <- Sampling.DP(AllDeaths=Msim$Nalive, cumtimes=Msim$cumtimes, len=16,tau=1)

no<- 5

tvec <- Wpreds[,1]
Wpred.hat <- Wpreds[,3]
Mpred.hat <- Mpreds[,3]
###### Fig 5 c???
svg("Figure5c.svg")
par(oma=c(2,2,2,2))
plot(tvec, Wpred.hat, type="l", col="blue", lwd=2, ylim=c(0,no), ylab="Number of mice alive", xlab="Days", cex.lab=1.5, main="Predicted deaths,1000 CFU replicates", xlim=c(0,17))
points(tvec, Wpreds[,2], type="l", col="blue",lty=2, lwd=1)
points(tvec, Wpreds[,4], type="l", col="blue",lty=2, lwd=1)
points(tvec, Mpred.hat, type="l", col="red", lwd=2)
points(tvec, Mpreds[,2], type="l", col="red",lty=2, lwd=1)
points(tvec, Mpreds[,4], type="l", col="red",lty=2, lwd=1)
legend("bottomleft", c("Wild Type", "Mutant"), lty=c(1,1), lwd=c(2,2), col=c("blue","red"), cex=1.25, bty="n")
dev.off()

svg("FigureS8a.svg")
par(oma=c(2,2,2,2))
plot(Wsim.samps$times,Wsim.samps$obs.ts, col="blue",type="s", lwd=2,ylim=c(0,no), 
     ylab="Number of mice alive", xlab="Days", cex.lab=1.5, main="Simulated deaths,1000 CFU replicates", xlim=c(0,17))
points(Wsim.samps$times,Msim.samps$obs.ts, col="red",type="s", lwd=2)
legend("bottomleft", c("Wild Type", "Mutant"), lty=c(1,1), lwd=c(2,2), col=c("blue","red"), cex=1.25, bty="n")
dev.off()

############  

#####  Example of simulating data set like the one observed, using the W and M types mles
Wdatasim <- dataset.sim(mles=Wmles[,2], no=5, nreps=3, id.reps=1:3, days=seq(0,15,by=3), Strain.name="W", Dose=1000)
Mdatasim <- dataset.sim(mles=Mmles[,2], no=5, nreps=3, id.reps=4:6,days=seq(0,15,by=3), Strain.name="M", Dose=1000)
#####  Example of testing the difference between the two simulated data sets
sim.diff.test <- joint.vs.sep.fit(dataframe1=Wdatasim, dataframe2=Mdatasim,Nreps=TRUE)


####### Out of 100 times, how many would we detect a difference when it is actually there
####### for initial number of mice of 5, 10, 20, 40, 80, 160 if we were to run the experiment for 15 days?

Nsims <- 100
success.counter <- rep(0, Nsims)

Strong <- "The evidence for separate death rates is conclusive"
Weak <- "Separate death rates is a better model but weakly so."

#A <- "The evidence for separate death rates is conclusive"
#A <- "Separate death rates is a better model but weakly so."
#A <- "Something else"

for(i in 1:Nsims){
	
	Wdatasim <- dataset.sim(mles=Wmles[,2], no=5, nreps=3, id.reps=1:3, days=seq(0,15,by=3), 
	Strain.name="W", Dose=1000)
	Mdatasim <- dataset.sim(mles=Mmles[,2], no=5, nreps=3, id.reps=4:6,days=seq(0,15,by=3), 
	Strain.name="M", Dose=1000)
	sim.diff.test <- joint.vs.sep.fit(dataframe1=Wdatasim, dataframe2=Mdatasim,Nreps=TRUE)
	A <- sim.diff.test$Best.model
	success.counter[i] <- (A==Strong)||(A==Weak)		
	
}

prop.success <- sum(success.counter)/Nsims # 0.13 (13%), days = 0  3  6  9 12 15

######  Now sampling every day

Nsims <- 100
success.counter <- rep(0, Nsims)

Strong <- "The evidence for separate death rates is conclusive"
Weak <- "Separate death rates is a better model but weakly so."
for(i in 1:Nsims){
	
	Wdatasim <- dataset.sim(mles=Wmles[,2], no=5, nreps=3, id.reps=1:3, days=0:15, 
	Strain.name="W", Dose=1000)
	Mdatasim <- dataset.sim(mles=Mmles[,2], no=5, nreps=3, id.reps=4:6,days=0:15, 
	Strain.name="M", Dose=1000)
	sim.diff.test <- joint.vs.sep.fit(dataframe1=Wdatasim, dataframe2=Mdatasim,Nreps=TRUE)
	A <- sim.diff.test$Best.model
	success.counter[i] <- (A==Strong)||(A==Weak)		
	
}

prop.success <- sum(success.counter)/Nsims # 0.12 (12%)
prop.success

##########  Now changing the initial number of mice

Nsims <- 100
Nmice <- c(5, 10, 20, 40, 80)
num.Nos <- length(Nmice)
power.vec <- rep(0,num.Nos)
Strong <- "The evidence for separate death rates is conclusive"
Weak <- "Separate death rates is a better model but weakly so."

for(j in 1:num.Nos){
	success.counter <- rep(0, Nsims)
	jth.no <- Nmice[j]
	for(i in 1:Nsims){
		Wdatasim <- dataset.sim(mles=Wmles[,2], no=jth.no, nreps=3, id.reps=1:3, days=0:40, 
		Strain.name="W", Dose=1000)
		Mdatasim <- dataset.sim(mles=Mmles[,2], no=jth.no, nreps=3, id.reps=4:6,days=0:40, 
		Strain.name="M", Dose=1000)
		sim.diff.test <- joint.vs.sep.fit(dataframe1=Wdatasim, dataframe2=Mdatasim,Nreps=TRUE)
		A <- sim.diff.test$Best.model
		success.counter[i] <- (A==Strong)||(A==Weak)		
	}
	power.vec[j] <- sum(success.counter)/Nsims 
}
	
Num.mice =    5,   10,   20,   40,   80
Power    = 0.02, 0.22, 0.73, 0.99, 1.00


#############  Now sampling only to 20 days:
##########  Now changing the initial number of mice

Nsims <- 100
Nmice <- c(5, 10, 20, 40, 80)
num.Nos <- length(Nmice)
power.vec <- rep(0,num.Nos)
Strong <- "The evidence for separate death rates is conclusive"
Weak <- "Separate death rates is a better model but weakly so."

for(j in 1:num.Nos){
	success.counter <- rep(0, Nsims)
	jth.no <- Nmice[j]
	for(i in 1:Nsims){
		Wdatasim <- dataset.sim(mles=Wmles[,2], no=jth.no, nreps=3, id.reps=1:3, days=0:20, 
		Strain.name="W", Dose=1000)
		Mdatasim <- dataset.sim(mles=Mmles[,2], no=jth.no, nreps=3, id.reps=4:6,days=0:20, 
		Strain.name="M", Dose=1000)
		sim.diff.test <- joint.vs.sep.fit(dataframe1=Wdatasim, dataframe2=Mdatasim,Nreps=TRUE)
		A <- sim.diff.test$Best.model
		success.counter[i] <- (A==Strong)||(A==Weak)		
	}
	power.vec[j] <- sum(success.counter)/Nsims 
}

Nmice
power.vec

#> Nmice
#[1]  5 10 20 40 80
#> power.vec
#[1] 0.06 0.13 0.57 0.92 1.00

# Draw power curve. 
par(oma=c(1,1,1,1), mar=c(5,5,3,2))
plot(Nmice, power.vec, pch=16, xlab="Number of mice per batch (3 replicates)", ylab="Pr(Correctly detecting a faster mutant death rate)", cex.lab=1.5)

save.image("MiceCalcs.RData")

######################  Galleria data analysis #####################

galleriaM.dat <- data.frame(read.csv("~/Documents/UFL/RESEARCH/JasonBlackburn/MikeNorris/Datasets/GalleriaDP-M.csv", header=TRUE))

galleriaW.dat <- data.frame(read.csv("~/Documents/UFL/RESEARCH/JasonBlackburn/MikeNorris/Datasets/GalleriaDP-W.csv", header=TRUE))

galleria.test <- joint.vs.sep.fit(dataframe1=galleriaW.dat, dataframe2=galleriaM.dat)


Wmles <- galleria.test$data1CIs.mat
Mmles <- galleria.test$data2CIs.mat

#######  Power simulation using the Galleria MLEs at different initial sample sizes:

Nsims <- 1000
Nmice <- c(5, 10, 20, 40, 80)
num.Nos <- length(Nmice)
power.vec <- rep(0,num.Nos)
Strong <- "The evidence for separate death rates is conclusive"
Weak <- "Separate death rates is a better model but weakly so."


for(j in 1:num.Nos){
  success.counter <- rep(0, Nsims)
  jth.no <- Nmice[j]
  for(i in 1:Nsims){
    Wdatasim <- dataset.sim(mles=Wmles[,2], no=jth.no, nreps=3, id.reps=1:3, days=0:72, 
                            Strain.name="W", Dose=100000)
    Mdatasim <- dataset.sim(mles=Mmles[,2], no=jth.no, nreps=3, id.reps=4:6,days=0:72, 
                            Strain.name="M", Dose=100000)
    sim.diff.test <- joint.vs.sep.fit(dataframe1=Wdatasim, dataframe2=Mdatasim,Nreps=TRUE)
    A <- sim.diff.test$Best.model
    success.counter[i] <- (A==Strong)||(A==Weak)		
  }
  power.vec[j] <- sum(success.counter)/Nsims 
}


cbind(Nmice,power.vec)

save.image("GalleriaCalcs.RData")


#### 100 sims
#> cbind(Nmice,power.vec)
#      Nmice power.vec
#[1,]     5      0.04
#[2,]    10      0.03
#[3,]    20      0.16
#[4,]    40      0.64
#[5,]    80      0.97

### 1000 sims
#> cbind(Nmice,power.vec)
#Nmice power.vec
#[1,]     5     0.012
#[2,]    10     0.048
#[3,]    20     0.207
#[4,]    40     0.611
#[5,]    80     0.981

###### Galleria test:
# > galleria.test
# $data1CIs.mat # W MLES
# 2.5%         MLE       97.5%
#   Beta0  0.00928471  0.01113821  0.01299171
# Beta1 -0.32579254 -0.30445290 -0.28311325
# 
# $data2CIs.mat # M MLES
# 2.5%          MLE         97.5%
#   Beta0 -0.00423441 -0.002506171 -0.0007779315
# Beta1 -0.25094424 -0.231046786 -0.2111493266
# 
# $data3CIs.mat
# 2.5%         MLE       97.5%
#   Beta0  0.01728562  0.01854582  0.01980601
# Beta1 -0.28816726 -0.27365843 -0.25914959
# 
# $data1.nll
# [1] 41.3517
# 
# $data2.nll
# [1] 24.04336
# 
# $data3.nll
# [1] 76.92255
# 
# $BIC.data1
# [1] 87.09784
# 
# $BIC.data2
# [1] 52.48117
# 
# $BIC.sep  ###  This is the best model
# [1] 142.3516
# 
# $BIC.joint ### This is the bic saying there's no difference
# [1] 159.6259
# 
# $Best.model
# [1] "The evidence for separate death rates is conclusive"


######  May 11th meeting:

####### And plot simulations side by side these estimations with the two sets of parameters:

####  Simulating the exact time of death
no<- 5
Wsim <- deaths.sims(mles=Wmles[,2], no=no,dose=100000)
Msim <- deaths.sims(mles=Mmles[,2], no=no, dose=100000)

####  Simulating the observations every day, at the same time
Wsim.samps <- Sampling.DP(AllDeaths=Wsim$Nalive, cumtimes=Wsim$cumtimes, len=72,tau=1)
Msim.samps <- Sampling.DP(AllDeaths=Msim$Nalive, cumtimes=Msim$cumtimes, len=72,tau=1)


##### Figure S8b  #######
svg(filename="FigureS8b.svg", width=8, height=8)
par(oma=c(2,2,2,2))
plot(Wsim.samps$times,Wsim.samps$obs.ts, col="blue",type="s", lwd=2,ylim=c(0,no), 
     ylab="Number of larvae alive", xlab="Hours", cex.lab=1.5, main="Simulated deaths,10000 CFU replicates", xlim=c(0,72))
points(Wsim.samps$times,Msim.samps$obs.ts, col="red",type="s", lwd=2)
legend("bottomleft", c("Wild Type", "Mutant"), lty=c(1,1), lwd=c(2,2), col=c("blue","red"), cex=1.25, bty="n")
dev.off()


####### Figure 5d #########
Wpreds <- Expected.pred(mlesmat=Wmles, no=80,max.days=72,dose=100000)
Mpreds <- Expected.pred(mlesmat=Mmles, no=80,max.days=72,dose=100000)

no<- 80
tvec <- Wpreds[,1]
Wpred.hat <- Wpreds[,3]
Mpred.hat <- Mpreds[,3]

svg(filename="Figure5D.svg", width=8, height=8)
par(oma=c(2,2,2,2))
plot(tvec, Wpred.hat, type="l", col="blue", lwd=2, ylim=c(0,no), ylab="Number of larvae alive", xlab="Hours", 
     cex.lab=1.5, main="Predicted deaths,100000 CFU replicates", xlim=c(0,72))
points(tvec, Wpreds[,2], type="l", col="blue",lty=2, lwd=1)
points(tvec, Wpreds[,4], type="l", col="blue",lty=2, lwd=1)
points(tvec, Mpred.hat, type="l", col="red", lwd=2)
points(tvec, Mpreds[,2], type="l", col="red",lty=2, lwd=1)
points(tvec, Mpreds[,4], type="l", col="red",lty=2, lwd=1)
legend("topright", c("Wild Type", "Mutant"), lty=c(1,1), lwd=c(2,2), col=c("blue","red"), cex=0.9, bty="n")
dev.off()

###########  Figure 5 e now:

#Nmice power.vec
#[1,]     5     0.012
#[2,]    10     0.048
#[3,]    20     0.207
#[4,]    40     0.611
#[5,]    80     0.981

svg(filename="Figure5e.svg", width=8, height=8)
N.start <- c(5,10,20,40,80)
mice.power <-c(0.06, 0.13, 0.57, 0.92, 1.00)
larvae.power <- c(0.012,0.048,0.207,0.611,0.981)
par(oma=c(2,2,2,2))
plot(N.start, mice.power, pch=19, xlab="Number of individuals per batch (3 replicates)", 
     ylab="Pr(Correctly detecting a faster mutant death rate)", cex.lab=1.05,
     bty="l", type="b", col="black", lwd=2, ylim=c(0,1))
points(N.start, larvae.power, pch=1, type="b", col="black", lwd=2)
legend("bottomright", legend=c("Mice", "Larvae"), pch=c(19,1), lwd=c(2,2), bty="n" )
dev.off()


