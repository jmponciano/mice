#library("extraDistr")
library("MASS")


negll.death<- function(guess, dataframe,nreps=FALSE){
	
	# Assumes labeled columns in data mat as follows: Day, Nexperiment, Strain, Nsurviv, Dose
	beta0 <- guess[1]
	beta1 <- guess[2]
	
	exper.id <- dataframe$Nexperiment
	unique.ids <- unique(exper.id)
	n.ids <- length(unique.ids)
	
	if(nreps==FALSE){
		unique.doses <- log(unique(dataframe$Dose))
		n.doses <- length(unique.doses)
		if(n.ids!=n.doses){print("STOP!!! N doses != N experiments!!")}
	}else{
		
		unique.doses <- rep(log(unique(dataframe$Dose)), n.ids)
		n.doses <- length(unique.doses)
	}
	
	nlls.perts <- rep(0,n.ids)  #vector of nll, one per experimental time series
	for(i in 1:n.ids){
		
		ith.id <- unique.ids[i]
		subdat <- dataframe[dataframe$Nexperiment==ith.id,]
		ntsteps<- nrow(subdat)
		nim1s  <- subdat$Nsurviv[1:(ntsteps-1)]
		nis    <- subdat$Nsurviv[2:ntsteps]
		#losses <- nim1s-nis
		Taus   <- subdat$Day[2:ntsteps]-subdat$Day[1:(ntsteps-1)]
		
		ith.dose   <- unique.doses[i]
		log.mu <- beta0 + beta1*ith.dose;
		mu     <- exp(log.mu)
		pvec   <- exp(-mu*Taus)
		pvec[pvec==0] <- .Machine$double.xmin
		
		binom.llikes <- dbinom(x=nis, size=nim1s, prob=pvec, log=TRUE)
		#print(binom.llikes)
		nlls.perts[i] <- sum(binom.llikes)
		
	}
	
	negll <- -sum(nlls.perts)
	return(negll)
} 

negll.joint.death<- function(guess, dataframe,nreps=FALSE){
	
	# Assumes labeled columns in data mat as follows: Day, Nexperiment, Strain, Nsurviv, Dose
	beta0 <- guess[1]
	beta1 <- guess[2]
	
	exper.id <- dataframe$Nexperiment
	unique.ids <- unique(exper.id)
	n.ids <- length(unique.ids)
	
	if(nreps==FALSE){
		unique.doses <- rep(log(unique(dataframe$Dose)),2) # assumes data entered by dose in same order!!
		n.doses <- length(unique.doses)
		if(n.ids!=n.doses){print("STOP!!! N doses != N experiments!!")}
	}else{
		unique.doses <- rep(log(unique(dataframe$Dose)), n.ids)
		n.doses <- length(unique.doses)
	}
	
	nlls.perts <- rep(0,n.ids)  #vector of nll, one per experimental time series
	for(i in 1:n.ids){
		
		ith.id <- unique.ids[i]
		subdat <- dataframe[dataframe$Nexperiment==ith.id,]
		ntsteps<- nrow(subdat)
		nim1s  <- subdat$Nsurviv[1:(ntsteps-1)]
		nis    <- subdat$Nsurviv[2:ntsteps]
		#losses <- nim1s-nis
		Taus   <- subdat$Day[2:ntsteps]-subdat$Day[1:(ntsteps-1)]
		
		ith.dose   <- unique.doses[i]
		log.mu <- beta0 + beta1*ith.dose;
		mu     <- exp(log.mu)
		pvec   <- exp(-mu*Taus)
		pvec[pvec==0] <- .Machine$double.xmin
		
		binom.llikes <- dbinom(x=nis, size=nim1s, prob=pvec, log=TRUE)
		#print(binom.llikes)
		nlls.perts[i] <- sum(binom.llikes)
		
	}
	
	negll <- -sum(nlls.perts)
	return(negll)
} 



joint.vs.sep.fit <- function(dataframe1, dataframe2, my.guess = c(0.0002,0.000004),my.method="BFGS",
							plot.data1=FALSE, plot.data2=FALSE,Nreps=FALSE){
	

	##### Data 1 fit
	exp.ids1 <- dataframe1$Nexperiment
	unique.ids1 <- unique(exp.ids1)
	n.ids1  <- length(unique.ids1)
	n.transitions <- rep(0,n.ids1)
	for(i in 1:n.ids1){
		ith.id <- unique.ids1[i]
		subdat <- dataframe1[dataframe1$Nexperiment==ith.id,]
		ntsteps<- nrow(subdat)
		n.transitions[i] <- ntsteps-1		
	}
	n1 <- sum(n.transitions)
	# negll.death(guess=c(0.0002,0.000004), dataframe=dataframe1)
	data1.fit <- optim(par=my.guess, fn=negll.death, method=my.method,dataframe=dataframe1,nreps=Nreps, hessian=TRUE)
	data1.mles <- data1.fit$par
	data1.nll  <- data1.fit$val
	my.hess <- data1.fit$hessian
	Fish.Inv <- ginv(my.hess)
	zalphahalf <- qnorm(p=0.975, mean=0, sd=1)
	st.errs <-zalphahalf*sqrt(diag(Fish.Inv))
	low.cis <- data1.mles - st.errs
	hi.cis <- data1.mles + st.errs
	data1CIs.mat <- cbind(low.cis, data1.mles, hi.cis)
	colnames(data1CIs.mat) <- c("2.5%", "MLE", "97.5%")
	row.names(data1CIs.mat) <- c("Beta0","Beta1")
	BIC.data1 <- 2*data1.nll + length(guess)*log(n1)	


	##### Data 2 fit
	exp.ids2 <- dataframe2$Nexperiment
	unique.ids2 <- unique(exp.ids2)
	n.ids2  <- length(unique.ids2)
	n.transitions2 <- rep(0,n.ids2)
	for(i in 1:n.ids2){
		ith.id <- unique.ids2[i]
		subdat <- dataframe2[dataframe2$Nexperiment==ith.id,]
		ntsteps<- nrow(subdat)
		n.transitions2[i] <- ntsteps-1		
	}
	n2 <- sum(n.transitions2)
	# negll.death(guess=c(0.00002,0.00004), dataframe=dataframe2)
	data2.fit <- optim(par=my.guess, fn=negll.death, method=my.method,dataframe=dataframe2,nreps=Nreps , hessian=TRUE)
	data2.mles <- data2.fit$par
	data2.nll  <- data2.fit$val
	my.hess <- data2.fit$hessian
	Fish.Inv <- ginv(my.hess)
	zalphahalf <- qnorm(p=0.975, mean=0, sd=1)
	st.errs <-zalphahalf*sqrt(diag(Fish.Inv))
	low.cis <- data2.mles - st.errs
	hi.cis <- data2.mles + st.errs
	data2CIs.mat <- cbind(low.cis, data2.mles, hi.cis)
	colnames(data2CIs.mat) <- c("2.5%", "MLE", "97.5%")
	row.names(data2CIs.mat) <- c("Beta0","Beta1")
	BIC.data2 <- 2*data2.nll + length(guess)*log(n2)	

	BIC.sep <- 2*(data1.nll+data2.nll) + 2*length(guess)*log(n1+n2)
	
	
	######  Joint data fit
	n3 <- n1+n2
	dataframe3 <- rbind(dataframe1,dataframe2)
	# negll.joint.death(guess=c(0.00002,0.00004), dataframe=dataframe3)
	data3.fit <- optim(par=c(0.0002,0.000004), fn=negll.joint.death, method=my.method,dataframe=dataframe3,nreps=Nreps , hessian=TRUE)
	data3.mles <- data3.fit$par
	data3.nll  <- data3.fit$val
	my.hess <- data3.fit$hessian
	Fish.Inv <- ginv(my.hess)
	zalphahalf <- qnorm(p=0.975, mean=0, sd=1)
	st.errs <-zalphahalf*sqrt(diag(Fish.Inv))
	low.cis <- data3.mles - st.errs
	hi.cis <- data3.mles + st.errs
	data3CIs.mat <- cbind(low.cis, data3.mles, hi.cis)
	colnames(data3CIs.mat) <- c("2.5%", "MLE", "97.5%")
	row.names(data3CIs.mat) <- c("Beta0","Beta1")
	
	BIC.joint <- 2*data3.nll + length(guess)*log(n3)	
	
	Delta.BIC <- BIC.joint-BIC.sep # >0 sep is better; <0 joint is better

	if(Delta.BIC>0){
		
		#sep is better
		if(Delta.BIC>3){Best.model <-"The evidence for separate death rates is conclusive"
			}else{
				Best.model <- "Separate death rates is a better model but weakly so."
			}
		
	}else if(Delta.BIC<0){
		
		#joint is better
		if(Delta.BIC< (-3)){Best.model <-"The evidence for a single death rate is conclusive"
			}else{
				Best.model <- "No difference in death rates is a better model but weakly so."
			}
		
		
		
	}

	
	out <- list(data1CIs.mat=data1CIs.mat,data2CIs.mat=data2CIs.mat,data3CIs.mat=data3CIs.mat, data1.nll=data1.nll,data2.nll=data2.nll, data3.nll=data3.nll, BIC.data1=BIC.data1, BIC.data2 = BIC.data2, BIC.sep=BIC.sep, BIC.joint = BIC.joint, Best.model=Best.model)
	
	return(out)
	
}


Expected.pred <- function(mlesmat,no,max.days,dose){
	
	tvec <- seq(0,max.days, by=0.05)
	b0.1 <- mlesmat[1,1];b0.2 <- mlesmat[1,3];b0.hat <- mlesmat[1,2];
	b1.1 <- mlesmat[2,1];b1.2 <- mlesmat[2,3];b1.hat <- mlesmat[2,2];
	p.ci1 <- exp(-exp(b0.1 + b1.1*log(dose))*tvec)
	p.mle <- exp(-exp(b0.hat + b1.hat*log(dose))*tvec)	
	p.ci2<- exp(-exp(b0.2 + b1.2*log(dose))*tvec)
	
	return(cbind(tvec,no*p.ci1,no*p.mle,no*p.ci2))
	
}

# Gillespie algorithm, death process from infection simulator:

deaths.sims <- function(mles, no,dose){
	
	len <- no+1 # This is the number of death events simulated
			  # In other words, this function simulates the death
			  # process until all animals in the experiment have died.

	b0 <- mles[1];
	b1 <- mles[2];
	log.mu <- b0 + b1*log(dose);
	mu     <- exp(log.mu)
	
	Nalive <- rep(0,len)
	times  <- rep(0,len)
	Nalive[1] <- no
	Uvec   <- runif(n=len)
	
	for(i in 2:len){
		
		N <- Nalive[(i-1)] # current number of alive mice
		death.rate <- mu*N
		
		Ui <- Uvec[(i-1)]
		wait.time <- -(1/death.rate)*log(Ui)
		times[i]  <- wait.time
		Nalive[i] <- N-1
		
	}
	cumtimes <- cumsum(times)
	
	return(list(Nalive=Nalive, times=times, cumtimes=cumtimes))	
	
} 

####  Death-process sampler, at arbitrary days:

# Perfect observations from a simulation of the death process: 
# Assuming that not all are observed, a time series of
# observations from the death process without sampling error is obtained by
# recording the states of the process at a series of times
# t(1) < t(2) <...<t(q).  The time intervals between recordings  
# tau.k = t(k+1)-t(k) can be regular or not.  If they are regular,  
# the likelihood calculations simplify enormously.
# The arguments are: the time series of births and deaths
# The vector of cumulative exponential times until an event occurs
# len is the length of the desired time series
# tau is the desired size of the time interval between observations,
# If tau is not specified, the function automatically determines the
# value of tau from the desired length of the time series and 
# the specified length len
Sampling.DP <- function(AllDeaths,cumtimes, len, start=0, tau = 0){
	
	nevents <- length(cumtimes)
	max.time <- cumtimes[nevents]
	if(tau==0){
		
		time.seq <- seq(start,tvec[nevents], by= (tvec[nevents]-start)/(len-1))		
		}else{
		
		time.ints <- rep(tau,(len-1))
		time.seq <- c(start,start + cumsum(time.ints))
	}
	tau <- time.seq[2]-time.seq[1]
	indices <- rep(0,len)
	for(i in 1:len){
	
		it  <- time.seq[i]
		it.m.cumt <- abs(it-cumtimes)
		pre.ind <- which(it.m.cumt==min(it.m.cumt), arr.ind=TRUE)
		t.diff <- it - cumtimes[pre.ind]
		index <- pre.ind
		if(t.diff<0){index <- pre.ind - 1}
		indices[i] <- index	
	
	}

	out <- list(obs.ts = AllDeaths[indices], indices = indices, times = time.seq,tau=tau) 
	return(out)
}

dataset.sim <- function(mles=Wmles[,2], no=5,nreps=3,id.reps = 1:3,days=0:14,Strain.name="W",Dose=1000){
	
	len <- length(days)
	Day <- rep(days,nreps)
	tau <- days[2]-days[1] #  Assuming all sampling occurs exactly at the same time every day!!!
	Strain <- rep(Strain.name,nreps*len)
	Nsurviv <- rep(0,nreps*len)
	Dose <- rep(Dose,nreps*len)
	Nexperiment <- rep(id.reps, each=len)
	
	for(i in 1:nreps){
		every.event.sim <- deaths.sims(mles=mles, no=no,dose=1000)
		Samp.sim <- Sampling.DP(AllDeaths=every.event.sim$Nalive, cumtimes=every.event.sim$cumtimes,len=len,tau=tau)
		Nsurviv[(len*(i-1) +1):(i*len)] <- Samp.sim$obs.ts
		
	}
	
	out.df <- data.frame(Day=Day, Nexperiment=Nexperiment, Strain=Strain, Nsurviv=Nsurviv, Dose=Dose)
	return(out.df)
}







