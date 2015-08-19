library(ape)
library(splines)
library(gplots)
library(RColorBrewer)

##internal functions, just read these in###########################
###########################
###########################
###########################
###Need this to get individual point at a time of profile
#this function is based on Townsend 2007
site.summer<-function(rate.vector,time)
{
	length(rate.vector)->calculation.length
		at.site<-matrix(nrow=calculation.length)
for(i in 1:calculation.length)
		

	{
		rate.vector[i]->current
		16*current*current*time*exp(-4*current*time)->at.site[i]
		}
		sum(at.site)->inform.at.time
		return(inform.at.time)
		}

##another internal function
###########################
###########################
###########################

get.ind.sites<-function(rate.output,breaks)
{
	
	rate.output->rates
	length(rates)->vector.length
	c(1:vector.length)->numbers
	cbind(numbers,rates)->unsorted.matrix
	length(breaks[,1])->n
	length(rates)->limit
	matrix(ncol=n, nrow=limit)->extracted.sites
	matrix(ncol=n)->names.of.columns
	for(i in 1:n)
	{
###this part looks through the breaks and extracts the site numbers for each user specified bin
		breaks[i,]->upper.lower
		upper.lower[1]->lower
		upper.lower[2]->upper
		which(rates>=lower)->lista
		which(rates<=upper)->listb
		####get the list of sites, which are bigger than lower bound but smaller than upper bound
		lista[(lista%in%listb)]->numbers
		length(numbers)->data.length
		limit-data.length->filler	
		rep("Na",filler)->fill
		c(numbers,fill)->output
		output-> extracted.sites[,i]
		 
	}
	###assign column names
	for(i in 1:n)
	{
		string1="Charset_"
		string2=paste(string1,i,sep="")
		string3=paste(string2,":",sep="")
		names.of.columns[,i]<-string3
		}
		colnames(extracted.sites)<-names.of.columns
		as.data.frame(extracted.sites)->ES
		return(ES)

}	

###internal function #3
###########################
###########################
###########################
###########################

inform.profile.generator2<-function(use.rates,tree)
	{
	branching.times(tree)->btimes
	c(0,btimes)->btimes2
	sort(btimes2)->sorted.btimes
	length(btimes2)->branching.points
	length(use.rates)->calculation.length
	
	inform.at.time<-matrix(ncol=branching.points)
	for(i in 1:branching.points)
	{
		sorted.btimes[i]->btime
		site.summer(use.rates,btime)->inform.at.time[i]

	
}
inform.at.time->close
return(close)
}

defined.multi.profile<-function(rate.vector,tree,breaks)
{
	
	length(rate.vector)->n
	branching.times(tree)->btimes
	c(0,btimes)->btimes2
	sort(btimes2)->sorted.btimes
		length(btimes2)->branching.points

	
	length(breaks[,1])->n.parts
	close<-matrix(ncol=branching.points,nrow=n.parts)
	for (i in 1:n.parts)
	{
	
	min(breaks[i,]):max(breaks[i,])->part
	as.matrix(part)->partx
	partx[partx%in%1:n]->part.check
	as.numeric(part.check)->part.check
	rate.vector->rates
	rates[part.check]->part.current
	inform.profile.generator2(part.current,tree)->close[i,]

	}

	

rbind(sorted.btimes,close)->closer
return(closer)
}


Approximator<-function(t,t0,rateVector,s)	
{	
rateVector->rv
	currentProbability<-matrix(nrow=length(rv), ncol=1)
	Expectationxinnersum1<-c(0)
	Expectationxinnersum2<-c(0)
	Expectationy<-c(0)
	Expectationy2<-c(0)
	ExpectationX1Y<-c(0)
	ExpectationSQROOTX1Y<-c(0)
	length(rv)->n

###Loop calculations and variance
for(i in 1:n)
{
	rv[i,]->rateVector2
	npnl<-pnl(rateVector2,t,t0,s)
	npro<-prother(rateVector2,t,t0,s)
	npsnr<-psnr(rateVector2,t,t0,s)
	
	Expectationy<-Expectationy+npsnr
	Expectationxinnersum1<-Expectationxinnersum1+npnl
	Expectationxinnersum2<-Expectationxinnersum2+npnl*npnl
	Expectationy2<-Expectationy2+ npsnr* npsnr
	ExpectationX1Y<-ExpectationX1Y+ npsnr* npnl
	ExpectationSQROOTX1Y<-ExpectationSQROOTX1Y+ npsnr*sqrt(npnl)
}


Expectationx<-Expectationxinnersum1+sqrt((Expectationxinnersum1/pi))


Expectation<- Expectationy-Expectationx
variancey<- Expectationy-Expectationy2
variancex<-((pi-1)/pi)*Expectationxinnersum1-Expectationxinnersum2
variance<-variancey+variancex-2*ExpectationX1Y-(2/sqrt(pi)) * ExpectationSQROOTX1Y	
	
rnorm(n, mean=Expectation, sd=sqrt(variance))->ndistr
princtree<-pnorm(-0.5,mean=Expectation, sd=sqrt(variance))
prpolytomy<-pnorm(0.5,mean=Expectation, sd=sqrt(variance))-pnorm(-0.5,mean=Expectation, sd=sqrt(variance))
prcortree=1-pnorm(0.5,mean=Expectation, sd=sqrt(variance))
c("Probabilty Correct", "Probability Polytomy", "Probability Incorrect" )->labels
c(prcortree,prpolytomy,princtree)->values
labels->names(values)
return(values)
}

psnr<-function(lambda,t,t0,s)
{
	(-1/s^3 + 1/s^2 + (3/s^3 -1/s^2 -2/s +1 + (-4/s^2 + 4/s -1)*exp(-t0*lambda)) *exp(-(4*s)/(s-1)*t*lambda) + (-8/s^3+4/s^2  +(8/s^2-4/s) *exp (-t0*lambda)) *exp((-3*s)/(s-1)*t*lambda) +(6/s^3 -4/s^2 +2/s -4/s^2 *exp(-t0*lambda) )*exp((-2*s)/(s-1)*t*lambda) )->psnr.value
	return(psnr.value)}
	
	
pnl<-function(lambda,t,t0,s)
{pnl.value<-( -1/s^3 +1/s^2 + (3/s^3-1/s^2 + (-4/s^2+2/s) * exp(-t0*lambda) ) *exp(((-4*s)/(s-1))*t*lambda) + (-8/s^3+4/s^2+ (8/s^2-4/s) * exp(-t0*lambda)) * exp((-3*s)/(s-1)*t*lambda) + (6/s^3 - 4/s^2 + (-4/s^2+2/s) * exp(-t0*lambda))*exp((-2*s/(s-1))*t*lambda) )
	return(pnl.value)
	}

pnL2<-function(lambda,t,t0,s)
{pnl.value<-(-1/s^3+1/s^2+ (3/s^3-1/s^2+ (-4/s^2+2/s) *exp(-t0*lambda) ) *exp(((-4*s)/(s-1))*t*lambda) + (-8/s^3+4/s^2+(8/s^2-4/s)*exp(-t0*lambda)) * exp((-3*s)/(s-1)*t*lambda) +(6/s^3 - 4/s^2 + (-4/s^2+2/s) *exp(-t0*lambda))*exp((-2*s/(s-1))*t*lambda) )
	return(pnl.value)
	}

	prother<-function(lambda,t,t0,s)
	{prother.value<-1-pnL2(lambda,t,t0,s)-pnl(lambda,t,t0,s)-psnr(lambda,t,t0,s)
		return(prother.value)}
		

Approximator.lite<-function(t,t0,rateVector,s)	
{	
rateVector->rv
	currentProbability<-matrix(nrow=length(rv), ncol=1)
	Expectationxinnersum1<-c(0)
	Expectationxinnersum2<-c(0)
	Expectationy<-c(0)
	Expectationy2<-c(0)
	ExpectationX1Y<-c(0)
	ExpectationSQROOTX1Y<-c(0)
	length(rv)->n

###Loop calculations and variance
for(i in 1:n)
{
	rv[i,]->rateVector2
	npnl<-pnl(rateVector2,t,t0,s)
	npro<-prother(rateVector2,t,t0,s)
	npsnr<-psnr(rateVector2,t,t0,s)
	
	Expectationy<-Expectationy+npsnr
	Expectationxinnersum1<-Expectationxinnersum1+npnl
	Expectationxinnersum2<-Expectationxinnersum2+npnl*npnl
	Expectationy2<-Expectationy2+ npsnr* npsnr
	ExpectationX1Y<-ExpectationX1Y+ npsnr* npnl
	ExpectationSQROOTX1Y<-ExpectationSQROOTX1Y+ npsnr*sqrt(npnl)
}


Expectationx<-Expectationxinnersum1+sqrt((Expectationxinnersum1/pi))


Expectation<- Expectationy-Expectationx
variancey<- Expectationy-Expectationy2
variancex<-((pi-1)/pi)*Expectationxinnersum1-Expectationxinnersum2
variance<-variancey+variancex-2*ExpectationX1Y-(2/sqrt(pi)) * ExpectationSQROOTX1Y	
	
rnorm(n, mean=Expectation, sd=sqrt(variance))->ndistr
princtree<-pnorm(-0.5,mean=Expectation, sd=sqrt(variance))
prpolytomy<-pnorm(0.5,mean=Expectation, sd=sqrt(variance))-pnorm(-0.5,mean=Expectation, sd=sqrt(variance))
prcortree=1-pnorm(0.5,mean=Expectation, sd=sqrt(variance))
c("Probabilty Correct", "Probability Polytomy", "Probability Incorrect" )->labels
c(prcortree,prpolytomy,princtree)->values
labels->names(values)
return(prcortree)
}

space.maker<-function(rateVector,t,s)
{
	t/20->by.this
	seq(by.this,t-0.0001,by=by.this)->lilts
	rowspace<-matrix(nrow=1,ncol=length(lilts))
	for (i in 1:length(lilts))
	{
	lilts[i]->to
	Approximator.lite(t,to,rateVector,s)->rowspace[i]
	}
	return(rowspace)
}

space.maker.narrow<-function(rateVector,t,s)
{
	t/2->halft
	halft/20->by.this
	seq(by.this, halft-0.0001,by=by.this)->lilts
	rowspace<-matrix(nrow=1,ncol=length(lilts))
	for (i in 1:length(lilts))
	{
	lilts[i]->to
	Approximator.lite(t,to,rateVector,s)->rowspace[i]
	}
	return(rowspace)
}


