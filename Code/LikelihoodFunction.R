##### Function that calculates likelihood

### Libraries
library(nloptr)
library(ucminf)
library(doMC)
library(Rcpp)
library(RcppEigen)
library(Rfast)
library(abind)



### Input data:
### 1. dat.in: five column dataframe: id1 id of tip 1, id2 id of tip2, loc1: Location of tip1, loc2: location of tip 2, TimeStart - time of MRCA

### 2. extTranMatDat.tmp. List with the following data:
#### pars: a vector of the parameter values
#### popByCell: a vector with the size of the population within each location
#### allcases: A vector with the total number of cases over all time in each location
#### Mosq: A vector with the mosquito density in each location
#### Popden: A vector with the population density in each location
#### RelSeroIncPriorYear: An array that has the historic incidence for the two year period going back in time for each time point per serotype and per location

#### 3. extProbSampleDat.tmp. A matrix which has the probability of sampling in each location at each time point

#### 4. probByGen.tmp. The probability that each pair of sequences is separated by a specified number of transmission events. An array with dimensions number of sequences, number of sequences, 200/the 200 represents the maximum number of transmission events.



LikFunc<-function(
	dat.in, # Input data on pairs of sequences
	extProbSampleDat.tmp,
	probByGen.tmp,
	maxGen,  	# Maximum number of generations
	extTranMatDat.tmp,
	cdr.mat, # Input mobility matrix for adults
	serotype=1, # Input serotype
	ncore=1, ### Number of cores available
	Eigencpp=TRUE){
	
	npairs=nrow(dat.in) ## Number of pairs of sequences
	if(npairs==0)return(0)
	nlocs=length(extTranMatDat.tmp$popbyCell) ## Number of locations
	
	lambda<-exp(extTranMatDat.tmp$pars$foi)/(1+exp(extTranMatDat.tmp$pars$foi))	## Lambda is provided on a logit scale
	tmp.pHome<-extTranMatDat.tmp$popbyCell*(extTranMatDat.tmp$allCases^extTranMatDat.tmp$pars$allInc)

	### The probability of the location of the MRCA is linked to population size
	mrcaVec<-extTranMatDat.tmp$popbyCell
	mrcaVec<-mrcaVec/sum(mrcaVec)

	### Create mobility matrix for susceptible individuals
	tmpbase<-cdr.mat
	tmppar<-exp(extTranMatDat.tmp$pars$homeSus)/(1+exp(extTranMatDat.tmp$pars$homeSus))
	tmpdiag<-diag(tmpbase)+tmppar
	tmpdiag[which(tmpdiag>0.999)]<-0.999
	diag(tmpbase)<-0
	tmpbase<-sweep(tmpbase,1,rowSums(tmpbase),"/")
	tmpbase<-sweep(tmpbase,1,(1-tmpdiag)/(1-diag(tmpbase)),"*")
	diag(tmpbase)<-tmpdiag
	mobility.sus<-tmpbase

	### Create mobility matrix for infected individuals
	tmpbase.sick<-cdr.mat
	tmppar<-exp(extTranMatDat.tmp$pars$homeSick)/(1+exp(extTranMatDat.tmp$pars$homeSick))
	tmpdiag<-diag(tmpbase.sick)+tmppar
	tmpdiag[which(tmpdiag>0.999)]<-0.999
	diag(tmpbase.sick)<-0
	tmpbase.sick<-sweep(tmpbase.sick,1,rowSums(tmpbase.sick),"/")
	tmpbase.sick<-sweep(tmpbase.sick,1,(1-tmpdiag)/(1-diag(tmpbase.sick)),"*")
	diag(tmpbase.sick)<-tmpdiag
	mobility.sick<-tmpbase.sick

	### Reltaive incidence for homotypic and heterotypic
	RelIncDat<-extTranMatDat.tmp$RelSeroIncPriorYear[,,serotype]
	RelIncDatHetSer<-apply(extTranMatDat.tmp$RelSeroIncPriorYear[,,-serotype],c(1,2),mean)

	### Contribution to likelihood for each pair
	llIndPair<-function(ii){

		if(ii>npairs){return(NA)}
					
		ind.series<-dat.in$TimeStart[ii]
		ind.series<-ind.series:(ind.series+maxGen-1)

		tmp.recent.ind<-RelIncDat[,ind.series[1]]
		tmp.recent.ind.het<-RelIncDatHetSer[,ind.series[1]]
		
		nGen.in.time.series<-length(ind.series[which(ind.series<=ncol(RelIncDat))])

		probInfec.tmp<-extTranMatDat.tmp$mosq^extTranMatDat.tmp$pars$mosq
		tmp.adj<-sweep(mobility.sus,2,probInfec.tmp,"*")
			
		tmp.pHome.tvar<-extTranMatDat.tmp$popbyCell*(extTranMatDat.tmp$allCases^extTranMatDat.tmp$pars$allInc)*(tmp.recent.ind)^extTranMatDat.tmp$pars$hom*(tmp.recent.ind.het^extTranMatDat.tmp$pars$het)		
		move3<-tcrossprod(mobility.sick,tmp.adj)
		move4<-sweep(move3,2,tmp.pHome.tvar,"*")
		TranMat.tmp<-sweep(move4,1,rowSums(move4),"/")

		probByGenA<-probByGen.tmp[dat.in[ii,1],dat.in[ii,2],1:maxGen]
		probByGenB<-probByGen.tmp[dat.in[ii,2],dat.in[ii,1],1:maxGen]

		gensA<-(1:maxGen)[which(probByGenA[1:nGen.in.time.series]>0)]
		gensB<-(1:maxGen)[which(probByGenB[1:nGen.in.time.series]>0)]

		if((length(gensA)<2|length(gensB)<2)){return(log(1e-16))}
		maxgen.tmp<-min(max(c(gensA,gensB)),nGen.in.time.series)

		TranMatArray<-array(NA,c(nlocs,nlocs,maxgen.tmp))
		TranMatArray[,,1]<-TranMat.tmp
		for (j in 2:maxgen.tmp){
			tmp.recent.ind<-RelIncDat[,ind.series[j]]
			tmp.recent.ind.het<-RelIncDatHetSer[,ind.series[j]]
			tmp.pHome.tvar<-extTranMatDat.tmp$popbyCell*(extTranMatDat.tmp$allCases^extTranMatDat.tmp$pars$allInc)*(tmp.recent.ind)^extTranMatDat.tmp$pars$hom*(tmp.recent.ind.het^extTranMatDat.tmp$pars$het)
			move4<-sweep(move3,2,tmp.pHome.tvar,"*")
			TranMat.tmp<-sweep(move4,1,rowSums(move4),"/")

			if(Eigencpp){TranMatArray[,,j]<-eigenMapMatMult(TranMatArray[,,j-1], TranMat.tmp)}else{
			TranMatArray[,,j]<-TranMatArray[,,j-1]%*%TranMat.tmp}
		}
		
		TranMatArrayA<-TranMatArray[,,gensA]
		TranMatArrayA<-sweep(TranMatArrayA,3,probByGenA[gensA],"*")
		TranMatArrayA<-sweep(TranMatArrayA,2,extProbSampleDat.tmp[dat.in[ii,1],],"*")

		TranMatArrayB<-TranMatArray[,,gensB]
		TranMatArrayB<-sweep(TranMatArrayB,3,probByGenB[gensB],"*")
		TranMatArrayB<-sweep(TranMatArrayB,2,extProbSampleDat.tmp[dat.in[ii,2],],"*")

		TranMatArrayA2<-sweep(TranMatArrayA,1,mrcaVec,"*")
		TranMatArrayA3<-matrix(TranMatArrayA2,nlocs,nlocs*length(gensA))
		TranMatArrayB2<-matrix(TranMatArrayB,nlocs*length(gensB),nlocs,byrow=T)
		probAllPrs<-eigenMapMatMult(TranMatArrayB2,TranMatArrayA3)
		den<-sum(probAllPrs)

		TranMatArrayA3<-matrix(TranMatArrayA2[,dat.in[ii,3],],nlocs,length(gensA))
		TranMatArrayB2<-matrix(TranMatArrayB[,dat.in[ii,4],],length(gensB),nlocs,byrow=T)
		num<-sum(eigenMapMatMult(TranMatArrayB2,TranMatArrayA3))

		lik=log(num/den)

		return(lik)
	}

	### Run in parallel
	n.per.core=ceiling(npairs/ncore)
	ntot=n.per.core*ncore
	lines.mat<-matrix(1:ntot,nc=ncore,nr=n.per.core)
	func<-function(lines){
		out<-rep(NaN,length(lines))
		for (iii in 1:length(lines)){out[iii]<-llIndPair(ii=lines[iii])}
		return(out)
	}

	if(ncore>1){
		aa<-foreach (jj=1:ncore)%dopar%func(lines=lines.mat[,jj])}else{
		aa<-list(func(lines=lines.mat[,1]))
	}


	out<-sum(do.call("c",aa),na.rm=T)

	return(out)
}
