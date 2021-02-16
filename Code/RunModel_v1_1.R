### Bring in data

### If want to run Province level analysis
output.prob.gDist.byGen<-readRDS(file = "output.prob.gDist.byGen_PROV.rds")
probSampSer<-readRDS(file = "probSampSer_PROV.rds")
output.cellindex<-readRDS(file = "output.cellindex_PROV.rds")
output.meanGens<-readRDS(file = "output.meanGens_PROV.rds")
outputTipDates<-readRDS(file = "outputTipDates_PROV.rds")
TimeMRCA<-readRDS(file = "TimeMRCA_PROV.rds")
extTranMatDat<-readRDS(file = "extTranMatDat_PROV.rds")
tranmat.cdr<-readRDS(file = "tranmat.cdr_PROV.rds")

#### If want to run bangkok data
output.prob.gDist.byGen<-readRDS(file = "output.prob.gDist.byGen_BKK.rds")
probSampSer<-readRDS(file = "probSampSer_BKK.rds")
output.cellindex<-readRDS(file = "output.cellindex_BKK.rds")
output.meanGens<-readRDS(file = "output.meanGens_BKK.rds")
outputTipDates<-readRDS(file = "outputTipDates_BKK.rds")
TimeMRCA<-readRDS(file = "TimeMRCA_BKK.rds")
extTranMatDat<-readRDS(file = "extTranMatDat_BKK.rds")
tranmat.cdr<-readRDS(file = "tranmat.cdr_BKK.rds")


### Functions
source("LikelihoodFunction.R")
sourceCpp("MatrixMultiplication.cpp")

### Set up for parallel processing
ncore=4
registerDoMC(cores=ncore)

extTranMatDat$pars$foi<-log(0.04/(1-0.04)) ### Assume FOI of 0.04 per serotype
maxTranGens<-25 ### maximum number of tranmission generations
starting.vals=c(
	-3, 	# mosquito density
	-0.05, 	# Recent homotypic immuity
	0.2, 	# Incidence over all time
	-2, 	# Increased probability of being home - sus
	-2, 	# Increased probability of being home - infected
	0.8 	# Recent Heterotypic immunity
	)

findLLfun<-function(pars,negLL=T){

	extTranMatDat$pars$mosq<-pars[1]
	extTranMatDat$pars$hom<-pars[2]
	extTranMatDat$pars$allInc<-pars[3]
	extTranMatDat$pars$homeSus<-pars[4]
	extTranMatDat$pars$homeSick<-pars[5]
	extTranMatDat$pars$het<-pars[6]

	LL<-rep(NaN,4)
	for (i in 1:4){
		LL[i]<-LikFunc(
			dat.in=dat.tmp.allser[[i]],
			extProbSampleDat.tmp=probSampSer.bs[[i]],
			probByGen.tmp=output.prob.gDist.byGen.bs[[i]],
			maxGen=maxTranGens,
			extTranMatDat.tmp=extTranMatDat,
			cdr.mat=tranmat.cdr,
			ncore=ncore,
			serotype=i)
	}
	LL<-sum(LL)

	print(paste0("pars: ",c(pars)))
	print(paste0("LL: ",LL))
	if(LL==0)return(NA)
	if(negLL){LL<--LL}
	
	return(LL)
}


optFunc<-function(start.pars){
	est<-ucminf(start.pars,fn=findLLfun,hessian=0,control=list(xtol=0.00001,grtol=0.00001))
	return(est)
}

### Bootstrap 
output.prob.gDist.byGen.bs<-output.prob.gDist.byGen
probSampSer.bs<-probSampSer
output.cellindex.bs<-output.cellindex
output.meanGens.bs<-output.meanGens
outputTipDates.bs<-outputTipDates
TimeMRCA.bs<-TimeMRCA

### Prepare all pairs of sequences
all.seqs<-matrix(NaN,0,2)
for (ser in 1:4){
	all.seqs<-rbind(all.seqs,cbind(which(output.cellindex[[ser]]>=0),ser))
}

### Run the bootstrap
nboot<-100
npars=6
out.par<-matrix(NaN,npars,nboot)
out.ll<-rep(NaN,nboot)
for (i in 1:nboot){
	all.seqs.bs<-all.seqs[sample(nrow(all.seqs),replace=T),]
	if(i==1){all.seqs.bs<-all.seqs[1:nrow(all.seqs),]}

	tmp<-NULL
	for (ser in 1:4){
		samp<-all.seqs.bs[,1][which(all.seqs.bs[,2]==ser)]
		output.prob.gDist.byGen.bs[[ser]]<-output.prob.gDist.byGen[[ser]][samp,samp,]
		probSampSer.bs[[ser]]<-probSampSer[[ser]][samp,]
		output.cellindex.bs[[ser]]<-output.cellindex[[ser]][samp]
		outputTipDates.bs[[ser]]<-outputTipDates[[ser]][samp]
		TimeMRCA.bs[[ser]]<-TimeMRCA[[ser]][samp,samp]
		output.meanGens.bs[[ser]]<-output.meanGens[[ser]][samp,samp]
	}

	dat.tmp.allser<-list()	
	for (ii in 1:4){
		dat.tmp.allser[[ii]]<-cbind(expand.grid(1:length(output.cellindex.bs[[ii]]),1:length(output.cellindex.bs[[ii]])),expand.grid(output.cellindex.bs[[ii]],output.cellindex.bs[[ii]]),as.numeric(output.meanGens.bs[[ii]]))
		colnames(dat.tmp.allser[[ii]])<-c("id1","id2","loc1","loc2","meanGen")
		dat.tmp.allser[[ii]]$y2<-1
		dat.tmp.allser[[ii]]<-cbind(dat.tmp.allser[[ii]],TimeStart=as.numeric(TimeMRCA.bs[[ii]]))

		a<-which(dat.tmp.allser[[ii]][,1]>dat.tmp.allser[[ii]][,2]&dat.tmp.allser[[ii]]$TimeStart>0&dat.tmp.allser[[ii]][,5]<=(2*maxTranGens)&is.na(rowSums(dat.tmp.allser[[ii]][,1:4]))==F)
		dat.tmp.allser[[ii]]<-dat.tmp.allser[[ii]][a,]

	}

	mle<-optFunc(start.pars=starting.vals)
	out.par[,i]<-mle$par
	out.ll[i]<-mle$value

	print(paste0("Interation: ",i,", mle: ",mle$value))
}


