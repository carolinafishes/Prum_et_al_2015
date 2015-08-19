identifyGoodSeqs<<-function(project,nInds,nLoci,nHomologs,minReadsMapped){
	nMapped<<-0
	for(i in 1:nHomologs){
		blah<<-read.table(paste("../Results/",project,"_AssemblySummary_nMappedReads.txt",collapse="",sep=""),skip=2+(i-1)*(nInds+3),nrow=nInds)
		nMapped<<-c(nMapped,as.numeric(t(as.matrix(blah[,2:(ncol(blah))]))))
	}
	hist(log10(nMapped),nclass=200)
	lines(c(log10(minReadsMapped),log10(minReadsMapped)),c(0,100000),col="red",lwd=2)
	firstWrite=T
	for(i in 1:nHomologs){
		for(j in 1:nInds){
			for(k in 1:nLoci){
				if(nMapped[1+(i-1)*(nInds)*(nLoci)+(j-1)*(nLoci)+k]>=minReadsMapped){
					write(paste(i,"\t",j,"\t",k,sep="",collapse=""),file="../Results/GoodConSeqIDs.txt",append=!firstWrite)
					firstWrite=F
				}
			}
		}
	}

}
