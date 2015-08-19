locus<<-1
plotMDS<<-function(project,loc=locus,maxDist=-1,minSep=0.01,xlm=0.5,ylm=0.5,plotPoints=TRUE){
	locus<<-loc+1
print(paste("Determing Orthology for Locus ",loc))
	temp=sqrt(length(scan(paste("../Distances/",project,"_L",loc,".txt",sep="",collapse=""))))
	if(temp<5){
		#write(c(loc),paste("../OrthoSets/",project,"_TaxaInOrthoSets.txt",sep="",collapse=""),append=loc!=1)
		#print(paste(loc,"(has fewer than 5 sequences)"),quote=FALSE)
		#write(rep(0,temp),paste("../OrthoSets/",project,"_L",loc,".txt",sep="",collapse=""),ncol=1)
	}else{
		distances<<-read.table(paste("../Distances/",project,"_L",loc,".txt",sep="",collapse=""))

		blah<<-scan(paste("../Homologs/",project,"_L",loc,".fasta",sep="",collapse=""),what=character(0))
		sequences<<-blah[seq(2,length(blah),2)]
		
		blah<<-blah[seq(1,length(blah),2)]

		#determine which of the sequences come from the reference
		reference<<-rep(FALSE,length(blah))
		for(i in 1:length(blah)){
			#print(i)
			#print(blah[i])
			#print(strsplit(blah[i],"REFERENCE"))
			#print(length(strsplit(blah[i],"REFERENCE")))
			reference[i]=length(strsplit(blah[i],"REFERENCE")[[1]])==2
		}

		blah<<-strsplit(blah,"Copy")
		taxa<<-rep(0,length(blah))
		copy<<-rep(0,length(blah))
		for(i in 1:length(blah)){taxa[i]=blah[[i]][1]}
		for(i in 1:length(blah)){copy[i]=blah[[i]][2]}
		sameTaxa<<-matrix(as.vector(t(outer(taxa,taxa,"=="))),nrow=length(blah))
		sameCopy<<-matrix(as.vector(t(outer(copy,copy,"=="))),nrow=length(blah))
		diagonal<<-matrix(as.vector(t(outer(1:length(blah),1:length(blah),"=="))),nrow=length(blah))

		#scores<<-rep(0,100)
		#for(i in 1:100){
		#	scores[i]=mean(distances[sameTaxa & !diagonal]>i/100)
		#}
		#plot(scores,pch=19,cex=0.5)
		#quartz()
		#clrs=c("red","blue")
		#par(mfrow=c(1,2))
		#plot(distances[!diagonal],ylim=c(0,1),cex=0.5,pch=19,col="blue")
#		points(distances[!diagonal & sameTaxa],cex=0.5,pch=19,col="red")
	#lines(c(0,length(blah)*length(blah)),c(maxDist,maxDist),col="red")

distsToTry=maxDist
if(maxDist<0){
	distsToTry=seq(0.001,1.001,0.001)
}

resultA=rep(FALSE,length(distsToTry))
resultB=rep(FALSE,length(distsToTry))
resultC=rep(FALSE,length(distsToTry))	#number of orthoGroups
resultD=rep(FALSE,length(distsToTry))	#number of taxa contained within a reasonably sized orthogroup

bestMaxDist=-1
bestResultB=5
bestResultD=0
bestOrthoSet=NULL
foundRuleA=FALSE
bestNTaxaInOrthoSet=NULL

allOrthoGroups<<-matrix(0,nrow=length(distsToTry),ncol=length(taxa))

for(d in 1:length(distsToTry)){	

		#automatically determine clusters based on distance threshold, which is chosen to optimize ortholgous sets
		size=nrow(distances)

		#SETUP AN EMPTY MATRIX OF SAME SIZE
		x<<-matrix(rep(0,size*size),nrow=size)
			#testA example x<<-matrix(c(1,0,1,0,0,1,0,1,1,0,1,0,0,1,0,1),nrow=4)
			#testB example x<<-matrix(c(1,1,0,0,1,0,0,1,1,0,0,1,0,0,0,0,1,1,0,1,0,0,0,1,1,0,1,0,1,1,0,0,1,0,0,0,0,1,1,0,1,0,0,0,0,0,0,0,1),nrow=7)

		#SET value equal to zero if above distance threshold, 0 otherwise
		x[distances<distsToTry[d]]=1



		#MAKE A CORRESPONDING MATRIX WITH VALUE EQUAL TO ROW INDEX IF VALUE IN X IS 1
		y1<<-x*(1:size)
	
		#set value equal to max of column if x>0 for that element
		y1prev=y1
		y1<<-x*matrix(rep(pmax(apply(y1,2,max),apply(y1,2,max)),size),nrow=size)
				
		#repeat until all connected samples have same group ID
		counter<<-0
		while(sum(y1prev!=y1)>0){
			counter<<-counter+1
			y1prev=y1
			y1<<-x*matrix(rep(pmax(apply(y1,2,max),apply(y1,2,max)),size),nrow=size)
		}

#if(d==283){write(y1,"test321.txt",ncol=size,sep="\t")}


		#renumber from 1:nClusters
		dy<<-diag(y1)
		uy<<-unique(dy)
		z<<-rep(0,max(dy))
		z[uy]<<-1:length(uy)

		yfinal<<-z[dy]


		yfinal2=yfinal		
		for(i in 1:length(taxa)){			
			if(sum(taxa[yfinal==yfinal[i]]==taxa[i])>1){
				yfinal2[i]=-yfinal2[i]
			}
		}
		allOrthoGroups[d,]=yfinal2
		
		#SUMMARIZE THE CURRENT ORTHO GROUPS AND DETERMINE IF THRESHOLD IS PROPER
		nTaxaInOrtho=rep(0,max(yfinal))
		nUniqueInOrtho=rep(0,max(yfinal))
		nRefsInOrtho=rep(0,max(yfinal))
		for(i in 1:max(yfinal)){
			nTaxaInOrtho[i]=length(taxa[yfinal==i])
			if(nTaxaInOrtho[i]>length(unique(taxa))*0.5 & nTaxaInOrtho[i]<length(unique(taxa))*1.5){
				resultD[d]=resultD[d]+nTaxaInOrtho[i]
			}
			nUniqueInOrtho[i]=length(unique(taxa[yfinal==i]))
			nRefsInOrtho[i]=sum(reference[yfinal==i])
		}
		
		#print(paste(nTaxaInOrtho,nUniqueInOrtho,nRefsInOrtho),sep="\t",collapse="")
		#FIRST RULE: two reference sequences cannot be in the same ortho group
		rule1=length(unique(yfinal[reference]))==sum(reference)

		#SECOND RULE: each reference seq in a ortho group of reasonable size
		rule2=nUniqueInOrtho[yfinal[reference]]>length(unique(taxa))/2

		rule2b=sum(rule2)==sum(reference)

		#THIRD RULE: DON'T OVERLUMP, number of sequences in each orthologus group should be less than 1.5 times the expected number
		score=0
		score2=0
		for(i in 1:max(yfinal)){
			if(nUniqueInOrtho[i]>0){
				if(nTaxaInOrtho[i]/nUniqueInOrtho[i]>score){
					score=nTaxaInOrtho[i]/nUniqueInOrtho[i]
					score2=length(unique(yfinal))
				}
			}
		}

		#OVERALL RULE
		resultA[d]=rule1 && rule2b
		resultB[d]=score
		resultC[d]=score2
		resultD[d]=resultD[d]/sum(nTaxaInOrtho)

		if(!foundRuleA && resultA[d]){
			#first time ruleA met, reset bestResulB to current resultB
			bestResultB=resultB[d]
		}
		foundRuleA = foundRuleA || resultA[d]

		#if((!foundRuleA && resultB[d]<=bestResultB) || (resultA[d] && resultB[d]<=bestResultB)){
		if(resultD[d]>=bestResultD){
			bestMaxDist=distsToTry[d]
			bestResultB=resultB[d]
			bestResultD=resultD[d]
			bestOrthoSet=yfinal
			bestNTaxaInOrthoSet=nTaxaInOrtho
		}
}

#allOrthoGroups[length(distsToTry),]=-1

		write(t(allOrthoGroups),"testing123.txt",ncol=length(taxa),sep="\t")

		#############
		#determine optimal ortho sets given results from all levels:
		#############
		#initialize finalOrthoIDs with the results form the largest d
		
		finalOrthoIDs<<- -1

		#work from greatest distance to least distance
		for(d in length(distsToTry):1){	


			if(max(finalOrthoIDs)<0){
#			if(max(finalOrthoIDs)<=0){


				finalOrthoIDs<<-allOrthoGroups[d,]
				finalOrthoIDs[finalOrthoIDs<0]=0		#disqualified taxa get a ortho group ID of 0 (means not in ortho group yet)				

			}else{

				#identify the IDs for the taxa that have not yet been joined, but are now eligible because they have a positive value at this distance threshold
				candidates=(1:length(taxa))[finalOrthoIDs==0 & allOrthoGroups[d,]>0]
				#compute distances to potential orthoGroup buddies
				avgDists<<-rep(0,length(taxa))
				for(can in candidates){
					if(sum(finalOrthoIDs!=0 & allOrthoGroups[d,]==allOrthoGroups[d,can])>0){ #if to be grouped to taxon already in final list, compute distance
						avgDists[can]=mean(as.numeric(distances[can,finalOrthoIDs!=0 & allOrthoGroups[d,]==allOrthoGroups[d,can]]))	#only include distances to taxa already included that are in the same ortho group for this threshold
					}else{
						avgDists[can]=0		#could not find valid taxa in already in final list, so use distance of 0 (ok to add since it cannot conflict)
					}
				}
				
				#in order of increasing distance, attempt to join candidates to pre-existing orthogroups
				distOrder=order(avgDists)		#this gets the indexes of the smallest to the largest distance
				for(i in distOrder){	
					if(sum(candidates==i)>0){	#only attempt if a candidate
						potentialBuddies<<-finalOrthoIDs>0 & allOrthoGroups[d,]==allOrthoGroups[d,i]
						foundBuddy<<-sum(potentialBuddies)>0
						potentialGroup<<-NULL
						if(foundBuddy){
							potentialGroup<<-unique(finalOrthoIDs[potentialBuddies])
							taxonAlreadyPresent<<-sum( taxa[finalOrthoIDs==potentialGroup]==taxa[i] )>0
						}

						if(length(potentialGroup)>1){
							print(paste("locus ",loc,"Trouble encountered, resetting finalOrthoIDs at d=",d))
							finalOrthoIDs[1:length(finalOrthoIDs)]=0
							break
						}else{

							if(foundBuddy && !taxonAlreadyPresent){	#join if new not needed and the candidate taxon is not already in the orthoGroup being considered (DON'T ALLOW DUPLICATES!)
								#check to make sure that the candidate does not match to two different groups in the final ortho id set
								finalOrthoIDs[i]=potentialGroup
								#print(paste(i,"joined old",finalOrthoIDs[i]," at d=",d))
							}else{
								finalOrthoIDs[i]=max(1,max(finalOrthoIDs)+1)	#since it could not be joined, put in its own group for now (others may join later)
								#print(paste(i,"joined new",finalOrthoIDs[i]," at d=",d))
							}
						}
					}
				}
			}
		}


		#check to see if any redundant sequences in each ortholog, if so, throw BOTH out to be safe
		for(i in unique(finalOrthoIDs)){
			for(j in 1:length(taxa)){
				if(finalOrthoIDs[j]==i){
					if(sum(taxa==taxa[j] & finalOrthoIDs==finalOrthoIDs[j])>1){
						finalOrthoIDs[taxa==taxa[j] & finalOrthoIDs==finalOrthoIDs[j]]=0	#set to 0 for now, will be assigned unique ID in next step
					}
				}
			}
		}

		#put each unassigned sequence into its own orthogroup
		for(i in 1:length(finalOrthoIDs)){
			if(finalOrthoIDs[i]<=0){
				finalOrthoIDs[i]=max(finalOrthoIDs)+1
			}
		}

		#write to file
		write(finalOrthoIDs,paste("../OrthoSets/",project,"_L",loc,".txt",sep="",collapse=""),ncol=1)

		nTaxaInFinalOrtho=rep(0,length(unique(finalOrthoIDs)))
		for(i in 1:max(finalOrthoIDs)){
			if(length(unique(taxa[finalOrthoIDs==i]))>length(taxa[finalOrthoIDs==i])){print(paste("WARNING! locus ",loc," Ortho set ",i," has redudant taxa!"))}
			nTaxaInFinalOrtho[i]=length(taxa[finalOrthoIDs==i])
		}
		write(c(loc,nTaxaInFinalOrtho),paste("../OrthoSets/",project,"_TaxaInOrthoSets.txt",sep="",collapse=""),ncol=length(nTaxaInFinalOrtho)+1,append=loc!=1)
		
		#make a plot
#		if(plotPoints){
#			par(mfrow=c(1,2))
#			plot(distsToTry,resultA,ylim=c(0,4),type="l",col="green",lwd=2)
#			lines(distsToTry,resultB,col="blue")
#			lines(distsToTry,resultC,col="gray")
#			lines(distsToTry,resultD,col="orange")
#			lines(c(bestMaxDist,bestMaxDist),c(0.9,1.1),col="red",lwd=2)
			

#			#print out cluster identity for each sequence
#			#print(paste("Cluster Assignments Chosen via Distance Threshold of bestMaxDist=",bestMaxDist))
#			#print(cbind(1:size,bestMaxDist))
			colrs=rainbow(length(unique(finalOrthoIDs))+2)
			blah<<-cmdscale(as.matrix(distances))*mean(as.matrix(distances))
			plot(blah[,1],blah[,2],cex=0,main=paste("Locus=",loc,"   nSeqs=",nrow(distances)),xlim=c(-xlm,xlm),ylim=c(-ylm,ylm))
			text(blah[,1],blah[,2],1:length(blah[,1]),cex=0.5,col=colrs[finalOrthoIDs])
#		}
		
		
		
#comment out line-drawing for now
if(FALSE){
		nSeqs=nrow(distances)
		
		print("Drawing clusters...")
		counter=0
		for(i in 1:nrow(distances)){
			for(j in i:nrow(distances)){
				if(distances[i,j]<maxDist && (abs(blah[i,1]-blah[j,1])>minSep || abs(blah[i,2]-blah[j,2])>minSep)){
					counter=counter+1
					lines(c(blah[i,1],blah[j,1]),c(blah[i,2],blah[j,2]),col=gray(distances[i,j]))
				}
			}
		}
		print(counter)
		text(blah[,1],blah[,2],1:length(blah[,1]),cex=0.5)

		clusterIDs <<-rep(0,nSeqs)
		clrs<<-rainbow(8)

		clusterIndex<<-1
		while(TRUE){
			print(paste("Draw a polygon around cluster ",clusterIndex,", then push esc. To end now, draw line, then push esc."))

			polygonA<<-locator(type="l",col=clrs[clusterIndex])
			len<<-len<<-length(polygonA$x)
			if(len==2){
				write(clusterIDs,paste("../OrthoSets/",project,"_L",loc,".txt",sep="",collapse=""),ncol=1)
				break
			}
			clusterIDs[in.out(cbind(c(polygonA$x,polygonA$x[len]),c(polygonA$y,polygonA$y[len])),cbind(blah[,1],blah[,2]))]=clusterIndex

			clusterIndex=clusterIndex+1
		}
}

}
}

plotSummary<<-function(nLoci){
	avgDists<<-rep(0,nLoci)
	nSeqs<<-rep(0,nLoci)
	for(loc in 1:nLoci){
		nSeqs[loc]=sqrt(length(scan(paste("../Distances/",project,"_L",loc,".txt",sep="",collapse=""))))
		
		if(nSeqs[loc]>2){
			distances<<-read.table(paste("../Distances/",project,"_L",loc,".txt",sep="",collapse=""))	
			avgDists[loc]=mean(as.matrix(distances))
		}
	}
	plot(nSeqs,avgDists,type="n")
	text(nSeqs,avgDists,1:nLoci,cex=0.5)
}