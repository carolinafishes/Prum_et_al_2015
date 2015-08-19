import java.io.*;
import java.util.concurrent.TimeUnit;

//version 2 is only different in that it tells the user that it will be creating the removeFile if not found
//there is also some code for skipping entire taxa, but it was disabled and so shoudl behave teh same as version 1

//version 3 has an additional parameter that allows for usage of alternative alignment folders for different ploidy specifications

///////////////////////////////////////////////////
//USAGE: java TrimAndMaskRawAlignments3 T45 403 96 20 14 0.5 5 ../T45_LociAndTaxaToRemove.txt
///////////////////////////////////////////////////
public class TrimAndMaskRawAlignments3{

  public static void main(String[] args){
      try{

		String project = args[0];								//e.g. T45
		int nLoci = Integer.parseInt(args[1]);					//e.g. 403
		int maxNTaxa = Integer.parseInt(args[2]);				//e.g. 96
		int windowSize = Integer.parseInt(args[3]);				//size of sliding window, e.g. 20
		int minGoodSites = Integer.parseInt(args[4]);			//e.g. 14 
		double propToBeGood = Double.parseDouble(args[5]); 		//e.g. 0.5
		int nTaxaMissingAllowed = Integer.parseInt(args[6]);	//e.g. 5
		String removeFile = args[7];	//e.g. T45_LociAndTaxaToRemove.txt    ...this file contains information on which loci x taxon combinations to remove, based on geneious visual inspection or other inspection
		String folderTail = "";
		System.out.println("args length = "+args.length);								//can be empty...
		if (args.length>8){
			folderTail=args[8];	//e.g. _2alleles
		}

		if(!removeFile.startsWith("../")){
			System.out.println("\n\n! ! ! ! WARNING: YOU PROBABLY WANT TO FIX THE PATH TO THE the LociAndTaxaToRemove file (e.g. add ../)\n"); 
			try{
				TimeUnit.MINUTES.sleep(10);
			}catch(InterruptedException ex){}		
		}

		if(!new File(removeFile).exists()){
			System.out.println("COULD NOT FIND REMOVE FILE, CREATING BLANK ONE...");
			BufferedWriter bwXXXX = new BufferedWriter ( new OutputStreamWriter(new FileOutputStream(   new File(removeFile) ) ));
			bwXXXX.close();
		}

		new File("../TrimmedAlignments"+folderTail).mkdir();


		boolean removeSequence[][] = new boolean[maxNTaxa][nLoci];	
		boolean skipTaxon[] = new boolean[maxNTaxa];

		BufferedReader br = new BufferedReader ( new InputStreamReader(new FileInputStream(   new File(removeFile) ) ));
		String tempS=br.readLine();
		while(tempS!=null){
			tempS=tempS.toUpperCase().replace(" ","\t");
			//get locus:
			int loc=-1;	//-1 indicates all loci for the individual
			if(!tempS.startsWith("ALL")){
				loc=Integer.parseInt(tempS.split("\t")[0]);
			}

			int ind=-1;	//-1 indicates all taxa for locus
			if(!tempS.endsWith("ALL")){
				ind=Integer.parseInt(tempS.split("\t")[1]);
			}

			if(ind==-1){System.out.println("Note that locus "+loc+" will be removed entirely because "+removeFile+" says to.");}
			if(loc==-1){System.out.println("Note that taxon "+ind+" will be removed entirely because "+removeFile+" says to.");}
			
			if(loc==-1){
//diable for now				skipTaxon[ind]=true;
//				System.out.println("Note that taxon "+ind+"will be removed entirely because "+removeFile+" says to.");
			}

			if(loc>0 && ind>0){	//a single sequence
				removeSequence[ind-1][loc-1]=true;
			}else if(loc>0 && ind<=0){	//remove an entire locus
				for(int i=0; i<maxNTaxa; i++){
					removeSequence[i][loc-1]=true;
				}
			}else if(loc<=0 && ind>0){
				for(int i=0; i<nLoci; i++){
					removeSequence[ind-1][i]=true;
				}				
			}else{
				System.out.println("INVALID CHOICE: "+tempS);
			}
			
			tempS=br.readLine();
		}
		br.close();


		BufferedWriter bwTrimmedLens = new BufferedWriter ( new OutputStreamWriter(new FileOutputStream(   new File("../Results/"+project+"_trimmedLens"+folderTail+".txt") ) ));
		BufferedWriter bwMissing = new BufferedWriter ( new OutputStreamWriter(new FileOutputStream(   new File("../Results/"+project+"_propMissing"+folderTail+".txt") ) ));

		int nChosenCharsByTaxon[][] = new int[maxNTaxa][maxNTaxa];

		int maxBaseTally[][]=new int[maxNTaxa+1][nLoci];

		int goodWindowTally[][]=new int[maxNTaxa][windowSize+1];
				
		int nTaxa=-1;
		
		for(int loc=1; loc<=nLoci; loc++){

			BufferedWriter bwTrimmed = new BufferedWriter ( new OutputStreamWriter(new FileOutputStream(   new File("../TrimmedAlignments"+folderTail+"/"+project+"_L"+loc+".fasta") ) ));
			
			System.out.println("L"+loc);
			char alignment[][]=getCharAlignment("../Alignments"+folderTail+"/"+project+"_L"+loc+".fasta",skipTaxon);
			String sequenceNames[]=getSequenceNames("../Alignments"+folderTail+"/"+project+"_L"+loc+".fasta",skipTaxon);
			nTaxa=alignment.length;
			if(nTaxa==0){
				alignment=new char[1][1];
				alignment[0][0] = 'N';
				sequenceNames=new String[1];
				sequenceNames[0]=">this locus was removed";
			}
			int nSites=alignment[0].length;
			
			//if(nChosenCharsByTaxon==null){ nChosenCharsByTaxon=new int[nTaxa][nTaxa]; }
						
			//tally the number of times each base is observed for each site
			int nA[] = new int[nSites];
			int nT[] = new int[nSites];
			int nC[] = new int[nSites];
			int nG[] = new int[nSites];
			int n_[] = new int[nSites];
			int nN[] = new int[nSites];	//also for other ambiguities
			int nPresent[] = new int[nSites];
			
			for(int site=0; site<nSites; site++){
				int maxTally=0;
				for(int seq=0; seq<nTaxa; seq++){
					switch(alignment[seq][site]){
						case 'A': nA[site]++; nPresent[site]++; maxTally=Math.max(nA[site],maxTally); break;
						case 'T': nT[site]++; nPresent[site]++; maxTally=Math.max(nT[site],maxTally); break;
						case 'C': nC[site]++; nPresent[site]++; maxTally=Math.max(nC[site],maxTally); break;
						case 'G': nG[site]++; nPresent[site]++; maxTally=Math.max(nG[site],maxTally); break;
						case '-': n_[site]++; break;
						default: nN[site]++;	nPresent[site]++; 
					}
				}
				maxBaseTally[maxTally][loc-1]++;
			}
			
			//trim the ends of each sequence until a windowSize consecuative ungapped bases contain at least minGoodSites bases that agree with at least propToBeGood% other sequences
			for(int seq=0; seq<nTaxa; seq++){
				//trim left end
				int firstGoodBase=nSites;	//begin by assuming there will be no good bases
				for(int site=0; site<=nSites-windowSize; site++){
					int nGood=0;
					int lastBadPos=-1;
					for(int i=0; i<windowSize; i++){
						switch(alignment[seq][site+i]){
							case 'A': if(nA[site+i]/(double)(nTaxa)>propToBeGood){nGood++;}else{lastBadPos=i;} break;
							case 'T': if(nT[site+i]/(double)(nTaxa)>propToBeGood){nGood++;}else{lastBadPos=i;} break;
							case 'C': if(nC[site+i]/(double)(nTaxa)>propToBeGood){nGood++;}else{lastBadPos=i;} break;
							case 'G': if(nG[site+i]/(double)(nTaxa)>propToBeGood){nGood++;}else{lastBadPos=i;} break;
							default: lastBadPos=i; break;
						}
					}
					goodWindowTally[seq][nGood]++;
					if(nGood>=minGoodSites){
						firstGoodBase=site+lastBadPos+1;
						break;
					}
				}
				
				//right trim end
				int lastGoodBase=-1;
				for(int site=nSites-1; site>=windowSize; site--){
					int nGood=0;
					int lastBadPos=-1;
					for(int i=0; i<windowSize; i++){
						switch(alignment[seq][site-i]){
							case 'A': if(nA[site-i]/(double)(nTaxa)>propToBeGood){nGood++;}else{lastBadPos=i;} break;
							case 'T': if(nT[site-i]/(double)(nTaxa)>propToBeGood){nGood++;}else{lastBadPos=i;} break;
							case 'C': if(nC[site-i]/(double)(nTaxa)>propToBeGood){nGood++;}else{lastBadPos=i;} break;
							case 'G': if(nG[site-i]/(double)(nTaxa)>propToBeGood){nGood++;}else{lastBadPos=i;} break;
							default: break;
						}
					}
					goodWindowTally[seq][nGood]++;
					if(nGood>=minGoodSites){
						lastGoodBase=site-(lastBadPos+1);
						break;
					}
				}
				
				//mask internally
				boolean needToMask[] = new boolean[nSites];
				for(int site=firstGoodBase; site<=lastGoodBase-windowSize; site++){
					int nGood=0;
				 	boolean containedGap=false;
					for(int i=0; i<windowSize; i++){
						switch(alignment[seq][site+i]){
							case 'A': if(nA[site+i]/(double)(nTaxa)>propToBeGood){nGood++;} break;
							case 'T': if(nT[site+i]/(double)(nTaxa)>propToBeGood){nGood++;} break;
							case 'C': if(nC[site+i]/(double)(nTaxa)>propToBeGood){nGood++;} break;
							case 'G': if(nG[site+i]/(double)(nTaxa)>propToBeGood){nGood++;} break;
							case '-': containedGap=true; break;
							default: break;
						}
					}
					goodWindowTally[seq][nGood]++;
					if(nGood<minGoodSites && !containedGap){
						for(int i=0; i<windowSize; i++){
							needToMask[site+i]=true;
						}
					}
				}
				
				//bwTrimmed.write(loc+"\t"+seq+"\t"+firstGoodBase+"\t"+lastGoodBase+"\n");
			

				//turn to - all sequences that should be removed based on removeFile
			
			
				//adjust the alignment accordingly
				//boolean oneGoodBaseFound=false;
				for(int site=0; site<nSites; site++){
					if(removeSequence[seq][loc-1]){		//remove entire sequence for those specified in removeFile
						alignment[seq][site]='-';
					}else if(!needToMask[site] && site>=firstGoodBase && site<=lastGoodBase){
						//oneGoodBaseFound=true;
						//leave alone
					}else{
						if( (needToMask[site] && alignment[seq][site]!='-')){// || (site==nSites-1 && !oneGoodBaseFound)){
							alignment[seq][site]='N';	//need to mask
						}else{
							alignment[seq][site]='-';	//need to trim
						}
					}
				}
			
/*			
				//write the new trimmed alignment
				bwTrimmed.write(">"+sequenceNames[seq]+"\n");
				boolean oneGoodBaseFound=false;
				int trimmedLen=0;
				int nCompleteEnough=0;
				for(int site=0; site<nSites; site++){
					if(!needToMask[site] && site>=firstGoodBase && site<=lastGoodBase){
						oneGoodBaseFound=true;
						bwTrimmed.write(alignment[seq][site]);
						trimmedLen++;
					}else{
						if( (needToMask[site] && alignment[seq][site]!='-') || (site==nSites-1 && !oneGoodBaseFound)){
							bwTrimmed.write("N");
						}else{
							bwTrimmed.write("-");
						}
					}
				}
				bwTrimmed.write("\n");
				bwTrimmedLens.write(trimmedLen+"\t");
				bwCompleteEnough.write(trimmedLen+"\t");
*/				
	
			}
			
			//determine which sites to keep based on n good bases in each site
			boolean goodSite[][] = new boolean[nTaxa][nSites];
			int nMissing[] = new int[nSites];
			for(int site=0; site<nSites; site++){
				int nTaxaMissing=0;
				for(int seq=0; seq<nTaxa; seq++){
					if(alignment[seq][site]=='-' ||  alignment[seq][site]=='N'){
						nTaxaMissing++;
					}
				}
//IS THERE AN ERROR HERE?
				nMissing[site]=nTaxaMissing;
				for(int x=0; x<nTaxa; x++){
					if(nTaxaMissing<=x){
						goodSite[x][site]=true;
					}
				}
			}
			
			//write the alignment
			for(int seq=0; seq<nTaxa; seq++){
				bwTrimmed.write(">"+sequenceNames[seq]+"\n");
				//count up the number of unambig bases, if less than 100, don't write
				int nUnambigCalls=0;
				for(int site=0; site<nSites; site++){
					if(goodSite[nTaxaMissingAllowed][site] && (alignment[seq][site]=='A' || alignment[seq][site]=='T' || alignment[seq][site]=='C' || alignment[seq][site]=='G')){
						nUnambigCalls++;
					}
				}
				boolean foundGoodSite=false;
				for(int site=0; site<nSites; site++){
					if(goodSite[nTaxaMissingAllowed][site]){
						if(nUnambigCalls>=100){
							bwTrimmed.write(""+alignment[seq][site]);
						}else{
							if(!foundGoodSite){
								bwTrimmed.write("n");
							}else{
								bwTrimmed.write("-");
							}
						}
						foundGoodSite=true;
					}
				}
				bwTrimmed.write("\n");
			}

			bwTrimmed.flush();
			bwTrimmed.close();

			//count up the number of sites that are good under each threshold			
			for(int x=0; x<nTaxa; x++){
				int nSitesGood=0;
				int nMissingChars=0;
				for(int site=0; site<nSites; site++){
					if(goodSite[x][site]){
						nSitesGood++;
						nMissingChars+=nMissing[site];	//add to tally if site is included
						for(int i=0; i<nTaxa; i++){
							if(alignment[i][site]!='-' && alignment[i][site]!='N'){nChosenCharsByTaxon[x][i]++;}
						}
					}
				}
				bwTrimmedLens.write(nSitesGood+"\t");
				bwMissing.write(nMissingChars/(double)(nSitesGood*nTaxa)+"\t");
			}
			bwTrimmedLens.write("\n");
			bwMissing.write("\n");

		}
		bwTrimmedLens.flush();
		bwTrimmedLens.close();
		bwMissing.flush();
		bwMissing.close();

		BufferedWriter bwMaxTally = new BufferedWriter ( new OutputStreamWriter(new FileOutputStream(   new File("../Results/"+project+"_maxTally"+folderTail+".txt") ) ));
		for(int i=0; i<nLoci; i++){
			for(int j=0; j<=maxNTaxa; j++){
				bwMaxTally.write(maxBaseTally[j][i]+"\t");
			}
			bwMaxTally.write("\n");
		}
		bwMaxTally.flush();
		bwMaxTally.close();

		BufferedWriter bwWindowTally = new BufferedWriter ( new OutputStreamWriter(new FileOutputStream(   new File("../Results/"+project+"_windowTally"+folderTail+".txt") ) ));
		for(int i=0; i<=windowSize; i++){
			for(int j=0; j<maxNTaxa; j++){
				bwWindowTally.write(goodWindowTally[j][i]+"\t");
			}
			bwWindowTally.write("\n");
		}
		bwWindowTally.flush();
		bwWindowTally.close();
		
		BufferedWriter bwCharsByTaxon = new BufferedWriter ( new OutputStreamWriter(new FileOutputStream(   new File("../Results/"+project+"_charsByTaxon"+folderTail+".txt") ) ));
		for(int i=0; i<maxNTaxa; i++){
			for(int j=0; j<maxNTaxa; j++){
				bwCharsByTaxon.write(nChosenCharsByTaxon[j][i]+"\t");	//each row is a different taxon
			}
			bwCharsByTaxon.write("\n");
		}
		bwCharsByTaxon.flush();
		bwCharsByTaxon.close();
		
      }catch(IOException ioe){System.out.println("<<!!ERROR main()!!>> MESSAGE:"+ioe.getMessage());}
  }

  static char[][] getCharAlignment(String filename, boolean skipTaxon[]){
   try{
	if(!new File(filename).exists()){return null;}

	//open the file and get the number of sequences and the length of the first sequences (alignment length)
	BufferedReader br = new BufferedReader ( new InputStreamReader(new FileInputStream(   new File(filename) ) ));
	String tempS=br.readLine();
	int nSeqs=0;
	int nSites=0;
	int count=0;
	while(tempS!=null){
		if(tempS.startsWith(">")){
			if(!skipTaxon[count]){
				nSeqs++;
			}
			count++;
		}else if(nSeqs==1){
			nSites+=tempS.length();
		}
		tempS=br.readLine();
	}
	br.close();
	char alignment[][] = new char[nSeqs][nSites];
	
	br = new BufferedReader ( new InputStreamReader(new FileInputStream(   new File(filename) ) ));
	tempS=br.readLine();
	int currSeq=0;
	int currSite=0;
	count=0;
	while(tempS!=null){
		if(tempS.startsWith(">")){
			if(!skipTaxon[count]){
				currSeq++;
				currSite=0;
			}
			count++;
		}else{
			if(!skipTaxon[count-1]){
				tempS=tempS.toUpperCase();
				for(int i=0; i<tempS.length(); i++){
					alignment[currSeq-1][currSite+i]=tempS.charAt(i);
				}
				currSite+=tempS.length();
			}
		}
		tempS=br.readLine();
	}
	br.close();
	return alignment;
    }catch(IOException ioe){System.out.println("<<!!ERROR main()!!>> MESSAGE:"+ioe.getMessage());}
	return null;
  }

  static String[] getSequenceNames(String filename, boolean skipTaxon[]){
   try{
	if(!new File(filename).exists()){return null;}
		//open the file and get the number of sequences and the length of the first sequences (alignment length)
		BufferedReader br = new BufferedReader ( new InputStreamReader(new FileInputStream(   new File(filename) ) ));
		String tempS=br.readLine();
		int nSeqs=0;
		int count=0;
		while(tempS!=null){
			if(tempS.startsWith(">")){
				if(!skipTaxon[count]){
					nSeqs++;
				}
				count++;
			}
			tempS=br.readLine();
		}
		br.close();
		String names[] = new String[count];

		br = new BufferedReader ( new InputStreamReader(new FileInputStream(   new File(filename) ) ));
		tempS=br.readLine();
		int currSeq=0;
		count=0;
		while(tempS!=null){
			if(tempS.startsWith(">")){
				if(!skipTaxon[count]){
					names[currSeq]=tempS.substring(1);
					currSeq++;
				}
				count++;
			}
			tempS=br.readLine();
		}
		br.close();
		return names;
    }catch(IOException ioe){System.out.println("<<!!ERROR main()!!>> MESSAGE:"+ioe.getMessage());}
	return null;
  }

}
