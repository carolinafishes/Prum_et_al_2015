import java.io.*;
import java.util.*;

////////////////////////////////////
//THIS PROGRAM COLLECTS THE SEQUENCES WITH THE GREATEST COVERAGE (NREADS MAPPED) TO USE IN ALIGNMENTS ACROSS SPECIES

public class GatherALLConSeqsWithOKCoverage{
  public static void main(String[] args){
      try{

		String project = args[0];							//e.g. P0040
		String taxonSetFile = args[1];						//e.g. ../P0040_TaxonSet.txt   //this file contains the taxon ids for the samples to be included in the alignment
		int nHomologs = Integer.parseInt(args[2]);			//e.g. 6
		int nLoci = Integer.parseInt(args[3]);				//e.g 263
		
		//count the number of taxa included
		BufferedReader br = new BufferedReader ( new InputStreamReader(new FileInputStream(   new File(taxonSetFile) ) ));	//in file							
		String tempS=br.readLine();
		int nInds=0;
		while(tempS!=null){
			nInds++;
			tempS=br.readLine();
		}
		
		br.close();
		//SETUP THE ARRAYS TO HOLD THE RESULTS
		String taxonSet[]=new String[nInds];
		String indSet[]=new String[nInds];
		
		//get the taxon identification
		br = new BufferedReader ( new InputStreamReader(new FileInputStream(   new File(taxonSetFile) ) ));	//in file
		for(int i=0; i<nInds; i++){
			taxonSet[i]=br.readLine();			//e.g. I4010_Cyperaceae_Carex_scoparia
			indSet[i]=taxonSet[i].split("_")[0];
		}
		br.close();

		//GET THE TAXON LIST FOUND IN THE MAPPEDREADS FILE
		br = new BufferedReader ( new InputStreamReader(new FileInputStream(   new File("../Results/"+project+"_AssemblySummary_nMappedReads.txt") ) ));	//in file
		tempS=br.readLine(); //skip supercontig label
		tempS=br.readLine(); //skip locus labels
		tempS=br.readLine(); //get first ind
		
		int nReadsMappedNTaxa=0;
		while(tempS.length()!=0){
			nReadsMappedNTaxa++;
			tempS=br.readLine();
		}
		br.close();
		
		int indexInTaxonSet[] = new int[nReadsMappedNTaxa];	//0 based
		
		br = new BufferedReader ( new InputStreamReader(new FileInputStream(   new File("../Results/"+project+"_AssemblySummary_nMappedReads.txt") ) ));	//in file
		tempS=br.readLine(); //skip supercontig label
		tempS=br.readLine(); //skip locus labels
		boolean foundTaxon[] = new boolean[nInds];	//keep track of which of the chosen taxa were found in the nReadsMappedFile, if not found, assume reference
		for(int i=0; i<nReadsMappedNTaxa; i++){
			indexInTaxonSet[i]=-1;
			tempS=br.readLine().split("\t")[0];
			for(int j=0; j<nInds; j++){	//find the taxon name in the list of chosen taxa
				if(indSet[j].equals(tempS)){
					indexInTaxonSet[i]=j;		//e.g.  for individual 21, i=20 corresponds to chosen taxon index 10
					foundTaxon[j]=true;
					break;
				}
			}
			if(indexInTaxonSet[i]==-1){System.out.println("Note: Taxon "+tempS+" found in ../Results/"+project+"_AssemblySummary_nMappedReads.txt was not found in "+taxonSetFile);}
		}
		br.close();
		
		//NOW MAKE A MATRIX TO KEEP TRACK OF THE CONSENSUS SEQUENCES THAT NEED TO BE KEPT
		boolean keepers[][][] = new boolean[nHomologs][nInds][nLoci]; //default false
		br = new BufferedReader ( new InputStreamReader(new FileInputStream(   new File("../Results/GoodConSeqIDs.txt") ) ));	//in file
		tempS=br.readLine();
		while(tempS!=null){
			if(indexInTaxonSet[Integer.parseInt(tempS.split("\t")[1])-1]==-1){tempS=br.readLine(); continue;}
			keepers[Integer.parseInt(tempS.split("\t")[0])-1][indexInTaxonSet[Integer.parseInt(tempS.split("\t")[1])-1]][Integer.parseInt(tempS.split("\t")[2])-1]=true;
			tempS=br.readLine();
		}
		br.close();

		//note the references as keepers
		for(int i=0; i<nInds; i++){
			//look for taxonID in 
			if(!foundTaxon[i]){
				System.out.println("Taxon "+taxonSet[i]+" is assumed to be a reference.");
				for(int j=0; j<nHomologs; j++){
					for(int k=0; k<nLoci; k++){
						keepers[j][i][k]=true;
					}
				}
			}
		}
		
		
		//CREATE THE DIRECTORY STRUCTURE FOR BOTH THE PREALIGNMENTS AND THE ALIGNMENTS
		new File("../Homologs").mkdir();
//		new File("../Alignments").mkdir();

		//now finally pull out the valid consensus sequences from each individual for each locus
//		BufferedWriter bwX=null;
		int numberOfCopies[][] = new int[nInds][nLoci];
		int scriptNumber=0;
		for(int loc=0; loc<nLoci; loc++){	//this will be slow...
			BufferedWriter bw	= new BufferedWriter ( new OutputStreamWriter(new FileOutputStream(   new File("../Homologs/"+project+"_L"+(loc+1)+".fasta") ) ));	//out file
			for(int i=0; i<nInds; i++){
				System.out.print("\rCollecting homologous sequences for locus "+(loc+1)+"..."+i+"\t"+i+"    ");


				if(!new File("../"+indSet[i]+"/"+indSet[i]+"_conSeqs.fasta").exists()){ 	//filter if no consensus sequences for this individual
//					bw.write(">"+taxonSet[i]+"\n");
//					bw.write("n\n"); 
					continue;	//note if file does not exist, the individual will not be represented
				}

				br = new BufferedReader ( new InputStreamReader(new FileInputStream(   new File("../"+indSet[i]+"/"+indSet[i]+"_conSeqs.fasta") ) ));	//in file
				tempS=br.readLine();	//e.g. >L1.2
				while(tempS!=null){
					int locFound=Integer.parseInt(tempS.substring(2,tempS.indexOf(".")));
					int copyFound=Integer.parseInt(tempS.substring(tempS.indexOf(".")+1));
					tempS=br.readLine();	//get sequence
					//System.out.println("!!!"+locFound+"!!!!");
					if(tempS.length()>20 && locFound-1==loc && copyFound<=nHomologs && keepers[copyFound-1][i][locFound-1]){	//relevant locus and keeper
						bw.write(">"+taxonSet[i]+"_Copy"+copyFound+"\n");
						bw.write(""+tempS+"\n");
						numberOfCopies[i][loc]++;
					}
					tempS=br.readLine();	//next header
				}
				br.close();
			}
			bw.flush();
			bw.close();

/*			
			if(loc%(nLoci/nAlignmentScripts)==0){
				scriptNumber++;
				if(bwX!=null){
					bwX.flush();
					bwX.close();
				}
				bwX = new BufferedWriter ( new OutputStreamWriter(new FileOutputStream(   new File("Align_AllViaMAFFT_"+scriptNumber+".sh") ) ));	//out file
				bwX.write("unset MAFFT_BINARIES\n");
			}
			bwX.write("mafft --genafpair --maxiterate 1000 --quiet ../Prealignments/"+project+"_L"+(loc+1)+".fasta > ../Alignments/"+project+"_L"+(loc+1)+".fasta\n");			
*/
		}
//		bwX.flush();
//		bwX.close();
		System.out.println();

		BufferedWriter bwCopy = new BufferedWriter ( new OutputStreamWriter(new FileOutputStream(   new File("../Results/"+project+"_CopyNumbers.txt") ) ));	//out file
		bwCopy.write("Sample\t");
		for(int j=0; j<nLoci; j++){bwCopy.write("L"+(j+1)+"\t");}
		bwCopy.write("\n");

		for(int i=0; i<nInds; i++){
			bwCopy.write(taxonSet[i]+"\t");
			for(int j=0; j<nLoci; j++){
				bwCopy.write(numberOfCopies[i][j]+"\t");
			}
			bwCopy.write("\n");
		}
		bwCopy.flush();
		bwCopy.close();
      }catch(IOException ioe){System.out.println("<<!!ERROR main()!!>>"+ioe.getMessage());}
  }

}
