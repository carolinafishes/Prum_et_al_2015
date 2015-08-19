import java.io.*;
import java.util.*;
import java.util.HashMap;

////////////////////////////////////
//THIS PROGRAM USES KMER MATCHES TO COMPUTE A MEASURE OF DISTANCE TO BE USED IN MULTIDIMENSIONAL SCALING FOR ORTHOLOGY ASSESSMENT

public class GetPairwiseDistanceMeasures{
  public static void main(String[] args){
      try{

		String project = args[0];							//e.g. P0040
		int nLoci = Integer.parseInt(args[1]);				//e.g 263
		int K = Integer.parseInt(args[2]);					//e.g. 20

		new File("../Distances").mkdir();

		for(int loc=1; loc<=nLoci; loc++){
			System.out.print("\rComputing distances for locus "+loc);
			//create a hashmap to store the observed kmers
			HashMap kmerMap = new HashMap();



			//open the file containing the sequences to be compared, count the number of sequences
			BufferedReader br = new BufferedReader ( new InputStreamReader(new FileInputStream(   new File("../Homologs/"+project+"_L"+loc+".fasta") ) ));	//in file
			String tempS=br.readLine();	//header e.g. >I2982_CER85.2_Gobiidae_Chlamydogobius_eremius_Copy1
			int nSeqs=0;
			while(tempS!=null){
				tempS=br.readLine();
				nSeqs++;
				tempS=br.readLine();
			}
			br.close();

			//reopen the file and process the sequences
			br = new BufferedReader ( new InputStreamReader(new FileInputStream(   new File("../Homologs/"+project+"_L"+loc+".fasta") ) ));	//in file
			tempS=br.readLine();	//header e.g. >I2982_CER85.2_Gobiidae_Chlamydogobius_eremius_Copy1
			int seqID=0;
			String kmer;
			boolean value[];
			int nKmers[]=new int[nSeqs];
			
			while(tempS!=null){
					seqID++;
					tempS=br.readLine().toUpperCase();	//e.g. gcagctagtagtagggcagtaggttacgtggatggatcatgGGAGATTGATCCTAGACTAGT

					//hash all kmers found...continuous kmers
					for(int i=0; i<=tempS.length()-K; i++){
						kmer=tempS.substring(i,i+K);
						if(kmer.indexOf("N")>=0){continue;}
						nKmers[seqID-1]++;						
						value = (boolean[])kmerMap.get(kmer);
						if(value==null){
							value = new boolean[nSeqs];
							kmerMap.put(kmer,value);
						}
						value[seqID-1]=true;
					}
					
					//hash all kmers found...spaced in three frames...
					for(int i=0; i<=tempS.length()-3*K; i++){
						kmer="";
						for(int j=0; j<K; j++){
							kmer+=tempS.charAt(i+j*3);
						}
						if(kmer.indexOf("N")>=0){continue;}
						nKmers[seqID-1]++;
						value = (boolean[])kmerMap.get(kmer);
						if(value==null){
							value = new boolean[nSeqs];
							kmerMap.put(kmer,value);
						}
						value[seqID-1]=true;
					}					

					tempS=br.readLine();
			}
			br.close();
			
			
			//now tally the matches to generate a pairwise similarity matrix
			int tallies[][] = new int[nSeqs][nSeqs];
			Iterator it2 = kmerMap.entrySet().iterator();
			while (it2.hasNext()) { 
				Map.Entry entry = (Map.Entry) it2.next();
				value = (boolean[])entry.getValue();
				for(int i=0; i<nSeqs; i++){
					if(!value[i]){continue;}
					for(int j=i; j<nSeqs; j++){
						if(value[j]){tallies[i][j]++; if(i!=j){tallies[j][i]++;}}
					}
				}
			}
		//System.out.println();
		//for(int i=0; i<nSeqs; i++){System.out.println(nKmers[i]);}
		
			//write the matrix to a file as a distance matrix
			BufferedWriter bw = new BufferedWriter ( new OutputStreamWriter(new FileOutputStream(   new File("../Distances/"+project+"_L"+loc+".txt") ) ));	//out file
			for(int i=0; i<nSeqs; i++){
				for(int j=0; j<nSeqs; j++){
					if(i==j){
						bw.write("0"+"\t");
					}else{
						bw.write(""+((Math.min(nKmers[i],nKmers[j])-tallies[i][j])/(double)Math.min(nKmers[i],nKmers[j]))+"\t");
					}
				}
				bw.write("\n");
			}
			bw.flush();
			bw.close();

		}
		System.out.println();

      }catch(IOException ioe){System.out.println("<<!!ERROR main()!!>>"+ioe.getMessage());}
  }

}

