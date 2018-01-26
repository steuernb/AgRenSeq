package agrenseqDataStructures;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.BitSet;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.zip.GZIPInputStream;

import support.BioSequence;
import support.FastaReader;


/**
 * 
 * Calculate association scores for kmers and project them onto contigs of an assembly. 
 * An association score is the sum of AgRenSeq phenotype scores of accessions where the kmer is present.
 * The AgRenSeq phenotype score is suposed to be more positive the more the accession is resistant and more negative the more susceptible an accession is.
 * 
 * 
 * 
 * 
 * @author steuernb
 *
 */
public class KmerProjection {

	Phenotype phenotype;
	HashMap<BitSet, Double> associationMatrix;
	int kmerSize;
	
	
	File assemblyFile;
	File nlrList;
	File presenceMatrix;
	
	
	
	/**
	 * 
	 * 
	 * @param phenotype
	 * 
	 * @param kmerSize
	 * 
	 * @param assemblyFile
	 * 		The assemblyFile contains the denovo assembly of the accession where association scores are projected on. This is fasta format.
	 * 
	 * @param nlrList
	 * 		The nlrList is a tsv file with contig ids in the first column. Only contigs from the denovo assembly that are listed in first column of this file are regarded for projection.
	 * 
	 * @param presenceMatrix
	 * 		The presenceMatrix is a text file. Gzip is supported. First line is assumed to start with a # followed by a comma separated list of accession names. Every line starting with a "#" is not regarded.
	 * 		An entry in the matrix is a kmer, then a tab, then a string of 0 and 1 for absence and presence of that kmer in accessions. Order is according to first line of the matrix. Trailing 0s are ommitted.
	 *  
	 * 
	 */
	public KmerProjection( Phenotype phenotype, File assemblyFile, File nlrList, File presenceMatrix )throws IOException{
		
		this.associationMatrix = new HashMap<BitSet, Double>();
		this.phenotype = phenotype;
		this.assemblyFile = assemblyFile;
		this.nlrList = nlrList;
		this.presenceMatrix = presenceMatrix;
		this.setKmerSize();
		
	}
	
	
	
	
	
	/**
	 * 
	 * This takes the first entry in the presence absense matrix and uses that as the kmer length.
	 * 
	 * 
	 * 
	 * @throws IOException
	 */
	public void setKmerSize()throws IOException{
		
		BufferedReader in;
		FileInputStream fis = new FileInputStream(presenceMatrix);
		byte[] bytes = new byte[2];
		fis.read(bytes);
		int head = ((int) bytes[0] & 0xff) | ((bytes[1] << 8) & 0xff00);
		boolean gzip = GZIPInputStream.GZIP_MAGIC == head;
		fis.close();
		if(gzip){
			in = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(presenceMatrix))));
		}else{
			in = new BufferedReader(new FileReader(presenceMatrix));
		}
		String inputline = in.readLine();
		while(inputline!= null && inputline.startsWith("#")){
			inputline = in.readLine();
		}
		int kmerSize = inputline.split("\t")[0].length();
		
		in.close();
		
		this.kmerSize = kmerSize;
	}
	
	
	
	/**
	 * 
	 * This method will extract kmers from a denovo assembly. 
	 * The file nlrList is assumed to be a tab separated table. The first column is assumed to contain contig ids. Only contigs where IDs are in that column will be regarded by this method.
	
	 * @throws IOException
	 */
	public void readAssembly()throws IOException{
		
		System.out.println("Reading assembly " + assemblyFile.getAbsolutePath() + " using kmer size of " + this.kmerSize);
		
		HashSet<String> nlrContigs = new HashSet<String>();
		BufferedReader in = new BufferedReader(new FileReader(nlrList));
		for (String inputline = in.readLine(); inputline != null; inputline = in.readLine()) {
			nlrContigs.add(inputline.split("\t")[0]);
		}
		in.close();

		
		FastaReader fastaReader = new FastaReader(assemblyFile);
		for (BioSequence seq = fastaReader.readEntry(); seq != null; seq = fastaReader.readEntry()) {
			if(!nlrContigs.contains(seq.getIdentifier())){
				continue;
			}
			String sequence = seq.getSequence();
			for( int i = 0; i< sequence.length()-kmerSize; i++){
				BitSet kmer = Kmer.convert(sequence.substring(i, i+kmerSize));
				this.associationMatrix.put(kmer, 0.0);
			}
		}
		fastaReader.close();
		System.out.println("...finished. Recorded " + this.associationMatrix.size() + " kmers.");
	}
	
	
	/**
	 * This method will read through the file containing the presence/absense matrix and  
	 * 
	 * 
	 * @param presenceMatrix
	 * @throws IOException
	 */
	public void readMatrix()throws IOException{
		
		System.out.println("Reading presence/absense matrix " + presenceMatrix.getAbsolutePath());
		
		BufferedReader in;
		FileInputStream fis = new FileInputStream(presenceMatrix);
		byte[] bytes = new byte[2];
		fis.read(bytes);
		int head = ((int) bytes[0] & 0xff) | ((bytes[1] << 8) & 0xff00);
		boolean gzip = GZIPInputStream.GZIP_MAGIC == head;
		fis.close();
		if(gzip){
			in = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(presenceMatrix))));
		}else{
			in = new BufferedReader(new FileReader(presenceMatrix));
		}
		
		HashMap<String, Integer> accessions = new HashMap<String, Integer>();
		String[] header = in.readLine().substring(1).split(",");  //get rid of the #.
		for( int i = 0; i< header.length; i++){
			accessions.put( header[i], i);
		}
		
		for( String inputline = in.readLine(); inputline != null; inputline = in.readLine()){
			if( inputline.startsWith("#")){
				continue;
			}
			String[] split = inputline.split("\t");
			
			BitSet kmerF = Kmer.convert(split[0]);
			BitSet kmerR = Kmer.convert(new BioSequence("", split[0]).getReverseComplementarySequence()); //jellyfish does canonical kmers. Just checking if the reverse complement is there.
			
			double associationScore = 0;
			if( associationMatrix.containsKey(kmerF)){
				associationScore = getAssociationScore(split[1], accessions);
				this.associationMatrix.put(kmerF, associationScore);
			}
			
			if(associationMatrix.containsKey(kmerR)){
				if(associationScore ==0){
					associationScore = getAssociationScore(split[1], accessions);
				}
				this.associationMatrix.put(kmerR, associationScore);
			}
		}
		in.close();
		
		
	}
	
	
	/**
	 * 
	 * calculate an association score. 
	 * 
	 * @param presenceString
	 * 		String of 0 and 1 according to presence and absence of a kmer in accessions. 
	 * @param accessions
	 * 		This points accessions to their accoring index in the presnceString.
	 * @return
	 */
	private double getAssociationScore(String presenceString, HashMap<String, Integer> accessions){
		
		
		
		char[] a = presenceString.toCharArray();
		BitSet presence = new BitSet();
		for( int i = 0; i< a.length; i++){
			if(a[i] =='1'){
				presence.set(i);
			}
		}
		
		double associationScore = 0;
		for( Iterator<String> iterator = phenotype.getPhenotypes().keySet().iterator(); iterator.hasNext();){
			String accession = iterator.next(); 
			
			try{
				int index = accessions.get(accession);
				if(presence.get(index) && this.phenotype.getPhenotypes().containsKey(accession)){
					double agRenSeqScore = this.phenotype.getPhenotypes().get(accession);
					associationScore = associationScore + agRenSeqScore;
				}
			}catch(NullPointerException e){}	
		}
		return associationScore;
	}
	
	
	
	
	/**
	 * 
	 * Writes the projection of association scores onto NLR contigs. This will read the denovo assembly, only regards the contigs in the nlrList, and prints scores.
	 * 
	 * 
	 * @param outputFile
	 * 		This is the output of the kmer projection. A TSV file.
	 * 		First column is contig ID. 
	 * 		Second column is a sequence along contigIDs, integers incremented with every new contig. This is for convenience to print in R.
	 * 		Third column is a unique kmer score for a kmer that ocurs in that contig.
	 * 		Fourth column is the number of kmers that have that score.
	 * 		
	 * @throws IOException
	 */
	public void writeAssociationScore(File outputFile)throws IOException{
		
		System.out.println("Writing association output");
		
		HashSet<String> nlrContigs = new HashSet<String>();
		BufferedReader in = new BufferedReader(new FileReader(nlrList));
		for (String inputline = in.readLine(); inputline != null; inputline = in.readLine()) {
			nlrContigs.add(inputline.split("\t")[0]);
		}
		in.close();

		
		FastaReader fastaReader = new FastaReader(assemblyFile);
		BufferedWriter out = new BufferedWriter(new FileWriter(outputFile));
		int contigCount =0;
		
		
		for (BioSequence seq = fastaReader.readEntry(); seq != null; seq = fastaReader.readEntry()) {
			if(!nlrContigs.contains(seq.getIdentifier())){
				continue;
			}
			contigCount++;
			HashMap<Double, Integer> h = new HashMap<Double, Integer>();
			String sequence = seq.getSequence();
			for( int i = 0; i< sequence.length()-kmerSize; i++){
				BitSet kmer = Kmer.convert(sequence.substring(i, i+kmerSize));
				Double associationScore = this.associationMatrix.get(kmer);
				if( associationScore!= null && associationScore.doubleValue() != 0.0){
					int num = 0;
					if(h.containsKey(associationScore)){
						num = h.get(associationScore);
					}
					num++;
					h.put(associationScore, num);
				}
			}
			
			for( Iterator<Double> iterator = h.keySet().iterator(); iterator.hasNext();){
				double associationScore = iterator.next();
				int num = h.get(associationScore);
				out.write(seq.getIdentifier() + "\t" + contigCount + "\t" + associationScore + "\t" + num);
				out.newLine();
			}
			
		}
		fastaReader.close();
		out.close();
	}
	
	
	
	/**
	 * 
	 * Writes the projection of association scores onto contigs of an assembly.
	 * 
	 * @param contigList
	 * 		This HashSet<String> contains the list of contigs kmers should be printed for
	 * @param outputFile
	 * 		This is the output of the kmer projection. A TSV file.
	 * 		First column is the contig id.
	 * 		Second column is the position in the contig where a kmer starts.
	 * 		Third column is the association score of the kmer.
	 * 
	 * @throws IOException
	 */
	public void writeAssociationScorePerPosition( HashSet<String> contigList, File outputFile)throws IOException{
		FastaReader fastaReader = new FastaReader(assemblyFile);
		BufferedWriter out = new BufferedWriter(new FileWriter(outputFile));
		
		
		for (BioSequence seq = fastaReader.readEntry(); seq != null; seq = fastaReader.readEntry()) {
			if(!contigList.contains(seq.getIdentifier())){
				continue;
			}
			String sequence = seq.getSequence();
			for( int i = 0; i< sequence.length()-kmerSize; i++){
				BitSet kmer = Kmer.convert(sequence.substring(i, i+kmerSize));
				Double associationScore = this.associationMatrix.get(kmer);
				if( associationScore!= null && associationScore.doubleValue() != 0.0){
					out.write(seq.getIdentifier() + "\t" + i + "\t" + associationScore);
					out.newLine();
				}
			}
		}
		fastaReader.close();
		out.close();
	}
	
	
	
	
	
}
