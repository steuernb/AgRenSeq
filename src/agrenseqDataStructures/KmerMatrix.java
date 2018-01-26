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



/**
 * This is a data structure to hold a kmer presence/absence matrix. For most applications except the creation or modification it is not efficient to hold the entire matrix in memory as this class is doing.
 * 
 * 
 * 
 * @author steuernb
 *
 */
public class KmerMatrix {

	
	HashMap<BitSet, BitSet> kmerMatrix;
	HashMap<String, Integer> accessions;
	int kmerSize;
	
	
	public static void main(String[] args) {
		try {
			KmerMatrix matrix = new KmerMatrix();
			matrix.addKmerSet(new File("/Users/steuernb/Documents/projects/AgRenSeq_CreateMatrix/kmer/BW_01000.dump.txt"), 10, "BW_01000");
			
		} catch (Exception e) {
			e.printStackTrace();
		}

	}

	
	public KmerMatrix(){
		this.kmerMatrix = new HashMap<BitSet, BitSet>();
		this.accessions = new HashMap<String, Integer>();
		this.kmerSize = 0;
	}
	
	/**
	 * 
	 * Add a kmer dump to the matrix. The input may be gzip. The format is tab separated. First colum is the kmer second colum is the count.
	 * 
	 * @param kmerDump
	 * 			input file with kmer dump.
	 * @param minCount
	 * 			minimum kmer count to consider a kmer present. This is to get rid of all the noise from sequencing errors.
	 * @param accession
	 * 			name of the accession
	 * 
	 * 
	 * @throws IOException
	 */
	public void addKmerSet(File kmerDump, int minCount, String accession)throws IOException{
		
		
		int maxAccessionIndex = -1;
		for(Iterator<String> iterator = this.accessions.keySet().iterator(); iterator.hasNext();){
			String key = iterator.next();
			int index = this.accessions.get(key);
			if(index > maxAccessionIndex){
				maxAccessionIndex = index;
			}
		}
		int accessionIndex = maxAccessionIndex + 1;
		accessions.put(accession, new Integer(accessionIndex));
		
		
		System.out.println("reading " + kmerDump.getAbsolutePath());
		
		int problemCounta = 0;
		
		BufferedReader in;
		FileInputStream fis = new FileInputStream(kmerDump);
		byte[] bytes = new byte[2];
		fis.read(bytes);
		int head = ((int) bytes[0] & 0xff) | ((bytes[1] << 8) & 0xff00);
		boolean gzip = GZIPInputStream.GZIP_MAGIC == head;
		fis.close();
		if(gzip){
			in = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(kmerDump))));
		}else{
			in = new BufferedReader((new FileReader(kmerDump)));
		}
	
		try{
			int kmerSize = in.readLine().split("\t")[0].length();
		
		if(this.kmerSize ==0){
			this.kmerSize = kmerSize;
		}else{
			if( this.kmerSize != kmerSize){
				System.err.println("Warning: found kmer size of "+kmerSize + ". Kmer size of the data set is "+this.kmerSize);
			}
		}
		}catch(NullPointerException e){
			System.err.println("Warning: empty data set " + kmerDump.getAbsolutePath());
		}
		in.close();
		
		
		
		if(gzip){
			in = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(kmerDump))));
		}else{
			in = new BufferedReader((new FileReader(kmerDump)));
		}
		
		for (String inputline = in.readLine(); inputline != null; inputline = in.readLine()) {
			String[] split  = inputline.split("\t");
			
			if( split[0].length() != kmerSize){
				problemCounta++;
				continue;
			}
			BitSet kmer = Kmer.convert(split[0]);
			
			int kmerCount = 0;
			try {
				
				kmerCount = Integer.parseInt(split[1]);
			} catch (NumberFormatException e) {
			
				problemCounta++;
				continue;
			} catch (ArrayIndexOutOfBoundsException e){
				problemCounta++;
				continue;
			}
			
			
			
			
			if(kmerCount >= minCount  ){
				
				if(!this.kmerMatrix.containsKey(kmer)){
					this.kmerMatrix.put(kmer, new BitSet());
				}
				kmerMatrix.get(kmer).set(accessionIndex);
				
			}
			
		}

		in.close();
		
		System.out.println("\tfinished reading "+kmerDump.getName()+". Found " + problemCounta + " problems. Size of data: " + this.kmerMatrix.size());
	}
	
	
	
	/**
	 * 
	 * Write a presence/absense matrix to a text file.
	 * First line starts with a # followed by a comma separated list of accession names. Every line starting with a "#" is not regarded.
	 * An entry in the matrix is a kmer, then a tab, then a string of 0 and 1 for absence and presence of that kmer in accessions. Order is according to first line of the matrix. Trailing 0s are ommitted.
	 *
	 * 
	 * 
	 * @param outputFile
	 * 		The output file
	 * 
	 * @throws IOException
	 */
	public void writePresenceMatrix(File outputFile)throws IOException{
		
		BufferedWriter out = new BufferedWriter(new FileWriter(outputFile));

		
		String[] accessionString = new String[this.accessions.size()];
		for(Iterator<String> iterator = this.accessions.keySet().iterator(); iterator.hasNext();){
			String accession = iterator.next();
			int index = this.accessions.get(accession);
			accessionString[index] = accession;
		}
		out.write("#");
		String s = "";
		for( int i = 0; i< accessionString.length; i++){
			s = s + "," + accessionString[i];
		}
		out.write(s.substring(1));
		out.newLine();
		
		for(Iterator<BitSet> iterator = this.kmerMatrix.keySet().iterator(); iterator.hasNext();){
			BitSet kmer = iterator.next();
			out.write(Kmer.convert(kmer, this.kmerSize) + "\t");
			
			BitSet presence = kmerMatrix.get(kmer);
			for(int i = 0; i< presence.length(); i++){
				if(presence.get(i)){
					out.write("1");
				}else{
					out.write("0");
				}
			}
			out.newLine();
		}
		
		
		out.close();
	}
	
	
	/**
	 * 
	 * Read in a presence/absense matrix.
	 * 
	 * @param inputFile
	 * 			The input file. Gzip is supported. First line is assumed to start with a # followed by a comma separated list of accession names. Every line starting with a "#" is not regarded.
	 * 			An entry in the matrix is a kmer, then a tab, then a string of 0 and 1 for absence and presence of that kmer in accessions. Order is according to first line of the matrix. Trailing 0s are ommitted.
	 * @throws IOException
	 */
	public void readPresenceMatrix(File inputFile)throws IOException{
		BufferedReader in;
		FileInputStream fis = new FileInputStream(inputFile);
		byte[] bytes = new byte[2];
		fis.read(bytes);
		int head = ((int) bytes[0] & 0xff) | ((bytes[1] << 8) & 0xff00);
		boolean gzip = GZIPInputStream.GZIP_MAGIC == head;
		fis.close();
		if(gzip){
			in = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(inputFile))));
		}else{
			in = new BufferedReader((new FileReader(inputFile)));
		}
	
		
		String[] header = in.readLine().substring(1).split(",");  //get rid of the #.
		for( int i = 0; i< header.length; i++){
			this.accessions.put( header[i], i);
		}
		
		for( String inputline = in.readLine(); inputline != null; inputline = in.readLine()){
			String[] split = inputline.split("\t");
			BitSet kmer = Kmer.convert(split[0]);
			char[] a = split[1].toCharArray();
			BitSet presence = new BitSet();
			for( int i = 0; i< a.length; i++){
				if(a[i] =='1'){
					presence.set(i);
				}
			}
			this.kmerMatrix.put(kmer, presence);
		}
		in.close();
	}
	
	
	/**
	 * 
	 * This removes kmers that are present in less accessions than threshold or all accessions but a number equal to the threshold.
	 * 
	 * @param threshold
	 * 			This will be the threshold then, I guess.	
	 * @throws IOException
	 */
	public void reduceMatrix(int threshold)throws IOException{
		HashSet<BitSet> removable = new HashSet<BitSet>();
		
		for( Iterator<BitSet> iterator = this.kmerMatrix.keySet().iterator(); iterator.hasNext();){
			BitSet kmer = iterator.next();
			BitSet presence = this.kmerMatrix.get(kmer);
			
			int count = 0;
			for( int i = 0; i< presence.length(); i++){
				if( presence.get(i)){
					count++;
				}
			}
			
			if(count < threshold  || this.accessions.size()-threshold < count){
				removable.add(kmer);
			}
			
		}
		
		for( Iterator<BitSet> iterator = removable.iterator(); iterator.hasNext();){
			this.kmerMatrix.remove(iterator.next());
		}
		
	}
	
	
	
}
