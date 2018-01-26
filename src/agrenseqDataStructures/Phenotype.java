package agrenseqDataStructures;


import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;


/**
 * 
 * This data structure holds the phenotypes of a diversity panel. 
 * The phenotype scores are assumed to be AgRenSeq scores, meaning the more positive the value is the more reistant is an accession and the more negative, the more susceptible is an accession. Intermediate should be around zero.
 * 
 * 
 * @author steuernb
 *
 */
public class Phenotype {

	
	HashMap<String,Double> phenoScores;
	
	/**
	 * Initialize an empty Phenotype data structure
	 */
	public Phenotype(){
		phenoScores = new HashMap<String, Double>();
	}
	
	
	/**
	 * 
	 * Initialize a phenotype data structure with given AgRenSeq scores
	 * 
	 * @param phenoScores
	 * 		AgRenSeq scores. A HashMap with accessions as keys and AgRenSeq scores as values.
	 * @throws IOException
	 */
	public Phenotype(HashMap<String,Double> phenoScores)throws IOException{
		
		this.phenoScores = phenoScores;
		
	}
	
	
	
	/**
	 * 
	 * Initialize a Phenotype data structure with a given phenotype file. 
	 * The File is assumed to be TSV and have accessions in the first column. The phenotype score in each line will be the average of all following columns.
	 * 
	 * @param phenotypeFile
	 * @throws IOException
	 */
	public Phenotype(File phenotypeFile)throws IOException{
		
		phenoScores = new HashMap<String,Double>();
		readScores(phenotypeFile);
		
	}
	
	/**
	 * 
	 * Get the HashMap that contains the AgRenSeq scores. 
	 * 
	 * @return
	 * 			A HashMap with accessions as keys and AgRenSeq scores as values.
	 */
	public HashMap<String, Double> getPhenotypes(){
		return this.phenoScores;
	}
	
	
	/**
	 * 
	 * read AgRenSeq scores from a phenotype file and add it to the data structure. First column has to be accessions. Following columns are AgRenSeq scores.
	 * The recorded score is the average from all scores in one line. Entries that are not numeric values are left out.
	 * 
	 * @param phenotypeFile
	 * @throws IOException
	 */
	private void readScores(File phenotypeFile)throws IOException{
		
		BufferedReader in = new BufferedReader(new FileReader(phenotypeFile));

		for (String inputline = in.readLine(); inputline != null; inputline = in.readLine()) {
			String[] split = inputline.split("\t");
			
			double d = 0;
			int count = 0;
			
			for( int i = 1; i< split.length; i++){
				try{
					d = d + Double.parseDouble(split[i]);
					count++;
					
				}catch(NumberFormatException e){}
			}
			
			phenoScores.put(split[0],  (d/count));
		}

		in.close();
		
	}
	
	
	
	/**
	 * provide a positive list of accessions. Every accession not in this list will be eliminated.
	 * @param usableFile
	 * 	A list of accession names.
	 * 
	 * @throws IOException
	 */
	public void selectAccessions(File usableFile)throws IOException{
		HashSet<String>  usable = new HashSet<String>();
		BufferedReader in = new BufferedReader(new FileReader(usableFile));

		for (String inputline = in.readLine(); inputline != null; inputline = in.readLine()) {
			String accession = inputline.trim().split("\t")[0];
			usable.add(accession);
		}

		in.close();
	
		HashSet<String> unusable =  new HashSet<String>();
		for(Iterator<String> iterator = this.phenoScores.keySet().iterator(); iterator.hasNext();){
			String accession = iterator.next();
			if( !usable.contains(accession)){
				unusable.add(accession);
			}
		}
		
		for(Iterator<String> iterator = unusable.iterator(); iterator.hasNext();){
			this.phenoScores.remove(iterator.next());
		}
		
		
	}
	
	/**
	 * provide a negative list of accessions. Every accession from this list will be eliminated
	 * @param unusableFile
	 * 	A list of accession names.
	 * @throws IOException
	 */
	public void removeAccessions(File unusableFile)throws IOException{
		BufferedReader in = new BufferedReader(new FileReader(unusableFile));

		for (String inputline = in.readLine(); inputline != null; inputline = in.readLine()) {
			String accession = inputline.trim().split("\t")[0];
			this.phenoScores.remove(accession);
		}

		in.close();
		
		
	}
	
	
	/**
	 * Get the maximum possible score for this phenotype. Essencially the sum of all positive entries.
	 * 
	 * @return
	 */
	public double getMaxScore(){
		double d = 0;
		for(Iterator<String> iterator = this.phenoScores.keySet().iterator(); iterator.hasNext();){
			String accession = iterator.next();
			double e = this.phenoScores.get(accession).doubleValue();
			if( e >0){
				d = d + e;
			}
		}
		return d;
	}
	
	
	/**
	 * 
	 * The total score is the score for a kmer that would be present in all accessions. 
	 * 
	 * @return
	 */
	public double getTotalScore(){
		double d = 0;
		for(Iterator<String> iterator = this.phenoScores.keySet().iterator(); iterator.hasNext();){
			String accession = iterator.next();
			double e = this.phenoScores.get(accession);
			d = d + e;
			
		}
		return d;
	}
	
	
	

	
	
}
