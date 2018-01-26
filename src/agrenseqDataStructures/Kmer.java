package agrenseqDataStructures;


import java.util.BitSet;

import support.BioSequence;

/**
 * 
 * A collection of methods to convert a kmer string into a kmer BitSet class and back. 
 * Having 4 nucleotides, each nucleotide can be encoded by 2 bits, which greatly reduces memory to store it.
 * 
 * This is how I encode them.
 * A = 11
 * T = 00
 * G = 01
 * C = 10
 * 
 * 
 * @author steuernb
 *
 */
public class Kmer {

	
	
	/**
	 * 
	 * Get the reverse complement of a kmer in BitSet representation. This method is implemented in a rather naive way. 
	 * Convert it back to a String, get the reverse complement and convert it back to BitSet.
	 * 
	 * @param set
	 * @param kmerSize
	 * @return
	 */
	public static BitSet complement(BitSet set, int kmerSize){
		BioSequence seq = new BioSequence("", convert(set, kmerSize));
		return convert(seq.getReverseComplementarySequence());
	}
	
	/**
	 * 
	 * Convert a kmer string into a BitSet
	 * 
	 * @param kmer
	 * @return
	 */
	public static  BitSet convert(String kmer){
		
		BitSet set = new BitSet();
		
		char[] a = kmer.toUpperCase().toCharArray();
		
		for( int i = 0; i< a.length; i++){
			int index = i*2;
			if( a[i] =='A' || a[i] =='C' ){
				set.set(index);
			}
			if( a[i] =='A' || a[i] =='G'){
				set.set(index+1);
			}
			
		}
		
		return set;
	}
	
	
	/**
	 * 
	 * Convert a BitSet back to a String. Since trailing "0"s are not recorded we need to know the actual kmer length to know how many "T"s are at the end.
	 * 
	 * @param set
	 * 		Kmer encoded as BitSet
	 * @param kmerSize
	 * 		length of the kmer
	 * @return
	 */
	public static String convert(BitSet set, int kmerSize){
		StringBuffer buffer = new StringBuffer();
		
		for( int i = 0; i < kmerSize *2 ; i = i+2){
			if(set.get(i)){
				if( set.get(i+1)){
					buffer = buffer.append("A");
				}else{
					buffer = buffer.append("C");
				}
				
			}else{
				if( set.get(i+1)){
					buffer = buffer.append("G");
				}else{
					buffer = buffer.append("T");
				}
			}
		}
		return buffer.toString();
	}
	
	
}
