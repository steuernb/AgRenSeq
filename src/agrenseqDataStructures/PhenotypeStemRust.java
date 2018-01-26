package agrenseqDataStructures;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;


public class PhenotypeStemRust extends Phenotype {

	
	
	public PhenotypeStemRust(File phenotypeFile)throws IOException{
		
		phenoScores = new HashMap<String,Double>();
		readScores(phenotypeFile);
		
	}
	
	
	
	
	public void readScores(File phenotypeFile)throws IOException{
		
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
	
	
	
	
}
