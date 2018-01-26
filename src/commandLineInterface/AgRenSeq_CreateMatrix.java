package commandLineInterface;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;

import agrenseqDataStructures.KmerMatrix;
import support.CLI;


public class AgRenSeq_CreateMatrix {

	public static final double version = 1.0;
	
	
	public static void main(String[] args){
		
		
		CLI cli = new CLI();
		
		String help = 	"AgRenSeq_CreateMatrix version " + AgRenSeq_CreateMatrix.version + "\n"+
						"-i <metainfofile>\n" +
						"-o <outputFile>\n"	+	
						"-t <threshold> kmers present in less than <t> or all but <t> will not be recorded.\n"+
						"-n <minimum kmer count>";
						
		cli.parseOptions(args);
		
		
		
		
		try {
			
			if( !cli.hasOption("i") || (!cli.hasOption("o"))){
				throw new Exception("parameters -i and -o must be set.");
			}
			
			KmerMatrix matrix = new KmerMatrix();

			
			int minCount = 10;
			if( cli.hasOption("n")){
				try{
				minCount = Integer.parseInt(cli.getArg("n"));
				}catch (Exception e) {System.out.println("-n argument should be and int. Continuing with default value (10)");}
			}
			
			int threshold = 3;
			if(cli.hasOption ("t")){
				threshold = Integer.parseInt(cli.getArg("t"));
			}
			
			
			
			BufferedReader in = new BufferedReader(new FileReader(cli.getArg("i")));

			for (String inputline = in.readLine(); inputline != null; inputline = in.readLine()) {
				if(inputline.trim().startsWith("#")){
					continue;
				}
				String[] split = inputline.split("\t");
				String accession = split[0];
				File file = new File(split[1]);
				if(file.exists()){
					matrix.addKmerSet(file, minCount, accession);
				}else{
					System.out.println("File " + file.getAbsolutePath() + " does not exist.");
				}
				
			}

			in.close();
			
			
			matrix.reduceMatrix(threshold);
			
			matrix.writePresenceMatrix(new File(cli.getArg("o")));
			
		} catch (Exception e) {
			e.printStackTrace();
			System.out.println(help);
		}
		
		
		
		
	}
	
	
	
	
	
	
}
