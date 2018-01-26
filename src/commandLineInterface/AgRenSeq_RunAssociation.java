package commandLineInterface;
import java.io.File;

import agrenseqDataStructures.KmerProjection;
import agrenseqDataStructures.Phenotype;
import support.CLI;


public class AgRenSeq_RunAssociation {

	public static final double version = 1.0;
	
	
	public static void main(String[] args){
		
		
		CLI cli = new CLI();
		
		String help = 	"AgRenSeq_CreateMatrix version " + AgRenSeq_CreateMatrix.version + "\n"+
						"-i <matrix.txt>\n"      +
						"-n <nlr.txt>: List of contigs associated with nlrs" +
						"-o <outputFile.txt>\n"	 +	
						"-a <assembly.fasta> \n" +
						"-p <phenotype.txt> \n"  +
						"-u <usable>: optional list of accessions. If this is set all other accessions are ommitted.";
						
		cli.parseOptions(args);
		
		
		
		
		try {
			
			if( !cli.hasOption("i") || !cli.hasOption("o") || !cli.hasOption("p") || !cli.hasOption("a")){
				throw new Exception ("Parameters -i, -o, -p and -a are mandatory");
			}
			
			File inputMatrix = new File(cli.getArg("i"));
			if(!inputMatrix.exists()){
				throw new Exception("File " + cli.getArg("i") + " does not exist.");
			}
			
			File outputFile = new File(cli.getArg("o"));
			if(outputFile.exists()){
				throw new Exception("File " + cli.getArg("o") + " exists. I won't overwrite it. ");
			}
			
			File phenotypeFile = new File(cli.getArg("p"));
			if(!phenotypeFile.exists()){
				throw new Exception("File " + cli.getArg("p") + " does not exist.");
			}
			
			File assemblyFile = new File(cli.getArg("a"));
			if(!assemblyFile.exists()){
				throw new Exception("File " + cli.getArg("a") + " does not exist.");
			}
			
			File nlrList = new File(cli.getArg("n"));
			if(!nlrList.exists()){
				throw new Exception("File " + cli.getArg("n") + " does not exist.");
			}
			
			
			Phenotype p = new Phenotype(phenotypeFile);
			if( cli.hasOption("u")){
				p.selectAccessions(new File(cli.getArg("u")));
			}
			
			
			
			
			KmerProjection projection = new KmerProjection( p, assemblyFile, nlrList, inputMatrix );
			
			projection.readAssembly();
			projection.readMatrix();
			projection.writeAssociationScore( outputFile );
			
		} catch (Exception e) {
			e.printStackTrace();
			System.out.println(help);
		}
		
		
		
		
	}
	
	
	
	
	

}
