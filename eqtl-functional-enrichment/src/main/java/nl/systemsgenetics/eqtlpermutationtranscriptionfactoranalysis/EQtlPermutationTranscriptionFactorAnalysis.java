/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.eqtlpermutationtranscriptionfactoranalysis;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.regex.Pattern;
//import nl.umcg.deelenp.regulomedb.RegulomeDbEntry;
//import nl.umcg.deelenp.regulomedb.RegulomeDbFiles;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.trityper.TriTyperGenotypeData;
import org.molgenis.genotype.util.Ld;
import org.molgenis.genotype.util.LdCalculatorException;
import org.molgenis.genotype.variant.GeneticVariant;

/**
 *
 * @author Matthieu
 */
public class EQtlPermutationTranscriptionFactorAnalysis {
	
	private static final Pattern TAB_PATTERN = Pattern.compile("\\t");
	private RandomAccessGenotypeData eQtlGenotypeData;
	
	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) {
		// TODO code application logic here
		
		/*
		 * TODO:
		 * 1. Lees de eQTLMappingPipeline resultaten voor Imputed (eQTLProbesFDR0.05-ProbeLevel.txt).
		 * 2. Lees de RegulomeDB data in mbv de RegulomeDbReader.
		 * 3. Lees de Bloed Genotype Data in mbv de GenotypeReader.
		 */
	}
	
	
	/*
	 * Constructor that calls the other methods. The current code will change later.
	 */
	public EQtlPermutationTranscriptionFactorAnalysis()throws IOException{
		ArrayList<EQtl> eQtlResultData = this.readEQtlResultData("");
		//ArrayList<RegulomeDbEntry> regulomeDbData = readRegulomeDbData("");
		eQtlGenotypeData = readEQtlGenotypeData("");
	}
	
	
	
	public ArrayList<EQtl> readEQtlResultData(String eQtlResultsFile)throws IOException{
		
		ArrayList<EQtl> eqtldata = new ArrayList<EQtl>();
		String fileLine;
		String[] fileLineData;
		int n = 0;
		
		
		BufferedReader br = new BufferedReader( new FileReader(eQtlResultsFile) );
		while( (fileLine = br.readLine())!=null ){
			if(n != 0){
				fileLineData = TAB_PATTERN.split(fileLine);
				eqtldata.add(new EQtl(Double.valueOf(fileLineData[0]), fileLineData[1],
						fileLineData[2], Integer.parseInt(fileLineData[3]), fileLineData[4],
						Double.valueOf(fileLineData[10])));
			}
			n++;
		}
		br.close();
		
		return eqtldata;
		
		
		/*
		 * fileLineData[0] = pvalue
		 * fileLineData[1] = eQtlName
		 * fileLineData[2] = eQtlChr
		 * fileLineData[3] = eQtlChrPos
		 * fileLineData[4] = probeName
		 * fileLineData[10] = ZScore
		 */
	}
	
	
	/*
	public ArrayList<RegulomeDbEntry> readRegulomeDbData(String regulomeDbFilesPath)throws IOException{
		
		ArrayList<File> regulomeDbFiles = new ArrayList<File>();
		regulomeDbFiles.add( new File(regulomeDbFilesPath + "\\RegulomeDB.dbSNP132.Category1.txt") );
		regulomeDbFiles.add( new File(regulomeDbFilesPath + "\\RegulomeDB.dbSNP132.Category2.txt") );
		regulomeDbFiles.add( new File(regulomeDbFilesPath + "\\RegulomeDB.dbSNP132.Category3.txt") );
		regulomeDbFiles.add( new File(regulomeDbFilesPath + "\\RegulomeDB.dbSNP132.Category4.txt") );
		regulomeDbFiles.add( new File(regulomeDbFilesPath + "\\RegulomeDB.dbSNP132.Category5.txt") );
		regulomeDbFiles.add( new File(regulomeDbFilesPath + "\\RegulomeDB.dbSNP132.Category6.txt") );
		regulomeDbFiles.add( new File(regulomeDbFilesPath + "\\RegulomeDB.dbSNP132.Category7.txt") );
		
		
		ArrayList<RegulomeDbEntry> regulomedbdata = new ArrayList<RegulomeDbEntry>();
		RegulomeDbFiles rdbfs = new RegulomeDbFiles(regulomeDbFiles);
		
		
		Iterator<RegulomeDbEntry> rdbfIterator = rdbfs.iterator();
		while(rdbfIterator.hasNext()){
			RegulomeDbEntry rdbe = rdbfIterator.next();
			regulomedbdata.add(rdbe);
		}
		
		return regulomedbdata;
	}
	*/
	
	
	public RandomAccessGenotypeData readEQtlGenotypeData(String genotypeData) throws IOException{
		
		//RandomAccessGenotypeData gonlImputedBloodGenotypeData = new TriTyperGenotypeData( new File(genotypeData) );
		//return gonlImputedBloodGenotypeData;
		return null;
	}
	
	
	
	public void calculateLd(){
		//Use a window size of 250k: eQTL pos - 250k and eQTL pos + 250k
		
		/*
		 * Loop through the list of genetic variants (the eQTLs) and 
		 */
		
		
		
		
	}
	
	private Iterable<GeneticVariant> getNearByVariants(String chr, int windowSize){
		return this.eQtlGenotypeData.getVariantsByRange(chr, 0, windowSize);
	}
	
	
	
	public ArrayList<GeneticVariant> getEQtlsAsGenotypeSnps(ArrayList<EQtl> eqtls, RandomAccessGenotypeData genotypeData){
		
		ArrayList<GeneticVariant> eqtlSnps = new ArrayList<GeneticVariant>();
		for( EQtl eqtl : eqtls ){
			GeneticVariant gv = genotypeData.getSnpVariantByPos(eqtl.getEQtlChr(), eqtl.getEQtlChrPos()) ;
			eqtlSnps.add(gv);
		}
		return eqtlSnps;
	}
	
	
	
	
	public void performAnalysisStep(ArrayList<EQtl> eqtlData, int windowSize){
		Ld ld = null;
		
		for(EQtl eqtl : eqtlData){
			GeneticVariant eQtlSnp = this.eQtlGenotypeData.getSnpVariantByPos(eqtl.getEQtlChr(), eqtl.getEQtlChrPos());
			
			for(GeneticVariant gv : this.eQtlGenotypeData.getVariantsByRange(eqtl.getEQtlChr(), 0, windowSize)){
				try {
					ld = eQtlSnp.calculateLd(gv);
				} catch (LdCalculatorException ex) {
					System.out.println("Error in LD calculation: " + ex.getMessage());
					System.exit(1);
				}
				ld.getR2();
				ld.getDPrime();
				
				//Write results to a file.
				//Place results in a convenient structure for later.
			}
		}
	}
}
