/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.eqtlpermutationtranscriptionfactoranalysis;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.regex.Pattern;
//import nl.umcg.deelenp.regulomedb.RegulomeDbEntry;
//import nl.umcg.deelenp.regulomedb.RegulomeDbFiles;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.trityper.TriTyperGenotypeData;
import org.molgenis.genotype.util.Ld;
import org.molgenis.genotype.util.LdCalculatorException;
import org.molgenis.genotype.variant.GeneticVariant;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.EQTL;
import umcg.genetica.io.trityper.eQTLTextFile;

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
	}
	
	
	/*
	 * Constructor that calls the other methods. The current code will change later.
	 */
	public EQtlPermutationTranscriptionFactorAnalysis()throws IOException{
		EQTL[] eQtlResultData = this.readEQtlResultData("");
		//ArrayList<RegulomeDbEntry> regulomeDbData = readRegulomeDbData("");
		eQtlGenotypeData = readEQtlGenotypeData("");
	}
	
	
	
	public EQTL[] readEQtlResultData(String eqtlFileLocation) throws IOException{
		eQTLTextFile eqtlData = new eQTLTextFile(eqtlFileLocation, false);
		return eqtlData.read();
	}
	
	
	public RandomAccessGenotypeData readEQtlGenotypeData(String genotypeData) throws IOException{
		RandomAccessGenotypeData gonlImputedBloodGenotypeData = new TriTyperGenotypeData( new File(genotypeData) );
		return gonlImputedBloodGenotypeData;
	}
	
	
	
	public void performAnalysisStep(EQTL[] eqtlData, int windowSize, String outputFilePath) throws IOException{
		HashMap<String, ArrayList<Ld>> calculateLd = calculateLd(eqtlData, windowSize, outputFilePath, 0.90);
	}
	
	
	public HashMap<String, ArrayList<Ld>> calculateLd(EQTL[] eqtlData, int windowSize, String outputFilePath, double r2CutOff)throws IOException{
		//Use a window size of 250k: eQTL pos - 250k and eQTL pos + 250k
		Ld ld = null;
		HashMap<String, ArrayList<Ld>> ldResultsMap = new HashMap<String, ArrayList<Ld>>();
		TextFile outputFile = new TextFile(outputFilePath, true);
		outputFile.write("SNP\tSNP2\tR2\tDPrime");
		
		
		for(EQTL eqtl : eqtlData){
			GeneticVariant eQtlSnp = this.eQtlGenotypeData.getSnpVariantByPos(eqtl.getRsChr().toString(), eqtl.getRsChrPos());
			
			for(GeneticVariant gv : this.eQtlGenotypeData.getVariantsByRange(eqtl.getRsChr().toString(), 0, windowSize)){
				try {
					ld = eQtlSnp.calculateLd(gv);
				} catch (LdCalculatorException ex) {
					System.out.println("Error in LD calculation: " + ex.getMessage());
					System.exit(1);
				}
				GeneticVariant variant1 = ld.getVariant1();
				GeneticVariant variant2 = ld.getVariant2();
				
				//Write results to a file.
				outputFile.write( variant1.getPrimaryVariantId() + "\t" + variant2.getPrimaryVariantId() + "\t" + ld.getR2());
				
				
				//Place results in a convenient structure for later.
				ArrayList<Ld> tmp;
				if(ldResultsMap.containsKey(variant1.getPrimaryVariantId())){
					tmp = ldResultsMap.get(variant1.getPrimaryVariantId());
					tmp.add(ld);
				}

				else{
					tmp = new ArrayList<Ld>();
					tmp.add(ld);
					ldResultsMap.put(variant1.getPrimaryVariantId(), tmp);
				}
			}
		}
		outputFile.close();
		
		return ldResultsMap;
	}
	
	
	public void findSnpsInRegulomeDb(ArrayList<String> regulomeDbFileLocations){
		
	}
	
	
	
	public ArrayList<GeneticVariant> getEQtlsAsGenotypeSnps(ArrayList<EQtl> eqtls, RandomAccessGenotypeData genotypeData){
		
		ArrayList<GeneticVariant> eqtlSnps = new ArrayList<GeneticVariant>();
		for( EQtl eqtl : eqtls ){
			GeneticVariant gv = genotypeData.getSnpVariantByPos(eqtl.getEQtlChr(), eqtl.getEQtlChrPos()) ;
			eqtlSnps.add(gv);
		}
		return eqtlSnps;
	}
	
	
	public RandomAccessGenotypeData getGenotypeData(){
		return this.eQtlGenotypeData;
	}
	
	private Iterable<GeneticVariant> getNearByVariants(String chr, int windowSize){
		return this.eQtlGenotypeData.getVariantsByRange(chr, 0, windowSize);
	} 
}
