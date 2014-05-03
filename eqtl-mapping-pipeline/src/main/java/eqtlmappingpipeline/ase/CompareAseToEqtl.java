/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.ase;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.HashMap;
import java.util.regex.Pattern;
import umcg.genetica.io.trityper.EQTL;
import umcg.genetica.io.trityper.eQTLTextFile;

/**
 *
 * @author Patrick Deelen
 */
public class CompareAseToEqtl {

	private static final Pattern TAB_PATTERN = Pattern.compile("\\t");
	
	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) throws Exception {
		
		eQTLTextFile eQTLsTextFile = new eQTLTextFile("D:\\UMCG\\Genetica\\Projects\\RnaSeqEqtl\\batch9_eQTLmapping\\result_geuvadis_maf0.05_call0.5_pcs100_normalizedPCA_meta\\eQTLSNPsFDR0.05-ProbeLevel.txt", false);
		
		HashMap<String, EQTL> eQtls = new HashMap<String, EQTL>();
		
		for(EQTL eQTL : eQTLsTextFile.read()){
			eQtls.put(eQTL.getRsChr() + ":" + eQTL.getRsChrPos(), eQTL);
		}
		
		int aseTotal = 0;
		int aseWithEQtl = 0;
		int sameDirection = 0;
		int oppositeDirection = 0;
				
		BufferedReader aseReader = new BufferedReader(new FileReader("D:\\UMCG\\Genetica\\Projects\\RnaSeqEqtl\\Ase\\test12\\ase_bonferroni_noNegativeCountR.txt"));
		aseReader.readLine();//header
		String line;
		String[] elements;
		while ((line = aseReader.readLine()) != null){
			
			elements = TAB_PATTERN.split(line);
						
			++aseTotal;
			
			EQTL eQtl = eQtls.get(elements[2] + ":" + elements[3]);
			if(eQtl != null){
				++aseWithEQtl;
				
				
				double aseZ = Double.parseDouble(elements[1]);
				
				if(elements[7].equals(eQtl.getAlleleAssessed())){
					if(aseZ < 0 && eQtl.getZscore() < 0 ){
						++sameDirection;
					} else {
						++oppositeDirection;
					}
				} else {
					if(aseZ < 0 && eQtl.getZscore() < 0 ){
						++oppositeDirection;
					} else {
						++sameDirection;
					}
				}
				
				
				
				
				
				
				
				
				
				
				
			}
			
			
			
			
		}
		
		System.out.println("Ase total: " + aseTotal);
		System.out.println("Ase SNP with eQTL effect: " + aseWithEQtl + " (" + aseWithEQtl / (double) aseTotal + ")");
		System.out.println(" - Same direction: " + sameDirection + " (" + sameDirection / (double) aseWithEQtl + ")");
		System.out.println(" - Opposite direction: " + oppositeDirection + " (" + oppositeDirection / (double) aseWithEQtl + ")");
	
	}
}
