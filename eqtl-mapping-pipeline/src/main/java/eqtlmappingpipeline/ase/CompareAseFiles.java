/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.ase;

import java.io.BufferedReader;
import java.io.FileReader;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.regex.Pattern;
import umcg.genetica.io.trityper.EQTL;
import umcg.genetica.io.trityper.eQTLTextFile;

/**
 *
 * @author Patrick Deelen
 */
public class CompareAseFiles {

	private static final Pattern TAB_PATTERN = Pattern.compile("\\t");
	private static final Pattern COMMA_PATTERN = Pattern.compile(",");
	private static final int ASE_ESTIMATE_COLUMN = 2;
	private static final int ASE_CHR_COLUMN = 5;
	private static final int ASE_POS_COLUMN = 6;
	private static final int ASE_GENES_COLUMN = 12;
	private static final int ASE_A1_COLUMN = 9;

	/**
	 * @param args the command line arguments
	 */
	@SuppressWarnings("ManualArrayToCollectionCopy")
	public static void main(String[] args) throws Exception {
        
		BufferedReader aseReader = new BufferedReader(new FileReader("D:\\ASEcheck\\geuvadis_maskAll4_r20_a10_p2_s1_rq17_m1_gatkGenoGq30\\ase_bh.txt"));
        BufferedReader aseReader2 = new BufferedReader(new FileReader("D:\\ASEcheck\\all_maskAll4_r20_a10_p2_s5_rq17_m1_gatkGenoGq30\\ase_bh.txt"));
        
		HashMap<String, ArrayList<AseVariant>> AseEffects = new HashMap<String, ArrayList<AseVariant>>();

        String line;
		while ((line = aseReader.readLine()) != null) {
			String eQtlKey = eQtl.getRsChr() + ":" + eQtl.getRsChrPos();
			ArrayList<EQTL> posEqtls = eQtls.get(eQtlKey);
			if (posEqtls == null) {
				posEqtls = new ArrayList<EQTL>(1);
				eQtls.put(eQtlKey, posEqtls);
			}
			posEqtls.add(eQtl);
		}

		int aseTotal = 0;
		int aseWithEQtl = 0;
		int sameDirection = 0;
		int oppositeDirection = 0;

		HashSet<String> countedGenes = new HashSet<String>();

		aseReader.readLine();//header
		String line;
		String[] elements;
		while ((line = aseReader.readLine()) != null) {

			elements = TAB_PATTERN.split(line);


			HashSet<String> aseGenes = new HashSet<String>();
			for (String gene : COMMA_PATTERN.split(elements[ASE_GENES_COLUMN])) {
				aseGenes.add(gene);
			}

			++aseTotal;

			ArrayList<EQTL> posEqtls = eQtls.get(elements[ASE_CHR_COLUMN] + ":" + elements[ASE_POS_COLUMN]);
			if (posEqtls != null) {
				for (EQTL eQtl : posEqtls) {
					if (eQtl != null && aseGenes.contains(eQtl.getProbe())) {

						if (countedGenes.contains(eQtl.getProbe())) {
							continue;
						}
						countedGenes.add(eQtl.getProbe());

						//System.out.println(eQtl.getProbe());

						//if(eQtl.getRsChr() == 6 && eQtl.getRsChrPos() > 20000000 && eQtl.getRsChrPos() < 40000000) { continue; }

						++aseWithEQtl;


						double aseEstimate = Double.parseDouble(elements[ASE_ESTIMATE_COLUMN]);
						double eQtlZ = elements[ASE_A1_COLUMN].equals(eQtl.getAlleleAssessed()) ? eQtl.getZscore() : eQtl.getZscore() * -1;

						if (aseEstimate > 0.5 && eQtlZ > 0 || aseEstimate < 0.5 && eQtlZ < 0) {
							++sameDirection;
						} else {
							System.out.println("Opposite: " + eQtl.getRsChr() + ":" + eQtl.getRsChrPos() + "\t" + elements[ASE_A1_COLUMN] + "\t" + eQtl.getAlleleAssessed() + "\t" + aseEstimate + "\t" + eQtl.getZscore());
							++oppositeDirection;
						}
					}
				}
			}
		}

		NumberFormat numberFormat = NumberFormat.getInstance();
		numberFormat.setMinimumFractionDigits(2);
		numberFormat.setMaximumFractionDigits(2);
		System.out.println("Ase total: " + aseTotal);
		System.out.println("Ase SNP with eQTL effect: " + aseWithEQtl + " (" + numberFormat.format(aseWithEQtl / (double) aseTotal) + ")");
		System.out.println(" - Same direction: " + sameDirection + " (" + numberFormat.format(sameDirection / (double) aseWithEQtl) + ")");
		System.out.println(" - Opposite direction: " + oppositeDirection + " (" + numberFormat.format(oppositeDirection / (double) aseWithEQtl) + ")");

	}
}
