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
import java.util.regex.Pattern;
import umcg.genetica.io.trityper.EQTL;

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
        
        HashMap<String, AseVariantBean> aseEffects = new HashMap<String, AseVariantBean>();
        HashSet<String> totalEffects = new HashSet<String>();

        String line;
        String[] elements;
        //Header
        aseReader.readLine();
		while ((line = aseReader.readLine()) != null) {
            elements = TAB_PATTERN.split(line);
			aseEffects.put(elements[7],new AseVariantBean(elements));
		}
        aseReader.close();
        
        
		int aseTotalInFile2 = 0;
		int asePresentInBoth = 0;
		int sameDirection = 0;
		int oppositeDirection = 0;

		aseReader2.readLine();//header
		line = null;
		elements = null;
		while ((line = aseReader2.readLine()) != null) {
            ++aseTotalInFile2;
			elements = TAB_PATTERN.split(line);	
            AseVariantBean currentEffect = new AseVariantBean(elements);
            totalEffects.add(currentEffect.getId().toString());
            
            if(aseEffects.containsKey(currentEffect.getId().toString())){
                AseVariantBean otherEffect =  aseEffects.get(currentEffect.getId().toString());
                asePresentInBoth++;
                if(otherEffect.getEffect()<0.5 && currentEffect.getEffect()<0.5){
                    sameDirection++;
                } else if(otherEffect.getEffect()>0.5 && currentEffect.getEffect()>0.5){
                    sameDirection++;
                } else {
                    oppositeDirection++;
                }
            }
		}
        totalEffects.addAll(aseEffects.keySet());
        aseReader2.close();
        
		NumberFormat numberFormat = NumberFormat.getInstance();
		numberFormat.setMinimumFractionDigits(2);
		numberFormat.setMaximumFractionDigits(2);
		System.out.println("Ase effects in file 1: " + aseEffects.size());
        System.out.println("Ase effects in file 2: " + aseTotalInFile2);
        System.out.println("Total unique Ase effects observed : " + totalEffects.size());
		System.out.println("Ase effects identified in both analyses: " + asePresentInBoth + " (" + numberFormat.format(asePresentInBoth / (double) totalEffects.size()) + ")");
		System.out.println(" - Same direction: " + sameDirection + " (" + numberFormat.format(sameDirection / (double) asePresentInBoth) + ")");
		System.out.println(" - Opposite direction: " + oppositeDirection + " (" + numberFormat.format(oppositeDirection / (double) asePresentInBoth) + ")");

	}
}
