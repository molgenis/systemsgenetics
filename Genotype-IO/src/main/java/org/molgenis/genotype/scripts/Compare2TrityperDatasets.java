/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.scripts;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.Sample;
import org.molgenis.genotype.trityper.TriTyperGenotypeData;
import org.molgenis.genotype.variant.GeneticVariant;

/**
 *
 * @author Patrick Deelen
 */
public class Compare2TrityperDatasets {

	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) throws IOException {
		
		RandomAccessGenotypeData dataset1;
		RandomAccessGenotypeData dataset2;
		
		HashSet<String> sharedSamples = new HashSet<String>();
		HashMap<String, Integer> dataset1SampleMap = new HashMap<String, Integer>();
		HashMap<String, Integer> dataset2SampleMap = new HashMap<String, Integer>();
		
		try {
			dataset1 = new TriTyperGenotypeData(new File(args[0]));
		} catch (IOException ex) {
			System.err.println("Error reading dataset 1: " + args[0]);
			throw ex;
		}
		
		try {
			dataset2 = new TriTyperGenotypeData(new File(args[1]));
		} catch (IOException ex) {
			System.err.println("Error reading dataset 2: " + args[1]);
			throw ex;
		}
		
		int i = 0;
		for(Sample sample : dataset1.getSamples()){
			dataset1SampleMap.put(sample.getId(), i);
			++i;
		}
		
		i = 0;
		for(Sample sample : dataset2.getSamples()){
			dataset2SampleMap.put(sample.getId(), i);
			++i;
		}
		
		sharedSamples.addAll(dataset1SampleMap.keySet());
		sharedSamples.addAll(dataset2SampleMap.keySet());
		
		System.out.println("variantId\tcountIdentical\tcountMismatch");
		
		for(GeneticVariant dataset1Variant : dataset1){
			
			GeneticVariant dataset2Variant = dataset2.getSnpVariantByPos(dataset1Variant.getSequenceName(), dataset1Variant.getStartPos());
			
			if(dataset2Variant == null){
				//skipping non shared variants
				continue;
			}
			
			List<Alleles> dataset1VariantAlleles = dataset1Variant.getSampleVariants();
			List<Alleles> dataset2VariantAlleles = dataset2Variant.getSampleVariants();
			
			int countIdentical = 0;
			int countMismatch = 0;
			
			for(String sharedSample : sharedSamples){
				
				//Get alleles for this shared sample
				Alleles x = dataset1VariantAlleles.get(dataset1SampleMap.get(sharedSample));
				Alleles y = dataset2VariantAlleles.get(dataset2SampleMap.get(sharedSample));
				
				//Compare if same alleles. In this case AT == TA
				if(x.sameAlleles(y)){
					++countIdentical;
				} else {
					++countMismatch;
				}
				
			}			
			
			System.out.println(dataset1Variant.getPrimaryVariantId() + "\t" + countIdentical + "\t" + countMismatch);
			
		}
		
	}
}
