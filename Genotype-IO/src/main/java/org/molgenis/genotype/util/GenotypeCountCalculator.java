/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.util;

import java.util.ArrayList;
import java.util.List;
import org.apache.log4j.Logger;
import org.molgenis.genotype.Allele;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.variant.GeneticVariant;

/**
 *
 * @author Patrick Deelen
 */
public class GenotypeCountCalculator {
	
	private static final Logger LOGGER = Logger.getLogger(GenotypeCountCalculator.class);
	
	public static ArrayList<GenotypeCount> countGenotypes(GeneticVariant variant){
		
		ArrayList<GenotypeCount> genotypeCounts = new ArrayList<GenotypeCount>();
		
		List<Allele> alleles = variant.getVariantAlleles().getAlleles();
		
		for(int i = 0 ; i < alleles.size() ; ++i){
			for(int j = i ; j < alleles.size() ; ++j ){
				genotypeCounts.add(new GenotypeCount(Alleles.createAlleles(alleles.get(i), alleles.get(j)), 0));
			}
		}
		
		samples:
		for(Alleles sampleAlleles : variant.getSampleVariants()){
			
			if(sampleAlleles.contains(Allele.ZERO)){
				//skip missing
				continue samples;
			}
			
			counts:
			for(GenotypeCount genotypeCount : genotypeCounts){
				if(genotypeCount.getGenotype().sameAlleles(sampleAlleles)){
					genotypeCount.increment();
					continue samples;
				}
			}
			
			genotypeCounts.add(new GenotypeCount(sampleAlleles, 1));
			
		}
		
		
		return genotypeCounts;
		
	}
	
	public static class GenotypeCount {
		
		private final Alleles genotype;
		private int count;

		public GenotypeCount(Alleles genotype, int count) {
			this.genotype = genotype;
			this.count = count;
		}

		public int getCount() {
			return count;
		}

		public Alleles getGenotype() {
			return genotype;
		}
		
		public void increment(){
			++count;
		}
			
	}
	

	
}
