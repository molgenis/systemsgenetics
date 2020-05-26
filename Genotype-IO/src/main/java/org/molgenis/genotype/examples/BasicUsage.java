/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.examples;

import java.io.IOException;
import java.util.HashMap;
import org.apache.log4j.Logger;
import org.molgenis.genotype.Allele;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.GenotypeData;
import org.molgenis.genotype.GenotypeDataException;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.RandomAccessGenotypeDataReaderFormats;
import org.molgenis.genotype.Sample;
import org.molgenis.genotype.plink.BedBimFamGenotypeData;
import org.molgenis.genotype.util.Ld;
import org.molgenis.genotype.util.LdCalculatorException;
import org.molgenis.genotype.variant.GeneticVariant;

/**
 * Untested usage examples
 * 
 * @author Patrick Deelen
 */
@edu.umd.cs.findbugs.annotations.SuppressWarnings(value="UCF_USELESS_CONTROL_FLOW_NEXT_LINE", justification="It is just an example")
public class BasicUsage {

	private static Logger LOGGER = Logger.getLogger(BasicUsage.class);
	
	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) {
		
		
		
		String datasetPath = "/some/path/file"; //Note omition of extentions
		
		//Create a binary plink genotype data object
		GenotypeData genotypeData = null;
		try {
			genotypeData = new BedBimFamGenotypeData(datasetPath, 1000);
		} catch (IOException ex) {
			LOGGER.fatal("IO error: " + ex.getMessage());
			System.exit(1);
		}
		
		for(GeneticVariant variant : genotypeData){
			//Iterate over all variants
			System.out.println(variant.getPrimaryVariantId());
		}
		
		for(Sample sample : genotypeData.getSamples()){
			//Note sex is an enum, to string outputs a human readable form of the gender
			System.out.println("Sample ID: " + sample.getId() + " sex: " + sample.getSex());
		}
		
		//Create a more advaced random access genotype 
		RandomAccessGenotypeData randomAccessGenotypeData = null;
		
		try {
			
			//Here we read the dataset of the specified type using the auto loader
			randomAccessGenotypeData = new BedBimFamGenotypeData(datasetPath);
			
			//Or
			randomAccessGenotypeData = RandomAccessGenotypeDataReaderFormats.PLINK_BED.createGenotypeData(datasetPath, 1000);
			
			//Or
			randomAccessGenotypeData = RandomAccessGenotypeDataReaderFormats.valueOf("PLINK_BED").createGenotypeData(datasetPath, 1000);
			
		} catch (IOException ex) {
			LOGGER.fatal("IO error: " + ex.getMessage());
			System.exit(1);
		} catch (GenotypeDataException ex) {
			LOGGER.fatal("Genotype data error: " + ex.getMessage());
			System.exit(1);
		}
		
		for(String sequenceNames : randomAccessGenotypeData.getSeqNames()){
			//No a sequence can be chr or contig
			System.out.println("Seq/Chr: " + sequenceNames);
		}
			
		for(GeneticVariant variant : randomAccessGenotypeData.getVariantsByRange("1", 1, 10)){
			//get variants from chr 1 pos 1 to 10 (exclusive)
			System.out.println("Variant ID: " + variant.getPrimaryVariantId());
		}
		
		//Create hashmap based on primary ID. 
		//Note this takes lots of memory. 
		//Variants without an ID will always be ignored.
		//If multiple variants have the same ID an arbritary variant is selected.
		HashMap<String, GeneticVariant> variantHashMap = randomAccessGenotypeData.getVariantIdMap();
		
		GeneticVariant snp = randomAccessGenotypeData.getSnpVariantByPos("1", 1);
		//Get snp variant at this position. null if not present
			
		snp.isSnp(); //must be true since we dit getSnpVariantByPos
		
		for(Alleles sampleAlleles : snp.getSampleVariants()){
			//Iterate over the alleles form all samples.
			System.out.println(sampleAlleles.getAllelesAsString());
		}
		
		for(byte dosage : snp.getSampleCalledDosages()){
			//Iterate over the dosage values of snps. (0,1,2)
			System.out.println(dosage);
		}
		
		for(float dosage : snp.getSampleDosages()){
			//Iterate over the dosage values of snps. range 0 - 2
			System.out.println(dosage);
		}
		
		snp.getMinorAllele();
		snp.getMinorAlleleFrequency();
		
		snp.getCallRate();
		
		snp.getHwePvalue();
		
		snp.isBiallelic();
		
		//chr != 0 && pos != 0
		snp.isMapped();
		
		snp.isAtOrGcSnp();
		
		GeneticVariant snp2 = null;
		//There can be multiple variants at a position. 
		for(GeneticVariant variantAtPos2 : randomAccessGenotypeData.getVariantsByPos("1", 2)){
			snp2 = variantAtPos2;
			break;
		}
		
		Ld ld = null;
		try {
			ld = snp.calculateLd(snp2);
		} catch (LdCalculatorException ex) {
			LOGGER.fatal("Error in LD calculation: " + ex.getMessage(), ex);
			System.exit(1);
		}
		ld.getR2();
		ld.getDPrime();
		
		Alleles snpAlleles =  snp.getVariantAlleles();
		
		for(Allele a : snpAlleles){
			//Print the alles found for this variant
			System.out.println("Allele: " + a);
		}
		
		Alleles alleles1 = Alleles.createAlleles(Allele.A, Allele.C);
		Alleles alleles2 = Alleles.createAlleles(Allele.C, Allele.A);
		Alleles alleles3 = Alleles.createAlleles(Allele.A, Allele.C);
		
		if(alleles1 == alleles2){}; // false because order is different
		if(alleles1 == alleles3){}; // true because it is guaranteed that if alleles and order is endentical it is the same object
		if(alleles1.sameAlleles(alleles2)){}; // true
		
		
		
	}
}
