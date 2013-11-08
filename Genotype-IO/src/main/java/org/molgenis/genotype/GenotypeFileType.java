/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype;

/**
 *
 * @author Patrick Deelen
 */
public enum GenotypeFileType {
	
	BIM,
	FAM,
	BED,
	PED,
	MAP,
	GEN,
	HAPS,
	SAMPLE,
	VCF;
	
	public static GenotypeFileType getTypeForPath(String path){
		
		if (path.toLowerCase().endsWith(".bim")) {
			return BIM;
		}
		
		if (path.toLowerCase().endsWith(".fam")) {
			return FAM;
		}
		
		if (path.toLowerCase().endsWith(".bed")) {
			return BED;
		}
		
		if (path.toLowerCase().endsWith(".ped")) {
			return PED;
		}
		
		if (path.toLowerCase().endsWith(".map")) {
			return MAP;
		}
		
		if (path.toLowerCase().endsWith(".gen")) {
			return GEN;
		}
		
		if (path.toLowerCase().endsWith(".haps")) {
			return HAPS;
		}
		
		if (path.toLowerCase().endsWith(".sample")) {
			return SAMPLE;
		}
		
		if (path.toLowerCase().endsWith(".vcf.gz")) {
			return VCF;
		}
		
		
		
		
		return null;
		
	}
	
}
