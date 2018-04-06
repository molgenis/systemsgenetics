package eqtlmappingpipeline.trans;

import eqtlmappingpipeline.metaqtl3.MetaQTL3;

public class TransQTL extends MetaQTL3 {
	
	public void run() {
		/* req:
		- thread for loading
		- thread for lz4+save
		- some calculation logic
		 - predefined chunks per chromosome
		 
		 - load a chunk of snps in memory
		 - run assoc in parallel for set of snps
		 - save results as shorts
		 
		 */
		
		// load expression data
		// load genotype data files
		// check whether SNPs are sorted
		// load chunk of data
		// divide snps over threads
		// calculate one snp per thread (i.e. 20k tests)
		// return:
		// - short for z-score table --> goes into matrix
		// - summary stats for SNP
		//    - sample size
		//    - MAF
		//    - callrate
		//    - hwe-p
		//    - alleles, allele-assessed
		// 
		// prepare output chunk (e.g. 100 chunks per chromsome, would be ~300mbytes, 5600 SNPs)
		// collect results
		// write to binary in chunks
		
		
	}
}
