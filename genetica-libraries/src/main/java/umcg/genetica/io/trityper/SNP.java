/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.io.trityper;

import umcg.genetica.math.stats.Descriptives;
import umcg.genetica.math.stats.HWE;

/**
 * @author harmjan
 */
public class SNP {
	
	private byte chr;
	private int chrPos;
	private int id;
	private String name;
	private byte[] dosage;
	private byte[] alleles1;
	private byte[] alleles2;
	private double HWEP = 1d;
	private double MAF = -1d;
	private double DOSAGEMAF = -1d;
	private double CR = 0d;
	private byte[] genotypes;
	private double[] alleleFreq;
	private byte minorAllele = 0;
	private boolean passesQC;
	private byte[] alleles;
	public int nrCalled;
	private int[] genotypeFreq;
	private byte alleleItr;
	private String[] alleleEncoding;
	
	
	public byte getAlleleItr() {
		return alleleItr;
	}
	
	public void setChr(byte chr) {
		this.chr = chr;
	}
	
	public void setChrPos(int chrPos) {
		this.chrPos = chrPos;
	}
	
	public byte getChr() {
		return this.chr;
	}
	
	public int getChrPos() {
		return this.chrPos;
	}
	
	public void setId(int snpId) {
		this.id = snpId;
	}
	
	public int getId() {
		return this.id;
	}
	
	public void setName(String snp) {
		this.name = snp;
	}
	
	public String getName() {
		return name;
	}
	
	public void setAlleles(byte[] allele1, byte[] allele2, Boolean[] indIncluded, Boolean[] indIsFemale) {
		this.alleles1 = allele1;
		this.alleles2 = allele2;
		genotypes = new byte[allele1.length];
		alleles = new byte[3];
		
		alleleItr = 0;
		int missingGenotypes = 0;
		genotypeFreq = new int[3];
		int nrTotal = 0;
		nrCalled = 0;
		boolean multiallelic = false;
		for (int ind = 0; ind < alleles1.length; ind++) {
			Boolean inc = indIncluded[ind];
			
			if ((inc != null && inc) && (chr != 23 || (indIsFemale[ind] == null || indIsFemale[ind]))) {
				nrTotal++;
			}
			
			if (allele1[ind] != 0 && allele2[ind] != 0 && (inc != null && inc)) {
				
				// X-chromosomal SNPs get different treatment, as they are heterozygous in all males
				if ((chr != 23 || (indIsFemale[ind] == null || indIsFemale[ind]))) {
					nrCalled++;
				}
				
				byte allelecode1 = -1;
				byte allelecode2 = -1;
				
				for (byte i = 0; i < 2; i++) {
					if (alleles[i] == allele1[ind]) {
						allelecode1 = i;
					}
				}
				
				if (allelecode1 == -1) {
					alleles[alleleItr] = allele1[ind];
					allelecode1 = alleleItr;
					alleleItr++;
					if (alleleItr > 2) {
						System.out.println("ERROR!!!: Number of different alleles for SNP\t" + name + "\t is more than two!");
						System.out.println("Allele 1:\t" + alleles[0] + " / " + new String(new byte[]{alleles[0]}));
						System.out.println("Allele 2:\t" + alleles[1] + " / " + new String(new byte[]{alleles[1]}));
						System.out.println("Allele 3:\t" + alleles[2] + " / " + new String(new byte[]{alleles[2]}));
						multiallelic = true;
						break;
					}
				}
				
				for (byte i = 0; i < 2; i++) {
					if (alleles[i] == allele2[ind]) {
						allelecode2 = i;
					}
				}
				if (allelecode2 == -1) {
					alleles[alleleItr] = allele2[ind];
					allelecode2 = alleleItr;
					alleleItr++;
					if (alleleItr > 2) {
						System.out.println("ERROR!!!: Number of different alleles for SNP\t" + name + "\t is more than two!");
						System.out.println("Allele 1:\t" + alleles[0] + " / " + new String(new byte[]{alleles[0]}));
						System.out.println("Allele 2:\t" + alleles[1] + " / " + new String(new byte[]{alleles[1]}));
						System.out.println("Allele 3:\t" + alleles[2] + " / " + new String(new byte[]{alleles[2]}));
						multiallelic = true;
						break;
					}
				}
				
				byte genotypeCode = (byte) (allelecode1 + allelecode2);
				genotypes[ind] = genotypeCode;
				
				// X-chromosomal SNPs get different treatment, as they are heterozygous in all males
				if ((chr != 23 || (indIsFemale[ind] == null || indIsFemale[ind]))) {
					genotypeFreq[genotypeCode]++;
				}
				
			} else {
				if (inc != null && inc) {
					missingGenotypes++;
				}
				
				genotypes[ind] = -1;
			}
		}
		
		
		if (!multiallelic) {
			
			// prune allele encoding, if the variant is multi-allelic, but only two alleles are present
			if (alleleEncoding != null && alleleEncoding.length > 2) {
				String[] tmpalleleEnc = new String[2];
				byte a1repl = 0;
				byte a2repl = 0;
				byte a1orig = alleles[0];
				byte a2orig = alleles[1];
				if (a1orig > a2orig) { // order of alleles in TriTyper is random; need to take this into account when recoding
					tmpalleleEnc[1] = alleleEncoding[alleles[0] - 100];
					tmpalleleEnc[0] = alleleEncoding[alleles[1] - 100];
					a1repl = 101;
					a2repl = 100;
				} else {
					tmpalleleEnc[0] = alleleEncoding[alleles[0] - 100];
					tmpalleleEnc[1] = alleleEncoding[alleles[1] - 100];
					a1repl = 100;
					a2repl = 101;
					
				}
				alleleEncoding = tmpalleleEnc;
				alleles[0] = a1repl;
				alleles[1] = a2repl;
				// replace alleles1 and alleles2
				for (int i = 0; i < alleles1.length; i++) {
					byte a1i = alleles1[i];
					byte a2i = alleles2[i];
					if (a1i == a1orig) {
						alleles[i] = a1repl;
					} else {
						alleles[i] = a2repl;
					}
					if (a2i == a2orig) {
						alleles2[i] = a2repl;
					} else {
						alleles2[i] = a1repl;
					}
				}
				
			}
			alleleFreq = new double[2];
			
			alleleFreq[0] = 2 * genotypeFreq[0] + genotypeFreq[1];
			alleleFreq[1] = 2 * genotypeFreq[2] + genotypeFreq[1];
			MAF = (alleleFreq[0]) / ((double) (nrCalled) * 2d);
			
			minorAllele = alleles[0];
			if (alleleFreq[0] > alleleFreq[1]) {
				minorAllele = alleles[1];
				MAF = 1 - MAF;
			}
			
			this.HWEP = HWE.calculateExactHWEPValue(genotypeFreq[1], genotypeFreq[0], genotypeFreq[2]);
			
			
			passesQC = false;
			
			CR = (double) nrCalled / nrTotal;
			
			if ((genotypeFreq[0] > 0 && genotypeFreq[1] > 0) || (genotypeFreq[1] > 0 && genotypeFreq[2] > 0) || (genotypeFreq[0] > 0 && genotypeFreq[2] > 0)) {
				passesQC = true;
			}
		}
		
	}
	
	public void setPassesQC(boolean b) {
		passesQC = b;
	}
	
	public boolean passesQC() {
		return passesQC;
	}
	
	/**
	 * @return the HWEP
	 */
	public double getHWEP() {
		return HWEP;
	}
	
	/**
	 * @param HWEP the HWEP to set
	 */
	public void setHWEP(Double HWEP) {
		this.HWEP = HWEP;
	}
	
	/**
	 * @return the MAF
	 */
	public double getMAF() {
		return MAF;
	}
	
	/**
	 * @return the MAF calculated on the basis of dosage values
	 */
	public double getDosageMAF() {
		return DOSAGEMAF;
	}
	
	/**
	 * @param MAF the MAF to set
	 */
	public void setMAF(Double MAF) {
		this.MAF = MAF;
	}
	
	/**
	 * @return the CR
	 */
	public double getCR() {
		return CR;
	}
	
	/**
	 * @param CR the CR to set
	 */
	public void setCR(Double CR) {
		this.CR = CR;
	}
	
	public void setAlleles(byte[] allele1, byte[] allele2) {
		this.alleles1 = allele1;
		this.alleles2 = allele2;
	}
	
	public void clearGenotypes() {
		this.alleles1 = null;
		this.alleles2 = null;
		this.genotypes = null;
		this.alleleFreq = null;
		this.dosage = null;
	}
	
	public byte[] getAllele1() {
		return alleles1;
	}
	
	public byte[] getAllele2() {
		return alleles2;
	}
	
	public byte[] getGenotypes() {
		return genotypes;
	}
	
	public void setDosage(byte[] dosageValues) {
		this.dosage = dosageValues;
		
		double[] dosages = getDosageValues();
		DOSAGEMAF = Descriptives.mean(dosages) / 2;
		if (DOSAGEMAF > 0.5) {
			DOSAGEMAF = 1 - DOSAGEMAF;
		}
	}
	
	public double[] getDosageValues() {
		if (dosage != null) {
			double[] dosagevalues = new double[dosage.length];
			for (int i = 0; i < dosage.length; i++) {
				dosagevalues[i] = ((double) (-Byte.MIN_VALUE + dosage[i])) / 100;
			}
			return dosagevalues;
		} else {
			return null;
		}
	}
	
	public boolean hasDosageInformation() {
		return (dosage != null);
	}
	
	public double[] selectGenotypes(int[] ids) {
		return selectGenotypes(ids, false, true);
	}
	
	public double[] selectGenotypes(short[] phenotypeToGenotypeId, boolean includeMissingGenotypes, boolean loadDosageWhenAvailable) {
		int[] idsInt = new int[phenotypeToGenotypeId.length];
		for (int i = 0; i < phenotypeToGenotypeId.length; i++) {
			idsInt[i] = phenotypeToGenotypeId[i];
		}
		return selectGenotypes(idsInt, includeMissingGenotypes, loadDosageWhenAvailable);
	}
	
	public double[] selectGenotypes(int[] phenotypeToGenotypeId, boolean includeMissingGenotypes, boolean loadDosageWhenAvailable) {
		int numReq = phenotypeToGenotypeId.length;
		int i;
		int numAvail = 0;
		
		for (i = 0; i < numReq; i++) {
			// if we're including missing genotypes, or the genotypes for the expression sample are unequal to -1
			if (includeMissingGenotypes || genotypes[phenotypeToGenotypeId[i]] != -1) {
				// if the expression sample has a genotype
				if (phenotypeToGenotypeId[i] != -1) {
					numAvail++;
				}
			}
		}
		
		double[] gtypes = new double[numAvail];
		int q = 0;
		for (i = 0; i < numReq; i++) {
			int l_id = phenotypeToGenotypeId[i];
			if (l_id != -1 && (includeMissingGenotypes || genotypes[l_id] != -1)) {
				if (dosage != null && loadDosageWhenAvailable) {
					double dosagevalue = ((double) (-Byte.MIN_VALUE + dosage[l_id])) / 100;
					gtypes[q] = dosagevalue;
				} else {
					gtypes[q] = genotypes[l_id];
				}
				q++;
			}
		}
		return gtypes;
	}
	
	public byte[] getAlleles() {
		return alleles;
	}
	
	public void setAlleleCodes(byte[] alleles) {
		this.alleles = alleles;
	}
	
	public byte getMinorAllele() {
		return minorAllele;
	}
	
	public double[] getAlleleFreq() {
		return alleleFreq;
	}
	
	/**
	 * @return the genotypeFreq
	 */
	public int[] getGenotypeFreq() {
		return genotypeFreq;
	}
	
	/**
	 * @param genotypeFreq the genotypeFreq to set
	 */
	public void setGenotypeFreq(int[] genotypeFreq) {
		this.genotypeFreq = genotypeFreq;
	}
	
	public void setMinorAllele(byte minorAllele) {
		this.minorAllele = minorAllele;
	}
	
	public void setNrCalled(int nrCalled) {
		this.nrCalled = nrCalled;
	}
	
	public int getNrCalled() {
		return nrCalled;
	}
	
	public void setAlleleEncoding(String[] alleles) {
		this.alleleEncoding = alleles;
	}
	
	public String[] getAlleleEncoding() {
		return alleleEncoding;
	}
	
	public boolean hasAlleleEncoding() {
		return (alleleEncoding != null);
	}
	
	public boolean isIndel() {
		if (alleleEncoding == null) {
			return false;
		}
		
		boolean indel = false;
		for (String d : alleleEncoding) {
			if (d.length() > 1) {
				indel = true;
			}
		}
		return indel;
	}
}
