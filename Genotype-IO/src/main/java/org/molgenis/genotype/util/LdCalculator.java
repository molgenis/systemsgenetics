package org.molgenis.genotype.util;

import java.util.LinkedHashMap;

import org.molgenis.genotype.variant.GeneticVariant;

public class LdCalculator
{

	/**
	 * LD calculator. Based on implementation of Harm-Jan Westra and Lude
	 * Franke.
	 * 
	 * 
	 * @param variant1
	 *            bi-allelic genetic variant
	 * @param variant2
	 *            bi-allelic genetic variant
	 * @return LD information
	 * @throws LdCalculatorException
	 */
	public static Ld calculateLd(GeneticVariant variant1, GeneticVariant variant2) throws LdCalculatorException
	{

		if (variant1 == null)
		{
			throw new IllegalArgumentException("Variant1 = null");
		}

		if (variant2 == null)
		{
			throw new IllegalArgumentException("Variant2 = null");
		}

		if (variant1.getAlleleCount() != 2 || variant2.getAlleleCount() != 2)
		{
			throw new UnsupportedOperationException("Ld calculator currently only supports biallelic variants");
		}

		final byte[] variant1Genotypes = variant1.getSampleCalledDosages();
		final byte[] variant2Genotypes = variant2.getSampleCalledDosages();

		if (variant1Genotypes.length != variant2Genotypes.length)
		{
			throw new LdCalculatorException("Error calculating LD: " + variant1.getPrimaryVariantId() + " contains "
					+ variant1Genotypes.length + " samples and " + variant2.getPrimaryVariantId() + " contains "
					+ variant2Genotypes.length + " samples. This should be identical");
		}

		// matrix with all combinations between variant 1 genotypes and variant
		// 2 genotypes
		int[][] genotypes = new int[3][3];

		int calledGenoypes = 0;

		for (int ind = 0; ind < variant1Genotypes.length; ++ind)
		{
			byte genotypeVariant1 = variant1Genotypes[ind];
			byte genotypeVariant2 = variant2Genotypes[ind];
			if (genotypeVariant1 != -1 && genotypeVariant2 != -1)
			{
				genotypes[genotypeVariant1][genotypeVariant2]++;
				++calledGenoypes;
			}
		}

		// matrix with freq for all combined genotypes
		double[][] genotypesFreq = new double[3][3];
		for (int x = 0; x < 3; x++)
		{
			for (int y = 0; y < 3; y++)
			{
				genotypesFreq[x][y] = (double) genotypes[x][y] / (double) calledGenoypes;
			}
		}

		// Determine allele freq of variants:
		double[][] alleleFreq = new double[2][2];
		// Variant 1:
		alleleFreq[0][0] = (genotypesFreq[0][0] + genotypesFreq[0][1] + genotypesFreq[0][2])
				+ (genotypesFreq[1][0] + genotypesFreq[1][1] + genotypesFreq[1][2]) / 2d;
		alleleFreq[0][1] = (genotypesFreq[2][0] + genotypesFreq[2][1] + genotypesFreq[2][2])
				+ (genotypesFreq[1][0] + genotypesFreq[1][1] + genotypesFreq[1][2]) / 2d;
		// Variant 2:
		alleleFreq[1][0] = (genotypesFreq[0][0] + genotypesFreq[1][0] + genotypesFreq[2][0])
				+ (genotypesFreq[0][1] + genotypesFreq[1][1] + genotypesFreq[2][1]) / 2d;
		alleleFreq[1][1] = (genotypesFreq[0][2] + genotypesFreq[1][2] + genotypesFreq[2][2])
				+ (genotypesFreq[0][1] + genotypesFreq[1][1] + genotypesFreq[2][1]) / 2d;

		// Precalculate triangles of non-double heterozygote:
		double[][] genotypesTriangleFreq = new double[3][3];
		genotypesTriangleFreq[0][0] = 2d * genotypesFreq[0][0] + genotypesFreq[1][0] + genotypesFreq[0][1];
		genotypesTriangleFreq[2][0] = 2d * genotypesFreq[2][0] + genotypesFreq[1][0] + genotypesFreq[2][1];
		genotypesTriangleFreq[0][2] = 2d * genotypesFreq[0][2] + genotypesFreq[1][2] + genotypesFreq[0][1];
		genotypesTriangleFreq[2][2] = 2d * genotypesFreq[2][2] + genotypesFreq[1][2] + genotypesFreq[2][1];

		// Calculate expected genotypes, assuming equilibrium, take this as
		// start:
		double h11 = alleleFreq[0][0] * alleleFreq[1][0];
		double h12 = alleleFreq[0][0] * alleleFreq[1][1];
		double h21 = alleleFreq[0][1] * alleleFreq[1][0];
		double h22 = alleleFreq[0][1] * alleleFreq[1][1];

		// Calculate the frequency of the two double heterozygotes:
		double x12y12 = (h11 * h22 / (h11 * h22 + h12 * h21)) * genotypesFreq[1][1];
		double x12y21 = (h12 * h21 / (h11 * h22 + h12 * h21)) * genotypesFreq[1][1];

		// Perform iterations using EM algorithm:
		for (int itr = 0; itr < 25; itr++)
		{

			h11 = (x12y12 + genotypesTriangleFreq[0][0]) / 2;
			h12 = (x12y21 + genotypesTriangleFreq[0][2]) / 2;
			h21 = (x12y21 + genotypesTriangleFreq[2][0]) / 2;
			h22 = (x12y12 + genotypesTriangleFreq[2][2]) / 2;

			x12y12 = (h11 * h22 / (h11 * h22 + h12 * h21)) * genotypesFreq[1][1];
			x12y21 = (h12 * h21 / (h11 * h22 + h12 * h21)) * genotypesFreq[1][1];

		}

		double d = h11 - (alleleFreq[0][0] * alleleFreq[1][0]);

		double rSquared = d * d / (alleleFreq[0][0] * alleleFreq[0][1] * alleleFreq[1][0] * alleleFreq[1][1]);

		double dMax = 0;
		if (d < 0)
		{
			double a = alleleFreq[0][1] * alleleFreq[1][1];
			if (alleleFreq[0][0] > alleleFreq[1][0])
			{
				a = alleleFreq[0][0] * alleleFreq[1][0];
			}
			double b = alleleFreq[0][0] * alleleFreq[1][0];
			if (alleleFreq[0][0] > alleleFreq[1][0])
			{
				b = alleleFreq[0][1] * alleleFreq[1][1];
			}
			dMax = Math.min(a, b);
		}
		else
		{
			double a = alleleFreq[0][1] * alleleFreq[1][0];
			if (alleleFreq[0][0] > alleleFreq[1][0])
			{
				a = alleleFreq[0][0] * alleleFreq[1][1];
			}
			double b = alleleFreq[0][0] * alleleFreq[1][1];
			if (alleleFreq[0][0] > alleleFreq[1][0])
			{
				b = alleleFreq[0][1] * alleleFreq[1][0];
			}
			dMax = Math.min(a, b);
		}
		double dPrime = Math.abs(d / dMax);

		// sometimes dPrime slightly larger then 1. Fixing this:
		dPrime = Math.min(1, dPrime);

		String variant1Alt = variant1.getVariantAlleles().get(1).getAlleleAsString();
		String variant1Ref = variant1.getVariantAlleles().get(0).getAlleleAsString();
		String variant2Alt = variant2.getVariantAlleles().get(1).getAlleleAsString();
		String variant2Ref = variant2.getVariantAlleles().get(0).getAlleleAsString();

		LinkedHashMap<String, Double> haplotypesFreq = new LinkedHashMap<String, Double>(4);
		haplotypesFreq.put(variant1Alt + "/" + variant2Alt, h11);
		haplotypesFreq.put(variant1Alt + "/"  + variant2Ref, h12);
		haplotypesFreq.put(variant1Ref + "/"  + variant2Alt, h21);
		haplotypesFreq.put(variant1Ref + "/"  + variant2Ref, h22);

		return new Ld(variant1, variant2, rSquared, dPrime, haplotypesFreq);

	}
    
    
    /**
	 * r2 calculator. Based on implementation of Harm-Jan Westra and Lude Franke.
     * 
     * Please be aware this comes with less tests and less EM rounds for speed.
	 * 
	 * @param variant1
	 *            bi-allelic genetic variant
	 * @param variant2
	 *            bi-allelic genetic variant
	 * @return LD information
	 */
    
	public static double calculateRsquare(GeneticVariant variant1, GeneticVariant variant2) throws LdCalculatorException{
        
        if (variant1.getAlleleCount() != 2 || variant2.getAlleleCount() != 2)	{
			throw new UnsupportedOperationException("Ld calculator currently only supports biallelic variants");
		}
        
		final byte[] variant1Genotypes = variant1.getSampleCalledDosages();
		final byte[] variant2Genotypes = variant2.getSampleCalledDosages();

		// matrix with all combinations between variant 1 genotypes and variant
		// 2 genotypes
		double[][] genotypesFreq = new double[3][3];

		double calledGenoypes = 0;

		for (int ind = 0; ind < variant1Genotypes.length; ++ind){
			byte genotypeVariant1 = variant1Genotypes[ind];
			byte genotypeVariant2 = variant2Genotypes[ind];
			if (genotypeVariant1 != -1 && genotypeVariant2 != -1){
				genotypesFreq[genotypeVariant1][genotypeVariant2]++;
				++calledGenoypes;
			}
		}

		for (int x = 0; x < 3; x++)	{
			for (int y = 0; y < 3; y++)	{
				genotypesFreq[x][y] /= calledGenoypes;
			}
		}

		// Determine allele freq of variants:
		double[][] alleleFreq = new double[2][2];
		// Variant 1:
		alleleFreq[0][0] = (genotypesFreq[0][0] + genotypesFreq[0][1] + genotypesFreq[0][2])
				+ (genotypesFreq[1][0] + genotypesFreq[1][1] + genotypesFreq[1][2]) / 2d;
		alleleFreq[0][1] = (genotypesFreq[2][0] + genotypesFreq[2][1] + genotypesFreq[2][2])
				+ (genotypesFreq[1][0] + genotypesFreq[1][1] + genotypesFreq[1][2]) / 2d;
		// Variant 2:
		alleleFreq[1][0] = (genotypesFreq[0][0] + genotypesFreq[1][0] + genotypesFreq[2][0])
				+ (genotypesFreq[0][1] + genotypesFreq[1][1] + genotypesFreq[2][1]) / 2d;
		alleleFreq[1][1] = (genotypesFreq[0][2] + genotypesFreq[1][2] + genotypesFreq[2][2])
				+ (genotypesFreq[0][1] + genotypesFreq[1][1] + genotypesFreq[2][1]) / 2d;

		// Precalculate triangles of non-double heterozygote:
		double[][] genotypesTriangleFreq = new double[2][2];
		genotypesTriangleFreq[0][0] = (2d * genotypesFreq[0][0] + genotypesFreq[1][0] + genotypesFreq[0][1])/2;
		genotypesTriangleFreq[1][0] = (2d * genotypesFreq[2][0] + genotypesFreq[1][0] + genotypesFreq[2][1])/2;
		genotypesTriangleFreq[0][1] = (2d * genotypesFreq[0][2] + genotypesFreq[1][2] + genotypesFreq[0][1])/2;
		genotypesTriangleFreq[1][1] = (2d * genotypesFreq[2][2] + genotypesFreq[1][2] + genotypesFreq[2][1])/2;

		// Calculate expected genotypes, assuming equilibrium, take this as start:
		double h11 = alleleFreq[0][0] * alleleFreq[1][0];
        double h11h22 = h11 * alleleFreq[0][1] * alleleFreq[1][1];
        double h12h21 = alleleFreq[0][0] * alleleFreq[1][1] * alleleFreq[0][1] * alleleFreq[1][0];
        double h11h12h21h22 = (genotypesFreq[1][1]/(h11h22+h12h21))/2;
        
		// Calculate the frequency of the two double heterozygotes:
		double x12y12 = (h11h22 * h11h12h21h22);
		double x12y21 = (h12h21 * h11h12h21h22);

		// Perform iterations using EM algorithm:
		for (int itr = 0; itr < 25; itr++){
			h11 = ((x12y12 + genotypesTriangleFreq[0][0]));
            h11h22 = h11 * ((x12y21 + genotypesTriangleFreq[1][1]));
            h12h21 = ((x12y21 + genotypesTriangleFreq[0][1])) * ((x12y12 + genotypesTriangleFreq[1][0]));
            h11h12h21h22 = (genotypesFreq[1][1]/(h11h22+h12h21))/2;
            
			x12y12 = (h11h22 * h11h12h21h22);
			x12y21 = (h12h21 * h11h12h21h22);
		}

		double d = h11 - (alleleFreq[0][0] * alleleFreq[1][0]);

        return (d * d / (alleleFreq[0][0] * alleleFreq[0][1] * alleleFreq[1][0] * alleleFreq[1][1]));
	}
}
