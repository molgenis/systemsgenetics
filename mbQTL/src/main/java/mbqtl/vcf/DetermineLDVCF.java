/*
 * DetermineLD.java
 *
 * Created on August 8, 2006, 4:32 PM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */
package mbqtl.vcf;


import umcg.genetica.containers.Pair;

/**
 * @author Lude Franke
 */
public class DetermineLDVCF {

	public final static int INCLUDE_CASES_AND_CONTROLS = 1;
	public final static int INCLUDE_CASES = 2;
	public final static int INCLUDE_CONTROLS = 3;

	public Pair<Double, Double> getLD(VCFVariant variant1, VCFVariant variant2) {
		// convert to byte[]
		if (variant1 == null || variant2 == null || variant1.getNrAlleles() > 2 || variant2.getNrAlleles() > 2) {
			return new Pair<>(Double.NaN, Double.NaN);
		}
		return getLD(variant1.getGenotypesAsByteVector(), variant2.getGenotypesAsByteVector(), null, INCLUDE_CASES_AND_CONTROLS, false);
	}

	public Pair<Double, Double> getLD(byte[] snp1Genotypes, byte[] snp2Genotypes, Boolean[] indIsCase, int individualsToInclude, boolean print) {

		double h11 = 0;
		double h12 = 0;
		double h21 = 0;
		double h22 = 0;
		int nrCalledGenotypes = 0;

		if (snp1Genotypes == null || snp2Genotypes == null) {
			return new Pair<Double, Double>(0d, 0d);
		}

		byte[] genotypesX = snp1Genotypes;
		byte[] genotypesY = snp2Genotypes;


		if (genotypesX.length != genotypesY.length) {
			return new Pair<Double, Double>(0d, 0d);
		}
		//Get genotypes:
		int[][] genotypes = new int[3][3];
		nrCalledGenotypes = 0;


		for (int ind = 0; ind < genotypesX.length; ind++) {
			if (individualsToInclude == INCLUDE_CASES_AND_CONTROLS || indIsCase == null
					|| (indIsCase[ind] != null && indIsCase[ind] && individualsToInclude == INCLUDE_CASES)
					|| (indIsCase[ind] != null && !indIsCase[ind] && individualsToInclude == INCLUDE_CONTROLS)) {
				int genotypeX = genotypesX[ind];
				int genotypeY = genotypesY[ind];
				if (genotypeX != -1 && genotypeY != -1) {
					genotypes[genotypeX][genotypeY]++;
					nrCalledGenotypes++;
				}
			}
		}
		//System.out.println("NrCalledGenotypes:\t" + nrCalledGenotypes);

		//Determine genotype frequencies:
		double[][] genotypesFreq = new double[3][3];
		for (int x = 0; x < 3; x++) {
			for (int y = 0; y < 3; y++) {
				genotypesFreq[x][y] = (double) genotypes[x][y] / (double) nrCalledGenotypes;
				if (print) {
					System.out.print(genotypes[x][y] + "\t");
				}
			}
			if (print) {
				System.out.println();
			}
		}

		//Determine alle frequencies:
		double[][] alleleFreq = new double[2][2];
		//SNP X:
		alleleFreq[0][0] = (genotypesFreq[0][0] + genotypesFreq[0][1] + genotypesFreq[0][2]) + (genotypesFreq[1][0] + genotypesFreq[1][1] + genotypesFreq[1][2]) / 2d;
		alleleFreq[0][1] = (genotypesFreq[2][0] + genotypesFreq[2][1] + genotypesFreq[2][2]) + (genotypesFreq[1][0] + genotypesFreq[1][1] + genotypesFreq[1][2]) / 2d;
		//SNP Y:
		alleleFreq[1][0] = (genotypesFreq[0][0] + genotypesFreq[1][0] + genotypesFreq[2][0]) + (genotypesFreq[0][1] + genotypesFreq[1][1] + genotypesFreq[2][1]) / 2d;
		alleleFreq[1][1] = (genotypesFreq[0][2] + genotypesFreq[1][2] + genotypesFreq[2][2]) + (genotypesFreq[0][1] + genotypesFreq[1][1] + genotypesFreq[2][1]) / 2d;

		if (print) {
			System.out.println("Allele freq:");
			System.out.println(alleleFreq[0][0] + "\tAA");
			System.out.println(alleleFreq[0][1] + "\tAB");
			System.out.println(alleleFreq[1][0] + "\tBA");
			System.out.println(alleleFreq[1][1] + "\tBB");
		}

		//Precalculate triangles of non-double heterozygote:
		double[][] genotypesTriangleFreq = new double[3][3];
		genotypesTriangleFreq[0][0] = 2d * genotypesFreq[0][0] + genotypesFreq[1][0] + genotypesFreq[0][1];
		genotypesTriangleFreq[2][0] = 2d * genotypesFreq[2][0] + genotypesFreq[1][0] + genotypesFreq[2][1];
		genotypesTriangleFreq[0][2] = 2d * genotypesFreq[0][2] + genotypesFreq[0][1] + genotypesFreq[1][2];
		genotypesTriangleFreq[2][2] = 2d * genotypesFreq[2][2] + genotypesFreq[1][2] + genotypesFreq[2][1];
		if (print) {
			System.out.println("Triangle freq:");
			System.out.println(genotypesTriangleFreq[0][0] + "\tAA");
			System.out.println(genotypesTriangleFreq[0][2] + "\tAB");
			System.out.println(genotypesTriangleFreq[2][0] + "\tBA");
			System.out.println(genotypesTriangleFreq[2][2] + "\tBB");
		}

		//Calculate expected genotypes, assuming equilibrium, take this as start:
		h11 = alleleFreq[0][0] * alleleFreq[1][0];
		h12 = alleleFreq[0][0] * alleleFreq[1][1];
		h21 = alleleFreq[0][1] * alleleFreq[1][0];
		h22 = alleleFreq[0][1] * alleleFreq[1][1];

		//Calculate the frequency of the two double heterozygotes:
		double x12y12 = h11 * h22 / (h11 * h22 + h12 * h21) * genotypesFreq[1][1];
		double x12y21 = h12 * h21 / (h11 * h22 + h12 * h21) * genotypesFreq[1][1];

		if (print) {
			System.out.println(h11 + "\t" + h12 + "\t" + h21 + "\t" + h22 + "\t\t" + x12y12 + "\t" + x12y21);
		}

		//Perform iterations using EM algorithm:
		for (int itr = 0; itr < 25; itr++) {

			h11 = (x12y12 + genotypesTriangleFreq[0][0]) / 2;
			h12 = (x12y21 + genotypesTriangleFreq[0][2]) / 2;
			h21 = (x12y21 + genotypesTriangleFreq[2][0]) / 2;
			h22 = (x12y12 + genotypesTriangleFreq[2][2]) / 2;

			x12y12 = h11 * h22 / (h11 * h22 + h12 * h21) * genotypesFreq[1][1];
			x12y21 = h12 * h21 / (h11 * h22 + h12 * h21) * genotypesFreq[1][1];

			if (print) {
				System.out.println(h11 + "\t" + h12 + "\t" + h21 + "\t" + h22 + "\t\t" + x12y12 + "\t" + x12y21);
			}

		}

		double d = h11 - (alleleFreq[0][0] * alleleFreq[1][0]);

//        if (returnType == RETURN_R_SQUARED) {
		double rSquared = d * d / (alleleFreq[0][0] * alleleFreq[0][1] * alleleFreq[1][0] * alleleFreq[1][1]);
		//return rSquared;
//        } else {
		double dMax = 0;
		if (d < 0) {
			double a = alleleFreq[0][1] * alleleFreq[1][1];
			if (alleleFreq[0][0] > alleleFreq[1][0]) {
				a = alleleFreq[0][0] * alleleFreq[1][0];
			}
			double b = alleleFreq[0][0] * alleleFreq[1][0];
			if (alleleFreq[0][0] > alleleFreq[1][0]) {
				b = alleleFreq[0][1] * alleleFreq[1][1];
			}
			dMax = Math.min(a, b);
		} else {
			double a = alleleFreq[0][1] * alleleFreq[1][0];
			if (alleleFreq[0][0] > alleleFreq[1][0]) {
				a = alleleFreq[0][0] * alleleFreq[1][1];
			}
			double b = alleleFreq[0][0] * alleleFreq[1][1];
			if (alleleFreq[0][0] > alleleFreq[1][0]) {
				b = alleleFreq[0][1] * alleleFreq[1][0];
			}
			dMax = Math.min(a, b);
		}
		double dPrime = Math.abs(d / dMax);
		/*
		 if (dPrime>1.01) {
         System.out.println("");
         System.out.println(genotypes[0][0] + "\t" + genotypes[0][1] + "\t" + genotypes[0][2]);
         System.out.println(genotypes[1][0] + "\t" + genotypes[1][1] + "\t" + genotypes[1][2]);
         System.out.println(genotypes[2][0] + "\t" + genotypes[2][1] + "\t" + genotypes[2][2]);
         System.out.println(alleleFreq[0][0] + "\t" + alleleFreq[0][1] + "\t" + alleleFreq[1][0] + "\t" + alleleFreq[1][1]);
         System.out.println(h11 + "\t" + h12 + "\t" + h21 + "\t" + h22);
         System.out.println(d + "\t" + dMax + "\t" + dPrime);
         }
         */
		return new Pair<Double, Double>(Math.min(1, dPrime), Math.min(1, rSquared));
//        }
	}

}
