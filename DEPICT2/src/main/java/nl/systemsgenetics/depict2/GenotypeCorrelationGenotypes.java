/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.depict2;

import java.util.ArrayList;
import jsat.math.FastMath;
import org.apache.commons.math3.stat.regression.SimpleRegression;
import org.apache.log4j.Logger;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.variant.GeneticVariant;

/**
 *
 * @author patri
 */
public class GenotypeCorrelationGenotypes implements GenotypeCorrelationSource {

	private static final Logger LOGGER = Logger.getLogger(GenotypeCorrelationGenotypes.class);
	
	private final RandomAccessGenotypeData referenceGenotypes;
	

	public GenotypeCorrelationGenotypes(RandomAccessGenotypeData referenceGenotypes) {
		this.referenceGenotypes = referenceGenotypes;
	}

	@Override
	public GenotypieCorrelationResult getCorrelationMatrixForRange(String chr, int start, int stop, double maxR) {

		final SimpleRegression regression = new SimpleRegression();

		final ArrayList<GeneticVariant> includedVariantsList = new ArrayList<>();

		start = start < 0 ? 0 : start;
		
		LOGGER.debug("Query genotype data: " + chr + ":" + start + "-" + stop);
		
		int variantsFoundInRegion = 0;
		int variantsExcludedDueToHighR = 0;
		
		newVariants:
		for (GeneticVariant newVariant : referenceGenotypes.getVariantsByRange(chr, start, stop)) {
			final float[] newVariantDosages = newVariant.getSampleDosages();
			variantsFoundInRegion++;
			for (GeneticVariant selectedVariant : includedVariantsList) {
				final float[] selectedVariantDosages = selectedVariant.getSampleDosages();

				//This loop is used to calculate correlation between variant dosages
				for (int i = 0; i < newVariantDosages.length; ++i) {
					regression.addData(newVariantDosages[i], selectedVariantDosages[i]);
				}

				//If correlation is too large stop with current newVariant and move to next variant
				if (Math.abs(regression.getR()) >= maxR) {
					variantsExcludedDueToHighR++;
					continue newVariants;
				}
				regression.clear();
			
			}
			includedVariantsList.add(newVariant);
		}
		
		LOGGER.debug(" * Variants found in region: " + variantsFoundInRegion);
		LOGGER.debug(" * Variants excluded due to high correlation: " + variantsExcludedDueToHighR);

		if (includedVariantsList.isEmpty()) {
			return new GenotypieCorrelationResult(new double[0][0], new String[0]);
		} else {
			
			final String[] includedVariants = new String[includedVariantsList.size()];

			final double[][] cor = new double[includedVariantsList.size()][includedVariantsList.size()];
			for (int p = 0; p < includedVariantsList.size(); p++) {
				
				GeneticVariant pVariant = includedVariantsList.get(p);
				
				includedVariants[p] = pVariant.getPrimaryVariantId();
				
				final float[] pVariantDosages = pVariant.getSampleDosages();
				cor[p][p] = 1;
				for (int q = p + 1; q < includedVariantsList.size(); q++) {
					final float[] qVariantDosages = includedVariantsList.get(q).getSampleDosages();

					//This loop is used to calculate correlation between variant dosages
					for (int i = 0; i < pVariantDosages.length; ++i) {
						regression.addData(pVariantDosages[i], qVariantDosages[i]);
					}

					cor[p][q] = regression.getR();
					cor[q][p] = cor[p][q];
					
					regression.clear();
				}
			}

			return new GenotypieCorrelationResult(cor, includedVariants);

		}

	}

}
