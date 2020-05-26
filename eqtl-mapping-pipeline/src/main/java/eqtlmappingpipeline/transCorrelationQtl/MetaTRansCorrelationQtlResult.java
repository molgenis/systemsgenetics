/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.transCorrelationQtl;

import cern.jet.stat.tdouble.Probability;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import java.util.ArrayList;

/**
 *
 * @author patri
 */
public class MetaTRansCorrelationQtlResult {

	private final String transSnp;
	private final double qtlZscore;
	private final String localGene;
	private final String transGene;
	private final String assessedAllele;
	private final ArrayList<String> cohorts;
	private final TDoubleList perCohortInteractionZscore;
	private final TIntArrayList perCohortSampleCount;

	public MetaTRansCorrelationQtlResult(String transSnp, double qtlZscore, String localGene, String transGene, String assessedAllele) {
		this.transSnp = transSnp;
		this.qtlZscore = qtlZscore;
		this.localGene = localGene;
		this.transGene = transGene;
		this.assessedAllele = assessedAllele;
		this.cohorts = new ArrayList<>();
		this.perCohortInteractionZscore = new TDoubleArrayList();
		this.perCohortSampleCount = new TIntArrayList();
	}

	public void addCohortRestult(String allele, int sampleCount, double zscore, String cohort) {

		if (!allele.equals(assessedAllele)) {
			zscore = zscore * -1;
		}
		cohorts.add(cohort);
		perCohortInteractionZscore.add(zscore);
		perCohortSampleCount.add(sampleCount);

	}

	public String getTransSnp() {
		return transSnp;
	}

	public String getLocalGene() {
		return localGene;
	}

	public String getTransGene() {
		return transGene;
	}

	public int getTotalSampleCount() {
		return perCohortSampleCount.sum();
	}

	public double getZscore() {
		return qtlZscore;
	}

	public double calculateMetaZscore() {

		double numerator = 0;
		double denominator = 0;

		for (int i = 0; i < perCohortInteractionZscore.size(); ++i) {
			int sampleCount = perCohortSampleCount.get(i);
			double zscore = perCohortInteractionZscore.get(i);

			numerator += sampleCount * zscore;
			denominator += sampleCount * sampleCount;

		}

		denominator = Math.sqrt(denominator);

		return numerator / denominator;

	}

	public double calculateMetaP() {
		return 2 * Probability.normal(-Math.abs(calculateMetaZscore()));
	}

	public String getCohorts() {
		return String.join(",", cohorts);
	}

	public String getCohortSizes() {

		StringBuilder cohortSizes = new StringBuilder();

		for (int i = 0; i < perCohortSampleCount.size(); ++i) {

			if (i > 0) {
				cohortSizes.append(",");
			}
			cohortSizes.append(String.valueOf(perCohortSampleCount.get(i)));

		}

		return cohortSizes.toString();

	}

	public String getCohortZscores() {

		StringBuilder cohortZscores = new StringBuilder();

		for (int i = 0; i < perCohortInteractionZscore.size(); ++i) {

			if (i > 0) {
				cohortZscores.append(",");
			}
			cohortZscores.append(String.valueOf(perCohortInteractionZscore.get(i)));

		}

		return cohortZscores.toString();

	}

}
