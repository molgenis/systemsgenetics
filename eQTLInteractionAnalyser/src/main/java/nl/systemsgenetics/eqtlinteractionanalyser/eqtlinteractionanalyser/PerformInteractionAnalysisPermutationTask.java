/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.eqtlinteractionanalyser.eqtlinteractionanalyser;

import gnu.trove.set.hash.TIntHashSet;
import java.util.concurrent.Callable;
import org.apache.commons.math3.linear.SingularMatrixException;
import org.apache.commons.math3.stat.regression.SimpleRegression;

/**
 *
 * @author lude
 */
public class PerformInteractionAnalysisPermutationTask implements Callable<DoubleArrayIntegerObject> {

	public ExpressionDataset datasetGenotypes;
	public ExpressionDataset datasetExpression;
	public ExpressionDataset datasetCovariates;
	ExpressionDataset datasetCovariatesPCAForceNormal;
	public int covToTest = -1;
	public int nrSamples = -1;
	public org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression regression = null;
	public cern.jet.random.tdouble.StudentT tDistColt = null;
	private final SkippedInteractionTracker skippedTracker;
	private final SkippedInteractionWriter skippedWriter;
	private final TIntHashSet snpsToTest;

	public PerformInteractionAnalysisPermutationTask(ExpressionDataset datasetGenotypes, ExpressionDataset datasetExpression, ExpressionDataset datasetCovariates, ExpressionDataset datasetCovariatesPCAForceNormal, int covToTest, SkippedInteractionWriter skippedWriter, final TIntHashSet snpsToTest) {
		this.datasetGenotypes = datasetGenotypes;
		this.datasetExpression = datasetExpression;
		this.datasetCovariates = datasetCovariates;
		this.datasetCovariatesPCAForceNormal = datasetCovariatesPCAForceNormal;
		this.covToTest = covToTest;
		this.nrSamples = datasetGenotypes.nrSamples;
		this.skippedTracker = new SkippedInteractionTracker(datasetCovariates.probeNames[covToTest]);
		this.skippedWriter = skippedWriter;
		this.snpsToTest = snpsToTest;

		this.regression = new org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression();
		cern.jet.random.tdouble.engine.DoubleRandomEngine randomEngine = new cern.jet.random.tdouble.engine.DRand();
		this.tDistColt = new cern.jet.random.tdouble.StudentT(this.nrSamples - 4, randomEngine);

	}

	@Override
	public DoubleArrayIntegerObject call() throws Exception {
		double corrPvalueThreshold = 0.0001;
		double[] zScores = new double[datasetGenotypes.nrProbes];
		for (int snp = 0; snp < datasetGenotypes.nrProbes; snp++) {

			if(snpsToTest != null && !snpsToTest.contains(snp)){
				continue;
			}
			
			double corrPvalue = correlateCovariateWithGenotype(snp);
			if (corrPvalue > corrPvalueThreshold) { // don't compute the interaction if the covariate expression is affected by this SNP
				try {

					double[][] valsX = new double[nrSamples][3];
					for (int s = 0; s < nrSamples; s++) {
						valsX[s][0] = datasetGenotypes.rawData[snp][s];
						valsX[s][1] = datasetCovariates.rawData[covToTest][s];
						valsX[s][2] = valsX[s][0] * valsX[s][1];
					}
					double[] valsY = datasetExpression.rawData[snp];
					regression.newSampleData(valsY, valsX);
					double betaInteraction = regression.estimateRegressionParameters()[3];
					double seInteraction = regression.estimateRegressionParametersStandardErrors()[3];
					double tInteraction = betaInteraction / seInteraction;
					double pValueInteraction = 1;
					double zScoreInteraction = 0;
					if (tInteraction < 0) {
						pValueInteraction = tDistColt.cdf(tInteraction);
						if (pValueInteraction < 2.0E-323) {
							pValueInteraction = 2.0E-323;
						}
						zScoreInteraction = cern.jet.stat.tdouble.Probability.normalInverse(pValueInteraction);
					} else {
						pValueInteraction = tDistColt.cdf(-tInteraction);
						if (pValueInteraction < 2.0E-323) {
							pValueInteraction = 2.0E-323;
						}
						zScoreInteraction = -cern.jet.stat.tdouble.Probability.normalInverse(pValueInteraction);
					}
					zScores[snp] = zScoreInteraction;
				} catch (SingularMatrixException e) {
					zScores[snp] = 0;
					skippedTracker.addSkipped(SkippedInteractionTracker.Reason.SINGULAR, datasetGenotypes.probeNames[snp]);
				}
			}
			else{
				//System.out.println("Removing covariate because of eQTL effect! " + datasetCovariatesPCAForceNormal.probeNames[covToTest] + " : " + datasetGenotypes.probeNames[snp]);
				skippedTracker.addSkipped(SkippedInteractionTracker.Reason.SHARED_QTL, datasetGenotypes.probeNames[snp]);
				zScores[snp] = 0;
			}

		}
		skippedWriter.add(skippedTracker);
		return new DoubleArrayIntegerObject(zScores, covToTest);
	}

	private double correlateCovariateWithGenotype(int snp){
		SimpleRegression simpleRegression = new SimpleRegression();
		double[] expression = datasetCovariatesPCAForceNormal.rawData[covToTest];
		double[] genotypes = datasetGenotypes.rawData[snp];
		for (int s = 0; s < expression.length; s++) {
			simpleRegression.addData(expression[s], genotypes[s]);
		}
		//This is not working now that we have the _rs next to the gene names
//		if (datasetGenotypes.probeNames[snp].equals(datasetCovariatesPCAForceNormal.probeNames[covToTest])){
//			System.out.println("Same gene! " + datasetGenotypes.probeNames[snp] + "\t" + datasetCovariatesPCAForceNormal.probeNames[covToTest] + "\t" + simpleRegression.getSignificance() + "\t" + simpleRegression.getR());
//		}
		return simpleRegression.getSignificance();
	}
}
