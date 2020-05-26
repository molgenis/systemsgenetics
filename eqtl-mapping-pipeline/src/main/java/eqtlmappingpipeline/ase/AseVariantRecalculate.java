package eqtlmappingpipeline.ase;

import cern.colt.list.tdouble.DoubleArrayList;
import cern.colt.list.tint.IntArrayList;
import cern.jet.stat.tdouble.Probability;
import java.util.ArrayList;
import org.apache.commons.math3.stat.inference.AlternativeHypothesis;
import org.apache.commons.math3.stat.inference.BinomialTest;
import org.apache.commons.math3.stat.regression.SimpleRegression;
import org.molgenis.genotype.Allele;
import org.molgenis.genotype.variant.id.GeneticVariantId;

/**
 *
 * @author Patrick Deelen
 */
public class AseVariantRecalculate implements AseVariant {

	private AseMleBeta mle = null;
	private static final double LARGEST_ZSCORE = Probability.normalInverse(Double.MIN_NORMAL);
	private final AseVariantBean originalAseVariant;

	public AseVariantRecalculate(AseVariantBean originalAseVariant) {
		this.originalAseVariant = originalAseVariant;
	}

	@Override
	public String getChr() {
		return originalAseVariant.getChr();
	}

	@Override
	public int getPos() {
		return originalAseVariant.getPos();
	}

	@Override
	public GeneticVariantId getId() {
		return originalAseVariant.getId();
	}

	@Override
	public Allele getA1() {
		return originalAseVariant.getA1();
	}

	@Override
	public Allele getA2() {
		return originalAseVariant.getA2();
	}

	@Override
	public IntArrayList getA1Counts() {
		return originalAseVariant.getA1Counts();
	}

	@Override
	public IntArrayList getA2Counts() {
		return originalAseVariant.getA2Counts();
	}

	@Override
	public void calculateStatistics() {

		mle = new AseMleBeta(originalAseVariant.getA1Counts(), originalAseVariant.getA2Counts());

	}

	@Override
	public double getMetaZscore() {
		return Double.NaN;
	}

	@Override
	public double getMetaPvalue() {
		return Double.NaN;
	}

	@Override
	public int compareTo(AseVariant o) {

		double thisRatioD = this.getMle().getRatioD();
		double otherRatioD = o.getMle().getRatioD();

		//Reverse compare. Largest first
		return Double.compare(otherRatioD, thisRatioD);

	}

	@Override
	public int getSampleCount() {
		return originalAseVariant.getSampleCount();
	}

	@Override
	public double getCountPearsonR() {
		return Double.NaN;
	}

	@Override
	public ArrayList<String> getSampleIds() {
		return originalAseVariant.getSampleIds();
	}

	@Override
	public DoubleArrayList getPValues() {
		return new DoubleArrayList();
	}

	@Override
	public AseMleBeta getMle() {
		if (mle == null) {
			calculateStatistics();
		}
		return mle;
	}

	@Override
	public double getLikelihoodRatioP() {
		return mle.getRatioP();
	}

	@Override
	public double getLikelihoodRatioD() {
		return mle.getRatioD();
	}

	@Override
	public double getEffect() {
		return mle.getMaxLikelihoodP();
	}
	
	public double getOriginalLikelihoodRatioP() {
		return originalAseVariant.getLikelihoodRatioP();
	}

	public double getOriginalLikelihoodRatioD() {
		return originalAseVariant.getLikelihoodRatioD();
	}
	
	public double getOriginalEffect() {
		return originalAseVariant.getEffect();
	}
	
	public String getGenes(){
		return originalAseVariant.getGenes();
	}
	
}
