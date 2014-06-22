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
public class AseVariant implements Comparable<AseVariant>{

	private final String chr;
	private final int pos;
	private final GeneticVariantId id;
	private final Allele a1;
	private final Allele a2;
	private final IntArrayList a1Counts;
	private final IntArrayList a2Counts;
	private final DoubleArrayList a1MeanBaseQualities;
	private final DoubleArrayList a2MeanBaseQualities;
	private final DoubleArrayList pValues;
	private final ArrayList<String> sampleIds;
	private double metaZscore;
	private double metaPvalue;
	private double countPearsonR;
	private static final BinomialTest btest = new BinomialTest();
	private static final double LARGEST_ZSCORE = Probability.normalInverse(Double.MIN_NORMAL);

	public AseVariant(String chr, int pos, GeneticVariantId id, Allele a1, Allele a2) {
		this.chr = chr;
		this.pos = pos;
		this.id = id;
		this.a1 = a1;
		this.a2 = a2;
		this.a1Counts = new IntArrayList();
		this.a2Counts = new IntArrayList();
		this.pValues = new DoubleArrayList();
		this.a1MeanBaseQualities = new DoubleArrayList();
		this.a2MeanBaseQualities = new DoubleArrayList();
		this.sampleIds = new ArrayList<String>();
		this.metaZscore = Double.NaN;
		this.metaPvalue = Double.NaN;
		this.countPearsonR = Double.NaN;
	}

	public String getChr() {
		return chr;
	}

	public int getPos() {
		return pos;
	}

	public GeneticVariantId getId() {
		return id;
	}

	public Allele getA1() {
		return a1;
	}

	public Allele getA2() {
		return a2;
	}

	public IntArrayList getA1Counts() {
		return a1Counts;
	}

	public IntArrayList getA2Counts() {
		return a2Counts;
	}

	public DoubleArrayList getA1MeanBaseQualities() {
		return a1MeanBaseQualities;
	}

	public DoubleArrayList getA2MeanBaseQualities() {
		return a2MeanBaseQualities;
	}
	
	public void calculateStatistics() {

		double zscoreSum = 0;

		SimpleRegression regression = new SimpleRegression();

		for (int i = 0 ; i < a1Counts.size() ; ++i){

			regression.addData(a1Counts.getQuick(i), a2Counts.getQuick(i));

			final double pvalue = pValues.getQuick(i);

			// we used 2 sided test so divide by 2
			//double zscore = normalDist.inverseCumulativeProbability(pvalue/2);
			final double pvalueDiv2 = pvalue / 2;
			final double zscore;
			if (pvalueDiv2 < Double.MIN_NORMAL){
				zscore = LARGEST_ZSCORE;	
			} else {
				zscore = Probability.normalInverse(pvalueDiv2);
			}
			// Min / plus might look counter intuative but i omit 1 - p/2 above so here I have to swap
			if(a1Counts.getQuick(i) < a2Counts.getQuick(i)){
				zscoreSum -= zscore;
			} else {
				zscoreSum += zscore;
			}
		}

		countPearsonR = regression.getR();
		metaZscore = zscoreSum / Math.sqrt(a1Counts.size());
		metaPvalue = 2 * Probability.normal(-Math.abs(metaZscore));


	}
	public double getMetaZscore() {
		if(Double.isNaN(metaZscore)){
			calculateStatistics();
		}
		return metaZscore;
	}

	public double getMetaPvalue() {
		if(Double.isNaN(metaZscore)){
			calculateStatistics();
		}
		return metaPvalue;
	}

	public synchronized void addCounts(int a1Count, int a2Count, String sampleId) {
		addCounts(a1Count, a2Count, sampleId, Double.NaN, Double.NaN);
	}
	
	public synchronized void addCounts(int a1Count, int a2Count, String sampleId, double a1MeanBaseQuality, double a2MeanBaseQuality) {

		this.metaZscore = Double.NaN;//Reset meta Z-score when adding new data
		this.metaPvalue = Double.NaN;
		this.countPearsonR = Double.NaN;

		a1Counts.add(a1Count);
		a2Counts.add(a2Count);
		
		pValues.add(btest.binomialTest(a1Count + a2Count, a1Count, 0.5, AlternativeHypothesis.TWO_SIDED));
		
		a1MeanBaseQualities.add(a1MeanBaseQuality);
		a2MeanBaseQualities.add(a2MeanBaseQuality);
		
		sampleIds.add(sampleId);

	}

	@Override
	public int compareTo(AseVariant o) {

		double thisZ = Math.abs(this.getMetaZscore());
		double otherZ = Math.abs(o.getMetaZscore());

		if(thisZ < otherZ){
			return 1;
		} else if(thisZ == otherZ){
			return 0;
		} else{
			return -1;
		}

	}

	public int getSampleCount() {
		return a1Counts.size();
	}

	public double getCountPearsonR() {
		if(Double.isNaN(countPearsonR)){
			calculateStatistics();
		}
		return countPearsonR;
	}

	public ArrayList<String> getSampleIds() {
		return sampleIds;
	}

	public DoubleArrayList getPValues() {
		return pValues;
	}
	
}