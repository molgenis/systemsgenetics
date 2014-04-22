package eqtlmappingpipeline.ase;

import cern.colt.list.tdouble.DoubleArrayList;
import cern.colt.list.tint.IntArrayList;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.stat.inference.AlternativeHypothesis;
import org.apache.commons.math3.stat.inference.BinomialTest;
import org.molgenis.genotype.Allele;
import org.molgenis.genotype.variant.id.GeneticVariantId;

/**
 *
 * @author Patrick Deelen
 */
public class AseVariant {

	private final String chr;
	private final int pos;
	private final GeneticVariantId id;
	private final Allele a1;
	private final Allele a2;
	private final IntArrayList a1Counts;
	private final IntArrayList a2Counts;
	private final DoubleArrayList binomPvalues;
	private double metaZscore;
	private static final BinomialTest btest = new BinomialTest();
	private static final NormalDistribution normalDist = new NormalDistribution();

	public AseVariant(String chr, int pos, GeneticVariantId id, Allele a1, Allele a2) {
		this.chr = chr;
		this.pos = pos;
		this.id = id;
		this.a1 = a1;
		this.a2 = a2;
		this.a1Counts = new IntArrayList();
		this.a2Counts = new IntArrayList();
		this.binomPvalues = new DoubleArrayList();
		this.metaZscore = Double.NaN;
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

	public DoubleArrayList getBinomPvalues() {
		return binomPvalues;
	}

	public void calculateMetaZscore() {
		
		double zscoreSum = 0;
		
		for (int i = 0 ; i < binomPvalues.size() ; ++i){
			// we used 2 sided test so divide by 2
			double absZscore = normalDist.inverseCumulativeProbability(binomPvalues.getQuick(i)/2);
			
			// Min / plus might look counter intuative but i omit 1 - p/2 above so here I have to swap
			if(a1Counts.getQuick(i) < a2Counts.getQuick(i)){
				zscoreSum -= absZscore;
			} else {
				zscoreSum += absZscore;
			}
		}
		
		metaZscore = zscoreSum / Math.sqrt(binomPvalues.size());
		
		
	}

	public double getMetaZscore() {
		if(Double.isNaN(metaZscore)){
			calculateMetaZscore();
		}
		return metaZscore;
	}

	public synchronized void addCounts(int a1Count, int a2Count) {
		
		this.metaZscore = Double.NaN;//Reset meta Z-score when adding new data
		
		a1Counts.add(a1Count);
		a2Counts.add(a2Count);

		binomPvalues.add(btest.binomialTest(a1Count + a2Count, a1Count, 0.5, AlternativeHypothesis.TWO_SIDED));

	}
}
