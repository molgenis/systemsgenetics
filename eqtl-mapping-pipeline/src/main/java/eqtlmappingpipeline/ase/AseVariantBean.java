package eqtlmappingpipeline.ase;

import cern.colt.list.tdouble.DoubleArrayList;
import cern.colt.list.tint.IntArrayList;
import java.util.ArrayList;
import org.molgenis.genotype.Allele;
import org.molgenis.genotype.variant.id.GeneticVariantId;

/**
 *
 * @author Patrick Deelen
 */
public class AseVariantBean implements AseVariant{

	private final String chr;
	private final int pos;
	private final GeneticVariantId id;
	private final Allele a1;
	private final Allele a2;
	private final IntArrayList a1Counts;
	private final IntArrayList a2Counts;
	private final DoubleArrayList pValues;
	private final ArrayList<String> sampleIds;
	private final double metaZscore;
	private final double metaPvalue;
	private final double countPearsonR;
	private AseMle mle;

	public AseVariantBean(String chr, int pos, GeneticVariantId id, Allele a1, Allele a2, IntArrayList a1Counts, IntArrayList a2Counts, DoubleArrayList pValues, ArrayList<String> sampleIds, double metaZscore, double metaPvalue, double countPearsonR) {
		this.chr = chr;
		this.pos = pos;
		this.id = id;
		this.a1 = a1;
		this.a2 = a2;
		this.a1Counts = a1Counts;
		this.a2Counts = a2Counts;
		this.pValues = pValues;
		this.sampleIds = sampleIds;
		this.metaZscore = metaZscore;
		this.metaPvalue = metaPvalue;
		this.countPearsonR = countPearsonR;
	}

	@Override
	public String getChr() {
		return chr;
	}

	@Override
	public int getPos() {
		return pos;
	}

	@Override
	public GeneticVariantId getId() {
		return id;
	}

	@Override
	public Allele getA1() {
		return a1;
	}

	@Override
	public Allele getA2() {
		return a2;
	}

	@Override
	public IntArrayList getA1Counts() {
		return a1Counts;
	}

	@Override
	public IntArrayList getA2Counts() {
		return a2Counts;
	}

	@Override
	public DoubleArrayList getPValues() {
		return pValues;
	}

	@Override
	public ArrayList<String> getSampleIds() {
		return sampleIds;
	}

	@Override
	public double getMetaZscore() {
		return metaZscore;
	}

	@Override
	public double getMetaPvalue() {
		return metaPvalue;
	}

	@Override
	public double getCountPearsonR() {
		return countPearsonR;
	}

	@Override
	public AseMle getMle() {
		return mle;
	}

	@Override
	public int getSampleCount() {
		return a1Counts.size();
	}

	@Override
	public int compareTo(AseVariant o) {
		double thisRatioD = this.getMle().getRatioD();
		double otherRatioD = o.getMle().getRatioD();
		
		//Reverse compare. Largest first
		return Double.compare(otherRatioD, thisRatioD);
	}
	
	
	
}
