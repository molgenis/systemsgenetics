package umcg.genetica.features;

import umcg.genetica.enums.Chromosome;

/**
 * Created by hwestra on 12/3/15.
 */
public class BedGraphFeature extends Feature {

	private double value;
	private double pos = Double.NaN;
	private double neg = Double.NaN;

	public BedGraphFeature(Chromosome chr, Integer pos, int i) {
		super(chr, pos, i);
	}

	public void setValue(double posStrand, double negStrand) {
		this.pos = posStrand;
		this.neg = negStrand;
	}

	public double getValue() {
		return value;
	}

	public void setValue(double v) {
		this.value = v;
	}

	public double getPos() {
		return pos;
	}

	public double getNeg() {
		return neg;
	}

	public boolean isStranded(){
		return (!Double.isNaN(neg) && !Double.isNaN(pos));
	}
}
