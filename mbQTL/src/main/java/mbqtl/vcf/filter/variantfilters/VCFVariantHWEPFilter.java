package mbqtl.vcf.filter.variantfilters;


import mbqtl.vcf.VCFVariant;

/**
 * Created by hwestra on 3/21/17.
 */
public class VCFVariantHWEPFilter implements VCFVariantFilter {
	private double threshold = 1E-4;

	public enum MODE {
		CASES,
		CONTROLS,
		OVERALL
	}

	private MODE m;

	public VCFVariantHWEPFilter(double threshold) {
		this(threshold, MODE.OVERALL);
	}

	public VCFVariantHWEPFilter(double threshold, MODE m) {
		this.m = m;
		this.threshold = threshold;
	}

	@Override
	public boolean passesThreshold(VCFVariant variant) {
		if (m.equals(MODE.CASES)) {
			return variant.getHwepCases() > threshold;
		} else if (m.equals(MODE.CONTROLS)) {
			return variant.getHwepControls() > threshold;
		} else {
			return variant.getHwep() > threshold;
		}
	}

	@Override
	public String toString() {
		return "VCFVariantHWEPFilter{" +
				"threshold=" + threshold +
				", m=" + m +
				'}';
	}
}
