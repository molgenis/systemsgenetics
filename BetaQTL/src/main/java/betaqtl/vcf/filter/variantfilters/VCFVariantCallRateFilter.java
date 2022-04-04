package betaqtl.vcf.filter.variantfilters;

import betaqtl.vcf.VCFVariant;

/**
 * Created by hwestra on 3/21/17.
 */
public class VCFVariantCallRateFilter implements VCFVariantFilter {

	private double threshold = 0.3;

	public VCFVariantCallRateFilter(double threshold) {
		this.threshold = threshold;
	}

	@Override
	public boolean passesThreshold(VCFVariant variant) {
		double cr = variant.getCallrate();
		return cr > threshold;
	}

	@Override
	public String toString() {
		return "VCFVariantCallRateFilter{" +
				"threshold=" + threshold +
				'}';
	}
}
