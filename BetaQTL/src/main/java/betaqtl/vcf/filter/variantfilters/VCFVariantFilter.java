package betaqtl.vcf.filter.variantfilters;


import betaqtl.vcf.VCFVariant;

/**
 * Created by hwestra on 3/21/17.
 */
public interface VCFVariantFilter {
	boolean passesThreshold(VCFVariant variant);
}
