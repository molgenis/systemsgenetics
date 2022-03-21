package mbqtl.vcf.filter.variantfilters;


import mbqtl.vcf.VCFVariant;

/**
 * Created by hwestra on 3/21/17.
 */
public interface VCFVariantFilter {
	boolean passesThreshold(VCFVariant variant);
}
