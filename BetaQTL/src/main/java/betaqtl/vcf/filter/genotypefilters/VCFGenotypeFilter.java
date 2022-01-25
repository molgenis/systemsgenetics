package betaqtl.vcf.filter.genotypefilters;

import betaqtl.vcf.VCFVariant;


/**
 * Created by hwestra on 2/8/16.
 */
public interface VCFGenotypeFilter {
	void filter(VCFVariant variant);
}
