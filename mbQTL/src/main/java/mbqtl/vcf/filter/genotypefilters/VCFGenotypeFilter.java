package mbqtl.vcf.filter.genotypefilters;

import mbqtl.vcf.VCFVariant;


/**
 * Created by hwestra on 2/8/16.
 */
public interface VCFGenotypeFilter {
	void filter(VCFVariant variant);
}
