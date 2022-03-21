package mbqtl.vcf.filter.variantfilters;

import mbqtl.vcf.VCFVariant;
import umcg.genetica.features.Feature;

import java.util.ArrayList;

/**
 * Created by hwestra on 4/26/17.
 */
public class VCFVariantRegionFilter implements VCFVariantFilter {

	private final ArrayList<Feature> regions;

	public VCFVariantRegionFilter(ArrayList<Feature> regions) {
		this.regions = regions;
	}

	@Override
	public String toString() {
		return "VCFVariantRegionFilter{" +
				"regions=" + regions.size() +
				'}';
	}

	@Override
	public boolean passesThreshold(VCFVariant variant) {
		for (Feature region : regions) {
			if (variant.asFeature().overlaps(region)) {
				return true;
			}
		}
		return false;
	}
}
