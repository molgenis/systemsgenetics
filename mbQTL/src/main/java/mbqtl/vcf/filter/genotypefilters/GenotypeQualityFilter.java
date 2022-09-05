package mbqtl.vcf.filter.genotypefilters;

import cern.colt.matrix.tdouble.DoubleMatrix2D;
import mbqtl.vcf.VCFVariant;


/**
 * Created by hwestra on 2/8/16.
 */
public class GenotypeQualityFilter implements VCFGenotypeFilter {

	private int minimalGenotypeQual = 30;

	public GenotypeQualityFilter() {
	}

	public GenotypeQualityFilter(int minimalGenotypeQual) {
		this.minimalGenotypeQual = minimalGenotypeQual;
	}

	public void filter(VCFVariant variant) {
		DoubleMatrix2D alleles = variant.getGenotypeAllelesAsMatrix2D();
		short[] genotypeQuals = variant.getGenotypeQuals();
		if (genotypeQuals != null) {
			for (int i = 0; i < genotypeQuals.length; i++) {
				int qual = genotypeQuals[i];
				if (qual < minimalGenotypeQual) {
					alleles.set(i,0,-1);
					alleles.set(i,1,-1);
				}
			}
		} else {
			// should everything be -1?
		}
	}
}
