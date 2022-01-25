package betaqtl.vcf.filter.genotypefilters;

import cern.colt.matrix.tdouble.DoubleMatrix2D;
import betaqtl.vcf.VCFVariant;


/**
 * Created by hwestra on 2/8/16.
 */
public class AllelicDepthFilter implements VCFGenotypeFilter {

	private int minimumReadDepth;
	private double abCutoff = 0.15;

	public AllelicDepthFilter() { }

	public AllelicDepthFilter(double allelicBalanceCutoff, int minimumReadDepth) {
		this.abCutoff = allelicBalanceCutoff;
		this.minimumReadDepth = minimumReadDepth;
	}

	public void filter(VCFVariant variant) {

		short[][] allelicDepth = variant.getAllelicDepth();
		if (allelicDepth != null) {
			DoubleMatrix2D alleles = variant.getGenotypeAllelesAsMatrix2D();
			for (int i = 0; i < variant.getNrSamples(); i++) {
				if (alleles.getQuick(i, 0) != -1) {
					int sum = 0;
					for (int j = 0; j < allelicDepth.length; j++) {
						sum += allelicDepth[j][i];
					}
					if (alleles.getQuick(i, 0) != alleles.getQuick(i, 1)) {

						double ab1 = (double) allelicDepth[(byte) alleles.getQuick(i, 0)][i] / sum;
						double ab2 = (double) allelicDepth[(byte) alleles.getQuick(i, 1)][i] / sum;
						boolean ok = true;
						if (ab1 < abCutoff || ab1 > (1 - abCutoff)) {
							ok = false;
						}
						if (ab2 < abCutoff || ab2 > (1 - abCutoff)) {
							ok = false;
						}

						if (!ok || sum < minimumReadDepth) {
							alleles.setQuick(i, 0, -1);
							alleles.setQuick(i, 1, -1);
						}
					} else {

						if (allelicDepth[(byte) alleles.getQuick(i, 0)][i] < minimumReadDepth) {
							alleles.setQuick(i, 0, -1);
							alleles.setQuick(i, 1, -1);
						}
					}
				}
			}
		}
	}
}
