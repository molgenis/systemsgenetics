package umcg.genetica.features;

import umcg.genetica.enums.Chromosome;

/**
 * Created by hwestra on 8/31/15.
 */
public class PeakFeature extends Feature {


	// 	abs_summit	pileup	-log10(pvalue)	fold_enrichment	-log10(qvalue)	name
	int summit;
	double pileup;
	double pval;
	double foldenrich;
	double qval;


	public PeakFeature(Chromosome chromosome, int alignmentStart, int alignmentEnd,
					   int summit, double pileup, double pval, double foldenrich, double qval) {
		this.chromosome = chromosome;
		this.start = alignmentStart;
		this.stop = alignmentEnd;

		this.summit = summit;
		this.pileup = pileup;

		this.pval = pval;

		this.foldenrich = foldenrich;
		this.qval = qval;
	}

	public int getSummit() {
		return summit;
	}

	public void setSummit(int summit) {
		this.summit = summit;
	}

	public double getPileup() {
		return pileup;
	}

	public void setPileup(double pileup) {
		this.pileup = pileup;
	}

	public double getPval() {
		return pval;
	}

	public void setPval(double pval) {
		this.pval = pval;
	}

	public double getFoldenrich() {
		return foldenrich;
	}

	public void setFoldenrich(double foldenrich) {
		this.foldenrich = foldenrich;
	}

	public double getQval() {
		return qval;
	}

	public void setQval(double qval) {
		this.qval = qval;
	}
}
