package umcg.genetica.math.stats;


import java.util.ArrayList;

/*
This is an adaptation of https://github.com/secastel/aFC/blob/master/aFC.py
originally written by Pejman Mohammadi and Stephane Castel
 */
public class AllelicFoldChange {


	// expect un-centered genotypes as input (i.e. variable between 0,1,2 and -1 as missing)
	public void calculateAllelicFoldChange(int[] genotypesInput, double[] expValsInput) {
		int[] alleleCounts = new int[2];
		ArrayList<Double> genotypes = new ArrayList<>();
		ArrayList<Double> expvals = new ArrayList<>();

		double[] meansPerGenotypeGroup = new double[3];
		int[] indsPerGenotypeGroup = new int[3];

		for (int i = 0; i < genotypesInput.length; i++) {
			double gt = genotypesInput[i];
			if (genotypesInput[i] > -1) {
				double exp = expValsInput[i];
				genotypes.add(gt);
				expvals.add(exp);
				if (genotypesInput[i] == 0) {
					alleleCounts[0] += 2;
					meansPerGenotypeGroup[0] += exp;
					indsPerGenotypeGroup[0]++;
				} else if (genotypesInput[i] == 1) {
					alleleCounts[0]++;
					alleleCounts[1]++;
					meansPerGenotypeGroup[1] += exp;
					indsPerGenotypeGroup[1]++;
				} else {
					alleleCounts[1] += 2;
					meansPerGenotypeGroup[2] += exp;
					indsPerGenotypeGroup[2]++;
				}
			}
		}

		// calculate effect size
		// calculate mean per genotype group

		for (int i = 0; i < 3; i++) {
			meansPerGenotypeGroup[i] /= indsPerGenotypeGroup[i];
		}

		// calculate log2tratio
		double[] log2ratios = new double[3];
		for (int i = 0; i < 3; i++) {
			switch (i) {
				case 0:
					// bound_basic(p_m[2] - p_m[0], -args.ecap, args.ecap)

					break;
				case 1:
					// bound_basic(p_m[1] - p_m[2], -0.9999999, args.ecap)
					break;
				case 2:
					// bound_basic(p_m[1] - p_m[0], -1, args.ecap)
					break;
			}
		}




		/*
		for each sample, if non missing genotype:
		     count allele A and B --> put in alleleCounts
		     store number of alt alleles --> put in alleleCounts; also grep any covariates that we may have

		if allele cts > min nr alleles:
		     make matrix, store genotypes, phenotypes and also the covariates, again?
		     correct for covariates:
		          // we have facilities for this; no need to implement here.
		     calculate effect size:


		 */


	}

	private double bound_basic(double x, double l, double h) {

		

		return 0;
	}

	public void bootstrapBCi() {
		/*
		// grabbed from: https://github.com/cgevans/scikits-bootstrap/blob/master/scikits/bootstrap/bootstrap.py

		// start of bootstrap routine
		# The value of the statistic function applied just to the actual data.
        ostat = statfunction(*tdata)

        # The bias correction value.
        z0 = nppf( ( 1.0*np.sum(stat < ostat, axis=0)  ) / n_samples )

        # Statistics of the jackknife distribution
        jackindexes = jackknife_indexes(tdata[0])
        jstat = [statfunction(*(x[indexes] for x in tdata)) for indexes in jackindexes]
        jmean = np.mean(jstat,axis=0)

        # Temporarily kill numpy warnings:
        oldnperr = np.seterr(invalid='ignore')
        # Acceleration value
        a = np.sum((jmean - jstat)**3, axis=0) / (
            6.0 * np.sum((jmean - jstat)**2, axis=0)**1.5)
        if np.any(np.isnan(a)):
            nanind = np.nonzero(np.isnan(a))
            warnings.warn("BCa acceleration values for indexes {} were undefined. \
Statistic values were likely all equal. Affected CI will \
be inaccurate.".format(nanind), InstabilityWarning, stacklevel=2)

        zs = z0 + nppf(alphas).reshape(alphas.shape+(1,)*z0.ndim)

        avals = ncdf(z0 + zs/(1-a*zs))
        np.seterr(**oldnperr)


        nvals = np.round((n_samples-1)*avals)
        // end of bootstrap routine

    oldnperr = np.seterr(invalid='ignore')
    if np.any(np.isnan(nvals)):
        warnings.warn("Some values were NaN; results are probably unstable " +
                      "(all values were probably equal)", InstabilityWarning,
                      stacklevel=2)
    if np.any(nvals == 0) or np.any(nvals == n_samples-1):
        warnings.warn("Some values used extremal samples; " +
                      "results are probably unstable.",
                      InstabilityWarning, stacklevel=2)
    elif np.any(nvals < 10) or np.any(nvals >= n_samples-10):
        warnings.warn("Some values used top 10 low/high samples; " +
                      "results may be unstable.",
                      InstabilityWarning, stacklevel=2)
    np.seterr(**oldnperr)

    nvals = np.nan_to_num(nvals).astype('int')

    if output == 'lowhigh':
        if nvals.ndim == 1:
            # All nvals are the same. Simple broadcasting
            return stat[nvals]
        else:
            # Nvals are different for each data point. Not simple broadcasting.
            # Each set of nvals along axis 0 corresponds to the data at the same
            # point in other axes.
            return stat[(nvals, np.indices(nvals.shape)[1:].squeeze())]
    elif output == 'errorbar':
        if nvals.ndim == 1:
            return abs(statfunction(data)-stat[nvals])[np.newaxis].T
        else:
            return abs(statfunction(data)-stat[(nvals, np.indices(nvals.shape)[1:])])[np.newaxis].T
    else:
        raise ValueError("Output option {0} is not supported.".format(output))
		 */
	}

}
