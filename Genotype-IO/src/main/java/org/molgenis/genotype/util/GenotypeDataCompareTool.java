package org.molgenis.genotype.util;

import java.util.Iterator;

import org.molgenis.genotype.GenotypeData;
import org.molgenis.genotype.variant.GeneticVariant;

/**
 * Check if two GenotypeData or GeneticVariant objects contain then same data.
 * Be careful in using this, only use for testing as this can be very expensive
 * methods.
 *
 * @author erwin
 *
 */
public class GenotypeDataCompareTool {

	public static boolean same(GenotypeData gd1, GenotypeData gd2) {
		if (gd1.getSamples().equals(gd2.getSamples())) {
			//System.out.println("Samples same");
			//if (gd1.getSampleAnnotations().equals(gd2.getSampleAnnotations())) {
				//System.out.println("Sample annotations same");
				if (gd1.getVariantAnnotations().equals(gd2.getVariantAnnotations())) {
					//System.out.println("Sample variants same");

					Iterator<GeneticVariant> it1 = gd1.iterator();
					Iterator<GeneticVariant> it2 = gd2.iterator();

					while (it1.hasNext()) {
						if (!it2.hasNext()) {
							//System.err.println("Number variants not same");
							return false;
						}

						GeneticVariant v1 = it1.next();
						GeneticVariant v2 = it2.next();

						//System.out.println(v1.getSampleVariants());
						//System.out.println(v2.getSampleVariants());
						
						if (!same(v1, v2)) {
							//System.err.println("Variants not the same");
							return false;
						}

					}
					if (it2.hasNext()) {
						//System.err.println("Number variants not same");
						return false;
					} else {
						//System.err.println("Number variants the same");
						return true;
					}
				}
			//}
		}
		//System.err.println("Not the same");
		return false;
	}

	public static boolean same(GeneticVariant v1, GeneticVariant v2) {
		if (!v1.getVariantId().isSameId(v2.getVariantId())) {
			//System.err.println("Id not the same");
			return false;
		}
		if (!v1.getSequenceName().equals(v2.getSequenceName())) {
			//System.err.println("Seq name not same");
			return false;
		}

		if (v1.getStartPos() != v2.getStartPos()) {
			//System.err.println("Start pos not the same");
			return false;
		}

		if (!v1.getVariantAlleles().equals(v2.getVariantAlleles())) {
			//System.err.println("Variant alleles not same");
			return false;
		}

		if (!v1.getSampleVariants().equals(v2.getSampleVariants())) {
			//System.err.println("Sample variants not same");
			return false;
		}

		if (!v1.getAnnotationValues().equals(v1.getAnnotationValues())) {
			//System.err.println("Variant annotations not the same");
			return false;
		}

		return true;
	}
}
