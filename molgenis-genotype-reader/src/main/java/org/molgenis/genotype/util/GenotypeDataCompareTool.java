package org.molgenis.genotype.util;

import java.util.Iterator;

import org.molgenis.genotype.GenotypeData;
import org.molgenis.genotype.variant.GeneticVariant;

/**
 * Check if two GenotypeData or GeneticVariant objects contain then same data. Be careful in using this, only use for
 * testing as this can be very expensive methods.
 * 
 * @author erwin
 * 
 */
public class GenotypeDataCompareTool
{
	public static boolean same(GenotypeData gd1, GenotypeData gd2)
	{
		if (gd1.getSamples().equals(gd2.getSamples()))
		{
			if (gd1.getSampleAnnotations().equals(gd2.getSampleAnnotations()))
			{
				if (gd1.getVariantAnnotations().equals(gd2.getVariantAnnotations()))
				{
					Iterator<GeneticVariant> it1 = gd1.iterator();
					Iterator<GeneticVariant> it2 = gd2.iterator();

					while (it1.hasNext())
					{
						if (!it2.hasNext())
						{
							return false;
						}

						GeneticVariant v1 = it1.next();
						GeneticVariant v2 = it2.next();

						if (!same(v1, v2))
						{
							return false;
						}

					}

					return !it2.hasNext();
				}
			}
		}

		return false;
	}

	public static boolean same(GeneticVariant v1, GeneticVariant v2)
	{
		if (!v1.getVariantId().isSameId(v2.getVariantId()))
		{
			return false;
		}

		if (!v1.getSequenceName().equals(v2.getSequenceName()))
		{
			return false;
		}

		if (v1.getStartPos() != v2.getStartPos())
		{
			return false;
		}

		if (!v1.getVariantAlleles().equals(v2.getVariantAlleles()))
		{
			return false;
		}

		if (!v1.getSampleVariants().equals(v2.getSampleVariants()))
		{
			return false;
		}

		if (!v1.getAnnotationValues().equals(v1.getAnnotationValues()))
		{
			return false;
		}

		return true;
	}
}
