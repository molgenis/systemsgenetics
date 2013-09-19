package org.molgenis.genotype.util;

import java.util.HashMap;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;

import org.molgenis.genotype.Allele;
import org.molgenis.genotype.Alleles;

public class MafCalculator
{

	public static MafResult calculateMaf(Alleles alleles, Allele reference, List<Alleles> samplesAlleles)
	{

		HashMap<Allele, AtomicInteger> alleleCounts = new HashMap<Allele, AtomicInteger>(alleles.getAlleles().size());
		for (Allele allele : alleles.getAlleles())
		{
			alleleCounts.put(allele, new AtomicInteger());
		}

		for (Alleles sampleAlleles : samplesAlleles)
		{
			if (sampleAlleles != null)
			{
				for (Allele sampleAllele : sampleAlleles.getAlleles())
				{
					if (sampleAllele != null && sampleAllele != Allele.ZERO)
					{
						if (!alleleCounts.containsKey(sampleAllele))
						{
							throw new NullPointerException("No counter for allele: " + sampleAllele);
						}
						alleleCounts.get(sampleAllele).incrementAndGet();
					}
				}
			}
		}

		// This does not really do something since the first allele should
		// always be the reference allele in our genetic variants.
		Allele provisionalMinorAllele = reference != null ? reference : alleles.getAlleles().get(0);

		int provisionalMinorAlleleCount = alleleCounts.get(provisionalMinorAllele).get();
		int totalAlleleCount = 0;

		for (Allele allele : alleles.getAlleles())
		{

			int alleleCount = alleleCounts.get(allele).get();

			totalAlleleCount += alleleCount;

			if (alleleCount < provisionalMinorAlleleCount)
			{
				provisionalMinorAlleleCount = alleleCount;
				provisionalMinorAllele = allele;
			}
		}

		return new MafResult(provisionalMinorAllele, provisionalMinorAlleleCount / (float) totalAlleleCount);

	}
}
