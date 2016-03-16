package org.molgenis.genotype.util;

import gnu.trove.map.hash.TObjectIntHashMap;
import java.util.List;

import org.molgenis.genotype.Allele;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.GenotypeDataException;

public class MafCalculator {

	public static MafResult calculateMaf(Alleles alleles, Allele reference, List<Alleles> samplesAlleles) {

		if (alleles.getAlleles().isEmpty() || alleles.getAlleles().size()==1) {
			return new MafResult(Allele.ZERO, 0);
		}

		TObjectIntHashMap<Allele> alleleCounts = new TObjectIntHashMap<Allele>(alleles.getAlleleCount());


		for (Allele allele : alleles.getAlleles()) {
			if(allele != Allele.ZERO){
				alleleCounts.put(allele, 0);
			}
		}
		
		if (alleleCounts.isEmpty()) {
			return new MafResult(Allele.ZERO, 0);
		}

		for (Alleles sampleAlleles : samplesAlleles) {
			if (sampleAlleles != null) {
				for (Allele sampleAllele : sampleAlleles.getAlleles()) {
					if (sampleAllele != null && sampleAllele != Allele.ZERO) {
						if (!alleleCounts.increment(sampleAllele)) {
							throw new GenotypeDataException("No counter for allele: " + sampleAllele + " expected one of the following alleles: " + alleles);
						}
					}
				}
			}
		}

		// This does not really do something since the first allele should
		// always be the reference allele in our genetic variants.
		Allele provisionalMinorAllele = reference != null ? reference : alleles.getAlleles().get(0);

		int provisionalMinorAlleleCount = alleleCounts.get(provisionalMinorAllele);
		int totalAlleleCount = 0;

		for (Allele allele : alleles.getAlleles()) {

			int alleleCount = alleleCounts.get(allele);

			totalAlleleCount += alleleCount;

			if (alleleCount < provisionalMinorAlleleCount) {
				provisionalMinorAlleleCount = alleleCount;
				provisionalMinorAllele = allele;
			}
		}

		return new MafResult(provisionalMinorAllele, provisionalMinorAlleleCount / (float) totalAlleleCount);

	}
}
