package org.molgenis.genotype.impute2;

import java.util.Arrays;
import java.util.List;

import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variant.ReadOnlyGeneticVariant;
import org.molgenis.genotype.variant.VariantLineMapper;

/**
 * Parses a line of a Impute2 haps file
 * 
 * @author erwin
 * 
 */
public class Impute2VariantLineMapper implements VariantLineMapper
{

	@Override
	public GeneticVariant mapLine(String line)
	{
		HapsEntry entry = HapsLineParser.parse(line);
		String sequence = getSequence(entry);
		List<String> alleles = Arrays.asList(entry.getFirstAllele(), entry.getSecondAllele());

		return ReadOnlyGeneticVariant.createVariant(entry.getSnpId(), entry.getPosition(), sequence,
				new HapsEntrySampleVariantsProvider(entry), alleles);

	}

	private String getSequence(HapsEntry entry)
	{
		switch (entry.getChromosomeNumber())
		{
			case 23:
				return "x";
			case 24:
				return "y";
			default:
				return Integer.toString(entry.getChromosomeNumber());
		}
	}
}
