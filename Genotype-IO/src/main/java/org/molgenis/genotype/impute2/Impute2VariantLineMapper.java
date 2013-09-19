package org.molgenis.genotype.impute2;

import java.util.Arrays;
import java.util.List;
import java.util.StringTokenizer;

import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variant.ReadOnlyGeneticVariant;
import org.molgenis.genotype.variant.VariantLineMapper;
import org.molgenis.genotype.variant.sampleProvider.SampleVariantsProvider;

/**
 * Parses a line of a Impute2 haps file that is tab separated for tabix
 * 
 * @author erwin fixed by Patrick
 * 
 */
public class Impute2VariantLineMapper implements VariantLineMapper
{

	private final SampleVariantsProvider sampleVariantProvider;

	public Impute2VariantLineMapper(SampleVariantsProvider sampleVariantProvider)
	{
		this.sampleVariantProvider = sampleVariantProvider;
	}

	@Override
	public GeneticVariant mapLine(String line)
	{

		StringTokenizer tokenizer = new StringTokenizer(line, "\t");
		String chrom = tokenizer.nextToken();
		String snpId = tokenizer.nextToken();
		int position = Integer.parseInt(tokenizer.nextToken());
		String firstAllele = tokenizer.nextToken();
		String secondAllele = tokenizer.nextToken();

		List<String> alleles = Arrays.asList(firstAllele, secondAllele);

		return ReadOnlyGeneticVariant.createVariant(snpId, position, chrom, sampleVariantProvider, alleles);

	}

}
