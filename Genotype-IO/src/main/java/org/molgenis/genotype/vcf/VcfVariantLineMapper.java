package org.molgenis.genotype.vcf;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.lang3.StringUtils;
import org.molgenis.genotype.Allele;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variant.GeneticVariantMeta;
import org.molgenis.genotype.variant.ReadOnlyGeneticVariant;
import org.molgenis.genotype.variant.VariantLineMapper;
import org.molgenis.genotype.variant.sampleProvider.SampleVariantsProvider;
import org.molgenis.vcf.VcfInfo;
import org.molgenis.vcf.VcfRecord;
import org.molgenis.vcf.meta.VcfMeta;

public class VcfVariantLineMapper implements VariantLineMapper
{
	private final VcfMeta vcfMeta;
	private final SampleVariantsProvider sampleVariantsProvider;
	
	public VcfVariantLineMapper(VcfMeta vcfMeta, SampleVariantsProvider sampleVariantsProvider) {
		if(vcfMeta == null) throw new IllegalArgumentException("vcfMeta is null");
		if(sampleVariantsProvider == null) throw new IllegalArgumentException("sampleVariantsProvider is null");
		this.vcfMeta = vcfMeta;
		this.sampleVariantsProvider = sampleVariantsProvider;
	}
	
	@Override
	public GeneticVariant mapLine(String line)
	{
		VcfRecord vcfRecord = new VcfRecord(vcfMeta, StringUtils.split(line, '\t'));
		return toGeneticVariant(vcfRecord);
	}

	private GeneticVariant toGeneticVariant(VcfRecord vcfRecord) {
		List<String> identifiers = vcfRecord.getIdentifiers();
		int pos = vcfRecord.getPosition();
		String sequenceName = vcfRecord.getChromosome();			
		Allele refAllele = vcfRecord.getReferenceAllele();
		List<Allele> altAlleles = vcfRecord.getAlternateAlleles();
		
		Map<String, Object> annotationMap = new HashMap<String, Object>();
		for(VcfInfo vcfInfo : vcfRecord.getInformation())
			annotationMap.put(vcfInfo.getKey(), vcfInfo.getVal());
		
		Alleles alleles;
		if (altAlleles == null || altAlleles.isEmpty()) {
			alleles = Alleles.createAlleles(refAllele);
		} else {
			ArrayList<Allele> allelesList = new ArrayList<Allele>(altAlleles.size() + 1);
			allelesList.add(refAllele);
			allelesList.addAll(altAlleles);
			alleles = Alleles.createAlleles(allelesList);
		}

		GeneticVariantMeta geneticVariantMeta = new VcfGeneticVariantMeta(vcfMeta, Arrays.asList(vcfRecord.getFormat()));
		GeneticVariant geneticVariant = ReadOnlyGeneticVariant.createVariant(geneticVariantMeta, identifiers, pos, sequenceName, annotationMap, sampleVariantsProvider, alleles, refAllele);
		
		return geneticVariant;
	}
}
