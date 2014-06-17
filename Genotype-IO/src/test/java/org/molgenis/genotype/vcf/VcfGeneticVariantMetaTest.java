package org.molgenis.genotype.vcf;

import static org.testng.Assert.assertEquals;
import static org.testng.Assert.assertNull;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.molgenis.genotype.variant.GeneticVariantMeta;
import org.molgenis.vcf.VcfRecord;
import org.molgenis.vcf.meta.VcfMeta;
import org.molgenis.vcf.meta.VcfMetaFormat;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

import com.google.common.collect.Lists;

public class VcfGeneticVariantMetaTest
{
	private VcfGeneticVariantMeta vcfGeneticVariantMeta;

	@BeforeMethod
	public void setUp()
	{
		VcfMeta vcfMeta = new VcfMeta();

		Map<String, String> propertiesGt = new HashMap<String, String>();
		propertiesGt.put("ID", "GT");
		propertiesGt.put("Number", "1");
		propertiesGt.put("Type", "String");
		propertiesGt.put("Description", "Genotype");
		vcfMeta.addFormatMeta(new VcfMetaFormat(propertiesGt));

		Map<String, String> propertiesDp = new HashMap<String, String>();
		propertiesDp.put("ID", "DP");
		propertiesDp.put("Number", "1");
		propertiesDp.put("Type", "Integer");
		propertiesDp.put("Description", "Read Depth");
		vcfMeta.addFormatMeta(new VcfMetaFormat(propertiesDp));

		Map<String, String> propertiesEc = new HashMap<String, String>();
		propertiesEc.put("ID", "EC");
		propertiesEc.put("Number", "A");
		propertiesEc.put("Type", "Integer");
		propertiesEc.put("Description", "alternate allele counts");
		vcfMeta.addFormatMeta(new VcfMetaFormat(propertiesEc));

		Map<String, String> propertiesConfs = new HashMap<String, String>();
		propertiesConfs.put("ID", "CONFS");
		propertiesConfs.put("Number", "7");
		propertiesConfs.put("Type", "Float");
		propertiesConfs.put("Description", "NextGENe Confidence Scores");
		vcfMeta.addFormatMeta(new VcfMetaFormat(propertiesConfs));

		String[] tokens = new String[]
		{ "1", "565286", "rs1578391", "C", "T", ".", "flt", "NS=1;DP=5;AF=1.000;ANNOT=INT;GI=LOC100131754",
				"GT:DP:EC:CONFS", "1/1:5:5:5.300,5.300,1.000,1.000,1.000,1.000,1.000" };
		VcfRecord vcfRecord = new VcfRecord(vcfMeta, tokens);
		vcfGeneticVariantMeta = new VcfGeneticVariantMeta(vcfMeta, Arrays.asList(vcfRecord.getFormat()));
	}

	@Test(expectedExceptions = IllegalArgumentException.class)
	public void VcfGeneticVariantMeta()
	{
		new VcfGeneticVariantMeta(null, null);
	}

	@Test
	public void getRecordIds()
	{
		List<String> ids = Lists.newArrayList(vcfGeneticVariantMeta.getRecordIds());
		assertEquals(ids, Arrays.asList("GT", "DP", "EC", "CONFS"));
	}

	@Test
	public void getRecordType()
	{
		assertEquals(GeneticVariantMeta.Type.STRING, vcfGeneticVariantMeta.getRecordType("GT"));
		assertEquals(GeneticVariantMeta.Type.INTEGER, vcfGeneticVariantMeta.getRecordType("DP"));
		assertEquals(GeneticVariantMeta.Type.INTEGER_LIST, vcfGeneticVariantMeta.getRecordType("EC"));
		assertEquals(GeneticVariantMeta.Type.FLOAT_LIST, vcfGeneticVariantMeta.getRecordType("CONFS"));
		assertNull(vcfGeneticVariantMeta.getRecordType("invalid"));
	}
}
