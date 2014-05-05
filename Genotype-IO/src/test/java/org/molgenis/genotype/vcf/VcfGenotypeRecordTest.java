package org.molgenis.genotype.vcf;

import static org.testng.Assert.assertEquals;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.molgenis.genotype.Allele;
import org.molgenis.vcf.VcfRecord;
import org.molgenis.vcf.VcfSample;
import org.molgenis.vcf.meta.VcfMeta;
import org.molgenis.vcf.meta.VcfMetaFormat;
import org.testng.annotations.Test;

public class VcfGenotypeRecordTest
{

	@Test(expectedExceptions = IllegalArgumentException.class)
	public void VcfGenotypeRecord()
	{
		new VcfGenotypeRecord(null, null, null);
	}
	
	@Test
	public void getGenotypeRecordData()
	{
		// ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
		Map<String, String> propertiesGt = new HashMap<String, String>();
		propertiesGt.put("ID", "GT");
		propertiesGt.put("Number", "1");
		propertiesGt.put("Type", "String");
		propertiesGt.put("Description", "Genotype");
		VcfMetaFormat vcfMetaFormatGt = new VcfMetaFormat(propertiesGt);
		
		// ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
		Map<String, String> propertiesDp = new HashMap<String, String>();
		propertiesDp.put("ID", "DP");
		propertiesDp.put("Number", "1");
		propertiesDp.put("Type", "Integer");
		propertiesDp.put("Description", "Read Depth");
		VcfMetaFormat vcfMetaFormatDp = new VcfMetaFormat(propertiesDp);
		
		// ##FORMAT=<ID=EC,Number=A,Type=Integer,Description="alternate allele counts">
		Map<String, String> propertiesEc = new HashMap<String, String>();
		propertiesEc.put("ID", "EC");
		propertiesEc.put("Number", "A");
		propertiesEc.put("Type", "Integer");
		propertiesEc.put("Description", "alternate allele counts");
		VcfMetaFormat vcfMetaFormatEc = new VcfMetaFormat(propertiesEc);
		
		// ##FORMAT=<ID=CONFS,Number=3,Type=Float,Description="confidence scores">
		Map<String, String> propertiesConfs = new HashMap<String, String>();
		propertiesConfs.put("ID", "CONFS");
		propertiesConfs.put("Number", "3");
		propertiesConfs.put("Type", "Float");
		propertiesConfs.put("Description", "confidence scores");
		VcfMetaFormat vcfMetaFormatConfs = new VcfMetaFormat(propertiesConfs);
		
		VcfMeta vcfMeta = new VcfMeta();
		vcfMeta.addFormatMeta(vcfMetaFormatGt);
		vcfMeta.addFormatMeta(vcfMetaFormatDp);
		vcfMeta.addFormatMeta(vcfMetaFormatEc);
		vcfMeta.addFormatMeta(vcfMetaFormatConfs);
		
		String[] recordTokens = new String[]{"x", "x", "x", "x", "x", "x", "x", "x", "GT:DP:EC:CONFS"};
		VcfRecord vcfRecord = new VcfRecord(vcfMeta, recordTokens);
		String[] sampleTokens = new String[]{"0/2", "5", "5", "5.300,5.300,1.000"};
		VcfSample vcfSample = new VcfSample(vcfRecord, sampleTokens );
		VcfGenotypeRecord vcfGenotypeRecord = new VcfGenotypeRecord(vcfMeta, vcfRecord, vcfSample);
		assertEquals(vcfGenotypeRecord.getGenotypeRecordData("GT"), "0/2");
		assertEquals(vcfGenotypeRecord.getGenotypeRecordData("DP"), Integer.valueOf(5));
		assertEquals(vcfGenotypeRecord.getGenotypeRecordData("EC"), Arrays.asList(Integer.valueOf(5)));
		assertEquals(vcfGenotypeRecord.getGenotypeRecordData("CONFS"), Arrays.asList(Double.valueOf(5.300), Double.valueOf(5.300), Double.valueOf(1.000)));
	}
	
	@Test
	public void getSampleAlleles()
	{
		// ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
		Map<String, String> propertiesGt = new HashMap<String, String>();
		propertiesGt.put("ID", "GT");
		propertiesGt.put("Number", "1");
		propertiesGt.put("Type", "String");
		propertiesGt.put("Description", "Genotype");
		VcfMetaFormat vcfMetaFormatGt = new VcfMetaFormat(propertiesGt);
		
		VcfMeta vcfMeta = new VcfMeta();
		vcfMeta.addFormatMeta(vcfMetaFormatGt);
		
		String[] recordTokens = new String[]{"x", "x", "x", "G", "A", "x", "x", "x", "GT:DP:EC:CONFS"};
		VcfRecord vcfRecord = new VcfRecord(vcfMeta, recordTokens);
		String[] sampleTokens = new String[]{"1/0", "5", "5", "5.300,5.300,1.000"};
		VcfSample vcfSample = new VcfSample(vcfRecord, sampleTokens );
		VcfGenotypeRecord vcfGenotypeRecord = new VcfGenotypeRecord(vcfMeta, vcfRecord, vcfSample);
		List<Allele> alleles = vcfGenotypeRecord.getSampleAlleles().getAlleles();
		
		assertEquals(alleles.get(0).getAlleleAsSnp(), 'A');
		assertEquals(alleles.get(1).getAlleleAsSnp(), 'G');
	}
}
