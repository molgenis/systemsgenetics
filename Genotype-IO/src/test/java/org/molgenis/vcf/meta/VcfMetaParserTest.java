package org.molgenis.vcf.meta;

import static org.testng.Assert.assertEquals;

import java.io.IOException;
import java.io.StringReader;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.Map;

import org.testng.annotations.Test;

import com.google.common.collect.Lists;

public class VcfMetaParserTest
{
	private static String vcfStr = "##fileformat=VCFv4.1\n" +
			"##mycustomheader=mycustomvalue\n" +
			"##contig=<ID=1,length=249240621>\n" +
			"##contig=<ID=2,length=259240621>\n" +
			"##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">\n" + 
			"##FILTER=<ID=flt,Description=\"Failing one of the filters\">\n" + 
			"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" + 
			"#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO\n" + 
			"1	565286	rs1578391	C	T	.	flt	NS=1;DP=5;AF=1.000;ANNOT=INT;GI=LOC100131754\n";
	
	@Test
	public void parse() throws IOException
	{
		VcfMetaParser vcfMetaParser = new VcfMetaParser(new StringReader(vcfStr));
		VcfMeta vcfMeta = vcfMetaParser.parse();
		assertEquals(vcfMeta.getFileFormat(), "VCFv4.1");
		assertEquals(vcfMeta.get("mycustomheader"), "mycustomvalue");
		
		Map<String, String> properties1 = new LinkedHashMap<String, String>();
		properties1.put("ID", "1");
		properties1.put("length", "249240621");
		VcfMetaContig vcfMetaContig1 = new VcfMetaContig(properties1);
		Map<String, String> properties2 = new LinkedHashMap<String, String>();
		properties2.put("ID", "2");
		properties2.put("length", "259240621");
		VcfMetaContig vcfMetaContig2 = new VcfMetaContig(properties2);
		assertEquals(vcfMeta.get("mycustomheader"), "mycustomvalue");
		assertEquals(Lists.newArrayList(vcfMeta.getContigMeta()), Arrays.asList(vcfMetaContig1, vcfMetaContig2));
		
		Map<String, String> infoProperties = new LinkedHashMap<String, String>();
		infoProperties.put("ID", "NS");
		infoProperties.put("Number", "1");
		infoProperties.put("Type", "Integer");
		infoProperties.put("Description", "Number of Samples With Data");
		VcfMetaInfo vcfMetaInfo = new VcfMetaInfo(infoProperties );
		assertEquals(Lists.newArrayList(vcfMeta.getInfoMeta()), Arrays.asList(vcfMetaInfo));
		
		Map<String, String> filterProperties = new LinkedHashMap<String, String>();
		filterProperties.put("ID", "flt");
		filterProperties.put("Description", "Failing one of the filters");
		VcfMetaFilter vcfMetaFilter = new VcfMetaFilter(filterProperties);
		assertEquals(Lists.newArrayList(vcfMeta.getFilterMeta()), Arrays.asList(vcfMetaFilter));
		
		Map<String, String> formatProperties = new LinkedHashMap<String, String>();
		formatProperties.put("ID", "GT");
		formatProperties.put("Number", "1");
		formatProperties.put("Type", "String");
		formatProperties.put("Description", "Genotype");
		VcfMetaFormat vcfMetaFormat = new VcfMetaFormat(formatProperties);
		assertEquals(Lists.newArrayList(vcfMeta.getFormatMeta()), Arrays.asList(vcfMetaFormat));
		
		assertEquals(vcfMeta.getColNames(), new String[]{"#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"});
	}
}
