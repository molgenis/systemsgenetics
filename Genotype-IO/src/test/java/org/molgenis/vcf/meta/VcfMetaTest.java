package org.molgenis.vcf.meta;

import static org.testng.Assert.assertEquals;
import static org.testng.Assert.assertNull;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import org.testng.annotations.Test;

import com.google.common.collect.Lists;

public class VcfMetaTest
{

	@Test
	public void add()
	{
		VcfMeta vcfMeta = new VcfMeta();
		vcfMeta.add("key", "val");
		assertEquals(vcfMeta.get("key"), "val");
	}

	@Test
	public void addAltMeta()
	{
		VcfMeta vcfMeta = new VcfMeta();
		VcfMetaAlt vcfMetaAlt1 = new VcfMetaAlt(Collections.singletonMap("ID", "val1"));
		VcfMetaAlt vcfMetaAlt2 = new VcfMetaAlt(Collections.singletonMap("ID", "val2"));
		vcfMeta.addAltMeta(vcfMetaAlt1);
		vcfMeta.addAltMeta(vcfMetaAlt1); // check for no duplicates
		vcfMeta.addAltMeta(vcfMetaAlt2);
		List<VcfMetaAlt> vcfMetaAlts = Lists.newArrayList(vcfMeta.getAltMeta());
		assertEquals(vcfMetaAlts, Arrays.asList(vcfMetaAlt1, vcfMetaAlt2));
	}

	@Test
	public void addContigMeta()
	{
		VcfMeta vcfMeta = new VcfMeta();
		VcfMetaContig vcfMetaContig1 = new VcfMetaContig(Collections.singletonMap("ID", "val1"));
		VcfMetaContig vcfMetaContig2 = new VcfMetaContig(Collections.singletonMap("ID", "val2"));
		vcfMeta.addContigMeta(vcfMetaContig1);
		vcfMeta.addContigMeta(vcfMetaContig1); // check for no duplicates
		vcfMeta.addContigMeta(vcfMetaContig2);
		List<VcfMetaContig> VcfMetaContigs = Lists.newArrayList(vcfMeta.getContigMeta());
		assertEquals(VcfMetaContigs, Arrays.asList(vcfMetaContig1, vcfMetaContig2));
	}

	@Test
	public void addFilterMeta()
	{
		VcfMeta vcfMeta = new VcfMeta();
		VcfMetaFilter vcfMetaFilter1 = new VcfMetaFilter(Collections.singletonMap("ID", "val1"));
		VcfMetaFilter vcfMetaFilter2 = new VcfMetaFilter(Collections.singletonMap("ID", "val2"));
		vcfMeta.addFilterMeta(vcfMetaFilter1);
		vcfMeta.addFilterMeta(vcfMetaFilter1); // check for no duplicates
		vcfMeta.addFilterMeta(vcfMetaFilter2);
		List<VcfMetaFilter> VcfMetaFilters = Lists.newArrayList(vcfMeta.getFilterMeta());
		assertEquals(VcfMetaFilters, Arrays.asList(vcfMetaFilter1, vcfMetaFilter2));
	}

	@Test
	public void addFormatMeta()
	{
		VcfMeta vcfMeta = new VcfMeta();
		VcfMetaFormat vcfMetaFormat1 = new VcfMetaFormat(Collections.singletonMap("ID", "val1"));
		VcfMetaFormat vcfMetaFormat2 = new VcfMetaFormat(Collections.singletonMap("ID", "val2"));
		vcfMeta.addFormatMeta(vcfMetaFormat1);
		vcfMeta.addFormatMeta(vcfMetaFormat1); // check for no duplicates
		vcfMeta.addFormatMeta(vcfMetaFormat2);
		List<VcfMetaFormat> VcfMetaFormats = Lists.newArrayList(vcfMeta.getFormatMeta());
		assertEquals(VcfMetaFormats, Arrays.asList(vcfMetaFormat1, vcfMetaFormat2));
	}

	@Test
	public void addInfoMeta()
	{
		VcfMeta vcfMeta = new VcfMeta();
		VcfMetaInfo vcfMetaInfo1 = new VcfMetaInfo(Collections.singletonMap("ID", "val1"));
		VcfMetaInfo vcfMetaInfo2 = new VcfMetaInfo(Collections.singletonMap("ID", "val2"));
		vcfMeta.addInfoMeta(vcfMetaInfo1);
		vcfMeta.addInfoMeta(vcfMetaInfo1); // check for no duplicates
		vcfMeta.addInfoMeta(vcfMetaInfo2);
		List<VcfMetaInfo> VcfMetaInfos = Lists.newArrayList(vcfMeta.getInfoMeta());
		assertEquals(VcfMetaInfos, Arrays.asList(vcfMetaInfo1, vcfMetaInfo2));
	}

	@Test
	public void addPedigreeMeta()
	{
		VcfMeta vcfMeta = new VcfMeta();
		VcfMetaPedigree vcfMetaPedigree1 = new VcfMetaPedigree(Collections.singletonMap("key", "val1"));
		VcfMetaPedigree vcfMetaPedigree2 = new VcfMetaPedigree(Collections.singletonMap("key", "val2"));
		vcfMeta.addPedigreeMeta(vcfMetaPedigree1);
		vcfMeta.addPedigreeMeta(vcfMetaPedigree2);
		List<VcfMetaPedigree> VcfMetaPedigrees = Lists.newArrayList(vcfMeta.getPedigreeMeta());
		assertEquals(VcfMetaPedigrees, Arrays.asList(vcfMetaPedigree1, vcfMetaPedigree2));
	}

	@Test
	public void addSampleMeta()
	{
		VcfMeta vcfMeta = new VcfMeta();
		VcfMetaSample vcfMetaSample1 = new VcfMetaSample(Collections.singletonMap("ID", "val1"));
		VcfMetaSample vcfMetaSample2 = new VcfMetaSample(Collections.singletonMap("ID", "val2"));
		vcfMeta.addSampleMeta(vcfMetaSample1);
		vcfMeta.addSampleMeta(vcfMetaSample1); // check for no duplicates
		vcfMeta.addSampleMeta(vcfMetaSample2);
		List<VcfMetaSample> VcfMetaSamples = Lists.newArrayList(vcfMeta.getSampleMeta());
		assertEquals(VcfMetaSamples, Arrays.asList(vcfMetaSample1, vcfMetaSample2));
	}

	@Test
	public void get()
	{
		VcfMeta vcfMeta = new VcfMeta();
		vcfMeta.add("key", "val");
		assertEquals(vcfMeta.get("key"), "val");
		assertEquals(vcfMeta.get("unknown-key"), null);
	}

	@Test
	public void getAltMetaString()
	{
		VcfMeta vcfMeta = new VcfMeta();
		VcfMetaAlt vcfMetaAlt1 = new VcfMetaAlt(Collections.singletonMap("ID", "val1"));
		vcfMeta.addAltMeta(vcfMetaAlt1);
		assertEquals(vcfMeta.getAltMeta("val1"), vcfMetaAlt1);
		assertNull(vcfMeta.getAltMeta("unknown-key"));
	}

	@Test
	public void getColNames()
	{
		String[] colNames = new String[]{"col1", "col2"};
		VcfMeta vcfMeta = new VcfMeta();
		vcfMeta.setColNames(colNames);
		assertEquals(vcfMeta.getColNames(), colNames);
	}

	@Test
	public void getContigMetaString()
	{
		VcfMeta vcfMeta = new VcfMeta();
		VcfMetaContig vcfMetaContig1 = new VcfMetaContig(Collections.singletonMap("ID", "val1"));
		vcfMeta.addContigMeta(vcfMetaContig1);
		assertEquals(vcfMeta.getContigMeta("val1"), vcfMetaContig1);
		assertNull(vcfMeta.getContigMeta("unknown-key"));
	}

	@Test
	public void getFileFormat()
	{
		VcfMeta vcfMeta = new VcfMeta();
		vcfMeta.add("fileformat", "format");
		assertEquals(vcfMeta.getFileFormat(), "format");
	}

	@Test
	public void getFilterMetaString()
	{
		VcfMeta vcfMeta = new VcfMeta();
		VcfMetaFilter vcfMetaFilter1 = new VcfMetaFilter(Collections.singletonMap("ID", "val1"));
		vcfMeta.addFilterMeta(vcfMetaFilter1);
		assertEquals(vcfMeta.getFilterMeta("val1"), vcfMetaFilter1);
		assertNull(vcfMeta.getFilterMeta("unknown-key"));
	}

	@Test
	public void getFormatMetaString()
	{
		VcfMeta vcfMeta = new VcfMeta();
		VcfMetaFormat vcfMetaFormat1 = new VcfMetaFormat(Collections.singletonMap("ID", "val1"));
		vcfMeta.addFormatMeta(vcfMetaFormat1);
		assertEquals(vcfMeta.getFormatMeta("val1"), vcfMetaFormat1);
		assertNull(vcfMeta.getFormatMeta("unknown-key"));
	}

	@Test
	public void getInfoMetaString()
	{
		VcfMeta vcfMeta = new VcfMeta();
		VcfMetaInfo vcfMetaInfo1 = new VcfMetaInfo(Collections.singletonMap("ID", "val1"));
		vcfMeta.addInfoMeta(vcfMetaInfo1);
		assertEquals(vcfMeta.getInfoMeta("val1"), vcfMetaInfo1);
		assertNull(vcfMeta.getInfoMeta("unknown-key"));
	}

	@Test
	public void getSampleMetaString()
	{
		VcfMeta vcfMeta = new VcfMeta();
		VcfMetaSample vcfMetaSample1 = new VcfMetaSample(Collections.singletonMap("ID", "val1"));
		vcfMeta.addSampleMeta(vcfMetaSample1);
		assertEquals(vcfMeta.getSampleMeta("val1"), vcfMetaSample1);
		assertNull(vcfMeta.getSampleMeta("unknown-key"));
	}

	@Test
	public void getSampleName()
	{
		String[] colNames = new String[]{"x", "x", "x", "x", "x", "x", "x", "x", "x", "s1", "s2"};
		VcfMeta vcfMeta = new VcfMeta();
		vcfMeta.setColNames(colNames);
		assertEquals(vcfMeta.getSampleName(0), "s1");
		assertEquals(vcfMeta.getSampleName(1), "s2");
	}

	@Test
	public void getSampleNames()
	{
		String[] colNames = new String[]{"x", "x", "x", "x", "x", "x", "x", "x", "x", "s1", "s2"};
		VcfMeta vcfMeta = new VcfMeta();
		vcfMeta.setColNames(colNames);
		assertEquals(Lists.newArrayList(vcfMeta.getSampleNames()), Arrays.asList("s1", "s2"));
	}
}
