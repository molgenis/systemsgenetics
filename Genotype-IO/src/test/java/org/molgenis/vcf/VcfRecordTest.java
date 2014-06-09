package org.molgenis.vcf;

import static org.mockito.Mockito.mock;
import static org.testng.Assert.assertEquals;
import static org.testng.Assert.assertFalse;
import static org.testng.Assert.assertTrue;

import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import org.molgenis.vcf.meta.VcfMeta;
import org.testng.annotations.Test;

import com.google.common.collect.Lists;
import org.molgenis.genotype.Allele;

public class VcfRecordTest
{

	@Test(expectedExceptions = IllegalArgumentException.class)
	public void VcfRecordVcfMeta()
	{
		new VcfRecord(null);
	}

	@Test(expectedExceptions = IllegalArgumentException.class)
	public void VcfRecordVcfMetaString()
	{
		new VcfRecord(null, null);
	}

	@Test
	public void createClone()
	{
		VcfMeta vcfMeta = mock(VcfMeta.class);
		String[] tokens = new String[0];
		VcfRecord vcfRecord = new VcfRecord(vcfMeta, tokens);
		assertEquals(vcfRecord.createClone(), vcfRecord);
	}

	@Test
	public void getAlternateAlleles()
	{
		VcfMeta vcfMeta = mock(VcfMeta.class);
		String[] tokens = new String[]{"x", "x", "x", "x", "A,C,TT"};
		List<Allele> alternateAlleles = new VcfRecord(vcfMeta, tokens).getAlternateAlleles();
		assertEquals(alternateAlleles, Arrays.asList(Allele.A, Allele.C, Allele.create("TT")));
	}

	@Test
	public void getAlternateAlleles_missingValue()
	{
		VcfMeta vcfMeta = mock(VcfMeta.class);
		String[] tokens = new String[]{"x", "x", "x", "x", "."};
		List<Allele> alternateAlleles = new VcfRecord(vcfMeta, tokens).getAlternateAlleles();
		assertEquals(alternateAlleles, Collections.emptyList());
	}
	
	@Test
	public void getChromosome()
	{
		VcfMeta vcfMeta = mock(VcfMeta.class);
		String[] tokens = new String[]{"1", "x", "x", "x"};
		String chromosome = new VcfRecord(vcfMeta, tokens).getChromosome();
		assertEquals(chromosome, "1");
	}

	@Test
	public void getFilterStatus()
	{
		VcfMeta vcfMeta = mock(VcfMeta.class);
		String[] tokens = new String[]{"x", "x", "x", "x", "x", "x", "status"};
		String filterStatus = new VcfRecord(vcfMeta, tokens).getFilterStatus();
		assertEquals(filterStatus, "status");
	}

	@Test
	public void getFilterStatus_missingValue()
	{
		VcfMeta vcfMeta = mock(VcfMeta.class);
		String[] tokens = new String[]{"x", "x", "x", "x", "x", "x", "."};
		String filterStatus = new VcfRecord(vcfMeta, tokens).getFilterStatus();
		assertEquals(filterStatus, null);
	}
	
	@Test
	public void getFormat()
	{
		VcfMeta vcfMeta = mock(VcfMeta.class);
		String[] tokens = new String[]{"x", "x", "x", "x", "x", "x", "x", "x", "F1:F2:F3"};
		String[] format = new VcfRecord(vcfMeta, tokens).getFormat();
		assertEquals(format, new String[]{"F1", "F2", "F3"});
	}

	@Test
	public void getFormatIndex()
	{
		VcfMeta vcfMeta = mock(VcfMeta.class);
		String[] tokens = new String[]{"x", "x", "x", "x", "x", "x", "x", "x", "F1:F2:F3"};
		assertEquals(0, new VcfRecord(vcfMeta, tokens).getFormatIndex("F1"));
		assertEquals(1, new VcfRecord(vcfMeta, tokens).getFormatIndex("F2"));
		assertEquals(2, new VcfRecord(vcfMeta, tokens).getFormatIndex("F3"));
		assertEquals(-1, new VcfRecord(vcfMeta, tokens).getFormatIndex("bad"));
	}

	@Test
	public void getIdentifiers()
	{
		VcfMeta vcfMeta = mock(VcfMeta.class);
		String[] tokens = new String[]{"x", "x", "id1;id2;id3"};
		List<String> identifiers = new VcfRecord(vcfMeta, tokens).getIdentifiers();
		assertEquals(identifiers, Arrays.asList("id1", "id2", "id3"));
	}

	@Test
	public void getIdentifiers_missingValue()
	{
		VcfMeta vcfMeta = mock(VcfMeta.class);
		String[] tokens = new String[]{"x", "x", "."};
		List<String> identifiers = new VcfRecord(vcfMeta, tokens).getIdentifiers();
		assertEquals(identifiers, Collections.emptyList());
	}
	
	@Test
	public void getInformation()
	{
		VcfMeta vcfMeta = mock(VcfMeta.class);
		String[] tokens = new String[]{"x", "x", "x", "x", "x", "x", "x", "K1=V1;K2=V2;K3", "x"};
		Iterator<VcfInfo> information = new VcfRecord(vcfMeta, tokens).getInformation().iterator();
		assertTrue(information.hasNext());
		assertEquals(information.next(), new VcfInfo(vcfMeta, "K1", "V1"));
		assertTrue(information.hasNext());
		assertEquals(information.next(), new VcfInfo(vcfMeta, "K2", "V2"));
		assertTrue(information.hasNext());
		assertEquals(information.next(), new VcfInfo(vcfMeta, "K3", null));
		assertFalse(information.hasNext());
	}

	@Test
	public void getNrSamples()
	{
		VcfMeta vcfMeta = mock(VcfMeta.class);
		String[] tokens = new String[]{"x", "x", "x", "x", "x", "x", "x", "x", "x", "s1", "s2", "s3"};
		assertEquals(new VcfRecord(vcfMeta, tokens).getNrSamples(), 3);
	}
	
	@Test
	public void getNrSamples_oneSample()
	{
		VcfMeta vcfMeta = mock(VcfMeta.class);
		String[] tokens = new String[]{"x", "x", "x", "x", "x", "x", "x", "x", "x", "s1"};
		assertEquals(new VcfRecord(vcfMeta, tokens).getNrSamples(), 1);
	}
	
	@Test
	public void getNrSamples_noSamples()
	{
		VcfMeta vcfMeta = mock(VcfMeta.class);
		String[] tokens = new String[]{"x", "x", "x"};
		assertEquals(new VcfRecord(vcfMeta, tokens).getNrSamples(), 0);
	}

	@Test
	public void getPosition()
	{
		VcfMeta vcfMeta = mock(VcfMeta.class);
		String[] tokens = new String[]{"x", "123"};
		assertEquals(new VcfRecord(vcfMeta, tokens).getPosition(), 123);
	}

	@Test
	public void getQuality()
	{
		VcfMeta vcfMeta = mock(VcfMeta.class);
		String[] tokens = new String[]{"x", "x", "x", "x", "x", "qual"};
		assertEquals(new VcfRecord(vcfMeta, tokens).getQuality(), "qual");
	}

	@Test
	public void getQuality_missingValue()
	{
		VcfMeta vcfMeta = mock(VcfMeta.class);
		String[] tokens = new String[]{"x", "x", "x", "x", "x", "."};
		assertEquals(new VcfRecord(vcfMeta, tokens).getQuality(), null);
	}
	
	@Test
	public void getReferenceAllele()
	{
		VcfMeta vcfMeta = mock(VcfMeta.class);
		String[] tokens = new String[]{"x", "x", "x", "A"};
		Allele referenceAllele = new VcfRecord(vcfMeta, tokens).getReferenceAllele();
		assertEquals(referenceAllele, Allele.A);		
	}

	@Test
	public void getSamples()
	{
		VcfMeta vcfMeta = mock(VcfMeta.class);
		String[] tokens = new String[]{"x", "x", "x", "x", "x", "x", "x", "x", "x", "s1", "s2", "s3"};
		Iterable<VcfSample> infoMap = new VcfRecord(vcfMeta, tokens).getSamples();
		assertEquals(Lists.newArrayList(infoMap).size(), 3);
	}

	@Test
	public void getSamples_noSamples()
	{
		VcfMeta vcfMeta = mock(VcfMeta.class);
		String[] tokens = new String[]{"x", "x", "x", "x", "x", "x", "x"};
		Iterable<VcfSample> infoMap = new VcfRecord(vcfMeta, tokens).getSamples();
		assertFalse(infoMap.iterator().hasNext());
	}
	
	@Test
	public void reset()
	{
		VcfMeta vcfMeta = mock(VcfMeta.class);
		String[] tokens = new String[]{"1", "x"};
		VcfRecord vcfRecord = new VcfRecord(vcfMeta, tokens);
		assertEquals(vcfRecord.getChromosome(), "1");
		String[] resetTokens = new String[]{"2", "x"};
		vcfRecord.reset(resetTokens);
		assertEquals(vcfRecord.getChromosome(), "2");
	}
}
