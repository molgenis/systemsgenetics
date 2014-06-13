package org.molgenis.vcf;

import static org.mockito.Mockito.mock;
import static org.mockito.Mockito.when;
import static org.testng.Assert.assertEquals;
import static org.testng.Assert.assertNull;

import java.util.Arrays;
import java.util.List;
import org.molgenis.genotype.Allele;

import org.testng.annotations.Test;

public class VcfSampleTest
{
	@Test(expectedExceptions = IllegalArgumentException.class)
	public void VcfSampleVcfMetaVcfRecord()
	{
		new VcfSample(null, null);
	}

	@Test
	public void getAlleles_noDataTypes()
	{
		VcfRecord vcfRecord = mock(VcfRecord.class);
		when(vcfRecord.getFormat()).thenReturn(new String[0]);
		assertNull(new VcfSample(vcfRecord, new String[]{}).getAlleles());
	}
	
	@Test
	public void getAlleles_noGt()
	{
		VcfRecord vcfRecord = mock(VcfRecord.class);
		when(vcfRecord.getFormat()).thenReturn(new String[]{"NOT-GT"});
		assertNull(new VcfSample(vcfRecord, new String[]{}).getAlleles());
	}
	
	@Test
	public void getAlleles_phased()
	{
		VcfRecord vcfRecord = mock(VcfRecord.class);
		when(vcfRecord.getReferenceAllele()).thenReturn(Allele.A);
		when(vcfRecord.getAlternateAlleles()).thenReturn(Arrays.asList(Allele.C));
		when(vcfRecord.getFormat()).thenReturn(new String[]{"GT"});
		String[] tokens = new String[]{"1|0"};
		List<Allele> alleles = new VcfSample(vcfRecord, tokens).getAlleles();
		assertEquals(Arrays.asList(Allele.C, Allele.A), alleles);
	}

	@Test
	public void getAlleles_unphased()
	{
		VcfRecord vcfRecord = mock(VcfRecord.class);
		when(vcfRecord.getReferenceAllele()).thenReturn(Allele.A);
		when(vcfRecord.getAlternateAlleles()).thenReturn(Arrays.asList(Allele.C));
		when(vcfRecord.getFormat()).thenReturn(new String[]{"GT"});
		String[] tokens = new String[]{"1/0"};
		List<Allele> alleles = new VcfSample(vcfRecord, tokens).getAlleles();
		assertEquals(Arrays.asList(Allele.C, Allele.A), alleles);
	}
	
	@Test
	public void getAlleles_triploid()
	{
		VcfRecord vcfRecord = mock(VcfRecord.class);
		when(vcfRecord.getReferenceAllele()).thenReturn(Allele.A);
		when(vcfRecord.getAlternateAlleles()).thenReturn(Arrays.asList(Allele.C, Allele.T, Allele.G));
		when(vcfRecord.getFormat()).thenReturn(new String[]{"GT"});
		String[] tokens = new String[]{"2/0/1"};
		List<Allele> alleles = new VcfSample(vcfRecord, tokens).getAlleles();
		assertEquals(Arrays.asList(Allele.T, Allele.A, Allele.C), alleles);
	}
	
	@Test
	public void getAlleles_multipleChars()
	{
		VcfRecord vcfRecord = mock(VcfRecord.class);
		when(vcfRecord.getReferenceAllele()).thenReturn(Allele.A);
		when(vcfRecord.getAlternateAlleles()).thenReturn(Arrays.asList(Allele.C, Allele.C, Allele.C, Allele.C, Allele.T, Allele.C, Allele.C, Allele.C, Allele.C, Allele.G));
		when(vcfRecord.getFormat()).thenReturn(new String[]{"GT"});
		String[] tokens = new String[]{"10/./5"};
		List<Allele> alleles = new VcfSample(vcfRecord, tokens).getAlleles();
		assertEquals(Arrays.asList(Allele.G, Allele.ZERO, Allele.T), alleles);
	}
	
	@Test
	public void getAlleles_noCall()
	{
		VcfRecord vcfRecord = mock(VcfRecord.class);
		when(vcfRecord.getReferenceAllele()).thenReturn(Allele.A);
		when(vcfRecord.getAlternateAlleles()).thenReturn(Arrays.asList(Allele.C));
		when(vcfRecord.getFormat()).thenReturn(new String[]{"GT"});
		String[] tokens = new String[]{"./0"};
		List<Allele> alleles = new VcfSample(vcfRecord, tokens).getAlleles();
		assertEquals(Arrays.asList(Allele.ZERO, Allele.A), alleles);
	}
	
	@Test
	public void getAlleles_noCall2()
	{
		VcfRecord vcfRecord = mock(VcfRecord.class);
		when(vcfRecord.getReferenceAllele()).thenReturn(Allele.A);
		when(vcfRecord.getAlternateAlleles()).thenReturn(Arrays.asList(Allele.A));
		when(vcfRecord.getFormat()).thenReturn(new String[]{"GT"});
		String[] tokens = new String[]{"0/."};
		List<Allele> alleles = new VcfSample(vcfRecord, tokens).getAlleles();
		assertEquals(Arrays.asList(Allele.A, Allele.ZERO), alleles);
	}

	@Test
	public void getData()
	{
		VcfRecord vcfRecord = mock(VcfRecord.class);
		String[] tokens = new String[]{"something", "0/."};
		String data = new VcfSample(vcfRecord, tokens).getData(1);
		assertEquals(data, "0/.");
	}

	@Test
	public void getData_missingValue()
	{
		VcfRecord vcfRecord = mock(VcfRecord.class);
		String[] tokens = new String[]{".", "0/."};
		String data = new VcfSample(vcfRecord, tokens).getData(0);
		assertEquals(data, null);
	}
	
	@Test
	public void getPhasings()
	{
		VcfRecord vcfRecord = mock(VcfRecord.class);
		when(vcfRecord.getReferenceAllele()).thenReturn(Allele.A);
		when(vcfRecord.getAlternateAlleles()).thenReturn(Arrays.asList(Allele.C, Allele.T, Allele.G));
		when(vcfRecord.getFormat()).thenReturn(new String[]{"GT"});
		String[] tokens = new String[]{"2/0|1"};
		List<Boolean> phasings = new VcfSample(vcfRecord, tokens).getPhasings();
		assertEquals(Arrays.asList(Boolean.FALSE, Boolean.TRUE), phasings);
	}

	@Test
	public void getPhasings_noGt()
	{
		VcfRecord vcfRecord = mock(VcfRecord.class);
		when(vcfRecord.getFormat()).thenReturn(new String[]{"NOT-GT"});
		assertNull(new VcfSample(vcfRecord, new String[]{}).getPhasings());
	}
	
	@Test
	public void getPhasings_noDataTypes()
	{
		VcfRecord vcfRecord = mock(VcfRecord.class);
		when(vcfRecord.getFormat()).thenReturn(new String[0]);
		assertNull(new VcfSample(vcfRecord, new String[]{}).getPhasings());
	}
	
	@Test
	public void reset()
	{
		VcfRecord vcfRecord = mock(VcfRecord.class);
		when(vcfRecord.getReferenceAllele()).thenReturn(Allele.A);
		when(vcfRecord.getAlternateAlleles()).thenReturn(Arrays.asList(Allele.C));
		when(vcfRecord.getFormat()).thenReturn(new String[]{"GT"});
		VcfSample vcfSample = new VcfSample(vcfRecord, new String[]{"1|0"});
		assertEquals(Arrays.asList(Allele.C, Allele.A), vcfSample.getAlleles());
		
		vcfSample.reset(new String[]{"0/1"});
		assertEquals(Arrays.asList(Allele.A, Allele.C), vcfSample.getAlleles());
	}
}
