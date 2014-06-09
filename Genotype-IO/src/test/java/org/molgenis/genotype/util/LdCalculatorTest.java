package org.molgenis.genotype.util;

import static org.mockito.Mockito.mock;
import static org.testng.Assert.assertEquals;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Iterator;

import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.DummySampleVariantsProvider;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.ResourceTest;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variant.GeneticVariantMeta;
import org.molgenis.genotype.variant.ReadOnlyGeneticVariant;
import org.molgenis.genotype.variant.sampleProvider.SampleVariantsProvider;
import org.molgenis.genotype.vcf.VcfGenotypeData;
import org.testng.annotations.Test;

public class LdCalculatorTest extends ResourceTest
{
	private GeneticVariantMeta variantMeta;
	
	@Test
	public void calculateLd() throws LdCalculatorException
	{

		ArrayList<Alleles> sampleAlleles = new ArrayList<Alleles>();
		sampleAlleles.add(Alleles.createBasedOnChars('A', 'C'));
		sampleAlleles.add(Alleles.createBasedOnChars('A', 'C'));
		sampleAlleles.add(Alleles.createBasedOnChars('C', 'C'));
		sampleAlleles.add(Alleles.createBasedOnChars('A', 'A'));
		sampleAlleles.add(Alleles.createBasedOnChars('A', 'A'));
		SampleVariantsProvider sampleAllelesProvider = new DummySampleVariantsProvider(sampleAlleles);

		variantMeta = mock(GeneticVariantMeta.class);
		GeneticVariant testInstance = ReadOnlyGeneticVariant.createSnp(variantMeta , "rs1", 1, "chr1", sampleAllelesProvider, 'A',
				'C');

		Ld ld = LdCalculator.calculateLd(testInstance, testInstance);

		assertEquals(ld.getR2(), 1, 0.1);
		assertEquals(ld.getDPrime(), 1, 0.1);

		ArrayList<Double> hapFreqExpect = new ArrayList<Double>(Arrays.asList(0.4d, 0.0d, 0.0d, 0.6d));
		ArrayList<String> haps = new ArrayList<String>(Arrays.asList("CC", "CA", "AC", "AA"));

		assertEquals(ld.getHaplotypesFreq().keySet(), haps);
		assertEqualsDoubleCollection(ld.getHaplotypesFreq().values(), hapFreqExpect, 0.0001);

	}

	@Test
	public void calculateLd2() throws LdCalculatorException
	{

		ArrayList<Alleles> sampleAlleles = new ArrayList<Alleles>();
		sampleAlleles.add(Alleles.createBasedOnChars('A', 'C'));
		sampleAlleles.add(Alleles.createBasedOnChars('A', 'C'));
		sampleAlleles.add(Alleles.createBasedOnChars('C', 'C'));
		sampleAlleles.add(Alleles.createBasedOnChars('A', 'A'));
		sampleAlleles.add(Alleles.createBasedOnChars('A', 'A'));
		SampleVariantsProvider sampleAllelesProvider = new DummySampleVariantsProvider(sampleAlleles);

		GeneticVariant testInstance = ReadOnlyGeneticVariant.createSnp(variantMeta, "rs1", 1, "chr1", sampleAllelesProvider, 'A',
				'C');

		sampleAlleles = new ArrayList<Alleles>();
		sampleAlleles.add(Alleles.createBasedOnChars('T', 'C'));
		sampleAlleles.add(Alleles.createBasedOnChars('T', 'T'));
		sampleAlleles.add(Alleles.createBasedOnChars('C', 'C'));
		sampleAlleles.add(Alleles.createBasedOnChars('T', 'T'));
		sampleAlleles.add(Alleles.createBasedOnChars('C', 'C'));
		sampleAllelesProvider = new DummySampleVariantsProvider(sampleAlleles);

		GeneticVariant testInstance2 = ReadOnlyGeneticVariant.createSnp(variantMeta, "rs2", 1, "chr1", sampleAllelesProvider, 'T',
				'C');

		Ld ld = LdCalculator.calculateLd(testInstance, testInstance2);

		assertEquals(ld.getR2(), 0.1, 0.1);
		assertEquals(ld.getDPrime(), 0.4, 0.1);

		ArrayList<Double> hapFreqExpect = new ArrayList<Double>(Arrays.asList(0.2801730369238241d, 0.1198269630761759d,
				0.21982696307617589d, 0.38017303692382415d));
		ArrayList<String> haps = new ArrayList<String>(Arrays.asList("CC", "CT", "AC", "AT"));
		assertEquals(ld.getHaplotypesFreq().keySet(), haps);
		assertEqualsDoubleCollection(ld.getHaplotypesFreq().values(), hapFreqExpect, 0.00001);

	}

	@Test
	public void calculateLd3() throws LdCalculatorException, Exception
	{

		RandomAccessGenotypeData testData = new VcfGenotypeData(getLdTestVcf(), getLdTestVcfTbi(), 0.8);
		try
		{
			GeneticVariant rs2073738 = testData.getSnpVariantByPos("22", 19170956);
			GeneticVariant rs2073738_B = testData.getSnpVariantByPos("22", 19170957);
			GeneticVariant rs1210711 = testData.getSnpVariantByPos("22", 19339635);
			GeneticVariant rs1210711_B = testData.getSnpVariantByPos("22", 19339636);
			GeneticVariant rs5748165 = testData.getSnpVariantByPos("22", 19368723);

			Ld ld;

			ld = rs1210711.calculateLd(rs2073738);
			ArrayList<Double> hapFreqExpect = new ArrayList<Double>(Arrays.asList(1.0736996276734832E-38d,
					0.9629258517034068d, 0.006012024048096192d, 0.031062124248496994d));
			ArrayList<String> haps = new ArrayList<String>(Arrays.asList("CT", "CC", "TT", "TC"));
			assertEquals(ld.getHaplotypesFreq().keySet(), haps);
			assertEqualsDoubleCollection(ld.getHaplotypesFreq().values(), hapFreqExpect, 0.00001);

			ld = rs1210711.calculateLd(rs2073738_B);
			hapFreqExpect = new ArrayList<Double>(Arrays.asList(0.9629258517034068d, 1.0736996276734832E-38d,
					0.031062124248496994d, 0.006012024048096192));
			haps = new ArrayList<String>(Arrays.asList("CC", "CT", "TC", "TT"));
			assertEquals(ld.getHaplotypesFreq().keySet(), haps);
			assertEqualsDoubleCollection(ld.getHaplotypesFreq().values(), hapFreqExpect, 0.00001);

			ld = rs1210711_B.calculateLd(rs2073738);
			hapFreqExpect = new ArrayList<Double>(Arrays.asList(0.006012024048096192d, 0.031062124248496994d,
					1.0736996276734832E-38d, 0.9629258517034068));
			haps = new ArrayList<String>(Arrays.asList("TT", "TC", "CT", "CC"));
			assertEquals(ld.getHaplotypesFreq().keySet(), haps);
			assertEqualsDoubleCollection(ld.getHaplotypesFreq().values(), hapFreqExpect, 0.00001);

			ld = rs1210711_B.calculateLd(rs2073738_B);
			hapFreqExpect = new ArrayList<Double>(Arrays.asList(0.031062124248496994d, 0.006012024048096192d,
					0.9629258517034068d, 1.0736996276734832E-38d));
			haps = new ArrayList<String>(Arrays.asList("TC", "TT", "CC", "CT"));
			assertEquals(ld.getHaplotypesFreq().keySet(), haps);
			assertEqualsDoubleCollection(ld.getHaplotypesFreq().values(), hapFreqExpect, 0.00001);

			ld = rs1210711.calculateLd(rs5748165);
			hapFreqExpect = new ArrayList<Double>(Arrays.asList(0.0d, 0.9629258517034068d, 0.03707414829659319d, 0d));
			haps = new ArrayList<String>(Arrays.asList("CC", "CT", "TC", "TT"));
			assertEquals(ld.getHaplotypesFreq().keySet(), haps);
			assertEqualsDoubleCollection(ld.getHaplotypesFreq().values(), hapFreqExpect, 0.00001);
		}
		finally
		{
			testData.close();
		}
	}

	public void assertEqualsDoubleCollection(Collection<Double> observed, Collection<Double> expected, double delta)
			throws AssertionError
	{
		if (observed.size() != expected.size())
		{
			throw new AssertionError("List size differs");
		}

		Iterator<Double> observedIterator = observed.iterator();
		Iterator<Double> expectedIterator = expected.iterator();

		while (observedIterator.hasNext())
		{
			assertEquals(observedIterator.next(), expectedIterator.next(), delta);
		}

	}
}
