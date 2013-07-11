package org.molgenis.genotype.util;

import static org.testng.Assert.assertEquals;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Iterator;

import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.DummySampleVariantsProvider;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variant.ReadOnlyGeneticVariant;
import org.molgenis.genotype.variant.sampleProvider.SampleVariantsProvider;
import org.testng.annotations.Test;

public class LdCalculatorTest
{

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

		GeneticVariant testInstance = ReadOnlyGeneticVariant.createSnp("rs1", 1, "chr1", sampleAllelesProvider, 'A',
				'C');

		Ld ld = LdCalculator.calculateLd(testInstance, testInstance);

		assertEquals(ld.getR2(), 1, 0.1);
		assertEquals(ld.getDPrime(), 1, 0.1);

		ArrayList<Double> hapFreqExpect = new ArrayList<Double>(Arrays.asList(0.4d, 0.0d, 0.0d, 0.6d));
		ArrayList<String> haps = new ArrayList<String>(Arrays.asList("CC", "CA", "AC", "AA"));

		assertEquals(ld.getHaplotypesFreq().keySet(), haps);
		// AssertJUnit.assertEquals(ld.getHaplotypesFreq().values(),
		// hapFreqExpect);

		assertEqualsDoubleCollection(ld.getHaplotypesFreq().values(), hapFreqExpect, 0.000001);
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

		GeneticVariant testInstance = ReadOnlyGeneticVariant.createSnp("rs1", 1, "chr1", sampleAllelesProvider, 'A',
				'C');

		sampleAlleles = new ArrayList<Alleles>();
		sampleAlleles.add(Alleles.createBasedOnChars('T', 'C'));
		sampleAlleles.add(Alleles.createBasedOnChars('T', 'T'));
		sampleAlleles.add(Alleles.createBasedOnChars('C', 'C'));
		sampleAlleles.add(Alleles.createBasedOnChars('T', 'T'));
		sampleAlleles.add(Alleles.createBasedOnChars('C', 'C'));
		sampleAllelesProvider = new DummySampleVariantsProvider(sampleAlleles);

		GeneticVariant testInstance2 = ReadOnlyGeneticVariant.createSnp("rs2", 1, "chr1", sampleAllelesProvider, 'T',
				'C');

		Ld ld = LdCalculator.calculateLd(testInstance, testInstance2);

		assertEquals(ld.getR2(), 0.1, 0.1);
		assertEquals(ld.getDPrime(), 0.4, 0.1);

		ArrayList<Double> hapFreqExpect = new ArrayList<Double>(Arrays.asList(0.2801730369238241d,
				0.21982696307617589d, 0.1198269630761759d, 0.38017303692382415d));
		ArrayList<String> haps = new ArrayList<String>(Arrays.asList("CC", "CT", "AC", "AT"));

		assertEquals(ld.getHaplotypesFreq().keySet(), haps);
		assertEquals(ld.getHaplotypesFreq().values(), hapFreqExpect);

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
