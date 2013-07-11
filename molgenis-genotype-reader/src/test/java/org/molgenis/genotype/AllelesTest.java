package org.molgenis.genotype;

import static org.testng.Assert.assertEquals;
import static org.testng.Assert.assertNotNull;

import java.util.Arrays;
import java.util.List;

import org.testng.annotations.Test;

public class AllelesTest
{
	@Test
	public void getAlleles()
	{
		Alleles alleles = Alleles.createBasedOnString(Arrays.asList("A", "T"));
		assertNotNull(alleles.getAlleles());
		assertEquals(alleles.getAlleles().size(), 2);
		assertEquals(alleles.getAlleles().get(0).getAlleleAsString(), "A");
		assertEquals(alleles.getAllelesAsString().get(1), "T");
		assertEquals(alleles.getAllelesAsChars()[0], 'A');
		assertEquals(alleles.getAllelesAsChars()[1], 'T');

		Alleles alleles2 = Alleles.createBasedOnChars('A', 'T');
		assertEquals(alleles2, alleles);
	}

	@Test
	public void getAsChar()
	{
		Alleles alleles = Alleles.createBasedOnString(Arrays.asList("A", "T"));
		assertNotNull(alleles.getAllelesAsChars());
		assertEquals(alleles.getAllelesAsChars().length, 2);
		assertEquals(alleles.getAllelesAsChars()[0], 'A');
		assertEquals(alleles.getAllelesAsChars()[1], 'T');

		assertEquals(alleles == Alleles.createBasedOnChars('A', 'T'), true);

	}

	@Test(expectedExceptions = RuntimeException.class)
	public void getNonSnpAsChar()
	{
		Alleles.createBasedOnString(Arrays.asList("A", "T", "CG")).getAllelesAsChars();
	}

	@Test
	public void createWithChars()
	{
		Alleles variantAlleles = Alleles.createBasedOnChars('C', 'A');
		assertEquals(variantAlleles.getAllelesAsChars(), new char[]
		{ 'C', 'A' });
	}

	@Test
	public void swap()
	{
		List<String> alleles = Arrays.asList("A", "T", "C", "G");
		Alleles variantAlleles = Alleles.createBasedOnString(alleles);
		Alleles swapped = variantAlleles.getComplement();
		assertEquals(swapped.getAllelesAsString(), Arrays.asList("T", "A", "G", "C"));
		assertEquals(swapped.getAllelesAsChars(), new char[]
		{ 'T', 'A', 'G', 'C' });
	}

	@Test
	public void swapSnp()
	{
		Alleles variantAlleles = Alleles.createBasedOnChars('A', 'G');
		Alleles swapped = variantAlleles.getComplement();
		assertEquals(swapped.getAllelesAsString(), Arrays.asList("T", "C"));
		assertEquals(swapped.getAllelesAsChars(), new char[]
		{ 'T', 'C' });
	}

	@Test
	public void sameAlleles()
	{
		Alleles variantAlleles = Alleles.createBasedOnChars('A', 'G');
		Alleles variantAlleles2 = Alleles.createBasedOnChars('T', 'G');
		Alleles variantAlleles3 = Alleles.createBasedOnChars('G', 'A');
		Alleles variantAlleles4 = Alleles.createBasedOnChars(new char[]
		{ 'A', 'G', 'T' });
		Alleles variantAlleles5 = Alleles.createBasedOnChars(new char[]
		{ 'A' });
		Alleles variantAlleles6 = Alleles.createBasedOnChars('A', 'G');

		assertEquals(variantAlleles.sameAlleles(variantAlleles2), false);
		assertEquals(variantAlleles.sameAlleles(variantAlleles3), true);
		assertEquals(variantAlleles.sameAlleles(variantAlleles4), false);
		assertEquals(variantAlleles.sameAlleles(variantAlleles5), false);
		assertEquals(variantAlleles.sameAlleles(variantAlleles6), true);
	}

	@Test
	public void isAtOrGcSnp()
	{

		Alleles variantAlleles;

		variantAlleles = Alleles.createBasedOnChars('A', 'G');
		assertEquals(variantAlleles.isAtOrGcSnp(), false);

		variantAlleles = Alleles.createBasedOnString(Arrays.asList("A", "G"));
		assertEquals(variantAlleles.isAtOrGcSnp(), false);

		variantAlleles = Alleles.createBasedOnString(Arrays.asList("G", "C"));
		assertEquals(variantAlleles.isAtOrGcSnp(), true);

		variantAlleles = Alleles.createBasedOnChars('A', 'T');
		assertEquals(variantAlleles.isAtOrGcSnp(), true);

		variantAlleles = Alleles.createBasedOnString(Arrays.asList("G", "C", "GC"));
		assertEquals(variantAlleles.isAtOrGcSnp(), false);

		variantAlleles = Alleles.createBasedOnString(Arrays.asList("G", "C", "G"));
		assertEquals(variantAlleles.isAtOrGcSnp(), true);

		variantAlleles = Alleles.createBasedOnString(Arrays.asList("G", "C", "T"));
		assertEquals(variantAlleles.isAtOrGcSnp(), false);
	}

	@Test
	public void testPooling()
	{

		Alleles alleles1 = Alleles.createBasedOnChars('A', 'T');
		Alleles alleles2 = Alleles.createBasedOnString("A", "T");
		Alleles alleles3 = Alleles.createBasedOnChars('T', 'A');
		Alleles alleles4 = Alleles.createBasedOnString("A", "TTTT");
		Alleles alleles5 = Alleles.createBasedOnString("TTTT", "C");
		Alleles alleles6 = Alleles.createBasedOnString("TTT", "C");

		assertEquals(alleles1 == alleles2, true);
		assertEquals(alleles1.getAllelesAsString() == alleles2.getAllelesAsString(), true);
		assertEquals(alleles1.getAllelesAsString() == alleles3.getAllelesAsString(), false);

		assertEquals(alleles1.getAllelesAsChars() == alleles2.getAllelesAsChars(), true);
		assertEquals(alleles1.getAllelesAsChars() == alleles3.getAllelesAsChars(), false);

		assertEquals(alleles1.getAllelesAsString() == alleles3.getAllelesAsString(), false);
		assertEquals(alleles1.getAllelesAsString() == alleles4.getAllelesAsString(), false);

		assertEquals(alleles1.getAlleles() == alleles2.getAlleles(), true);
		assertEquals(alleles2.getAlleles() == alleles3.getAlleles(), false);

		assertEquals(alleles1.getAlleles().get(0) == alleles2.getAlleles().get(0), true);
		assertEquals(alleles1.getAlleles().get(1) == alleles2.getAlleles().get(1), true);
		assertEquals(alleles1.getAlleles().get(1) == alleles3.getAlleles().get(0), true);

		assertEquals(alleles4.getAlleles().get(1) == alleles5.getAlleles().get(0), true);
		assertEquals(alleles6.getAlleles().get(0) == alleles5.getAlleles().get(0), false);

		assertEquals(alleles4.getAllelesAsString().get(1) == alleles5.getAllelesAsString().get(0), true);
		assertEquals(alleles6.getAllelesAsString().get(0) == alleles5.getAllelesAsString().get(0), false);

	}

	@Test
	public void contains()
	{
		Alleles alleles = Alleles.createBasedOnString(Arrays.asList("A", "T"));
		assertEquals(alleles.contains(Allele.A), true);
		assertEquals(alleles.contains(Allele.T), true);
		assertEquals(alleles.contains(Allele.C), false);
		assertEquals(alleles.contains(Allele.create("AA")), false);
	}
}
