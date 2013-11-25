package org.molgenis.genotype.util;

import static org.testng.Assert.assertEquals;

import java.util.Arrays;
import java.util.Comparator;

import org.testng.annotations.Test;

/**
 * This comparator can be used to sort chromosome names. In this order: 1-9
 * 10-19 etc ending with X Y MT
 * 
 * @author Patrick Deelen
 * 
 */
public class ChromosomeComparatorTest
{

	@Test
	public void compare()
	{

		Comparator<String> testComparator = new ChromosomeComparator();

		assertEquals(testComparator.compare("1", "2") < 0, true);
		assertEquals(testComparator.compare("1", "1") == 0, true);
		assertEquals(testComparator.compare("3", "1") > 0, true);

		assertEquals(testComparator.compare("3", "13") < 0, true);
		assertEquals(testComparator.compare("13", "3") > 0, true);

		assertEquals(testComparator.compare("9", "10") < 0, true);
		assertEquals(testComparator.compare("10", "11") < 0, true);

		assertEquals(testComparator.compare("9", "21") < 0, true);
		assertEquals(testComparator.compare("19", "20") < 0, true);

		assertEquals(testComparator.compare("22", "X") < 0, true);
		assertEquals(testComparator.compare("10", "X") < 0, true);
		assertEquals(testComparator.compare("1", "X") < 0, true);

		assertEquals(testComparator.compare("X", "5") > 0, true);
		assertEquals(testComparator.compare("5", "X") < 0, true);

		assertEquals(testComparator.compare("X", "13") > 0, true);
		assertEquals(testComparator.compare("13", "X") < 0, true);

		assertEquals(testComparator.compare("X", "X") == 0, true);
		assertEquals(testComparator.compare("Y", "X") > 0, true);

		assertEquals(testComparator.compare("MT", "X") > 0, true);
		assertEquals(testComparator.compare("X", "MT") < 0, true);

		assertEquals(testComparator.compare("MT", "Y") > 0, true);
		assertEquals(testComparator.compare("Y", "MT") < 0, true);

		assertEquals(testComparator.compare("MT", "1") > 0, true);
		assertEquals(testComparator.compare("1", "MT") < 0, true);

		assertEquals(testComparator.compare("MT", "10") > 0, true);
		assertEquals(testComparator.compare("10", "MT") < 0, true);
		
		assertEquals(testComparator.compare("22", "23") < 0, true);
		assertEquals(testComparator.compare("23", "22") > 0, true);

	}

	@Test
	public void sortTest()
	{

		String[] chrs = new String[]
		{ "X", "1", "16", "2", "MT", "Y", "22", "3", "1" };

		String[] chrsSorted = new String[]
		{ "1", "1", "2", "3", "16", "22", "X", "Y", "MT" };

		Arrays.sort(chrs, new ChromosomeComparator());

		assertEquals(chrs, chrsSorted);

	}
}
