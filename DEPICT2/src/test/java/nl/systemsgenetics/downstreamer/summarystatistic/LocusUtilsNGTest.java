/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.downstreamer.summarystatistic;

import nl.systemsgenetics.downstreamer.summarystatistic.LocusUtils;
import nl.systemsgenetics.downstreamer.summarystatistic.OverlappableGenomicRange;
import java.util.List;
import java.util.Map;
import java.util.Set;
import nl.systemsgenetics.downstreamer.gene.Gene;
import nl.systemsgenetics.downstreamer.summarystatistic.filters.PvalueFilterSmaller;
import org.molgenis.genotype.RandomAccessGenotypeData;
import static org.testng.Assert.*;
import org.testng.annotations.AfterClass;
import org.testng.annotations.AfterMethod;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

/**
 *
 * @author patri
 */
public class LocusUtilsNGTest {
	
	public LocusUtilsNGTest() {
	}

	@BeforeClass
	public static void setUpClass() throws Exception {
	}

	@AfterClass
	public static void tearDownClass() throws Exception {
	}

	@BeforeMethod
	public void setUpMethod() throws Exception {
	}

	@AfterMethod
	public void tearDownMethod() throws Exception {
	}


	/**
	 * Test of partialGenomicRangeOverlapWindow method, of class LocusUtils.
	 */
	@Test
	public void testPartialGenomicRangeOverlapWindow() {
		OverlappableGenomicRange a = new Gene("A", "chr1", 10, 20, "");
		OverlappableGenomicRange b = new Gene("B", "chr2", 10, 20, "");
		int window = 0;
		boolean expResult = false;
		boolean result = LocusUtils.partialGenomicRangeOverlapWindow(a, b, window);
		assertEquals(result, expResult);
		
		a = new Gene("A", "chr1", 10, 20, "");
		b = new Gene("B", "1", 10, 20, "");
		window = 0;
		expResult = true;
		result = LocusUtils.partialGenomicRangeOverlapWindow(a, b, window);
		assertEquals(result, expResult);
		
		a = new Gene("A", "chr1", 10, 20, "");
		b = new Gene("B", "chr1", 12, 14, "");
		window = 0;
		expResult = true;
		result = LocusUtils.partialGenomicRangeOverlapWindow(a, b, window);
		assertEquals(result, expResult);
	
		a = new Gene("A", "chr1", 20, 10, "");
		b = new Gene("B", "chr1", 12, 14, "");
		window = 0;
		expResult = true;
		result = LocusUtils.partialGenomicRangeOverlapWindow(a, b, window);
		assertEquals(result, expResult);

		a = new Gene("A", "chr1", 12, 22, "");
		b = new Gene("B", "chr1", 10, 20, "");
		window = 0;
		expResult = true;
		result = LocusUtils.partialGenomicRangeOverlapWindow(a, b, window);
		assertEquals(result, expResult);

		a = new Gene("A", "chr1", 12, 19, "");
		b = new Gene("B", "chr1", 10, 20, "");
		window = 0;
		expResult = true;
		result = LocusUtils.partialGenomicRangeOverlapWindow(a, b, window);
		assertEquals(result, expResult);
		
		
		
		
		a = new Gene("A", "chr1", 10, 20, "");
		b = new Gene("B", "1", 10, 20, "");
		window = 1000;
		expResult = true;
		result = LocusUtils.partialGenomicRangeOverlapWindow(a, b, window);
		assertEquals(result, expResult);
		
		a = new Gene("A", "chr1", 10, 20, "");
		b = new Gene("B", "chr1", 12, 14, "");
		window = 1000;
		expResult = true;
		result = LocusUtils.partialGenomicRangeOverlapWindow(a, b, window);
		assertEquals(result, expResult);
	
		a = new Gene("A", "chr1", 20, 10, "");
		b = new Gene("B", "chr1", 12, 14, "");
		window = 100;
		expResult = true;
		result = LocusUtils.partialGenomicRangeOverlapWindow(a, b, window);
		assertEquals(result, expResult);

		a = new Gene("A", "chr1", 12, 22, "");
		b = new Gene("B", "chr1", 10, 20, "");
		window = 1000;
		expResult = true;
		result = LocusUtils.partialGenomicRangeOverlapWindow(a, b, window);
		assertEquals(result, expResult);

		a = new Gene("A", "chr1", 12, 19, "");
		b = new Gene("B", "chr1", 10, 20, "");
		window = 100;
		expResult = true;
		result = LocusUtils.partialGenomicRangeOverlapWindow(a, b, window);
		assertEquals(result, expResult);

		a = new Gene("A", "chr1", 10, 20, "");
		b = new Gene("B", "chr1", 21, 30, "");
		window = 0;
		expResult = false;
		result = LocusUtils.partialGenomicRangeOverlapWindow(a, b, window);
		assertEquals(result, expResult);
		
		a = new Gene("A", "chr1", 10, 20, "");
		b = new Gene("B", "chr1", 21, 30, "");
		window = 1;
		expResult = true;
		result = LocusUtils.partialGenomicRangeOverlapWindow(a, b, window);
		assertEquals(result, expResult);
		
		a = new Gene("A", "chr1", 10, 20, "");
		b = new Gene("B", "chr1", 15, 16, "");
		window = 100;
		expResult = true;
		result = LocusUtils.partialGenomicRangeOverlapWindow(a, b, window);
		assertEquals(result, expResult);
		
		a = new Gene("A", "chr1", 13, 15, "");
		b = new Gene("B", "chr1", 10, 20, "");
		window = 1;
		expResult = true;
		result = LocusUtils.partialGenomicRangeOverlapWindow(a, b, window);
		assertEquals(result, expResult);
		
		a = new Gene("A", "chr1", 13, 15, "");
		b = new Gene("B", "chr1", 10, 20, "");
		window = 100;
		expResult = true;
		result = LocusUtils.partialGenomicRangeOverlapWindow(a, b, window);
		assertEquals(result, expResult);
		
		a = new Gene("A", "chr1", 10, 15, "");
		b = new Gene("B", "chr1", 21, 22, "");
		window = 5;
		expResult = false;
		result = LocusUtils.partialGenomicRangeOverlapWindow(a, b, window);
		assertEquals(result, expResult);
		
		a = new Gene("A", "chr1", 10, 15, "");
		b = new Gene("B", "chr1", 21, 20, "");
		window = 5;
		expResult = true;
		result = LocusUtils.partialGenomicRangeOverlapWindow(a, b, window);
		assertEquals(result, expResult);
		
		a = new Gene("A", "chr1", 15, 10, "");
		b = new Gene("B", "chr1", 21, 20, "");
		window = 5;
		expResult = true;
		result = LocusUtils.partialGenomicRangeOverlapWindow(a, b, window);
		assertEquals(result, expResult);
		
		a = new Gene("A", "chr1", 15, 10, "");
		b = new Gene("B", "chr1", 0, 2, "");
		window = 5;
		expResult = false;
		result = LocusUtils.partialGenomicRangeOverlapWindow(a, b, window);
		assertEquals(result, expResult);
		
		
		a = new Gene("A", "chr1", 15, 10, "");
		b = new Gene("B", "chr1", 0, 2, "");
		window = 100;
		expResult = true;
		result = LocusUtils.partialGenomicRangeOverlapWindow(a, b, window);
		assertEquals(result, expResult);
		
	}
	
	
}
