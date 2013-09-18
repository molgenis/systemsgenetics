/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.genomicboundaries;


import static org.testng.Assert.*;
import org.testng.annotations.Test;

/**
 *
 * @author Patrick Deelen
 */
public class GenomicBoundaryTest {
	
	public GenomicBoundaryTest() {
	}
	
	

	/**
	 * Test of getChromosome method, of class GenomicBoundary.
	 */
	@Test
	public void testGetChromosome() {
		//System.out.println("getChromosome");
		GenomicBoundary instance = new GenomicBoundary("1", 1, 2);
		String expResult = "1";
		String result = instance.getChromosome();
		assertEquals(result, expResult);
	}

	/**
	 * Test of getStart method, of class GenomicBoundary.
	 */
	@Test
	public void testGetStart() {
		//System.out.println("getStart");
		GenomicBoundary instance = new GenomicBoundary("1", 1, 2);
		Integer expResult = 1;
		Integer result = instance.getStart();
		assertEquals(result, expResult);
	}

	/**
	 * Test of getStop method, of class GenomicBoundary.
	 */
	@Test
	public void testGetStop() {
		//System.out.println("getStop");
		GenomicBoundary instance = new GenomicBoundary("1", 1, 2);
		int expResult = 2;
		int result = instance.getStop();
		assertEquals(result, expResult);
	}

	/**
	 * Test of isInBoundarie method, of class GenomicBoundary.
	 */
	@Test
	public void testIsInBoundarie_int() {
		//System.out.println("isInBoundarie");
		int position = 1;
		GenomicBoundary instance = new GenomicBoundary("1", 1, 2);
		boolean expResult = true;
		boolean result = instance.isInBoundarie(position);
		assertEquals(result, expResult);
		
		position = 2;
		expResult = true;
		result = instance.isInBoundarie(position);
		assertEquals(result, expResult);
		
		position = 3;
		expResult = false;
		result = instance.isInBoundarie(position);
		assertEquals(result, expResult);
		
		position = 0;
		expResult = false;
		result = instance.isInBoundarie(position);
		assertEquals(result, expResult);
	}

	/**
	 * Test of isInBoundarie method, of class GenomicBoundary.
	 */
	@Test
	public void testIsInBoundarie_int_int() {
		//System.out.println("isInBoundarie");
		int margin = 2;
		GenomicBoundary instance= new GenomicBoundary("1", 10, 12);
		
		int position = 7;
		boolean expResult = false;
		boolean result = instance.isInBoundarie(position, margin);
		assertEquals(result, expResult);
		
		position = 8;
		expResult = true;
		result = instance.isInBoundarie(position, margin);
		assertEquals(result, expResult);
		
		position = 9;
		expResult = true;
		result = instance.isInBoundarie(position, margin);
		assertEquals(result, expResult);

		position = 14;
		expResult = true;
		result = instance.isInBoundarie(position, margin);
		assertEquals(result, expResult);
		
		position = 15;
		expResult = false;
		result = instance.isInBoundarie(position, margin);
		assertEquals(result, expResult);

		
	}

	/**
	 * Test of isPartOfBoundary method, of class GenomicBoundary.
	 */
	@Test
	public void testIsPartOfBoundary() {
		//System.out.println("isPartOfBoundary");
		GenomicBoundary other = new GenomicBoundary("1", 1, 5);
		GenomicBoundary instance = new GenomicBoundary("1", 1, 5);
		boolean expResult = true;
		boolean result = instance.isPartOfBoundary(other);
		assertEquals(result, expResult);

		
		other = new GenomicBoundary("1", 1, 6);
		expResult = true;
		result = instance.isPartOfBoundary(other);
		assertEquals(result, expResult);
		
		other = new GenomicBoundary("1", 1, 4);
		expResult = false;
		result = instance.isPartOfBoundary(other);
		assertEquals(result, expResult);
		
		other = new GenomicBoundary("1", 2, 6);
		expResult = false;
		result = instance.isPartOfBoundary(other);
		assertEquals(result, expResult);
		
		other = new GenomicBoundary("1", 2, 7);
		expResult = false;
		result = instance.isPartOfBoundary(other);
		assertEquals(result, expResult);
		
	}
	
	/**
	 * Test of isOverlaping method, of class GenomicBoundary.
	 */
	@Test
	public void isOverlaping() {
		
		//System.out.println("isOverlaping");
		GenomicBoundary other = new GenomicBoundary("1", 3, 5);
		GenomicBoundary instance = new GenomicBoundary("1", 3, 5);
		boolean expResult = true;
		boolean result = instance.isOverlaping(other);
		assertEquals(result, expResult);

		
		other = new GenomicBoundary("1", 1, 6);
		expResult = true;
		result = instance.isOverlaping(other);
		assertEquals(result, expResult);
		
		other = new GenomicBoundary("1", 1, 4);
		expResult = true;
		result = instance.isOverlaping(other);
		assertEquals(result, expResult);
		
		other = new GenomicBoundary("1", 2, 6);
		expResult = true;
		result = instance.isOverlaping(other);
		assertEquals(result, expResult);
		
		other = new GenomicBoundary("1", 2, 7);
		expResult = true;
		result = instance.isOverlaping(other);
		assertEquals(result, expResult);
		
		other = new GenomicBoundary("1", 5, 7);
		expResult = true;
		result = instance.isOverlaping(other);
		assertEquals(result, expResult);
		
		other = new GenomicBoundary("1", 6, 7);
		expResult = false;
		result = instance.isOverlaping(other);
		assertEquals(result, expResult);
		
		other = new GenomicBoundary("1", 7, 7);
		expResult = false;
		result = instance.isOverlaping(other);
		assertEquals(result, expResult);
		
		other = new GenomicBoundary("1", 7, 8);
		expResult = false;
		result = instance.isOverlaping(other);
		assertEquals(result, expResult);
		
		other = new GenomicBoundary("1", 2, 3);
		expResult = true;
		result = instance.isOverlaping(other);
		assertEquals(result, expResult);
		
		other = new GenomicBoundary("1", 1, 2);
		expResult = false;
		result = instance.isOverlaping(other);
		assertEquals(result, expResult);
		
		other = new GenomicBoundary("1", 1, 1);
		expResult = false;
		result = instance.isOverlaping(other);
		assertEquals(result, expResult);
		
	}
	
}
