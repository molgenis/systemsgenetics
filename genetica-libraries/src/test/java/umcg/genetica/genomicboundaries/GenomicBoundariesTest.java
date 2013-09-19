/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.genomicboundaries;


import java.util.Iterator;
import static org.testng.Assert.*;
import org.testng.annotations.Test;

/**
 *
 * @author Patrick Deelen
 */
public class GenomicBoundariesTest {
	
	public GenomicBoundariesTest() {
	}
	
	

	/**
	 * Test of addBoundary method, of class GenomicBoundaries.
	 */
	@Test
	public void testAddBoundary() {
		//System.out.println("addBoundary");
		String chromosome = "1";
		Integer beginPoint = 1;
		int endPoint = 3;
		GenomicBoundaries instance = new GenomicBoundaries();
		instance.addBoundary(chromosome, beginPoint, endPoint);
		boolean expResult = true;
		boolean result = instance.isInBoundary(chromosome, beginPoint, 0);
		assertEquals(result, expResult);
		
	}
	
	
	/**
	 * Test of removeSubBoundaries method, of class GenomicBoundaries.
	 */
	@Test
	public void testRemoveSubBoundaries() {
		//System.out.println("removeSubBoundaries");
		GenomicBoundaries instance = new GenomicBoundaries();
		
		//instance.addBoundary("1", 1, 20);
		//instance.addBoundary("1", 3, 15);
		instance.addBoundary("1", 2, 3);
		instance.addBoundary("1", 1, 4);
		instance.removeSubBoundaries();
		
		boolean expResult = true;
		boolean result = instance.isInBoundary("1", 1, 0);
		assertEquals(result, expResult);
		
		expResult = true;
		result = instance.isInBoundary("1", 2, 0);
		assertEquals(result, expResult);
		
		expResult = true;
		result = instance.isInBoundary("1", 3, 0);
		assertEquals(result, expResult);
		
		expResult = false;
		result = instance.isInBoundary("1", 5, 0);
		assertEquals(result, expResult);
		
		//instance.getBoundary("1", 2, 0);

		int expBoundaryCount = 1;
		int boundaryCount = instance.getBoundaryCountChromosome("1");
		assertEquals(expBoundaryCount, boundaryCount);
	}

	/**
	 * Test of mergeOverlappingBoundaries method, of class GenomicBoundaries.
	 */
	@Test
	public void testMergeOverlappingBoundaries() {
		//System.out.println("mergeOverlappingBoundaries");
		GenomicBoundaries instance = new GenomicBoundaries();
		instance.mergeOverlappingBoundaries();

		instance.addBoundary("1", 1, 3);
		instance.addBoundary("1", 1, 5);
		instance.addBoundary("1", 2, 3);
		instance.addBoundary("1", 3, 4);
		instance.addBoundary("1", 4, 5);
		instance.addBoundary("1", 7, 8);
		instance.addBoundary("1", 9, 10);
		
		int expBoundaryCount = 6;
		int boundaryCount = instance.getBoundaryCountChromosome("1");
		assertEquals(expBoundaryCount, boundaryCount);
		
		instance.mergeOverlappingBoundaries();
		
		expBoundaryCount = 3;
		boundaryCount = instance.getBoundaryCountChromosome("1");
		assertEquals(expBoundaryCount, boundaryCount);
		
		boolean expResult = true;
		boolean result = instance.isInBoundary("1", 1, 0);
		assertEquals(result, expResult);
		
		expResult = true;
		result = instance.isInBoundary("1", 4, 0);
		assertEquals(result, expResult);
		
		expResult = true;
		result = instance.isInBoundary("1", 5, 0);
		assertEquals(result, expResult);
		
		expResult = false;
		result = instance.isInBoundary("1", 6, 0);
		assertEquals(result, expResult);
		
		expResult = false;
		result = instance.isInBoundary("1", 11, 0);
		assertEquals(result, expResult);
		
	}
	
	
	@Test
	public void testIterator(){
		//System.out.println("Testing Iterator");
		String chromosome = "1";
		Integer beginPoint = 1;
		int endPoint = 3;
		
		GenomicBoundaries<Object> instance = new GenomicBoundaries();
		instance.addBoundary(chromosome, 1, 8);
		instance.addBoundary(chromosome, 2, 11);
		instance.addBoundary(chromosome, 3, 14);
		
		int[] expectedResults = new int[3];
		expectedResults[0] = 1;
		expectedResults[1] = 2;
		expectedResults[2] = 3;
		
		int n=0;
		/*for(GenomicBoundary<Object> boundary : instance){
			int k = expectedResults[n];
			int realResult = boundary.getStart();
			System.out.println(k + " :: " + realResult);
			assertEquals(k, realResult);
			n++;
		}*/
	}
	
	
	@Test
	public void testIteratorRemoveBoundaryAndList(){
		//System.out.println("Testing Iterator Remove");
		String chromosome = "1";
		int beginPoint = 1;
		int endPoint = 4;
		
		GenomicBoundaries<Object> instance = new GenomicBoundaries();
		instance.addBoundary(chromosome, beginPoint, endPoint);
		
		
		//boolean expResults = true;
		//boolean result = instance.isInBoundary(chromosome, beginPoint, 0);
		//assertEquals(expResults, result);
		
		
		Iterator<GenomicBoundary<Object>> gboit = instance.iterator();
		while(gboit.hasNext()){
			gboit.next();
			gboit.remove();
		}
		
		
		int expectedSize = 0;
		int actualSize = instance.getBoudaryCount();
		assertEquals(expectedSize, actualSize);
	}
	
	
	@Test
	public void testIteratorRemoveBoundaryNotList(){
		//System.out.println("Testing Iterator Remove only boundary");
		String chrom = "1";
		int start = 1;
		
		
		GenomicBoundaries instance = new GenomicBoundaries();
		instance.addBoundary(chrom, start, 5);
		instance.addBoundary(chrom, start, 64);
		
		/*
		int n=0;
		Iterator<GenomicBoundary<Object>> gboit = instance.iterator();
		while(gboit.hasNext()){
			gboit.next();
			
			if(n == 0){
				gboit.remove();
			}
			n++;
		}
		*/
		
		//Check if only one of two Boundary objects have bee removed from the GenomicBoundaries instance.
		int expectedBoundaryCount = 1;
		int actualBoundaryCount = instance.getBoudaryCount();
		assertEquals(expectedBoundaryCount, actualBoundaryCount);
	}
	
	
	/*
	@Test
	public void testAddBoundaryMore(){
		System.out.println("addBoundary more");
		String chromosome = "1";
		Integer beginPoint = 1;
		int endPoint = 3;
		GenomicBoundaries instance = new GenomicBoundaries();
		instance.addBoundary(chromosome, beginPoint, endPoint);
		instance.addBoundary(chromosome, 2, 5);
		
		boolean expResult = true;
		boolean result = instance.isInBoundary(chromosome, beginPoint, 0);
		assertEquals(result, expResult);
		
		expResult = true;
		result = instance.isInBoundary(chromosome, 2, 0);
		assertEquals(result, expResult);
	}
	
	
	@Test
	public void testAddBoundaryUnderSameStart(){
		System.out.println("addBoundary under the same start");
		String chromosome = "1";
		Integer beginPoint = 1;
		int endPoint = 3;
		GenomicBoundaries instance = new GenomicBoundaries();
		instance.addBoundary(chromosome, beginPoint, endPoint);
		instance.addBoundary(chromosome, beginPoint, 5);
		
		boolean expResult = true;
		boolean result = instance.isInBoundary(chromosome, beginPoint, 0);
		assertEquals(result, expResult);
		
		expResult = true;
		result = instance.isInBoundary(chromosome, beginPoint, 1);
		assertEquals(result, expResult);
	}
	*/
	
	@Test
	public void testIteratorNew(){
		//System.out.println("Testing Iterator New");
		String chromosome = "1";
		Integer beginPoint = 1;
		int endPoint = 3;
		
		GenomicBoundaries<Object> instance = new GenomicBoundaries();
		instance.addBoundary(chromosome, 1, 8);
		instance.addBoundary(chromosome, 1, 12);
		instance.addBoundary(chromosome, 1, 15);
		instance.addBoundary(chromosome, 1, 20);
		instance.addBoundary(chromosome, 2, 11);
		instance.addBoundary(chromosome, 3, 14);
		instance.addBoundary("2", 64, 128);
		
		GenomicBoundary boundary;
		/*
		Iterator<GenomicBoundary<Object>> gboit = instance.iterator();
		while(gboit.hasNext()){
			boundary = gboit.next();
			System.out.println("bd: " + boundary.getStart() + " - " + boundary.getStop());
		}
		
		System.out.println( instance.getBoundaryCountChromosome(chromosome) );*/
	}
}
