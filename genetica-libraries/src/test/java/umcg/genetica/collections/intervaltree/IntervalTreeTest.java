/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.collections.intervaltree;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import static org.testng.Assert.*;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

/**
 *
 * @author Patrick Deelen
 */
public class IntervalTreeTest {
	
	private final IntervalTree<StringRange> testIntervalTree;
	private final ArrayList<StringRange> elements;
	
	public IntervalTreeTest() {
		
		elements = new ArrayList<StringRange>();
		
		//elements.add(new StringRange(, , ""));
		elements.add(new StringRange(0, 10, "0-10a"));
		elements.add(new StringRange(5, 10, "5-10"));
		elements.add(new StringRange(10, 11, "10-11"));
		elements.add(new StringRange(10, 10, "10-10"));
		elements.add(new StringRange(0, 20, "0-20"));
		elements.add(new StringRange(5, 10, "5-10b"));
		
		testIntervalTree = new IntervalTree<StringRange>(elements, StringRange.class);
		
	}

	/**
	 * Test of getElementsOverlappingQuery method, of class IntervalTree.
	 */
	@Test
	public void testGetElementsOverlappingQuery() {
		
		ArrayList result = testIntervalTree.getElementsOverlappingQuery(0);
		assertEquals(result.size(), 2);
		assertTrue(result.contains(elements.get(0)));
		assertTrue(result.contains(elements.get(4)));
		
		result = testIntervalTree.getElementsOverlappingQuery(2);
		assertEquals(result.size(), 2);
		assertTrue(result.contains(elements.get(0)));
		assertTrue(result.contains(elements.get(4)));
		
		result = testIntervalTree.getElementsOverlappingQuery(5);
		assertEquals(result.size(), 4);
		assertTrue(result.contains(elements.get(0)));
		assertTrue(result.contains(elements.get(4)));
		assertTrue(result.contains(elements.get(1)));
		assertTrue(result.contains(elements.get(5)));
		
		result = testIntervalTree.getElementsOverlappingQuery(10);
		assertEquals(result.size(), 6);
		assertTrue(result.contains(elements.get(0)));
		assertTrue(result.contains(elements.get(4)));
		assertTrue(result.contains(elements.get(1)));
		assertTrue(result.contains(elements.get(5)));
		assertTrue(result.contains(elements.get(2)));
		assertTrue(result.contains(elements.get(3)));
		
		result = testIntervalTree.getElementsOverlappingQuery(11);
		assertEquals(result.size(), 2);
		assertTrue(result.contains(elements.get(4)));
		assertTrue(result.contains(elements.get(2)));
		
		result = testIntervalTree.getElementsOverlappingQuery(25);
		assertEquals(result.size(), 0);
		
		result = testIntervalTree.getElementsOverlappingQuery(20);
		assertEquals(result.size(), 1);
		assertTrue(result.contains(elements.get(4)));
		
	}


}