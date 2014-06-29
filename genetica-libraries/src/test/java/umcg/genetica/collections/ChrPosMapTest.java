/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.collections;

import static org.testng.Assert.*;
import org.testng.annotations.Test;

/**
 *
 * @author Patrick Deelen
 */
public class ChrPosMapTest {
	
	public ChrPosMapTest() {
	}

	/**
	 * Test of put method, of class ChrPosMap.
	 */
	@Test
	public void test() {
		ChrPosMap<String> map = new ChrPosMap<String>();
		
		assertNull(map.get("1", 0));
		
		map.put("1", 1, "Test1");
		assertEquals(map.get("1", 1), "Test1");
		assertNull(map.get("1", 0));
		assertNull(map.get("2", 1));
		
		map.put("1", 0, "Test2");
		assertEquals(map.get("1", 1), "Test1");
		assertEquals(map.get("1", 0), "Test2");
		assertNull(map.get("2", 0));
		
		map.put("2", 1, "Test3");
		assertEquals(map.get("1", 1), "Test1");
		assertEquals(map.get("1", 0), "Test2");
		assertEquals(map.get("2", 1), "Test3");
		
		map.put("1", 0, "Test4");
		assertEquals(map.get("1", 1), "Test1");
		assertEquals(map.get("1", 0), "Test4");
		assertEquals(map.get("2", 1), "Test3");
		
	}

}