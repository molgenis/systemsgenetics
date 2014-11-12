/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.io.fasta;

import static org.testng.Assert.*;
import org.testng.annotations.Test;

/**
 *
 * @author Patrick Deelen
 */
public class LargeByteArrayTest {
	
	public LargeByteArrayTest() {
	}

	/**
	 * Test of setQuick method, of class LargeByteArray.
	 */
	@Test
	public void testSetGetQuick() {

		long index;
		byte value;

		
		LargeByteArray instance = new LargeByteArray(10);		
		
		index = 0L;
		value = 1;
		instance.setQuick(index, value);
		assertEquals(instance.getQuick(index), value);
		
		index = 0L;
		value = 2;
		instance.setQuick(index, value);
		assertEquals(instance.getQuick(index), value);

		index = 1L;
		value = 2;
		instance.setQuick(index, value);
		assertEquals(instance.getQuick(index), value);
		
		index = 9L;
		value = 10;
		instance.setQuick(index, value);
		assertEquals(instance.getQuick(index), value);
		
		
		
		
	}


}