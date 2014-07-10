/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.io.bedgraph;

import java.io.File;
import java.util.Iterator;
import java.util.List;
import static org.testng.Assert.*;
import org.testng.annotations.Test;
import umcg.genetica.collections.intervaltree.PerChrIntervalTree;

/**
 *
 * @author Patrick Deelen
 */
public class BedGraphFileTest {
	
	public BedGraphFileTest() {
		
		
		
	}

	/**
	 * Test of createIntervalTree method, of class BedGraphFile.
	 */
	@Test
	public void testCreateIntervalTree() throws Exception {
		
		BedGraphFile testBedGraphFile = new BedGraphFile(new File(this.getClass().getResource("/wgEncodeCrgMapabilityAlign50mer.bedGraph").toURI()), true, true);
		
		PerChrIntervalTree<BedGraphEntry> bedGraphIntervalTree = testBedGraphFile.createIntervalTree();
		
		testIntervalForMapability(bedGraphIntervalTree.searchPosition("1", 10001), 0.00172712);
		testIntervalForMapability(bedGraphIntervalTree.searchPosition("1", 10062), 0.00172712);
		testIntervalForMapability(bedGraphIntervalTree.searchPosition("1", 10064), 0.00172712);
		testIntervalForMapability(bedGraphIntervalTree.searchPosition("1", 10102), 0.0769231);
		testIntervalForMapability(bedGraphIntervalTree.searchPosition("1", 10102), 0.0769231);
		testIntervalForMapability(bedGraphIntervalTree.searchPosition("1", 10102), 0.0769231);
		testIntervalForMapability(bedGraphIntervalTree.searchPosition("2", 1), 0.0384615);
		testIntervalForMapability(bedGraphIntervalTree.searchPosition("2", 2), 0.0384615);
		testIntervalForMapability(bedGraphIntervalTree.searchPosition("2", 10), 0.0384615);
		testIntervalForMapability(bedGraphIntervalTree.searchPosition("2", 10108), 0.047619);
		testIntervalForMapability(bedGraphIntervalTree.searchPosition("2", 10121), 0.047619);
		
		assertTrue(bedGraphIntervalTree.searchPosition("1", 10000).isEmpty());
		assertTrue(bedGraphIntervalTree.searchPosition("1", 10103).isEmpty());
		assertTrue(bedGraphIntervalTree.searchPosition("2", 0).isEmpty());
		assertTrue(bedGraphIntervalTree.searchPosition("2", 11).isEmpty());
		assertTrue(bedGraphIntervalTree.searchPosition("3", 1).isEmpty());
		
	}
	
	private void testIntervalForMapability(List<BedGraphEntry> entries, double expectedMapability){
		assertEquals(entries.size(), 1);
		assertEquals(entries.get(0).getValue(), expectedMapability);
	}
}