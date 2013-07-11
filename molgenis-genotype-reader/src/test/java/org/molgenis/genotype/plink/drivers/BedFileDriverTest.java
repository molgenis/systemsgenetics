package org.molgenis.genotype.plink.drivers;

import static org.testng.AssertJUnit.assertEquals;

import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

/**
 * Looks like the PED driver test, except gets the alleles only, and in binary encoding
 * 
 * PLEASE NOTE THAT:
 * The driver DOES NOT know about the padding bits, because the BED file alone does
 * not know about the amount of individuals or SNPs. Therefore, although the test file
 * contains 90 SNP elements, the BED driver will return 120. If you want to retrieve 
 * ONLY the REAL genotypes, an additional layer of logic is needed to leave out the
 * elements that are padding and not real data.
 * 
 * HOWEVER: the elements are retrieved in the 'correct' order: in every byte, read
 * from 'right to left' starting at byte nr. 4 (the first 3 are reserved, see plink docs)
 * 
 * ALSO NOTE: this file is "SNP-major" oriented, meaning:
 * it lists ALL individuals for 1 SNP, then moves on to the next SNP, and so on
 * 
 * see: http://pngu.mgh.harvard.edu/~purcell/plink/binary.shtml
 * 
 * binary dump of the test file:
 * 
	0000000: 01101100 00011011 00000001 11001011 10111111 00000010  l.....
	0000006: 11111011 11111111 00000011 11111111 11111111 00000011  ......
	000000c: 11111111 11111111 00000011 11111011 11111011 00000010  ......
	0000012: 11111111 10111111 00000011 11111011 11111111 00000011  ......
	0000018: 11111111 11111111 00000011 11111111 11111111 00000010  ......
	000001e: 11111111 11111111 00000011                             ...
 * 
 * NOTE: the bytes with '000000' as first 6 values have the padding
 * 
 * @author jvelde
 *
 */
public class BedFileDriverTest extends AbstractResourceTest
{
	private BedFileDriver bedfd;

	@BeforeClass
	public void setup() throws Exception
	{
		bedfd = new BedFileDriver(getTestResource("/test.bed"));
	}

	@Test
	public void BED_construct() throws Exception
	{
		assertEquals(1, bedfd.getMode());
		assertEquals(120, bedfd.getNrOfElements());
	}

	@Test
	public void BED_getElement() throws Exception
	{
		//the first SNP
		assertEquals("11", bedfd.getElement(0));
		assertEquals("01", bedfd.getElement(1));
		assertEquals("00", bedfd.getElement(2));
		assertEquals("11", bedfd.getElement(3));
		assertEquals("11", bedfd.getElement(4));
		assertEquals("11", bedfd.getElement(5));
		assertEquals("11", bedfd.getElement(6));
		assertEquals("01", bedfd.getElement(7));
		assertEquals("01", bedfd.getElement(8));
		
		//padding bits
		assertEquals("00", bedfd.getElement(9));
		assertEquals("00", bedfd.getElement(10));
		assertEquals("00", bedfd.getElement(11));
		
		//the second SNP
		assertEquals("11", bedfd.getElement(12));
		assertEquals("01", bedfd.getElement(13));
		assertEquals("11", bedfd.getElement(14));
		assertEquals("11", bedfd.getElement(15));

	}
	
	@Test
	public void BED_getElements() throws Exception
	{
		assertEquals(9, bedfd.getSNPs(0, 9).length);
		assertEquals(9, bedfd.getSNPs(1, 9).length);
		assertEquals(9, bedfd.getSNPs(2, 9).length);
		
		//indv 1, first 3
		assertEquals("11", bedfd.getSNPs(0, 9)[0]);
		assertEquals("01", bedfd.getSNPs(0, 9)[1]);
		assertEquals("00", bedfd.getSNPs(0, 9)[2]);
		
		//indv 1, last 3
		assertEquals("11", bedfd.getSNPs(0, 9)[6]);
		assertEquals("01", bedfd.getSNPs(0, 9)[7]);
		assertEquals("01", bedfd.getSNPs(0, 9)[8]);
		
		//indv 2, first 3
		assertEquals("11", bedfd.getSNPs(1, 9)[0]);
		assertEquals("01", bedfd.getSNPs(1, 9)[1]);
		assertEquals("11", bedfd.getSNPs(1, 9)[2]);
		
		//indv 2, last 3
		assertEquals("11", bedfd.getSNPs(1, 9)[6]);
		assertEquals("11", bedfd.getSNPs(1, 9)[7]);
		assertEquals("11", bedfd.getSNPs(1, 9)[8]);
		
		//indv 10, first 3
		assertEquals("11", bedfd.getSNPs(9, 9)[0]);
		assertEquals("11", bedfd.getSNPs(9, 9)[1]);
		assertEquals("11", bedfd.getSNPs(9, 9)[2]);
		
		//indv 10, last 3
		assertEquals("11", bedfd.getSNPs(9, 9)[6]);
		assertEquals("11", bedfd.getSNPs(9, 9)[7]);
		assertEquals("11", bedfd.getSNPs(9, 9)[8]);
		
		
	}
}
