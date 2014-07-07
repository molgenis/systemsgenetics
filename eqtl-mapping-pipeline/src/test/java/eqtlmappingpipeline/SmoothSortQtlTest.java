/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline;

import eqtlmappingpipeline.metaqtl3.containers.QTL;
import java.util.Random;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;
import umcg.genetica.util.SmoothSort;
import static org.testng.Assert.*;

/**
 *
 * @author MarcJan
 */
public class SmoothSortQtlTest {

	private static final Random random = new Random(1);
	
	public SmoothSortQtlTest() {
	}

	@BeforeMethod
	public void setUpMethod() throws Exception {
	}

	/**
	 * Test of sort method, of class SmoothSort.
	 */
	@Test
	public void testSortBig() {

		System.out.println("Testing sort");
		
		QTL[] eqtls = new QTL[100000];
		
		int numberToSort = 90000;

		for (int i = 0; i < eqtls.length ; ++i) {
			
			double pvalue = random.nextDouble();
			double zscore = random.nextDouble();
			
			eqtls[i] = new QTL(pvalue, 0, 0, (byte) 0, zscore, null, null, null, null, null, null, null, 0d, 0d);
		}
		
		
		SmoothSort.sort(eqtls, 0, numberToSort);
		
		double lastPvalue = Double.MIN_VALUE;
		
		for (int i = 0; i < numberToSort ; ++i) {
			
			assertFalse(eqtls[i].getPvalue() < lastPvalue, "Index " + i + " last p-value " + lastPvalue + " current p-value " + eqtls[i].getPvalue());
			lastPvalue = eqtls[i].getPvalue();
			
		}
		
		
		
	}
}