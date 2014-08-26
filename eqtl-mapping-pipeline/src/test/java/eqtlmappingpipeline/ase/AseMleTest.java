/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.ase;

import cern.colt.list.tint.IntArrayList;
import java.util.Random;
import static org.testng.Assert.*;
import org.testng.annotations.Test;

/**
 *
 * @author Patrick Deelen
 */
public class AseMleTest {

	public AseMleTest() {
	}

	@Test
	public void test() {

		Random random = new Random(1);
		
		final IntArrayList a1Counts = new IntArrayList();
		final IntArrayList a2Counts = new IntArrayList();

		int nrSamples = 100;
		for (int s = 0; s < nrSamples; s++) {
			int nrReads = 1000 + (int) (random.nextDouble() * 1000d); //Simulate random number of reads, assume this is unrelated to tissue
			int readsWT = 0;
			int readsAlt = 0;
			for (int r = 0; r < nrReads; r++) {
				if (random.nextDouble() < 0.6d) {
					readsWT++;
				} else {
					readsAlt++;
				}
			}

			a1Counts.add(readsWT);
			a2Counts.add(readsAlt);

		}

		AseMle mle = new AseMle(a1Counts, a2Counts);

		assertEquals(mle.getMaxLikelihood(), -428.3175944807302, 0.00001);
		assertEquals(mle.getMaxLikelihoodP(), 0.599, 0.00001);
		assertEquals(mle.getRatioD(), 5966.782310979877, 0.00001);
		assertEquals(mle.getRatioP(), 0, 0.00001);

	}

	@Test
	public void testLnBico() {

		assertEquals(AseMle.lnbico(1369, 689), 945.052036055064, 0.000001);
		assertEquals(AseMle.lnbico(12345, 6789), 8490.2927640914, 0.000001);
		assertEquals(AseMle.lnbico(123456, 67890), 84950.948114749, 0.000001);
		assertEquals(AseMle.lnbico(20, 10), 12.126791314602454, 0.000001);
		assertEquals(AseMle.lnbico(200, 10), 37.6501117434664266, 0.000001);
	}
	
	@Test
	public void preCalculateP(){
		
		int expectedLenght = 999;
		assertEquals(AseMle.probabilities.length, expectedLenght);
		assertEquals(AseMle.logProbabilities.length, expectedLenght);
		assertEquals(AseMle.log1minProbabilities.length, expectedLenght);
		
		assertEquals(AseMle.probabilities[0], 0.001);
		assertEquals(AseMle.probabilities[499], 0.5);
		assertEquals(AseMle.probabilities[998], 0.999);
		
		assertEquals(AseMle.logProbabilities[0], -6.907755, 0.00001);
		assertEquals(AseMle.logProbabilities[499], -0.693147, 0.00001);
		assertEquals(AseMle.logProbabilities[998], -0.00100050, 0.00001);
		
		assertEquals(AseMle.log1minProbabilities[0], -0.00100050, 0.00001);
		assertEquals(AseMle.log1minProbabilities[499], -0.693147, 0.00001);
		assertEquals(AseMle.log1minProbabilities[998], -6.907755, 0.00001);
		
	}
}