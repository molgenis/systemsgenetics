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
public class AseMleNegativeTest {

	public AseMleNegativeTest() {
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

		AseMleNegative mle = new AseMleNegative(a1Counts, a2Counts);

		//assertEquals(mle.getMaxLikelihood(), -428.3175944807302, 0.00001);
		assertEquals(mle.getMaxLikelihoodP(), 0.599, 0.00001);
		assertEquals(mle.getRatioD(), 5966.782310979877, 0.00001);
		assertEquals(mle.getRatioP(), 0, 0.00001);

	}
	
	@Test
	public void test2() {
		
		final IntArrayList a1Counts = new IntArrayList();
		final IntArrayList a2Counts = new IntArrayList();
		
		a1Counts.add(10);
		a2Counts.add(12);
		
		AseMle mle = new AseMle(a1Counts, a2Counts);
		
		
		System.out.println(mle.getRatioP());
		System.out.println(mle.getMaxLikelihoodP());
		System.out.println(mle.getMaxLikelihood());
		System.out.println(mle.getRatioD());
		
		System.out.println("---");
		AseMleNegative mleNegative = new AseMleNegative(a1Counts, a2Counts);
		
		System.out.println(mleNegative.getRatioP());
		System.out.println(mleNegative.getMaxLikelihoodP());
		System.out.println(mleNegative.getMaxLikelihood());
		System.out.println(mleNegative.getRatioD());
		
	}

}