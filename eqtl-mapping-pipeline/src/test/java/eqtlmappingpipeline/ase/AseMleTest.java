/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.ase;

import cern.colt.list.tint.IntArrayList;
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

		final IntArrayList a1Counts = new IntArrayList();
		final IntArrayList a2Counts = new IntArrayList();
		
		int nrSamples = 100;
		for (int s = 0; s < nrSamples; s++) {
			int nrReads = 10 + (int) (Math.random() * 1000d); //Simulate random number of reads, assume this is unrelated to tissue
			int readsWT = 0;
			int readsAlt = 0;
			for (int r = 0; r < nrReads; r++) {
				if (Math.random() < 0.6d) {
					readsWT++;
				} else {
					readsAlt++;
				}
			}
			
			a1Counts.add(readsWT);
			a2Counts.add(readsAlt);
			
		}
		
		AseMle mle = new AseMle(a1Counts, a2Counts);
		
		System.out.println(mle.getMaxLikelihood());
		System.out.println(mle.getMaxLikelihoodP());
		System.out.println(mle.getRatioD());
		System.out.println(mle.getRatioP());

	}

	
}