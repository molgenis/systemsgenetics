/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.sampleFilter;

import org.molgenis.genotype.Sample;

/**
 *
 * @author Patrick Deelen
 */
public interface SampleFilter {
	
	public boolean doesSamplePassFilter(Sample sample);
	
}
