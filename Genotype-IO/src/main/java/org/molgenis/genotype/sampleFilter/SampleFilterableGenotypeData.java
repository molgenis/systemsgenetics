/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.sampleFilter;

import java.util.List;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.Sample;

/**
 *
 * @author Patrick Deelen
 */
public interface SampleFilterableGenotypeData extends RandomAccessGenotypeData {

	public List<Sample> getOriginalSampleList();

	public SampleFilter getSampleFilter();

	public int getIncludedSampleCount();
}
