/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.sampleFilter;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import org.molgenis.genotype.Sample;

/**
 *
 * @author Patrick Deelen
 */
public class SampleCombinedFilter implements SampleFilter {

	private final List<SampleFilter> sampleFilters;

	public SampleCombinedFilter(SampleFilter... sampleFilters) {
		this.sampleFilters = Arrays.asList(sampleFilters);
	}

	public SampleCombinedFilter(List<SampleFilter> sampleFilters) {
		this.sampleFilters = sampleFilters == null ? new ArrayList<SampleFilter>() : sampleFilters;
	}

	@Override
	public boolean doesSamplePassFilter(Sample sample) {

		for (SampleFilter filter : sampleFilters) {
			if (!filter.doesSamplePassFilter(sample)) {
				return false;
			}
		}
		return true;

	}
}
