/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.sampleFilter;

import java.util.Collection;
import java.util.HashSet;
import org.molgenis.genotype.Sample;

/**
 *
 * @author Patrick Deelen
 */
public class SampleIdIncludeFilter implements SampleFilter {

	private final HashSet<String> includedSampleIds;

	public SampleIdIncludeFilter(Collection<String> includedSampleIds) {
		this.includedSampleIds = new HashSet<String>(includedSampleIds);
	}
    
    public SampleIdIncludeFilter(String... ids){
		this.includedSampleIds = new HashSet<String>();
		for(String id : ids){
			includedSampleIds.add(id);
		}
	}
	
	public SampleIdIncludeFilter(String sample){
		this.includedSampleIds = new HashSet<String>(1);
		includedSampleIds.add(sample);
	}
    
	public SampleIdIncludeFilter(HashSet<String> includedSampleIds) {
		this.includedSampleIds = includedSampleIds;
	}

	@Override
	public boolean doesSamplePassFilter(Sample sample) {
		return includedSampleIds.contains(sample.getId());
	}
}
