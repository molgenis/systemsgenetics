package betaqtl.vcf.filter.variantfilters;


import betaqtl.vcf.VCFVariant;

import java.util.ArrayList;

/**
 * Created by hwestra on 3/21/17.
 */
public class VCFVariantFilters {
	
	ArrayList<VCFVariantFilter> filters;
	public boolean hasRegionOrSetFilter = false;
	
	boolean filterschecked = false;
	
	public VCFVariantFilters() {
		filters = new ArrayList<VCFVariantFilter>();
		
	}
	
	public void add(VCFVariantFilter filter) {
		addFilter(filter);
		checkRegionFilter(filter);
	}
	
	public void addFilter(VCFVariantFilter filter) {
		filters.add(filter);
		checkRegionFilter(filter);
	}
	
	private void checkRegionFilter(VCFVariantFilter f) {
		
		if (f instanceof VCFVariantSetFilter || f instanceof VCFVariantRegionFilter) {
			hasRegionOrSetFilter = true;
		}
		
	}
	
	VCFVariantFilter failedFilter = null;
	
	public boolean passesFilters(VCFVariant v) {
		for (VCFVariantFilter f : filters) {
			if (!f.passesThreshold(v)) {
				failedFilter = f;
				return false;
			}
		}
		failedFilter = null;
		return true;
	}
	
	public VCFVariantFilter getFailedFilter() {
		return failedFilter;
	}
	
	public String toString() {
		String output = "VCF Variant filters:\n";
		for (VCFVariantFilter filter : filters) {
			output += filter.toString() + "\n";
		}
		return output;
	}
	
	public int size() {
		return filters.size();
	}
	
	public boolean passesRegionOrVariantFilter(VCFVariant v) {
		boolean passesfilter = true;
		for (VCFVariantFilter f : filters) {
			if (f instanceof VCFVariantSetFilter) {
				passesfilter = f.passesThreshold(v);
			}
		}
		return passesfilter;
		
	}
	
	
	public boolean hasRegionOrVariantSetFilter() {
		
		return hasRegionOrSetFilter;
	}
	
	public ArrayList<VCFVariantFilter> getFilters() {
		return filters;
	}
}

