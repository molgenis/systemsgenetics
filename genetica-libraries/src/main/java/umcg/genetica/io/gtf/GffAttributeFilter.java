/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.io.gtf;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;

/**
 *
 * @author Patrick Deelen
 */
public class GffAttributeFilter {
	
	protected final Map<String, ArrayList<String>> filter;

	public GffAttributeFilter() {
		this.filter = new HashMap<String, ArrayList<String>>();
	}
	
	public void addFilter(String attribute, String value){
		
		ArrayList<String> attributeFilterValues;
		
		if(filter.containsKey(attribute)){
			attributeFilterValues = filter.get(attribute);
		} else {
			attributeFilterValues = new ArrayList<String>();
			filter.put(attribute, attributeFilterValues);
		}
		
		attributeFilterValues.add(value);
		
	}
	
	public boolean doesMatchFilter(GffElement gffElement){
		
		attributesFilters:
		for(Entry<String, ArrayList<String>> attributeFilterEntry : filter.entrySet()){
			
			String attributeName = attributeFilterEntry.getKey();
			ArrayList<String> attributeFilterValues = attributeFilterEntry.getValue();
			
			if(!gffElement.hasAttribute(attributeName)){
				return false;
			}
			
			String gtfElementAttributeValue = gffElement.getAttributeValue(attributeName);
			
			for(String attributeFilterValue : attributeFilterValues){
				if(attributeFilterValue.equals(gtfElementAttributeValue)){
					continue attributesFilters;
				}
			}
			
			//The programs only reaches this when there as been not value match for this attribute
			return false;
			
		}
		
		return true;
		
	}
	
}
