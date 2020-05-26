/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.io.regulomedb;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Pattern;

/**
 *
 * @author Patrick Deelen
 */
public class RegulomeDbEntry {
	
	private static final Pattern TAB_PATTERN = Pattern.compile("\\t");
	private static final Pattern SUPPORT_SEP_PATTERN = Pattern.compile(", ");
	
	private final String chr;
	private final int chrPos;
	private final String variantId;
	private final String regulomeDbScore;
	private final Map<String, List<RegulomeDbSupportingData>> supportData;

	public RegulomeDbEntry(String regulomeDbFileLine) throws Exception {
		
		String[] lineElements = TAB_PATTERN.split(regulomeDbFileLine);
		
		if(lineElements.length != 5){
			throw new Exception("Error in RegulomeDB file. Expected 5 columns but found this: " + regulomeDbFileLine);
		}
		
		chr = lineElements[0];
		chrPos = Integer.parseInt(lineElements[1]);
		variantId = lineElements[2];
		regulomeDbScore = lineElements[4];
		
		HashMap<String, ArrayList<RegulomeDbSupportingData>> supportDataLoader = new HashMap<String, ArrayList<RegulomeDbSupportingData>>();
		
		if(!lineElements[3].equals(".")){

			for(String supportDataElementString : SUPPORT_SEP_PATTERN.split(lineElements[3])){

				RegulomeDbSupportingData supportDataElement = new RegulomeDbSupportingData(supportDataElementString);

				ArrayList<RegulomeDbSupportingData> supportDataClassElements;
				if(supportDataLoader.containsKey(supportDataElement.getSupportClass())){
					supportDataClassElements = supportDataLoader.get(supportDataElement.getSupportClass());
				} else {
					supportDataClassElements = new ArrayList<RegulomeDbSupportingData>();
					supportDataLoader.put(supportDataElement.getSupportClass(), supportDataClassElements);
				}
				supportDataClassElements.add(supportDataElement);

			}
		}
			
		HashMap<String, List<RegulomeDbSupportingData>> supportDataLoaderUnmodifiableLists = new HashMap<String, List<RegulomeDbSupportingData>>();
		
		for(Map.Entry<String, ArrayList<RegulomeDbSupportingData>> supportDataLoaderEntry : supportDataLoader.entrySet()){
			supportDataLoaderUnmodifiableLists.put(supportDataLoaderEntry.getKey(), Collections.unmodifiableList(supportDataLoaderEntry.getValue()));
		}
		
		supportData = Collections.unmodifiableMap(supportDataLoaderUnmodifiableLists);
		
	}

	public String getChr() {
		return chr;
	}

	public int getChrPos() {
		return chrPos;
	}

	public String getVariantId() {
		return variantId;
	}

	public String getRegulomeDbScore() {
		return regulomeDbScore;
	}

	public Map<String, List<RegulomeDbSupportingData>> getSupportData() {
		return supportData;
	}
	
	public boolean hasSupportDataClass(String supportDataClass){
		return supportData.containsKey(supportDataClass);
	}

	@Override
	public int hashCode() {
		int hash = 7;
		hash = 61 * hash + (this.chr != null ? this.chr.hashCode() : 0);
		hash = 61 * hash + this.chrPos;
		hash = 61 * hash + (this.variantId != null ? this.variantId.hashCode() : 0);
		hash = 61 * hash + (this.regulomeDbScore != null ? this.regulomeDbScore.hashCode() : 0);
		hash = 61 * hash + (this.supportData != null ? this.supportData.hashCode() : 0);
		return hash;
	}

	@Override
	public boolean equals(Object obj) {
		if (obj == null) {
			return false;
		}
		if (getClass() != obj.getClass()) {
			return false;
		}
		final RegulomeDbEntry other = (RegulomeDbEntry) obj;
		if ((this.chr == null) ? (other.chr != null) : !this.chr.equals(other.chr)) {
			return false;
		}
		if (this.chrPos != other.chrPos) {
			return false;
		}
		if ((this.variantId == null) ? (other.variantId != null) : !this.variantId.equals(other.variantId)) {
			return false;
		}
		if ((this.regulomeDbScore == null) ? (other.regulomeDbScore != null) : !this.regulomeDbScore.equals(other.regulomeDbScore)) {
			return false;
		}
		if (this.supportData != other.supportData && (this.supportData == null || !this.supportData.equals(other.supportData))) {
			return false;
		}
		return true;
	}
	
	
	
	
}
