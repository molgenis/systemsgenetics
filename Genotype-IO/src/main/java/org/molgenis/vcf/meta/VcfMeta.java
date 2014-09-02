package org.molgenis.vcf.meta;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

public class VcfMeta
{	
	private static final String KEY_FILEFORMAT = "fileformat";
	
	public static final int COL_CHROM_IDX = 0;
	public static final int COL_POS_IDX = 1;
	public static final int COL_ID_IDX = 2;
	public static final int COL_REF_IDX = 3;
	public static final int COL_ALT_IDX = 4;
	public static final int COL_QUAL_IDX = 5;
	public static final int COL_FILTER_IDX = 6;
	public static final int COL_INFO_IDX = 7;
	public static final int COL_FORMAT_IDX = 8;
	
	private final Map<String, String> meta;
	private Map<String, VcfMetaAlt> vcfMetaAlts;
	private Map<String, VcfMetaContig> vcfMetaContigs;
	private Map<String, VcfMetaFilter> vcfMetaFilters;
	private Map<String, VcfMetaFormat> vcfMetaFormats;
	private Map<String, VcfMetaInfo> vcfMetaInfos;
	private Map<String, VcfMetaSample> vcfMetaSamples;
	private List<VcfMetaPedigree> vcfMetaPedigrees;
	private String[] colNames;
	
	public VcfMeta() {
		this.meta = new HashMap<String, String>();
	}
	
	public String getFileFormat() {
		return meta.get(KEY_FILEFORMAT);
	}
	
	public void add(String key, String value) {
		meta.put(key, value);
	}
	
	public String get(String key) {
		return meta.get(key);
	}
	
	public void setColNames(String[] colNames) {
		this.colNames = colNames;
	}
	
	public String[] getColNames() {
		return colNames;
	}
	
	public String getSampleName(int nr) {
		return colNames[COL_FORMAT_IDX + 1 + nr];
	}
	
	public Iterable<String> getSampleNames() {
		if(colNames.length <= COL_FORMAT_IDX ){
			return Collections.emptyList();
		}
		return Arrays.asList(colNames).subList(COL_FORMAT_IDX + 1, colNames.length);
	}
	
	public void addAltMeta(VcfMetaAlt vcfMetaAlt) {
		if(vcfMetaAlts == null) vcfMetaAlts = new LinkedHashMap<String, VcfMetaAlt>();
		vcfMetaAlts.put(vcfMetaAlt.getId(), vcfMetaAlt);
	}
	
	public Iterable<VcfMetaAlt> getAltMeta() {
		return vcfMetaAlts != null ? vcfMetaAlts.values() : Collections.<VcfMetaAlt>emptyList();
	}
	
	public VcfMetaAlt getAltMeta(String id) {
		return vcfMetaAlts != null ? vcfMetaAlts.get(id) : null;
	}
	
	public void addContigMeta(VcfMetaContig vcfMetaContig) {
		if(vcfMetaContigs == null) vcfMetaContigs = new LinkedHashMap<String, VcfMetaContig>();
		vcfMetaContigs.put(vcfMetaContig.getId(), vcfMetaContig);
	}
	
	public Iterable<VcfMetaContig> getContigMeta() {
		return vcfMetaContigs != null ? vcfMetaContigs.values() : Collections.<VcfMetaContig>emptyList();
	}
	
	public VcfMetaContig getContigMeta(String id) {
		return vcfMetaContigs != null ? vcfMetaContigs.get(id) : null;
	}
	
	public void addFilterMeta(VcfMetaFilter vcfMetaFilter) {
		if(vcfMetaFilters == null) vcfMetaFilters = new LinkedHashMap<String, VcfMetaFilter>();
		vcfMetaFilters.put(vcfMetaFilter.getId(), vcfMetaFilter);
	}
	
	public Iterable<VcfMetaFilter> getFilterMeta() {
		return vcfMetaFilters != null ? vcfMetaFilters.values() : Collections.<VcfMetaFilter>emptyList();
	}
	
	public VcfMetaFilter getFilterMeta(String id) {
		return vcfMetaFilters != null ? vcfMetaFilters.get(id) : null;
	}
	
	public void addFormatMeta(VcfMetaFormat vcfMetaFormat) {
		if(vcfMetaFormats == null) vcfMetaFormats = new LinkedHashMap<String, VcfMetaFormat>();
		vcfMetaFormats.put(vcfMetaFormat.getId(), vcfMetaFormat);
	}
	
	public Iterable<VcfMetaFormat> getFormatMeta() {
		return vcfMetaFormats != null ? vcfMetaFormats.values() : Collections.<VcfMetaFormat>emptyList();
	}
	
	public VcfMetaFormat getFormatMeta(String id) {
		return vcfMetaFormats != null ? vcfMetaFormats.get(id) : null;
	}
	
	public void addInfoMeta(VcfMetaInfo vcfMetaInfo) {
		if(vcfMetaInfos == null) vcfMetaInfos = new LinkedHashMap<String, VcfMetaInfo>();
		vcfMetaInfos.put(vcfMetaInfo.getId(), vcfMetaInfo);
	}
	
	public Iterable<VcfMetaInfo> getInfoMeta() {
		return vcfMetaInfos != null ? vcfMetaInfos.values() : Collections.<VcfMetaInfo>emptyList();
	}
	
	public VcfMetaInfo getInfoMeta(String id) {
		return vcfMetaInfos != null ? vcfMetaInfos.get(id) : null;
	}
	
	public void addPedigreeMeta(VcfMetaPedigree vcfMetaPedigree) {
		if(vcfMetaPedigrees == null) vcfMetaPedigrees = new ArrayList<VcfMetaPedigree>();
		vcfMetaPedigrees.add(vcfMetaPedigree);
	}
	
	public Iterable<VcfMetaPedigree> getPedigreeMeta() {
		return vcfMetaPedigrees != null ? vcfMetaPedigrees : Collections.<VcfMetaPedigree>emptyList();
	}
	
	public void addSampleMeta(VcfMetaSample vcfMetaSample) {
		if(vcfMetaSamples == null) vcfMetaSamples = new LinkedHashMap<String, VcfMetaSample>();
		vcfMetaSamples.put(vcfMetaSample.getId(), vcfMetaSample);
	}
	
	public Iterable<VcfMetaSample> getSampleMeta() {
		return vcfMetaSamples != null ? vcfMetaSamples.values() : Collections.<VcfMetaSample>emptyList();
	}
	
	public VcfMetaSample getSampleMeta(String id) {
		return vcfMetaSamples != null ? vcfMetaSamples.get(id) : null;
	}
	
	@Override
	public String toString() {
		StringBuilder strBuilder = new StringBuilder();
		for(Map.Entry<String, String> entry : meta.entrySet()) {
			strBuilder.append(entry.getKey()).append('=').append(entry.getValue()).append('\n');
		}
		if(vcfMetaAlts != null) {
			for(VcfMetaAlt vcfMetaAlt : vcfMetaAlts.values())
				strBuilder.append(vcfMetaAlt);
		}
		if(vcfMetaContigs != null) {
			for(VcfMetaContig vcfMetaContig : vcfMetaContigs.values())
				strBuilder.append(vcfMetaContig);
		}
		if(vcfMetaFilters != null) {
			for(VcfMetaFilter vcfMetaFilter : vcfMetaFilters.values())
				strBuilder.append(vcfMetaFilter);
		}
		if(vcfMetaFormats != null) {
			for(VcfMetaFormat vcfMetaFormat : vcfMetaFormats.values())
				strBuilder.append(vcfMetaFormat);
		}
		if(vcfMetaInfos != null) {
			for(VcfMetaInfo vcfMetaInfo : vcfMetaInfos.values())
				strBuilder.append(vcfMetaInfo);
		}
		if(vcfMetaSamples != null) {
			for(VcfMetaSample vcfMetaSample : vcfMetaSamples.values())
				strBuilder.append(vcfMetaSample);
		}
		if(vcfMetaPedigrees != null) {
			for(VcfMetaPedigree vcfMetaPedigree : vcfMetaPedigrees)
				strBuilder.append(vcfMetaPedigree);
		}
		return strBuilder.toString();
	}
}
