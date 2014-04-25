package org.molgenis.vcf;

import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import org.apache.commons.lang3.StringUtils;
import org.molgenis.vcf.meta.VcfMeta;

public class VcfRecord
{	
	private final VcfMeta vcfMeta;
	private String[] tokens;

	private transient List<String> cachedIdentifiers; 		
	private transient List<String> cachedAlternateAlleles;
	private transient String[] cachedSampleDataTypes;
	
	public VcfRecord(VcfMeta vcfMeta) {
		this(vcfMeta, null);
	}
	
	public VcfRecord(VcfMeta vcfMeta, String[] tokens)
	{
		if(vcfMeta == null) throw new IllegalArgumentException("colNames is null");
		this.vcfMeta = vcfMeta;
		this.tokens = tokens;
	}
	
	public String getChromosome() {
		return tokens[VcfMeta.COL_CHROM_IDX];
	}
	
	public int getPosition() {
		return Integer.valueOf(tokens[VcfMeta.COL_POS_IDX]);
	}
	
	public List<String> getIdentifiers() {
		if(cachedIdentifiers == null) {
			String identifiersStr = tokens[VcfMeta.COL_ID_IDX];
			
			if (identifiersStr == null || identifiersStr.equals(".")) {
				cachedIdentifiers = Collections.emptyList();
			} else {
				cachedIdentifiers = Arrays.asList(StringUtils.split(identifiersStr, ';'));
			}
		}
		return cachedIdentifiers;
	}
	
	public String getReferenceAllele() {
		return tokens[VcfMeta.COL_REF_IDX];
	}
	
	public List<String> getAlternateAlleles() {
		if(cachedAlternateAlleles == null) {
			String alternateBasesStr = tokens[VcfMeta.COL_ALT_IDX];
			if(alternateBasesStr.length() == 0) {
				cachedAlternateAlleles = Collections.emptyList();
			} else {
				cachedAlternateAlleles = Arrays.asList(StringUtils.split(alternateBasesStr, ','));
			}
		}
		return cachedAlternateAlleles;
	}
	
	public String getQuality() {
		String quality = tokens[VcfMeta.COL_QUAL_IDX];
		return quality != null && quality.equals(".") ? null : quality;
	}
	
	public String getFilterStatus() {
		return tokens[VcfMeta.COL_FILTER_IDX];
	}
	
	public Iterable<VcfInfo> getInformation() {
		final String[] infoTokens = StringUtils.split(tokens[VcfMeta.COL_INFO_IDX], ';');
		return new Iterable<VcfInfo>() {
			
			@Override
			public Iterator<VcfInfo> iterator()
			{
				return new Iterator<VcfInfo>(){
					private final VcfInfo recycableVcfInfo = new VcfInfo(vcfMeta);
					private int nrToken = 0; 
					
					@Override
					public boolean hasNext()
					{
						return nrToken < infoTokens.length;
					}

					@Override
					public VcfInfo next()
					{
						String infoToken = infoTokens[nrToken++];
						int idx = infoToken.indexOf('=');
						String key = idx != -1 ? infoToken.substring(0, idx) : infoToken;
						String val = idx != -1 ? infoToken.substring(idx + 1) : null;
						recycableVcfInfo.reset(key, val);
						return recycableVcfInfo;
					}

					@Override
					public void remove()
					{
						throw new UnsupportedOperationException();
					}};
			}
		};
	}
	
	public String[] getFormat() {
		if(cachedSampleDataTypes == null) {
			cachedSampleDataTypes = StringUtils.split(tokens[VcfMeta.COL_FORMAT_IDX], ':'); 
		}
		return cachedSampleDataTypes;
	}
	
	public int getFormatIndex(String dataType) {
		String[] dataTypes = getFormat();
		for(int i = 0; i < dataTypes.length; ++i)
			if(dataTypes[i].equals(dataType))
				return i;
		return -1;
	}
	
	public int getNrSamples() {
		return tokens.length > VcfMeta.COL_FORMAT_IDX + 1 ? tokens.length - (VcfMeta.COL_FORMAT_IDX + 1) : 0;
	}
	
	public Iterable<VcfSample> getSamples() {
		final VcfRecord vcfRecord = this;
		return new Iterable<VcfSample>() {

			@Override
			public Iterator<VcfSample> iterator()
			{
				return new Iterator<VcfSample>(){
					private final VcfSample recycableVcfSample = new VcfSample(vcfRecord);
					
					private final int nrSamples = getNrSamples();
					private int nrSample = 0;
					
					@Override
					public boolean hasNext()
					{
						return nrSample < nrSamples;
					}

					@Override
					public VcfSample next()
					{
						recycableVcfSample.reset(StringUtils.split(tokens[VcfMeta.COL_FORMAT_IDX + 1 + nrSample], ':'));
						++nrSample;
						return recycableVcfSample;
					}

					@Override
					public void remove()
					{
						throw new UnsupportedOperationException();
					}};
			}
		};
	}

	public void reset(String[] tokens) {
		this.tokens = tokens;
		this.cachedIdentifiers = null;
		this.cachedAlternateAlleles = null;
		this.cachedSampleDataTypes = null;
	}
	
	public VcfRecord createClone() {
		return new VcfRecord(vcfMeta, tokens);
	}
	
	@Override
	public String toString() {
		return StringUtils.join(tokens, '\t');
	}

	@Override
	public int hashCode()
	{
		final int prime = 31;
		int result = 1;
		result = prime * result + Arrays.hashCode(tokens);
		result = prime * result + ((vcfMeta == null) ? 0 : vcfMeta.hashCode());
		return result;
	}

	@Override
	public boolean equals(Object obj)
	{
		if (this == obj) return true;
		if (obj == null) return false;
		if (getClass() != obj.getClass()) return false;
		VcfRecord other = (VcfRecord) obj;
		if (!Arrays.equals(tokens, other.tokens)) return false;
		if (vcfMeta == null)
		{
			if (other.vcfMeta != null) return false;
		}
		else if (!vcfMeta.equals(other.vcfMeta)) return false;
		return true;
	}
}
