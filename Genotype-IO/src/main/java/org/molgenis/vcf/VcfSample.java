package org.molgenis.vcf;

import com.google.common.base.CharMatcher;
import java.util.ArrayList;
import java.util.List;
import org.molgenis.genotype.Allele;

public class VcfSample
{
	private static final String FIELD_GT = "GT";
	private static final char GENOTYPE_UNPHASED = '/';
	private static final char GENOTYPE_PHASED = '|';
	
	private final VcfRecord vcfRecord;
	private String[] tokens;

	private transient List<Allele> cachedAlleles;
	
	public VcfSample(VcfRecord vcfRecord) {
		this(vcfRecord, null);
	}
	
	public VcfSample(VcfRecord vcfRecord, String[] tokens) {
		if(vcfRecord == null) throw new IllegalArgumentException("vcfRecord is null");
		this.vcfRecord = vcfRecord;
		this.tokens = tokens;
	}
	
	/**
	 * @param idx
	 * @return data for the sample record at the given position or null if data is set to the missing value
	 */
	public String getData(int idx) {
		String data = idx == 0 || idx < tokens.length ? tokens[idx] : null; // return null for trailing values that are not specified 
		return data != null && data.equals(VcfRecord.MISSING_VALUE) ? null : data;
	}
	
	public List<Boolean> getPhasings() {
		// the first sub-field must always be the genotype if it is present
		String[] dataTypes = vcfRecord.getFormat();
		if(dataTypes.length == 0 || !dataTypes[0].equals(FIELD_GT)) return null;
		String genotype = tokens[0];
		
		// parse phasing
		List<Boolean> phasings = new ArrayList<Boolean>(1);
		final int nrChars = genotype.length();
		for (int i = 0; i < nrChars; ++i){
			switch(genotype.charAt(i)) {
				case GENOTYPE_PHASED:
					phasings.add(Boolean.TRUE);
					break;
				case GENOTYPE_UNPHASED:
					phasings.add(Boolean.FALSE);
					break;
				default:
					break;
			}
		}
		return phasings;
	}
	
	public List<Allele> getAlleles() {
		if(cachedAlleles == null) {
			// the first sub-field must always be the genotype if it is present
			String[] dataTypes = vcfRecord.getFormat();
			if(dataTypes.length == 0 || !dataTypes[0].equals(FIELD_GT)) return null;
			String genotype = tokens[0];
			
			Allele referenceAllele = vcfRecord.getReferenceAllele();
			List<Allele> alternateAlleles = vcfRecord.getAlternateAlleles();
			
			// performance optimization for the common case that a sample consists of two alleles
			cachedAlleles = new ArrayList<Allele>(2);
			final int nrGenotypeChars = genotype.length();
			for (int j = 0, start = 0; j < nrGenotypeChars; ++j)
			{
				char c = genotype.charAt(j);
				if (c == GENOTYPE_PHASED || c == GENOTYPE_UNPHASED || j == nrGenotypeChars - 1)
				{
					if(j - start == 1) {
						// performance optimization for the common case that an allele is described by one char
						char alleleChar = j == nrGenotypeChars - 1 ?  c : genotype.charAt(j - 1);
						if(alleleChar != '.') {
							int alleleIndex = Character.digit(alleleChar, 10);
							if(alleleIndex == 0)
								cachedAlleles.add(referenceAllele);
							else
								cachedAlleles.add(alternateAlleles.get(alleleIndex - 1));
						} else {
							cachedAlleles.add(Allele.ZERO);
						}
					} else {
						String alleleIndexStr = j == nrGenotypeChars - 1 ? genotype.substring(start) : genotype.substring(start, j);
						if (!alleleIndexStr.equals(".")) {
							int alleleIndex;
							try {
								alleleIndex = Integer.parseInt(alleleIndexStr);
							} catch (NumberFormatException ex){
								alleleIndex = Integer.parseInt(CharMatcher.JAVA_ISO_CONTROL.removeFrom(alleleIndexStr));
							}
							if(alleleIndex == 0)
								cachedAlleles.add(referenceAllele);
							else
								cachedAlleles.add(alternateAlleles.get(alleleIndex - 1));
						} else {
							cachedAlleles.add(Allele.ZERO);
						}
					}
					start = j + 1;
				}
			}
		}
		return cachedAlleles;
	}
	
	public VcfSample createClone() {
		return new VcfSample(vcfRecord, tokens);
	}
	
	public void reset(String[] tokens) {
		this.tokens = tokens;
		this.cachedAlleles = null;
	}
}
