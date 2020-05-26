package org.molgenis.genotype.vcf;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.lang3.StringUtils;
import org.molgenis.genotype.Allele;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.GenotypeDataException;
import org.molgenis.genotype.variant.GenotypeRecord;
import org.molgenis.vcf.VcfRecord;
import org.molgenis.vcf.VcfSample;
import org.molgenis.vcf.meta.VcfMeta;
import org.molgenis.vcf.meta.VcfMetaFormat;

public class VcfGenotypeRecord implements GenotypeRecord
{
	private final VcfMeta vcfMeta;
	private final VcfRecord vcfRecord;
	private final VcfSample vcfSample;

	public VcfGenotypeRecord(VcfMeta vcfMeta, VcfRecord vcfRecord, VcfSample vcfSample) {
		if(vcfMeta == null) throw new IllegalArgumentException("vcfMeta is null");
		if(vcfRecord == null) throw new IllegalArgumentException("vcfRecord is null");
		if(vcfSample == null) throw new IllegalArgumentException("vcfSample is null");
		this.vcfMeta = vcfMeta;
		this.vcfRecord = vcfRecord.createClone();
		this.vcfSample = vcfSample.createClone();
	}
	
	@Override
	public Object getGenotypeRecordData(String recordId)
	{
		int idx = vcfRecord.getFormatIndex(recordId);
		if(idx == -1) throw new GenotypeDataException("unknown record id [" + recordId + "]");
		
		String data = vcfSample.getData(idx);
		VcfMetaFormat format = vcfMeta.getFormatMeta(recordId);
		if(format == null) throw new GenotypeDataException("missing vcf format data for record id [" + recordId + "]");
		
		Object value;
		switch(format.getType()) {
			case CHARACTER:
				return Character.valueOf(data.charAt(0));
			case FLOAT:
				if (isListValue(format.getNumber()))
				{
					String[] tokens = StringUtils.split(data, ',');
					List<Double> doubles = new ArrayList<Double>(tokens.length);
					for(String token : tokens)
						doubles.add(Double.valueOf(token));
					value = doubles;
				}
				else
				{
					value = Double.valueOf(data);
				}
				break;
			case INTEGER:
				if (isListValue(format.getNumber()))
				{
					String[] tokens = StringUtils.split(data, ',');
					List<Integer> integers = new ArrayList<Integer>(tokens.length);
					for(String token : tokens)
						integers.add(Integer.valueOf(token));
					value = integers;
				}
				else
				{
					value = Integer.valueOf(data);
				}
				break;
			case STRING:
				if (isListValue(format.getNumber()))
				{
					value = Arrays.asList(StringUtils.split(data, ','));
				}
				else
				{
					value = data;
				}
				break;
			default:
				throw new IllegalArgumentException("invalid vcf format type [" + format.getType() + "]");
		}
		return value;
	}

	@Override
	public Alleles getSampleAlleles()
	{
		List<Allele> alleles = vcfSample.getAlleles();
		if(vcfSample.getAlleles() != null){
			return Alleles.createAlleles(alleles);
		} else {
			return null;
		}
		
	}
	
	private boolean isListValue(String number) {
		// A: one value per alternate allele
		// R: one value for each possible allele (including the reference)
		// G: one value for each possible genotype
		// .: number of possible values varies, is unknown, or is unbounded
		return number.equals("A") || number.equals("R") || number.equals("G") || number.equals(".") || Integer.valueOf(number) > 1;
	}

	@Override
	public float[] getSampleProbs() {
		
		List<Double> probs = (List<Double>) getGenotypeRecordData("GP");
		
		if(probs != null){
			float[] probsArray = new float[probs.size()];
			
			for(int i = 0 ; i < probs.size() ; ++i){
				probsArray[i] = probs.get(i).floatValue();
			}
			
			return probsArray;
			
		} else {
			return null;
		}
		
		
	}

	@Override
	public float getSampleDosage() {
		//TODO
		throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
	}

	@Override
	public boolean containsGenotypeRecord(String recordId) {
		return vcfRecord.getFormatIndex(recordId) != -1;
	}
}
