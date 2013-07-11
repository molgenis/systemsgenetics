package org.molgenis.genotype.impute2;

import java.util.List;

/**
 * Represents one line in a haps file
 * 
 * @author erwin
 * 
 */
public class HapsEntry
{
	private int chromosomeNumber; // 23:x,24:y
	private String snpId;
	private int position;
	private String firstAllele;
	private String secondAllele;
	private List<String[]> sampleAlelles;// 0,1 or ?

	public HapsEntry(int chromosomeNumber, String snpId, int position, String firstAllele, String secondAllele,
			List<String[]> sampleAlelles)
	{
		this.chromosomeNumber = chromosomeNumber;
		this.snpId = snpId;
		this.position = position;
		this.firstAllele = firstAllele;
		this.secondAllele = secondAllele;
		this.sampleAlelles = sampleAlelles;
	}

	public int getChromosomeNumber()
	{
		return chromosomeNumber;
	}

	public String getSnpId()
	{
		return snpId;
	}

	public int getPosition()
	{
		return position;
	}

	public String getFirstAllele()
	{
		return firstAllele;
	}

	public String getSecondAllele()
	{
		return secondAllele;
	}

	public List<String[]> getSampleAlleles()
	{
		return sampleAlelles;
	}

}
