package org.molgenis.genotype.util;

import org.molgenis.genotype.Allele;

public class MafResult
{

	private final Allele minorAllele;
	private final double freq;

	public MafResult(Allele minorAllele, double freq)
	{
		super();
		this.minorAllele = minorAllele;
		this.freq = freq;
	}

	/**
	 * @return the minorAllele
	 */
	public Allele getMinorAllele()
	{
		return minorAllele;
	}

	/**
	 * @return the freq
	 */
	public double getFreq()
	{
		return freq;
	}

}
