package org.molgenis.genotype.util;

import java.util.HashMap;

import org.molgenis.genotype.variant.GeneticVariant;

public class Ld
{

	private final GeneticVariant variant1;
	private final GeneticVariant variant2;
	private final double r2;
	private final double dPrime;
	private final HashMap<String, Double> haplotypesFreq;

	public Ld(GeneticVariant variant1, GeneticVariant variant2, double r2, double dPrime,
			HashMap<String, Double> haplotypesFreq)
	{
		this.variant1 = variant1;
		this.variant2 = variant2;
		this.r2 = r2;
		this.dPrime = dPrime;
		this.haplotypesFreq = haplotypesFreq;
	}

	public GeneticVariant getVariant1()
	{
		return variant1;
	}

	public GeneticVariant getVariant2()
	{
		return variant2;
	}

	public double getR2()
	{
		return r2;
	}

	public double getDPrime()
	{
		return dPrime;
	}

	public HashMap<String, Double> getHaplotypesFreq()
	{
		return haplotypesFreq;
	}

}
