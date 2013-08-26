package org.molgenis.genotype;

import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

import org.molgenis.genotype.variant.GeneticVariant;

public class SimpleSequence implements Sequence
{
	private static final List<String> CHROMOSOMES = Arrays.asList("1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
			"11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "x", "y", "xy", "mt", "23", "24",
			"25", "26");
	private static final List<String> AUTOSOMES = Arrays.asList("1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
			"11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22");

	private String name;
	private Integer length;
	private RandomAccessGenotypeData genotypeData;

	public SimpleSequence(String name, Integer length, RandomAccessGenotypeData genotypeData)
	{
		this.name = name.intern();
		this.length = length;
		this.genotypeData = genotypeData;
	}

	@Override
	public String getName()
	{
		return name;
	}

	@Override
	public Integer getLength()
	{
		return length;
	}

	@Override
	public boolean isChromosome()
	{
		return CHROMOSOMES.contains(name.toLowerCase());
	}

	@Override
	public boolean isAutosome()
	{
		return AUTOSOMES.contains(name.toLowerCase());
	}

	@Override
	public int hashCode()
	{
		final int prime = 31;
		int result = 1;
		result = prime * result + ((genotypeData == null) ? 0 : genotypeData.hashCode());
		result = prime * result + ((name == null) ? 0 : name.hashCode());
		return result;
	}

	@Override
	public boolean equals(Object obj)
	{
		if (this == obj) return true;
		if (obj == null) return false;
		if (getClass() != obj.getClass()) return false;
		SimpleSequence other = (SimpleSequence) obj;
		if (genotypeData == null)
		{
			if (other.genotypeData != null) return false;
		}
		else if (!genotypeData.equals(other.genotypeData)) return false;
		if (name == null)
		{
			if (other.name != null) return false;
		}
		else if (!name.equals(other.name)) return false;
		return true;
	}

}
