package org.molgenis.genotype;

import org.molgenis.genotype.impute2.Impute2GenotypeWriter;
import org.molgenis.genotype.plink.PedMapGenotypeWriter;

public enum GenotypedDataWriterFormats
{

	PED_MAP("ped / map file"), SHAPEIT2("Impute2 haplotypes haps / sample files");

	private final String name;

	GenotypedDataWriterFormats(String name)
	{
		this.name = name;
	}

	public String getName()
	{
		return name;
	}

	public GenotypeWriter createGenotypeWriter(GenotypeData genotypeData)
	{

		switch (this)
		{
			case PED_MAP:
				return new PedMapGenotypeWriter(genotypeData);
			case SHAPEIT2:
				return new Impute2GenotypeWriter(genotypeData);
			default:
				throw new RuntimeException("This should not be reachable. Please contact the authors");
		}
	}
}
