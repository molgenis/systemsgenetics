package org.molgenis.genotype;

import org.molgenis.genotype.oxford.GenGenotypeWriter;
import org.molgenis.genotype.oxford.HapsGenotypeWriter;
import org.molgenis.genotype.plink.BedBimFamGenotypeWriter;
import org.molgenis.genotype.plink.PedMapGenotypeWriter;

public enum GenotypedDataWriterFormats
{

	PED_MAP("PED / MAP files"), 
	SHAPEIT2("Impute2 haplotypes haps / sample files"), 
	PLINK_BED("Plink BED / BIM / FAM files"),
	GEN("Oxford gen / sample files");

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
				return new HapsGenotypeWriter(genotypeData);
			case PLINK_BED:
				return new BedBimFamGenotypeWriter(genotypeData);
			case GEN:
				return new GenGenotypeWriter(genotypeData);
			default:
				throw new RuntimeException("This should not be reachable. Please contact the authors");
		}
		
	}
}
