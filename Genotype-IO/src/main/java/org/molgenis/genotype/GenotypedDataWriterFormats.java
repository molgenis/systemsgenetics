package org.molgenis.genotype;

import org.molgenis.genotype.bgen.BgenGenotypeData;
import org.molgenis.genotype.bgen.BgenGenotypeWriter;
import org.molgenis.genotype.oxford.GenGenotypeWriter;
import org.molgenis.genotype.oxford.HapsGenotypeWriter;
import org.molgenis.genotype.plink.BedBimFamGenotypeWriter;
import org.molgenis.genotype.plink.PedMapGenotypeWriter;
import org.molgenis.genotype.table.TableGenotypeWriter;
import org.molgenis.genotype.trityper.TriTyperGenotypeWriter;

public enum GenotypedDataWriterFormats
{

	PED_MAP("PED / MAP files"), 
	SHAPEIT2("Impute2 haplotypes haps / sample files"), 
	PLINK_BED("Plink BED / BIM / FAM files"),
	GEN("Oxford GEN / SAMPLE files"),
	BGEN("Oxford Binary GEN / SAMPLE files"),
	TRITYPER("Trityper folder"),
	TABLE("Simple tab separated files with dosage and genotypes");

	private final String name;

	GenotypedDataWriterFormats(String name)
	{
		this.name = name;
	}

	public String getName()
	{
		return name;
	}

	public GenotypeWriter createGenotypeWriter(GenotypeData genotypeData, Integer bgenBitRepresentation)
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
			case BGEN:
				BgenGenotypeWriter bgenGenotypeWriter = new BgenGenotypeWriter(genotypeData);
				bgenGenotypeWriter.setProbabilityPrecisionInBits(bgenBitRepresentation);
				return bgenGenotypeWriter;
			case TRITYPER:
				return new TriTyperGenotypeWriter(genotypeData);
			case TABLE:
				return new TableGenotypeWriter(genotypeData);
			default:
				throw new RuntimeException("This should not be reachable. Please contact the authors");
		}
		
	}
}
