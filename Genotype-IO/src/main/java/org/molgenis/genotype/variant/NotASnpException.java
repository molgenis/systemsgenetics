package org.molgenis.genotype.variant;

import org.molgenis.genotype.GenotypeDataException;

public class NotASnpException extends GenotypeDataException
{

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	private final GeneticVariant geneticVariant;

	@SuppressWarnings("unused")
	private NotASnpException()
	{
		geneticVariant = null;
	}

	public NotASnpException(GeneticVariant geneticVariant)
	{
		super("This variant is not a SNP " + geneticVariant.getSequenceName() + ":" + geneticVariant.getStartPos()
				+ " " + geneticVariant.getPrimaryVariantId());
		this.geneticVariant = geneticVariant;
	}

	/**
	 * @return the geneticVariant
	 */
	public GeneticVariant getGeneticVariant()
	{
		return geneticVariant;
	}

}
