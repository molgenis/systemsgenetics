package org.molgenis.genotype.variant;

import org.molgenis.genotype.Alleles;

public class IllegalReferenceAlleleException extends Exception
{

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	@SuppressWarnings("unused")
	private IllegalReferenceAlleleException()
	{
	}

	public IllegalReferenceAlleleException(Alleles variantAlleles)
	{
		super("Illigal reference allele. Only 1 allele can be ref but found " + variantAlleles.getAlleleCount()
				+ " alleles");
	}

}
