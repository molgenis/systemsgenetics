package org.molgenis.genotype;

public class GenotypeDataException extends RuntimeException
{

	private static final long serialVersionUID = -4244460220718452547L;

	public GenotypeDataException()
	{
	}

	public GenotypeDataException(String message)
	{
		super(message);
	}

	public GenotypeDataException(Throwable cause)
	{
		super(cause);
	}

	public GenotypeDataException(String message, Throwable cause)
	{
		super(message, cause);
	}

}
