package org.molgenis.genotype.annotation;

public enum SexAnnotation
{

	MALE((byte) 1), FEMALE((byte) 2), UNKNOWN((byte) 0);

	private final byte plinkSex;

	private SexAnnotation(byte plinkSex)
	{
		this.plinkSex = plinkSex;
	}

	/**
	 * @return the plinkSex
	 */
	public byte getPlinkSex()
	{
		return plinkSex;
	}

	public static SexAnnotation getSexAnnotationForPlink(byte plinkSex)
	{
		switch (plinkSex)
		{
			case 1:
				return SexAnnotation.MALE;
			case 2:
				return SexAnnotation.FEMALE;
			default:
				return SexAnnotation.UNKNOWN;
		}
	}

}
