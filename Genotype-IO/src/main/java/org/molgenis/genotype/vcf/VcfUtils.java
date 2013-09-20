package org.molgenis.genotype.vcf;

public class VcfUtils
{
	private static final String NULL_VALUE = ".";

	private VcfUtils()
	{
		// Do not instanciate
	}

	/**
	 * Filter null values ('.')
	 */
	public static String checkNullValue(String s)
	{
		if (s == null)
		{
			return null;
		}

		if (s.equals(NULL_VALUE))
		{
			return null;
		}

		return s;
	}

}
