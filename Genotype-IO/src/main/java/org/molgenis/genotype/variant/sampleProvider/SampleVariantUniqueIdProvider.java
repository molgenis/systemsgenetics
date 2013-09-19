package org.molgenis.genotype.variant.sampleProvider;

public class SampleVariantUniqueIdProvider
{

	private static int counter = 0;

	public static synchronized int getNextUniqueId()
	{
		int nextId = counter;
		counter++;
		return nextId;
	}
}
