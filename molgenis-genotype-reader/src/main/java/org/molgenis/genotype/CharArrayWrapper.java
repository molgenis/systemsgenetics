package org.molgenis.genotype;

import java.util.Arrays;

public final class CharArrayWrapper
{
	private final char[] data;

	public CharArrayWrapper(char[] data)
	{
		if (data == null)
		{
			throw new NullPointerException();
		}
		this.data = data;
	}

	@Override
	public boolean equals(Object other)
	{
		if (!(other instanceof CharArrayWrapper))
		{
			return false;
		}
		return Arrays.equals(data, ((CharArrayWrapper) other).data);
	}

	@Override
	public int hashCode()
	{
		return Arrays.hashCode(data);
	}
}