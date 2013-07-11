package org.molgenis.genotype;


public enum GenotypedDataWriterFormats
{

	PED_MAP("PED / MAP file");

	private final String name;

	GenotypedDataWriterFormats(String name)
	{
		this.name = name;
	}

	public String getName()
	{
		return name;
	}

}
