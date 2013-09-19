package org.molgenis.genotype.vcf;

import org.molgenis.util.tuple.Tuple;

public class VcfContig
{
	private String id;
	private Integer length;

	public VcfContig(Tuple t)
	{
		id = t.getString("ID");
		length = t.getInt("length");
	}

	public VcfContig(String id, Integer length)
	{
		this.id = id;
		this.length = length;
	}

	public String getId()
	{
		return id;
	}

	public Integer getLength()
	{
		return length;
	}

}
