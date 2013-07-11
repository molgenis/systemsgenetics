package org.molgenis.genotype.vcf;

import org.molgenis.util.tuple.Tuple;

/**
 * Info metadata atttached to VcfReader
 * 
 * Example: #INFO=<ID=NS,Number=1,Type=Integer,Description=
 * "Number of Samples With Data">
 */
public class VcfInfo
{
	public enum Type
	{
		INTEGER, FLOAT, FLAG, CHARACTER, STRING, UNKNOWN;

		public static Type getValue(String value)
		{
			if ("Integer".equals(value)) return INTEGER;
			if ("Float".equals(value)) return FLOAT;
			if ("Flag".equals(value)) return FLAG;
			if ("Character".equals(value)) return CHARACTER;
			if ("String".equals(value)) return STRING;
			return UNKNOWN;
		}
	};

	private final String id;
	private final Type type;
	private final String number;
	private final String description;

	public VcfInfo(Tuple settings)
	{
		id = settings.getString("ID");
		type = Type.getValue(settings.getString("Type"));
		number = settings.getString("Number");
		description = settings.getString("Description");
	}

	public VcfInfo(String id, Type type, String number, String description)
	{
		this.id = id;
		this.type = type;
		this.number = number;
		this.description = description;
	}

	public String getId()
	{
		return id;
	}

	public Type getType()
	{
		return type;
	}

	public String getNumber()
	{
		return number;
	}

	public String getDescription()
	{
		return description;
	}

	@Override
	public String toString()
	{
		return String.format("VcfInfo(ID=%s,Number=%s,Type=%s,Description=\"%s\")", getId(), getNumber(), getType(),
				getDescription());
	}
}
