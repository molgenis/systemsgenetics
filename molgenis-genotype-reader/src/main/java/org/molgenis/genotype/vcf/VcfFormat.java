package org.molgenis.genotype.vcf;

import org.molgenis.util.tuple.Tuple;

public class VcfFormat
{
	public enum Type
	{
		INTEGER, FLOAT, CHARACTER, STRING, UNKNOWN;

		public static Type getValue(String value)
		{
			if ("Integer".equals(value)) return INTEGER;
			if ("Float".equals(value)) return FLOAT;
			if ("Character".equals(value)) return CHARACTER;
			if ("String".equals(value)) return STRING;
			return UNKNOWN;
		}
	};

	public VcfFormat(Tuple settings)
	{
		id = settings.getString("ID");
		number = settings.getString("Number");
		type = Type.getValue(settings.getString("Type"));
		description = settings.getString("Description");
	}

	public VcfFormat(String id, Type type, String number, String description)
	{
		this.id = id;
		this.type = type;
		this.number = number;
		this.description = description;
	}

	private String id;
	private Type type;
	private String number;
	private String description;

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
		return String.format("VcfFormat(ID=%s,Number=%s,Type=%s,Description=\"%s\")", getId(), getNumber(), getType(),
				getDescription());
	}
}
