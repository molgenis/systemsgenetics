package org.molgenis.genotype.annotation;

import org.apache.commons.lang3.builder.ToStringBuilder;

/**
 * Generic Annotation class. Describes samples, sequences etc.
 */
public abstract class Annotation
{
	public enum Type
	{
		INTEGER, BOOLEAN, FLOAT, STRING, CHAR, UNKOWN
	}

	private String id;
	private String name;
	private String description;
	private Annotation.Type type;

	public Annotation(String id, String name, String description, Type type)
	{
		super();
		this.id = id;
		this.name = name;
		this.description = description;
		this.type = type;
	}

	public String getId()
	{
		return id;
	}

	public String getName()
	{
		return name;
	}

	public String getDescription()
	{
		return description;
	}

	public Annotation.Type getType()
	{
		return type;
	}

	public abstract boolean isList();

	@Override
	public int hashCode()
	{
		final int prime = 31;
		int result = 1;
		result = prime * result + ((id == null) ? 0 : id.hashCode());
		return result;
	}

	@Override
	public boolean equals(Object obj)
	{
		if (this == obj) return true;
		if (obj == null) return false;
		if (getClass() != obj.getClass()) return false;
		Annotation other = (Annotation) obj;
		if (id == null)
		{
			if (other.id != null) return false;
		}
		else if (!id.equals(other.id)) return false;
		return true;
	}

	@Override
	public String toString()
	{
		return ToStringBuilder.reflectionToString(this);
	}

}
