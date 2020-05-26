package org.molgenis.genotype.annotation;

import org.apache.commons.lang3.math.NumberUtils;
import org.molgenis.vcf.meta.VcfMetaInfo;

public class VcfAnnotation extends Annotation
{
	public static final int NUMBER_UNKNOWN = -1;

	private final Integer number;
	boolean unbounded;
	private final boolean perAltAllele;
	private final boolean perGenotype;

	public static VcfAnnotation fromVcfInfo(VcfMetaInfo info)
	{
		Annotation.Type type = VcfAnnotation.toAnnotationType(info.getType());

		Integer number = null;
		boolean unbounded = false;
		boolean perAltAllele = false;
		boolean perGenotype = false;

		// Number can be A -> field has one value per alternate allele
		// Or G -> the field has one value for each possible genotype
		// Or . -> is unknown, or is unbounded
		// Or an Integer
		if (info.getNumber() != null)
		{
			if (info.getNumber().equalsIgnoreCase("A"))
			{
				perAltAllele = true;
			}
			else if (info.getNumber().equalsIgnoreCase("G"))
			{
				perGenotype = true;
			}
			else if (info.getNumber().equals("."))
			{
				unbounded = true;
			}
			else if (NumberUtils.isDigits(info.getNumber()))
			{
				number = Integer.parseInt(info.getNumber());
			}
		}

		return new VcfAnnotation(info.getId(), info.getDescription(), type, number, unbounded, perAltAllele,
				perGenotype);
	}

	public VcfAnnotation(String id, String description, Annotation.Type type, Integer number, boolean unbounded,
			boolean perAltAllele, boolean perGenotype)
	{
		super(id, id, description, type);
		this.number = number;
		this.perAltAllele = perAltAllele;
		this.perGenotype = perGenotype;
		this.unbounded = unbounded;
	}

	public Integer getNumber()
	{
		return number;
	}

	public boolean isPerAltAllele()
	{
		return perAltAllele;
	}

	public boolean isPerGenotype()
	{
		return perGenotype;
	}

	public boolean isUnbounded()
	{
		return unbounded;
	}

	@Override
	public boolean isList()
	{
		return (number != null) && (number > 1) || perAltAllele || perGenotype || unbounded;
	}

	private static Annotation.Type toAnnotationType(VcfMetaInfo.Type infoType)
	{
		if (infoType == VcfMetaInfo.Type.CHARACTER) return Type.CHAR;
		if (infoType == VcfMetaInfo.Type.STRING) return Type.STRING;
		if (infoType == VcfMetaInfo.Type.FLAG) return Type.BOOLEAN;
		if (infoType == VcfMetaInfo.Type.FLOAT) return Type.FLOAT;
		if (infoType == VcfMetaInfo.Type.INTEGER) return Type.INTEGER;
		return Type.UNKOWN;
	}
}
