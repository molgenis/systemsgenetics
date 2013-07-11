package org.molgenis.genotype.vcf;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.molgenis.genotype.annotation.Annotation;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variant.ReadOnlyGeneticVariant;
import org.molgenis.genotype.variant.VariantLineMapper;
import org.molgenis.genotype.variant.sampleProvider.SampleVariantsProvider;

public class VcfVariantLineMapper implements VariantLineMapper
{
	private final List<String> colNames;
	private final List<Annotation> infoAnnotations;
	private final Map<String, String> altDescriptions;
	private final SampleVariantsProvider sampleVariantsProvider;

	public VcfVariantLineMapper(List<String> colNames, List<Annotation> infoAnnotations,
			Map<String, String> altDescriptions, SampleVariantsProvider sampleVariantsProvider)
	{
		this.colNames = colNames;
		this.infoAnnotations = infoAnnotations;
		this.altDescriptions = altDescriptions;
		this.sampleVariantsProvider = sampleVariantsProvider;
	}

	@Override
	public GeneticVariant mapLine(String line)
	{
		VcfRecord record = new VcfRecord(line, colNames);

		List<String> ids = record.getId();
		String sequenceName = record.getChrom();
		int startPos = record.getPos();
		List<String> alleles = record.getAlleles();
		String refAllele = record.getRef();

		Map<String, Object> annotationValues = getAnnotationValues(record, infoAnnotations);

		// Check if the alt alleles contain references to alt annotations
		List<String> altTypes = new ArrayList<String>();
		List<String> altDescriptions = new ArrayList<String>();

		if (alleles.size() > 0)
		{
			// First allele is ref
			for (int i = 1; i < alleles.size(); i++)
			{
				String alt = alleles.get(i);
				if ((alt != null) && alt.startsWith("<") && alt.endsWith(">"))
				{
					String altType = alt.substring(1, alt.length() - 1);
					altTypes.add(altType);
					String altDescription = this.altDescriptions.get(altType);
					if (altDescription != null)
					{
						altDescriptions.add(altDescription);
					}
				}
			}
		}

		GeneticVariant variant = ReadOnlyGeneticVariant.createVariant(ids, startPos, sequenceName, annotationValues,
				sampleVariantsProvider, alleles, refAllele);

		return variant;
	}

	private Map<String, Object> getAnnotationValues(VcfRecord record, List<Annotation> annotations)
	{
		Map<String, Object> annotationValues = new HashMap<String, Object>();

		for (Annotation annotation : annotations)
		{
			String annoId = annotation.getId();
			Object annoValue = null;

			List<String> values = record.getInfo(annotation.getId());
			if ((values != null) && !values.isEmpty())
			{
				switch (annotation.getType())
				{
					case INTEGER:
						if (annotation.isList())
						{
							List<Integer> ints = new ArrayList<Integer>();
							for (String value : values)
							{

								ints.add(Integer.valueOf(value));
							}
							annoValue = ints;
						}
						else
						{
							annoValue = Integer.valueOf(values.get(0));
						}
						break;
					case BOOLEAN:
						if (annotation.isList())
						{
							List<Boolean> bools = new ArrayList<Boolean>();
							for (String value : values)
							{
								bools.add(Boolean.parseBoolean(value));
							}
							annoValue = bools;
						}
						else
						{
							annoValue = Boolean.parseBoolean(values.get(0));
						}
						break;
					case FLOAT:
						if (annotation.isList())
						{
							List<Float> floats = new ArrayList<Float>();
							for (String value : values)
							{
								floats.add(Float.parseFloat(value));
							}
							annoValue = floats;
						}
						else
						{
							annoValue = Float.parseFloat(values.get(0));
						}
						break;
					case CHAR:
						if (annotation.isList())
						{
							List<Character> chars = new ArrayList<Character>();
							for (String value : values)
							{
								chars.add(value.charAt(0));
							}
							annoValue = chars;
						}
						else
						{
							annoValue = Character.valueOf(values.get(0).charAt(0));
						}
						break;
					default:
						if (annotation.isList())
						{
							annoValue = values;
						}
						else
						{
							annoValue = values.get(0);
						}
				}

				if (annoValue != null)
				{
					annotationValues.put(annoId, annoValue);
				}
			}
		}

		return annotationValues;
	}

}
