package org.molgenis.genotype.vcf;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import org.molgenis.util.tuple.Tuple;

/**
 * Sample meta header
 * 
 * ##SAMPLE=<ID=S_ID,Genomes=G1_ID;G2_ID; ...;GK_ID,Mixture=N1;N2;
 * ...;NK,Description=S1;S2; ...; SK >
 * 
 * @author erwin
 * 
 */
public class VcfSample
{
	private final Tuple tuple;

	public VcfSample(Tuple tuple)
	{
		this.tuple = tuple;
	}

	public String getId()
	{
		return tuple.getString("ID");
	}

	public List<String> getGenomes()
	{
		String genomesValue = tuple.getString("Genomes");
		if (genomesValue == null)
		{
			return Collections.emptyList();
		}

		return Collections.unmodifiableList(Arrays.asList(genomesValue.split(";")));
	}

	public List<String> getMixture()
	{
		String mixtureValue = tuple.getString("Mixture");
		if (mixtureValue == null)
		{
			return Collections.emptyList();
		}

		return Collections.unmodifiableList(Arrays.asList(mixtureValue.split(";")));
	}

	public List<String> getDescription()
	{
		String descriptionValue = tuple.getString("Description");
		if (descriptionValue == null)
		{
			return Collections.emptyList();
		}

		return Collections.unmodifiableList(Arrays.asList(descriptionValue.split(";")));
	}
}
