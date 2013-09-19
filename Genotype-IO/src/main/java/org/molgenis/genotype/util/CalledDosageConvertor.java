package org.molgenis.genotype.util;

import java.util.List;

import org.molgenis.genotype.Allele;
import org.molgenis.genotype.Alleles;

public class CalledDosageConvertor
{

	public static float[] convertCalledAllelesToDosage(List<Alleles> sampleAlleles, Alleles alleles, Allele refAllele)
	{
		Allele dosageRef = refAllele == null ? alleles.getAlleles().get(0) : refAllele;

		float[] dosages = new float[sampleAlleles.size()];

		for (int i = 0; i < dosages.length; ++i)
		{
			Alleles sampleVariant = sampleAlleles.get(i);
			boolean missing = false;
			float dosage = 0;

			sampleAlleles: for (Allele allele : sampleVariant)
			{
				if (allele == null || allele == Allele.ZERO)
				{
					missing = true;
					break sampleAlleles;
				}
				else if (allele == dosageRef)
				{
					++dosage;
				}
			}

			dosages[i] = missing ? -1 : dosage;
		}

		return dosages;
	}

	public static byte[] convertCalledAllelesToCalledDosage(List<Alleles> sampleAlleles, Alleles alleles,
			Allele refAllele)
	{
		Allele dosageRef = refAllele == null ? alleles.getAlleles().get(0) : refAllele;

		byte[] dosages = new byte[sampleAlleles.size()];

		for (int i = 0; i < dosages.length; ++i)
		{
			Alleles sampleVariant = sampleAlleles.get(i);
			boolean missing = false;
			byte dosage = 0;

			for (Allele allele : sampleVariant)
			{

				if (allele == null || allele == Allele.ZERO)
				{
					missing = true;
				}
				else if (allele == dosageRef)
				{
					++dosage;
				}
			}

			dosages[i] = missing ? -1 : dosage;
		}

		return dosages;
	}

	public static float[] convertCalledDosageToDosage(byte[] calledDosage)
	{
		float[] dosage = new float[calledDosage.length];

		for (int i = 0; i < calledDosage.length; ++i)
		{
			dosage[i] = calledDosage[i];
		}

		return dosage;
	}

}
