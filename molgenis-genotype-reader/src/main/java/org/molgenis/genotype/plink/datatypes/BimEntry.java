package org.molgenis.genotype.plink.datatypes;

import org.molgenis.genotype.Alleles;
import org.molgenis.util.tuple.KeyValueTuple;
import org.molgenis.util.tuple.Tuple;
import org.molgenis.util.tuple.WritableTuple;

public class BimEntry extends MapEntry
{

	private Alleles biallele;

	public BimEntry(String chromosome, String SNP, double cM, long bpPos, Alleles biallele)
	{
		super(chromosome, SNP, cM, bpPos);
		this.biallele = biallele;
	}

	public Alleles getBiallele()
	{
		return biallele;
	}

	public static String[] bimHeader()
	{
		return new String[]
		{ "chr", "snp", "cm", "bp", "al1", "al2" };
	}

	public static Tuple bimToTuple(BimEntry bim)
	{
		WritableTuple tuple = new KeyValueTuple();
		tuple.set("chr", bim.getChromosome());
		tuple.set("snp", bim.getSNP());
		tuple.set("cm", bim.getcM());
		tuple.set("bp", bim.getBpPos());
		tuple.set("al1", bim.getBiallele().get(0));
		tuple.set("al2", bim.getBiallele().get(1));
		return tuple;
	}
}
