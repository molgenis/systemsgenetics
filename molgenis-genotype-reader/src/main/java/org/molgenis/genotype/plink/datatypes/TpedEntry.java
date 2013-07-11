package org.molgenis.genotype.plink.datatypes;

import java.util.List;

import org.molgenis.genotype.Alleles;
import org.molgenis.util.tuple.KeyValueTuple;
import org.molgenis.util.tuple.Tuple;
import org.molgenis.util.tuple.WritableTuple;

public class TpedEntry extends MapEntry
{

	// list iterates individuals, so 1 list per SNP
	// NOTE: this is the inverse of PED format!
	List<Alleles> bialleles;

	public TpedEntry(String chromosome, String SNP, double cM, long bpPos, List<Alleles> bialleles)
	{
		super(chromosome, SNP, cM, bpPos);
		this.bialleles = bialleles;
	}

	public static String[] tpedHeader()
	{
		return new String[]
		{ "chr", "snp", "cm", "bp", "bial" };
	}

	public static Tuple tpedToTuple(TpedEntry tped)
	{
		WritableTuple tuple = new KeyValueTuple();
		tuple.set("chr", tped.getChromosome());
		tuple.set("snp", tped.getSNP());
		tuple.set("cm", tped.getcM());
		tuple.set("bp", tped.getBpPos());
		tuple.set("bial", tped.getBialleles());
		return tuple;
	}

	public List<Alleles> getBialleles()
	{
		return bialleles;
	}

}
