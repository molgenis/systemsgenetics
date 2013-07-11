package org.molgenis.genotype.plink.datatypes;

import org.molgenis.util.tuple.KeyValueTuple;
import org.molgenis.util.tuple.Tuple;
import org.molgenis.util.tuple.WritableTuple;

/**
 * See: http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#map
 * 
 */
public class MapEntry
{
	private String chromosome;
	private String SNP;
	private double cM;
	private long bpPos;

	public MapEntry(String chromosome, String SNP, double cM, long bpPos)
	{
		this.chromosome = chromosome;
		this.SNP = SNP;
		this.cM = cM;
		this.bpPos = bpPos;
	}

	public static String[] mapHeader()
	{
		return new String[]
		{ "chr", "snp", "cm", "bp" };
	}

	public static Tuple mapToTuple(MapEntry map)
	{
		WritableTuple tuple = new KeyValueTuple();
		tuple.set("chr", map.getChromosome());
		tuple.set("snp", map.getSNP());
		tuple.set("cm", map.getcM());
		tuple.set("bp", map.getBpPos());
		return tuple;
	}

	public String getChromosome()
	{
		return chromosome;
	}

	public String getSNP()
	{
		return SNP;
	}

	public double getcM()
	{
		return cM;
	}

	public long getBpPos()
	{
		return bpPos;
	}

}
