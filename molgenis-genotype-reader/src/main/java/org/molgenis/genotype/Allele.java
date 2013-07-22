package org.molgenis.genotype;

import java.util.HashMap;
import java.util.Map;

import org.molgenis.genotype.util.Utils;

public class Allele implements Comparable<Allele>
{

	private static Map<String, Allele> pool = new HashMap<String, Allele>();
	private static Map<Character, Allele> snpPool = new HashMap<Character, Allele>();

	public final static Allele A = create('A');
	public final static Allele C = create('C');
	public final static Allele G = create('G');
	public final static Allele T = create('T');
	public final static Allele ZERO = create('0');

	private final String allele;
	private final char snpAllele;
	private Allele complement;
	private final int hashCode;

	private Allele(String allele)
	{
		if (allele.length() == 1)
		{
			if (allele.charAt(0) == 'I' || allele.charAt(0) == 'D')
			{
				// Some times these are found in plink files. These are not SNPs
				this.snpAllele = (char) -1;
			}
			else
			{
				this.snpAllele = allele.charAt(0);
			}

		}
		else
		{
			this.snpAllele = (char) -1;
		}
		this.allele = allele;
		this.hashCode = allele.hashCode();
	}

	private Allele(char allele)
	{
		this(String.valueOf(allele));
	}

	public boolean isSnpAllele()
	{
		return (byte) snpAllele != -1;
	}

	/**
	 * @return the allele
	 */
	public String getAlleleAsString()
	{
		return allele;
	}

	/**
	 * @return the snpAllele
	 */
	public char getAlleleAsSnp()
	{
		return snpAllele;
	}

	private void addComplement(Allele complement)
	{
		this.complement = complement;
	}

	public Allele getComplement()
	{
		if (isSnpAllele())
		{
			return complement;
		}
		else
		{
			throw new RuntimeException("Complement currenlty only supported for SNPs");
		}
	}

	public static Allele create(String alleleString)
	{
		
		if(alleleString == null){
			return ZERO;
		}
		
		if(alleleString.isEmpty()){
			return ZERO;
		}

		if (pool.containsKey(alleleString))
		{
			return pool.get(alleleString);
		}
		else
		{

			Allele newAllele = new Allele(alleleString);
			pool.put(alleleString, newAllele);
			if (newAllele.isSnpAllele())
			{
				snpPool.put(newAllele.getAlleleAsSnp(), newAllele);
				newAllele.addComplement(Allele.create(Utils.getComplementNucleotide(newAllele.getAlleleAsSnp())));
			}
			return newAllele;

		}
	}

	public static Allele create(char alleleChar)
	{

		if(alleleChar == '\0'){
			return Allele.ZERO;
		}
		
		if (snpPool.containsKey(alleleChar))
		{

			return snpPool.get(alleleChar);
		}
		else
		{
			Allele newAllele = new Allele(alleleChar);
			snpPool.put(alleleChar, newAllele);
			pool.put(newAllele.getAlleleAsString(), newAllele);
			newAllele.addComplement(Allele.create(Utils.getComplementNucleotide(alleleChar)));
			return newAllele;
		}

	}

	@Override
	public int hashCode()
	{
		return hashCode;
	}

	@Override
	public boolean equals(Object obj)
	{
		if (this == obj) return true;
		if (obj == null) return false;
		if (getClass() != obj.getClass()) return false;
		Allele other = (Allele) obj;
		if (allele == null)
		{
			if (other.allele != null) return false;
		}
		else if (!allele.equals(other.allele)) return false;
		return true;
	}

	@Override
	public String toString()
	{
		return this.getAlleleAsString();
	}

	@Override
	public int compareTo(Allele other)
	{
		return this.allele.compareTo(other.allele);
	}
}
