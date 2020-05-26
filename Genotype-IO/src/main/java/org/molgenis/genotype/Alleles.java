package org.molgenis.genotype;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;

public final class Alleles implements Iterable<Allele>, Comparable<Alleles> {

	private static final Map<List<Allele>, Alleles> pool;
	private final List<Allele> alleles;
	private final boolean snp;
	private Alleles complement;
	private final boolean isAtOrGcSnp;
	private final List<String> allelesAsString;
	private final char[] allelesAsChar;
	private final int hashCode;
	public static final Alleles BI_ALLELIC_MISSING;

	static {
		pool = new HashMap<List<Allele>, Alleles>();
		BI_ALLELIC_MISSING = createAlleles(Allele.ZERO, Allele.ZERO);
	}

	private Alleles(List<Allele> alleles) {
		this.alleles = Collections.unmodifiableList(alleles);
		hashCode = alleles.hashCode();// result;

		boolean isSnp = true;
		ArrayList<String> allelesAsStringBuilder = new ArrayList<String>(alleles.size());
		for (Allele allele : alleles) {
			if (!allele.isSnpAllele()) {
				isSnp = false;
			}
			allelesAsStringBuilder.add(allele.getAlleleAsString());
		}
		this.allelesAsString = Collections.unmodifiableList(allelesAsStringBuilder);

		this.snp = isSnp;
		if (snp) {
			allelesAsChar = new char[alleles.size()];
			int i = 0;
			for (Allele allele : alleles) {
				allelesAsChar[i] = allele.getAlleleAsSnp();
				++i;
			}
		} else {
			allelesAsChar = null;
		}

		this.isAtOrGcSnp = areAlleleCharsAtOrGc(alleles);

	}

	private static boolean areAlleleCharsAtOrGc(List<Allele> alleles) {

		if (alleles.isEmpty()) {
			return false;
		}

		boolean onlyAt = true;
		boolean onlyGc = true;
		for (Allele allele : alleles) {
			if (!allele.isSnpAllele()) {
				return false;
			}
			if (allele == Allele.A || allele == Allele.T) {
				onlyGc = false;
			}
			if (allele == Allele.C || allele == Allele.G) {
				onlyAt = false;
			}
		}
		return onlyAt || onlyGc;

	}

	public static Alleles createAlleles(List<Allele> alleleList) {
		Alleles alleles = pool.get(alleleList);
		if (alleles == null) {
			alleles = new Alleles(alleleList);
			pool.put(alleleList, alleles);
			alleles.addComplement();
		}
		return alleles;
	}

	public static Alleles createAlleles(Allele... allele) {
		return createAlleles(Arrays.asList(allele));
	}

	public static Alleles createBasedOnString(List<String> stringAlleles) {
		ArrayList<Allele> alleles = new ArrayList<Allele>(stringAlleles.size());

		for (String stringAllele : stringAlleles) {
			alleles.add(Allele.create(stringAllele));
		}

		return createAlleles(alleles);
	}

	public static Alleles createBasedOnString(String allele1, String allele2) {

		return createAlleles(Allele.create(allele1), Allele.create(allele2));

	}

	public static Alleles createBasedOnChars(char allele1, char allele2) {

		return createAlleles(Allele.create(allele1), Allele.create(allele2));

	}

	public static Alleles createBasedOnChars(char[] charAlleles) {
		ArrayList<Allele> alleles = new ArrayList<Allele>(charAlleles.length);
		for (char charAllele : charAlleles) {
			alleles.add(Allele.create(charAllele));
		}
		return createAlleles(alleles);
	}

	/**
	 * Add complement. Not done in constructor to prevent infinite loop. Pool
	 * must be up to date before this is called
	 */
	private void addComplement() {
		if (snp) {
			ArrayList<Allele> complementAlleles = new ArrayList<Allele>(alleles.size());
			for (Allele allele : alleles) {
				complementAlleles.add(allele.getComplement());
			}
			this.complement = Alleles.createAlleles(complementAlleles);
		} else {
			this.complement = null;
		}

	}

	/**
	 * List of the possible alleles, can contain null if not known!!!!!
	 *
	 * @return
	 */
	public List<Allele> getAlleles() {
		return alleles;
	}

	public List<String> getAllelesAsString() {
		return allelesAsString;
	}

	public int getAlleleCount() {
		return alleles.size();
	}

	public boolean isSnp() {
		return snp;
	}

	public char[] getAllelesAsChars() {
		if (!isSnp()) {
			throw new RuntimeException("Not a snp");
		}

		return allelesAsChar.clone();
	}

	@Override
	public String toString() {

		if (alleles.isEmpty()) {
			return "";
		} else {
			StringBuilder s = new StringBuilder(3);

			s.append(allelesAsString.get(0));
			for (int i = 1; i < alleles.size(); ++i) {
				s.append('\\');
				s.append(allelesAsString.get(i));
			}

			return s.toString();
		}


	}

	/**
	 * Returns the complements of this variant alleles. Currently only works for
	 * SNPs
	 *
	 * @return complement of current variant alleles
	 */
	public Alleles getComplement() {

		if (!isSnp()) {
			throw new RuntimeException("Complement currenlty only supported for SNPs");
		}

		return complement;
	}

	/**
	 * Assess if two variantAllele instances have same alleles regardless of
	 * order. Only true if also identical number of alleles
	 *
	 * @param other
	 * @return
	 */
	public boolean sameAlleles(Alleles other) {
		if (this == other) {
			return true;
		}
		if (this.alleles.size() != other.alleles.size()) {
			return false;
		}
		return this.alleles.containsAll(other.alleles) && other.alleles.containsAll(this.alleles);
	}

	public boolean isAtOrGcSnp() {
		return isAtOrGcSnp;
	}

	@Override
	public Iterator<Allele> iterator() {
		return alleles.iterator();
	}

	public Allele get(int alleleIndex) {
		return alleles.get(alleleIndex);
	}

	public boolean contains(Allele queryAllele) {
		return (alleles.contains(queryAllele));
	}

	public boolean containsAll(Alleles queryAlleles) {
		for (Allele queryAllele : queryAlleles) {
			if (!contains(queryAllele)) {
				return false;
			}
		}
		return true;
	}

	@Override
	public int compareTo(Alleles other) {
		if (this == other) {
			return 0;
		}

		Iterator<Allele> thisAlleleIterator = this.alleles.iterator();
		Iterator<Allele> otherAlleleIterator = other.alleles.iterator();

		while (thisAlleleIterator.hasNext() && otherAlleleIterator.hasNext()) {
			Allele thisCurrentAllele = thisAlleleIterator.next();
			Allele otherCurrentAllele = otherAlleleIterator.next();

			if (thisCurrentAllele != otherCurrentAllele) {
				return thisCurrentAllele.compareTo(otherCurrentAllele);
			}
		}

		if (thisAlleleIterator.hasNext()) {
			return 1;
		}

		if (otherAlleleIterator.hasNext()) {
			return -1;
		}

		return 0;

	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#hashCode()
	 */
	@Override
	public int hashCode() {
		return hashCode;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#equals(java.lang.Object)
	 */
	@Override
	public boolean equals(Object obj) {
		if (this == obj) {
			return true;
		}
		if (obj == null) {
			return false;
		}
		if (getClass() != obj.getClass()) {
			return false;
		}
		Alleles other = (Alleles) obj;
		if (alleles == null) {
			if (other.alleles != null) {
				return false;
			}
		} else if (!alleles.equals(other.alleles)) {
			return false;
		}
		return true;
	}

	public Alleles createCopyWithoutDuplicates() {

		LinkedHashSet<Allele> uniqueAlleles = new LinkedHashSet<Allele>(alleles.size());

		for (Allele allele : alleles) {
			if (!uniqueAlleles.contains(allele)) {
				uniqueAlleles.add(allele);
			}
		}

		return Alleles.createAlleles(new ArrayList<Allele>(uniqueAlleles));

	}
}
