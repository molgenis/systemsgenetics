package org.molgenis.genotype.plink.datatypes;

public class FamEntry {
	// see: http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml

	String family;
	String individual;
	String father;
	String mother;
	byte sex;
	double phenotype;

	public FamEntry(String family, String individual, String father, String mother, byte sex, double phenotype) {
		super();
		this.family = family;
		this.individual = individual;
		this.father = father;
		this.mother = mother;
		this.sex = sex;
		this.phenotype = phenotype;
	}

	public static String[] famHeader() {
		return new String[]{"fam", "ind", "fa", "mo", "sex", "phen"};
	}

	public String getFamily() {
		return family;
	}

	public String getIndividual() {
		return individual;
	}

	public String getFather() {
		return father;
	}

	public String getMother() {
		return mother;
	}

	public byte getSex() {
		return sex;
	}

	public double getPhenotype() {
		return phenotype;
	}
}
