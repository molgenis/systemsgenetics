package org.molgenis.genotype.plink.datatypes;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import org.molgenis.genotype.Alleles;

public class PedEntry implements Iterable<Alleles> {

	// list iterates SNP's, so 1 list per individual
	private final Iterator<Alleles> bialleles;
	private final String family;
	private final String individual;
	private final String father;
	private final String mother;
	private final byte sex;
	private final double phenotype;
	
	public PedEntry(String family, String individual, String father, String mother, byte sex, double phenotype,
			Iterator<Alleles> bialleles) {
		this.family = family;
		this.individual = individual;
		this.father = father;
		this.mother = mother;
		this.sex = sex;
		this.phenotype = phenotype;
		this.bialleles = bialleles;
	}

	@Override
	public Iterator<Alleles> iterator() {
		return bialleles;
	}

	public List<Alleles> getBialleles() {
		List<Alleles> bialleleList = new ArrayList<Alleles>();
		Iterator<Alleles> it = iterator();
		while (it.hasNext()) {
			bialleleList.add(it.next());
		}

		return bialleleList;
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
