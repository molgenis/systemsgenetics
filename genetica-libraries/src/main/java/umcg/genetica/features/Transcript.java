/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.features;

import umcg.genetica.enums.Chromosome;
import umcg.genetica.enums.Strand;

import java.util.ArrayList;
import java.util.HashMap;

/**
 * @author Harm-Jan
 */
public class Transcript extends Feature {

	private final Gene gene;
	private ArrayList<Exon> exons;
	private HashMap<Exon, Integer> exonRanks;
	private ArrayList<UTR> UTRs;

	public Transcript(String name, Chromosome chromosome, Strand strand, Gene gene) {
		this.name = name;
		this.chromosome = chromosome;
		this.strand = strand;
		this.gene = gene;
	}

	public Transcript(String name, Chromosome chromosome, Strand strand, Gene gene, int start, int stop) {
		this.name = name;
		this.chromosome = chromosome;
		this.strand = strand;
		this.gene = gene;
		this.start = start;
		this.stop = stop;
	}

	public Gene getGene() {
		return gene;
	}

	public ArrayList<Exon> getExons() {
		return exons;
	}

	public void addExon(Exon e) {
		if (exons == null) {
			exons = new ArrayList<Exon>();
		}
		int estart = e.getStart();
		int estop = e.getStop();
		if (estart < start) {
			start = estart;
		}
		if (estop > stop) {
			stop = estop;
		}
		exons.add(e);
	}

	public void addUTR(UTR u) {
		if (UTRs == null) {
			UTRs = new ArrayList<>();
		}
		UTRs.add(u);
	}

	public ArrayList<UTR> getUTRs() {
		return UTRs;
	}

	@Override
	public String toString() {
		return "Transcript{" + "name=" + name + ", chromosome=" + chromosome + ", strand=" + strand + ", gene=" + gene.getName() + ", start=" + start + ", stop=" + stop + '}';
	}


	public void setExonRank(Exon currExo, Integer exonrankintranscript) {
		if (exonRanks == null) {
			exonRanks = new HashMap<>();
		}
		exonRanks.put(currExo, exonrankintranscript);

	}

	public Exon[] getExonsRanked() {
		if (exonRanks == null) {
			return null;
		}

//		Exon[] sortedexons = new Exon[exons.size()];
//		Vector<StringIntegerObject> exonnames = new Vector<StringIntegerObject>();
//
//		for (Exon key : exons) {
////            System.out.println(key+""+exons.get(key).getRank());
//			Integer exonrank = exonRanks.get(key);
//			exonnames.add(new StringIntegerObject(key.toString(), exonrank));
//		}
//		StringIntegerObjectSorterSortOnInteger sorter = new StringIntegerObjectSorterSortOnInteger();
//		sorter.sort(exonnames);
//		for (int i = 0; i < exonnames.size(); i++) {
//			sortedexons[i] = exons.get(exonnames.get(i).stringValue);
//		}
		return null;
	}
}
