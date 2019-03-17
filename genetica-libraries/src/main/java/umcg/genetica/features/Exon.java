/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.features;

import umcg.genetica.enums.Chromosome;
import umcg.genetica.enums.Strand;

import java.util.ArrayList;

/**
 * @author Harm-Jan
 */
public class Exon extends Feature {
	
	private ArrayList<Transcript> transcripts;
	private final Gene gene;
	private FeatureType type;
	
	
	public Exon(String name, Chromosome chr, Strand strand, Gene gene, int start, int stop) {
		this.name = name;
		
		if (this.name != null) {
			type = FeatureType.parse(name);
		} else {
			type = FeatureType.EXON;
		}
		
		this.chromosome = chr;
		this.strand = strand;
		
		this.gene = gene;
		this.start = start;
		this.stop = stop;
		
	}
	
	public ArrayList<Transcript> getTranscripts() {
		return transcripts;
	}
	
	public Gene getGene() {
		return gene;
	}
	
	public void addTranscript(Transcript t) {
		if (this.transcripts == null) {
			this.transcripts = new ArrayList<Transcript>();
		}
		this.transcripts.add(t);
	}
	
	public FeatureType getType() {
		return type;
	}
	
	@Override
	public String toString() {
		return "Exon{" +
				"type=" + type +
				", chromosome=" + chromosome +
				", name='" + name + '\'' +
				", strand=" + strand +
				", start=" + start +
				", stop=" + stop +
				'}';
	}
}
