/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.downstreamer.containers;

import htsjdk.samtools.util.Locatable;
import java.util.ArrayList;
import nl.systemsgenetics.downstreamer.gene.Gene;

/**
 *
 * @author patri
 */
public class GwasLocus implements Locatable {

	private final LeadVariant leadVariant;
	private final String contig;
	private final int start;
	private final int end;
	private final ArrayList<Gene> overlappingGenes = new ArrayList<>();

	public GwasLocus(LeadVariant leadVariant, String chr, int startPos, int stopPos) {
		this.leadVariant = leadVariant;
		this.contig = chr;
		this.start = startPos;
		this.end = stopPos;
	}
	
	public void addGene(Gene gene){
		overlappingGenes.add(gene);
	}

	@Override
	public int getStart() {
		return start;
	}

	@Override
	public int getEnd() {
		return end;
	}

	public LeadVariant getLeadVariant() {
		return leadVariant;
	}

	public ArrayList<Gene> getOverlappingGenes() {
		return overlappingGenes;
	}

	@Override
	public String getContig() {
		return contig;
	}
	
	public double getPvalue(){
		return leadVariant.getpValue();
	}
	
}
