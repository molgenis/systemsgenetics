/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.variantFilter;

import java.util.HashSet;
import org.molgenis.genotype.util.ChrPos;
import org.molgenis.genotype.variant.GeneticVariant;

/**
 *
 * @author Patrick Deelen
 */
public class VariantFilterSeqPos implements VariantFilter{

	private static final String SEP = ":";
	
	private final HashSet<String> toInclude;

	public VariantFilterSeqPos() {
		this.toInclude = new HashSet<String>();
	}
	
	public void addSeqPos(String seq, int pos){
		toInclude.add(seq + SEP + pos);
	}
	
	public void addSeqPos(ChrPos chrPos){
		addSeqPos(chrPos.getChr(), chrPos.getPos());
	}
	
	public void addSeqPos(GeneticVariant variant){
		toInclude.add(variant.getSequenceName()	+ SEP + variant.getStartPos());
	}
	
	@Override
	public boolean doesVariantPassFilter(GeneticVariant variant) {
		return toInclude.contains(variant.getSequenceName() + SEP + variant.getStartPos());
	}

	@Override
	public boolean doesIdPassFilter(String id) {
		return true;
	}
	
	
	
}
