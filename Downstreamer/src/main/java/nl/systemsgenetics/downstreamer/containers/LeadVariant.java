/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.downstreamer.containers;

import htsjdk.samtools.util.Locatable;

/**
 *
 * @author patri
 */
public class LeadVariant implements Locatable{

	final String variantId;
	final String chr;
	final int pos;
	final double pValue;

	public LeadVariant(String variantId, String chr, int pos, double pValue) {
		this.variantId = variantId;
		this.chr = chr;
		this.pos = pos;
		this.pValue = pValue;
	}

	public String getVariantId() {
		return variantId;
	}

	@Override
	public String getContig() {
		return chr;
	}

	public int getPos() {
		return pos;
	}

	public double getpValue() {
		return pValue;
	}

	@Override
	public int getStart() {
		return pos;
	}

	@Override
	public int getEnd() {
		return pos;
	}
	
	
	
}
