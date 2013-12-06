/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.util;

/**
 *
 * @author Patrick Deelen
 */
public class ChrPos {
	
	private final String chr;
	private final int pos;

	public ChrPos(String chr, int pos) {
		this.chr = chr;
		this.pos = pos;
	}

	public String getChr() {
		return chr;
	}

	public int getPos() {
		return pos;
	}
	
	
	
}
