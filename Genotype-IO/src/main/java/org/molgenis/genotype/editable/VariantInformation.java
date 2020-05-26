/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.editable;

import org.molgenis.genotype.Alleles;

/**
 *
 * @author Patrick Deelen
 */
public class VariantInformation {
	
	private final String variantId;
	private final int startPos;
	private final String sequenceName;
	private final Alleles alleles;

	public VariantInformation(String variantId, int startPos, String sequenceName, Alleles alleles) {
		this.variantId = variantId;
		this.startPos = startPos;
		this.sequenceName = sequenceName.intern();
		this.alleles = alleles;
	}

	public String getVariantId() {
		return variantId;
	}

	public int getStartPos() {
		return startPos;
	}

	public String getSequenceName() {
		return sequenceName;
	}

	public Alleles getAlleles() {
		return alleles;
	}

	@Override
	public int hashCode() {
		int hash = 3;
		hash = 83 * hash + this.startPos;
		hash = 83 * hash + (this.sequenceName != null ? this.sequenceName.hashCode() : 0);
		return hash;
	}

	@Override
	public boolean equals(Object obj) {
		if (obj == null) {
			return false;
		}
		if (getClass() != obj.getClass()) {
			return false;
		}
		final VariantInformation other = (VariantInformation) obj;
		if (this.startPos != other.startPos) {
			return false;
		}
		if ((this.sequenceName == null) ? (other.sequenceName != null) : !this.sequenceName.equals(other.sequenceName)) {
			return false;
		}
		if (this.alleles != other.alleles && (this.alleles == null || !this.alleles.sameAlleles(other.alleles))) {
			return false;
		}
		return true;
	}

	
	
	
}
