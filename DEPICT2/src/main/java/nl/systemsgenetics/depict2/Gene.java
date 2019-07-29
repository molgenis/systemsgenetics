/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.depict2;

import java.util.Objects;

/**
 *
 * @author patri
 */
public class Gene {

	private final String gene;
	private final String chr;
	private final int start;
	private final int stop;
	private final String band;

	public Gene(String gene, String chr, int start, int stop, String band) {
		this.gene = gene;
		this.chr = chr.intern();
		this.start = start;
		this.stop = stop;
		this.band = band;
	}

	public String getGene() {
		return gene;
	}

	public String getChr() {
		return chr;
	}

	public int getStart() {
		return start;
	}

	public int getStop() {
		return stop;
	}

	public String getBand() {
		return band;
	}

	public String getChrAndArm() {

		StringBuilder sb = new StringBuilder(chr);
		if (!this.band.equals("")) {
			sb.append('_');
			sb.append(band.charAt(0));
		}
		return sb.toString().intern();

	}

	@Override
	public int hashCode() {
		int hash = 7;
		hash = 67 * hash + Objects.hashCode(this.gene);
		return hash;
	}

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
		final Gene other = (Gene) obj;
		if (this.start != other.start) {
			return false;
		}
		if (this.stop != other.stop) {
			return false;
		}
		if (!Objects.equals(this.gene, other.gene)) {
			return false;
		}
		if (!Objects.equals(this.chr, other.chr)) {
			return false;
		}
		if (!Objects.equals(this.band, other.band)) {
			return false;
		}
		return true;
	}

}
