/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.depict2;

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

}
