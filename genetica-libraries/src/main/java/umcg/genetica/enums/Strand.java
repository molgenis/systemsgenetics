/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.enums;

/**
 * @author Harm-Jan
 */
public enum Strand {

	POS("+", 0), NEG("-", 1), NA("NA", 2);

	public static Strand parseStr(String string) {
		if (string.equals("+") || string.equals("1")) {
			return Strand.POS;
		}
		if (string.equals("-") || string.equals("-1")) {
			return Strand.NEG;
		}
		return Strand.NA;
	}

	int num;
	String strand;

	private Strand(String s, int num) {
		this.strand = s;
		this.num = num;
	}

	@Override
	public String toString() {
		return this.getName();
	}

	public int getNumber() {
		return num;
	}

	public String getName() {
		if (this == Strand.NEG) {
			return "-";
		} else if (this == Strand.POS) {
			return "+";
		} else {
			return "NA";
		}
	}

}
