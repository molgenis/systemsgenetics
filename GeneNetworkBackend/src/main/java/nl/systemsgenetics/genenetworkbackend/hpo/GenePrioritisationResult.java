/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.genenetworkbackend.hpo;

/**
 *
 * @author patri
 */
public class GenePrioritisationResult implements Comparable<GenePrioritisationResult> {

	private final String ensg;
	private final String symbol;
	private final double geneScore;

	public GenePrioritisationResult(String ensg, String symbol, double geneScore) {
		this.ensg = ensg;
		this.symbol = symbol;
		this.geneScore = geneScore;
	}

	public String getEnsg() {
		return ensg;
	}

	public String getSymbol() {
		return symbol;
	}

	public double getGeneScore() {
		return geneScore;
	}

	@Override
	public int compareTo(GenePrioritisationResult o) {
		//Code from double comapre but reversed
		if (o.geneScore < this.geneScore)
            return -1;           // Neither val is NaN, thisVal is larger
        if (o.geneScore > this.geneScore)
            return 1;            // Neither val is NaN, thisVal is smaler

        // Cannot use doubleToRawLongBits because of possibility of NaNs.
        long thisBits    = Double.doubleToLongBits(this.geneScore);
        long anotherBits = Double.doubleToLongBits(o.geneScore);

        return (thisBits == anotherBits ?  0 : // Values are equal
                (anotherBits < thisBits ? -1 : // (-0.0, 0.0) or (!NaN, NaN)
                 1));     
	}
	
}
