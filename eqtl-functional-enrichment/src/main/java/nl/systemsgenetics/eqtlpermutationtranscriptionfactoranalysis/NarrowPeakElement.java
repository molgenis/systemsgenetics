/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.eqtlpermutationtranscriptionfactoranalysis;

/**
 *
 * @author Matthieu
 */
public class NarrowPeakElement {
	private String chromosome;
	private int startPos;
	private int stopPos;
	private double peakValue;
	private double log10PValue;
	private String infoLine;
	
	public NarrowPeakElement(String chr, int start, int stop, double peakValue, double pvalue, String info){
		this.chromosome = chr;
		this.startPos = start;
		this.stopPos = stop;
		this.peakValue = peakValue;
		this.log10PValue = pvalue;
		this.infoLine = info;
	}
	
	public String getChromosome(){
		return this.chromosome;
	}
	
	public int getStartPosition(){
		return this.startPos;
	}
	
	public int getStopPosition(){
		return this.stopPos;
	}
	
	public double getPeakValue(){
		return this.peakValue;
	}
	
	public double getLog10PValue(){
		return this.log10PValue;
	}
	
	public double getPValue(){
		return Math.pow(10, this.log10PValue);
	}
	
	public String getInfoLine(){
		return this.infoLine;
	}
}
