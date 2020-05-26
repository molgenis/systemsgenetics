/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.eqtlpermutationtranscriptionfactoranalysis;

/**
 *
 * @author Matthieu
 */
public class EQtl {
	
	private double pValue;
	private String eQtlName;
	private String eQtlChr;
	private int eQtlChrPos;
	private String probeName;
	private double zScore;
	
	
	
	public EQtl(double pv, String name, String chr, int pos, String probe, double score){
		this.pValue = pv;
		this.eQtlName = name;
		this.eQtlChr = chr;
		this.eQtlChrPos = pos;
		this.probeName = probe;
		this.zScore = score;
	}
	
	
	public double getPValue(){
		return this.pValue;
	}
	
	public String getEQtlName(){
		return this.eQtlName;
	}
	
	public String getEQtlChr(){
		return this.eQtlChr;
	}
	
	public int getEQtlChrPos(){
		return this.eQtlChrPos;
	}
	
	public String getProbeName(){
		return this.probeName;
	}
	
	public double getZScore(){
		return this.zScore;
	}
	
}
