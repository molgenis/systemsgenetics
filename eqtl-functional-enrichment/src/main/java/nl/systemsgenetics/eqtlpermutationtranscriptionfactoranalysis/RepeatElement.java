/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.eqtlpermutationtranscriptionfactoranalysis;

/**
 *
 * @author Matthieu
 */
public class RepeatElement {
	private String chromosome;
	private int startPos;
	private int stopPos;
	private String elementClass;
	private String elementFamily;
	private String elementSubFamily;
	
	public RepeatElement(String chr, int start, int stop, String eclass, String efam, String esubfam){
		this.chromosome = chr;
		this.startPos = start;
		this.stopPos = stop;
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
	
	public String getElementClass(){
		return this.elementClass;
	}
	
	public String getElementFamily(){
		return this.elementFamily;
	}
	
	public String getElementSubFamily(){
		return this.elementSubFamily;
	}
	
	public int getLength(){
		return(this.stopPos - this.startPos);
	}
	
	public boolean snpIsInBoundary(int snpPos){
		return(snpPos >= this.startPos && snpPos <= this.stopPos);
	}
	
	public boolean regionIsInBoundary(int regionStart, int regionStop){
		return( (this.startPos >= regionStart && this.stopPos <= regionStart)
				|| (this.startPos >= regionStop && this.stopPos >= regionStop) );
	}
}