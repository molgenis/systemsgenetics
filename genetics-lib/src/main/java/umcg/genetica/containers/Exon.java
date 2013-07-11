/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.containers;

import java.util.HashSet;

/**
 *
 * @author harmjan
 */
public class Exon {
    private HashSet<Transcript> parentTranscripts;
    private Chromosome parentChromosome;
    private int start;
    private int end;
    private int strand;
    private String name;

    /**
     * @return the parentTranscript
     */
    public Transcript[] getParentTranscript() {
        return parentTranscripts.toArray(new Transcript[0]);
    }

    /**
     * @param parentTranscript the parentTranscript to set
     */
    public void setParentTranscript(Transcript parentTranscript) {
        if(parentTranscripts == null){
	    parentTranscripts = new HashSet<Transcript>();
	}
	parentTranscripts.add(parentTranscript);
    }

    /**
     * @return the start
     */
    public int getStart() {
        return start;
    }

    /**
     * @param start the start to set
     */
    public void setStart(int start) {
        this.start = start;
    }

    /**
     * @return the end
     */
    public int getEnd() {
       return end;
    }

    /**
     * @param end the end to set
     */
    public void setEnd(int end) {
        this.end = end;
    }

    /**
     * @return the strand
     */
    public int getStrand() {
        return strand;
    }

    /**
     * @param strand the strand to set
     */
    public void setStrand(int strand) {
        this.strand = strand;

    }

    /**
     * @return the name
     */
    public String getName() {
        return name;
    }

    /**
     * @param name the name to set
     */
    public void setName(String name) {
        this.name = name;
    }

    public int getLength() {
        return (end-start) + 1;
    }

    public void setParentChromosome(Chromosome currChr) {
	if(parentChromosome!=null){
	    if(parentChromosome != currChr){
		System.err.println("WARNING: transcipt maps to multiple chromosomes: "+currChr.getName());
	    }
	}
	parentChromosome = currChr;
    }

    public Chromosome getParentChromosome() {
	return parentChromosome;
    }

}
