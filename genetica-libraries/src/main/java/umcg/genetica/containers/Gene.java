/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.containers;

import java.util.HashMap;

/**
 *
 * @author harmjan
 */
public class Gene {
    private Chromosome parentChromosome;
    private HashMap<String, Transcript> transcripts;
    private int start;
    private int end;
    private int strand;
    private String name;
    private String annotation;

    /**
     * @return the parentChromosome
     */
    public Chromosome getParentChromosome() {
        return parentChromosome;
    }

    /**
     * @param parentChromosome the parentChromosome to set
     */
    public void setParentChromosome(Chromosome parentChromosome) {
        this.parentChromosome = parentChromosome;
    }

    /**
     * @return the transcripts
     */
    public HashMap<String, Transcript> getTranscripts() {
        return transcripts;
    }

    /**
     * @param transcripts the transcripts to set
     */
    public void setTranscripts(HashMap<String, Transcript> transcripts) {
        this.transcripts = transcripts;
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

    public void addTranscript(Transcript currTra) {
        if(transcripts == null){
            transcripts = new HashMap<String, Transcript>();
        }
        transcripts.put(currTra.getName(), currTra);
    }

    public void setAnnotation(String annotation){
	this.annotation = annotation;
    }

    public String getAnnotation(){
	return this.annotation;
    }
}
