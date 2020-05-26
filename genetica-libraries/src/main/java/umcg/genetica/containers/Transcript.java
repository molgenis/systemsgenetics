/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.containers;

import java.util.HashMap;
import java.util.Set;
import java.util.Vector;
import umcg.genetica.util.StringIntegerObjectSorterSortOnInteger;

/**
 *
 * @author harmjan
 */
public class Transcript {
    private Gene parentGene;
    private HashMap<String, Exon> exons;
    private HashMap<Exon, Integer> exonRanks;
    private Chromosome parentChromosome;
    private int start;
    private int end;
    private int strand;
    private String name;
    private String protein;

    /**
     * @return the parentGene
     */
    public Gene getParentGene() {
        return parentGene;
    }

    /**
     * @param parentGene the parentGene to set
     */
    public void setParentGene(Gene parentGene) {
        this.parentGene = parentGene;
    }

    /**
     * @return the exons
     */
    public HashMap<String, Exon> getExons() {
        return exons;
    }

    /**
     * @param exons the exons to set
     */
    public void setExons(HashMap<String, Exon> exons) {
        this.exons = exons;
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

    public void addExon(Exon currExo) {
        if(exons == null){
            exons = new HashMap<String, Exon>();
        }
        if(currExo == null){
            System.out.println("ERROR!");
            System.exit(0);
        }
        exons.put(currExo.getName(), currExo);
    }

    public void setProtein(String protein) {
        this.protein = protein;
    }

    public String getProtein(){
        return protein;
    }

    public Exon[] getExonsRanked() {
        Exon[] sortedexons = new Exon[exons.size()];
        Vector<StringIntegerObject> exonnames = new Vector<StringIntegerObject>();
        Set<String> keys = exons.keySet();
        for(String key: keys){
//            System.out.println(key+""+exons.get(key).getRank());
	    Integer exonrank = exonRanks.get(exons.get(key));
            exonnames.add(new StringIntegerObject(key, exonrank));
        }
        StringIntegerObjectSorterSortOnInteger sorter = new StringIntegerObjectSorterSortOnInteger();
        sorter.sort(exonnames);
        for(int i=0; i<exonnames.size(); i++){
            sortedexons[i] = exons.get( exonnames.get(i).stringValue );
        }
        return sortedexons;
    }

    public int getLength() {
        return end - start;
    }

    public void setParentChromosome(Chromosome currChr) {
	if(parentChromosome!=null){
	    if(parentChromosome != currChr){
		System.err.println("WARNING: transcipt maps to multiple chromosomes: "+currChr.getName());
	    }
	}
	parentChromosome = currChr;
    }

    public void setExonRank(Exon currExo, Integer exonrankintranscript) {
	if(exonRanks == null){
	    exonRanks = new HashMap<Exon, Integer>();
	}
	exonRanks.put(currExo, exonrankintranscript);
    }
}
