/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.containers;

import java.util.ArrayList;
import java.util.HashMap;

/**
 *
 * @author harmjan
 */
public class Chromosome {
    private String name;
    private HashMap<String, Gene> geneHash;
    private ArrayList<Gene> genes;
    
    public Chromosome(String name){
        
        this.name = name;
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

    /**
     * @return the genes
     */
    public HashMap<String, Gene> getGenesHash() {
        return geneHash;
    }

    
    
    /**
     * @param genes the genes to set
     */
    public void setGenes(HashMap<String, Gene> genes) {
        this.geneHash = genes;
    }

    public void addGene(Gene currGen) {
        if(geneHash == null){
            genes = new ArrayList<Gene>();
            geneHash = new HashMap<String, Gene>();
        }
        genes.add(currGen);
        geneHash.put(currGen.getName(), currGen);
    }

    public ArrayList<Gene> getGenes() {
        return genes;
    }

    public void setGenes(ArrayList<Gene> genes) {
        this.genes = genes;
    }

}
