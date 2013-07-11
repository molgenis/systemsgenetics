/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.io.gwascatalog;

import java.util.HashMap;
import java.util.HashSet;

/**
 *
 * @author harmjan
 */
public class GWASSNP {

    private HashSet<GWASTrait> associatedTraits = new HashSet<GWASTrait>();
    private HashSet<GWASPublication> publishedIn = new HashSet<GWASPublication>();
    private HashMap<GWASTrait, String> riskAllele = new HashMap<GWASTrait, String>();
    private int id;
    private String name;
    private byte chr;
    private int position;
    private String locus;
    private HashMap<GWASTrait, Double> pvalPerTrait;

    public String getName() {
        return name;
    }

    public String getRiskAllele(GWASTrait t) {
        return riskAllele.get(t);
    }

    /**
     * @return the associatedTraits
     */
    public HashSet<GWASTrait> getAssociatedTraits() {
        return associatedTraits;
    }



    public GWASTrait[] getAssociatedTraitsArray() {
        GWASTrait[] r = new GWASTrait[associatedTraits.size()];
        r = associatedTraits.toArray(r);
        return r;
    }

    /**
     * @param associatedTraits the associatedTraits to set
     */
    public void setAssociatedTraits(HashSet<GWASTrait> associatedTraits) {
        this.associatedTraits = associatedTraits;
    }

    /**
     * @return the publishedIn
     */
    public HashSet<GWASPublication> getPublishedIn() {
        return publishedIn;
    }

    /**
     * @param publishedIn the publishedIn to set
     */
    public void setPublishedIn(HashSet<GWASPublication> publishedIn) {
        this.publishedIn = publishedIn;
    }

    /**
     * @return the riskAllele
     */
    public HashMap<GWASTrait, String> getRiskAllele() {
        return riskAllele;
    }

    /**
     * @param riskAllele the riskAllele to set
     */
    public void setRiskAllele(HashMap<GWASTrait, String> riskAllele) {
        this.riskAllele = riskAllele;
    }

    /**
     * @return the id
     */
    public int getId() {
        return id;
    }

    /**
     * @param id the id to set
     */
    public void setId(int id) {
        this.id = id;
    }

    /**
     * @param name the name to set
     */
    public void setName(String name) {
        this.name = name;
    }

    /**
     * @return the chr
     */
    public byte getChr() {
        return chr;
    }

    /**
     * @param chr the chr to set
     */
    public void setChr(byte chr) {
        this.chr = chr;
    }

    /**
     * @return the position
     */
    public int getPosition() {
        return position;
    }

    /**
     * @param position the position to set
     */
    public void setPosition(int position) {
        this.position = position;
    }

    /**
     * @return the locus
     */
    public String getLocus() {
        return locus;
    }

    /**
     * @param locus the locus to set
     */
    public void setLocus(String locus) {
        this.locus = locus;
    }

    public void setPValueAssociatedWithTrait(GWASTrait gwasTraitObj, Double pval) {
	if(pvalPerTrait == null){
	    pvalPerTrait = new HashMap<GWASTrait, Double>();
	}
	pvalPerTrait.put(gwasTraitObj, pval);
    }

    public Double getPValueAssociatedWithTrait(GWASTrait t){
	return pvalPerTrait.get(t);
    }
    
    @Override
    public String toString() {
        return name;
    }
}
