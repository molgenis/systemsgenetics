/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.cellTypeSpecificAlleleSpecificExpression;

import java.util.ArrayList;
import org.jdom.IllegalDataException;

/**
 *
 * @author adriaan
 * small class that contains genomic region.
 */
public class GenomicRegion {
    private String annotation;
    private String sequence;
    private int startPosition = -1;
    private int endPosition = -1;
    
    private ArrayList<String> snpInRegions;

    
    public void addSnpToRegion(String snpName){
        getSnpInRegions().add(snpName);
    }
    
    
    /**
     * @return the sequence
     */
    public String getSequence() {
        return sequence;
    }

    /**
     * @param sequence the sequence to set
     */
    public void setSequence(String sequence) {
        this.sequence = sequence;
    }

    /**
     * @return the startPosition
     */
    public int getStartPosition() {
        return startPosition;
    }

    /**
     * @param startPosition the startPosition to set
     */
    public void setStartPosition(int startPosition) throws IllegalDataException {
        
        if(startPosition < 0){
            throw new IllegalDataException("startPosition cannot be negative");
        }
        
        this.startPosition = startPosition;
    }

    /**
     * @return the endPosition
     */
    public int getEndPosition() {
        return endPosition;
    }

    /**
     * @param endPosition the endPosition to set
     */
    public void setEndPosition(int endPosition) {
      
        if(endPosition < 0){
            throw new IllegalDataException("endPosition cannot be negative");
        }
        
        this.endPosition = endPosition;
    }

    /**
     * @return the annotation
     */
    public String getAnnotation() {
        return annotation;
    }

    /**
     * @param annotation the annotation to set
     */
    public void setAnnotation(String annotation) {
        this.annotation = annotation;
    }

    /**
     * @return the snpInRegions
     */
    public ArrayList<String> getSnpInRegions() {
        return snpInRegions;
    }

    /**
     * @param snpInRegions the snpInRegions to set
     */
    public void setSnpInRegions(ArrayList<String> snpInRegions) {
        this.snpInRegions = snpInRegions;
    }
    
    
    
}
