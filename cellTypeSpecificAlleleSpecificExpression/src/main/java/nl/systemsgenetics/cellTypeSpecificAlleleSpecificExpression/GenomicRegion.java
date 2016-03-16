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
 * annotation is the name of the region (gene name or whatever.)
 * sequence is the name of the sequence or chromosome or construct
 * startPosition is the start of the gene region
 * endPosition is the end of the gene region
 * 
 * testStart is the start of the testing region
 * testEnd is the start of the testing region
 */

public class GenomicRegion {
    private String annotation;
    private String sequence;
    private int startPosition = -1;
    private int endPosition   = -1;
    private int testStart     = -1;
    private int testEnd       = -1;
    
    private boolean hasTestRegion = false;
    
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
        
        if(startPosition < 0){
            throw new IllegalDataException("startPosition was negative, is wrong.");
        }
        
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
        
        if(endPosition < 0){
            throw new IllegalDataException("endPosition was negative, invalid value");
        }
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
     * @return the testStart
     */
    public int getTestStart() {
        
        if(testStart < 0){
            throw new IllegalDataException("testStart was negative, is wrong.");
        }
        
        return testStart;
    }

    /**
     * @param testStart the testStart to set
     */
    public void setTestStart(int testStart) throws IllegalDataException {
        
        if(testStart < 0){
            throw new IllegalDataException("startPosition cannot be negative");
        }
        
        this.testStart = testStart;
    }
    
        /**
     * @return the testStart
     */
    public int getTestEnd() {
        
        if(testEnd < 0){
            throw new IllegalDataException("testEnd was negative, is wrong.");
        }
        
        return testEnd;
    }

    /**
     * @param testEnd the testStart to set
     */
    public void setTestEnd(int testEnd) throws IllegalDataException {
        
        if(testEnd < 0){
            throw new IllegalDataException("testEnd cannot be negative");
        }
        
        this.testEnd = testEnd;
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

    /**
     * @return the hasTestRegion
     */
    public boolean HasTestRegion() {
        return hasTestRegion;
    }

    /**
     * @param hasTestRegion the hasTestRegion to set
     */
    public void setHasTestRegion(boolean hasTestRegion) {
        this.hasTestRegion = hasTestRegion;
    }
    
    
    
}
