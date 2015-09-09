/*
 * Copyright (C) 2015 Adriaan van der Graaf
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
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
