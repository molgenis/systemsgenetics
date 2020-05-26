/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.metaqtl3.containers;

/**
 *
 * @author harm-jan
 */
public class TextResult implements Comparable<TextResult> {
    public double pvalue;
    public int pid;
    public String description;
    public double zscore;

    @Override
     public int compareTo(TextResult o) {
        if(pvalue == o.pvalue){
            if(zscore == o.zscore){
                return 0;
            } else if(zscore > o.zscore){
                return 1;
            } else {
                return -1;
            }
        } else if(pvalue > o.pvalue){
            return 1;
        } else {
            return -1;
        }
        
    }
    
    public boolean equals(TextResult o) {
        if(pvalue == o.pvalue){
            if(zscore == o.zscore){
                return true;
            } else {
                return false;
            } 
        } else {
            return false;
        }
    }  
    
    @Override
    public String toString(){
        return pvalue+"\t"+description;
    }
    
            
}
