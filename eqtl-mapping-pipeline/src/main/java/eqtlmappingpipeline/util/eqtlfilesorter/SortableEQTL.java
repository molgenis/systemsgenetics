/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.util.eqtlfilesorter;

/**
 *
 * @author harmjan
 */
public class SortableEQTL implements Comparable<SortableEQTL> {


    double pvalue    = 1;
    double absZScore = 0;
    String line = "";

    @Override
    public int compareTo(SortableEQTL t) {
	if(t.pvalue < pvalue){
	    return 1;
	} else if(t.pvalue == pvalue){
	    if(absZScore < t.absZScore){
		return 1;
	    } else {
		return -1;
	    }
	} else {
	    return -1;
	}
    }

    public boolean equals(SortableEQTL e){
	if(e.absZScore == absZScore){
	    return true;
	} else {
	    return false;
	}
    }

    @Override
    public String toString(){
	return line;
    }

}
