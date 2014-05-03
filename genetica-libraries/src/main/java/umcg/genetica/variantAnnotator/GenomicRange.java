/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.variantAnnotator;

/**
 *
 * @author Patrick Deelen
 */
public interface GenomicRange {
	
	String getSeqname();
	int getStart();
	int getEnd();
	
}
