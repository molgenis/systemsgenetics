/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.variantAnnotator;

import umcg.genetica.collections.intervaltree.Range;

/**
 *
 * @author Patrick Deelen
 */
public interface GenomicRange extends Range {
	
	String getSeqname();
	
}
