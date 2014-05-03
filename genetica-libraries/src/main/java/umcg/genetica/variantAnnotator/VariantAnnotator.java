/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.variantAnnotator;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;

/**
 *
 * @author Patrick Deelen
 */
public class VariantAnnotator<E extends GenomicRange> {

	HashMap<String, ArrayList<E>> annotations;
	
	public VariantAnnotator() {
		annotations = new HashMap<String, ArrayList<E>>();
	}
	
	public void addAnnotaion(E annotation) throws Exception{
		ArrayList<E> chrAnnotations = annotations.get(annotation.getSeqname());
		if(chrAnnotations == null){
			chrAnnotations = new ArrayList<E>();
			annotations.put(annotation.getSeqname(), chrAnnotations);
			chrAnnotations.add(annotation);
		} else {
			E last = chrAnnotations.get(chrAnnotations.size() - 1);
			if(last.getStart() > annotation.getStart()){
				throw new Exception("Annotations not sorted");
			} 
			chrAnnotations.add(annotation);
		}
	}
	

	
	
}
