/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.variant;

import java.util.Comparator;
import org.molgenis.genotype.util.ChromosomeComparator;

/**
 *
 * @author Patrick Deelen
 */
public class GeneticVariantChrPosComparator implements Comparator<GeneticVariant>{

	@Override
	public int compare(GeneticVariant o1, GeneticVariant o2) {
		
		if(o1 == o2){
			return 0;
		} else {
			int chrCompareRes = ChromosomeComparator.comparator.compare(o1.getSequenceName(), o2.getSequenceName());
			if(chrCompareRes == 0){
				return o1.getStartPos() - o2.getStartPos();
			} else {
				return chrCompareRes;
			}
		}
		
	}
	
}
