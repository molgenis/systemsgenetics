/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.features;

import java.util.Comparator;

/**
 * @author hwestra
 */
public class FeatureComparator implements Comparator<Feature> {

	boolean allowOverlap = false;

	public FeatureComparator() {
	}

	public FeatureComparator(boolean allowOverlap) {
		this.allowOverlap = allowOverlap;
	}

	public void setAllowOverlap(boolean b) {
		this.allowOverlap = b;
	}

	@Override
	public int compare(Feature obj1, Feature obj2) {
		if (obj1.equals(obj2)) {
			return 0;
		}
		if (allowOverlap) {
			if (obj1.overlaps(obj2)) {
				return 0;
			}
		}

		if (obj1.getChromosome().getNumber() > obj2.getChromosome().getNumber()) {
			return 1;
		} else if (obj1.getChromosome().getNumber() < obj2.getChromosome().getNumber()) {
			return -1;
		}

		if (obj1.getStart() > obj2.getStart()) {
			return 1;
		} else if (obj1.getStart() < obj2.getStart()) {
			return -1;
		}

//		if (obj1.getStrand() == null || obj2.getStrand() == null) {
//			return 0;
//		} else if (obj1.getStrand().getNumber() > obj2.getStrand().getNumber()) {
//			return 1;
//		} else if (obj1.getStrand().getNumber() < obj2.getStrand().getNumber()) {
//			return -1;
//		}

		return 0;



	}

}
