package umcg.genetica.features;


import umcg.genetica.containers.Pair;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

/**
 * Created by hwestra on 5/13/15.
 */
public class FeatureMerger {

	public static ArrayList<Feature> merge(ArrayList<Feature> f, boolean returnCopies) {
		// make copies of the features first
		Feature[] featureArr = new Feature[f.size()];
		for (int i = 0; i < f.size(); i++) {
			if (returnCopies) {
				featureArr[i] = new Feature(f.get(i));
			} else {
				featureArr[i] = f.get(i);
			}
		}

		Arrays.sort(featureArr, new FeatureComparator(false));
//		System.out.println(featureArr.length);
		int overlap = 0;
		while (overlap != 0) {
			overlap = 0;
			for (int i = 0; i < featureArr.length; i++) {
				Feature f1 = featureArr[i];
				if (f1 != null) {
					for (int j = i + 1; j < featureArr.length; j++) {
						Feature f2 = featureArr[j];
						if (f2 != null) {
							if (!f1.equals(f2.getChromosome())) {
								break;// list is sorted: no need to continue merging
							}
							if (f1.overlaps(f2)) {
								if (f2.getStart() < f1.getStart()) {
									f1.setStart(f2.getStart());
								}
								if (f2.getStop() > f1.getStop()) {
									f1.setStop(f2.getStop());
								}

								featureArr[j] = null;
								overlap++;
							}
						}

					}
				}
			}

		}

		ArrayList<Feature> output = new ArrayList<Feature>();
		for (int i = 0; i < featureArr.length; i++) {
			if (featureArr[i] != null) {
				output.add(featureArr[i]);
			}
		}
//		System.out.println(output.size());
		return output;
	}


	// make subregions from areas in features1 that overlap with features2
	// this method assumes that there is no overlap WITHIN features1 or features2
	public static ArrayList<Pair<Feature, ArrayList<Feature>>> makeSubregionsUsingOverlap(ArrayList<Feature> features1,
																						  ArrayList<Feature> features2) {


		Collections.sort(features1, new FeatureComparator(false));
		Collections.sort(features2, new FeatureComparator(false));

		int lastj = 0;
		ArrayList<Pair<Feature, ArrayList<Feature>>> regionsOfAnnotation1OverlappingAnnotation2 = new ArrayList<Pair<Feature, ArrayList<Feature>>>();
		for (int i = 0; i < features2.size(); i++) {
			Feature f2 = features2.get(i);
			ArrayList<Feature> overlappingSubregions = new ArrayList<Feature>();

			for (int j = 0; j < features1.size(); j++) {

				Feature f1 = features1.get(j);
				if (f1.overlaps(f2)) {
					// determine which part of f1 overlaps f2
					// create new region from those coordinates

					int f1start = f1.getStart();
					int f1stop = f1.getStop();
					int f2start = f2.getStart();
					int f2stop = f2.getStop();

					Feature out = null;
//					if (f1start >= f2start && f1stop <= f2stop) {
//						// perfect overlap
//						//       |-----|
//						//     |--------|
//						out = new Feature(f1);
//					} else if (f1start <= f2start && f1stop >= f2start) {
//						// region between f2start and f1stop is overlapping
//						//      |-----|
//						//         |----|
//						// AND
//						//      |----------|
//						//         |----|
//						int sto = f1stop;
//						if (sto > f2stop) {
//							sto = f2stop;
//						}
//						int sta = f1start;
//						if (sta < f2start) {
//							sta = f2start;
//						}
//						out = new Feature(f2.getChromosome(), sta, sto);
//					} else if (f1start < f2stop && f1stop > f2stop) {
//						// region between f1start and f2stop is overlapping
//						//      |-----|
//						//    |----|
//
//						int sto = f1stop;
//						if (sto > f2stop) {
//							sto = f2stop;
//						}
//						int sta = f1start;
//						if (sta < f2start) {
//							sta = f2start;
//						}
//
//						out = new Feature(f2.getChromosome(), sta, sto);
//					}
					int sto = f1stop;
					if (sto > f2stop) {
						sto = f2stop;
					}
					int sta = f1start;
					if (sta < f2start) {
						sta = f2start;
					}
					out = new Feature(f2.getChromosome(), sta, sto);
					overlappingSubregions.add(out);


				}
			}
			regionsOfAnnotation1OverlappingAnnotation2.add(new Pair<Feature, ArrayList<Feature>>(f2, overlappingSubregions));
		}

		return regionsOfAnnotation1OverlappingAnnotation2;
	}
}
