/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.features;


import umcg.genetica.enums.Chromosome;
import umcg.genetica.enums.Strand;

import java.util.ArrayList;
import java.util.NavigableSet;
import java.util.TreeSet;

/**
 * @author Harm-Jan
 *         <p>
 *         I've found this class to be very unsafe and confusing.
 */
public class Track extends Feature {
	//
	private TreeSet<Feature> features;
	private ArrayList<Feature> allFeatures;

	public Track(String name) {
//		this(name, 0, Integer.MAX_VALUE);
	}

	public Track(String name, boolean removeOverlappingFeaturesFromTree) {
		this(name, 0, Integer.MAX_VALUE, null, removeOverlappingFeaturesFromTree);
	}

	public Track(String name, int start, int stop) {
		this(name, start, stop, null, false);
	}

	public Track(String name, int start, int stop, ArrayList<Feature> features) {
		this(name, start, stop, features, false);
	}

	public Track(ArrayList<Feature> features) {
		this(null, 0, Integer.MAX_VALUE, features, false);
		this.addFeatures(features);
	}

	public Track(ArrayList<Feature> features, boolean removeOverlappingFeaturesFromTree) {
		this(null, 0, Integer.MAX_VALUE, features, removeOverlappingFeaturesFromTree);
	}

	public Track(String name, int start, int stop, ArrayList<Feature> features, boolean removeOverlappingFeaturesFromTree) {
		this.name = name;
		this.features = new TreeSet<Feature>(new FeatureComparator(removeOverlappingFeaturesFromTree));
		this.allFeatures = new ArrayList<Feature>();
		this.start = start;
		this.stop = stop;
		if (features != null) {
			this.addFeatures(features);
		}
	}


	public void addFeature(Feature f) {
		features.add(f);
	}

	public Iterable<Feature> getFeatures() {
		return features;
	}

	public ArrayList<Feature> getAllFeatures() {
		return allFeatures;
	}

	public int getNrReads() {
		return features.size();
	}

	public NavigableSet<Feature> getFeatureSet(Chromosome chr, int start, int end) {
		Feature left = new Feature();
		Feature right = new Feature();
		left.setChromosome(chr);
		left.setStart(start);
		left.setStop(start);
		left.setStrand(Strand.POS);
		right.setChromosome(chr);
		right.setStart(end);
		right.setStop(end);
		right.setStrand(Strand.POS);

		NavigableSet<Feature> set = this.features.subSet(left, true, right, true);


		left.setStrand(Strand.NEG);
		right.setStrand(Strand.NEG);
		set.addAll(this.features.subSet(left, true, right, true));

		left.setStrand(Strand.NA);
		right.setStrand(Strand.NA);
		set.addAll(this.features.subSet(left, true, right, true));
		return set;
	}

	public void printNrFeatures() {
		System.out.println(this.features.size() + " features in track.");
	}

	public void setFeatures(NavigableSet<Feature> f) {
		this.features = new TreeSet<Feature>(new FeatureComparator(false));
		this.allFeatures = new ArrayList<Feature>();
		this.features.addAll(f);
		this.allFeatures.addAll(f);
	}

	public Track getSubset(Chromosome chr, int start, int stop) {
		Track t = new Track(this.name, start, stop);
		NavigableSet<Feature> set = getFeatureSet(chr, start, stop);

		t.setFeatures(set);
		return t;
	}

	public int getNrUniqueFeatures() {
		return this.features.size();
	}

	public int getNrFeatures() {
		return this.allFeatures.size();
	}

	public boolean containsFeature(Feature f) {
		return features.contains(f);
	}

	public void addFeatures(Track t) {
		this.features.addAll(t.features);
		this.allFeatures.addAll(t.features);
	}

	public void addFeatures(ArrayList<Feature> features) {
		// TODO: I should really kill this class.......
		if (this.features == null) {
			FeatureTree tree = new FeatureTree(features);
			this.features = tree.getFeatureTree();
			this.allFeatures = new ArrayList<>();
			this.allFeatures.addAll(features);
		} else {
			this.features.addAll(features);
			this.allFeatures.addAll(features);
		}
		System.out.println(this.features.size() + " features loaded into track");
	}

	public NavigableSet<Feature> getFeatureSet(Feature feat) {
		if (feat.getStop() == feat.getStart()) {
			return getFeatureSet(feat.getChromosome(), feat.getStart() - 1, feat.getStop() + 1);
		} else {
			return getFeatureSet(feat.getChromosome(), feat.getStart(), feat.getStop());
		}

	}

	@Override
	public String toString() {
		return "Track{" +
				"features=" + this.features.size() +
				", allFeatures=" + this.allFeatures.size() +
				'}';
	}
}
