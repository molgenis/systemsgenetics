package umcg.genetica.features;

import umcg.genetica.enums.Chromosome;
import umcg.genetica.enums.Strand;

import java.util.ArrayList;
import java.util.NavigableSet;
import java.util.TreeSet;

/**
 * Created by hwestra on 5/3/16.
 */
public class FeatureTree {

	// setting comparator to true removes all overlapping elements from the tree
	// also, when getting subsets from the tree, overlapping elements are not retrieved
	FeatureComparator comparator = new FeatureComparator(false);
	TreeSet<Feature> treeSet = new TreeSet<Feature>(comparator);
	ArrayList<Feature> features = new ArrayList<>();
	private String name;

	public FeatureTree(ArrayList<Feature> features, boolean mergedups) {
		if (mergedups) {
			FeatureMerger merger = new FeatureMerger();
			this.features = merger.merge(features, true);
			treeSet.addAll(this.features);
		} else {
			this.features = features;
			treeSet.addAll(features);
		}

		// switch the comparator for retrieval
		comparator.setAllowOverlap(true);
	}

	public FeatureTree(ArrayList<Feature> features) {
		this(features, false);
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

		NavigableSet<Feature> set = this.treeSet.subSet(left, true, right, true);
		left.setStrand(Strand.NEG);
		right.setStrand(Strand.NEG);
		set.addAll(this.treeSet.subSet(left, true, right, true));

		left.setStrand(Strand.NA);
		right.setStrand(Strand.NA);
		set.addAll(this.treeSet.subSet(left, true, right, true));
		return set;
	}

	public NavigableSet<Feature> getFeatureSet(Feature feat) {
		if (feat.getStop() == feat.getStart()) {
			return getFeatureSet(feat.getChromosome(), feat.getStart() - 1, feat.getStop() + 1);
		} else {
			return getFeatureSet(feat.getChromosome(), feat.getStart(), feat.getStop());
		}
	}

	public int getTreeSize() {
		return treeSet.size();
	}

	public int getListSize() {
		return features.size();
	}

	public String getName() {

		return name;
	}

	public TreeSet<Feature> getFeatureTree() {
		return treeSet;
	}

	public void setName(String name) {
		this.name = name;
	}
}
