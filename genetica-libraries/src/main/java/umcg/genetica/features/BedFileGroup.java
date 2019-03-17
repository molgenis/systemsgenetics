package umcg.genetica.features;


import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

/**
 * Created by hwestra on 2/24/17.
 */
public class BedFileGroup {

//	HashMap<String, Integer> groupIndex = new HashMap<String, Integer>();
//	ArrayList<String> groupNames = new ArrayList<String>();
//	ArrayList<Annotation> annotations = new ArrayList<>();
//
//	public void loadGroupsFromFile(String groupfile, ArrayList<Feature> regions) throws IOException {
//		// read annotations in the regions
//		groupIndex = new HashMap<String, Integer>();
//		groupNames = new ArrayList<String>();
//		annotations = new ArrayList<>();
//		BedFileReader reader = new BedFileReader();
//
//		int groupctr = 0;
//		if (groupfile != null) {
//			System.out.println("Getting annotations from: " + groupfile);
//			String defaultgroup = "default";
//			TextFile tf = new TextFile(groupfile, TextFile.R);
//			String[] elems = tf.readLineElems(TextFile.tab);
//			while (elems != null) {
//				String file = null;
//				String name = null;
//				String group = null;
//				if (elems.length == 1) {
//					// there is only one group 'default'
//					file = elems[0];
//					group = defaultgroup;
//				} else if (elems.length == 3) {
//					group = elems[2];
//					name = elems[1];
//					file = elems[0];
//				}
//
//
//				if (file != null && Gpio.exists(file)) {
//					Integer groupId = groupIndex.get(group);
//					if (groupId == null) {
//						groupIndex.put(group, groupctr);
//						groupNames.add(group);
//						groupId = groupctr;
//						groupctr++;
//					}
//
//					Annotation a = new Annotation();
//					a.features = reader.readAsList(file, regions);
//					a.name = name;
//					a.filename = file;
//					a.groupId = groupId;
//					annotations.add(a);
//				}
//				elems = tf.readLineElems(TextFile.tab);
//			}
//			tf.close();
//			System.out.println(annotations.size() + " annotations, " + groupNames.size() + " groups of annotations.");
//		}
//	}
//
//	public Integer getGroupId(String groupName) {
//		return groupIndex.get(groupName);
//	}
//
//	public ArrayList<String> getGroups() {
//		return groupNames;
//	}
//
//
//	public ArrayList<Feature> getJointFeaturesForGroup(Integer groupId, Feature region) {
//		ArrayList<Feature> features = new ArrayList<>();
//		for (Annotation a : annotations) {
//			if (groupId == null || a.groupId.equals(groupId)) {
//				for (Feature f : a.features) {
//					if (region == null || region.overlaps(f)) {
//						features.add(f);
//					}
//				}
//			}
//		}
//		return features;
//	}
//
//	public ArrayList<ArrayList<Feature>> getFeaturesForGroup(int groupId, Feature region) {
//		ArrayList<ArrayList<Feature>> features = new ArrayList<>();
//		for (Annotation a : annotations) {
//			if (a.groupId.equals(groupId)) {
//				ArrayList<Feature> ann = new ArrayList<>();
//				for (Feature f : a.features) {
//					if (region == null || region.overlaps(f)) {
//						ann.add(f);
//					}
//				}
//				features.add(ann);
//			}
//		}
//		return features;
//	}
//
//	public int size() {
//		return groupNames.size();
//	}
//
//	public String getGroupName(int i) {
//		return groupNames.get(i);
//	}
//
//	private class Annotation {
//		Integer groupId;
//		ArrayList<Feature> features;
//		String name;
//		String filename;
//
//	}
}
