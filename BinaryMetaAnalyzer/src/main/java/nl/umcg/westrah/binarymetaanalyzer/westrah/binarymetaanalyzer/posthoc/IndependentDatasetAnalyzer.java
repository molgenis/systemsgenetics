package nl.umcg.westrah.binarymetaanalyzer.westrah.binarymetaanalyzer.posthoc;

import nl.umcg.westrah.binarymetaanalyzer.BinaryMetaAnalysis;
import nl.umcg.westrah.binarymetaanalyzer.BinaryMetaAnalysisSettings;

import java.io.IOException;
import java.util.ArrayList;

public class IndependentDatasetAnalyzer {


	public IndependentDatasetAnalyzer(String settingsFile, String textToReplace, String replaceTextWith, boolean usetmp) {
		BinaryMetaAnalysisSettings originalSettings = new BinaryMetaAnalysisSettings();
		originalSettings.parse(settingsFile, textToReplace, replaceTextWith);

		BinaryMetaAnalysisSettings settingsToUse = new BinaryMetaAnalysisSettings();
		settingsToUse.parse(settingsFile, textToReplace, replaceTextWith);


		ArrayList<String> annotations = originalSettings.getDatasetannotations();
		ArrayList<String> locations = originalSettings.getDatasetlocations();
		ArrayList<String> names = originalSettings.getDatasetnames();
		ArrayList<String> prefixes = originalSettings.getDatasetPrefix();


		for (int i = 0; i < annotations.size(); i++) {
			ArrayList<String> newannotations = new ArrayList<>();
			ArrayList<String> newprefixes = new ArrayList<>();
			ArrayList<String> newnames = new ArrayList<>();
			ArrayList<String> newlocations = new ArrayList<>();

			newannotations.add(annotations.get(i));
			newprefixes.add(prefixes.get(i));
			newnames.add(names.get(i));
			newlocations.add(locations.get(i));

			settingsToUse.setDatasetannotations(newnames);
			settingsToUse.setDatasetnames(newnames);
			settingsToUse.setDatasetlocations(newlocations);
			settingsToUse.setDatasetPrefix(newprefixes);

			settingsToUse.setOutput(originalSettings.getOutput() + "/" + names.get(i)+"/");

			BinaryMetaAnalysis b = new BinaryMetaAnalysis();
			settingsToUse.usetmp = usetmp;
			try {
				b.initialize(settingsToUse);
				b.run();
			} catch (IOException e) {
				e.printStackTrace();
			} catch (Exception e) {
				e.printStackTrace();
			}


		}
	}

}
