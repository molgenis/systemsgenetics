/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.depict2;

import com.opencsv.CSVParser;
import com.opencsv.CSVParserBuilder;
import com.opencsv.CSVReader;
import com.opencsv.CSVReaderBuilder;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import org.apache.log4j.Logger;

/**
 *
 * @author patri
 */
public class PathwayAnnotations {

	private final HashMap<String, ArrayList<String>> annotations;
	private final int maxNumberOfAnnotations;
	
	private static final Logger LOGGER = Logger.getLogger(Depict2Options.class);

	public PathwayAnnotations(final File annotationFile) throws FileNotFoundException, IOException {

		if (annotationFile.exists()) {

			annotations = new HashMap<>();
			int maxNumberOfAnnotationsTmp = 0;
			final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
			final CSVReader reader = new CSVReaderBuilder(new BufferedReader(new FileReader(annotationFile))).withCSVParser(parser).withSkipLines(0).build();

			String[] nextLine;
			while ((nextLine = reader.readNext()) != null) {

				ArrayList<String> thisAnnotions = new ArrayList();

				for (int i = 1; i < nextLine.length; ++i) {
					thisAnnotions.add(nextLine[i]);
				}

				annotations.put(nextLine[0], thisAnnotions);

				if (thisAnnotions.size() > maxNumberOfAnnotationsTmp) {
					maxNumberOfAnnotationsTmp = thisAnnotions.size();
				}

			}

			maxNumberOfAnnotations = maxNumberOfAnnotationsTmp;

			LOGGER.debug("Found annotations at: " + annotationFile.getAbsolutePath() + ". Max annotations is: " + maxNumberOfAnnotations);

		} else {
			LOGGER.debug("Did not found annotations at: " + annotationFile.getAbsolutePath());
			annotations = null;
			maxNumberOfAnnotations = 0;
		}

	}

	public ArrayList<String> getAnnotationsForPathway(String pathway) {
		return annotations.get(pathway);
	}

	public int getMaxNumberOfAnnotations() {
		return maxNumberOfAnnotations;
	}

}
