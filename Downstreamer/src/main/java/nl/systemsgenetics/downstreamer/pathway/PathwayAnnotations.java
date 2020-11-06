/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.downstreamer.pathway;

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

import nl.systemsgenetics.downstreamer.DownstreamerOptions;
import org.apache.log4j.Logger;

/**
 *
 * @author patri
 */
public class PathwayAnnotations {

	private final HashMap<String, ArrayList<String>> annotations;
	private final ArrayList<String> annotationHeaders;
	private final String setName;

	private static final Logger LOGGER = Logger.getLogger(DownstreamerOptions.class);

	public PathwayAnnotations(final File annotationFile) throws FileNotFoundException, IOException {

		if (annotationFile.exists()) {

			annotations = new HashMap<>();
			final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
			final CSVReader reader = new CSVReaderBuilder(new BufferedReader(new FileReader(annotationFile))).withCSVParser(parser).withSkipLines(0).build();

			String[] nextLine = reader.readNext();
			if (nextLine == null) {
				throw new IOException("Empty annotation file: " + annotationFile.getAbsolutePath());
			}

			annotationHeaders = new ArrayList<>(nextLine.length - 1);
			setName=nextLine[0];
			for (int i = 1; i < nextLine.length; ++i) {
				annotationHeaders.add(nextLine[i]);
			}
			
			final int numberOfAnnotations = annotationHeaders.size();

			while ((nextLine = reader.readNext()) != null) {

				ArrayList<String> thisAnnotions = new ArrayList<>(numberOfAnnotations);

				for (int i = 1; i < nextLine.length; ++i) {
					thisAnnotions.add(nextLine[i]);
				}

				annotations.put(nextLine[0], thisAnnotions);

			}

			LOGGER.debug("Found annotations at: " + annotationFile.getAbsolutePath() + ". Number of annotations is: " + annotationHeaders.size());

		} else {
			LOGGER.debug("Did not found annotations at: " + annotationFile.getAbsolutePath());
			annotations = null;
			setName = null;
			annotationHeaders = null;

		}

	}

	public ArrayList<String> getAnnotationsForPathway(String pathway) {
		return annotations.get(pathway);
	}

	public int getMaxNumberOfAnnotations() {
		if (annotationHeaders == null) {
			return 0;
		}
		return annotationHeaders.size();
	}

	public String getSetName() {
		return setName;
	}

	public ArrayList<String> getAnnotationHeaders() {
		return annotationHeaders;
	}

}
