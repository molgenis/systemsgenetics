/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.oxford;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.io.IOUtils;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.LogManager;

import static org.molgenis.genotype.GenotypeData.SAMPLE_MISSING_RATE_FLOAT;

import org.molgenis.genotype.Sample;
import org.molgenis.genotype.annotation.Annotation;

import org.molgenis.genotype.annotation.SampleAnnotation;

import au.com.bytecode.opencsv.CSVReader;

/**
 *
 * @author Patrick Deelen
 */
public class OxfordSampleFile {

	private static final Logger LOGGER = LogManager.getLogger(OxfordSampleFile.class);
	private final File sampleFile;
	private final Map<String, SampleAnnotation> sampleAnnotations;
	private final List<Sample> samples;

	public OxfordSampleFile(File sampleFile) throws FileNotFoundException, IOException {

		if (sampleFile == null) {
			throw new IllegalArgumentException("sampleFile is null");
		}
		if (!sampleFile.isFile()) {
			throw new FileNotFoundException("sample file file not found at "
					+ sampleFile.getAbsolutePath());
		}
		if (!sampleFile.canRead()) {
			throw new IOException("Can not read sample file " + sampleFile.getAbsolutePath());
		}

		this.sampleFile = sampleFile;

		sampleAnnotations = new LinkedHashMap<String, SampleAnnotation>();
		samples = new ArrayList<Sample>();

		loadAnnotations();
		loadSamples();

	}

	private void loadAnnotations() throws IOException {

		SampleAnnotation missingAnnotation = new SampleAnnotation(SAMPLE_MISSING_RATE_FLOAT, "missing", "Missing data proportion of each individual", Annotation.Type.FLOAT, SampleAnnotation.SampleAnnotationType.OTHER, false);
		sampleAnnotations.put(missingAnnotation.getId(), missingAnnotation);

		CSVReader reader = new CSVReader(new InputStreamReader(new FileInputStream(sampleFile), Charset.forName("UTF-8")), ' ');
		try {
			
			List<String> colNames = Arrays.asList(reader.readNext());
			List<String> dataTypes = Arrays.asList(reader.readNext());
			
			for (int i = 3; i < colNames.size(); i++) {
				SampleAnnotation annotation = null;
				if (dataTypes.get(i).equalsIgnoreCase("D")) {
					annotation = new SampleAnnotation(colNames.get(i), colNames.get(i), "", Annotation.Type.STRING,
							SampleAnnotation.SampleAnnotationType.COVARIATE, false);

				} else if (dataTypes.get(i).equalsIgnoreCase("C")) {
					annotation = new SampleAnnotation(colNames.get(i), colNames.get(i), "", Annotation.Type.FLOAT,
							SampleAnnotation.SampleAnnotationType.COVARIATE, false);
				} else if (dataTypes.get(i).equalsIgnoreCase("P")) {
					annotation = new SampleAnnotation(colNames.get(i), colNames.get(i), "", Annotation.Type.FLOAT,
							SampleAnnotation.SampleAnnotationType.PHENOTYPE, false);
				} else if (dataTypes.get(i).equalsIgnoreCase("B")) {
					annotation = new SampleAnnotation(colNames.get(i), colNames.get(i), "", Annotation.Type.BOOLEAN,
							SampleAnnotation.SampleAnnotationType.PHENOTYPE, false);
				} else {
					LOGGER.warn("Unknown datatype [" + dataTypes.get(i) + "]");
				}

				if (annotation != null) {
					sampleAnnotations.put(annotation.getId(), annotation);
				}
			}
		} finally {
			IOUtils.closeQuietly(reader);
		}
	}

	private void loadSamples() throws IOException {
		CSVReader reader = new CSVReader(new InputStreamReader(new FileInputStream(sampleFile), Charset.forName("UTF-8")), ' ');
		try {
			// create col names index
			Map<String, Integer> colNamesIndex = new HashMap<String, Integer>();
			String[] colNames = reader.readNext();
			for(int i = 0; i < colNames.length; ++i) {
				colNamesIndex.put(colNames[i], i);
			}
			reader.readNext(); // skip datatypes row
			
			String[] tokens;
			while((tokens = reader.readNext()) != null) {

				String familyId = tokens[0];
				String sampleId = tokens[1];

				Map<String, Object> annotationValues = new LinkedHashMap<String, Object>();
				annotationValues.put(SAMPLE_MISSING_RATE_FLOAT, tokens[2].equals("NA") ? Float.NaN : Float.parseFloat(tokens[2]));


				for (String colName : sampleAnnotations.keySet()) {

					if (colName.equals(SAMPLE_MISSING_RATE_FLOAT)) {
						continue;//already done
					}

					SampleAnnotation annotation = sampleAnnotations.get(colName);

					Object value = null;
					String token = tokens[colNamesIndex.get(colName)];
					switch (annotation.getType()) {
						case STRING:
							value = token;
							break;
						case INTEGER:
							value = token.equals("NA") ? null : Integer.valueOf(token);
							break;
						case BOOLEAN:
							if (token.equals("-9") || token.equals("NA")) {
								value = null;
							} else {
								value = token.equals("1") ? true : false;
							}
							break;
						case FLOAT:
							value = token.equals("NA") || token.equals("-9") ? Float.NaN : Float.parseFloat(token);
							break;
						default:
							LOGGER.warn("Unsupported data type encountered for column [" + colName + "]");
					}


					annotationValues.put(colName, value);
				}

				samples.add(new Sample(sampleId, familyId, annotationValues));
			}
		} finally {
			IOUtils.closeQuietly(reader);
		}
	}

	public List<Sample> getSamples() {
		return Collections.unmodifiableList(samples);
	}

	public Map<String, SampleAnnotation> getSampleAnnotations() {
		return sampleAnnotations;
	}
}
