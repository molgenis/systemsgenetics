/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.oxford;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import org.apache.commons.io.IOUtils;
import org.apache.log4j.Logger;
import static org.molgenis.genotype.GenotypeData.SAMPLE_MISSING_RATE_DOUBLE;
import org.molgenis.genotype.Sample;
import org.molgenis.genotype.annotation.Annotation;
import static org.molgenis.genotype.annotation.Annotation.Type.BOOLEAN;
import static org.molgenis.genotype.annotation.Annotation.Type.FLOAT;
import static org.molgenis.genotype.annotation.Annotation.Type.INTEGER;
import static org.molgenis.genotype.annotation.Annotation.Type.STRING;
import org.molgenis.genotype.annotation.SampleAnnotation;
import org.molgenis.genotype.util.Utils;
import org.molgenis.io.csv.CsvReader;
import org.molgenis.util.tuple.Tuple;

/**
 *
 * @author Patrick Deelen
 */
public class OxfordSampleFile {

	private static final Logger LOGGER = Logger.getLogger(OxfordSampleFile.class);
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

		SampleAnnotation missingAnnotation = new SampleAnnotation(SAMPLE_MISSING_RATE_DOUBLE, "missing",
				"Missing data proportion of each individual", Annotation.Type.FLOAT, SampleAnnotation.SampleAnnotationType.OTHER, false);
		sampleAnnotations.put(missingAnnotation.getId(), missingAnnotation);

		CsvReader reader = null;
		try {
			reader = new CsvReader(sampleFile, ' ');

			List<String> colNames = Utils.iteratorToList(reader.colNamesIterator());
			Tuple dataTypes = reader.iterator().next();
			for (int i = 3; i < colNames.size(); i++) {
				SampleAnnotation annotation = null;
				if (dataTypes.getString(i).equalsIgnoreCase("D")) {
					annotation = new SampleAnnotation(colNames.get(i), colNames.get(i), "", Annotation.Type.STRING,
							SampleAnnotation.SampleAnnotationType.COVARIATE, false);

				} else if (dataTypes.getString(i).equalsIgnoreCase("C")) {
					annotation = new SampleAnnotation(colNames.get(i), colNames.get(i), "", Annotation.Type.FLOAT,
							SampleAnnotation.SampleAnnotationType.COVARIATE, false);
				} else if (dataTypes.getString(i).equalsIgnoreCase("P")) {
					annotation = new SampleAnnotation(colNames.get(i), colNames.get(i), "", Annotation.Type.FLOAT,
							SampleAnnotation.SampleAnnotationType.PHENOTYPE, false);
				} else if (dataTypes.getString(i).equalsIgnoreCase("B")) {
					annotation = new SampleAnnotation(colNames.get(i), colNames.get(i), "", Annotation.Type.BOOLEAN,
							SampleAnnotation.SampleAnnotationType.PHENOTYPE, false);
				} else {
					LOGGER.warn("Unknown datatype [" + dataTypes.getString(i) + "]");
				}

				if (annotation != null) {
					sampleAnnotations.put(annotation.getId(), annotation);
				}
			}
		} finally {
			IOUtils.closeQuietly(reader);
		}

	}

	private void loadSamples() {

		CsvReader reader = null;
		try {
			reader = new CsvReader(sampleFile, ' ');
			Iterator<Tuple> it = reader.iterator();

			it.next();// Datatype row

			while (it.hasNext()) {
				Tuple tuple = it.next();

				String familyId = tuple.getString(0);
				String sampleId = tuple.getString(1);

				Map<String, Object> annotationValues = new LinkedHashMap<String, Object>();
				annotationValues.put(SAMPLE_MISSING_RATE_DOUBLE, tuple.getString(2).equals("NA") ? Double.NaN : tuple.getDouble(2));


				for (String colName : sampleAnnotations.keySet()) {

					if (colName.equals(SAMPLE_MISSING_RATE_DOUBLE)) {
						continue;//already done
					}

					SampleAnnotation annotation = sampleAnnotations.get(colName);

					Object value = null;
					if (!tuple.getString(annotation.getName()).equalsIgnoreCase("NA")) {
						switch (annotation.getType()) {
							case STRING:
								value = tuple.getString(colName);
								break;
							case INTEGER:
								value = tuple.getInt(colName);
								break;
							case BOOLEAN:
								if (tuple.getString(colName).equals("-9")) {
									value = null;
								} else {
									value = tuple.getBoolean(colName);
								}
								break;
							case FLOAT:
								value = tuple.getDouble(colName);
								break;
							default:
								LOGGER.warn("Unsupported data type encountered for column [" + colName + "]");
						}
					}

					annotationValues.put(colName, value);
				}

				samples.add(new Sample(sampleId, familyId, annotationValues));
			}

		} catch (FileNotFoundException e) {
			throw new RuntimeException("File [" + sampleFile.getAbsolutePath() + "] does not exists", e);
		} finally {
			IOUtils.closeQuietly(reader);
		}

	}

	public List<Sample> getSamples() {
		return Collections.unmodifiableList(samples);
	}

	public Map<String, SampleAnnotation> getSampleAnnotations() {
		return Collections.unmodifiableMap(sampleAnnotations);
	}
}
