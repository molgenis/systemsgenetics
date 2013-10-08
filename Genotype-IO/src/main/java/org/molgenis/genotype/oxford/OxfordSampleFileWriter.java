/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.oxford;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import org.apache.log4j.Logger;
import org.molgenis.genotype.GenotypeData;
import org.molgenis.genotype.Sample;
import org.molgenis.genotype.annotation.SampleAnnotation;
import org.molgenis.genotype.annotation.SampleAnnotation.SampleAnnotationType;
import org.molgenis.genotype.annotation.SexAnnotation;

/**
 *
 * @author Patrick Deelen
 */
public class OxfordSampleFileWriter {

	public static final Charset FILE_ENCODING = Charset.forName("UTF-8");
	public static final String LINE_ENDING = "\n";
	private static final char SEPARATOR = ' ';
	private static final Logger LOG = Logger.getLogger(OxfordSampleFileWriter.class);

	public static void writeSampleFile(File sampleFile, GenotypeData genotypeData, HashMap<Sample, Float> sampleMissingness) throws IOException {

		Writer sampleWriter = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(sampleFile),
				FILE_ENCODING));

		// Write headers
		sampleWriter.append("ID_1");
		sampleWriter.append(SEPARATOR);
		sampleWriter.append("ID_2");
		sampleWriter.append(SEPARATOR);
		sampleWriter.append("missing");

		List<String> colNames = new ArrayList<String>();
		List<String> dataTypes = new ArrayList<String>();

		for (SampleAnnotation annotation : genotypeData.getSampleAnnotations()) {
			if (annotation.getId().equals(GenotypeData.SAMPLE_MISSING_RATE_FLOAT)) {
				continue;
			}

			if (annotation.getSampleAnnotationType() == SampleAnnotationType.COVARIATE) {
				switch (annotation.getType()) {
					case INTEGER:
						colNames.add(annotation.getId());
						dataTypes.add("D");
						break;
					case STRING:
						colNames.add(annotation.getId());
						dataTypes.add("D");
						break;
					case SEX:
						colNames.add(annotation.getId());
						dataTypes.add("D");
						break;
					case FLOAT:
						colNames.add(annotation.getId());
						dataTypes.add("C");
						break;
					default:
						LOG.warn("Unsupported covariate datatype [" + annotation.getType() + "]");
						break;
				}
			} else if (annotation.getSampleAnnotationType() == SampleAnnotationType.PHENOTYPE) {
				switch (annotation.getType()) {
					case BOOLEAN:
						colNames.add(annotation.getId());
						dataTypes.add("B");
						break;
					case FLOAT:
						colNames.add(annotation.getId());
						dataTypes.add("P");
						break;
					default:
						LOG.warn("Unsupported phenotype datatype [" + annotation.getType() + "]");
						break;
				}
			} else {
				LOG.warn("'OTHER' sample annotation type not supported by oxford sample file");
			}
		}

		for (String colName : colNames) {
			if(colName.equals(GenotypeData.SEX_SAMPLE_ANNOTATION_NAME)){
				colName = "sex";
			}
			sampleWriter.append(SEPARATOR);
			sampleWriter.append(colName);
		}

		sampleWriter.append(LINE_ENDING);

		// Write datatypes
		sampleWriter.append("0");
		sampleWriter.append(SEPARATOR);
		sampleWriter.append("0");
		sampleWriter.append(SEPARATOR);
		sampleWriter.append("0");

		for (String dataType : dataTypes) {
			sampleWriter.append(SEPARATOR);
			sampleWriter.append(dataType);
		}

		sampleWriter.append(LINE_ENDING);

		// Write values
		for (Sample sample : genotypeData.getSamples()) {
			sampleWriter.append(sample.getFamilyId() == null ? "NA" : sample.getFamilyId());
			sampleWriter.append(SEPARATOR);
			sampleWriter.append(sample.getId() == null ? "NA" : sample.getId());
			sampleWriter.append(SEPARATOR);
			sampleWriter.append(String.valueOf(sampleMissingness.get(sample)));

			for (String colName : colNames) {
				sampleWriter.append(SEPARATOR);
				sampleWriter.append(getValue(colName, sample, "NA"));
			}

			sampleWriter.append(LINE_ENDING);
		}
		sampleWriter.close();
	}

	private static String getValue(String colName, Sample sample, String nullValue) {
		Object value = sample.getAnnotationValues().get(colName);
		if (value == null) {
			return nullValue;
		}

		if (value instanceof Boolean) {
			return value.equals(true) ? "1" : "0";
		}

		if (value instanceof Double) {

			if (Double.isNaN((Double) value)) {
				return nullValue;
			}

			String result = value.toString();
			if (result.equals("0.0")) {
				result = "0";
			}

			return result;

		}

		if (value instanceof Float) {
			if (Float.isNaN((Float) value)) {
				return nullValue;
			}

			String result = value.toString();
			if (result.equals("0.0")) {
				result = "0";
			}

			return result;
		}

		if (value instanceof SexAnnotation) {
			SexAnnotation sexValue = (SexAnnotation) value;
			switch(sexValue){
				case UNKNOWN:
					return nullValue;
				default:
					return Byte.toString(sexValue.getPlinkSex());
			}
		}

		return value.toString();
	}
}
