package org.molgenis.genotype.impute2;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.List;
import org.apache.commons.io.IOUtils;
import org.apache.log4j.Logger;
import org.molgenis.genotype.Allele;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.GenotypeData;
import org.molgenis.genotype.GenotypeWriter;
import org.molgenis.genotype.Sample;
import org.molgenis.genotype.annotation.SampleAnnotation;
import org.molgenis.genotype.annotation.SampleAnnotation.SampleAnnotationType;
import org.molgenis.genotype.annotation.SexAnnotation;
import org.molgenis.genotype.util.Utils;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variant.NotASnpException;

/**
 * Export a GenotypeData object to an impute2 haps/sample files
 * 
 * Missing values in the samplefile are written as 'NA'
 * 
 * @author erwin
 * 
 */
public class Impute2GenotypeWriter implements GenotypeWriter
{
	public static final Charset FILE_ENCODING = Charset.forName("UTF-8");
	public static final String LINE_ENDING = "\n";
	private static final char SEPARATOR = ' ';
	private static final char UNPHASED_INDICATOR = '*';
	private static final char MISSING_INDICATOR = '?';
	private static final Logger LOG = Logger.getLogger(Impute2GenotypeWriter.class);
	private GenotypeData genotypeData;

	public Impute2GenotypeWriter(GenotypeData genotypeData)
	{
		this.genotypeData = genotypeData;
	}

    @Override
	public void write(String basePath) throws IOException
	{
		write(new File(basePath + ".haps"), new File(basePath + ".sample"));
	}

	public void write(File hapsFile, File sampleFile) throws IOException
	{
		LOG.info("Writing haps file [" + hapsFile.getAbsolutePath() + "] and sample file ["
				+ sampleFile.getAbsolutePath() + "]");

		Utils.createEmptyFile(hapsFile, "haps");
		Utils.createEmptyFile(sampleFile, "sample");
		
		Writer hapsFileWriter = null;
		Writer sampleFileWriter = null;
		try
		{
			hapsFileWriter = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(hapsFile), FILE_ENCODING));
			sampleFileWriter = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(sampleFile),
					FILE_ENCODING));
			write(hapsFileWriter, sampleFileWriter);
		}
		finally
		{
			IOUtils.closeQuietly(hapsFileWriter);
			IOUtils.closeQuietly(sampleFileWriter);
		}
	}

	public void write(Writer hapsFileWriter, Writer sampleFileWriter) throws IOException
	{
		writeHapsFile(hapsFileWriter);
		writeSampleFile(sampleFileWriter);
	}

	private void writeSampleFile(Writer sampleWriter) throws IOException
	{
		// Write headers
		StringBuilder sb = new StringBuilder();
		sb.append("ID_1");
		sb.append(SEPARATOR);
		sb.append("ID_2");
		sb.append(SEPARATOR);
		sb.append("missing");

		List<String> colNames = new ArrayList<String>();
		List<String> dataTypes = new ArrayList<String>();

		for (SampleAnnotation annotation : genotypeData.getSampleAnnotations())
		{
			if(annotation.getId().equals(GenotypeData.SAMPLE_MISSING_RATE_DOUBLE)){
				continue;
			}
			
			
			if (annotation.getSampleAnnotationType() == SampleAnnotationType.COVARIATE)
			{
				switch (annotation.getType())
				{
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
			}
			else if (annotation.getSampleAnnotationType() == SampleAnnotationType.PHENOTYPE)
			{
				switch (annotation.getType())
				{
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
			}
			else
			{
				LOG.warn("'OTHER' sample annotation type not supported by impute2");
			}
		}

		for (String colName : colNames)
		{
			sb.append(SEPARATOR);
			sb.append(colName);
		}

		sb.append(LINE_ENDING);

		// Write datatypes
		sb.append("0");
		sb.append(SEPARATOR);
		sb.append("0");
		sb.append(SEPARATOR);
		sb.append("0");

		for (String dataType : dataTypes)
		{
			sb.append(SEPARATOR);
			sb.append(dataType);
		}

		sb.append(LINE_ENDING);
		sampleWriter.write(sb.toString());

		// Write values
		for (Sample sample : genotypeData.getSamples())
		{
			sb = new StringBuilder();
			sb.append(sample.getFamilyId() == null ? "NA" : sample.getFamilyId());
			sb.append(SEPARATOR);
			sb.append(sample.getId() == null ? "NA" : sample.getId());
			sb.append(SEPARATOR);
			sb.append(getValue(GenotypeData.SAMPLE_MISSING_RATE_DOUBLE, sample, "NA"));

			for (String colName : colNames)
			{
				sb.append(SEPARATOR);
				sb.append(getValue(colName, sample, "NA"));
			}

			sb.append(LINE_ENDING);
			sampleWriter.write(sb.toString());
		}
	}

	private String getValue(String colName, Sample sample, String nullValue)
	{
		Object value = sample.getAnnotationValues().get(colName);
		if (value == null)
		{
			return nullValue;
		}

		if (value instanceof Boolean)
		{
			return value.equals(true) ? "1" : "0";
		}

		if (value instanceof Double || value instanceof Float)
		{
			String result = value.toString();
			if (result.equals("0.0"))
			{
				result = "0";
			}

			return result;
		}
		
		if (value instanceof SexAnnotation){
			return Byte.toString(((SexAnnotation) value).getPlinkSex());
		}

		return value.toString();
	}

	private void writeHapsFile(Writer hapsFileWriter) throws IOException
	{
		for (GeneticVariant variant : genotypeData)
		{
			if (!variant.isSnp())
			{
				throw new NotASnpException(variant);
			}

			Allele allele0 = variant.getVariantAlleles().get(0);
			Allele allele1 = variant.getVariantAlleles().get(1);

			StringBuilder sb = new StringBuilder();
			sb.append(variant.getSequenceName());
			sb.append(SEPARATOR);
			sb.append(variant.getPrimaryVariantId());
			sb.append(SEPARATOR);
			sb.append(variant.getStartPos());
			sb.append(SEPARATOR);
			sb.append(allele0);
			sb.append(SEPARATOR);
			sb.append(allele1);

			List<Alleles> sampleAlleles = variant.getSampleVariants();
			List<Boolean> phasing = variant.getSamplePhasing();

			if ((sampleAlleles != null) && !sampleAlleles.isEmpty())
			{
				for (int i = 0; i < sampleAlleles.size(); i++)
				{
					sb.append(SEPARATOR);

					Alleles alleles = sampleAlleles.get(i);
					if ((alleles == null) || alleles.getAllelesAsString().isEmpty() || (alleles.get(0) == Allele.ZERO)
							|| (alleles.get(1) == Allele.ZERO))
					{
						sb.append(MISSING_INDICATOR);
						sb.append(SEPARATOR);
						sb.append(MISSING_INDICATOR);
					}
					else
					{
						if (alleles.get(0).equals(allele0))
						{
							sb.append("0");
						}
						else if (alleles.get(0).equals(allele1))
						{
							sb.append("1");
						}
						else
						{
							throw new RuntimeException("SampleAllele [" + alleles.get(0) + "] for SNP ["
									+ variant.getPrimaryVariantId() + "] does not match one of the variant alleles");
						}

						if (!phasing.get(i))
						{
							sb.append(UNPHASED_INDICATOR);
						}

						sb.append(SEPARATOR);

						if (alleles.get(1).equals(allele0))
						{
							sb.append("0");
						}
						else if (alleles.get(1).equals(allele1))
						{
							sb.append("1");
						}
						else
						{
							throw new RuntimeException("SampleAllele [" + alleles.get(1) + "] for SNP ["
									+ variant.getPrimaryVariantId() + "] does not match one of the variant alleles");
						}

						if (!phasing.get(i))
						{
							sb.append(UNPHASED_INDICATOR);
						}
					}
				}
			}

			sb.append(LINE_ENDING);

			hapsFileWriter.write(sb.toString());
		}
	}
}
