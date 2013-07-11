package org.molgenis.genotype.impute2;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.io.IOUtils;
import org.apache.log4j.Logger;
import org.molgenis.genotype.GenotypeDataIndex;
import org.molgenis.genotype.IndexedGenotypeData;
import org.molgenis.genotype.Sample;
import org.molgenis.genotype.Sequence;
import org.molgenis.genotype.SimpleSequence;
import org.molgenis.genotype.annotation.Annotation;
import org.molgenis.genotype.annotation.SampleAnnotation;
import org.molgenis.genotype.annotation.SampleAnnotation.SampleAnnotationType;
import org.molgenis.genotype.tabix.TabixIndex;
import org.molgenis.genotype.util.Utils;
import org.molgenis.io.csv.CsvReader;
import org.molgenis.util.tuple.Tuple;

/**
 * GenotypeData for haps/sample files see http://www.shapeit.fr/
 * 
 * First run <code>index-haps.sh yourfile.haps<code> to create the tabix index file
 * 
 * The two character string 'NA' is treated as missing when encountered in the sample file
 * 
 * @author erwin
 * 
 */
public class Impute2GenotypeData extends IndexedGenotypeData
{
	private GenotypeDataIndex index;
	private File sampleFile;
	private Map<String, SampleAnnotation> sampleAnnotations = new LinkedHashMap<String, SampleAnnotation>();

	private static final Logger LOG = Logger.getLogger(Impute2GenotypeData.class);

	public Impute2GenotypeData(File bzipHapsFile, File tabixIndexFile, File sampleFile) throws IOException
	{
		if (bzipHapsFile == null) throw new IllegalArgumentException("bzipHapsFile is null");
		if (!bzipHapsFile.isFile()) throw new FileNotFoundException("bzipHapsFile file not found at "
				+ bzipHapsFile.getAbsolutePath());
		if (!bzipHapsFile.canRead()) throw new IOException("bzipHapsFile file not found at "
				+ bzipHapsFile.getAbsolutePath());

		if (tabixIndexFile == null) throw new IllegalArgumentException("tabixIndexFile is null");
		if (!tabixIndexFile.isFile()) throw new FileNotFoundException("tabixIndexFile file not found at "
				+ tabixIndexFile.getAbsolutePath());
		if (!tabixIndexFile.canRead()) throw new IOException("tabixIndexFile file not found at "
				+ tabixIndexFile.getAbsolutePath());

		if (sampleFile == null) throw new IllegalArgumentException("sampleFile is null");
		if (!sampleFile.isFile()) throw new FileNotFoundException("sampleFile file not found at "
				+ sampleFile.getAbsolutePath());
		if (!sampleFile.canRead()) throw new IOException("sampleFile file not found at " + sampleFile.getAbsolutePath());

		index = new TabixIndex(tabixIndexFile, bzipHapsFile, new Impute2VariantLineMapper());
		LOG.info("Read tabix index");

		this.sampleFile = sampleFile;

		loadAnnotations();
		LOG.info("Annotations loaded");
	}

	@Override
	public Iterable<Sequence> getSequences()
	{
		List<Sequence> sequences = new ArrayList<Sequence>();
		for (String seqName : getSeqNames())
		{
			sequences.add(new SimpleSequence(seqName, null, this));
		}

		return sequences;
	}

	@Override
	public List<Sample> getSamples()
	{
		List<Sample> samples = new ArrayList<Sample>();

		CsvReader reader = null;
		try
		{
			reader = new CsvReader(sampleFile, ' ');
			Iterator<Tuple> it = reader.iterator();

			it.next();// Datatype row

			while (it.hasNext())
			{
				Tuple tuple = it.next();

				String familyId = tuple.getString(0);
				String sampleId = tuple.getString(1);

				Map<String, Object> annotationValues = new LinkedHashMap<String, Object>();
				annotationValues.put("missing", tuple.getDouble(2));

				for (String colName : sampleAnnotations.keySet())
				{
					SampleAnnotation annotation = sampleAnnotations.get(colName);

					Object value = null;
					if (!tuple.getString(colName).equalsIgnoreCase("NA"))
					{
						switch (annotation.getType())
						{
							case INTEGER:
								value = tuple.getInt(colName);
								break;
							case BOOLEAN:
								value = tuple.getBoolean(colName);
								break;
							case FLOAT:
								value = tuple.getDouble(colName);
								break;
							default:
								LOG.warn("Unsupported data type encountered for column [" + colName + "]");
						}
					}

					annotationValues.put(colName, value);
				}

				samples.add(new Sample(sampleId, familyId, annotationValues));
			}

		}
		catch (FileNotFoundException e)
		{
			throw new RuntimeException("File [" + sampleFile.getAbsolutePath() + "] does not exists", e);
		}
		finally
		{
			IOUtils.closeQuietly(reader);
		}
		return samples;
	}

	@Override
	protected GenotypeDataIndex getIndex()
	{
		return index;
	}

	@Override
	protected Map<String, Annotation> getVariantAnnotationsMap()
	{
		return Collections.emptyMap();
	}

	@Override
	protected Map<String, SampleAnnotation> getSampleAnnotationsMap()
	{
		return sampleAnnotations;
	}

	private void loadAnnotations() throws IOException
	{
		sampleAnnotations.clear();

		SampleAnnotation missingAnnotation = new SampleAnnotation("missing", "missing",
				"Missing data proportion of each individual", Annotation.Type.FLOAT, SampleAnnotationType.OTHER, false);
		sampleAnnotations.put(missingAnnotation.getId(), missingAnnotation);

		CsvReader reader = null;
		try
		{
			reader = new CsvReader(sampleFile, ' ');

			List<String> colNames = Utils.iteratorToList(reader.colNamesIterator());
			Tuple dataTypes = reader.iterator().next();
			for (int i = 3; i < colNames.size(); i++)
			{
				SampleAnnotation annotation = null;
				if (dataTypes.getString(i).equalsIgnoreCase("D"))
				{
					annotation = new SampleAnnotation(colNames.get(i), colNames.get(i), "", Annotation.Type.INTEGER,
							SampleAnnotationType.COVARIATE, false);

				}
				else if (dataTypes.getString(i).equalsIgnoreCase("C"))
				{
					annotation = new SampleAnnotation(colNames.get(i), colNames.get(i), "", Annotation.Type.FLOAT,
							SampleAnnotationType.COVARIATE, false);
				}
				else if (dataTypes.getString(i).equalsIgnoreCase("P"))
				{
					annotation = new SampleAnnotation(colNames.get(i), colNames.get(i), "", Annotation.Type.FLOAT,
							SampleAnnotationType.PHENOTYPE, false);
				}
				else if (dataTypes.getString(i).equalsIgnoreCase("B"))
				{
					annotation = new SampleAnnotation(colNames.get(i), colNames.get(i), "", Annotation.Type.BOOLEAN,
							SampleAnnotationType.PHENOTYPE, false);
				}
				else
				{
					LOG.warn("Unknown datatype [" + dataTypes.getString(i) + "]");
				}

				if (annotation != null)
				{
					sampleAnnotations.put(annotation.getId(), annotation);
				}
			}
		}
		finally
		{
			IOUtils.closeQuietly(reader);
		}

	}
}
