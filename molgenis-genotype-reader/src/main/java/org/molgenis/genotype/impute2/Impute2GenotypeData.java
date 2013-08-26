package org.molgenis.genotype.impute2;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.StringTokenizer;

import org.apache.commons.io.IOUtils;
import org.apache.log4j.Logger;
import org.molgenis.genotype.Allele;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.GenotypeDataException;
import org.molgenis.genotype.GenotypeDataIndex;
import org.molgenis.genotype.IndexedGenotypeData;
import org.molgenis.genotype.RawLineQueryResult;
import org.molgenis.genotype.Sample;
import org.molgenis.genotype.Sequence;
import org.molgenis.genotype.SimpleSequence;
import org.molgenis.genotype.annotation.Annotation;
import org.molgenis.genotype.annotation.SampleAnnotation;
import org.molgenis.genotype.annotation.SampleAnnotation.SampleAnnotationType;
import org.molgenis.genotype.tabix.TabixIndex;
import org.molgenis.genotype.util.CalledDosageConvertor;
import org.molgenis.genotype.util.Utils;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variant.ReadOnlyGeneticVariant;
import org.molgenis.genotype.variant.VariantLineMapper;
import org.molgenis.genotype.variant.sampleProvider.CachedSampleVariantProvider;
import org.molgenis.genotype.variant.sampleProvider.SampleVariantUniqueIdProvider;
import org.molgenis.genotype.variant.sampleProvider.SampleVariantsProvider;
import org.molgenis.io.csv.CsvReader;
import org.molgenis.util.tuple.Tuple;

/**
 * GenotypeData for haps/sample files see http://www.shapeit.fr/
 * 
 * First run
 * <code>index-haps.sh yourfile.haps<code> to create the tabix index file
 * 
 * The two character string 'NA' is treated as missing when encountered in the sample file
 * 
 * @author erwin fixed by patrick
 * 
 */
public class Impute2GenotypeData extends IndexedGenotypeData implements SampleVariantsProvider
{
	private GenotypeDataIndex index;
	private File sampleFile;
	private Map<String, SampleAnnotation> sampleAnnotations = new LinkedHashMap<String, SampleAnnotation>();
	private final int sampleVariantProviderUniqueId;
	private final VariantLineMapper lineMapper;
	private final SampleVariantsProvider sampleVariantProvider;

	private static final Logger LOG = Logger.getLogger(Impute2GenotypeData.class);

	public Impute2GenotypeData(File bzipHapsFile, File tabixIndexFile, File sampleFile) throws IOException
	{
		this(bzipHapsFile, tabixIndexFile, sampleFile, 0);
	}

	public Impute2GenotypeData(File bzipHapsFile, File tabixIndexFile, File sampleFile, int cacheSize)
			throws IOException
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

		if (cacheSize > 0)
		{
			sampleVariantProvider = new CachedSampleVariantProvider(this, cacheSize);
		}
		else
		{
			sampleVariantProvider = this;
		}
		lineMapper = new Impute2VariantLineMapper(sampleVariantProvider);

		index = new TabixIndex(tabixIndexFile, bzipHapsFile, lineMapper);
		LOG.info("Read tabix index");

		this.sampleFile = sampleFile;

		loadAnnotations();
		LOG.info("Annotations loaded");

		sampleVariantProviderUniqueId = SampleVariantUniqueIdProvider.getNextUniqueId();

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
	public Map<String, Annotation> getVariantAnnotationsMap()
	{
		return Collections.emptyMap();
	}

	@Override
	public Map<String, SampleAnnotation> getSampleAnnotationsMap()
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

	@Override
	public List<Alleles> getSampleVariants(GeneticVariant variant)
	{
		RawLineQueryResult queryResult = index.createRawLineQuery().executeQuery(variant.getSequenceName(),
				variant.getStartPos());

		ArrayList<Alleles> genotypes = new ArrayList<Alleles>();
		try
		{
			for (String line : queryResult)
			{

				StringTokenizer tokenizer = new StringTokenizer(line, "\t");
				String chrom = tokenizer.nextToken();
				String snpId = tokenizer.nextToken();
				int position = Integer.parseInt(tokenizer.nextToken());
				String firstAllele = tokenizer.nextToken();
				String secondAllele = tokenizer.nextToken();

				List<String> alleles = Arrays.asList(firstAllele, secondAllele);

				GeneticVariant variantInFile = ReadOnlyGeneticVariant.createVariant(snpId, position, chrom,
						sampleVariantProvider, alleles);

				if (variantInFile.equals(variant)
						&& variantInFile.getPrimaryVariantId().equals(variant.getPrimaryVariantId()))
				{
					// same variant now load sample alleles
					while (tokenizer.hasMoreTokens())
					{

						Allele allele1 = createAllele(tokenizer.nextToken(), variant);
						Allele allele2 = createAllele(tokenizer.nextToken(), variant);
						genotypes.add(Alleles.createAlleles(allele1, allele2));

					}

				}

			}
		}
		finally
		{
			IOUtils.closeQuietly(queryResult);
		}

		return Collections.unmodifiableList(genotypes);
	}

	@Override
	public List<Boolean> getSamplePhasing(GeneticVariant variant)
	{
		RawLineQueryResult queryResult = index.createRawLineQuery().executeQuery(variant.getSequenceName(),
				variant.getStartPos());

		ArrayList<Boolean> phasing = new ArrayList<Boolean>();
		try
		{
			for (String line : queryResult)
			{

				StringTokenizer tokenizer = new StringTokenizer(line, "\t");
				String chrom = tokenizer.nextToken();
				String snpId = tokenizer.nextToken();
				int position = Integer.parseInt(tokenizer.nextToken());
				String firstAllele = tokenizer.nextToken();
				String secondAllele = tokenizer.nextToken();

				List<String> alleles = Arrays.asList(firstAllele, secondAllele);

				GeneticVariant variantInFile = ReadOnlyGeneticVariant.createVariant(snpId, position, chrom,
						sampleVariantProvider, alleles);

				if (variantInFile.equals(variant))
				{
					// same variant now load sample alleles
					while (tokenizer.hasMoreTokens())
					{

						boolean phased = !tokenizer.nextToken().endsWith("*");
						phasing.add(phased);

						// only use first allele of sample to determine phasing
						tokenizer.nextToken();

					}

				}

			}
		}
		finally
		{
			IOUtils.closeQuietly(queryResult);
		}

		return Collections.unmodifiableList(phasing);

	}

	@Override
	public int cacheSize()
	{
		return 0;
	}

	@Override
	public int getSampleVariantProviderUniqueId()
	{
		return sampleVariantProviderUniqueId;
	}

	@Override
	public byte[] getSampleCalledDosage(GeneticVariant variant)
	{
		return CalledDosageConvertor.convertCalledAllelesToCalledDosage(getSampleVariants(variant),
				variant.getVariantAlleles(), variant.getRefAllele());
	}

	@Override
	public float[] getSampleDosage(GeneticVariant variant)
	{
		return CalledDosageConvertor.convertCalledAllelesToDosage(getSampleVariants(variant),
				variant.getVariantAlleles(), variant.getRefAllele());
	}

	/**
	 * Sample can have an asterisk directly after the allele indicating that it
	 * is unphased
	 */
	private Allele createAllele(String sample, GeneticVariant variant)
	{
		char allele = sample.charAt(0);
		switch (allele)
		{
			case '?':
				return Allele.ZERO;
			case '0':
				return variant.getVariantAlleles().get(0);
			case '1':
				return variant.getVariantAlleles().get(1);
			default:
				throw new GenotypeDataException("[" + sample + "] is an invalid value for a haps sample value");
		}

	}

	@Override
	public void close() throws IOException {

	}
}