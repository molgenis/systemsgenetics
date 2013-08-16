package org.molgenis.genotype.plink;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import org.apache.commons.io.IOUtils;
import org.apache.log4j.Logger;
import org.molgenis.genotype.AbstractRandomAccessGenotypeData;
import org.molgenis.genotype.Allele;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.Sample;
import org.molgenis.genotype.Sequence;
import org.molgenis.genotype.SimpleSequence;
import org.molgenis.genotype.annotation.Annotation;
import org.molgenis.genotype.annotation.SampleAnnotation;
import org.molgenis.genotype.annotation.SexAnnotation;
import org.molgenis.genotype.plink.datatypes.MapEntry;
import org.molgenis.genotype.plink.datatypes.PedEntry;
import org.molgenis.genotype.plink.drivers.PedFileDriver;
import org.molgenis.genotype.plink.readers.MapFileReader;
import org.molgenis.genotype.util.Cache;
import org.molgenis.genotype.util.CalledDosageConvertor;
import org.molgenis.genotype.util.GeneticVariantTreeSet;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variant.ReadOnlyGeneticVariant;
import org.molgenis.genotype.variant.sampleProvider.SampleVariantUniqueIdProvider;
import org.molgenis.genotype.variant.sampleProvider.SampleVariantsProvider;

public class PedMapGenotypeData extends AbstractRandomAccessGenotypeData implements SampleVariantsProvider
{
	private static final Logger LOG = Logger.getLogger(PedMapGenotypeData.class);
	private final int sampleVariantProviderUniqueId;

	private final File pedFile;
	private Map<Integer, List<Alleles>> sampleAllelesBySnpIndex = new HashMap<Integer, List<Alleles>>();

	private GeneticVariantTreeSet<GeneticVariant> snps = new GeneticVariantTreeSet<GeneticVariant>();
	private Map<String, Integer> snpIndexById = new HashMap<String, Integer>(1000);
	private Map<String, List<GeneticVariant>> snpBySequence = new TreeMap<String, List<GeneticVariant>>();
	private Map<String, SampleAnnotation> sampleAnnotations;

	private final Cache<GeneticVariant, byte[]> calledDosageCache;
	private final Cache<GeneticVariant, float[]> dosageCache;

	public PedMapGenotypeData(String basePath) throws FileNotFoundException, IOException
	{
		this(new File(basePath + ".ped"), new File(basePath + ".map"));
	}

	public PedMapGenotypeData(File pedFile, File mapFile) throws FileNotFoundException, IOException
	{
		if (pedFile == null) throw new IllegalArgumentException("PedFile is null");
		if (mapFile == null) throw new IllegalArgumentException("MapFile is null");
		if (!mapFile.isFile()) throw new FileNotFoundException("MAP index file not found at "
				+ mapFile.getAbsolutePath());
		if (!mapFile.canRead()) throw new IOException("MAP index file not found at " + mapFile.getAbsolutePath());
		if (!pedFile.isFile()) throw new FileNotFoundException("PED file not found at " + pedFile.getAbsolutePath());
		if (!pedFile.canRead()) throw new IOException("PED file not found at " + pedFile.getAbsolutePath());

		this.pedFile = pedFile;

		MapFileReader mapFileReader = null;
		PedFileDriver pedFileDriver = null;
		try
		{
			pedFileDriver = new PedFileDriver(pedFile);
			loadSampleBialleles(pedFileDriver);

			mapFileReader = new MapFileReader(new FileInputStream(mapFile));
			loadSnps(mapFileReader);

		}
		finally
		{
			IOUtils.closeQuietly(pedFileDriver);
			IOUtils.closeQuietly(mapFileReader);
		}

		sampleVariantProviderUniqueId = SampleVariantUniqueIdProvider.getNextUniqueId();

		sampleAnnotations = PlinkSampleAnnotations.getSampleAnnotations();

		this.calledDosageCache = new Cache<GeneticVariant, byte[]>(100);
		this.dosageCache = new Cache<GeneticVariant, float[]>(100);

	}

	private void loadSampleBialleles(PedFileDriver pedFileDriver)
	{
		int count = 0;
		for (PedEntry entry : pedFileDriver)
		{
			int index = 0;
			for (Alleles biallele : entry)
			{
				List<Alleles> biallelesForSnp = sampleAllelesBySnpIndex.get(index);
				if (biallelesForSnp == null)
				{
					biallelesForSnp = new ArrayList<Alleles>();
					sampleAllelesBySnpIndex.put(index, biallelesForSnp);
				}

				biallelesForSnp.add(biallele);
				index++;
			}

			++count;
			if ((count % 100) == 0)
			{
				LOG.info("Loaded [" + (count) + "] samples");
			}
		}

		LOG.info("Total [" + count + "] samples");
	}

	private void loadSnps(MapFileReader reader)
	{
		int index = 0;
		for (MapEntry entry : reader)
		{
			String id = entry.getSNP();
			String sequenceName = entry.getChromosome();
			int startPos = (int) entry.getBpPos();

			List<Alleles> sampleAlleles = sampleAllelesBySnpIndex.get(index);
			List<String> alleles = new ArrayList<String>(2);
			for (Alleles biallele : sampleAlleles)
			{

				String allele1 = biallele.get(0) == Allele.ZERO ? null : biallele.get(0).toString();
				if ((allele1 != null) && !alleles.contains(allele1))
				{
					alleles.add(allele1);
				}

				String allele2 = biallele.get(1) == Allele.ZERO ? null : biallele.get(1).toString();
				if ((allele2 != null) && !alleles.contains(allele2))
				{
					alleles.add(allele2);
				}
			}

			GeneticVariant snp = ReadOnlyGeneticVariant.createVariant(id, startPos, sequenceName, this, alleles);

			snps.add(snp);
			snpIndexById.put(snp.getPrimaryVariantId(), index);

			List<GeneticVariant> seqGeneticVariants = snpBySequence.get(sequenceName);
			if (seqGeneticVariants == null)
			{
				seqGeneticVariants = new ArrayList<GeneticVariant>();
				snpBySequence.put(sequenceName, seqGeneticVariants);
			}
			seqGeneticVariants.add(snp);

			index++;

			if ((index % 100000) == 0)
			{
				LOG.info("Loaded [" + index + "] snps");
			}
		}

		LOG.info("Total [" + index + "] snps");
	}

	@Override
	public List<Sequence> getSequences()
	{
		List<String> seqNames = getSeqNames();

		List<Sequence> sequences = new ArrayList<Sequence>(seqNames.size());
		for (String seqName : seqNames)
		{
			sequences.add(new SimpleSequence(seqName, null, this));
		}

		return sequences;
	}

	@Override
	public List<Sample> getSamples()
	{
		PedFileDriver pedFileDriver = null;

		try
		{
			pedFileDriver = new PedFileDriver(pedFile);
			List<Sample> samples = new ArrayList<Sample>();
			for (PedEntry pedEntry : pedFileDriver)
			{
				Map<String, Object> annotationValues = new LinkedHashMap<String, Object>();
				annotationValues.put(FATHER_SAMPLE_ANNOTATION_NAME, pedEntry.getFather());
				annotationValues.put(MOTHER_SAMPLE_ANNOTATION_NAME, pedEntry.getMother());
				annotationValues.put(SEX_SAMPLE_ANNOTATION_NAME,
						SexAnnotation.getSexAnnotationForPlink(pedEntry.getSex()));
				annotationValues.put(DOUBLE_PHENOTYPE_SAMPLE_ANNOTATION_NAME, pedEntry.getPhenotype());

				samples.add(new Sample(pedEntry.getIndividual(), pedEntry.getFamily(), annotationValues));
			}

			return samples;
		}
		finally
		{
			IOUtils.closeQuietly(pedFileDriver);
		}
	}

	@Override
	public List<Alleles> getSampleVariants(GeneticVariant variant)

	{
		if (variant.getPrimaryVariantId() == null)
		{
			throw new IllegalArgumentException("Not a snp, missing primaryVariantId");
		}

		Integer index = snpIndexById.get(variant.getPrimaryVariantId());

		if (index == null)
		{
			throw new IllegalArgumentException("Unknown primaryVariantId [" + variant.getPrimaryVariantId() + "]");
		}

		List<Alleles> bialleles = sampleAllelesBySnpIndex.get(index);

		return Collections.unmodifiableList(bialleles);
	}

	@Override
	protected Map<String, Annotation> getVariantAnnotationsMap()
	{
		return Collections.emptyMap();
	}

	@Override
	public List<String> getSeqNames()
	{
		return new ArrayList<String>(snpBySequence.keySet());
	}

	@Override
	public List<GeneticVariant> getVariantsByPos(String seqName, int startPos)
	{
		List<GeneticVariant> variants = new ArrayList<GeneticVariant>(snps.getSequencePosVariants(seqName, startPos));

		return variants;
	}

	@Override
	public Iterator<GeneticVariant> iterator()
	{
		return snps.iterator();
	}

	@Override
	public Iterable<GeneticVariant> getSequenceGeneticVariants(String seqName)
	{
		// TODO remove snpBySequence this makes no sense now that we have the
		// treset
		List<GeneticVariant> variants = snpBySequence.get(seqName);
		if (variants == null)
		{
			throw new IllegalArgumentException("Unknown sequence [" + seqName + "]");
		}

		return variants;
	}

	@Override
	public int cacheSize()
	{
		return 0;
	}

	/**
	 * Ped/Map daoes not support phasing, always return false
	 * 
	 */
	@Override
	public List<Boolean> getSamplePhasing(GeneticVariant variant)
	{

		List<Boolean> phasing = Collections.nCopies(getSampleVariants(variant).size(), false);

		return phasing;
	}

	@Override
	public int getSampleVariantProviderUniqueId()
	{
		return sampleVariantProviderUniqueId;
	}

	@Override
	protected Map<String, SampleAnnotation> getSampleAnnotationsMap()
	{
		return sampleAnnotations;
	}

	@Override
	public Iterable<GeneticVariant> getVariantsByRange(String seqName, int rangeStart, int rangeEnd)
	{
		return snps.getSequenceRangeVariants(seqName, rangeStart, rangeEnd);
	}

	@Override
	public byte[] getSampleCalledDosage(GeneticVariant variant)
	{

		if (calledDosageCache.containsKey(variant))
		{
			return calledDosageCache.get(variant);
		}

		byte[] calledDosage = CalledDosageConvertor.convertCalledAllelesToCalledDosage(getSampleVariants(variant),
				variant.getVariantAlleles(), variant.getRefAllele());
		calledDosageCache.put(variant, calledDosage);
		return calledDosage;

	}

	@Override
	public float[] getSampleDosage(GeneticVariant variant)
	{

		if (dosageCache.containsKey(variant))
		{
			return dosageCache.get(variant);
		}

		float[] dosage = CalledDosageConvertor.convertCalledAllelesToDosage(getSampleVariants(variant),
				variant.getVariantAlleles(), variant.getRefAllele());
		dosageCache.put(variant, dosage);
		return dosage;

	}

	@Override
	public void close() throws IOException {
	}
}
