package org.molgenis.genotype.plink;

import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.io.IOUtils;
import org.apache.log4j.Logger;
import org.molgenis.genotype.*;
import org.molgenis.genotype.annotation.Annotation;
import org.molgenis.genotype.annotation.SampleAnnotation;
import org.molgenis.genotype.annotation.SexAnnotation;
import org.molgenis.genotype.plink.datatypes.MapEntry;
import org.molgenis.genotype.plink.datatypes.PedEntry;
import org.molgenis.genotype.plink.drivers.PedFileDriver;
import org.molgenis.genotype.plink.readers.MapFileReader;
import org.molgenis.genotype.util.Cache;
import org.molgenis.genotype.util.CalledDosageConvertor;
import org.molgenis.genotype.util.FixedSizeIterable;
import org.molgenis.genotype.util.ProbabilitiesConvertor;
import org.molgenis.genotype.util.RecordIteratorCreators;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variant.GeneticVariantMeta;
import org.molgenis.genotype.variant.GeneticVariantMetaMap;
import org.molgenis.genotype.variant.GenotypeRecord;
import org.molgenis.genotype.variant.ReadOnlyGeneticVariant;
import org.molgenis.genotype.variant.range.GeneticVariantRange;
import org.molgenis.genotype.variant.sampleProvider.SampleVariantUniqueIdProvider;
import org.molgenis.genotype.variant.sampleProvider.SampleVariantsProvider;

public class PedMapGenotypeData extends AbstractRandomAccessGenotypeData implements SampleVariantsProvider {

	private static final Logger LOG = Logger.getLogger(PedMapGenotypeData.class);
	private final int sampleVariantProviderUniqueId;
	private TIntObjectMap<List<Alleles>> sampleAllelesBySnpIndex = new TIntObjectHashMap<List<Alleles>>(1000, 0.75f);
	private GeneticVariantRange snps;
	private TObjectIntMap<String> snpIndexById = new TObjectIntHashMap<String>(1000, 0.75f, -1);
	private Map<String, SampleAnnotation> sampleAnnotations;
	private final Cache<GeneticVariant, byte[]> calledDosageCache;
	private final Cache<GeneticVariant, float[]> dosageCache;
	private ArrayList<Sample> samples = new ArrayList<Sample>();
	private Set<String> seqNames = new LinkedHashSet<String>();
	private GeneticVariantMeta geneticVariantMeta = GeneticVariantMetaMap.getGeneticVariantMetaGt();

	public PedMapGenotypeData(String basePath) throws FileNotFoundException, IOException {
		this(new File(basePath + ".ped"), new File(basePath + ".map"));
	}

	public PedMapGenotypeData(File pedFile, File mapFile) throws FileNotFoundException, IOException {
		if (pedFile == null) {
			throw new IllegalArgumentException("PedFile is null");
		}
		if (mapFile == null) {
			throw new IllegalArgumentException("MapFile is null");
		}
		if (!mapFile.isFile()) {
			throw new FileNotFoundException("MAP index file not found at "
					+ mapFile.getAbsolutePath());
		}
		if (!mapFile.canRead()) {
			throw new IOException("MAP index file not found at " + mapFile.getAbsolutePath());
		}
		if (!pedFile.isFile()) {
			throw new FileNotFoundException("PED file not found at " + pedFile.getAbsolutePath());
		}
		if (!pedFile.canRead()) {
			throw new IOException("PED file not found at " + pedFile.getAbsolutePath());
		}

		MapFileReader mapFileReader = null;
		PedFileDriver pedFileDriver = null;
		try {
			pedFileDriver = new PedFileDriver(pedFile);
			loadSampleBialleles(pedFileDriver);

			mapFileReader = new MapFileReader(new FileInputStream(mapFile));
			loadSnps(mapFileReader);

		} finally {
			IOUtils.closeQuietly(pedFileDriver);
			IOUtils.closeQuietly(mapFileReader);
		}

		sampleVariantProviderUniqueId = SampleVariantUniqueIdProvider.getNextUniqueId();

		sampleAnnotations = PlinkSampleAnnotations.getSampleAnnotations();

		this.calledDosageCache = new Cache<GeneticVariant, byte[]>(100);
		this.dosageCache = new Cache<GeneticVariant, float[]>(100);

	}

	private void loadSampleBialleles(PedFileDriver pedFileDriver) {
		int count = 0;
		for (PedEntry entry : pedFileDriver) {

			Map<String, Object> annotationValues = new LinkedHashMap<String, Object>();
			annotationValues.put(FATHER_SAMPLE_ANNOTATION_NAME, entry.getFather());
			annotationValues.put(MOTHER_SAMPLE_ANNOTATION_NAME, entry.getMother());
			annotationValues.put(SEX_SAMPLE_ANNOTATION_NAME,
					SexAnnotation.getSexAnnotationForPlink(entry.getSex()));
			annotationValues.put(DOUBLE_PHENOTYPE_SAMPLE_ANNOTATION_NAME, entry.getPhenotype());

			samples.add(new Sample(entry.getIndividual(), entry.getFamily(), annotationValues));

			int index = 0;
			for (Alleles biallele : entry) {
				List<Alleles> biallelesForSnp = sampleAllelesBySnpIndex.get(index);
				if (biallelesForSnp == null) {
					biallelesForSnp = new ArrayList<Alleles>();
					sampleAllelesBySnpIndex.put(index, biallelesForSnp);
				}

				biallelesForSnp.add(biallele);
				index++;
			}

			++count;
			if ((count % 100) == 0) {
				LOG.info("Loaded [" + (count) + "] samples");
			}
		}

		LOG.info("Total [" + count + "] samples");
	}

	private void loadSnps(MapFileReader reader) {
		int index = 0;
		
		GeneticVariantRange.GeneticVariantRangeCreate rangeFactory = GeneticVariantRange.createRangeFactory();
		
		for (MapEntry entry : reader) {
			String id = entry.getSNP();
			String sequenceName = entry.getChromosome();

			seqNames.add(sequenceName);

			int startPos = (int) entry.getBpPos();

			List<Alleles> sampleAlleles = sampleAllelesBySnpIndex.get(index);
			List<String> alleles = new ArrayList<String>(2);
			for (Alleles biallele : sampleAlleles) {

				String allele1 = biallele.get(0) == Allele.ZERO ? null : biallele.get(0).toString();
				if ((allele1 != null) && !alleles.contains(allele1)) {
					alleles.add(allele1);
				}

				String allele2 = biallele.get(1) == Allele.ZERO ? null : biallele.get(1).toString();
				if ((allele2 != null) && !alleles.contains(allele2)) {
					alleles.add(allele2);
				}
			}

			GeneticVariant snp = ReadOnlyGeneticVariant.createVariant(geneticVariantMeta, id, startPos, sequenceName, this, alleles);

			rangeFactory.addVariant(snp);
			snpIndexById.put(snp.getPrimaryVariantId(), index);

			index++;

			if ((index % 100000) == 0) {
				LOG.info("Loaded [" + index + "] snps");
			}
		}
		snps = rangeFactory.createRange();
		LOG.info("Total [" + index + "] snps");
	}

	@Override
	public List<Sequence> getSequences() {
		List<Sequence> sequences = new ArrayList<Sequence>(seqNames.size());
		for (String seqName : seqNames) {
			sequences.add(new SimpleSequence(seqName, null, this));
		}

		return sequences;
	}

	@Override
	public List<Sample> getSamples() {
		return samples;
	}

	@Override
	public List<Alleles> getSampleVariants(GeneticVariant variant) {
		if (variant.getPrimaryVariantId() == null) {
			throw new IllegalArgumentException("Not a snp, missing primaryVariantId");
		}

		int index = snpIndexById.get(variant.getPrimaryVariantId());

		if (index == -1) {
			throw new IllegalArgumentException("Unknown primaryVariantId [" + variant.getPrimaryVariantId() + "]");
		}

		List<Alleles> bialleles = sampleAllelesBySnpIndex.get(index);

		return Collections.unmodifiableList(bialleles);
	}

	@Override
	public Map<String, Annotation> getVariantAnnotationsMap() {
		return Collections.emptyMap();
	}

	@Override
	public List<String> getSeqNames() {
		return new ArrayList<String>(seqNames);
	}

	@Override
	public Iterable<GeneticVariant> getVariantsByPos(String seqName, int startPos) {
		return snps.getVariantAtPos(seqName, startPos);
	}

	@Override
	public Iterator<GeneticVariant> iterator() {
		return snps.iterator();
	}

	@Override
	public Iterable<GeneticVariant> getSequenceGeneticVariants(String seqName) {
		return snps.getVariantsBySequence(seqName);
	}

	@Override
	public int cacheSize() {
		return 0;
	}

	/**
	 * Ped/Map daoes not support phasing, always return false
	 *
	 */
	@Override
	public List<Boolean> getSamplePhasing(GeneticVariant variant) {

		List<Boolean> phasing = Collections.nCopies(getSampleVariants(variant).size(), false);

		return phasing;
	}

	@Override
	public int getSampleVariantProviderUniqueId() {
		return sampleVariantProviderUniqueId;
	}

	@Override
	public Map<String, SampleAnnotation> getSampleAnnotationsMap() {
		return sampleAnnotations;
	}

	@Override
	public Iterable<GeneticVariant> getVariantsByRange(String seqName, int rangeStart, int rangeEnd) {
		return snps.getVariantsByRange(seqName, rangeStart, rangeEnd);
	}

	@Override
	public byte[] getSampleCalledDosage(GeneticVariant variant) {

		if (calledDosageCache.containsKey(variant)) {
			return calledDosageCache.get(variant);
		}

		byte[] calledDosage = CalledDosageConvertor.convertCalledAllelesToCalledDosage(getSampleVariants(variant),
				variant.getVariantAlleles(), variant.getRefAllele());
		calledDosageCache.put(variant, calledDosage);
		return calledDosage;

	}

	@Override
	public float[] getSampleDosage(GeneticVariant variant) {

		if (dosageCache.containsKey(variant)) {
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

	@Override
	public boolean isOnlyContaingSaveProbabilityGenotypes() {
		return true;
	}

	@Override
	public float[][] getSampleProbilities(GeneticVariant variant) {
		return ProbabilitiesConvertor.convertCalledAllelesToProbability(variant.getSampleVariants(), variant.getVariantAlleles());
	}

	@Override
	public double[][] getSampleProbabilitiesComplex(GeneticVariant variant) {
		return ProbabilitiesConvertor.convertProbabilitiesToComplexProbabilities(getSampleProbilities(variant));
	}

	@Override
	public double[][][] getSampleProbabilitiesPhased(GeneticVariant variant) {
		throw new GenotypeDataException("Phased data not available");
	}

	@Override
	public FixedSizeIterable<GenotypeRecord> getSampleGenotypeRecords(GeneticVariant variant) {
		
		return RecordIteratorCreators.createIteratorFromAlleles(variant.getSampleVariants());
		
	}
}
