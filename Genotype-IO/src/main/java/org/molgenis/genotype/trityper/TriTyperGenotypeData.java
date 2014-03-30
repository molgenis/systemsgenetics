/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.trityper;

import gnu.trove.map.hash.TObjectIntHashMap;
import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import org.apache.log4j.Logger;
import org.molgenis.genotype.AbstractRandomAccessGenotypeData;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.GenotypeData;
import static org.molgenis.genotype.GenotypeData.BOOL_INCLUDE_SAMPLE;
import org.molgenis.genotype.GenotypeDataException;
import org.molgenis.genotype.Sample;
import org.molgenis.genotype.Sequence;
import org.molgenis.genotype.SimpleSequence;
import org.molgenis.genotype.annotation.Annotation;
import org.molgenis.genotype.annotation.CaseControlAnnotation;
import org.molgenis.genotype.annotation.SampleAnnotation;
import org.molgenis.genotype.annotation.SexAnnotation;
import org.molgenis.genotype.sampleFilter.SampleFilter;
import org.molgenis.genotype.sampleFilter.SampleIncludedFilter;
import org.molgenis.genotype.util.CalledDosageConvertor;
import org.molgenis.genotype.util.ProbabilitiesConvertor;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variant.ReadOnlyGeneticVariantTriTyper;
import org.molgenis.genotype.variant.range.GeneticVariantRange;
import org.molgenis.genotype.variant.sampleProvider.CachedSampleVariantProvider;
import org.molgenis.genotype.variant.sampleProvider.SampleVariantUniqueIdProvider;
import org.molgenis.genotype.variant.sampleProvider.SampleVariantsProvider;
import org.molgenis.genotype.variantFilter.VariantFilter;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author harmjan
 */
public class TriTyperGenotypeData extends AbstractRandomAccessGenotypeData implements SampleVariantsProvider {

	private final List<Boolean> samplePhasing;
	private final GeneticVariantRange snps;
	private final SampleVariantsProvider variantProvider;
	private final File genotypeDataFile;
	private final File imputedDosageDataFile;
	private final File snpFile;
	private final File snpMapFile;
	private final File individualFile;
	private final File phenotypeAnnotationFile;
	private static final Logger LOG = Logger.getLogger(TriTyperGenotypeData.class);
	private final int cacheSize;
	private final RandomAccessFile dosageHandle;
	private final RandomAccessFile genotypeHandle;
	private final FileChannel dosageChannel;
	private final int sampleVariantProviderUniqueId;
	private HashMap<String, SampleAnnotation> sampleAnnotationMap;
	private HashMap<String, Sequence> sequences;
	private final VariantFilter variantFilter;
	private final SampleFilter sampleFilter;
	private int unfilteredSnpCount;
	/**
	 * These are the samples as visible to the outside. if sample filter is used
	 * then a subset of all samples in dataset otherwise ref to all samples
	 * arraylist.
	 */
	private ArrayList<Sample> includedSamples;
	/**
	 * These are samples present in the dataset. If sample filters are used then
	 * the it could be that there are fewer samples returned
	 */
	private ArrayList<Sample> samples;

	public TriTyperGenotypeData(String location) throws IOException {
		this(new File(location), 1024, null, null);
	}

	public TriTyperGenotypeData(String location, int cacheSize) throws IOException {
		this(new File(location), cacheSize, null, null);
	}

	public TriTyperGenotypeData(String location, int cacheSize, VariantFilter variantFilter) throws IOException {
		this(new File(location), cacheSize, variantFilter, null);
	}

	public TriTyperGenotypeData(String location, int cacheSize, VariantFilter variantFilter, boolean readOnlyIncludedIndividuals) throws IOException {
		this(new File(location), cacheSize, variantFilter, readOnlyIncludedIndividuals ? new SampleIncludedFilter() : null);
	}

	public TriTyperGenotypeData(File location) throws IOException {
		this(location, 1024, null, null);
	}

	public TriTyperGenotypeData(File location, int cacheSize, VariantFilter variantFilter, boolean readOnlyIncludedIndividuals) throws IOException {
		this(location, cacheSize, variantFilter, readOnlyIncludedIndividuals ? new SampleIncludedFilter() : null);
	}

	public TriTyperGenotypeData(File location, int cacheSize, VariantFilter variantFilter, SampleFilter sampleFilter) throws IOException {
		this(new File(location, "GenotypeMatrix.dat"), new File(location, "ImputedDosageMatrix.dat").exists() ? new File(location, "ImputedDosageMatrix.dat") : null, new File(location, "SNPs.txt.gz").exists() ? new File(location, "SNPs.txt.gz") : new File(location, "SNPs.txt"), new File(location, "SNPMappings.txt.gz").exists() ? new File(location, "SNPMappings.txt.gz") : new File(location, "SNPMappings.txt"), new File(location, "Individuals.txt.gz").exists() ? new File(location, "Individuals.txt.gz") : new File(location, "Individuals.txt"), new File(location, "PhenotypeInformation.txt.gz").exists() ? new File(location, "PhenotypeInformation.txt.gz") : new File(location, "PhenotypeInformation.txt"), cacheSize, variantFilter, sampleFilter);
	}

	public TriTyperGenotypeData(File genotypeDataFile, File imputedDosageDataFile, File snpFile, File snpMapFile, File individualFile, File phenotypeAnnotationFile, int cacheSize, VariantFilter variantFilter, SampleFilter sampleFilter) throws IOException {

		this.variantFilter = variantFilter;
		this.sampleFilter = sampleFilter;
		this.sampleVariantProviderUniqueId = SampleVariantUniqueIdProvider.getNextUniqueId();

		if (cacheSize <= 0) {
			variantProvider = this;
		} else {
			variantProvider = new CachedSampleVariantProvider(this, cacheSize);
		}
		this.cacheSize = cacheSize;

		this.genotypeDataFile = genotypeDataFile;
		if (!genotypeDataFile.exists()) {
			throw new GenotypeDataException("GenotypeMatrix.dat not found at: " + genotypeDataFile.getAbsolutePath());
		}

		this.imputedDosageDataFile = imputedDosageDataFile;
		if (this.imputedDosageDataFile != null && !this.imputedDosageDataFile.exists()) {
			//ofcourse this file is optional but if it is explicitly specified we check if it exsists. This error will not occure during normal operations
			throw new GenotypeDataException("ImputedDosageMatrix.dat not found at:" + this.imputedDosageDataFile.getAbsolutePath());
		}

		this.snpFile = snpFile;
		if (!this.snpFile.exists()) {
			throw new GenotypeDataException("SNPs.txt or SNPs.txt.gz at:" + this.snpFile.getAbsolutePath());
		}

		this.snpMapFile = snpMapFile;
		if (!this.snpMapFile.exists()) {
			throw new GenotypeDataException("SNPMappings.txt or SNPMappings.txt.gz at:" + this.snpMapFile.getAbsolutePath());
		}

		this.individualFile = individualFile;
		if (!this.individualFile.exists()) {
			throw new GenotypeDataException("Individuals.txt or Individuals.txt.gz at:" + this.individualFile.getAbsolutePath());
		}

		this.phenotypeAnnotationFile = phenotypeAnnotationFile;
		if (!this.phenotypeAnnotationFile.exists()) {
			throw new GenotypeDataException("PhenotypeInformation.txt or PhenotypeInformation.txt.gz at:" + this.phenotypeAnnotationFile.getAbsolutePath());
		}

		// create file handles //
		if (imputedDosageDataFile != null) {
			dosageHandle = new RandomAccessFile(imputedDosageDataFile, "r");
			dosageChannel = dosageHandle.getChannel();
		} else {
			dosageHandle = null;
			dosageChannel = null;
		}

		genotypeHandle = new RandomAccessFile(genotypeDataFile, "r");

		loadSamples();
		samplePhasing = Collections.nCopies(includedSamples.size(), false);

		GeneticVariantRange.ClassGeneticVariantRangeCreate snpsFactory = GeneticVariantRange.createRangeFactory();
		loadSNPAnnotation(snpsFactory);
		snps = snpsFactory.createRange();

		checkFileSize();



	}

	@Override
	public GeneticVariant getSnpVariantByPos(String seqName, int startPos) {
		Iterable<GeneticVariant> variants = getVariantsByPos(seqName, startPos);
		for (GeneticVariant variant : variants) {
			return variant; // this is a speed hack for TriTyper data: all variants in TriTyper data are SNPs or equivalents
		}
		return null;
	}

	private void loadSamples() throws IOException {
		// load sample list
		int i = 0;
		TextFile t = new TextFile(individualFile, TextFile.R);
		String[] lineElems = t.readLineElemsReturnReference(TextFile.tab);
		samples = new ArrayList<Sample>();
		HashMap<String, Sample> sampleNameToSampleObj = new HashMap<String, Sample>();
		while (lineElems != null) {
			String individual = new String(lineElems[0].getBytes("UTF-8"));
			Sample sample = new Sample(individual, null, null);

			sampleNameToSampleObj.put(individual, sample);
			samples.add(sample);
			i++;
			lineElems = t.readLineElemsReturnReference(TextFile.tab);
		}
		t.close();

		// load annotation
		t = new TextFile(phenotypeAnnotationFile, TextFile.R);

		int numIncluded = 0;
		int numFemale = 0;
		int numMale = 0;
		int numUnknownSex = 0;
		int numCase = 0;
		int numControl = 0;
		int numUnknown = 0;
		int numAnnotated = 0;

		lineElems = t.readLineElemsReturnReference(TextFile.tab);

		HashSet<Sample> visitedSamples = new HashSet<Sample>();
		while (lineElems != null) {
			String individual = lineElems[0];
			Sample sampleObj = sampleNameToSampleObj.get(individual);
			if (sampleObj != null) {
				if (visitedSamples.contains(sampleObj)) {
					LOG.warn("Sample " + sampleObj.getId() + " may have duplicate annotation in PhenotypeInformation.txt.");
				} else {
					Boolean includeSample = parseIncludeExcludeStatus(lineElems[2]);
					if (!includeSample) {
						numIncluded++;
					}
					CaseControlAnnotation caseControlStatus = CaseControlAnnotation.getCaseAnnotationForTriTyper(lineElems[1]);
					if (caseControlStatus == CaseControlAnnotation.CASE) {
						numCase++;
					} else if (caseControlStatus == CaseControlAnnotation.CONTROL) {
						numControl++;
					} else {
						numUnknown++;
					}
					SexAnnotation sex = SexAnnotation.getSexAnnotationForTriTyper(lineElems[3]);
					if (sex == SexAnnotation.FEMALE) {
						numFemale++;
					} else if (sex == SexAnnotation.MALE) {
						numMale++;
					} else {
						numUnknownSex++;
					}
					sampleObj.putAnnotationValues(GenotypeData.BOOL_INCLUDE_SAMPLE, includeSample);
					sampleObj.putAnnotationValues(GenotypeData.CASE_CONTROL_SAMPLE_ANNOTATION_NAME, caseControlStatus);
					sampleObj.putAnnotationValues(GenotypeData.SEX_SAMPLE_ANNOTATION_NAME, sex);
					visitedSamples.add(sampleObj);
				}
			}
			lineElems = t.readLineElemsReturnReference(TextFile.tab);
		}
		t.close();

		if (sampleFilter != null) {
			includedSamples = new ArrayList<Sample>(numIncluded);
			for (Sample sample : samples) {
				if (sampleFilter.doesSamplePassFilter(sample)) {
					includedSamples.add(sample);
				}
			}
		} else {
			includedSamples = samples;
		}

		LOG.info("Loaded " + includedSamples.size() + " out of " + samples.size() + " samples.");

		sampleAnnotationMap = new HashMap<String, SampleAnnotation>(3);
		sampleAnnotationMap.put(GenotypeData.BOOL_INCLUDE_SAMPLE, new SampleAnnotation(BOOL_INCLUDE_SAMPLE, BOOL_INCLUDE_SAMPLE, null, Annotation.Type.BOOLEAN, SampleAnnotation.SampleAnnotationType.OTHER, false));
		sampleAnnotationMap.put(GenotypeData.CASE_CONTROL_SAMPLE_ANNOTATION_NAME, new SampleAnnotation(CASE_CONTROL_SAMPLE_ANNOTATION_NAME, CASE_CONTROL_SAMPLE_ANNOTATION_NAME, null, Annotation.Type.CASECONTROL, SampleAnnotation.SampleAnnotationType.PHENOTYPE, false));
		sampleAnnotationMap.put(GenotypeData.SEX_SAMPLE_ANNOTATION_NAME, new SampleAnnotation(SEX_SAMPLE_ANNOTATION_NAME, SEX_SAMPLE_ANNOTATION_NAME, null, Annotation.Type.SEX, SampleAnnotation.SampleAnnotationType.COVARIATE, false));
	}

	private Boolean parseIncludeExcludeStatus(String status) {
		if (status == null) {
			return false;
		} else if (status.toLowerCase().equals("exclude")) {
			return false;
		} else if (status.toLowerCase().equals("include")) {
			return true;
		} else {
			return false;
		}
	}

	private void loadSNPAnnotation(GeneticVariantRange.ClassGeneticVariantRangeCreate snpsFactory) throws IOException {
		TextFile tf = new TextFile(snpFile, TextFile.R);

		final TObjectIntHashMap<String> allSNPHash = new TObjectIntHashMap<String>();

		unfilteredSnpCount = 0;
		for (String line : tf) {
			if (variantFilter == null || variantFilter.doesIdPassFilter(line)) {
				allSNPHash.put(line, unfilteredSnpCount);
			}
			++unfilteredSnpCount;
		}
		tf.close();


		TextFile tfSNPMap = new TextFile(snpMapFile, TextFile.R);

		int numberOfIncludedSNPsWithAnnotation = 0;

		sequences = new HashMap<String, Sequence>();

		int lineCount = 0;
		for (String[] chrPosId : tfSNPMap.readLineElemsIterable(TextFile.tab)) {
			++lineCount;
			if (chrPosId.length != 3) {
				throw new GenotypeDataException("Error in Trityper SNPMappings.txt. Line number " + lineCount + " does not contain 3 elements: ");
			}

			if (allSNPHash.containsKey(chrPosId[2])) {
				String snp = chrPosId[2];

				int pos = 0;

				String chr = chrPosId[0].intern();
				if (!sequences.containsKey(chr)) {
					if (!chr.equals("0")) {
						sequences.put(chr, new SimpleSequence(chr, 0, this));
					}
				}

				try {
					pos = Integer.parseInt(chrPosId[1]);
				} catch (NumberFormatException e) {
					throw new GenotypeDataException("Position defined for " + snp + " on chromosome " + chr + " is not an integer: " + chrPosId[1]);
				}

				//Index will be removed from snp hash so we can later iterator over remaing without a mapping
				GeneticVariant variant = new ReadOnlyGeneticVariantTriTyper(snp, pos, chr, variantProvider, allSNPHash.remove(snp));

				if (variantFilter == null || variantFilter.doesVariantPassFilter(variant)) {
					snpsFactory.addVariant(variant);
					numberOfIncludedSNPsWithAnnotation++;
				}

			}

		}

		tfSNPMap.close();

		//loop over reaming variant without annotation
		for (String variantId : allSNPHash.keySet()) {
			GeneticVariant variant = new ReadOnlyGeneticVariantTriTyper(variantId, 0, "0", variantProvider, allSNPHash.get(variantId));

			if (variantFilter == null || variantFilter.doesVariantPassFilter(variant)) {
				snpsFactory.addVariant(variant);
			}

		}

		LOG.info("Loaded " + snpsFactory.size() + " out of " + unfilteredSnpCount + " SNPs, " + numberOfIncludedSNPsWithAnnotation + " of loaded SNPs have annotation.");
	}

	private void checkFileSize() {

//		System.out.println(unfilteredSnpCount);

		long expectedfilesize = (long) (unfilteredSnpCount * 2) * (long) samples.size();
		long detectedsize = genotypeDataFile.length();
		if (expectedfilesize != detectedsize) {
			throw new GenotypeDataException("Size of GenotypeMatrix.dat does not match size defined by Indivuals.txt and SNPs.txt. Expected size: " + expectedfilesize + " (" + Gpio.humanizeFileSize(expectedfilesize) + ")\tDetected size: " + detectedsize + " (" + Gpio.humanizeFileSize(detectedsize) + ")\tDiff: " + Math.abs(expectedfilesize - detectedsize));
		}

		if (imputedDosageDataFile != null) {
			expectedfilesize = (long) (unfilteredSnpCount) * (long) samples.size();
			detectedsize = imputedDosageDataFile.length();

			if (expectedfilesize != detectedsize) {
				throw new GenotypeDataException("Size of ImputedDosageMatrix.dat does not match size defined by Indivuals.txt and SNPs.txt. Expected size: " + expectedfilesize + " (" + Gpio.humanizeFileSize(expectedfilesize) + ")\tDetected size: " + detectedsize + " (" + Gpio.humanizeFileSize(detectedsize) + ")\tDiff: " + Math.abs(expectedfilesize - detectedsize));
			}
		}
	}

	@Override
	public boolean isOnlyContaingSaveProbabilityGenotypes() {
		return imputedDosageDataFile == null;
	}

	@Override
	public float[][] getSampleProbilities(GeneticVariant variant) {
		return ProbabilitiesConvertor.convertDosageToProbabilityHeuristic(variant.getSampleDosages());
	}

	@Override
	public List<Alleles> getSampleVariants(GeneticVariant variant) {

		//This is save to do because it would not make sence that a non trityper variant would call this functioon. Unless someone is hacking the api (which they should not do) :)
		int index = ((ReadOnlyGeneticVariantTriTyper) variant).getIndexOfVariantInTriTyperData();

		int numIndividuals = samples.size();
		long indexLong = (long) (index) * (numIndividuals * 2);

		byte[] buffer = new byte[2 * numIndividuals];
		try {
			genotypeHandle.seek(indexLong);
			if (genotypeHandle.read(buffer) != buffer.length) {
				throw new GenotypeDataException("Could not read bytes from: " + indexLong + " in genotype file " + genotypeDataFile.getAbsolutePath() + " (size: " + genotypeDataFile.length() + ")");
			}

		} catch (IOException e) {
			throw new GenotypeDataException("Could not read bytes from: " + indexLong + " in genotype file " + genotypeDataFile.getAbsolutePath() + " (size: " + genotypeDataFile.length() + ")");
		}

		List<Alleles> alleles = new ArrayList<Alleles>(includedSamples.size());

		for (int i = 0; i < numIndividuals; i++) {
			if (sampleFilter == null || sampleFilter.doesSamplePassFilter(samples.get(i))) {
				int allele2Pos = numIndividuals + i;
				Alleles a = Alleles.createAlleles(TriTyperAlleleAnnotation.convertByteToAllele(buffer[i]), TriTyperAlleleAnnotation.convertByteToAllele(buffer[allele2Pos]));
				alleles.add(a);
			}
		}

		return alleles;
	}

	@Override
	public float[] getSampleDosage(GeneticVariant variant) {
		// if there is a dosage file, read from there.. if not, conver genotypes.
		// now transcode into dosage..

		// TODO: optimize this step: no need to get ALL alleles.
		float[] genotypes = CalledDosageConvertor.convertCalledAllelesToDosage(variant.getSampleVariants(), variant.getVariantAlleles(), variant.getRefAllele());
		if (imputedDosageDataFile != null) {

			//This is save to do because it would not make sence that a non trityper variant would call this functioon. Unless someone is hacking the api (which they should not do) :)
			int index = ((ReadOnlyGeneticVariantTriTyper) variant).getIndexOfVariantInTriTyperData();

			int numIndividuals = samples.size();
			long indexLong = (long) index * (long) numIndividuals;
			ByteBuffer buffer = ByteBuffer.allocate(numIndividuals);
			try {
				dosageChannel.read(buffer, indexLong);
			} catch (IOException e) {
				throw new GenotypeDataException("Could not read bytes from: " + indexLong
						+ " in genotype file " + genotypeDataFile.getAbsolutePath()
						+ " (size: " + genotypeDataFile.length() + ")");
			}

			byte[] dosageValuesAll = buffer.array();
			byte[] dosageValues;

			//Filter on included samples
			if (sampleFilter == null) {
				dosageValues = dosageValuesAll;
			} else {
				dosageValues = new byte[includedSamples.size()];
				for (int i = 0, j = 0; i < dosageValuesAll.length; i++) {

					if (sampleFilter.doesSamplePassFilter(samples.get(i))) {
						dosageValues[j] = dosageValuesAll[i];
						++j;
					}

				}
			}


			boolean takeComplement = false;
			for (int ind = 0; ind < dosageValues.length; ind++) {

				if (dosageValues[ind] != 127) {
					double dosagevalue = ((double) (-Byte.MIN_VALUE + dosageValues[ind])) / 100;
					if (genotypes[ind] == 0 && dosagevalue > 1) {
						takeComplement = true;
						break;
					}
					if (genotypes[ind] == 2 && dosagevalue < 1) {
						takeComplement = true;
						break;
					}
				}
			}
			if (takeComplement) {
				for (int ind = 0; ind < dosageValues.length; ind++) {
					if (dosageValues[ind] != 127) {
						byte dosageValue = (byte) (200 - (-Byte.MIN_VALUE + dosageValues[ind]) + Byte.MIN_VALUE);
						dosageValues[ind] = dosageValue;
					}
				}
			}
			float[] dosageValuesFloat = new float[includedSamples.size()];
			for (int i = 0; i < dosageValues.length; i++) {

				if (dosageValues[i] == 127) {
					dosageValuesFloat[i] = -1;
				} else {
					dosageValuesFloat[i] = ((float) (-Byte.MIN_VALUE + dosageValues[ i])) / 100;
				}


			}
			return dosageValuesFloat;
		} else {
			return genotypes;
		}
	}

	@Override
	public byte[] getSampleCalledDosage(GeneticVariant variant) {
		byte[] genotypes = CalledDosageConvertor.convertCalledAllelesToCalledDosage(variantProvider.getSampleVariants(variant),
				variant.getVariantAlleles(), variant.getRefAllele());

		return genotypes;
	}

	@Override
	public List<Sample> getSamples() {
		return includedSamples;
	}

	@Override
	public List<Boolean> getSamplePhasing(GeneticVariant variant) {
		return samplePhasing;
	}

	@Override
	public int cacheSize() {
		return cacheSize;
	}

	@Override
	public void close() throws IOException {
		try {
			if (imputedDosageDataFile != null) {
				dosageChannel.close();
				dosageHandle.close();
			}
			genotypeHandle.close();
		} catch (IOException e) {
			throw new GenotypeDataException("Could not close file handle to TriTyper file: " + genotypeDataFile);
		}
	}

	@Override
	public int getSampleVariantProviderUniqueId() {
		return sampleVariantProviderUniqueId;
	}

	@Override
	public Map<String, SampleAnnotation> getSampleAnnotationsMap() {
		return sampleAnnotationMap;
	}

	@Override
	public Map<String, Annotation> getVariantAnnotationsMap() {
		return Collections.emptyMap();
	}

	@Override
	public List<String> getSeqNames() {
		return new ArrayList<String>(sequences.keySet());
	}

	@Override
	public Iterable<Sequence> getSequences() {
		return sequences.values();
	}

	@Override
	public Iterable<GeneticVariant> getVariantsByPos(String seqName, int startPos) {
		return snps.getVariantAtPos(seqName, startPos);
	}

	@Override
	public Iterable<GeneticVariant> getSequenceGeneticVariants(String seqName) {
		return snps.getVariantsBySequence(seqName);
	}

	@Override
	public Iterable<GeneticVariant> getVariantsByRange(String seqName, int rangeStart, int rangeEnd) {
		return snps.getVariantsByRange(seqName, rangeStart, rangeEnd);
	}

	@Override
	public Iterator<GeneticVariant> iterator() {
		return snps.iterator();
	}
}
