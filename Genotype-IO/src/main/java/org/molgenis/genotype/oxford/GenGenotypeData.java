/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.oxford;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import org.apache.log4j.Logger;
import org.molgenis.genotype.AbstractRandomAccessGenotypeData;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.GenotypeDataException;
import org.molgenis.genotype.Sample;
import org.molgenis.genotype.Sequence;
import org.molgenis.genotype.SimpleSequence;
import org.molgenis.genotype.annotation.Annotation;
import org.molgenis.genotype.annotation.SampleAnnotation;
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
import org.molgenis.genotype.variant.range.GeneticVariantRange.GeneticVariantRangeCreate;
import org.molgenis.genotype.variant.sampleProvider.CachedSampleVariantProvider;
import org.molgenis.genotype.variant.sampleProvider.SampleVariantUniqueIdProvider;
import org.molgenis.genotype.variant.sampleProvider.SampleVariantsProvider;

/**
 *
 * @author Patrick Deelen
 */
public class GenGenotypeData extends AbstractRandomAccessGenotypeData implements SampleVariantsProvider {

	private final RandomAccessFile genFileReader;
	private Map<String, SampleAnnotation> sampleAnnotations;
	private final int sampleVariantProviderUniqueId;
	private final SampleVariantsProvider sampleVariantProvider;
	private final GeneticVariantRange variants;
	private final LinkedHashMap<GeneticVariant, Long> variantSampleAllelesIndex;
	private final List<Sample> samples;
	private final LinkedHashSet<String> sequenceNames;
	private final int byteToReadForSampleAlleles;
	private static final Logger LOGGER = Logger.getLogger(GenGenotypeData.class);
	private final double minimumPosteriorProbabilityToCall;
	private final List<Boolean> phasing;
	private static final double DEFAULT_MINIMUM_POSTERIOR_PROBABILITY_TO_CALL = 0.4f;
	private GeneticVariantMeta geneticVariantMeta = GeneticVariantMetaMap.getGeneticVariantMetaGp();

	public GenGenotypeData(String path) throws IOException {
		this(new File(path + ".gen"), new File(path + ".sample"));
	}

	public GenGenotypeData(String path, double minimumPosteriorProbabilityToCall) throws IOException {
		this(new File(path + ".gen"), new File(path + ".sample"), minimumPosteriorProbabilityToCall);
	}

	public GenGenotypeData(File genFile, File sampleFile) throws IOException {
		this(genFile, sampleFile, 1000);
	}

	public GenGenotypeData(File genFile, File sampleFile, double minimumPosteriorProbabilityToCall) throws IOException {
		this(genFile, sampleFile, 1000, minimumPosteriorProbabilityToCall);
	}

	public GenGenotypeData(File genFile, File sampleFile, int cacheSize) throws IOException {
		this(genFile, sampleFile, cacheSize, null, DEFAULT_MINIMUM_POSTERIOR_PROBABILITY_TO_CALL);
	}

	public GenGenotypeData(File genFile, File sampleFile, int cacheSize, double minimumPosteriorProbabilityToCall) throws IOException {
		this(genFile, sampleFile, cacheSize, null, minimumPosteriorProbabilityToCall);
	}

	public GenGenotypeData(File genFile, File sampleFile, String forceSeqName) throws IOException {
		this(genFile, sampleFile, 1000, forceSeqName, DEFAULT_MINIMUM_POSTERIOR_PROBABILITY_TO_CALL);
	}

	public GenGenotypeData(File genFile, File sampleFile, int cacheSize, String forceSeqName, double minimumPosteriorProbabilityToCall)
			throws IOException {

		if (genFile == null) {
			throw new IllegalArgumentException("genFile is null");
		}
		if (!genFile.isFile()) {
			throw new FileNotFoundException("gen file file not found at "
					+ genFile.getAbsolutePath());
		}
		if (!genFile.canRead()) {
			throw new IOException("cannot read gen file at "
					+ genFile.getAbsolutePath());
		}

		this.minimumPosteriorProbabilityToCall = minimumPosteriorProbabilityToCall;

		LOGGER.debug("Using " + minimumPosteriorProbabilityToCall + " as cutoff to call genotypes from: " + genFile.getAbsolutePath());

		sampleVariantProviderUniqueId = SampleVariantUniqueIdProvider.getNextUniqueId();
		if (cacheSize > 0) {
			sampleVariantProvider = new CachedSampleVariantProvider(this, cacheSize);
		} else {
			sampleVariantProvider = this;
		}

		OxfordSampleFile oxfordSampleFile = new OxfordSampleFile(sampleFile);

		sampleAnnotations = oxfordSampleFile.getSampleAnnotations();
		samples = oxfordSampleFile.getSamples();

		LOGGER.info("Loaded " + samples.size() + " samples from " + sampleFile.getAbsolutePath());

		phasing = Collections.unmodifiableList(Collections.nCopies((int) samples.size(), false));

		GeneticVariantRangeCreate variantRangeFactory = GeneticVariantRange.createRangeFactory();
		variantSampleAllelesIndex = new LinkedHashMap<GeneticVariant, Long>();
		sequenceNames = new LinkedHashSet<String>();
		genFileReader = new RandomAccessFile(genFile, "r");

		byteToReadForSampleAlleles = loadVariants(forceSeqName, variantRangeFactory);

		variants = variantRangeFactory.createRange();
		
		LOGGER.info("Loaded " + variants.size() + " variants from " + genFile.getAbsolutePath());

	}

	@Override
	public List<Sequence> getSequences() {
		List<Sequence> sequences = new ArrayList<Sequence>();
		for (String seqName : getSeqNames()) {
			sequences.add(new SimpleSequence(seqName, null, this));
		}

		return sequences;
	}

	@Override
	public List<Sample> getSamples() {
		return Collections.unmodifiableList(samples);
	}

	@Override
	public Map<String, Annotation> getVariantAnnotationsMap() {
		return Collections.emptyMap();
	}

	@Override
	public Map<String, SampleAnnotation> getSampleAnnotationsMap() {
		return sampleAnnotations;
	}

	@Override
	public List<Alleles> getSampleVariants(GeneticVariant variant) {

		return ProbabilitiesConvertor.convertProbabilitiesToAlleles(variant.getSampleGenotypeProbilities(), variant.getVariantAlleles(), minimumPosteriorProbabilityToCall);

	}

	@Override
	public List<Boolean> getSamplePhasing(GeneticVariant variant) {
		return phasing;
	}

	@Override
	public int cacheSize() {
		return 0;
	}

	@Override
	public int getSampleVariantProviderUniqueId() {
		return sampleVariantProviderUniqueId;
	}

	@Override
	public byte[] getSampleCalledDosage(GeneticVariant variant) {
		return CalledDosageConvertor.convertCalledAllelesToCalledDosage(variant.getSampleVariants(), variant.getVariantAlleles(), null);
	}

	@Override
	public float[] getSampleDosage(GeneticVariant variant) {
		return ProbabilitiesConvertor.convertProbabilitiesToDosage(variant.getSampleGenotypeProbilities(), minimumPosteriorProbabilityToCall);
	}

	@Override
	public void close() throws IOException {
		genFileReader.close();
	}

	@Override
	public Iterator<GeneticVariant> iterator() {
		return variants.iterator();
	}

	@Override
	public List<String> getSeqNames() {
		return new ArrayList<String>(sequenceNames);
	}

	@Override
	public Iterable<GeneticVariant> getVariantsByPos(String seqName, int startPos) {
		return variants.getVariantAtPos(seqName, startPos);
	}

	@Override
	public Iterable<GeneticVariant> getSequenceGeneticVariants(String seqName) {
		return variants.getVariantsBySequence(seqName);
	}

	@Override
	public Iterable<GeneticVariant> getVariantsByRange(String seqName, int rangeStart, int rangeEnd) {
		return variants.getVariantsByRange(seqName, rangeStart, rangeEnd);
	}

	/**
	 * Very complicated low level reading of the variants. Also stores for each
	 * variant the start pos of the sample alleles in the file
	 *
	 * @return the number of bytes of the longest chunk of sample alleles
	 * @throws IOException
	 */
	private int loadVariants(String forceSeqName, GeneticVariantRangeCreate variantRangeFactory) throws IOException {
		StringBuilder stringBuilder = new StringBuilder();
		byte[] buffer = new byte[8192];
		boolean eol = false;

		String seqName = null;
		String variantId = null;
		int position = 0;
		String allele1 = null;
		String allele2;

		int longestedChunk = 0;
		int currentChunk = 0;

		int column = 0;
		long posBeforeBufferRead;

		while (!eol) {

			posBeforeBufferRead = genFileReader.getFilePointer();
			int bytesRead = genFileReader.read(buffer);

			if (bytesRead == -1) {
				eol = true;
			} else {
				for (int i = 0; i < bytesRead; ++i) {
					switch (buffer[i]) {
						case '\n':
						case '\r':
							if (currentChunk == 0) {
								//Ignore empty lines or second line break
								currentChunk = -1;
								continue;
							}
							if (column != 5) {
								LOGGER.fatal("Error reading gen file, did not detect first 5 columns with variant information \n"
										+ "current column is:" + column + "\n"
										+ "content in current column: " + stringBuilder.toString());
								throw new GenotypeDataException("Error reading gen file, did not detect first 5 columns with variant information. Note gen files must be space separted");
							}
							longestedChunk = longestedChunk < currentChunk ? currentChunk : longestedChunk;
							column = 0;
							currentChunk = -1;
							break;
						case ' ':
							switch (column) {
								case 0:
									seqName = forceSeqName == null ? stringBuilder.toString().intern() : forceSeqName;
									sequenceNames.add(seqName);
									++column;
									stringBuilder = new StringBuilder();
									break;
								case 1:
									variantId = stringBuilder.toString();
									++column;
									stringBuilder = new StringBuilder();
									break;
								case 2:
									try {
										position = Integer.parseInt(stringBuilder.toString());
									} catch (NumberFormatException ex) {
										throw new GenotypeDataException("Error parsing \"" + stringBuilder.toString() + "\" as position in haps file");
									}
									++column;
									stringBuilder = new StringBuilder();
									break;
								case 3:
									allele1 = stringBuilder.toString();
									++column;
									stringBuilder = new StringBuilder();
									break;
								case 4:
									allele2 = stringBuilder.toString();
									GeneticVariant variant = ReadOnlyGeneticVariant.createVariant(geneticVariantMeta, variantId, position, seqName, sampleVariantProvider, allele1, allele2);
									variantRangeFactory.addVariant(variant);
									variantSampleAllelesIndex.put(variant, posBeforeBufferRead + i + 1);
									++column;
									stringBuilder = new StringBuilder();
									currentChunk = -1; // because otherwise we would count this space
									break;
								default:
									//We do not care about the rest of the information. We simply seak the end off line
									break;
							}
							break;
						default:
							if (column < 5) {
								//We do not care about the rest of the information. We simply seak the end off line
								stringBuilder.append((char) buffer[i]);
							}
							break;
					}
					++currentChunk;
				}
			}
		}
		return longestedChunk;

	}

	@Override
	public boolean isOnlyContaingSaveProbabilityGenotypes() {
		return true;
	}

	@Override
	public float[][] getSampleProbilities(GeneticVariant variant) {


		float[][] probs = new float[getSamples().size()][3];

		long start = variantSampleAllelesIndex.get(variant);

		byte[] buffer = new byte[byteToReadForSampleAlleles];

		int bytesRead;

		try {
			synchronized (genFileReader) {
				genFileReader.seek(start);
				bytesRead = genFileReader.read(buffer);
			}
		} catch (IOException ex) {
			throw new GenotypeDataException("Error loading probs for variant: " + variant.getPrimaryVariantId() + " from gen file");
		}

		if (bytesRead == -1) {
			throw new GenotypeDataException("Error loading probs for variant: " + variant.getPrimaryVariantId() + " from gen file");
		}



		int i = 0;
		for (int s = 0; s < samples.size(); ++s) {

			float[] sampleProbs = new float[3];
			int currentProbIndex = 0;
			StringBuilder currentProb = new StringBuilder();

			while (currentProbIndex < 3) {

				if (i == bytesRead) {

					if (currentProb.length() == 0) {
						throw new GenotypeDataException("Error parsing gen file: " + variant.getPrimaryVariantId() + " probability " + (currentProbIndex + 1) + " of sample " + samples.get(s).getId() + " sample index: " + s + " empty probability");
					}

					try {
						sampleProbs[currentProbIndex] = Float.parseFloat(currentProb.toString());
					} catch (NumberFormatException e) {
						throw new GenotypeDataException("Error parsing gen file: " + variant.getPrimaryVariantId() + " probability " + (currentProbIndex + 1) + " of sample " + samples.get(s).getId() + " sample index: " + s + " error: " + e.getMessage());
					}
					++currentProbIndex;
					break;
				}
				
				switch ((char) buffer[i]) {
					case '\n':
					case '\r':
						//if not at last probability of line give error
						if(s < samples.size() - 1 && currentProbIndex == 2){
							throw new GenotypeDataException("Error parsing gen file: variant ID: " + variant.getPrimaryVariantId() + " unexpected new line when parsing sample " + samples.get(s).getId() + " sample index: " + s);
						}
						//if not at premature line ending just parse the last probability
					case ' ':
						try {
							sampleProbs[currentProbIndex] = Float.parseFloat(currentProb.toString());
						} catch (NumberFormatException e) {
							throw new GenotypeDataException("Error parsing gen file: variant ID: " + variant.getPrimaryVariantId() + " genotype probability value " + (currentProbIndex + 1) + " of sample " + samples.get(s).getId() + " sample index: " + s + " problem parsing probability with value \"" + currentProb.toString() + "\": " + e.getMessage());
						}
						currentProb = new StringBuilder();
						++currentProbIndex;
						break;
					
					default:
						currentProb.append((char) buffer[i]);
				}

				++i;

			}

			probs[s] = sampleProbs;

		}


		return probs;



	}

	@Override
	public double[][] getSampleGenotypeProbabilitiesBgen(GeneticVariant variant) {
		return ProbabilitiesConvertor.convertProbabilitiesToBgenProbabilities(getSampleProbilities(variant));
	}

	@Override
	public double[][][] getSampleGenotypeProbabilitiesBgenPhased(GeneticVariant variant) {
		throw new UnsupportedOperationException("Not supported yet.");
	}

	@Override
	public FixedSizeIterable<GenotypeRecord> getSampleGenotypeRecords(GeneticVariant variant) {
		
		return RecordIteratorCreators.createIteratorFromProbs(variant.getSampleGenotypeProbilities());
		
	}
}
