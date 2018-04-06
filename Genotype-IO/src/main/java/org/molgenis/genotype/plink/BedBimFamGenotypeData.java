/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.plink;

import gnu.trove.map.hash.TObjectIntHashMap;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.RandomAccessFile;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Pattern;
import org.molgenis.genotype.AbstractRandomAccessGenotypeData;
import org.molgenis.genotype.Allele;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.GenotypeDataException;
import org.molgenis.genotype.Sample;
import org.molgenis.genotype.Sequence;
import org.molgenis.genotype.SimpleSequence;
import org.molgenis.genotype.annotation.Annotation;
import org.molgenis.genotype.annotation.SampleAnnotation;
import org.molgenis.genotype.annotation.SexAnnotation;
import org.molgenis.genotype.plink.readers.FamFileReader;
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
import org.molgenis.genotype.variant.sampleProvider.CachedSampleVariantProvider;
import org.molgenis.genotype.variant.sampleProvider.SampleVariantUniqueIdProvider;
import org.molgenis.genotype.variant.sampleProvider.SampleVariantsProvider;

/**
 *
 * @author Patrick Deelen
 */
public class BedBimFamGenotypeData extends AbstractRandomAccessGenotypeData implements SampleVariantsProvider {

	private static final Pattern SEPARATOR_PATTERN = Pattern.compile("[ \\t]+");
	private static final byte MAGIC_NUMBER_1 = 108;
	private static final byte MAGIC_NUMBER_2 = 27;
	private static final byte MODE = 1; //We only write snp major mode
	private static final int READER_MASK = 3;
	private static final int HOMOZYGOTE_FIRST = 0;
	private static final int HOMOZYGOTE_SECOND = 3;
	private static final int HETEROZYGOTE = 2;
	private static final int MISSING = 1;
	private static final Alleles BI_ALLELIC_MISSING = Alleles.createAlleles(Allele.ZERO, Allele.ZERO);
	private static final org.apache.log4j.Logger LOGGER = org.apache.log4j.Logger.getLogger(BedBimFamGenotypeWriter.class);
	private static final Charset FILE_ENCODING = Charset.forName("UTF-8");
	private final ArrayList<Sample> samples;
	private final Map<String, SampleAnnotation> sampleAnnotations;
	private final LinkedHashMap<String, Sequence> sequences;
	private final GeneticVariantRange snps;
	private final TObjectIntHashMap<GeneticVariant> snpIndexces;
	private final RandomAccessFile bedFileReader;
	private final SampleVariantsProvider sampleVariantProvider;
	private final int sampleVariantProviderUniqueId;
	private final int cacheSize;
	private final List<Boolean> phasing;
	private final long bytesPerVariant;
	/**
	 * The original SNP count in the data regardless of number of read SNPs
	 */
	private final int originalSnpCount;
	private GeneticVariantMeta geneticVariantMeta = GeneticVariantMetaMap.getGeneticVariantMetaGt();

	public BedBimFamGenotypeData(String basePath) throws IOException {
		this(basePath, 100);
	}

	public BedBimFamGenotypeData(String basePath, int cacheSize) throws IOException {
		this(new File(basePath + ".bed"), new File(basePath + ".bim"), new File(basePath + ".fam"), cacheSize);
	}

	public BedBimFamGenotypeData(File bedFile, File bimFile, File famFile, int cacheSize) throws IOException {

		if (bedFile == null) {
			throw new IllegalArgumentException("BedFile is null");
		}
		if (bimFile == null) {
			throw new IllegalArgumentException("BimFile is null");
		}
		if (famFile == null) {
			throw new IllegalArgumentException("FamFile is null");
		}

		if (!bedFile.isFile()) {
			throw new FileNotFoundException("BED file not found at "
					+ bedFile.getAbsolutePath());
		}
		if (!bedFile.canRead()) {
			throw new IOException("BED file not found at " + bedFile.getAbsolutePath());
		}

		if (!bimFile.isFile()) {
			throw new FileNotFoundException("BIM file not found at " + bimFile.getAbsolutePath());
		}
		if (!bimFile.canRead()) {
			throw new IOException("BIM file not found at " + bimFile.getAbsolutePath());
		}

		if (!famFile.isFile()) {
			throw new FileNotFoundException("FAM file not found at " + famFile.getAbsolutePath());
		}
		if (!famFile.canRead()) {
			throw new IOException("FAM file not found at " + famFile.getAbsolutePath());
		}

		sampleVariantProviderUniqueId = SampleVariantUniqueIdProvider.getNextUniqueId();

		if (cacheSize <= 0) {
			sampleVariantProvider = this;
		} else {
			sampleVariantProvider = new CachedSampleVariantProvider(this, cacheSize);
		}
		this.cacheSize = cacheSize;

		sampleAnnotations = PlinkSampleAnnotations.getSampleAnnotations();
		samples = FamFileReader.readFamFile(famFile);

		phasing = Collections.unmodifiableList(Collections.nCopies((int) samples.size(), false));

		snpIndexces = new TObjectIntHashMap<GeneticVariant>(10000, 0.75f, -1);
		GeneticVariantRange.GeneticVariantRangeCreate snpsFactory = GeneticVariantRange.createRangeFactory();
		sequences = new LinkedHashMap<String, Sequence>();
		originalSnpCount = readBimFile(bimFile, snpsFactory);
		snps = snpsFactory.createRange();

		bytesPerVariant = samples.size() % 4 == 0 ? samples.size() / 4 : (samples.size() / 4 + 1);

		//Check file size of bed file
		if (bedFile.length() != ((long) bytesPerVariant * (long) originalSnpCount + 3)) {
			throw new GenotypeDataException("Invalid plink BED file not the expected file size. " + bedFile.getAbsolutePath() + " expected: " + ((long) bytesPerVariant * (long) snpIndexces.size() + 3) + " found: " + bedFile.length() + " bytes per variant: " + bytesPerVariant + " snps: " + originalSnpCount);
		}

		//Check first two bytes for magic number
		bedFileReader = new RandomAccessFile(bedFile, "r");
		if (bedFileReader.read() != MAGIC_NUMBER_1 || bedFileReader.read() != MAGIC_NUMBER_2) {
			throw new GenotypeDataException("Error reading plink BED file, magic number not found. " + bedFile.getAbsolutePath());
		}

		int bedFileMode = bedFileReader.read();
		if (bedFileMode != MODE) {
			if (bedFileMode == 0) {
				throw new GenotypeDataException("Error reading BED file, only SNP major mode is supported. " + bedFile.getAbsolutePath());
			} else {
				throw new GenotypeDataException("Error reading BED file, ivalid mode byte detected. " + bedFile.getAbsolutePath());
			}
		}

	}

	@Override
	public Map<String, SampleAnnotation> getSampleAnnotationsMap() {
		return sampleAnnotations;
	}

	@Override
	public Map<String, Annotation> getVariantAnnotationsMap() {
		return Collections.emptyMap();
	}

	@Override
	public List<Sample> getSamples() {
		return Collections.unmodifiableList(samples);
	}

	@Override
	public void close() throws IOException {
		bedFileReader.close();
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
	public List<Alleles> getSampleVariants(GeneticVariant variant) {

		int index = snpIndexces.get(variant);
		
		if(index == -1){
			throw new GenotypeDataException("Error reading variant from bed file. ID: " + variant.getPrimaryVariantId() + " chr: " + variant.getSequenceName() + " pos: " + variant.getStartPos() + " alleles" + variant.getVariantAlleles().toString());
		}

		long startByte = (index * bytesPerVariant) + 3;

		long stopByte = startByte + bytesPerVariant;

		byte[] variantBytes = new byte[(int) (stopByte - startByte)];
		try {
			bedFileReader.seek(startByte);
			if (bedFileReader.read(variantBytes) != variantBytes.length) {
				throw new GenotypeDataException("Error reading bed file");
			}
		} catch (IOException ex) {
			throw new GenotypeDataException("Error reading bed file", ex);
		}

		ArrayList<Alleles> alleles = new ArrayList<Alleles>(samples.size());

		Alleles heterozygote = variant.getVariantAlleles();
		Alleles homozygoteFirst = Alleles.createAlleles(heterozygote.get(0), heterozygote.get(0));
		Alleles homozygoteSecond = Alleles.createAlleles(heterozygote.get(1), heterozygote.get(1));

		int sampleCounter = 0;

		for (int variantByte : variantBytes) {

			for (int i = 0; i < 4; ++i) {

				if (sampleCounter < samples.size()) {
					switch (variantByte & READER_MASK) {
						case HOMOZYGOTE_FIRST:
							alleles.add(homozygoteFirst);
							break;
						case HOMOZYGOTE_SECOND:
							alleles.add(homozygoteSecond);
							break;
						case HETEROZYGOTE:
							alleles.add(heterozygote);
							break;
						case MISSING:
							alleles.add(BI_ALLELIC_MISSING);
							break;
						default:
							throw new GenotypeDataException("Error reading BED, this should not be reachable");
					}
				} else {
					if ((variantByte & READER_MASK) != 0) {
						throw new GenotypeDataException("Error reading BED file, found data in padding bits of variant: " + variant.getPrimaryVariantId());
					}
				}
				variantByte = variantByte >>> 2;
				++sampleCounter;
			}

		}


		return Collections.unmodifiableList(alleles);

	}

	@Override
	public List<Boolean> getSamplePhasing(GeneticVariant variant) {
		return phasing;
	}

	@Override
	public int cacheSize() {
		return cacheSize;
	}

	@Override
	public int getSampleVariantProviderUniqueId() {
		return sampleVariantProviderUniqueId;
	}

	@Override
	public byte[] getSampleCalledDosage(GeneticVariant variant) {
		return CalledDosageConvertor.convertCalledAllelesToCalledDosage(getSampleVariants(variant),
				variant.getVariantAlleles(), variant.getRefAllele());
	}

	@Override
	public float[] getSampleDosage(GeneticVariant variant) {
		return CalledDosageConvertor.convertCalledAllelesToDosage(getSampleVariants(variant),
				variant.getVariantAlleles(), variant.getRefAllele());
	}

	private void readFamFile(File famFile) throws FileNotFoundException, IOException {

		BufferedReader famFileReader = new BufferedReader(new InputStreamReader(new FileInputStream(famFile), FILE_ENCODING));

		String line;
		while ((line = famFileReader.readLine()) != null) {

			String[] elements = SEPARATOR_PATTERN.split(line);

			Map<String, Object> annotationValues = new LinkedHashMap<String, Object>();
			annotationValues.put(FATHER_SAMPLE_ANNOTATION_NAME, elements[2]);
			annotationValues.put(MOTHER_SAMPLE_ANNOTATION_NAME, elements[3]);
			annotationValues.put(SEX_SAMPLE_ANNOTATION_NAME, SexAnnotation.getSexAnnotationForPlink(Byte.parseByte(elements[4])));
			annotationValues.put(DOUBLE_PHENOTYPE_SAMPLE_ANNOTATION_NAME, Double.parseDouble(elements[5]));

			samples.add(new Sample(elements[1], elements[0], annotationValues));

		}

		LOGGER.info("Read " + samples.size() + " samples from " + famFile.getAbsolutePath());

		famFileReader.close();

	}

	private int readBimFile(File bimFile, GeneticVariantRange.GeneticVariantRangeCreate snpsFactory) throws FileNotFoundException, IOException {

		BufferedReader bimFileReader = new BufferedReader(new InputStreamReader(new FileInputStream(bimFile), FILE_ENCODING));

		String line;
		int snpIndex = 0;
		while ((line = bimFileReader.readLine()) != null) {

			String[] elements = SEPARATOR_PATTERN.split(line);

			String sequenceName = elements[0].intern();

			if (!sequences.containsKey(sequenceName)) {
				sequences.put(sequenceName, new SimpleSequence(sequenceName, 0, this));
			}

			//Create new strign to make sure it is not backed by the whole line
			GeneticVariant variant = ReadOnlyGeneticVariant.createVariant(geneticVariantMeta, new String(elements[1]), Integer.parseInt(elements[3]), sequenceName, sampleVariantProvider, Allele.create(elements[4]), Allele.create(elements[5]));

			if (snpIndexces.containsKey(variant)) {
				LOGGER.warn("Found two SNPs at " + sequenceName + ":" + variant.getStartPos() + " Only first is read!");
			} else {
				snpsFactory.addVariant(variant);
				snpIndexces.put(variant, snpIndex);
			}



			++snpIndex;

		}

		LOGGER.info("Read " + snpIndex + " SNPs from " + bimFile.getAbsolutePath());

		bimFileReader.close();

		return snpIndex;

	}

	@Override
	public Iterator<GeneticVariant> iterator() {
		return snps.iterator();
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
	public FixedSizeIterable<GenotypeRecord> getSampleGenotypeRecords(GeneticVariant variant) {
		
		return RecordIteratorCreators.createIteratorFromAlleles(variant.getSampleVariants());
		
	}
}