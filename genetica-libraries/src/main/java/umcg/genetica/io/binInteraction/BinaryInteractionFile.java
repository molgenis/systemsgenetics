package umcg.genetica.io.binInteraction;

import umcg.genetica.io.binInteraction.gene.BinaryInteractionGene;
import umcg.genetica.io.binInteraction.variant.BinaryInteractionVariantStatic;
import com.google.common.io.CountingInputStream;
import gnu.trove.map.hash.TObjectIntHashMap;
import java.io.BufferedInputStream;
import java.io.Closeable;
import java.io.DataInputStream;
import java.io.EOFException;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;
import java.text.SimpleDateFormat;
import java.util.Arrays;
import java.util.Collections;
import java.util.Date;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import org.molgenis.genotype.Allele;
import umcg.genetica.collections.ChrPosMap;
import umcg.genetica.io.binInteraction.gene.BinaryInteractionGeneStatic;
import umcg.genetica.io.binInteraction.variant.BinaryInteractionVariant;

/**
 *
 * @author Patrick Deelen
 */
public class BinaryInteractionFile implements Closeable {

	//Static variables
	protected static final byte MAGIC_1 = 81;
	protected static final byte MAGIC_2 = 73;
	private static final long POINTER_TO_CLOSED_BOOLEAN = 4;
	private static final byte MAJOR_VERSION = 1;
	private static final byte MINOR_VERSION = 0;
	private static final int NO_ENTRY_INT_MAP = -1;
	private static final SimpleDateFormat DEFAULT_DATE_FORMAT = new java.text.SimpleDateFormat("dd-MM-yyyy HH:mm:ss");
	private static final int BUFFER_SIZE = 8192;
	//Instance variables
	private final File interactionFile;
	private boolean readOnly;
	private final BinaryInteractionCohort[] cohorts;
	private final BinaryInteractionGene[] genes;
	private final BinaryInteractionVariant[] variants;
	private ChrPosMap<BinaryInteractionVariant> variantsByPos = null;
	protected final String[] covariates;
	private final int[][] covariatesTested;
	private final long timeStamp;
	private final boolean allCovariants;
	private final boolean metaAnalysis;
	private final boolean normalQtlStored;
	private final boolean flippedZscoreStored;
	private final String fileDescription;
	private final long interactions;
	private final long startQtlBlock;
	private final long startInteractionBlock;
	private final long sizeQtlBlock;
	protected final long sizeInteractionBlock;
	/**
	 * Number of variant gene combinations upto a variant. Length = variants + 1
	 */
	private final int[] cummulativeGeneCountUpToVariant;
	/**
	 * Number of interactions upto a variant-gene combinations. Length = total
	 * number of variant-gene combinations + 1
	 */
	private final long[] cummalitiveInteractionCountUptoVariantGene;
	private final TObjectIntHashMap<String> variantMap;
	private final TObjectIntHashMap<String> genesMap;
	private final TObjectIntHashMap<String> covariatesMap;
	private RandomAccessFile randomAccess;
	private FileChannel channel;
	private final ByteBuffer qtlBuffer = ByteBuffer.allocate(BUFFER_SIZE);
	private final ByteBuffer interactionBuffer = ByteBuffer.allocate(BUFFER_SIZE);
	private boolean qtlBufferWriting = false;
	private boolean interactionBufferWriting = false;
	private long qtlBufferStart = Long.MIN_VALUE;
	private long interactionBufferStart = Long.MIN_VALUE;
	private long qtlZscoresSet = 0;
	private long interactionZscoresSet = 0;
	private long qtlZscoresRead = 0;
	private long interactionZscoresRead = 0;
	private long interactionWriteBufferFlushed = 0;
	private long qtlWriteBufferFlushed = 0;
	private long interactionReadBufferLoaded = 0;
	private long qtlReadBufferLoaded = 0;

	protected BinaryInteractionFile(File interactionFile, boolean readOnly, BinaryInteractionCohort[] cohorts, BinaryInteractionGene[] genes, BinaryInteractionVariant[] variants, String[] covariates, int[][] covariatesTested, long timeStamp, boolean allCovariants, boolean metaAnalysis, boolean normalQtlStored, boolean flippedZscoreStored, String fileDescription, long interactions, long startQtlBlock, long startInteractionBlock) throws BinaryInteractionFileException, FileNotFoundException, IOException {
		this.interactionFile = interactionFile;
		this.readOnly = readOnly;
		this.cohorts = cohorts;
		this.genes = genes;
		this.variants = variants;
		this.covariates = covariates;
		this.covariatesTested = covariatesTested;
		this.timeStamp = timeStamp;
		this.allCovariants = allCovariants;
		this.metaAnalysis = metaAnalysis;
		this.normalQtlStored = normalQtlStored;
		this.flippedZscoreStored = flippedZscoreStored;
		this.fileDescription = fileDescription;
		this.interactions = interactions;
		this.startQtlBlock = startQtlBlock;
		this.startInteractionBlock = startInteractionBlock;

		this.sizeQtlBlock = calculateSizeNormalQtlBlock(cohorts.length, metaAnalysis);
		this.sizeInteractionBlock = calculateSizeInteractionResultBlock(cohorts.length, flippedZscoreStored, metaAnalysis);

		this.cummulativeGeneCountUpToVariant = new int[variants.length + 1];
		for (int i = 0; i < variants.length; ++i) {
			this.cummulativeGeneCountUpToVariant[i + 1] = this.cummulativeGeneCountUpToVariant[i] + variants[i].getGeneCount();
		}

		this.cummalitiveInteractionCountUptoVariantGene = new long[cummulativeGeneCountUpToVariant[variants.length] + 1];
		{
			int i = 1;
			for (int v = 0; v < variants.length; ++v) {
				int variantGeneCount = variants[v].getGeneCount();
				for (int g = 0; g < variantGeneCount; ++g) {
					cummalitiveInteractionCountUptoVariantGene[i] = cummalitiveInteractionCountUptoVariantGene[i - 1] + (allCovariants ? covariates.length : this.covariatesTested[i - 1].length);
					++i;
				}
			}
		}

		if (cummalitiveInteractionCountUptoVariantGene[cummulativeGeneCountUpToVariant[variants.length]] != interactions) {
			throw new BinaryInteractionFileException("Something went wrong");
		}

		variantMap = new TObjectIntHashMap<String>(variants.length, 0.75f, -1);
		genesMap = new TObjectIntHashMap<String>(genes.length, 0.75f, -1);
		covariatesMap = new TObjectIntHashMap<String>(covariates.length, 0.75f, -1);

		for (int i = 0; i < variants.length; ++i) {
			if (variantMap.put(variants[i].getName(), i) != NO_ENTRY_INT_MAP) {
				throw new BinaryInteractionFileException("Cannot store the same variant twice (" + variants[i].getName() + ")");
			}
		}

		for (int i = 0; i < genes.length; ++i) {
			if (genesMap.put(genes[i].getName(), i) != NO_ENTRY_INT_MAP) {
				throw new BinaryInteractionFileException("Cannot store the same gene twice (" + genes[i].getName() + ")");
			}
		}

		for (int i = 0; i < covariates.length; ++i) {
			if (covariatesMap.put(covariates[i], i) != NO_ENTRY_INT_MAP) {
				throw new BinaryInteractionFileException("Cannot store the same covariate twice (" + covariates[i] + ")");
			}
		}

		if (this.sizeQtlBlock >= BUFFER_SIZE) {
			throw new BinaryInteractionFileException("QTL block size larger than buffer");
		}

		if (this.sizeInteractionBlock >= BUFFER_SIZE) {
			throw new BinaryInteractionFileException("Interaction block size larger than buffer");
		}

		this.open();

	}

	public static BinaryInteractionFile load(File interactionFile) throws FileNotFoundException, IOException, BinaryInteractionFileException {
		return load(interactionFile, true);
	}

	public static BinaryInteractionFile load(File interactionFile, boolean readOnly) throws FileNotFoundException, IOException, BinaryInteractionFileException {

		final BinaryInteractionFileConstructorBuilder builder = new BinaryInteractionFileConstructorBuilder();

		builder.setInteractionFile(interactionFile);
		builder.setReadOnly(readOnly);

		final CountingInputStream inputStreamCounted = new CountingInputStream(new BufferedInputStream(new FileInputStream(interactionFile)));
		final DataInputStream inputStream = new DataInputStream(inputStreamCounted);

		try {

			//First parse heaeder		
			if (inputStream.readByte() != MAGIC_1 || inputStream.readByte() != MAGIC_2) {
				throw new BinaryInteractionFileException("No a valid interaction file");
			}

			if (inputStream.readByte() != MAJOR_VERSION || inputStream.readByte() != MINOR_VERSION) {
				throw new BinaryInteractionFileException("Interaction file version not supported");
			}

			//Test if closed proparly
			if (!inputStream.readBoolean()) {
				throw new BinaryInteractionFileException("Interaction file not properly closed and might be corrupt");
			}

			final boolean allCovariants = inputStream.readBoolean();
			builder.setAllCovariants(allCovariants);

			final boolean metaAnalysis = inputStream.readBoolean();
			builder.setMetaAnalysis(metaAnalysis);

			final boolean normalQtlStored = inputStream.readBoolean();
			builder.setNormalQtlStored(normalQtlStored);

			final boolean flippedZscoreStored = inputStream.readBoolean();
			builder.setFlippedZscoreStored(flippedZscoreStored);

			inputStream.skip(4);//Skip reserved

			builder.setTimeStamp(inputStream.readLong());

			builder.setFileDescription(readString(inputStream));

			final int chortsCount = inputStream.readInt();
			final int chrsCount = inputStream.readInt();
			final int allelesCount = inputStream.readInt();
			final int genesCount = inputStream.readInt();
			final int variantCount = inputStream.readInt();
			final int covariatesCount = inputStream.readInt();

			long totalInteractions = inputStream.readLong();
			builder.setInteractions(totalInteractions);

			builder.setCohorts(readCohorts(inputStream, chortsCount));
			final String[] chrDictionary = readStringArray(inputStream, chrsCount, "Chromosomes");
			final Allele[] alleleDictionary = readAlleles(inputStream, allelesCount);

			final BinaryInteractionVariant[] variants = readVariants(inputStream, variantCount, chrDictionary, alleleDictionary);

			builder.setVariants(variants);
			builder.setGenes(readGenes(inputStream, genesCount, chrDictionary));
			builder.setCovariates(readStringArray(inputStream, covariatesCount, "Chromosomes"));

			final int totalVariantGeneCombinations = getTotalVariantGeneCombinations(variants);

			if (!allCovariants) {
				builder.setCovariatesTested(readCovariatesData(inputStream, totalVariantGeneCombinations, totalInteractions));
			}

			final long startData = inputStreamCounted.getCount();

			final long sizeNormalQtlSection;
			final long startNormalQtlSection;

			if (normalQtlStored) {
				startNormalQtlSection = startData;
				final long sizeQtlBlock = calculateSizeNormalQtlBlock(chortsCount, metaAnalysis);
				sizeNormalQtlSection = sizeQtlBlock * totalVariantGeneCombinations;
			} else {
				sizeNormalQtlSection = 0;
				startNormalQtlSection = -1;
			}

			builder.setStartQtlBlock(startNormalQtlSection);

			final long startInteractionSection = startData + sizeNormalQtlSection;

			builder.setStartInteractionBlock(startInteractionSection);

			final long sizeInteractionSection = calculateSizeInteractionResultBlock(chortsCount, flippedZscoreStored, metaAnalysis) * totalInteractions;

			if (startData + sizeNormalQtlSection + sizeInteractionSection != interactionFile.length()) {
				throw new BinaryInteractionFileException("Incorrect file size. Expected: " + (startData + sizeNormalQtlSection + sizeInteractionSection) + " found: " + interactionFile.length() + " diff: " + (startData + sizeNormalQtlSection + sizeInteractionSection - interactionFile.length()));
			}

			inputStream.close();
			inputStreamCounted.close();

			return builder.createBinaryInteractionFile();

		} catch (EOFException ex) {
			throw new BinaryInteractionFileException("Error parsing header, unexpected EOF", ex);
		}

	}

	private static BinaryInteractionCohort[] readCohorts(DataInputStream inputStream, int cohortsCount) throws EOFException, IOException, BinaryInteractionFileException {

		final BinaryInteractionCohort[] cohorts = new BinaryInteractionCohort[cohortsCount];
		final HashSet<String> cohortNames = new HashSet<String>(cohortsCount);

		for (int i = 0; i < cohortsCount; ++i) {

			String name = readString(inputStream);

			if (!cohortNames.add(name)) {
				throw new BinaryInteractionFileException("Cohort names must be unique. " + name + " has been found multiple times");
			}

			cohorts[i] = new BinaryInteractionCohort(name, inputStream.readInt());

		}
		return cohorts;
	}

	private static String[] readStringArray(DataInputStream inputStream, int size, String name) throws BinaryInteractionFileException, EOFException, IOException {

		final String[] array = new String[size];
		final HashSet<String> names = new HashSet<String>(size);

		for (int i = 0; i < size; ++i) {

			array[i] = readString(inputStream);
			if (!names.add(array[i])) {
				throw new BinaryInteractionFileException(name + " entries must be unique. " + array[i] + " has been found multiple times");
			}

		}
		return array;

	}

	private static Allele[] readAlleles(DataInputStream inputStream, int size) throws BinaryInteractionFileException, EOFException, IOException {

		final Allele[] alleles = new Allele[size];
		final HashSet<String> allelesCheck = new HashSet<String>(size);

		for (int i = 0; i < size; ++i) {

			String alleleString = readString(inputStream);
			alleles[i] = Allele.create(alleleString);
			if (!allelesCheck.add(alleleString)) {
				throw new BinaryInteractionFileException("Allele entries must be unique. " + alleleString + " has been found multiple times");
			}

		}
		return alleles;

	}

	private static BinaryInteractionVariant[] readVariants(DataInputStream inputStream, int variantCount, String[] chrDictionary, Allele[] alleleDictionary) throws EOFException, IOException, BinaryInteractionFileException {

		BinaryInteractionVariant[] variants = new BinaryInteractionVariant[variantCount];

		for (int i = 0; i < variantCount; ++i) {

			String name = readString(inputStream);
			String chr = chrDictionary[inputStream.readInt()];
			int pos = inputStream.readInt();
			Allele refAllele = alleleDictionary[inputStream.readInt()];
			Allele altAllele = alleleDictionary[inputStream.readInt()];
			int geneCount = inputStream.readInt();
			int[] genes = readIntArray(inputStream, geneCount);

			variants[i] = new BinaryInteractionVariantStatic(name, chr, pos, refAllele, altAllele, genes);

		}

		return variants;

	}

	private static BinaryInteractionGene[] readGenes(DataInputStream inputStream, int geneCount, String[] chrDictionary) throws EOFException, IOException {

		BinaryInteractionGene[] genes = new BinaryInteractionGene[geneCount];

		for (int i = 0; i < geneCount; ++i) {

			String name = readString(inputStream);
			String chr = chrDictionary[inputStream.readInt()];
			int start = inputStream.readInt();
			int end = inputStream.readInt();
			int variantCount = inputStream.readInt();
			int[] variants = readIntArray(inputStream, variantCount);

			genes[i] = new BinaryInteractionGeneStatic(name, chr, start, end, variants);

		}

		return genes;

	}

	private static int[][] readCovariatesData(DataInputStream inputStream, int totalSnpGeneCombinations, long totalInteractions) throws BinaryInteractionFileException, EOFException, IOException {

		int[][] testedCovariates = new int[totalSnpGeneCombinations][];

		long interactionSumCovariates = 0;

		for (int i = 0; i < totalSnpGeneCombinations; ++i) {
			int covariatesCount = inputStream.readInt();
			testedCovariates[i] = readIntArray(inputStream, covariatesCount);
			interactionSumCovariates += testedCovariates[i].length;
		}

		if (interactionSumCovariates != totalInteractions) {
			throw new BinaryInteractionFileException("Interactions in header differs from interactions in covariate block. Header: " + totalInteractions + " covariate sum: " + interactionSumCovariates);
		}

		return testedCovariates;

	}

	/**
	 * Will first read int to determine length. Will then read length chars and
	 * convert to string.
	 *
	 * @param inputStream
	 * @return
	 */
	private static String readString(DataInputStream inputStream) throws EOFException, IOException {

		int lenght = inputStream.readInt();
		char[] chars = new char[lenght];

		for (int i = 0; i < lenght; ++i) {
			chars[i] = inputStream.readChar();
		}

		return new String(chars);

	}

	private static int[] readIntArray(DataInputStream inputStream, int arraySize) throws EOFException, IOException {

		int[] array = new int[arraySize];
		for (int i = 0; i < arraySize; ++i) {
			array[i] = inputStream.readInt();
		}
		return array;

	}

	protected static long calculateSizeNormalQtlBlock(int cohorts, boolean metaAnalysis) {
		long size = (cohorts * 12);
		if (metaAnalysis) {
			size += 8;
		}
		return size;
	}

	protected static long calculateSizeInteractionResultBlock(int cohorts, boolean flippedZscoreStored, boolean metaAnalysis) {
		long size = cohorts * 36;
		if (flippedZscoreStored) {
			size += (cohorts * 8);
			if (metaAnalysis) {
				size += 8;
			}
		}
		if (metaAnalysis) {
			size += 24;
		}
		return size;
	}

	private static int getTotalVariantGeneCombinations(BinaryInteractionVariant[] variants) throws BinaryInteractionFileException {
		long totalSnpGeneCombinationsLong = 0;
		for (BinaryInteractionVariant variant : variants) {
			totalSnpGeneCombinationsLong += variant.getGeneCount();
		}

		if (totalSnpGeneCombinationsLong > Integer.MAX_VALUE) {
			throw new BinaryInteractionFileException("Cannot handle variant-gene combinations bigger than 2^31-1");
		}

		return (int) totalSnpGeneCombinationsLong;
	}

	public File getInteractionFile() {
		return interactionFile;
	}

	public boolean isReadOnly() {
		return readOnly;
	}

	public long getCreationDataEpoch() {
		return timeStamp;
	}

	public String getCreationDataTimeString() {
		return DEFAULT_DATE_FORMAT.format(new Date(timeStamp * 1000));
	}

	public boolean isMetaAnalysis() {
		return metaAnalysis;
	}

	public boolean isNormalQtlStored() {
		return normalQtlStored;
	}

	public boolean isFlippedZscoreStored() {
		return flippedZscoreStored;
	}

	public String getFileDescription() {
		return fileDescription;
	}

	public boolean areAllCovariatesTestedForAllVariantGenes() {
		return allCovariants;
	}

	public long getTotalNumberInteractions() {
		return interactions;
	}

	public List<BinaryInteractionCohort> getCohorts() {
		return Collections.unmodifiableList(Arrays.asList(cohorts));
	}

	public List<BinaryInteractionGene> getGenes() {
		return Collections.unmodifiableList(Arrays.asList(genes));
	}
	
	public BinaryInteractionGene getGene(String name) throws BinaryInteractionFileException{
		int index = genesMap.get(name);
		if(index == NO_ENTRY_INT_MAP){
			throw new BinaryInteractionFileException("Gene not found: " + name);
		}
		return genes[index];
	}
	
	public BinaryInteractionGene getGene(int pointer) throws BinaryInteractionFileException{
		return genes[pointer];
	}


	public List<BinaryInteractionVariant> getVariants() {
		return Collections.unmodifiableList(Arrays.asList(variants));
	}
	
	public BinaryInteractionVariant getVariant(String name) throws BinaryInteractionFileException{
		int index = variantMap.get(name);
		if(index == NO_ENTRY_INT_MAP){
			throw new BinaryInteractionFileException("Variant not found: " + name);
		}
		return variants[index];
	}
	
	public BinaryInteractionVariant getVariant(int pointer) throws BinaryInteractionFileException{
		return variants[pointer];
	}

	public List<String> getCovariates() {
		return Collections.unmodifiableList(Arrays.asList(covariates));
	}

	public void finalizeWriting() throws IOException, BinaryInteractionFileException {
		if (!readOnly) {

			if (qtlBufferWriting) {
				writeQtlBuffer();
			}

			if (interactionBufferWriting) {
				writeInteractionBuffer();
			}

			randomAccess.seek(POINTER_TO_CLOSED_BOOLEAN);
			randomAccess.writeBoolean(true);

			this.close();
			readOnly = true;
			try {
				this.open();
			} catch (IOException ex) {
				throw new RuntimeException(ex);
			}
		}
	}

	public long getQtlZscoresSet() {
		return qtlZscoresSet;
	}

	public long getInteractionZscoresSet() {
		return interactionZscoresSet;
	}

	public long getQtlZscoresRead() {
		return qtlZscoresRead;
	}

	public long getInteractionZscoresRead() {
		return interactionZscoresRead;
	}

	public long getInteractionWriteBufferFlushed() {
		return interactionWriteBufferFlushed;
	}

	public long getQtlWriteBufferFlushed() {
		return qtlWriteBufferFlushed;
	}

	public long getInteractionReadBufferLoaded() {
		return interactionReadBufferLoaded;
	}

	public long getQtlReadBufferLoaded() {
		return qtlReadBufferLoaded;
	}

	@Override
	public void close() throws IOException {
		channel.close();
		randomAccess.close();
	}

	private void open() throws FileNotFoundException, IOException {
		randomAccess = new RandomAccessFile(interactionFile, readOnly ? "r" : "rw");
		if (!readOnly) {
			randomAccess.seek(POINTER_TO_CLOSED_BOOLEAN);
			randomAccess.writeBoolean(false);
		}
		channel = randomAccess.getChannel();
	}

	private long getQtlPointer(String variantName, String geneName) throws BinaryInteractionFileException {

		int variantIndex = variantMap.get(variantName);
		int geneIndex = genesMap.get(geneName);

		if (variantIndex < 0) {
			throw new BinaryInteractionFileException("Variant not found: " + variantName);
		}

		if (geneIndex < 0) {
			throw new BinaryInteractionFileException("Gene not found: " + geneName);
		}

		int geneIndexInVariant = variants[variantIndex].getIndexOfGenePointer(geneIndex);

		if (geneIndexInVariant < 0) {
			throw new BinaryInteractionFileException("Cannot find QTL for: " + variantName + "-" + geneName);
		}

		return startQtlBlock + ((cummulativeGeneCountUpToVariant[variantIndex] + geneIndexInVariant) * sizeQtlBlock);

	}

	private long getInteractionPointer(String variantName, String geneName, String covariateName) throws BinaryInteractionFileException {

		int variantIndex = variantMap.get(variantName);
		int geneIndex = genesMap.get(geneName);
		int covariateIndex = covariatesMap.get(covariateName);

		if (variantIndex < 0) {
			throw new BinaryInteractionFileException("Variant not found: " + variantName);
		}

		if (geneIndex < 0) {
			throw new BinaryInteractionFileException("Gene not found: " + geneName);
		}

		if (covariateIndex < 0) {
			throw new BinaryInteractionFileException("Covariate not found: " + covariateName);
		}

		int geneIndexInVariant = variants[variantIndex].getIndexOfGenePointer(geneIndex);

		if (geneIndexInVariant < 0) {
			throw new BinaryInteractionFileException("Cannot find variant gene combination for: " + variantName + "-" + geneName);
		}

		int variantGeneIndex = cummulativeGeneCountUpToVariant[variantIndex] + geneIndexInVariant;

		int variantGeneCovariateIndex;
		if (allCovariants) {
			variantGeneCovariateIndex = covariateIndex;
		} else {
			variantGeneCovariateIndex = Arrays.binarySearch(covariatesTested[variantGeneIndex], covariateIndex);
			if (variantGeneCovariateIndex < 0) {
				throw new BinaryInteractionFileException("Cannot find covariate " + covariateName + " for: " + variantName + "-" + geneName);
			}
		}

		return startInteractionBlock + ((cummalitiveInteractionCountUptoVariantGene[variantGeneIndex] + variantGeneCovariateIndex) * sizeInteractionBlock);

	}

	public BinaryInteractionQtlZscores readQtlResults(String variantName, String geneName) throws BinaryInteractionFileException, IOException {

		if(!normalQtlStored){
			throw new BinaryInteractionFileException("This file does not store normal QTL restuls");
		}
		
		//Check will be done in get pointer
		long qtlPointer = getQtlPointer(variantName, geneName);

		setQtlBuffer(qtlPointer, false);

		double[] zscore = new double[cohorts.length];
		for (int i = 0; i < cohorts.length; i++) {
			zscore[i] = qtlBuffer.getDouble();
		}

		int[] sampleCounts = new int[cohorts.length];
		for (int i = 0; i < cohorts.length; i++) {
			sampleCounts[i] = qtlBuffer.getInt();
		}
		
		++qtlZscoresRead;

		if (metaAnalysis) {
			double metaZscore = qtlBuffer.getDouble();
			return new BinaryInteractionQtlZscores(zscore, sampleCounts, metaZscore);
		} else {
			return new BinaryInteractionQtlZscores(zscore, sampleCounts);
		}

	}

	public void setQtlResults(String variantName, String geneName, BinaryInteractionQtlZscores zscores) throws BinaryInteractionFileException, IOException {
		
		if(!normalQtlStored){
			throw new BinaryInteractionFileException("This file does not store normal QTL restuls");
		}

		if (zscores.getZscores().length != cohorts.length) {
			throw new BinaryInteractionFileException("Error setting qtl " + variantName + "-" + geneName + " expected " + cohorts.length + " but found " + zscores.getZscores().length + " Z-scores");
		}

		if (zscores.getSampleCounts().length != cohorts.length) {
			throw new BinaryInteractionFileException("Error setting qtl " + variantName + "-" + geneName + " expected " + cohorts.length + " but found " + zscores.getSampleCounts().length + " samples counts");
		}

		//Check will be done in get pointer
		long qtlPointer = getQtlPointer(variantName, geneName);

		setQtlBuffer(qtlPointer, true);

		for (int i = 0; i < cohorts.length; i++) {
			qtlBuffer.putDouble(zscores.getZscores()[i]);
		}

		for (int i = 0; i < cohorts.length; i++) {
			qtlBuffer.putInt(zscores.getSampleCounts()[i]);
		}

		if (metaAnalysis) {
			qtlBuffer.putDouble(zscores.getMetaZscore());
		}
		
		++qtlZscoresSet;

	}

	/**
	 * Set the QTL buffer at the position to start reading from qtlPointer.
	 *
	 * @param qtlPointer
	 * @param writing
	 * @throws BinaryInteractionFileException
	 * @throws IOException
	 */
	private void setQtlBuffer(long qtlPointer, boolean writing) throws BinaryInteractionFileException, IOException {

		if (writing) {

			if (qtlBufferWriting && qtlPointer == qtlBufferStart + qtlBuffer.position() && qtlBuffer.remaining() >= sizeQtlBlock) {
				//Current write buffer;
				//return;
			} else {
				if (qtlBufferWriting) {
					writeQtlBuffer();
				}
				qtlBuffer.clear();
				qtlBufferWriting = true;
				qtlBufferStart = qtlPointer;
				//return;
			}

		} else { //reading 

			if (qtlBufferWriting) {
				writeQtlBuffer();
			}

			if (qtlPointer >= qtlBufferStart && (qtlPointer + sizeQtlBlock) <= qtlBufferStart + qtlBuffer.limit()) {
				int positionInBuffer = (int) (qtlPointer - qtlBufferStart);
				qtlBuffer.position(positionInBuffer);
				//return;
			} else {
				channel.position(qtlPointer);
				qtlBuffer.clear();
				channel.read(qtlBuffer);
				qtlBuffer.flip();
				qtlBufferStart = qtlPointer;
				++qtlReadBufferLoaded;
				//return;
			}

		}

	}

	private void writeQtlBuffer() throws BinaryInteractionFileException, IOException {
		if (readOnly) {
			throw new BinaryInteractionFileException("Interaction file is in read only mode");
		}
		if (!qtlBufferWriting) {
			throw new BinaryInteractionFileException("Interaction file has no QTLs to write");
		}

		channel.position(qtlBufferStart);
		qtlBuffer.flip();
		channel.write(qtlBuffer);
		qtlBufferWriting = false;
		++qtlWriteBufferFlushed;
	}
	
	public BinaryInteractionZscores readInteractionResults(String variantName, String geneName, String covariateName) throws BinaryInteractionFileException, IOException {
		//Check will be done in get pointer
		long interactionPointer = getInteractionPointer(variantName, geneName, covariateName);
		return readInteractionResults(interactionPointer);
	}
	
	public BinaryInteractionQueryResult readVariantGeneCovariateResults(String variantName, String geneName, String covariateName) throws BinaryInteractionFileException, IOException{
		
		BinaryInteractionQtlZscores qtlZscores;
		if(isNormalQtlStored()){
			qtlZscores = readQtlResults(variantName, geneName);
		} else {
			qtlZscores = null;
		}
		
		BinaryInteractionZscores interactionRestuls = readInteractionResults(variantName, geneName, covariateName);
		
		return new BinaryInteractionQueryResult(variantName, geneName, covariateName, qtlZscores, interactionRestuls);
		
	}
	
	public Iterator<BinaryInteractionQueryResult> readVariantGeneResults(String variantName, String geneName) throws BinaryInteractionFileException, IOException{
		
		int variantIndex = variantMap.get(variantName);
		int geneIndex = genesMap.get(geneName);

		if (variantIndex < 0) {
			throw new BinaryInteractionFileException("Variant not found: " + variantName);
		}

		if (geneIndex < 0) {
			throw new BinaryInteractionFileException("Gene not found: " + geneName);
		}

		int geneIndexInVariant = variants[variantIndex].getIndexOfGenePointer(geneIndex);

		if (geneIndexInVariant < 0) {
			throw new BinaryInteractionFileException("Cannot find variant gene combination for: " + variantName + "-" + geneName);
		}

		int variantGeneIndex = cummulativeGeneCountUpToVariant[variantIndex] + geneIndexInVariant;
		
		BinaryInteractionQtlZscores qtlZscore;
		if(isNormalQtlStored()){
			qtlZscore = readQtlResults(variantName, geneName);
		} else {
			qtlZscore = null;
		}
		
		long startVariantGeneBlock = startInteractionBlock + (cummalitiveInteractionCountUptoVariantGene[variantGeneIndex] * sizeInteractionBlock);

		return new BinaryInteractionCovariateIterator(variantName, geneName, allCovariants ? null : covariatesTested[variantGeneIndex] , this, startVariantGeneBlock, qtlZscore);

		
	}

	protected BinaryInteractionZscores readInteractionResults(long interactionPointer) throws BinaryInteractionFileException, IOException {

		setInteactionBuffer(interactionPointer, false);

		++interactionZscoresRead;

		final int[] samplesInteractionCohort = readIntArrayFromInteractionBuffer(cohorts.length);
		final double[] zscoreSnpCohort = readDoubleArrayFromInteractionBuffer(cohorts.length);
		final double[] zscoreCovariateCohort = readDoubleArrayFromInteractionBuffer(cohorts.length);
		final double[] zscoreInteractionCohort = readDoubleArrayFromInteractionBuffer(cohorts.length);
		final double[] rSquaredCohort = readDoubleArrayFromInteractionBuffer(cohorts.length);
		final double[] zscoreInteractionFlippedCohort;
		if (flippedZscoreStored) {
			zscoreInteractionFlippedCohort = readDoubleArrayFromInteractionBuffer(cohorts.length);
		} else {
			zscoreInteractionFlippedCohort = new double[cohorts.length];
			Arrays.fill(zscoreInteractionFlippedCohort, Double.NaN);
		}
		if (metaAnalysis) {
			final double zscoreSnpMeta = interactionBuffer.getDouble();
			final double zscoreCovariateMeta = interactionBuffer.getDouble();
			final double zscoreInteractionMeta = interactionBuffer.getDouble();
			if (flippedZscoreStored) {
				final double zscoreInteractionFlippedMeta = interactionBuffer.getDouble();
				return new BinaryInteractionZscores(samplesInteractionCohort, zscoreSnpCohort, zscoreCovariateCohort, zscoreInteractionCohort, rSquaredCohort, zscoreInteractionFlippedCohort, zscoreSnpMeta, zscoreCovariateMeta, zscoreInteractionMeta, zscoreInteractionFlippedMeta);
			} else {
				return new BinaryInteractionZscores(samplesInteractionCohort, zscoreSnpCohort, zscoreCovariateCohort, zscoreInteractionCohort, rSquaredCohort, zscoreSnpMeta, zscoreCovariateMeta, zscoreInteractionMeta);
			}
		} else {
			return new BinaryInteractionZscores(samplesInteractionCohort, zscoreSnpCohort, zscoreCovariateCohort, zscoreInteractionCohort, rSquaredCohort, zscoreInteractionFlippedCohort);
		}
	}

	public void setInteractionResults(String variantName, String geneName, String covariateName, BinaryInteractionZscores zscores) throws BinaryInteractionFileException, IOException {

		if (zscores.getSamplesInteractionCohort().length != cohorts.length) {
			throw new BinaryInteractionFileException("Error setting interaction " + variantName + "-" + geneName + " expected " + cohorts.length + " but found " + zscores.getSamplesInteractionCohort().length + " cohorts");
		}

		//Check will be done in get pointer
		long interactionPointer = getInteractionPointer(variantName, geneName, covariateName);
		
		setInteactionBuffer(interactionPointer, true);

		writeIntArrayToInteractionBuffer(zscores.getSamplesInteractionCohort());
		writeDoubleArrayToInteractionBuffer(zscores.getZscoreSnpCohort());
		writeDoubleArrayToInteractionBuffer(zscores.getZscoreCovariateCohort());
		writeDoubleArrayToInteractionBuffer(zscores.getZscoreInteractionCohort());
		writeDoubleArrayToInteractionBuffer(zscores.getrSquaredCohort());
		if (flippedZscoreStored) {
			writeDoubleArrayToInteractionBuffer(zscores.getZscoreInteractionFlippedCohort());
		}

		if (metaAnalysis) {
			interactionBuffer.putDouble(zscores.getZscoreSnpMeta());
			interactionBuffer.putDouble(zscores.getZscoreCovariateMeta());
			interactionBuffer.putDouble(zscores.getZscoreInteractionMeta());
			if (flippedZscoreStored) {
				interactionBuffer.putDouble(zscores.getZscoreInteractionFlippedMeta());
			}
		}
		
		++interactionZscoresSet;

	}

	/**
	 * Set the interaction buffer at the position to start reading or writing
	 * from interactionPointer.
	 *
	 * @param interactionPointer
	 * @param writing
	 * @throws BinaryInteractionFileException
	 * @throws IOException
	 */
	private void setInteactionBuffer(long interactionPointer, boolean writing) throws BinaryInteractionFileException, IOException {

		if (writing) {

			if (interactionBufferWriting && interactionPointer == interactionBufferStart + interactionBuffer.position() && interactionBuffer.remaining() >= sizeInteractionBlock) {
				//Current write buffer;
				//return;
			} else {
				if (interactionBufferWriting) {
					writeInteractionBuffer();
				}
				interactionBuffer.clear();
				interactionBufferWriting = true;
				interactionBufferStart = interactionPointer;
				//return;
			}

		} else { //reading 

			if (interactionBufferWriting) {
				writeInteractionBuffer();
			}

			if (interactionPointer >= interactionBufferStart && (interactionPointer + sizeInteractionBlock) <= interactionBufferStart + interactionBuffer.limit()) {
				int positionInBuffer = (int) (interactionPointer - interactionBufferStart);
				interactionBuffer.position(positionInBuffer);
				//return;
			} else {
				channel.position(interactionPointer);
				interactionBuffer.clear();
				channel.read(interactionBuffer);
				interactionBuffer.flip();
				interactionBufferStart = interactionPointer;
				++interactionReadBufferLoaded;
				//return;
			}

		}



	}

	private void writeInteractionBuffer() throws BinaryInteractionFileException, IOException {
		if (readOnly) {
			throw new BinaryInteractionFileException("Interaction file is in read only mode");
		}
		if (!interactionBufferWriting) {
			throw new BinaryInteractionFileException("Interaction file has no interactions to write");
		}

		channel.position(interactionBufferStart);
		interactionBuffer.flip();
		channel.write(interactionBuffer);
		interactionBufferWriting = false;
		++interactionWriteBufferFlushed;

	}

	private double[] readDoubleArrayFromInteractionBuffer(int length) {
		if (length == 0) {
			return BinaryInteractionZscores.emptyDoubleArray;
		}
		double[] array = new double[length];
		for (int i = 0; i < length; ++i) {
			array[i] = interactionBuffer.getDouble();
		}
		return array;
	}

	private int[] readIntArrayFromInteractionBuffer(int length) {
		if (length == 0) {
			return BinaryInteractionZscores.emptyIntArray;
		}
		int[] array = new int[length];
		for (int i = 0; i < length; ++i) {
			array[i] = interactionBuffer.getInt();
		}
		return array;
	}

	private void writeDoubleArrayToInteractionBuffer(double[] array) {
		for (int i = 0; i < array.length; ++i) {
			interactionBuffer.putDouble(array[i]);
		}
	}

	private void writeIntArrayToInteractionBuffer(int[] array) {
		for (int i = 0; i < array.length; ++i) {
			interactionBuffer.putInt(array[i]);
		}
	}
	
	public int getGeneCount(){
		return genes.length;
	}
	
	public int getCovariateCount(){
		return covariates.length;
	}
	
	public int getVariantCount(){
		return variants.length;
	}
	
	public int getCohortCount(){
		return cohorts.length;
	}
	
	public int getVariantGeneCombinations(){
		return cummalitiveInteractionCountUptoVariantGene.length - 1;
	}
	
	public boolean containsGene(String geneName){
		return genesMap.containsKey(geneName);
	}
	
	public boolean containsVariant(String variantName){
		return variantMap.containsKey(variantName);
	}
	
	public boolean containsCovariant(String covariateName){
		return covariatesMap.containsKey(covariateName);
	}
	
	public boolean containsVariantGene(String variantName, String geneName){
		
		int variantIndex = variantMap.get(variantName);
		int geneIndex = genesMap.get(geneName);

		if (variantIndex < 0) {
			return false;
		}

		if (geneIndex < 0) {
			return false;
		}

		int geneIndexInVariant = variants[variantIndex].getIndexOfGenePointer(geneIndex);

		if (geneIndexInVariant < 0) {
			return false;
		} else {
			return true;
		}
		
	}
	
	public boolean containsInteraction(String variantName, String geneName, String covariateName) throws BinaryInteractionFileException {

		int variantIndex = variantMap.get(variantName);
		int geneIndex = genesMap.get(geneName);
		int covariateIndex = covariatesMap.get(covariateName);

		if (variantIndex < 0) {
			return false;
		}

		if (geneIndex < 0) {
			return false;
		}

		if (covariateIndex < 0) {
			return false;
		}

		int geneIndexInVariant = variants[variantIndex].getIndexOfGenePointer(geneIndex);

		if (geneIndexInVariant < 0) {
			return false;
		}

		if (allCovariants) {
			return true;
		} else {
			int variantGeneIndex = cummulativeGeneCountUpToVariant[variantIndex] + geneIndexInVariant;
			int variantGeneCovariateIndex = Arrays.binarySearch(covariatesTested[variantGeneIndex], covariateIndex);
			if (variantGeneCovariateIndex < 0) {
				return false;
			} else {
				return true;
			}
		}
	}
	
	private void createVariantChrPosMap(){
		variantsByPos = new ChrPosMap<BinaryInteractionVariant>();
		for(BinaryInteractionVariant variant : variants){
			variantsByPos.put(variant.getChr(), variant.getPos(), variant);
		}
	}
	
	/**
	 * 
	 * @param chr
	 * @param pos
	 * @return null if not found
	 */
	public BinaryInteractionVariant getVariant(String chr, int pos){
		if(variantsByPos == null){
			createVariantChrPosMap();
		}
		return variantsByPos.get(chr, pos);
	}
}
