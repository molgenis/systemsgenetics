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
import java.text.SimpleDateFormat;
import java.util.Arrays;
import java.util.Collections;
import java.util.Date;
import java.util.HashSet;
import java.util.List;
import org.molgenis.genotype.Allele;
import umcg.genetica.io.binInteraction.gene.BinaryInteractionGeneStatic;
import umcg.genetica.io.binInteraction.variant.BinaryInteractionVariant;

/**
 *
 * @author Patrick Deelen
 */
public class BinaryInteractionFile implements Closeable {

	protected static final byte MAGIC_1 = 81;
	protected static final byte MAGIC_2 = 73;
	private static final long POINTER_TO_CLOSED_BOOLEAN = 4;
	private static final byte MAJOR_VERSION = 1;
	private static final byte MINOR_VERSION = 0;
	private static final int NO_ENTRY_INT_MAP = -1;
	private static final SimpleDateFormat DEFAULT_DATE_FORMAT = new java.text.SimpleDateFormat("dd-MM-yyyy HH:mm:ss");
	private final File interactionFile;
	private boolean readOnly;
	private final BinaryInteractionCohort[] cohorts;
	private final BinaryInteractionGene[] genes;
	private final BinaryInteractionVariant[] variants;
	private final String[] covariats;
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
	private final long sizeInteractionBlock;
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

	public BinaryInteractionFile(File interactionFile, boolean readOnly, BinaryInteractionCohort[] cohorts, BinaryInteractionGene[] genes, BinaryInteractionVariant[] variants, String[] covariats, int[][] covariatesTested, long timeStamp, boolean allCovariants, boolean metaAnalysis, boolean normalQtlStored, boolean flippedZscoreStored, String fileDescription, long interactions, long startQtlBlock, long startInteractionBlock) throws BinaryInteractionFileException, FileNotFoundException, IOException {
		this.interactionFile = interactionFile;
		this.readOnly = readOnly;
		this.cohorts = cohorts;
		this.genes = genes;
		this.variants = variants;
		this.covariats = covariats;
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
					cummalitiveInteractionCountUptoVariantGene[i] = cummalitiveInteractionCountUptoVariantGene[i - 1] + (allCovariants ? covariats.length : this.covariatesTested[i - 1].length);
				}
			}
		}

		if (cummalitiveInteractionCountUptoVariantGene[cummulativeGeneCountUpToVariant[variants.length]] != interactions) {
			System.out.println(cummalitiveInteractionCountUptoVariantGene[cummulativeGeneCountUpToVariant[variants.length]]);
			System.out.println(interactions);
			throw new BinaryInteractionFileException("Something went wrong");
		}

		variantMap = new TObjectIntHashMap<String>(variants.length, 0.75f, -1);
		genesMap = new TObjectIntHashMap<String>(genes.length, 0.75f, -1);
		covariatesMap = new TObjectIntHashMap<String>(covariats.length, 0.75f, -1);

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

		for (int i = 0; i < covariats.length; ++i) {
			if (covariatesMap.put(covariats[i], i) != NO_ENTRY_INT_MAP) {
				throw new BinaryInteractionFileException("Cannot store the same covariate twice (" + covariats[i] + ")");
			}
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
			final int covariatsCount = inputStream.readInt();

			long totalInteractions = inputStream.readLong();
			builder.setInteractions(totalInteractions);

			builder.setCohorts(readCohorts(inputStream, chortsCount));
			final String[] chrDictionary = readStringArray(inputStream, chrsCount, "Chromosomes");
			final Allele[] alleleDictionary = readAlleles(inputStream, allelesCount);

			final BinaryInteractionVariant[] variants = readVariants(inputStream, variantCount, chrDictionary, alleleDictionary);

			builder.setVariants(variants);
			builder.setGenes(readGenes(inputStream, genesCount, chrDictionary));
			builder.setCovariats(readStringArray(inputStream, covariatsCount, "Chromosomes"));

			final int totalVariantGeneCombinations = getTotalVariantGeneCombinations(variants);

			if (!allCovariants) {
				readCovariantsData(inputStream, totalVariantGeneCombinations, totalInteractions);
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

			final long sizeInteractionBlock = calculateSizeInteractionResultBlock(chortsCount, flippedZscoreStored, metaAnalysis);

			if (startData + sizeNormalQtlSection + sizeInteractionBlock != interactionFile.length()) {
				throw new BinaryInteractionFileException("Incorrect file size. Expected: " + (startData + sizeNormalQtlSection + sizeInteractionBlock) + " found: " + interactionFile.length() + " diff: " + (startData + sizeNormalQtlSection + sizeInteractionBlock - interactionFile.length()));
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

	private static int[][] readCovariantsData(DataInputStream inputStream, int totalSnpGeneCombinations, long totalInteractions) throws BinaryInteractionFileException, EOFException, IOException {

		int[][] testedCovariats = new int[totalSnpGeneCombinations][];

		long interactionSumCovariats = 0;

		for (int i = 0; i < totalSnpGeneCombinations; ++i) {
			int covariatsCount = inputStream.readInt();
			testedCovariats[i] = readIntArray(inputStream, covariatsCount);
			interactionSumCovariats += testedCovariats[i].length;
		}

		if (interactionSumCovariats != totalInteractions) {
			throw new BinaryInteractionFileException("Interactions in header differs from interactions in covariate block. Header: " + totalInteractions + " covariate sum: " + interactionSumCovariats);
		}

		return testedCovariats;

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

	public boolean areAllCovariatsTestedForAllVariantGenes() {
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

	public List<BinaryInteractionVariant> getVariants() {
		return Collections.unmodifiableList(Arrays.asList(variants));
	}

	public List<String> getCovariats() {
		return Collections.unmodifiableList(Arrays.asList(covariats));
	}

	public void makeReadOnly() {
		if (!readOnly) {
			this.close();
			readOnly = true;
			try {
				this.open();
			} catch (IOException ex) {
				throw new RuntimeException(ex);
			}
		}
	}

	@Override
	public void close() {
		try {
			randomAccess.seek(POINTER_TO_CLOSED_BOOLEAN);
			randomAccess.writeBoolean(true);
			randomAccess.close();
		} catch (IOException ex) {
			throw new RuntimeException(ex);
		}
	}

	private void open() throws FileNotFoundException, IOException {
		randomAccess = new RandomAccessFile(interactionFile, readOnly ? "r" : "rw");
		if (!readOnly) {
			randomAccess.seek(POINTER_TO_CLOSED_BOOLEAN);
			randomAccess.writeBoolean(false);
		}
	}

	public long getQtlPointer(String variantName, String geneName) throws BinaryInteractionFileException {

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

	public long getInteractionPointer(String variantName, String geneName, String covariateName) throws BinaryInteractionFileException {

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
		if(allCovariants){
			variantGeneCovariateIndex = covariateIndex;
		} else {
			variantGeneCovariateIndex = Arrays.binarySearch(covariatesTested[variantGeneIndex], covariateIndex);
			if(variantGeneCovariateIndex < 0 ){
				throw new BinaryInteractionFileException("Cannot find covariate " + covariateName + " for: " + variantName + "-" + geneName);
			}
		}

		return startInteractionBlock + (cummalitiveInteractionCountUptoVariantGene[variantGeneIndex] + variantGeneCovariateIndex * sizeInteractionBlock);
		
	}
}
