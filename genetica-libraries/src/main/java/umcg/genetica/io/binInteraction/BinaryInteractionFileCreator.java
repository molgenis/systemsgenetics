package umcg.genetica.io.binInteraction;

import gnu.trove.map.hash.TObjectIntHashMap;
import java.io.BufferedOutputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.Arrays;
import java.util.HashSet;
import org.molgenis.genotype.Allele;
import static umcg.genetica.io.binInteraction.BinaryInteractionFile.calculateInteractionResultBlock;
import static umcg.genetica.io.binInteraction.BinaryInteractionFile.calculateSizeNormalQtlBlock;
import umcg.genetica.io.binInteraction.gene.BinaryInteractionGene;
import umcg.genetica.io.binInteraction.gene.BinaryInteractionGeneCreator;
import umcg.genetica.io.binInteraction.variant.BinaryInteractionVariantCreator;

/**
 *
 * @author Patrick Deelen
 */
public class BinaryInteractionFileCreator {

	private static final int NO_ENTRY_INT_MAP = -1;
	private final BinaryInteractionVariantCreator[] variants;
	private final BinaryInteractionGeneCreator[] genes;
	private final BinaryInteractionCohort[] cohorts;
	private final String[] covariats;
	private final boolean allCovariants;
	private final boolean metaAnalysis;
	private final boolean normalQtlStored;
	private final boolean flippedZscoreStored;
	private String description = "";
	private int[][] covariatesTested;
	private long interactions = 0;
	private long startQtlBlock;
	private long startInteractionBlock;
	private boolean startedAddingCovariates = false;
	private int countVariantGeneCombinations = 0;
	private final TObjectIntHashMap<String> variantMap;
	private final TObjectIntHashMap<String> genesMap;
	private final TObjectIntHashMap<String> covariatesMap;
	private int[] variantCummulativeGeneCounts;
	private int variantGenesWithCovariatesAdded = 0;
	private boolean sortedIndices = false;
	private final File file;
	private boolean created = false;

	/**
	 * 
	 * @param file
	 * @param variants
	 * @param genes
	 * @param cohorts
	 * @param covariats
	 * @param allCovariants
	 * @param metaAnalysis
	 * @param normalQtlStored
	 * @param flippedZscoreStored
	 * @throws BinaryInteractionFileException 
	 */
	public BinaryInteractionFileCreator(File file, BinaryInteractionVariantCreator[] variants, BinaryInteractionGeneCreator[] genes, BinaryInteractionCohort[] cohorts, String[] covariats, boolean allCovariants, boolean metaAnalysis, boolean normalQtlStored, boolean flippedZscoreStored) throws BinaryInteractionFileException {
		this.file = file;
		this.variants = variants;
		this.genes = genes;
		this.cohorts = cohorts;
		this.covariats = covariats;
		this.allCovariants = allCovariants;
		this.metaAnalysis = metaAnalysis;
		this.normalQtlStored = normalQtlStored;
		this.flippedZscoreStored = flippedZscoreStored;

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

		HashSet<String> tmp = new HashSet<String>(cohorts.length);
		for (int i = 0; i < cohorts.length; ++i) {
			if (!tmp.add(cohorts[i].getName())) {
				throw new BinaryInteractionFileException("Cannot store the same cohort twice (" + cohorts[i] + ")");
			}
		}

	}

	public synchronized void addTestedVariantGene(String variantName, String geneName) throws BinaryInteractionFileException {

		sortedIndices = false;

		if (created) {
			throw new BinaryInteractionFileException("You already created this file. Adding variant-gene combinations no longer possible");
		}

		if (startedAddingCovariates) {
			throw new BinaryInteractionFileException("All variant-gene combinations must be added before setting covariates. (sorry)");
		}

		int variantIndex = variantMap.get(variantName);
		int geneIndex = genesMap.get(geneName);

		if (variantIndex == NO_ENTRY_INT_MAP) {
			throw new BinaryInteractionFileException("Unable to add variant-gene combination. Variant " + variantName + " not found");
		}

		if (geneIndex == NO_ENTRY_INT_MAP) {
			throw new BinaryInteractionFileException("Unable to add variant-gene combination. Gene " + geneName + " not found");
		}

		BinaryInteractionVariantCreator variant = variants[variantIndex];
		BinaryInteractionGeneCreator gene = genes[geneIndex];

		//Below can not be done using the search method since the pointers are not yet sorted
		for (int existinVariantGenePointer : variant.getGenePointers()) {
			if (existinVariantGenePointer == geneIndex) {
				throw new BinaryInteractionFileException("Cannot add a variant-gene combination twice: " + variantName + "-" + geneName);
			}
		}
		++countVariantGeneCombinations;
		variant.addGene(geneIndex);
		gene.addVariant(variantIndex);

	}

	public synchronized void addTestedInteraction(String variantName, String geneName, String[] covariateNames) throws BinaryInteractionFileException {

		if (created) {
			throw new BinaryInteractionFileException("You already created this file. Adding covariates no longer possible");
		}

		if (allCovariants) {
			throw new BinaryInteractionFileException("Cannot set specific covariates for a file with all covariates tested");
		}

		if (!sortedIndices) {
			sortIndices();
		}

		if (!startedAddingCovariates) {
			startedAddingCovariates = true;
			covariatesTested = new int[countVariantGeneCombinations][];
			variantCummulativeGeneCounts = new int[variants.length];
			variantCummulativeGeneCounts[0] = variants[0].getGeneCount();
			for (int i = 1; i < variants.length; ++i) {
				variantCummulativeGeneCounts[i] = variantCummulativeGeneCounts[i - 1] + variants[i].getGeneCount();
			}
		}

		int variantIndex = variantMap.get(variantName);
		int geneIndex = genesMap.get(geneName);

		if (variantIndex == NO_ENTRY_INT_MAP) {
			throw new BinaryInteractionFileException("Unable to add interaction. Variant " + variantName + " not found");
		}

		if (geneIndex == NO_ENTRY_INT_MAP) {
			throw new BinaryInteractionFileException("Unable to add interaction. Gene " + geneName + " not found");
		}

		BinaryInteractionVariantCreator variant = variants[variantIndex];

		int variantGenePointerIndex = variant.getIndexOfGenePointer(geneIndex);

		if (variantGenePointerIndex < 0) {
			throw new BinaryInteractionFileException("Cannot add a interaction variant-gene combination does not exist: " + variantName + "-" + geneName);
		}

		int[] variantGeneCovariateArray = new int[covariateNames.length];

		for (int i = 0; i < covariateNames.length; ++i) {

			int covariteIndex = covariatesMap.get(covariateNames[i]);
			if (covariteIndex == NO_ENTRY_INT_MAP) {
				throw new BinaryInteractionFileException("Unable to add interaction combination. Covariate " + covariateNames[i] + " not found");
			}

			variantGeneCovariateArray[i] = covariteIndex;

			++interactions;

		}

		Arrays.sort(variantGeneCovariateArray);

		int indexInCovariatesTested = variantCummulativeGeneCounts[variantIndex - 1] + variantGenePointerIndex;

		if (indexInCovariatesTested >= covariatesTested.length) {
			throw new BinaryInteractionFileException("Something has gone wrong :(");
		}

		covariatesTested[indexInCovariatesTested] = variantGeneCovariateArray;

		++variantGenesWithCovariatesAdded;

	}

	public BinaryInteractionFile create() throws FileNotFoundException, IOException, BinaryInteractionFileException {

		if (created) {
			throw new BinaryInteractionFileException("You already created this file.");
		}

		created = true;

		if (!sortedIndices) {
			sortIndices();
		}

		if (allCovariants) {
			interactions = (long) countVariantGeneCombinations * (long) covariats.length;
		}

		if (!allCovariants && variantGenesWithCovariatesAdded != countVariantGeneCombinations) {
			throw new BinaryInteractionFileException("Trying to create binary file without setting covariates for all variant-gene combinations");
		}

		HashSet<Allele> alleles = createAlleleDictionary(variants);
		HashSet<String> chrs = createChrDictionary(variants, genes);

		DataOutputStream dataOutputStream = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(file)));

		dataOutputStream.writeByte(BinaryInteractionFile.MAGIC_1);
		dataOutputStream.writeByte(BinaryInteractionFile.MAGIC_2);

		dataOutputStream.writeByte(1);
		dataOutputStream.writeByte(0);

		long timeStamp = System.currentTimeMillis() / 1000L;
		dataOutputStream.writeLong(timeStamp);

		dataOutputStream.writeBoolean(allCovariants);
		dataOutputStream.writeBoolean(metaAnalysis);
		dataOutputStream.writeBoolean(normalQtlStored);
		dataOutputStream.writeBoolean(flippedZscoreStored);

		dataOutputStream.writeByte(0);
		dataOutputStream.writeByte(0);
		dataOutputStream.writeByte(0);
		dataOutputStream.writeByte(0);

		writeString(dataOutputStream, description);

		dataOutputStream.writeInt(cohorts.length);
		dataOutputStream.writeInt(chrs.size());
		dataOutputStream.writeInt(alleles.size());
		dataOutputStream.writeInt(genes.length);
		dataOutputStream.writeInt(variants.length);
		dataOutputStream.writeInt(covariats.length);
		dataOutputStream.writeLong(interactions);

		for (BinaryInteractionCohort cohort : cohorts) {
			writeString(dataOutputStream, cohort.getName());
			dataOutputStream.writeInt(cohort.getSampleCount());
		}

		TObjectIntHashMap<String> chrDictionary = new TObjectIntHashMap<String>(chrs.size(), 0.75f, NO_ENTRY_INT_MAP);
		for (String chr : chrs) {
			chrDictionary.put(chr, chrDictionary.size());
			writeString(dataOutputStream, chr);
		}

		TObjectIntHashMap<Allele> alleleDictionary = new TObjectIntHashMap<Allele>(alleles.size(), 0.75f, NO_ENTRY_INT_MAP);
		for (Allele alelle : alleles) {
			alleleDictionary.put(alelle, alleleDictionary.size());
			writeString(dataOutputStream, alelle.getAlleleAsString());
		}

		for (BinaryInteractionVariantCreator variant : variants) {

			writeString(dataOutputStream, variant.getName());
			dataOutputStream.writeInt(chrDictionary.get(variant.getChr()));
			dataOutputStream.writeInt(variant.getPos());
			dataOutputStream.writeInt(alleleDictionary.get(variant.getRefAllele()));
			dataOutputStream.writeInt(alleleDictionary.get(variant.getAltAllele()));
			dataOutputStream.writeInt(variant.getGeneCount());
			writeIntArray(dataOutputStream, variant.getGenePointers());

		}

		for (BinaryInteractionGeneCreator gene : genes) {

			writeString(dataOutputStream, gene.getName());
			dataOutputStream.writeInt(chrDictionary.get(gene.getChr()));
			dataOutputStream.writeInt(gene.getStart());
			dataOutputStream.writeInt(gene.getEnd());
			dataOutputStream.writeInt(gene.getVariantCount());
			writeIntArray(dataOutputStream, gene.getVariantPointers());

		}

		for (String covariate : covariats) {
			writeString(dataOutputStream, covariate);
		}

		if (!allCovariants) {
			for (int i = 0; i < covariatesTested.length; ++i) {
				dataOutputStream.writeInt(covariatesTested[i].length);
				writeIntArray(dataOutputStream, covariatesTested[i]);
			}
		}

		long startData = dataOutputStream.size();
		dataOutputStream.close();

		final long sizeNormalQtlSection;
		final long startNormalQtlSection;

		if (normalQtlStored) {
			startNormalQtlSection = startData;
			final long sizeQtlBlock = calculateSizeNormalQtlBlock(cohorts.length);
			sizeNormalQtlSection = sizeQtlBlock * countVariantGeneCombinations;
		} else {
			sizeNormalQtlSection = 0;
			startNormalQtlSection = -1;
		}


		final long startInteractionSection = startData + sizeNormalQtlSection;
		final long sizeInteractionBlock = calculateInteractionResultBlock(cohorts.length, flippedZscoreStored);

		RandomAccessFile fileRandomAccess = new RandomAccessFile(file, "rw"); //rw stands for open in read/write mode.
		fileRandomAccess.setLength(startData + sizeNormalQtlSection + sizeInteractionBlock);

		BinaryInteractionFileConstructorBuilder constructorBuilder = new BinaryInteractionFileConstructorBuilder();
		constructorBuilder.setAllCovariants(allCovariants);
		constructorBuilder.setCohorts(cohorts);
		constructorBuilder.setCovariatesTested(covariatesTested);
		constructorBuilder.setCovariats(covariats);
		constructorBuilder.setFileDescription(description);
		constructorBuilder.setFlippedZscoreStored(flippedZscoreStored);
		constructorBuilder.setGenes(genes);
		constructorBuilder.setInteractionFile(file);
		constructorBuilder.setMetaAnalysis(metaAnalysis);
		constructorBuilder.setNormalQtlStored(normalQtlStored);
		constructorBuilder.setReadOnly(false);
		constructorBuilder.setStartInteractionBlock(startInteractionSection);//TODO
		constructorBuilder.setStartQtlBlock(startNormalQtlSection);//TODO
		constructorBuilder.setTimeStamp(timeStamp);
		constructorBuilder.setVariants(variants);

		return constructorBuilder.createBinaryInteractionFile();

	}

	private void sortIndices() {
		for (BinaryInteractionGeneCreator gene : genes) {
			gene.sortGenePointers();
		}
		for (BinaryInteractionVariantCreator variant : variants) {
			variant.sortGenePointers();
		}
		sortedIndices = true;
	}

	public void setDescription(String description) {
		this.description = description;
	}

	private static void writeString(DataOutputStream dataOutputStream, String string) throws IOException {

		char[] chars = string.toCharArray();

		dataOutputStream.writeInt(chars.length);
		for (char c : chars) {
			dataOutputStream.writeChar(c);
		}

	}

	private static void writeIntArray(DataOutputStream dataOutputStream, int[] array) throws IOException {

		for (int e : array) {
			dataOutputStream.writeInt(e);
		}

	}

	private static HashSet<Allele> createAlleleDictionary(BinaryInteractionVariantCreator[] variants) {

		HashSet<Allele> alleles = new HashSet<Allele>();

		for (BinaryInteractionVariantCreator variant : variants) {
			alleles.add(variant.getRefAllele());
			alleles.add(variant.getAltAllele());
		}

		return alleles;

	}

	private static HashSet<String> createChrDictionary(BinaryInteractionVariantCreator[] variants, BinaryInteractionGene[] genes) {

		HashSet<String> chrs = new HashSet<String>();

		for (BinaryInteractionVariantCreator variant : variants) {
			chrs.add(variant.getChr());
		}

		for (BinaryInteractionGene gene : genes) {
			chrs.add(gene.getChr());
		}

		return chrs;

	}
}