package umcg.genetica.io.binInteraction;

import gnu.trove.map.hash.TObjectIntHashMap;
import java.util.Arrays;
import java.util.HashSet;
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

	public BinaryInteractionFileCreator(BinaryInteractionVariantCreator[] variants, BinaryInteractionGeneCreator[] genes, BinaryInteractionCohort[] cohorts, String[] covariats, boolean allCovariants, boolean metaAnalysis, boolean normalQtlStored, boolean flippedZscoreStored) throws BinaryInteractionFileException {
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
		BinaryInteractionGeneCreator gene = genes[geneIndex];

		int variantGenePointerIndex = variant.getIndexOfGenePointer(geneIndex);

		if (variantGenePointerIndex < 0) {
			throw new BinaryInteractionFileException("Cannot add a interaction variant-gene combination does not exist: " + variantName + "-" + geneName);
		}

		int[] variantGeneCovariateArray = new int[covariateNames.length];

		for (int i = 0 ; i < covariateNames.length ; ++i) {

			int covariteIndex = covariatesMap.get(covariateNames[i]);
			if (covariteIndex == NO_ENTRY_INT_MAP) {
				throw new BinaryInteractionFileException("Unable to add interaction combination. Covariate " + covariateNames[i] + " not found");
			}

			variantGeneCovariateArray[i] = covariteIndex;

		}
		
		Arrays.sort(variantGeneCovariateArray);

		int indexInCovariatesTested = variantCummulativeGeneCounts[variantIndex - 1] + variantGenePointerIndex;

		if (indexInCovariatesTested >= covariatesTested.length) {
			throw new BinaryInteractionFileException("Something has gone wrong :(");
		}
		
		covariatesTested[indexInCovariatesTested] = variantGeneCovariateArray;

		++variantGenesWithCovariatesAdded;

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
}
