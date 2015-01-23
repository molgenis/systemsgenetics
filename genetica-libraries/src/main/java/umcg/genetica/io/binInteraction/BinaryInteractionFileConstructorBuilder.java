package umcg.genetica.io.binInteraction;

import umcg.genetica.io.binInteraction.gene.BinaryInteractionGene;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import umcg.genetica.io.binInteraction.variant.BinaryInteractionVariant;


public class BinaryInteractionFileConstructorBuilder {
	private File interactionFile;
	private boolean readOnly;
	private BinaryInteractionCohort[] cohorts;
	private BinaryInteractionGene[] genes;
	private BinaryInteractionVariant[] variants;
	private int[][] covariatesTested = null;
	private String[] covariates;
	private long timeStamp;
	private boolean allCovariants;
	private boolean metaAnalysis;
	private boolean normalQtlStored;
	private boolean flippedZscoreStored;
	private String fileDescription;
	private long interactions;
	private long startQtlBlock;
	private long startInteractionBlock;

	public BinaryInteractionFileConstructorBuilder() {
	}

	public BinaryInteractionFileConstructorBuilder setInteractionFile(File interactionFile) {
		this.interactionFile = interactionFile;
		return this;
	}

	public BinaryInteractionFileConstructorBuilder setReadOnly(boolean readOnly) {
		this.readOnly = readOnly;
		return this;
	}

	public BinaryInteractionFileConstructorBuilder setCohorts(BinaryInteractionCohort[] cohorts) {
		this.cohorts = cohorts;
		return this;
	}

	public BinaryInteractionFileConstructorBuilder setGenes(BinaryInteractionGene[] genes) {
		this.genes = genes;
		return this;
	}

	public BinaryInteractionFileConstructorBuilder setVariants(BinaryInteractionVariant[] variants) {
		this.variants = variants;
		return this;
	}

	public BinaryInteractionFileConstructorBuilder setCovariatesTested(int[][] covariatesTested) {
		this.covariatesTested = covariatesTested;
		return this;
	}

	public BinaryInteractionFileConstructorBuilder setCovariates(String[] covariates) {
		this.covariates = covariates;
		return this;
	}
	
	public BinaryInteractionFileConstructorBuilder setTimeStamp(long timeStamp) {
		this.timeStamp = timeStamp;
		return this;
	}

	public BinaryInteractionFileConstructorBuilder setAllCovariants(boolean allCovariants) {
		this.allCovariants = allCovariants;
		return this;
	}

	public BinaryInteractionFileConstructorBuilder setMetaAnalysis(boolean metaAnalysis) {
		this.metaAnalysis = metaAnalysis;
		return this;
	}

	public BinaryInteractionFileConstructorBuilder setNormalQtlStored(boolean normalQtlStored) {
		this.normalQtlStored = normalQtlStored;
		return this;
	}
	
	public BinaryInteractionFileConstructorBuilder setFlippedZscoreStored(boolean flippedZscoreStored) {
		this.flippedZscoreStored = flippedZscoreStored;
		return this;
	}

	public BinaryInteractionFileConstructorBuilder setFileDescription(String fileDescription) {
		this.fileDescription = fileDescription;
		return this;
	}

	public BinaryInteractionFileConstructorBuilder setInteractions(long interactions) {
		this.interactions = interactions;
		return this;
	}

	public BinaryInteractionFileConstructorBuilder setStartQtlBlock(long startQtlBlock) {
		this.startQtlBlock = startQtlBlock;
		return this;
	}

	public BinaryInteractionFileConstructorBuilder setStartInteractionBlock(long startInteractionBlock) {
		this.startInteractionBlock = startInteractionBlock;
		return this;
	}

	public BinaryInteractionFile createBinaryInteractionFile() throws BinaryInteractionFileException, FileNotFoundException, IOException {
		return new BinaryInteractionFile(interactionFile, readOnly, cohorts, genes, variants, covariates, covariatesTested, timeStamp, allCovariants, metaAnalysis, normalQtlStored, flippedZscoreStored, fileDescription, interactions, startQtlBlock, startInteractionBlock);
	}

}
