package umcg.genetica.io.binInteraction;

import java.io.File;


public class BinaryInteractionFileConstructorBuilder {
	private File interactionFile;
	private boolean readOnly;
	private BinaryInteractionCohort[] cohorts;
	private BinaryInteractionGene[] genes;
	private BinaryInteractionVariant[] variants;
	private String[] covariats;
	private String[] chrDictionary;
	private String[] alleleDictionary;
	private long timeStamp;
	private boolean allCovariants;
	private boolean metaAnalysis;
	private boolean normalQtlStored;
	private String fileDescription;
	private int interactions;
	private int startQtlBlock;
	private int startInteractionBlock;

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

	public BinaryInteractionFileConstructorBuilder setCovariats(String[] covariats) {
		this.covariats = covariats;
		return this;
	}

	public BinaryInteractionFileConstructorBuilder setChrDictionary(String[] chrDictionary) {
		this.chrDictionary = chrDictionary;
		return this;
	}

	public BinaryInteractionFileConstructorBuilder setAlleleDictionary(String[] alleleDictionary) {
		this.alleleDictionary = alleleDictionary;
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

	public BinaryInteractionFileConstructorBuilder setFileDescription(String fileDescription) {
		this.fileDescription = fileDescription;
		return this;
	}

	public BinaryInteractionFileConstructorBuilder setInteractions(int interactions) {
		this.interactions = interactions;
		return this;
	}

	public BinaryInteractionFileConstructorBuilder setStartQtlBlock(int startQtlBlock) {
		this.startQtlBlock = startQtlBlock;
		return this;
	}

	public BinaryInteractionFileConstructorBuilder setStartInteractionBlock(int startInteractionBlock) {
		this.startInteractionBlock = startInteractionBlock;
		return this;
	}

	public BinaryInteractionFile createBinaryInteractionFile() {
		return new BinaryInteractionFile(interactionFile, readOnly, cohorts, genes, variants, covariats, chrDictionary, alleleDictionary, timeStamp, allCovariants, metaAnalysis, normalQtlStored, fileDescription, interactions, startQtlBlock, startInteractionBlock);
	}

}
