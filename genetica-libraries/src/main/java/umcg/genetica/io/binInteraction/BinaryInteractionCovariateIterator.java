package umcg.genetica.io.binInteraction;

import java.io.IOException;
import java.util.Iterator;
import java.util.NoSuchElementException;

/**
 *
 * @author Patrick Deelen
 */
public class BinaryInteractionCovariateIterator implements Iterator<BinaryInteractionQueryResult>{

	private final String variantName;
	private final String geneName;
	private final int[] covariatsTested;
	private final BinaryInteractionFile file;
	private final long startVariantGeneBlock;
	private final BinaryInteractionQtlZscores qtlZscores;
	private int variantGeneCovariateIndex;

	public BinaryInteractionCovariateIterator(String variantName, String geneName, int[] covariatsTested, BinaryInteractionFile file, long startVariantGeneBlock, BinaryInteractionQtlZscores qtlZscores) {
		this.variantName = variantName;
		this.geneName = geneName;
		this.covariatsTested = covariatsTested;
		this.file = file;
		this.startVariantGeneBlock = startVariantGeneBlock;
		this.qtlZscores = qtlZscores;
		variantGeneCovariateIndex = -1;
	}
	
	@Override
	public boolean hasNext() {
		if(covariatsTested == null){
			return variantGeneCovariateIndex + 1 < file.getCovariateCount();
		} else {
			return variantGeneCovariateIndex + 1 < covariatsTested.length;
		}
	}

	@Override
	public BinaryInteractionQueryResult next() {
		if(!hasNext()){
			throw new NoSuchElementException();
		}
		
		variantGeneCovariateIndex++;
		BinaryInteractionZscores interactionRestuls;
		try {
			interactionRestuls = file.readInteractionResults( startVariantGeneBlock + (variantGeneCovariateIndex * file.sizeInteractionBlock));
		} catch (BinaryInteractionFileException | IOException ex) {
			throw new RuntimeException(ex);
		}
		String covariateName = covariatsTested == null ? file.covariates[variantGeneCovariateIndex] : file.covariates[covariatsTested[variantGeneCovariateIndex]];
		return new BinaryInteractionQueryResult(variantName, geneName, covariateName, qtlZscores, interactionRestuls);
	}

	@Override
	public void remove() {
		throw new UnsupportedOperationException("Not supported ever.");
	}

}
