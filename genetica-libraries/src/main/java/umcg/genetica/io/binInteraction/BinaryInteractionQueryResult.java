package umcg.genetica.io.binInteraction;

/**
 *
 * @author Patrick Deelen
 */
public class BinaryInteractionQueryResult {

	private final String variantName;
	private final String geneName;
	private final String covariateName;
	private final BinaryInteractionQtlZscores qtlZscores;
	private final BinaryInteractionZscores interactionZscores;

	public BinaryInteractionQueryResult(String variantName, String geneName, String covariateName, BinaryInteractionQtlZscores qtlZscores, BinaryInteractionZscores interactionZscores) {
		this.variantName = variantName;
		this.geneName = geneName;
		this.covariateName = covariateName;
		this.qtlZscores = qtlZscores;
		this.interactionZscores = interactionZscores;
	}

	public String getVariantName() {
		return variantName;
	}

	public String getGeneName() {
		return geneName;
	}

	public String getCovariateName() {
		return covariateName;
	}

	public BinaryInteractionQtlZscores getQtlZscores() {
		return qtlZscores;
	}

	public BinaryInteractionZscores getInteractionZscores() {
		return interactionZscores;
	}
	
}
