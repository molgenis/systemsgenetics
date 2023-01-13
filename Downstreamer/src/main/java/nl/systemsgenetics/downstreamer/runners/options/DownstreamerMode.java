/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.downstreamer.runners.options;

/**
 *
 * @author patri
 */
public enum DownstreamerMode {

	// Runners for the updated workflow
	REGRESS("Run the Downstreamer regression model without any pre-proccessing using preformatted data"),

	// Converters, use OptionsModeConverters
	CONVERT_TXT("Convert a txt z-score matrix to binary. Use --gwas, --output"),
	CONVERT_BIN("Convert a binary matrix to a txt. Use --gwas and --output optionally --columnsToExtract"),
	CONVERT_MERGE_BIN("Merge multiple .dat matrix files on overlapping rows"),
	CONVERT_DAT_TO_DATG("Convert a .dat binary matrix to a row compressed .datg"),
	CONVERT_EQTL("Convert binary matrix with eQTL z-scores from our pipeline. Use --gwas and --output"),
	CONVERT_GTEX("Convert GTEx median tissue GCT file. Use --gwas for the GCT file and --output"),
	CONVERT_EXP("Convert a tab seperated expression matrix and normalize genes. Use --gwas (for exp data) and --output optionally --columnsToExtract"),
	CONVERT_PTOZSCORE("Convert a .txt matrix of p-values to Z scores (two-tailed)"),

	MATRIX_TRANSPOSE("Transposes a binary matrix. Use --gwas and --output"),
	MATRIX_SUBSET("Subset either a binary (.dat, .datg) or .txt matrix on specified rows or columns. Is saved in same format as input."),

	PREPARE_GWAS("Convert a txt z-score matrix to binary. Use --gwas, --output and optionally --pvalueToZscore if the matrix contains p-values instead of z-scores."),
	PREPARE_GWAS_MERGE("Merge multiple txt pvalue files into one matrix containing only overlapping snps"),

	// Utilities relating to calculating co-regulation, use OptionsModeCoreg
	COREG_RTOZSCORE("Convert correlation matrix with r values to Z-score matrix. Must be used together with -ns. --gwas Input can be .txt or binary."),
	COREG_ZSCORETOR("Convert Z-score matrix to correlation matrix with r values. Must be used together with -ns. --gwas Input can be .txt or binary."),
	COREG_CORRELATE_GENES("Create gene correlation matrix with 0 on diagonal. Use --gwas as input matrix (genes on row, tab separated), --output and --genes. Optionally use --corZscore to create Z-score matrix"),
	COREG_REMOVE_CIS_COEXP("Set cis gene-gene correlation in the co-expression / co-regulation matrix to zero [not recommended]."),
	COREG_PCA("Run a PCA on a matrix. Matrix must be binary."),
	COREG_INVESTIGATE_NETWORK("Calculate degree statistics on a co-regulation matrix. Expects Z-scores, so run -m COREG_RTOZSCORE first."),

	// Testing modes
	TEST_DECOMP("Test entrypoint for eigen decomposition. Remove upon release."),

	// Runners that expect the old main workflow, use DownstreamerOptionsDeprecated
	STEP1("[deprecated] Run the Downstreamer prioritization"),
	STEP2("[deprecated] Run the Downstreamer prioritization starting at stage 2"),
	CREATE_EXCEL("[deprecated] Convert old implementation DS output to an Excel sheet"),
	PRIO_GENE_ENRICH("[deprecated] Perform over-representation enrichment of Downstreamer prioritized genes."),
	EXPAND_PATHWAYS("[to update] Expand known pathway annotations based on predicted memberships in a gene network."),

	// To update in future if we need them
	GET_NORMALIZED_GENEP("[to update] [currently not avail]"),
	GET_PATHWAY_LOADINGS("[to update] [currently not avail]"),
	GET_MARKER_GENES("[to update] [currently not avail]"),
	TOP_HITS("[to update] [currently not avail]"),
	PREPARE_GENE_PVALUES("[to update] [currently not avail]"),
	SPECIAL("[deprecated][remove?]"),
	FIRST1000("[deprecated][remove?]");

	private final String description;

	DownstreamerMode(String description) {
		this.description = description;
	}

	public String getDescription() {
		return  description;
	}

	public static String getFullDescriptionString() {
		StringBuilder output = new StringBuilder();
		for (DownstreamerMode mode: DownstreamerMode.values()) {
			output.append(mode.toString());
			output.append(" - ");
			output.append(mode.getDescription());

			output.append("\n");
			for (int i = 0; i < 40; i++) {
				output.append("-");
			}
			output.append("\n");
		}

		return output.toString();
	}
}
