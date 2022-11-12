/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.downstreamer;

/**
 *
 * @author patri
 */
public enum DownstreamerMode {
	STEP1,
	STEP2,
	CONVERT_TXT,
	CONVERT_BIN,
	CONVERT_DAT_TO_DATG,
	CONVERT_EQTL,
	CONVERT_GTEX,
	CONVERT_TXT_MERGE,
	GET_NORMALIZED_GENEP,
	GET_PATHWAY_LOADINGS,
	GET_MARKER_GENES,
	GET_NETWORK_DEGREE,
	PTOZSCORE,
	PCA,
	FIRST1000,
	CORRELATE_GENES,
	TRANSPOSE,
	CONVERT_EXP,
	MERGE_BIN,
	PRIO_GENE_ENRICH,
	INVESTIGATE_NETWORK,
	R_2_Z_SCORE,
	CREATE_EXCEL,
	TOP_HITS,
	REMOVE_CIS_COEXP,
	SPECIAL,
	SUBSET_MATRIX,
	EXPAND_PATHWAYS,
	PREPARE_GENE_PVALUES;
}
