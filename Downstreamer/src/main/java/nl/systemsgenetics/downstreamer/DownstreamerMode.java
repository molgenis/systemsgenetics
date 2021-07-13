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
	CONVERT_EQTL,
	CONVERT_GTEX,
	CONVERT_TXT_MERGE,
	GET_NORMALIZED_GENEP,
	GET_PATHWAY_LOADINGS,
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
	SPECIAL;
}
