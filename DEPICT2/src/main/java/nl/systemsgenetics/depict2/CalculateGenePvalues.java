/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.depict2;

import java.util.ArrayList;
import java.util.HashMap;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

/**
 *
 * @author patri
 */
public class CalculateGenePvalues {

	/**
	 * 
	 * @param variantPhenotypeMatrix rows: variants, cols: phenotypes
	 * @param genotypeCovarianceSource Source data used to calculate correlations between SNPs
	 * @param genes genes for which to calculate p-values. 
	 * @param snpMapping genomic positions variants 
	 * @return gene p-value matrix for each phenotype. rows: genes in same order as genes parameter, cols: phenotypes
	 */
	public static DoubleMatrixDataset<String, String> calculatorGenePvalues(DoubleMatrixDataset<String, String> variantPhenotypeMatrix, GenotypeCovarianceSource genotypeCovarianceSource, ArrayList<Genes> genes, HashMap<String, ChrPos> snpMapping){
		
		return null;
		
	}
	
}
