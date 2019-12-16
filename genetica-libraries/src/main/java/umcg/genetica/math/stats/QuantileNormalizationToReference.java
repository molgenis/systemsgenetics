/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.math.stats;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.algo.DoubleSorting;

/**
 *
 * @author patri
 */
public class QuantileNormalizationToReference {

	/**
	 * The columns of the matrix will q-normed to the reference.
	 * 
	 * @param reference
	 * @param matrix 
	 */
	public static void inplaceQuantileNormalizationToReference(DoubleMatrix1D reference, DoubleMatrix2D matrix){
		
		if(reference.size() != matrix.rows()){
			throw new RuntimeException("Reference not equalt to rows in matrix");
		}
		
		DoubleMatrix1D referenceSorted = DoubleSorting.quickSort.sort(reference);
		
		for(int c = 0 ; c < matrix.columns() ; ++c ){
			
			DoubleMatrix1D colSorted = DoubleSorting.quickSort.sort(matrix.viewColumn(c));
			colSorted.assign(referenceSorted);
			
		}
		
	}
	
}
