/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.depict2;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.Map;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

/**
 *
 * @author patri
 */
public class ConvertGtexGct {
	
	public static void convertGct(String gctFile, String outputMatrixPath) throws Exception{
		
		final HashSet<String> colToExclude = new HashSet<>();
		colToExclude.add("Description");
		
		DoubleMatrixDataset<String, String> gtexMedianExp = DoubleMatrixDataset.loadDoubleTextDoubleDataExlcudeCols(gctFile, '\t', colToExclude, 2);
		
		LinkedHashMap<String, Integer> newHashRow = new LinkedHashMap<>(gtexMedianExp.rows());
		
		for(Map.Entry<String, Integer> rowEntry : gtexMedianExp.getHashRows().entrySet()){
			
			String gene = rowEntry.getKey();
			
			int indexOfPoint = gene.indexOf('.');
			if(indexOfPoint > 0){
				gene = gene.substring(0, indexOfPoint);
			}
			
			newHashRow.put(gene, rowEntry.getValue());
			
		}
		
		gtexMedianExp.setHashRows(newHashRow);
		
		final DoubleMatrix2D gtexMatrix = gtexMedianExp.getMatrix();
		final int numberOfTissues = gtexMatrix.columns();
		
		final double numberOfTissuesDouble = gtexMatrix.columns();
		
		for(int g = 0; g < gtexMatrix.rows() ; ++g){
			
			DoubleMatrix1D geneRow = gtexMatrix.viewRow(g);
			
			double geneMean = geneRow.zSum() / numberOfTissuesDouble;
			
			for(int t = 0 ; t < numberOfTissues ; ++t){
				geneRow.setQuick(t, geneRow.getQuick(t) - geneMean);
			}
			
		}
		
		gtexMedianExp.saveBinary(outputMatrixPath);
		
	}
	
}
