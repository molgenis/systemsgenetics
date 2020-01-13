/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.depict2;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.Map;
import org.apache.log4j.Logger;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

/**
 *
 * @author patri
 */
public class ConvertGtexGct {
	
	private static final Logger LOGGER = Logger.getLogger(ConvertGtexGct.class);
	
	public static void convertGct(String gctFile, String outputMatrixPath) throws Exception{
		
		final HashSet<String> colToExclude = new HashSet<>();
		colToExclude.add("Description");
		
		DoubleMatrixDataset<String, String> gtexMedianExp = DoubleMatrixDataset.loadDoubleTextDoubleDataExlcudeCols(gctFile, '\t', colToExclude, 2);
		
		LOGGER.info("Loaded gtext data of " + gtexMedianExp.rows() + " genes in " + gtexMedianExp.columns() + " tissues");
		
		HashSet<String> includeGtextGenes = new HashSet<>(gtexMedianExp.getHashRows().keySet());
		HashMap<String, String> ensgToGtexMapping = new HashMap<>();
		
		for(String gtexGene : gtexMedianExp.getHashRows().keySet()){
			
			int indexOfPoint = gtexGene.indexOf('.');
			String gene;
			if(indexOfPoint > 0){
				gene = gtexGene.substring(0, indexOfPoint);
			} else {
				gene = gtexGene;
			}
			
			if(ensgToGtexMapping.containsKey(gene)){
				includeGtextGenes.remove(gtexGene);
				includeGtextGenes.remove(ensgToGtexMapping.get(gene));
				LOGGER.debug("Excluding: " + gene);
			}
			ensgToGtexMapping.put(gene, gtexGene);
			
		}
		
		if(includeGtextGenes.size() < gtexMedianExp.rows()){
			LOGGER.info("Excluding " + (gtexMedianExp.rows() - includeGtextGenes.size()) + " duplicated genes");
			gtexMedianExp = gtexMedianExp.viewRowSelection(includeGtextGenes);
		}
		
		
		
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

		gtexMedianExp.normalizeRows();
		
		ArrayList<String> rowNames = gtexMedianExp.getRowObjects();
		ArrayList<String> nonNanRowNames = new ArrayList<>(gtexMedianExp.rows());
		
		rows:
		for(int r = 0; r < gtexMedianExp.rows() ; ++r){
			for(int c = 0 ; c < gtexMedianExp.columns() ; ++c){
				if(Double.isNaN(gtexMedianExp.getElementQuick(r, c))){
					LOGGER.debug("Excluding NaN: " + rowNames.get(r));
					continue rows;
				}
			}
			nonNanRowNames.add(rowNames.get(r));
		}
		
		if(nonNanRowNames.size() < rowNames.size()){
			gtexMedianExp = gtexMedianExp.viewRowSelection(nonNanRowNames);
			LOGGER.info("Removing " + (rowNames.size() - nonNanRowNames.size()) + " rows with NaN after normalizing");
		}

		
		gtexMedianExp.saveBinary(outputMatrixPath);
		
	}
	
}
