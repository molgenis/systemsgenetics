
/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.genenetworkbackend.div;

import com.opencsv.CSVParser;
import com.opencsv.CSVParserBuilder;
import com.opencsv.CSVReader;
import com.opencsv.CSVReaderBuilder;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import org.apache.commons.lang.StringUtils;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

/**
 *
 * @author patri
 */
public class CreateCgdMatrix {

	/**
	 * @param args the command line arguments
	 * @throws java.io.FileNotFoundException
	 */
	public static void main(String[] args) throws FileNotFoundException, IOException {
		
		final File cgdFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\CGD.txt");
		final File cgdMatrixFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\CGD_matrix.txt");
		
		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		final CSVReader reader = new CSVReaderBuilder(new BufferedReader(new FileReader(cgdFile))).withSkipLines(1).withCSVParser(parser).build();
		
		HashSet<String> cgdInfo = new HashSet<>(30);
		cgdInfo.add("AD");
		cgdInfo.add("AD");
		cgdInfo.add("XL");
		cgdInfo.add("Pediatric");
		
		HashMap<String, HashSet<String>> geneInfoMap = new HashMap<>();
		
		String[] nextLine;
		while ((nextLine = reader.readNext()) != null) {
			
			HashSet<String> geneInfo = new HashSet<>();
			geneInfoMap.put(nextLine[0], geneInfo);
			
			String diseaseManifistationField = nextLine[7];
			String[] diseaseManifistations = StringUtils.split(diseaseManifistationField, ';');
			
			geneInfo.add(nextLine[4].trim());
			geneInfo.add(nextLine[5].trim());
			
			for( int i = 0 ; i < diseaseManifistations.length ; ++i){
				
				diseaseManifistations[i] = diseaseManifistations[i].trim();
				cgdInfo.add(diseaseManifistations[i]);
				geneInfo.add(diseaseManifistations[i]);
				
			}
			
			
		}
		
		cgdInfo.remove("General");
		
		DoubleMatrixDataset<String, String> cgdInfoTable = new DoubleMatrixDataset<>(geneInfoMap.keySet(), cgdInfo);
		
		for(Map.Entry<String, HashSet<String>> geneInfoEntry : geneInfoMap.entrySet()){
			
			String gene = geneInfoEntry.getKey();
			HashSet<String> geneInfo = geneInfoEntry.getValue();
			
			for(String geneProperty : geneInfo){
				
				if(cgdInfoTable.containsCol(geneProperty)){
					cgdInfoTable.setElement(gene, geneProperty, 1);
				}
				
			}
			
			
		}
		
		cgdInfoTable.save(cgdMatrixFile);
		
		
	}
	
}
