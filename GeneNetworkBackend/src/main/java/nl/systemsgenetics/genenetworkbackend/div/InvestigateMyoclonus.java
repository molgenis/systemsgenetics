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
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

/**
 *
 * @author patri
 */
public class InvestigateMyoclonus {

	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) throws Exception {

		final File ensgSymbolMappingFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\ensgHgnc.txt");
		final Map<String, String> symbolEnsgMapping = loadHgncToEnsg(ensgSymbolMappingFile);

		final DoubleMatrixDataset<String, String> predictions = DoubleMatrixDataset.loadDoubleTextData("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Myoclonus\\2019_12_17_myoclonusClusters_predictions.txt", '\t');

		final ArrayList<String> genes = predictions.getRowObjects();
		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		

		ArrayList<String> dnaNumbers = new ArrayList<>();
		ArrayList<File> vcfFiles = new ArrayList<>();

		final Pattern pattern = Pattern.compile("DNA[^_]+");

		for (File vcfFile : new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Myoclonus\\gavinRes").listFiles()) {

			Matcher matcher = pattern.matcher(vcfFile.getName());
			matcher.find();
			String dnaNumber = matcher.group(0);

			dnaNumbers.add(dnaNumber);
			vcfFiles.add(vcfFile);

		}
		
		DoubleMatrixDataset<String, String> results = DoubleMatrixDataset.concatColumns(predictions, new DoubleMatrixDataset<>(genes, dnaNumbers));

		String[] nextLine;
		for (int i = 0; i < dnaNumbers.size(); ++i) {
			
			String dnaNumber = dnaNumbers.get(i);
			
			final CSVReader gavinFileReader = new CSVReaderBuilder(new BufferedReader(new FileReader(vcfFiles.get(i)))).withSkipLines(1).withCSVParser(parser).build();
			while ((nextLine = gavinFileReader.readNext()) != null) {
				
				String ensgGene = symbolEnsgMapping.get(nextLine[8]);
				
				if(ensgGene != null && results.containsRow(ensgGene)){
					System.out.println(dnaNumber + "\t" + ensgGene);
					results.setElement(ensgGene, dnaNumber, 1);
				}
				
			}

		}
		
		results.save("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Myoclonus\\results.txt");

	}

	private static Map<String, String> loadHgncToEnsg(File mappingFile) throws IOException {

		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		final CSVReader reader = new CSVReaderBuilder(new BufferedReader(new FileReader(mappingFile))).withSkipLines(1).withCSVParser(parser).build();

		HashMap<String, String> mapping = new HashMap<>();

		String[] nextLine;
		while ((nextLine = reader.readNext()) != null) {

			mapping.put(nextLine[1], nextLine[0]);

		}

		return Collections.unmodifiableMap(mapping);

	}

}
