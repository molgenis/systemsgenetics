/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.genenetworkbackend.div;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import com.opencsv.CSVParser;
import com.opencsv.CSVParserBuilder;
import com.opencsv.CSVReader;
import com.opencsv.CSVReaderBuilder;
import com.opencsv.CSVWriter;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

/**
 *
 * @author patri
 */
public class ExctractAnnotatedGenes {

	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) throws Exception {

		final File annotationMatrixFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Data31995Genes05-12-2017\\PCA_01_02_2018\\PathwayMatrix\\ALL_SOURCES_ALL_FREQUENCIES_phenotype_to_genes.txt_matrix.txt.gz");
		final File outputFolder = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Data31995Genes05-12-2017\\PCA_01_02_2018\\PathwayMatrix\\Extracts");
		final File ensgSymbolMappingFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\ensgHgnc.txt");

		Map<String, String> ensgSymbolMapping = loadEnsgToHgnc(ensgSymbolMappingFile);

		HashSet<String> termsToExtract = new HashSet<>();
		termsToExtract.add("HP:0001644");
		termsToExtract.add("HP:0001638");

		DoubleMatrixDataset<String, String> annotationMatrix = DoubleMatrixDataset.loadSubsetOfTextDoubleData(annotationMatrixFile.getAbsolutePath(), '\t', null, termsToExtract);

		ArrayList<String> rowNames = annotationMatrix.getRowObjects();

		for (String term : termsToExtract) {

			System.out.println(term);

			DoubleMatrix1D termAnnoations = annotationMatrix.getCol(term);
			ArrayList<String> termGenes = new ArrayList<>();

			for (int r = 0; r < termAnnoations.size(); ++r) {
				if (termAnnoations.getQuick(r) > 0) {
					termGenes.add(rowNames.get(r));
				}
			}

			CSVWriter writer = new CSVWriter(new FileWriter(new File(outputFolder, term.replace(':', '_') + ".txt")), '\t', '\0', '\0', "\n");
			String[] outputLine = new String[1];
			int c = 0;
			outputLine[c++] = "Gene";
			writer.writeNext(outputLine);

			for (String gene : termGenes) {
				c = 0;
				outputLine[c++] = gene;
				writer.writeNext(outputLine);
			}

			writer.close();

			CSVWriter writer2 = new CSVWriter(new FileWriter(new File(outputFolder, term.replace(':', '_') + "_symbol.txt")), '\t', '\0', '\0', "\n");
			outputLine = new String[1];
			c = 0;
			outputLine[c++] = "Gene";
			writer2.writeNext(outputLine);

			for (String gene : termGenes) {
				c = 0;
				outputLine[c++] = ensgSymbolMapping.get(gene);
				writer2.writeNext(outputLine);
			}

			writer2.close();

		}

	}

	private static Map<String, String> loadEnsgToHgnc(File mappingFile) throws IOException {

		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		final CSVReader reader = new CSVReaderBuilder(new BufferedReader(new FileReader(mappingFile))).withSkipLines(1).withCSVParser(parser).build();

		HashMap<String, String> mapping = new HashMap<>();

		String[] nextLine;
		while ((nextLine = reader.readNext()) != null) {

			mapping.put(nextLine[0], nextLine[1]);

		}

		return Collections.unmodifiableMap(mapping);

	}
}
