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
import java.util.LinkedHashSet;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

/**
 *
 * @author patri
 */
public class ParseLog {

	/**
	 * @param args the command line arguments
	 * @throws java.io.FileNotFoundException
	 */
	public static void main(String[] args) throws FileNotFoundException, IOException {

		final File logFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Data31995Genes05-12-2017\\PCA_01_02_2018\\predictions\\hpo_predictions.log");
		final int pcCount = 1588;
		final File outputMatrixFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Data31995Genes05-12-2017\\PCA_01_02_2018\\predictions\\hpo_compPathwayZscoreMatrix.txt.gz");

		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		CSVReader reader = new CSVReaderBuilder(new BufferedReader(new FileReader(logFile))).withSkipLines(0).withCSVParser(parser).build();

		LinkedHashSet<String> pathways = new LinkedHashSet<>();

		String[] nextLine;
		while ((nextLine = reader.readNext()) != null) {

			if (nextLine[0].equals("Processing geneset:")) {

				if (!pathways.add(nextLine[1])) {

					System.err.println("Duplicate pathway: " + nextLine[1]);
				}

			}

		}

		System.out.println("Number of pathways: " + pathways.size());

		LinkedHashSet<String> pcs = new LinkedHashSet<>();
		for (int pc = 1; pc <= pcCount; pc++) {
			pcs.add("PC" + pc);
		}

		DoubleMatrixDataset zscoreMatrix = new DoubleMatrixDataset(pathways, pcs);

		reader.close();

		reader = new CSVReaderBuilder(new BufferedReader(new FileReader(logFile))).withSkipLines(0).withCSVParser(parser).build();

		while ((nextLine = reader.readNext()) != null) {

			if (nextLine[0].equals("Processing geneset:")) {

				String pathway = nextLine[1];

				reader.readNext();//skip header

				for (int pc = 1; pc <= pcCount; pc++) {

					nextLine = reader.readNext();
					
					if(!nextLine[0].equals("PC" + pc)){
						System.err.println("Error parsing z-scores for: " + pathway);
					}
					
					zscoreMatrix.setElement(pathway, nextLine[0], Double.valueOf(nextLine[4]));
					
				}

			}

		}

		reader.close();
		zscoreMatrix.save(outputMatrixFile);

	}

}
