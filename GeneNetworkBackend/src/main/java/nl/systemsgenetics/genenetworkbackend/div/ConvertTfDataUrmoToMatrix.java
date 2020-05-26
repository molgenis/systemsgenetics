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
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import static nl.systemsgenetics.genenetworkbackend.ConvertGoToMatrix.readGenes;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

/**
 *
 * @author patri
 */
public class ConvertTfDataUrmoToMatrix {

	/**
	 * @param args the command line arguments
	 * @throws java.lang.Exception
	 */
	public static void main(String[] args) throws Exception {

		final File tfFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\TF\\TF_miR_databases_combined_ENSG_added.txt");
		final File geneOrderFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Data31995Genes05-12-2017\\PCA_01_02_2018\\genes.txt");
		final File outputFolder = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Data31995Genes05-12-2017\\PCA_01_02_2018\\PathwayMatrixTf\\");

		HashMap<String, HashMap<String, HashSet<String>>> tfdatabasesPathwayToGenes = loadTfData(tfFile);

		ArrayList<String> geneOrder = readGenes(geneOrderFile);

		System.out.println("Genes in order file: " + geneOrder.size());

		for (String tfDatabase : tfdatabasesPathwayToGenes.keySet()) {

			final File outputFile = new File(outputFolder, tfDatabase + "_matrix.txt.gz");
			final File outputFile2 = new File(outputFolder, tfDatabase + "_genesInPathways.txt");

			HashMap<String, HashSet<String>> pathwayToGenes = tfdatabasesPathwayToGenes.get(tfDatabase);

			System.out.println("Total genesets of " + tfDatabase + ": " + pathwayToGenes.size());

			DoubleMatrixDataset<String, String> pathwayMatrix = new DoubleMatrixDataset(geneOrder, pathwayToGenes.keySet());

			HashSet<String> genesWithPathway = new HashSet<>(10000);
			BufferedWriter geneWriter = new BufferedWriter(new FileWriter(outputFile2));

			for (Map.Entry<String, HashSet<String>> pathwayToGenesEntry : pathwayToGenes.entrySet()) {

				String pathway = pathwayToGenesEntry.getKey();

				for (String gene : pathwayToGenesEntry.getValue()) {
					
					if (pathwayMatrix.containsRow(gene)) {

						if (genesWithPathway.add(gene)) {
							//add to genes file if not already done
							geneWriter.write(gene);
							geneWriter.write('\n');
						}

						pathwayMatrix.setElement(gene, pathway, 1);
					}

				}

			}

			pathwayMatrix.save(outputFile);
			geneWriter.close();

			System.out.println("Genes in " + tfDatabase + " pathway: " + genesWithPathway.size());

		}

	}

	private static HashMap<String, HashMap<String, HashSet<String>>> loadTfData(File tfFile) throws FileNotFoundException, IOException {

		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		final CSVReader reader = new CSVReaderBuilder(new BufferedReader(new FileReader(tfFile))).withSkipLines(1).withCSVParser(parser).build();

		HashMap<String, HashMap<String, HashSet<String>>> tfdatabasesPathwayToGenes = new HashMap<>();

		String[] nextLine;
		while ((nextLine = reader.readNext()) != null) {

			if (nextLine[0].charAt(0) == '!') {
				continue;
			}

			String database = nextLine[0];
			String pathway = nextLine[1];
			String ensgId = nextLine[3];

			HashMap<String, HashSet<String>> pathwayToGenes = tfdatabasesPathwayToGenes.get(database);
			if (pathwayToGenes == null) {
				pathwayToGenes = new HashMap<>();
				tfdatabasesPathwayToGenes.put(database, pathwayToGenes);
			}

			HashSet<String> pathwayGenes = pathwayToGenes.get(pathway);
			if (pathwayGenes == null) {
				pathwayGenes = new HashSet<>();
				pathwayToGenes.put(pathway, pathwayGenes);
			}

			pathwayGenes.add(ensgId);

		}

		return tfdatabasesPathwayToGenes;

	}

}
