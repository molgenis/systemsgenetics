package nl.systemsgenetics.genenetworkbackend;

import com.opencsv.CSVParser;
import com.opencsv.CSVParserBuilder;
import com.opencsv.CSVReader;
import com.opencsv.CSVReaderBuilder;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author patri
 */
public class ConvertReactomeToMatrix {

	/**
	 * @param args the command line arguments
	 * @throws java.io.IOException
	 * @throws java.lang.Exception
	 */
	public static void main(String[] args) throws IOException, Exception {

		final File pathwayFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Reactome\\Ensembl2Reactome_All_Levels.txt");
		final File geneOrderFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Data31995Genes05-12-2017\\PCA_01_02_2018\\genes.txt");
		final File outputFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Data31995Genes05-12-2017\\PCA_01_02_2018\\PathwayMatrix\\" + pathwayFile.getName() + "_matrix.txt");
		final File outputFile2 = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Data31995Genes05-12-2017\\PCA_01_02_2018\\PathwayMatrix\\" + pathwayFile.getName() + "_genesInPathways.txt");
		
		HashMap<String, HashSet<String>> pathwayToGenes = readPathwayFile(pathwayFile);

		ArrayList<String> geneOrder = readGenes(geneOrderFile);

		System.out.println("Total genesets: " + pathwayToGenes.size());
		System.out.println("Genes in order file: " + geneOrder.size());

		DoubleMatrixDataset<String, String> pathwayMatrix = new DoubleMatrixDataset(geneOrder, pathwayToGenes.keySet());
		
		HashSet<String> genesWithPathway = new HashSet<>(10000);
		BufferedWriter geneWriter = new BufferedWriter(new FileWriter(outputFile2));

		for (Map.Entry<String, HashSet<String>> pathwayToGenesEntry : pathwayToGenes.entrySet()) {

			String pathway = pathwayToGenesEntry.getKey();

			for (String gene : pathwayToGenesEntry.getValue()) {

				if (pathwayMatrix.containsRow(gene)) {
					
					if(genesWithPathway.add(gene)){
						//add to genes file if not already done
						geneWriter.write(gene);
						geneWriter.write('\n');
					}
					
					pathwayMatrix.setElement(gene, pathway, 1);
				}

			}

		}

		geneWriter.close();
		pathwayMatrix.save(outputFile);

		System.out.println("Genes in pathway: " + genesWithPathway.size());
		
	}

	private static HashMap<String, HashSet<String>> readPathwayFile(File pathwayFile) throws Exception {

		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		final CSVReader reader = new CSVReaderBuilder(new BufferedReader(new FileReader(pathwayFile))).withSkipLines(0).withCSVParser(parser).build();

		HashMap<String, HashSet<String>> pathwayToGenes = new HashMap<>();

		String[] nextLine;
		while ((nextLine = reader.readNext()) != null) {

			if (nextLine[5].equals("Homo sapiens")) {

				String pathway = nextLine[1];
				String ensgId = nextLine[0];

				HashSet<String> pathwayGenes = pathwayToGenes.get(pathway);
				if (pathwayGenes == null) {
					pathwayGenes = new HashSet<>();
					pathwayToGenes.put(pathway, pathwayGenes);
				}

				pathwayGenes.add(ensgId);

			}
		}

		return pathwayToGenes;

	}

	private static ArrayList<String> readGenes(File geneOrderFile) throws IOException {

		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		final CSVReader reader = new CSVReaderBuilder(new BufferedReader(new FileReader(geneOrderFile))).withSkipLines(0).withCSVParser(parser).build();

		String[] nextLine;
		ArrayList<String> geneOrder = new ArrayList<>();

		while ((nextLine = reader.readNext()) != null) {

			geneOrder.add(nextLine[0]);

		}

		return geneOrder;

	}

}
