package nl.systemsgenetics.genenetworkbackend;

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
public class ConvertMyoclonusClustersToMatrix {

	/**
	 * @param args the command line arguments
	 * @throws java.io.IOException
	 * @throws java.lang.Exception
	 */
	public static void main(String[] args) throws IOException, Exception {

		//final File hpoFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\HPO\\135\\ALL_SOURCES_ALL_FREQUENCIES_diseases_to_genes_to_phenotypes.txt");
		
		final File geneOrderFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Data31995Genes05-12-2017\\PCA_01_02_2018\\genes.txt");
		
		final File clusterFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\myoclonusClusters.txt");
		final File outputFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Data31995Genes05-12-2017\\PCA_01_02_2018\\PathwayMatrix\\" + clusterFile.getName() + "_matrix.txt.gz");
		
		
		HashMap<String, HashSet<String>> clusterToGenes = readClusterFile(clusterFile);

		ArrayList<String> geneOrder = readGenes(geneOrderFile);

		System.out.println("Total HPO terms: " + clusterToGenes.size());
		System.out.println("Genes in order file: " + geneOrder.size());

		DoubleMatrixDataset<String, String> hpoMatrix = new DoubleMatrixDataset(geneOrder, clusterToGenes.keySet());

		for (Map.Entry<String, HashSet<String>> clusterToGenesEntry : clusterToGenes.entrySet()) {

			String cluster = clusterToGenesEntry.getKey();

			for (String gene : clusterToGenesEntry.getValue()) {

				if (hpoMatrix.containsRow(gene)) {
					System.out.println(gene + " - " + cluster);
					hpoMatrix.setElement(gene, cluster, 1);
				}

			}

		}

		hpoMatrix.save(outputFile);


	}

	
	private static HashMap<String, HashSet<String>> readClusterFile(File hpoFile) throws Exception {

		final CSVParser hpoParser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		final CSVReader hpoReader = new CSVReaderBuilder(new BufferedReader(new FileReader(hpoFile))).withSkipLines(1).withCSVParser(hpoParser).build();

		HashMap<String, HashSet<String>> hpoToGenes = new HashMap<>();

		String[] nextLine;
		while ((nextLine = hpoReader.readNext()) != null) {
			String gene = nextLine[0];
			String cluster = nextLine[1];

				HashSet<String> hpoGenes = hpoToGenes.get(cluster);
				if (hpoGenes == null) {
					hpoGenes = new HashSet<>();
					hpoToGenes.put(cluster, hpoGenes);
				}

				
					hpoGenes.add(gene);
				


		}

		return hpoToGenes;

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

	public static HashMap<String, ArrayList<String>> loadHgncToEnsgMap(File map) throws FileNotFoundException, IOException, Exception {

		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		final CSVReader reader = new CSVReaderBuilder(new BufferedReader(new FileReader(map))).withSkipLines(0).withCSVParser(parser).build();

		String[] nextLine = reader.readNext();

//		if (!nextLine[0].equals("Gene stable ID") || !nextLine[1].equals("Transcript stable ID") || !nextLine[2].equals("UniProtKB Gene Name ID") || !nextLine[3].equals("UniProtKB/Swiss-Prot ID") || !nextLine[4].equals("UniProtKB/TrEMBL ID")) {
//			throw new Exception("Header of uniprot to ensg map should be: \"Gene stable ID[tab]Transcript stable ID[tab]UniProtKB Gene Name ID[tab]UniProtKB/Swiss-Prot ID[tab]UniProtKB/TrEMBL ID\"");
//		}
		if (!nextLine[0].equals("Gene stable ID") || !nextLine[1].equals("HGNC symbol")) {
			throw new Exception("Header of hgnc to ensg map should be: \"Gene stable ID[tab]HGNC symbol\"");
		}

		HashMap<String, ArrayList<String>> hgncToEnsgMap = new HashMap<>(70000);

		while ((nextLine = reader.readNext()) != null) {

			String hgncId = nextLine[1];

			ArrayList<String> hgncEnsgIds = hgncToEnsgMap.get(hgncId);
			if (hgncEnsgIds == null) {
				hgncEnsgIds = new ArrayList<>();
				hgncToEnsgMap.put(hgncId, hgncEnsgIds);
			}

			hgncEnsgIds.add(nextLine[0]);

		}

		return hgncToEnsgMap;

	}

}
