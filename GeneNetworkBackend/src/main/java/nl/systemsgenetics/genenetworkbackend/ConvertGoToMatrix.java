package nl.systemsgenetics.genenetworkbackend;

import com.opencsv.CSVParser;
import com.opencsv.CSVParserBuilder;
import com.opencsv.CSVReader;
import com.opencsv.CSVReaderBuilder;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.EnumMap;
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
public class ConvertGoToMatrix {

	private static enum GoType {
		F, P, C
	};

	/**
	 * @param args the command line arguments
	 * @throws java.io.IOException
	 * @throws java.lang.Exception
	 */
	public static void main(String[] args) throws IOException, Exception {

		final File pathwayFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\HPO\\135\\ALL_SOURCES_ALL_FREQUENCIES_diseases_to_genes_to_phenotypes.txt");
		final File geneOrderFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\geneOrder.txt");
		final File uniPortToEnsgMapFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\ensgUniProtId.txt");
		

		HashMap<String, ArrayList<String>> uniProtToEnsgMap = loadUniProtToEnsgMap(uniPortToEnsgMapFile);

		EnumMap<GoType, HashMap<String, HashSet<String>>> goTypePathwayToGenes = readPathwayFile(pathwayFile, uniProtToEnsgMap);

		ArrayList<String> geneOrder = readGenes(geneOrderFile);

		
		System.out.println("Genes in order file: " + geneOrder.size());

		for (GoType goType : GoType.values()) {

			final File outputFile = new File(pathwayFile.getAbsolutePath() + "_" + goType.toString() + "_matrix.txt");
			
			HashMap<String, HashSet<String>> pathwayToGenes = goTypePathwayToGenes.get(goType);
			
			System.out.println("Total genesets of" + goType.toString() + ": " + pathwayToGenes.size());
			
			DoubleMatrixDataset<String, String> pathwayMatrix = new DoubleMatrixDataset(geneOrder, pathwayToGenes.keySet());

			for (Map.Entry<String, HashSet<String>> pathwayToGenesEntry : pathwayToGenes.entrySet()) {

				String pathway = pathwayToGenesEntry.getKey();

				for (String gene : pathwayToGenesEntry.getValue()) {

					pathwayMatrix.setElement(gene, pathway, 1);

				}

			}

			pathwayMatrix.save(outputFile);
			
		}

	}

	private static EnumMap<GoType, HashMap<String, HashSet<String>>> readPathwayFile(File pathwayFile, HashMap<String, ArrayList<String>> uniProtToEnsgMap) throws Exception {

		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		final CSVReader reader = new CSVReaderBuilder(new BufferedReader(new FileReader(pathwayFile))).withSkipLines(0).withCSVParser(parser).build();

		EnumMap<GoType, HashMap<String, HashSet<String>>> goTypePathwayToGenes = new EnumMap<>(GoType.class);

		for (GoType goType : GoType.values()) {
			goTypePathwayToGenes.put(goType, new HashMap<>());
		}

		String[] nextLine;
		while ((nextLine = reader.readNext()) != null) {

			if (nextLine[0].charAt(0) == '!') {
				continue;
			}

			//Exclude special tuples like empty
			if (nextLine[4].isEmpty()) {

				String uniProtId = nextLine[1];
				String pathway = nextLine[0];

				ArrayList<String> ensgIds = uniProtToEnsgMap.get(uniProtId);

				if (ensgIds == null) {
					System.err.println("Missing mapping for gene: " + uniProtId + " " + nextLine[1]);
				} else {
					
					HashMap<String, HashSet<String>> pathwayToGenes = goTypePathwayToGenes.get(GoType.valueOf(nextLine[3]));
					
					HashSet<String> pathwayGenes = pathwayToGenes.get(pathway);
					if (pathwayGenes == null) {
						pathwayGenes = new HashSet<>();
						pathwayToGenes.put(pathway, pathwayGenes);
					}

					for (String ensgId : ensgIds) {
						pathwayGenes.add(ensgId);
					}

				}

			}

		}

		return goTypePathwayToGenes;

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

	private static HashMap<String, ArrayList<String>> loadUniProtToEnsgMap(File map) throws FileNotFoundException, IOException, Exception {

		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		final CSVReader reader = new CSVReaderBuilder(new BufferedReader(new FileReader(map))).withSkipLines(0).withCSVParser(parser).build();

		String[] nextLine = reader.readNext();

		if (!nextLine[0].equals("Gene stable ID") || !nextLine[1].equals("UniProtKB Gene Name ID")) {
			throw new Exception("Header of ncbi to ensg map should be: \"Gene stable ID[tab]UniProtKB Gene Name ID\"");
		}

		HashMap<String, ArrayList<String>> uniProtToEnsgMap = new HashMap<>(70000);

		while ((nextLine = reader.readNext()) != null) {

			String uniProtId = nextLine[1];

			ArrayList<String> uniProtEnsgIds = uniProtToEnsgMap.get(uniProtId);
			if (uniProtEnsgIds == null) {
				uniProtEnsgIds = new ArrayList<>();
			}

			uniProtEnsgIds.add(nextLine[0]);

		}

		return uniProtToEnsgMap;

	}

}
