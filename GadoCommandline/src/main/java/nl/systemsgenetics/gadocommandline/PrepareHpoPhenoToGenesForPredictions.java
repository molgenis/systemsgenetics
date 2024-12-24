/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.gadocommandline;

import com.opencsv.CSVParser;
import com.opencsv.CSVParserBuilder;
import com.opencsv.CSVReader;
import com.opencsv.CSVReaderBuilder;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.zip.GZIPInputStream;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

/**
 * Converts the HPO file "phenotype_to_genes.txt" to a matrix that can be used
 * in for gene network predictions. Also converts the gene symbols to ensg IDs.
 *
 *
 * @author Patrick Deelen
 */
public class PrepareHpoPhenoToGenesForPredictions {

	public static void prioritize(GadoOptions options) throws IOException, ParseException, Exception {

		Map<String, String> symbolGeneMap = loadHgncToEnsg(options.getGenesFile());
		
		ArrayList<String> geneIds = new ArrayList();
		for(String ensgGene : symbolGeneMap.values()){
			geneIds.add(ensgGene);
		}
		
		HashMap<String, HashSet<String>> hpoToGenes = readHpoFile(options.getPhenotypeToGenes(), symbolGeneMap);


		System.out.println("Total HPO terms: " + hpoToGenes.size());

		DoubleMatrixDataset<String, String> hpoMatrix = new DoubleMatrixDataset(geneIds, hpoToGenes.keySet());

		HashSet<String> genesWithHpo = new HashSet<>(10000);
		BufferedWriter geneWriter = new BufferedWriter(new FileWriter(options.getOutputBasePath() + "_genesInPathways.txt"));

		for (Map.Entry<String, HashSet<String>> hpoToGenesEntry : hpoToGenes.entrySet()) {

			String hpo = hpoToGenesEntry.getKey();

			for (String gene : hpoToGenesEntry.getValue()) {

				if (hpoMatrix.containsRow(gene)) {
					if (genesWithHpo.add(gene)) {
						//add to genes file if not already done
						geneWriter.write(gene);
						geneWriter.write('\n');
					}
					hpoMatrix.setElement(gene, hpo, 1);
				}

			}

		}

		geneWriter.close();
		hpoMatrix.save(options.getOutputBasePath() + ".txt.gz");

		System.out.println("Genes in pathway: " + genesWithHpo.size());

		
	}

	private static HashMap<String, HashSet<String>> readHpoFile(File hpoFile, Map<String, String> symbolGeneMap) throws Exception {

		final CSVParser hpoParser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		CSVReader hpoReader = null;
		if (hpoFile.getName().endsWith(".gz")) {
			hpoReader = new CSVReaderBuilder(new BufferedReader(new InputStreamReader((new GZIPInputStream(new FileInputStream(hpoFile)))))).withSkipLines(1).withCSVParser(hpoParser).build();
		} else {
			hpoReader = new CSVReaderBuilder(new BufferedReader(new FileReader(hpoFile))).withSkipLines(1).withCSVParser(hpoParser).build();
		}


		HashMap<String, HashSet<String>> hpoToGenes = new HashMap<>();

		HashSet<String> genesNotMapped = new HashSet<>();

		String[] nextLine;
		while ((nextLine = hpoReader.readNext()) != null) {
			String hpo = nextLine[0];
			
			if(hpo.equals("HP:0000001")){
				continue;
			}
			
			String hgcnId = nextLine[3];
			String ensgId = symbolGeneMap.get(hgcnId);
			
			if (ensgId == null) {
				genesNotMapped.add(hgcnId);
			} else {

				HashSet<String> hpoGenes = hpoToGenes.get(hpo);
				if (hpoGenes == null) {	
					hpoGenes = new HashSet<>();
					hpoToGenes.put(hpo, hpoGenes);
				}
				hpoGenes.add(ensgId);
				

			}

		}
		
		System.err.println("Number of genes that could not be mapped to ensg IDs: " + genesNotMapped.size());
		if(hpoToGenes.size() > 0){
			System.err.println("A few 100 are expected when using only protein coding genes.");
		}

		return hpoToGenes;

	}

	
	private static Map<String, String> loadHgncToEnsg(File mappingFile) throws IOException {

		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		CSVReader reader = null;
		if (mappingFile.getName().endsWith(".gz")) {
			reader = new CSVReaderBuilder(new BufferedReader(
					new InputStreamReader((new GZIPInputStream(new FileInputStream(mappingFile))))
			)).withSkipLines(1).withCSVParser(parser).build();
		} else {
			reader = new CSVReaderBuilder(new BufferedReader(new FileReader(mappingFile))).withSkipLines(1).withCSVParser(parser).build();
		}

		HashMap<String, String> mapping = new HashMap<>();

		String[] nextLine;
		while ((nextLine = reader.readNext()) != null) {

			mapping.put(nextLine[1], nextLine[0]);

		}

		return Collections.unmodifiableMap(mapping);

	}

}
