/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
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
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.Map;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

/**
 *
 * @author patri
 */
public class ConvertMgiToMatrix {

	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) throws Exception {

		final File hgncMapFile = new File("D:\\UMCG\\Genetica\\Projects\\GeneNetwork\\ensgHgncIdV98.txt");
		final File geneOrderFile = new File("D:\\UMCG\\Genetica\\Projects\\GeneNetwork\\ensgNcbiIdV98.txt");
		final File mgiMappingFile = new File("D:\\UMCG\\Genetica\\Projects\\GeneNetwork\\MGI\\2020_10_20\\HGNC_homologene.rpt");
		final File mgiPhenoFile = new File("D:\\UMCG\\Genetica\\Projects\\GeneNetwork\\MGI\\2020_10_20\\MGI_PhenoGenoMP.rpt");
		final File outputFile = new File("D:\\UMCG\\Genetica\\Projects\\GeneNetwork\\MGI\\2020_10_20\\MGI_PhenoGenoMP_matrix.txt.gz");
		final File outputFile2 = new File("D:\\UMCG\\Genetica\\Projects\\GeneNetwork\\MGI\\2020_10_20\\MGI_PhenoGenoMP_genesInPathways.txt");

		HashMap<String, String> hgncIdToEnsgMap = hgncMapReader(hgncMapFile);

		HashMap<String, String> mgiGeneHgncIdMap = mgiHgncReader(mgiMappingFile, 15);

		HashMap<String, String> mgiToEnsgMap = new HashMap<>();

		for (Map.Entry<String, String> mgiGeneHgncIdEntry : mgiGeneHgncIdMap.entrySet()) {
			String hgncId = mgiGeneHgncIdEntry.getValue();
			if (hgncId == null || hgncId.equals("")) {
				continue;
			}
			if (hgncIdToEnsgMap.containsKey(hgncId)) {
				mgiToEnsgMap.put(mgiGeneHgncIdEntry.getKey(), hgncIdToEnsgMap.get(hgncId));
			}
		}

		ArrayList<String> geneOrder = readGenes(geneOrderFile);
		LinkedHashSet<String> geneOrder2 = new LinkedHashSet<>(geneOrder);

		System.out.println("Genes with mapping from MGI to ENSG: " + mgiToEnsgMap.size());

		HashMap<String, HashSet<String>> mgiPathwayToGenes = readMgiPhenoFile(mgiPhenoFile, mgiToEnsgMap);

		DoubleMatrixDataset<String, String> gmiMatrix = new DoubleMatrixDataset(geneOrder2, mgiPathwayToGenes.keySet());
		
		HashSet<String> genesWithPathway = new HashSet<>(10000);
		BufferedWriter geneWriter = new BufferedWriter(new FileWriter(outputFile2));

		
		for (Map.Entry<String, HashSet<String>> mgiPathwayToGenesEntry : mgiPathwayToGenes.entrySet()) {
			String gmtPathway = mgiPathwayToGenesEntry.getKey();

			for (String gene : mgiPathwayToGenesEntry.getValue()) {

				if (gmiMatrix.containsRow(gene)) {

					if (genesWithPathway.add(gene)) {
						//add to genes file if not already done
						geneWriter.write(gene);
						geneWriter.write('\n');
					}

					if(!gmiMatrix.containsCol(gmtPathway)){
						System.out.println("missing pathway: " + gmtPathway);
					}
					
					if(!gmiMatrix.containsRow(gene)){
						System.out.println("Missing gene: " + gene);
					}
					
					gmiMatrix.setElement(gene, gmtPathway, 1);
					
				}
			}

		}

		gmiMatrix.save(outputFile);
		geneWriter.close();

		System.out.println("Genes in pathway: " + genesWithPathway.size());


	}

	private static HashMap<String, String> mgiHgncReader(File hgncMapFile, int hgncCol) throws FileNotFoundException, IOException, Exception {

		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		final CSVReader reader = new CSVReaderBuilder(new BufferedReader(new FileReader(hgncMapFile))).withSkipLines(0).withCSVParser(parser).build();

		String[] nextLine = reader.readNext();

		if (!nextLine[0].equals("MGI Accession ID") || !nextLine[hgncCol].equals("HGNC ID")) {
			throw new Exception("Did not find expect columns");
			
		}

		HashMap<String, String> hgncIdToEnsgMap = new HashMap<>(70000);

		while ((nextLine = reader.readNext()) != null) {
			hgncIdToEnsgMap.put(nextLine[0], nextLine[hgncCol]);
		}

		return hgncIdToEnsgMap;

	}

	private static HashMap<String, HashSet<String>> readMgiPhenoFile(File mgiPhenoFile, HashMap<String, String> mgiToEnsgMap) throws Exception {

		final CSVParser mgiParser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		final CSVReader mgiReader = new CSVReaderBuilder(new BufferedReader(new FileReader(mgiPhenoFile))).withSkipLines(1).withCSVParser(mgiParser).build();

		HashMap<String, HashSet<String>> mgiPathwayToGenes = new HashMap<>();

		String[] nextLine;
		while ((nextLine = mgiReader.readNext()) != null) {

			String gmtPathway = nextLine[3];
			String gmtGene = nextLine[5];
			String ensgGene = mgiToEnsgMap.get(gmtGene);
			if (ensgGene == null) {
				continue;
			}

			HashSet<String> hpoGenes = mgiPathwayToGenes.get(gmtPathway);
			if (hpoGenes == null) {
				hpoGenes = new HashSet<>();
				mgiPathwayToGenes.put(gmtPathway, hpoGenes);
			}
			hpoGenes.add(ensgGene);

		}

		return mgiPathwayToGenes;

	}

	private static HashMap<String, String> hgncMapReader(File hgncMapFile) throws FileNotFoundException, IOException, Exception {

		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		final CSVReader reader = new CSVReaderBuilder(new BufferedReader(new FileReader(hgncMapFile))).withSkipLines(0).withCSVParser(parser).build();

		String[] nextLine = reader.readNext();

		if (!nextLine[0].equals("Gene stable ID") || !nextLine[1].equals("HGNC ID")) {
			throw new Exception("Header of ncbi to ensg map should be: \"Gene stable HGNC ID\"");
		}

		HashMap<String, String> hgncIdToEnsgMap = new HashMap<>(70000);

		while ((nextLine = reader.readNext()) != null) {
			hgncIdToEnsgMap.put(nextLine[1], nextLine[0]);
		}

		return hgncIdToEnsgMap;

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
