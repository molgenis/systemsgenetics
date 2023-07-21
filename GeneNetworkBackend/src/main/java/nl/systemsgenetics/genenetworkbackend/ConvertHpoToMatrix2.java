package nl.systemsgenetics.genenetworkbackend;

import com.opencsv.CSVParser;
import com.opencsv.CSVParserBuilder;
import com.opencsv.CSVReader;
import com.opencsv.CSVReaderBuilder;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.zip.GZIPInputStream;
import static nl.systemsgenetics.genenetworkbackend.ConvertGeneList.readGenes;

import umcg.genetica.math.matrix2.DoubleMatrixDataset;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 * @author patri
 */
public class ConvertHpoToMatrix2 {

	/**
	 * @param args the command line arguments
	 * @throws java.io.IOException
	 * @throws java.lang.Exception
	 */
	public static void main(String[] args) throws IOException, Exception {

		//final File hpoFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\HPO\\135\\ALL_SOURCES_ALL_FREQUENCIES_diseases_to_genes_to_phenotypes.txt");

		final File genesFile = new File("D:\\UMCG\\Genetica\\Projects\\GeneNetwork\\genes_Ensembl94_protein_coding.txt");
//		final File geneOrderFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Data31995Genes05-12-2017\\PCA_01_02_2018\\genes.txt");
//		
		final File hpoFile = new File("D:\\UMCG\\Genetica\\Projects\\GeneNetwork\\HPO\\2023_06_17\\phenotype_to_genes.txt");
//		final File outputFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Data31995Genes05-12-2017\\PCA_01_02_2018\\PathwayMatrix\\" + hpoFile.getName() + "_matrix.txt.gz");
//		final File outputFile2 = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Data31995Genes05-12-2017\\PCA_01_02_2018\\PathwayMatrix\\" + hpoFile.getName() + "_genesInPathways.txt");

		//final File geneOrderFile = new File("C:\\UMCG\\Genetica\\Projects\\Depict2Pgs\\testPredictions\\mergedGeneNetworkCoregBoth.rows.txt");
		final File outputFile = new File("D:\\UMCG\\Genetica\\Projects\\GeneNetwork\\HPO\\2023_06_17\\HPO_2023_06_17.txt.gz");
		final File outputFile2 = new File("D:\\UMCG\\Genetica\\Projects\\GeneNetwork\\HPO\\2023_06_17\\HPO_2023_06_17_genesInPathways.txt");


//		final File hpoFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\HPO\\135\\bavWithoutWilliams.txt");
//		final File outputFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Data31995Genes05-12-2017\\PCA_01_02_2018\\PathwayMatrix\\" + hpoFile.getName() + "_matrix.txt.gz");
//		final File outputFile2 = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Data31995Genes05-12-2017\\PCA_01_02_2018\\PathwayMatrix\\" + hpoFile.getName() + "_genesInPathways.txt");


		List<Gene> genes = readGenes(genesFile);
		
		ArrayList<String> geneIds = new ArrayList();
		HashMap<String, String> symbolGeneMap = new HashMap<>();
		for(Gene gene : genes){
			geneIds.add(gene.getGene());
			symbolGeneMap.put(gene.getGeneSymbol(), gene.getGene());
		}

		HashMap<String, HashSet<String>> hpoToGenes = readHpoFile(hpoFile, symbolGeneMap);


		System.out.println("Total HPO terms: " + hpoToGenes.size());

		DoubleMatrixDataset<String, String> hpoMatrix = new DoubleMatrixDataset(geneIds, hpoToGenes.keySet());

		HashSet<String> genesWithHpo = new HashSet<>(10000);
		BufferedWriter geneWriter = new BufferedWriter(new FileWriter(outputFile2));

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
		hpoMatrix.save(outputFile);

		System.out.println("Genes in pathway: " + genesWithHpo.size());

	}

	public static HashMap<String, ArrayList<String>> loadNcbiToEnsgMap(File ncbiToEnsgMapFile) throws FileNotFoundException, IOException, Exception {

		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		CSVReader reader = null;

		if (ncbiToEnsgMapFile.getName().endsWith(".gz")) {
			reader = new CSVReaderBuilder(new BufferedReader(new InputStreamReader((new GZIPInputStream(new FileInputStream(ncbiToEnsgMapFile)))))).withSkipLines(0).withCSVParser(parser).build();
		} else {
			reader = new CSVReaderBuilder(new BufferedReader(new FileReader(ncbiToEnsgMapFile))).withSkipLines(0).withCSVParser(parser).build();
		}

		String[] nextLine = reader.readNext();

		if (!nextLine[0].equals("Gene stable ID") || !nextLine[1].equals("NCBI gene ID")) {
			throw new Exception("Header of ncbi to ensg map should be: \"Gene stable ID	NCBI gene ID\"");
		}

		HashMap<String, ArrayList<String>> ncbiToEnsgMap = new HashMap<>(70000);

		while ((nextLine = reader.readNext()) != null) {

			String ncbiId = nextLine[1];

			ArrayList<String> ncbiEnsgIds = ncbiToEnsgMap.get(ncbiId);
			if (ncbiEnsgIds == null) {
				ncbiEnsgIds = new ArrayList<>();
				ncbiToEnsgMap.put(ncbiId, ncbiEnsgIds);
			}

			ncbiEnsgIds.add(nextLine[0]);

		}

		return ncbiToEnsgMap;

	}

	private static HashMap<String, HashSet<String>> readHpoFile(File hpoFile, HashMap<String, String> symbolGeneMap) throws Exception {

		final CSVParser hpoParser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		CSVReader hpoReader = null;
		if (hpoFile.getName().endsWith(".gz")) {
			hpoReader = new CSVReaderBuilder(new BufferedReader(new InputStreamReader((new GZIPInputStream(new FileInputStream(hpoFile)))))).withSkipLines(1).withCSVParser(hpoParser).build();
		} else {
			hpoReader = new CSVReaderBuilder(new BufferedReader(new FileReader(hpoFile))).withSkipLines(1).withCSVParser(hpoParser).build();
		}


		HashMap<String, HashSet<String>> hpoToGenes = new HashMap<>();

		String[] nextLine;
		while ((nextLine = hpoReader.readNext()) != null) {
			String hpo = nextLine[0];
			
			if(hpo.equals("HP:0000001")){
				continue;
			}
			
			String hgcnId = nextLine[3];
			String ensgId = symbolGeneMap.get(hgcnId);
			
			if (ensgId == null) {
				System.err.println("Missing mapping for gene: " + hgcnId);
			} else {

				HashSet<String> hpoGenes = hpoToGenes.get(hpo);
				if (hpoGenes == null) {
					hpoGenes = new HashSet<>();
					hpoToGenes.put(hpo, hpoGenes);
				}
				hpoGenes.add(ensgId);
				

			}

		}

		return hpoToGenes;

	}

//	private static ArrayList<String> readGenes(File geneOrderFile) throws IOException {
//
//		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
//		CSVReader reader = null;
//		if (geneOrderFile.getName().endsWith(".gz")) {
//			reader = new CSVReaderBuilder(new BufferedReader(new InputStreamReader((new GZIPInputStream(new FileInputStream(geneOrderFile)))))).withSkipLines(0).withCSVParser(parser).build();
//		} else {
//			reader = new CSVReaderBuilder(new BufferedReader(new FileReader(geneOrderFile))).withSkipLines(0).withCSVParser(parser).build();
//		}
//
//
//		String[] nextLine;
//		ArrayList<String> geneOrder = new ArrayList<>();
//
//		while ((nextLine = reader.readNext()) != null) {
//
//			geneOrder.add(nextLine[0]);
//
//		}
//
//		return geneOrder;
//
//	}

	public static HashMap<String, ArrayList<String>> loadHgncToEnsgMap(File map) throws FileNotFoundException, IOException, Exception {

		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		CSVReader reader = null;
		if (map.getName().endsWith(".gz")) {
			reader = new CSVReaderBuilder(new BufferedReader(
					new InputStreamReader((new GZIPInputStream(new FileInputStream(map))))
			)).withSkipLines(0).withCSVParser(parser).build();
		} else {
			reader = new CSVReaderBuilder(new BufferedReader(new FileReader(map))).withSkipLines(0).withCSVParser(parser).build();
		}


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
