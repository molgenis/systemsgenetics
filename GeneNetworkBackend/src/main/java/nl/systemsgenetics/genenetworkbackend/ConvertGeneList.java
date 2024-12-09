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

import umcg.genetica.math.matrix2.DoubleMatrixDataset;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 * @author patri
 */
public class ConvertGeneList {

	/**
	 * @param args the command line arguments
	 * @throws java.io.IOException
	 * @throws java.lang.Exception
	 */
	public static void main(String[] args) throws IOException, Exception {

		//final File hpoFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\HPO\\135\\ALL_SOURCES_ALL_FREQUENCIES_diseases_to_genes_to_phenotypes.txt");
		final File genesFile = new File("D:\\UMCG\\Genetica\\Projects\\GeneNetwork\\genes_Ensembl94_protein_coding.txt");

		final File geneFile = new File("D:\\UMCG\\Genetica\\Projects\\GeneNetwork\\UCU\\HLHS.txt");
		final File outputFile = new File("D:\\UMCG\\Genetica\\Projects\\GeneNetwork\\UCU\\" + geneFile.getName() + "_matrix.txt.gz");

		String setName = geneFile.getName();

		HashSet<String> clusterToGenes = readGeneList(geneFile);

		List<Gene> genes = readGenes(genesFile);
		
		ArrayList<String> geneIds = new ArrayList();
		HashMap<String, String> symbolGeneMap = new HashMap<>();
		for(Gene gene : genes){
			geneIds.add(gene.getGene());
			symbolGeneMap.put(gene.getGeneSymbol(), gene.getGene());
		}

		System.out.println("Total HPO terms: " + clusterToGenes.size());
		System.out.println("Genes: " + genes.size());

		HashSet<String> colnames = new HashSet<>(1);
		colnames.add(setName);

		DoubleMatrixDataset<String, String> hpoMatrix = new DoubleMatrixDataset(geneIds, colnames);

		for (String geneSymbol : clusterToGenes) {

			
			
			//String geneId = symbolGeneMap.get(geneSymbol);
			String geneId = geneSymbol;//in case file contains ensg ids
			
			System.out.println(geneSymbol + " - " + geneId);
			

			if (hpoMatrix.containsRow(geneId)) {
				hpoMatrix.setElement(geneId, setName, 1);
			}

		}

		System.out.println(hpoMatrix.getCol(0).cardinality());
		
		hpoMatrix.save(outputFile);

	}

	public static List<Gene> readGenes(File geneFile) throws IOException {

		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		CSVReader reader = null;
		if (geneFile.getName().endsWith(".gz")) {
			reader = new CSVReaderBuilder((new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(geneFile)))))).withSkipLines(1).withCSVParser(parser).build();
		} else {
			reader = new CSVReaderBuilder(new BufferedReader(new FileReader(geneFile))).withCSVParser(parser).withSkipLines(1).build();
		}
		final ArrayList<Gene> genes = new ArrayList<>();

		String[] nextLine;
		while ((nextLine = reader.readNext()) != null) {

			genes.add(new Gene(nextLine[0], nextLine[1], Integer.parseInt(nextLine[2]), Integer.parseInt(nextLine[3]), nextLine[5], nextLine[6]));

		}

		return genes;

	}

	private static HashSet<String> readGeneList(File hpoFile) throws Exception {

		final CSVParser hpoParser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		CSVReader hpoReader = null;
		if (hpoFile.getName().endsWith(".gz")) {
			hpoReader = new CSVReaderBuilder(new BufferedReader(
					new InputStreamReader((new GZIPInputStream(new FileInputStream(hpoFile))))
			)).withSkipLines(0).withCSVParser(hpoParser).build();
		} else {
			hpoReader = new CSVReaderBuilder(new BufferedReader(new FileReader(hpoFile))).withSkipLines(0).withCSVParser(hpoParser).build();
		}

		HashSet<String> hpoGenes = new HashSet<>();

		String[] nextLine;
		while ((nextLine = hpoReader.readNext()) != null) {
			String gene = nextLine[0];

			hpoGenes.add(gene);

		}

		return hpoGenes;

	}

	private static ArrayList<String> readGenes2(File geneOrderFile) throws IOException {

		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		CSVReader reader = null;
		if (geneOrderFile.getName().endsWith(".gz")) {
			reader = new CSVReaderBuilder(new BufferedReader(new InputStreamReader((new GZIPInputStream(new FileInputStream(geneOrderFile)))))).withSkipLines(0).withCSVParser(parser).build();
		} else {
			reader = new CSVReaderBuilder(new BufferedReader(new FileReader(geneOrderFile))).withSkipLines(0).withCSVParser(parser).build();
		}
		String[] nextLine;
		ArrayList<String> geneOrder = new ArrayList<>();

		while ((nextLine = reader.readNext()) != null) {

			geneOrder.add(nextLine[0]);

		}

		return geneOrder;

	}

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
