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
import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Objects;

/**
 *
 * @author patri
 */
public class InvestigateString {

	/**
	 * @param args the command line arguments
	 * @throws java.io.FileNotFoundException
	 */
	public static void main(String[] args) throws FileNotFoundException, IOException, Exception {

		final File stringFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\String\\9606.protein.links.full.v10.5.txt");
		final File ensgEnspFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\ensgGeneProteinV83.txt");

		HashMap<String, String> enspToEnsgMap = loadEnspToEnsgMap(ensgEnspFile);

		final CSVParser parser = new CSVParserBuilder().withSeparator(' ').withIgnoreQuotations(true).build();
		final CSVReader sampleFileReader = new CSVReaderBuilder(new BufferedReader(new FileReader(stringFile))).withSkipLines(1).withCSVParser(parser).build();

		HashSet<GenePair> observedGenePairs = new HashSet<>();
		TObjectIntMap<String> geneCombinedConfidentCount = new TObjectIntHashMap<>();
		TObjectIntMap<String> geneCombinedCoexpressionCount = new TObjectIntHashMap<>();

		String[] nextLine;
		while ((nextLine = sampleFileReader.readNext()) != null) {

			final String protein1 = nextLine[0].substring(5);
			final String protein2 = nextLine[1].substring(5);

			final String gene1 = enspToEnsgMap.get(protein1);
			final String gene2 = enspToEnsgMap.get(protein2);
		

			if (gene1 == null || gene2 == null) {
				continue;
			}

			GenePair gene1Gene2Pair = new GenePair(gene1, gene2);
			if (observedGenePairs.add(gene1Gene2Pair)) {

				if (Integer.parseInt(nextLine[15]) >= 700) {
					geneCombinedConfidentCount.adjustOrPutValue(gene1, 1, 1);
					geneCombinedConfidentCount.adjustOrPutValue(gene2, 1, 1);
				}
				
				if (Integer.parseInt(nextLine[8]) >= 700) {
					geneCombinedCoexpressionCount.adjustOrPutValue(gene1, 1, 1);
					geneCombinedCoexpressionCount.adjustOrPutValue(gene2, 1, 1);
				}

			} else {
//				System.out.println("Observed duplicate pair");
//				System.out.println(protein1);
//				System.out.println(protein2);
//				System.out.println(gene1);
//				System.out.println(gene2);

			}

		}
		
		System.out.println("Gene\tCoexpressionCount\tCombinedCount");
		for(String gene : geneCombinedConfidentCount.keySet()){
			System.out.println(gene + "\t" + geneCombinedCoexpressionCount.get(gene) + "\t" + geneCombinedConfidentCount.get(gene)) ;
		}

		System.out.println("Total unique interactions: " + observedGenePairs.size());

	}

	private static HashMap<String, String> loadEnspToEnsgMap(File mapFile) throws FileNotFoundException, IOException, Exception {

		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		final CSVReader reader = new CSVReaderBuilder(new BufferedReader(new FileReader(mapFile))).withSkipLines(0).withCSVParser(parser).build();

		String[] nextLine = reader.readNext();

		if (!nextLine[0].equals("Ensembl Gene ID") || !nextLine[1].equals("Ensembl Protein ID")) {
			throw new Exception("Header should be: \"Gene stable ID	NCBI gene ID\"");
		}

		HashMap<String, String> enspToEnsgMap = new HashMap<>(70000);

		while ((nextLine = reader.readNext()) != null) {

			if (!nextLine[1].equals("")) {

				enspToEnsgMap.put(nextLine[1], nextLine[0]);

			}
		}

		return enspToEnsgMap;

	}

	private static final class GenePair {
		
		final String g1;
		final String g2;

		public GenePair(String g1, String g2) {
			this.g1 = g1;
			this.g2 = g2;
		}

		@Override
		public int hashCode() {
			int hash = Objects.hashCode(this.g1) * Objects.hashCode(this.g2);
			return hash;
		}

		@Override
		public boolean equals(Object obj) {
			if (this == obj) {
				return true;
			}
			if (obj == null) {
				return false;
			}
			if (getClass() != obj.getClass()) {
				return false;
			}
			final GenePair other = (GenePair) obj;
			if (Objects.equals(this.g1, other.g1) && Objects.equals(this.g2, other.g2)) {
				return true;
			}
			if (Objects.equals(this.g1, other.g2) && Objects.equals(this.g2, other.g1)) {
				return true;
			}
			return false;
		}
		
		
		
	}

}
