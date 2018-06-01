/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.genenetworkbackend.hpo;

import com.opencsv.CSVParser;
import com.opencsv.CSVParserBuilder;
import com.opencsv.CSVReader;
import com.opencsv.CSVReaderBuilder;
import gnu.trove.map.hash.TObjectDoubleHashMap;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;

/**
 *
 * @author patri
 */
public class SkewnessInfo {

	private final static String[] EXPECTED_HEADER = {"gene", "maxSkewnessExHpo", "meanSkewnessExHpo", "hpoSkewness"};
	
	private final TObjectDoubleHashMap<String> hpoSkewnessMap = new TObjectDoubleHashMap<>();
	private final TObjectDoubleHashMap<String> maxSkewnessExHpoMap = new TObjectDoubleHashMap<>();
	private final TObjectDoubleHashMap<String> meanSkewnessExHpoMap = new TObjectDoubleHashMap<>();

	public SkewnessInfo(File skewnessFile) throws FileNotFoundException, IOException {

		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		final CSVReader reader = new CSVReaderBuilder(new BufferedReader(new FileReader(skewnessFile))).withSkipLines(0).withCSVParser(parser).build();

		String[] nextLine = reader.readNext();
		
		//check header
		if(!Arrays.equals(nextLine, EXPECTED_HEADER)){
			throw new RuntimeException("skewness header different than expected, is: " + String.join("\t", nextLine)); 
		}
		
		
		while ((nextLine = reader.readNext()) != null) {
			String gene = nextLine[0];
			
			hpoSkewnessMap.put(gene, Double.valueOf(nextLine[1]));
			maxSkewnessExHpoMap.put(gene, Double.valueOf(nextLine[2]));
			meanSkewnessExHpoMap.put(gene, Double.valueOf(nextLine[3]));
			
		}
		
		reader.close();

		
	}

	public double getHpoSkewness(String gene) {
		return hpoSkewnessMap.get(gene);
	}

	public double getMaxSkewnessExHpo(String gene) {
		return maxSkewnessExHpoMap.get(gene);
	}

	public double getMeanSkewnessExHpo(String gene) {
		return meanSkewnessExHpoMap.get(gene);
	}

}
