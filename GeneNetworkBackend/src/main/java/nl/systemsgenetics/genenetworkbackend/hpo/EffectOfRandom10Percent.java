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
import gnu.trove.map.TObjectDoubleMap;
import gnu.trove.map.hash.TObjectDoubleHashMap;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.text.ParseException;
import java.util.HashSet;
import java.util.LinkedHashSet;
import org.biojava.nbio.ontology.Ontology;
import org.biojava.nbio.ontology.Term;
import org.biojava.nbio.ontology.Triple;

/**
 *
 * @author patri
 */
public class EffectOfRandom10Percent {

	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) throws IOException, FileNotFoundException, ParseException {

		final File predictedHpoTermFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Data31995Genes05-12-2017\\PCA_01_02_2018\\predictions\\hpo_predictions_auc_bonferroni.txt");

		TObjectDoubleMap<String> hpoAuc = readPredictedHpoTermFile(predictedHpoTermFile);

		

	}

	public static TObjectDoubleMap<String> readPredictedHpoTermFile(File predictedHpoTermFile) throws FileNotFoundException, IOException {

		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		final CSVReader reader = new CSVReaderBuilder(new BufferedReader(new FileReader(predictedHpoTermFile))).withSkipLines(1).withCSVParser(parser).build();

		TObjectDoubleMap<String> hpos = new TObjectDoubleHashMap<>();

		String[] nextLine;
		while ((nextLine = reader.readNext()) != null) {

			hpos.put(nextLine[0], Double.parseDouble(nextLine[3]));

		}

		reader.close();

		return hpos;

	}

}
