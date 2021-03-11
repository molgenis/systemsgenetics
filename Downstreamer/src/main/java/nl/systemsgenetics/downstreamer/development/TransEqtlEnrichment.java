/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.downstreamer.development;

import cern.jet.math.tdouble.DoubleFunctions;
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
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.RandomAccessGenotypeDataReaderFormats;
import org.molgenis.genotype.util.LdCalculatorException;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variantFilter.VariantIdIncludeFilter;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import static umcg.genetica.math.matrix2.DoubleMatrixDataset.loadDoubleBinaryData;

/**
 *
 * @author patri
 */
public class TransEqtlEnrichment {

	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) throws Exception {

		String transEqtlMatrixPath = "D:\\UMCG\\Genetica\\Projects\\Depict2Pgs\\eqtlgen\\ZScoreMatrix";
		File gwasFile = new File("D:\\UMCG\\Genetica\\Projects\\Depict2Pgs\\eqtlgen\\gwas.txt");
		String outputPath = "D:\\UMCG\\Genetica\\Projects\\Depict2Pgs\\eqtlgen\\SumChi2PerGwas.txt";

		if (false) {
			HashSet<String> colsToExclude = new HashSet<>();
			colsToExclude.add("Alleles");
			colsToExclude.add("AlleleAssessed");

			DoubleMatrixDataset<String, String> transEqtlDataset = DoubleMatrixDataset.loadDoubleTextDoubleDataExlcudeCols(transEqtlMatrixPath + ".txt.gz", '\t', colsToExclude);

			fixHeader(transEqtlDataset);
			//dont do this all rows contain some NaN values
			//transEqtlDataset = removeNaNRows(transEqtlDataset);
			transEqtlDataset.saveBinary(transEqtlMatrixPath);
		}

		DoubleMatrixDataset<String, String> transEqtlDataset = loadDoubleBinaryData(transEqtlMatrixPath);
		System.out.println("trans matrix " + transEqtlDataset.rows() + " x " + transEqtlDataset.columns());

		for(int r = 0 ; r < transEqtlDataset.rows() ; ++r){
			for(int c = 0 ; c < transEqtlDataset.columns(); ++c){
				if(Double.isNaN(transEqtlDataset.getElementQuick(r, c))){
					transEqtlDataset.setElementQuick(r, c, 0);
				}
			}
		}	
		
		HashMap<String, HashSet<String>> gwasVariants = readGwasFile(gwasFile);

		HashSet<String> allGwasVariants = new HashSet<>();
		for (HashSet<String> variants : gwasVariants.values()) {
			allGwasVariants.addAll(variants);
		}

		RandomAccessGenotypeData referenceGenotypeData = RandomAccessGenotypeDataReaderFormats.VCF_FOLDER.createFilteredGenotypeData("D:\\UMCG\\Genetica\\Projects\\Depict2Pgs\\reference_datasets\\human_b37\\VCF_1Kgp1v3_EUR\\", 20000, new VariantIdIncludeFilter(allGwasVariants), null);

		System.out.println(gwasVariants.get("Inflammatory bowel disease").size());
		prune(gwasVariants, referenceGenotypeData.getVariantIdMap());
		System.out.println(gwasVariants.get("Inflammatory bowel disease").size());

		DoubleMatrixDataset<String, String> sumChi2Dataset = new DoubleMatrixDataset<String, String>(transEqtlDataset.getHashCols().keySet(), gwasVariants.keySet());

		for (String trait : gwasVariants.keySet()) {
			
			System.out.println(trait);
			
			HashSet<String> variants = gwasVariants.get(trait);

			variants.retainAll(transEqtlDataset.getRowObjects());

			DoubleMatrixDataset<String, String> transEffectOfGwasVariants = transEqtlDataset.viewRowSelection(variants);

			for (String gene : transEffectOfGwasVariants.getHashCols().keySet()) {
				double sumchi2 = transEffectOfGwasVariants.getCol(gene).aggregate(DoubleFunctions.plus, DoubleFunctions.square);
				sumChi2Dataset.setElement(gene, trait, sumchi2);
			}

		}

		sumChi2Dataset.save(outputPath);

	}

	private static HashMap<String, HashSet<String>> readGwasFile(File gwas) throws FileNotFoundException, IOException {

		final CSVParser gwasParser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		final CSVReader gwasReader = new CSVReaderBuilder(new BufferedReader(new FileReader(gwas))).withSkipLines(1).withCSVParser(gwasParser).build();

		HashMap<String, HashSet<String>> gwasVariants = new HashMap<>();

		String[] inputLine;
		while ((inputLine = gwasReader.readNext()) != null) {

			String trait = inputLine[8];
			String snp = inputLine[0];

			HashSet<String> variants = gwasVariants.get(trait);
			if (variants == null) {
				variants = new HashSet<>();
				gwasVariants.put(trait, variants);
			}
			variants.add(snp);
		}

		for (Map.Entry<String, HashSet<String>> e : gwasVariants.entrySet()) {
			System.out.println(e.getKey() + "\t" + e.getValue().size());
		}

		return gwasVariants;

	}

	public static void fixHeader(DoubleMatrixDataset dataset) throws Exception {

		LinkedHashMap<String, Integer> oldHash = dataset.getHashCols();
		LinkedHashMap<String, Integer> newHash = new LinkedHashMap<>(oldHash.size());

		for (Map.Entry<String, Integer> oldEntry : oldHash.entrySet()) {

			String oldGeneName = oldEntry.getKey();
			int indexOfPoint = oldGeneName.indexOf('_');

			if (indexOfPoint < 0) {
				if (newHash.put(oldGeneName, oldEntry.getValue()) != null) {
					throw new Exception("Can't trim gene names if this causes duplicate genes: " + oldGeneName);
				}
			} else {
				if (newHash.put(oldGeneName.substring(0, indexOfPoint), oldEntry.getValue()) != null) {
					throw new Exception("Can't trim gene names if this causes duplicate genes: " + oldGeneName);
				}
			}

		}

		dataset.setHashCols(newHash);

	}

	private static DoubleMatrixDataset removeNaNRows(DoubleMatrixDataset dataset) {

		ArrayList<String> rowNames = dataset.getRowObjects();
		ArrayList<String> nonNanRowNames = new ArrayList<>(dataset.rows());

		rows:
		for (int r = 0; r < dataset.rows(); ++r) {
			boolean nonZeroValue = false;
			double e;
			for (int c = 0; c < dataset.columns(); ++c) {
				e = dataset.getElementQuick(r, c);
				if (Double.isNaN(e)) {
					continue rows;
				} else if (e != 0) {
					nonZeroValue = true;
				}
			}
			//if here not NaN values;
			if (nonZeroValue) {
				nonNanRowNames.add(rowNames.get(r));
			}

		}

		if (nonNanRowNames.size() < rowNames.size()) {
			dataset = dataset.viewRowSelection(nonNanRowNames);
			System.out.println("Removing " + (rowNames.size() - nonNanRowNames.size()) + " rows with only zero's or NaN");
		}

		return dataset;

	}

	public static void prune(final HashMap<String, HashSet<String>> gwasVariants, final HashMap<String, GeneticVariant> varMap) throws IOException {

		final double variantLdThreshold = 0.1;

		gwasVariants.values().parallelStream().forEach((HashSet<String> variants) -> {

			ArrayList<String> variantsList = new ArrayList<>(variants);

			variants:
			for (int i = 0; i < variantsList.size(); ++i) {
				String var1Name = variantsList.get(i);
				GeneticVariant var1 = varMap.get(var1Name);
				if (var1 == null) {
					//System.out.println(var1Name);
					variants.remove(var1Name);
					continue variants;
				}
				if (!var1.isBiallelic()) {
					variants.remove(var1Name);
					continue variants;
				}

				for (int j = i + i; j < variantsList.size(); ++j) {
					String var2Name = variantsList.get(i);
					GeneticVariant var2 = varMap.get(var2Name);

					try {
						if (var2.isBiallelic() && var1.calculateLd(var2).getR2() >= variantLdThreshold) {
							//now remove var1
							variants.remove(var1Name);
							continue variants;
						}
					} catch (LdCalculatorException ex) {
						throw new RuntimeException(ex);
					}

				}

			}

		});

	}

}
