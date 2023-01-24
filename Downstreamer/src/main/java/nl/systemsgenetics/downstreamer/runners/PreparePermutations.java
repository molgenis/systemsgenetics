package nl.systemsgenetics.downstreamer.runners;

import cern.colt.matrix.tdouble.DoubleMatrix2D;
import com.opencsv.CSVParserBuilder;
import com.opencsv.CSVReader;
import com.opencsv.CSVReaderBuilder;
import gnu.trove.map.hash.TObjectIntHashMap;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.zip.GZIPInputStream;
import nl.systemsgenetics.downstreamer.gene.Gene;
import nl.systemsgenetics.downstreamer.io.IoUtils;
import nl.systemsgenetics.downstreamer.runners.options.OptionsModePerparePermutations;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.matrix2.DoubleMatrixDatasetRowCompressedWriter;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author patri
 */
public class PreparePermutations {

	private static final Logger LOGGER = LogManager.getLogger(PreparePermutations.class);

	public static void run(OptionsModePerparePermutations options) throws IOException {

		// Load the genes to run the analysis on
		final LinkedHashMap<String, List<Gene>> genes = IoUtils.readGenesAsChrMap(options.getGeneInfoFile());

		final File permutationFolder = options.getPermutationFolder();

		final TObjectIntHashMap<String> chunkInfo = readChunks(options.getChunkFile());

		final List<String> permutationNames = readRowColNames(new File(permutationFolder, "NullGwasX_1_1.cols.txt"));
		final int nrPermutations = permutationNames.size();

		//gwasVariants.values().parallelStream().forEach((HashSet<String> variants) -> {
		chunkInfo.keySet().parallelStream().forEach((String chr) -> {
			//for (String chr : chunkInfo.keySet()) {

			try {

				LinkedHashSet<String> chrGeneNames = new LinkedHashSet<>();
				int nonNanGeneCount = 0;
				HashMap<String, ArrayList<String>> nonNanGenesPerArm = new HashMap<>();
				HashMap<String, String> geneArmMap = new HashMap<>();

				for (Gene gene : genes.get(chr)) {
					chrGeneNames.add(gene.getGene());
					String arm = gene.getChrAndArm();
					if (!nonNanGenesPerArm.containsKey(arm)) {
						nonNanGenesPerArm.put(arm, new ArrayList<>());
					}
					geneArmMap.put(gene.getGene(), gene.getChrAndArm());
				}

				int totalGenes = chrGeneNames.size();

				DoubleMatrixDataset<String, String> chrGenePvalues = new DoubleMatrixDataset<>(permutationNames, chrGeneNames);

				final int nrChunks = chunkInfo.get(chr);

				for (int c = 1; c <= nrChunks; ++c) {

					final List<String> chunkGenes = readRowColNames(new File(permutationFolder, "NullGwasX_" + chr + "_" + String.valueOf(c) + ".rows.txt"));

					if (!chrGeneNames.containsAll(chunkGenes)) {
						for (String gene : chunkGenes) {
							if (!chrGeneNames.contains(gene)) {
								System.out.println(gene);
							}
						}
						throw new RuntimeException("Found genes in permutation not in gene file");
					}

					chrGeneNames.removeAll(chunkGenes);

					DoubleMatrix2D chunkGenePvalues = chrGenePvalues.viewColSelection(chunkGenes).getMatrix();

					File chunkPvalueFile = new File(permutationFolder, "genePvaluesNullGwasX_" + chr + "_" + String.valueOf(c) + ".tsv");

					final CSVReader reader
							= new CSVReaderBuilder((new BufferedReader(
									chunkPvalueFile.getName().endsWith(".gz")
									? new InputStreamReader(new GZIPInputStream(new FileInputStream(chunkPvalueFile)))
									: new FileReader(chunkPvalueFile)
							))).withSkipLines(0).withCSVParser(
									new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build()
							).build();

					String[] nextLine;
					genes:
					for (int g = 0; g < chunkGenes.size(); ++g) {

						if ((nextLine = reader.readNext()) == null) {
							throw new RuntimeException("Not enough lines in: " + chunkPvalueFile.getAbsolutePath());
						}
						if (nextLine.length != nrPermutations) {
							throw new RuntimeException("Not enough columns in: " + chunkPvalueFile.getAbsolutePath());
						}

						for (int p = 0; p < nrPermutations; ++p) {

							if (nextLine[p].equalsIgnoreCase("nan")) {
								continue genes;
							}

							chunkGenePvalues.setQuick(p, g, Double.parseDouble(nextLine[p]));
						}

						//if here there are no nan values
						String arm = geneArmMap.get(chunkGenes.get(g));
						nonNanGenesPerArm.get(arm).add(chunkGenes.get(g));

						nonNanGeneCount++;

					}

				}

				if (!chrGeneNames.isEmpty()) {
					throw new RuntimeException("Not all genes of chr " + chr + " are found in the permuation files");
				}

				LOGGER.info("Chr " + chr + " total genes: " + totalGenes + " genes without NaN values: " + nonNanGeneCount);

				for (Map.Entry<String, ArrayList<String>> armSet : nonNanGenesPerArm.entrySet()) {
					
					if(armSet.getValue().isEmpty()){
						LOGGER.info("No genes in: " + armSet.getKey());
					}
					
					DoubleMatrixDataset<String, String> chrGenePvaluesArm = chrGenePvalues.viewColSelection(armSet.getValue());

					if (options.isForceNormalGenePvalues()) {
						chrGenePvaluesArm.createColumnForceNormalInplace();
					}

					chrGenePvaluesArm.normalizeColumns();

					DoubleMatrixDataset<String, String> chrCorrelationMatrix = chrGenePvaluesArm.calculateCorrelationMatrixOnNormalizedColumns();

					DoubleMatrixDatasetRowCompressedWriter.saveDataset(
							options.getOutputBasePath() + "chr_" + armSet.getKey() + "_correlations",
							chrCorrelationMatrix,
							"Gene-Gene correlation matrix of chr " + armSet.getKey(), "Genes", "Genes");
					
					LOGGER.info("Done calculating correlation and writing restuls for " + armSet.getValue().size() + " genes of " + armSet.getKey());
				}

			} catch (Exception ex) {
				throw new RuntimeException(ex);
			}

		});

	}

	private static TObjectIntHashMap<String> readChunks(File chuckFile) throws IOException {
		final CSVReader reader
				= new CSVReaderBuilder((new BufferedReader(
						chuckFile.getName().endsWith(".gz")
						? new InputStreamReader(new GZIPInputStream(new FileInputStream(chuckFile)))
						: new FileReader(chuckFile)
				))).withSkipLines(0).withCSVParser(
						new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build()
				).build();

		TObjectIntHashMap<String> chunkInfo = new TObjectIntHashMap<>();

		String[] nextLine;
		while ((nextLine = reader.readNext()) != null) {
			chunkInfo.put(nextLine[0], Integer.parseInt(nextLine[1]));
		}
		reader.close();

		return chunkInfo;
	}

	private static List<String> readRowColNames(File file) throws IOException {
		final CSVReader reader
				= new CSVReaderBuilder((new BufferedReader(
						file.getName().endsWith(".gz")
						? new InputStreamReader(new GZIPInputStream(new FileInputStream(file)))
						: new FileReader(file)
				))).withSkipLines(0).withCSVParser(
						new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build()
				).build();

		ArrayList<String> names = new ArrayList<>();

		String[] nextLine;
		while ((nextLine = reader.readNext()) != null) {
			names.add(nextLine[0]);
		}
		reader.close();

		return Collections.unmodifiableList(names);
	}

}
