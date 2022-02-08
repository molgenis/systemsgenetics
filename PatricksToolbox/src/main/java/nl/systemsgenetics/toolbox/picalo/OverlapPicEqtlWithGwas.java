/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.toolbox.picalo;

import com.opencsv.CSVParser;
import com.opencsv.CSVParserBuilder;
import com.opencsv.CSVReader;
import com.opencsv.CSVReaderBuilder;
import com.opencsv.CSVWriter;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.Set;
import java.util.StringJoiner;
import java.util.zip.GZIPInputStream;
import nl.systemsgenetics.toolbox.Options;
import nl.systemsgenetics.toolbox.Utils;
import org.apache.log4j.Logger;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.util.LdCalculator;
import org.molgenis.genotype.variant.GeneticVariant;
import umcg.genetica.io.gwascatalog.GWASCatalog;
import umcg.genetica.io.gwascatalog.GWASSNP;
import umcg.genetica.io.gwascatalog.GWASTrait;

/**
 *
 * @author patri
 */
public class OverlapPicEqtlWithGwas {

	private static final Logger LOGGER = Logger.getLogger(OverlapPicEqtlWithGwas.class);

	public static void overlap(Options options) throws IOException, Exception {

		GWASCatalog gwasCatalog = new GWASCatalog(options.getGwasCatalogFile(), 5e-8);
		HashMap<String, GWASSNP> gwasCatalogSnpMap = gwasCatalog.getSnpToObj();

		Set<String> gwasVariants = gwasCatalog.getAllSnpsIds();

		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		final CSVReader reader = new CSVReaderBuilder((new InputStreamReader(new GZIPInputStream(new FileInputStream(options.getInputPath()))))).withSkipLines(1).withCSVParser(parser).build();

		LinkedHashSet<String> interactionSnps = new LinkedHashSet<>();

		LOGGER.info("Assuming last col is the FDR");
		
		String[] nextLine;
		while ((nextLine = reader.readNext()) != null) {

			if(Double.parseDouble(nextLine[nextLine.length - 1]) <= 0.05){
				interactionSnps.add(nextLine[0]);
			}
			
			

		}

		LOGGER.info("Loaded interaction SNPs: " + interactionSnps.size());

		HashSet<String> combinedVariants = new HashSet<>();
		combinedVariants.addAll(gwasVariants);
		combinedVariants.addAll(interactionSnps);

		RandomAccessGenotypeData genotypeReference = Utils.loadGenotypes(options, combinedVariants);

		HashMap<String, GeneticVariant> variantMap = genotypeReference.getVariantIdMap();

		CSVWriter writer = new CSVWriter(new FileWriter(options.getOutputBasePath() + "overlap.txt"), '\t', '\0', '\0', "\n");

		String[] outputLine = new String[4	];
		int c = 0;
		outputLine[c++] = "interactionSNP";
		outputLine[c++] = "LD";
		outputLine[c++] = "gwasSNP";
		outputLine[c++] = "Traits";
		writer.writeNext(outputLine);

		for (String snp : interactionSnps) {

			if (variantMap.containsKey(snp)) {
				GeneticVariant interactionVariant = variantMap.get(snp);

				Iterable<GeneticVariant> otherVariants = genotypeReference.getVariantsByRange(interactionVariant.getSequenceName(), Integer.max(0, interactionVariant.getStartPos() - 100000), interactionVariant.getStartPos() + 100000);
				for (GeneticVariant otherVar : otherVariants) {
					String otherVariantId = otherVar.getVariantId().getPrimairyId();
					if (gwasVariants.contains(otherVariantId)) {
						double ld = LdCalculator.calculateLd(interactionVariant, otherVar).getR2();
						if (ld >= 0.8) {

							HashSet<GWASTrait> traits = gwasCatalogSnpMap.get(otherVariantId).getAssociatedTraits();
							StringJoiner joiner = new StringJoiner(";");
							for (GWASTrait trait : traits) {
								joiner.add(trait.getName());
							}
							String traitsString = joiner.toString();

							outputLine[0] = snp;
							outputLine[1] = String.valueOf(ld);
							outputLine[2] = otherVariantId;
							outputLine[3] = traitsString;
							writer.writeNext(outputLine);

						}
					}
				}

			}

		}

		writer.close();

	}

}
