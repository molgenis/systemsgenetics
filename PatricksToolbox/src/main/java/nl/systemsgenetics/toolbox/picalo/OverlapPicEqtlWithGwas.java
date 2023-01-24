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
import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;
import nl.systemsgenetics.toolbox.Options;
import nl.systemsgenetics.toolbox.Utils;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.LogManager;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.util.LdCalculator;
import org.molgenis.genotype.variant.GeneticVariant;
import umcg.genetica.io.gwascatalog.GWASCatalog;
import umcg.genetica.io.gwascatalog.GWASSNP;
import umcg.genetica.io.gwascatalog.GWASTrait;
import umcg.genetica.math.stats.FisherExactTest;

/**
 *
 * @author patri
 */
public class OverlapPicEqtlWithGwas {

	private static final Logger LOGGER = LogManager.getLogger(OverlapPicEqtlWithGwas.class);

	public static void overlap(Options options) throws IOException, Exception {

		GWASCatalog gwasCatalog = new GWASCatalog(options.getGwasCatalogFile(), 5e-8);
		HashMap<String, GWASSNP> gwasCatalogSnpMap = gwasCatalog.getSnpToObj();

		Set<String> gwasVariants = gwasCatalog.getAllSnpsIds();

		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		final CSVReader reader = new CSVReaderBuilder((new InputStreamReader(new GZIPInputStream(new FileInputStream(options.getInputPath()))))).withSkipLines(1).withCSVParser(parser).build();

		LinkedHashSet<String> interactionSnps = new LinkedHashSet<>();
		LinkedHashSet<String> nonInteractionSnps = new LinkedHashSet<>();
		LinkedHashSet<String> eqtlSnps = new LinkedHashSet<>();

		LOGGER.info("Assuming last col is the FDR");

		Pattern pattern = Pattern.compile("(rs\\d+)", Pattern.CASE_INSENSITIVE);

		String[] nextLine;
		while ((nextLine = reader.readNext()) != null) {

			String snpString = nextLine[0];

			Matcher matcher = pattern.matcher(snpString);
			boolean matchFound = matcher.find();
			if (matchFound) {
				String snp = matcher.group(1);
				eqtlSnps.add(snp);
				if (Double.parseDouble(nextLine[nextLine.length - 1]) <= 0.05) {
					interactionSnps.add(snp);
				} else {
					nonInteractionSnps.add(snp);
				}
			}

		}

		LOGGER.info("Loaded interaction SNPs: " + interactionSnps.size());
		LOGGER.info("Loaded non interaction SNPs: " + nonInteractionSnps.size());

		HashSet<String> combinedVariants = new HashSet<>();
		combinedVariants.addAll(gwasVariants);
		combinedVariants.addAll(eqtlSnps);

		RandomAccessGenotypeData genotypeReference = Utils.loadGenotypes(options, combinedVariants);

		HashMap<String, GeneticVariant> variantMap = genotypeReference.getVariantIdMap();

		//LOGGER.debug("VariantMap contains: " + variantMap.size());
		//count per trait how many variants are found in the genotype data
		TObjectIntMap<String> traitVariantCounts = new TObjectIntHashMap<>();
		HashMap<String, HashSet<String>> traitToVariantMapOverlappingWithPic = new HashMap<>();
		HashMap<String, HashSet<String>> traitToPicVariantMapOverlappingWithPic = new HashMap<>();
		HashMap<String, HashSet<String>> traitToVariantMapNotOverlappingWithPic = new HashMap<>();
		HashMap<String, HashSet<String>> traitToPicVariantMapNotOverlappingWithPic = new HashMap<>();
		for (GWASTrait trait : gwasCatalog.getTraitToObj().values()) {
			int traitVariantCount = 0;
			Set<GWASSNP> traitSnps = trait.getSnps();
			for (GWASSNP traitSnp : traitSnps) {

				String traitSnpId = traitSnp.getName();
				//LOGGER.debug(trait.getName() + " " + traitSnpId);
				if (variantMap.containsKey(traitSnpId)) {
					traitVariantCount++;
				}
			}
			traitToVariantMapOverlappingWithPic.put(trait.getName(), new HashSet<>());
			traitToPicVariantMapOverlappingWithPic.put(trait.getName(), new HashSet<>());
			traitToVariantMapNotOverlappingWithPic.put(trait.getName(), new HashSet<>());
			traitToPicVariantMapNotOverlappingWithPic.put(trait.getName(), new HashSet<>());
			traitVariantCounts.put(trait.getName(), traitVariantCount);

			//LOGGER.debug("Fase1: " + trait.getName() + " " + traitVariantCounts.get(trait.getName()));
		}

		CSVWriter writer = new CSVWriter(new FileWriter(options.getOutputBasePath() + "overlap.txt"), '\t', '\0', '\0', "\n");

		String[] outputLine = new String[2];
		int c = 0;
		outputLine[c++] = "interactionSNP";
		//outputLine[c++] = "LD";
		//outputLine[c++] = "gwasSNP";
		outputLine[c++] = "Traits";
		writer.writeNext(outputLine);

		int interactionsSnpsInGenotypeData = 0;
		int nonInteractionsSnpsInGenotypeData = 0;

		for (String snp : eqtlSnps) {

			if (variantMap.containsKey(snp)) {
				if (interactionSnps.contains(snp)) {
					interactionsSnpsInGenotypeData++;
				} else {
					nonInteractionsSnpsInGenotypeData++;
				}

				GeneticVariant eqtlVariant = variantMap.get(snp);

				HashSet<String> combinedTraits = new HashSet<>();

				Iterable<GeneticVariant> otherVariants = genotypeReference.getVariantsByRange(eqtlVariant.getSequenceName(), Integer.max(0, eqtlVariant.getStartPos() - 250000), eqtlVariant.getStartPos() + 250000);
				for (GeneticVariant otherVar : otherVariants) {
					String otherVariantId = otherVar.getVariantId().getPrimairyId();
					if (gwasVariants.contains(otherVariantId)) {
						double ld = LdCalculator.calculateLd(eqtlVariant, otherVar).getR2();
						if (ld >= 0.8) {

							HashSet<GWASTrait> traits = gwasCatalogSnpMap.get(otherVariantId).getAssociatedTraits();
							for (GWASTrait trait : traits) {
								combinedTraits.add(trait.getName());

								if (interactionSnps.contains(snp)) {
									traitToVariantMapOverlappingWithPic.get(trait.getName()).add(otherVariantId);
									traitToPicVariantMapOverlappingWithPic.get(trait.getName()).add(snp);
									if (trait.getName().equals("Asthma")) {
										LOGGER.info(snp + "\t" + otherVariantId + "\t" + ld);
									}

								} else {
									traitToVariantMapNotOverlappingWithPic.get(trait.getName()).add(otherVariantId);
									traitToPicVariantMapNotOverlappingWithPic.get(trait.getName()).add(snp);
								}

							}

						}
					}
				}

				outputLine[0] = snp;
				outputLine[1] = String.join(";", combinedTraits);
				writer.writeNext(outputLine);

			}

		}

		LOGGER.info(String.join(";", traitToPicVariantMapOverlappingWithPic.get("Asthma")));

		writer.close();

		CSVWriter enrichmentWriter = new CSVWriter(new FileWriter(options.getOutputBasePath() + "gwasEnrichment.txt"), '\t', '\0', '\0', "\n");

		outputLine = new String[7];
		c = 0;
		outputLine[c++] = "trait";
		outputLine[c++] = "traitVariants";
		outputLine[c++] = "nonPicVariantsOverlappingWithTrait";
		outputLine[c++] = "nonPicVariantsNotOverlappingWithTrait";
		outputLine[c++] = "picVariantsOverlappingWithTrait";
		outputLine[c++] = "picVariantsNotOverlappingWithTrait";
		outputLine[c++] = "FisherPvalue";
		enrichmentWriter.writeNext(outputLine);

		FisherExactTest ft = new FisherExactTest();

		for (GWASTrait trait : gwasCatalog.getTraitToObj().values()) {

			//int traitSnpsOverlappingWithPic = traitToVariantMapOverlappingWithPic.get(trait.getName()).size();
			int totalTraitSnps = traitVariantCounts.get(trait.getName());
			//int traitSnpsNotOverlappingWithPic = totalTraitSnps - traitSnpsOverlappingWithPic;

			int nonPicSnpsOverlapping = traitToPicVariantMapNotOverlappingWithPic.get(trait.getName()).size();
			int nonPicSnpsNotOverllaping = nonInteractionsSnpsInGenotypeData - nonPicSnpsOverlapping;

			int picSnpsOverlapping = traitToPicVariantMapOverlappingWithPic.get(trait.getName()).size();
			int picSnpsNotOverllaping = interactionsSnpsInGenotypeData - picSnpsOverlapping;

			double pvalue = ft.getFisherPValue(nonPicSnpsOverlapping, nonPicSnpsNotOverllaping, picSnpsOverlapping, picSnpsNotOverllaping);

			LOGGER.debug("Fase2: " + trait.getName() + " " + totalTraitSnps);

			c = 0;
			outputLine[c++] = trait.getName();
			outputLine[c++] = String.valueOf(totalTraitSnps);
			outputLine[c++] = String.valueOf(nonPicSnpsOverlapping);
			outputLine[c++] = String.valueOf(nonPicSnpsNotOverllaping);
			outputLine[c++] = String.valueOf(picSnpsOverlapping);
			outputLine[c++] = String.valueOf(picSnpsNotOverllaping);
			outputLine[c++] = String.valueOf(pvalue);
			enrichmentWriter.writeNext(outputLine);

		}

		enrichmentWriter.close();

	}

}
