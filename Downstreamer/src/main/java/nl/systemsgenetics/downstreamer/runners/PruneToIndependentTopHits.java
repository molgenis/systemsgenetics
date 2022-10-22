/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.downstreamer.runners;

import com.opencsv.CSVWriter;

import java.io.*;
import java.util.*;
import java.util.stream.IntStream;
import java.util.zip.GZIPOutputStream;

import nl.systemsgenetics.downstreamer.DownstreamerOptions;
import nl.systemsgenetics.downstreamer.io.IoUtils;
import nl.systemsgenetics.downstreamer.summarystatistic.Locus;
import nl.systemsgenetics.downstreamer.summarystatistic.LocusUtils;
import nl.systemsgenetics.downstreamer.summarystatistic.SummaryStatisticRecord;
import org.apache.log4j.Logger;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.util.LdCalculatorException;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variantFilter.VariantIdIncludeFilter;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.matrix2.DoubleMatrixDatasetRowIterable;
import umcg.genetica.math.stats.ZScores;

/**
 *
 * @author patri
 */
public class PruneToIndependentTopHits {

	private static final Logger LOGGER = Logger.getLogger(PruneToIndependentTopHits.class);

	public static void prune(final DownstreamerOptions options) throws Exception {
		prune(options, null);
	}

	public static void prune(final DownstreamerOptions options, RandomAccessGenotypeData referenceGenotypeData) throws Exception {

		options.getLeadSnpsFolder().mkdirs();

		// GWAS Z-scores
		final String variantPhenotypeZscoreMatrixPath = options.getGwasZscoreMatrixPath();
		final DoubleMatrixDatasetRowIterable gwasStreamer = new DoubleMatrixDatasetRowIterable(variantPhenotypeZscoreMatrixPath);

		// Column names of GWAS Z-scores
		final ArrayList<String> phenotypeNames = new ArrayList<>(gwasStreamer.getCols());
		final int numberPheno = gwasStreamer.getNrCols();

		// Variant identifiers
		final ArrayList<String> allGwasVariants = new ArrayList<>(gwasStreamer.getRows());

		// Clumping parameters
		final double variantLdThreshold = 0.1;
		final int window=1000000;
		final double pvalueThreshold = 5e-8;

		// Determine p-value threshold as a z-score
		final double zscoreThresholdNeg = ZScores.pToZTwoTailed(pvalueThreshold);
		final double zscoreThresholdPos = zscoreThresholdNeg * -1;

		// Init storage
		final ArrayList<Map<String, SummaryStatisticRecord>> significantVariantsPerTrait = new ArrayList<>(numberPheno);
		for (int p = 0; p < numberPheno; ++p) {
			significantVariantsPerTrait.add(new HashMap<>());
		}

		// Determine significant variants
		final HashSet<String> allSignificantVariants = new HashSet<>();
		int variantI = -1; //keeps track of current row in gwasStreamer
		for (double[] variantZscores : gwasStreamer) {
			//variantZscores contains a single row from the gwas zscore matrix, one element per studied trait
			variantI++;
			for (int p = 0; p < numberPheno; ++p) {
				if (variantZscores[p] <= zscoreThresholdNeg || variantZscores[p] >= zscoreThresholdPos) {
					// Significant variantId
					allSignificantVariants.add(allGwasVariants.get(variantI));
				}
			}
		}

		LOGGER.info("Found " + allSignificantVariants.size() + " variants at p<" + pvalueThreshold);
		LOGGER.info("Loading genotypes");

		// Load genotype data
		Map<String, GeneticVariant> variantMap;

		// Only significant variants are loaded so make map of all
		if (referenceGenotypeData == null) {
			referenceGenotypeData = IoUtils.loadGenotypes(options);
		}
		LOGGER.info("Debug");

		variantMap = referenceGenotypeData.getVariantIdMap(new VariantIdIncludeFilter(allSignificantVariants));

		LOGGER.info("Done loading genotypes");

		// Construct SummaryStatisticsRecords for clumping
		DoubleMatrixDataset<String, String> gwasZscores = DoubleMatrixDataset.loadSubsetOfRowsBinaryDoubleData(variantPhenotypeZscoreMatrixPath, allSignificantVariants);

		for (String variant : allSignificantVariants) {
			GeneticVariant curVariant = variantMap.get(variant);
			for (int p = 0; p < numberPheno; ++p) {
				if (curVariant != null && curVariant.isBiallelic()) {
					SummaryStatisticRecord curRecord = new SummaryStatisticRecord(curVariant, ZScores.zToP(gwasZscores.getElement(variant, phenotypeNames.get(p))));
					significantVariantsPerTrait.get(p).put(variant, curRecord);
				}
			}
		}

		// Args for clumping: Upper p-value threshold, lower p-value threshold, window, LD threshold
		// TODO: Update this uglyness in selectIndepTopHits
		String[] args = new String[] {Double.toString(pvalueThreshold), Double.toString(pvalueThreshold), Integer.toString(window), Double.toString(variantLdThreshold)};

		LOGGER.info("Starting clumping");

		// Loop over phenotypes
		for (int p = 0; p < numberPheno; ++p) {

			CSVWriter writer = new CSVWriter(new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(new File(options.getLeadSnpsFolder(),phenotypeNames.get(p) + "_leads.txt.gz"))))), '\t', '\0', '\0', "\n");

			// Clump
			List<Locus> curLoci = LocusUtils.selectIndepTopHits(args, referenceGenotypeData, significantVariantsPerTrait.get(p));

			int nIndepVariants = 0;
			String[] outputLine = new String[4];
			int c = 0;
			outputLine[c++] = "Variant";
			outputLine[c++] = "Chr";
			outputLine[c++] = "Pos";
			outputLine[c++] = "Pvalue";
			writer.writeNext(outputLine);

			for (Locus curLocus : curLoci) {
				for (SummaryStatisticRecord topVariant : curLocus.getIndepSummaryStatisticRecords()) {
					c = 0;
					outputLine[c++] = topVariant.getPrimaryVariantId();
					outputLine[c++] = topVariant.getContig();
					outputLine[c++] = String.valueOf(topVariant.getStart());
					outputLine[c++] = String.valueOf(topVariant.getPvalue());
					writer.writeNext(outputLine);
					nIndepVariants++;
				}
			}
			writer.close();

			LOGGER.info("Identified: " + curLoci.size() + " independent loci with : " + nIndepVariants + " variants for: " + phenotypeNames.get(p));
		}

		LOGGER.info("Done");
	}

	@Deprecated
	public static void pruneOld(final DownstreamerOptions options) throws IOException {
		pruneOld(options, null);
	}

	@Deprecated
	public static void pruneOld(final DownstreamerOptions options, final RandomAccessGenotypeData referenceGenotypeData) throws IOException {

		options.getLeadSnpsFolder().mkdirs();

		final String variantPhenotypeZscoreMatrixPath = options.getGwasZscoreMatrixPath();
		final DoubleMatrixDatasetRowIterable gwasStreamer = new DoubleMatrixDatasetRowIterable(variantPhenotypeZscoreMatrixPath);
		final ArrayList<String> phenotypeNames = new ArrayList<>(gwasStreamer.getCols());
		final int numberPheno = gwasStreamer.getNrCols();
		final ArrayList<String> allGwasVariants = new ArrayList<>(gwasStreamer.getRows());
		final double pvalueThreshold = 5e-8;
		final double variantLdThreshold = 0.1;

		final double zscoreThresholdNeg = ZScores.pToZTwoTailed(pvalueThreshold);
		final double zscoreThresholdPos = zscoreThresholdNeg * -1;

		final ArrayList<TreeSet<VariantZscore>> significantVariantsPerTrait = new ArrayList<>(numberPheno);
		final ArrayList<ArrayList<VariantZscore>> topVariantsPerTrait = new ArrayList<>(numberPheno);
		for (int p = 0; p < numberPheno; ++p) {
			significantVariantsPerTrait.add(new TreeSet<>());
			topVariantsPerTrait.add(new ArrayList<>());
		}

		final HashSet<String> allSignificantVariants = new HashSet<>();

		int variantI = -1; //keeps track of current row in gwasStreamer
		for (double[] variantZscores : gwasStreamer) {
			//variantZscores contains a single row from the gwas zscore matrix, one element per studied trait
			variantI++;

			for (int p = 0; p < numberPheno; ++p) {

				if (variantZscores[p] <= zscoreThresholdNeg || variantZscores[p] >= zscoreThresholdPos) {

					//Significant variantId
					String variant = allGwasVariants.get(variantI);
					significantVariantsPerTrait.get(p).add(new VariantZscore(variant, variantZscores[p]));
					allSignificantVariants.add(variant);

				}

			}

		}

		final HashMap<String, GeneticVariant> genotypeIdVariantMap;

		if (referenceGenotypeData == null) {
			RandomAccessGenotypeData referenceGenotypeData2 = IoUtils.loadGenotypes(options, allSignificantVariants);
			genotypeIdVariantMap = referenceGenotypeData2.getVariantIdMap();//only significant variants are loaded so make map of all
		} else {
			//only make map of significant variants, we don't need any others.
			genotypeIdVariantMap = referenceGenotypeData.getVariantIdMap(new VariantIdIncludeFilter(allSignificantVariants));
		}

		IntStream.range(0, numberPheno).parallel().forEach(p -> {

			TreeSet<VariantZscore> signficantVariants = significantVariantsPerTrait.get(p);

			LOGGER.debug("Significant variants: " + signficantVariants.size());
			
			final ArrayList<VariantZscore> topVariants = topVariantsPerTrait.get(p);

			//While loop allows save delete within the loop
			VariantZscore topVariant;
			while ((topVariant = signficantVariants.pollFirst()) != null) {
			
				final GeneticVariant topVariantGeno = genotypeIdVariantMap.get(topVariant.getVariantId());

				if (topVariantGeno == null || !topVariantGeno.isBiallelic()) {
					continue;
				}

				LOGGER.debug("Top variant: " + topVariant.getVariantId()+ " - " + topVariant.getNegZscore());

				topVariant.setVariant(topVariantGeno);
				topVariants.add(topVariant);

				Iterator<VariantZscore> otherVariantIterator = signficantVariants.iterator();

				while (otherVariantIterator.hasNext()) {
					VariantZscore otherVariant = otherVariantIterator.next();
					GeneticVariant otherVariantGeno = genotypeIdVariantMap.get(otherVariant.getVariantId());

					if (otherVariantGeno == null || !otherVariantGeno.isBiallelic()) {
						otherVariantIterator.remove();
						continue;
					}

					if (topVariantGeno.getSequenceName().equals(otherVariantGeno.getSequenceName())) {
						try {
							if (topVariantGeno.calculateLd(otherVariantGeno).getR2() >= variantLdThreshold) {
								//LOGGER.info(" - Variant in LD: " + otherVariant.getVariant());
								otherVariantIterator.remove();
							}
						} catch (LdCalculatorException | java.lang.UnsupportedOperationException ex) {
							LOGGER.fatal(topVariantGeno.getVariantId() + "\t" + topVariantGeno.getMinorAlleleFrequency() + "\t" + topVariantGeno.getVariantAlleles().toString());
							LOGGER.fatal(otherVariantGeno.getVariantId() + "\t" + otherVariantGeno.getMinorAlleleFrequency() + "\t" + otherVariantGeno.getVariantAlleles().toString());
							throw new RuntimeException(ex);
						}
					}

				}

				//LOGGER.info("Significant variants left: " + signficantVariants.size());
			}

		});

		for (int p = 0; p < numberPheno; ++p) {

			CSVWriter writer = new CSVWriter(new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(new File(options.getLeadSnpsFolder(), phenotypeNames.get(p) + "_leads.txt.gz"))))), '\t', '\0', '\0', "\n");

			final ArrayList<VariantZscore> topVariants = topVariantsPerTrait.get(p);

			LOGGER.info("Identified: " + topVariants.size() + " independent loci for: " + phenotypeNames.get(p));

			String[] outputLine = new String[4];
			int c = 0;
			outputLine[c++] = "Variant";
			outputLine[c++] = "Chr";
			outputLine[c++] = "Pos";
			outputLine[c++] = "Pvalue";
			writer.writeNext(outputLine);

			for (VariantZscore topVariant : topVariants) {
				c = 0;
				outputLine[c++] = topVariant.getVariantId();
				outputLine[c++] = topVariant.getVariant().getSequenceName();
				outputLine[c++] = String.valueOf(topVariant.getVariant().getStartPos());
				outputLine[c++] = String.valueOf(ZScores.zToP(topVariant.getNegZscore()));
				writer.writeNext(outputLine);
			}

			writer.close();

		}

	}

	@Deprecated
	private static class VariantZscore implements Comparable<VariantZscore> {

		final String variantId;
		final double negZscore;
		GeneticVariant variant;

		public VariantZscore(String variantId, double zscore) {
			this.variantId = variantId;
			this.negZscore = zscore <= 0 ? zscore : (zscore * -1);
		}

		public String getVariantId() {
			return variantId;
		}

		public double getNegZscore() {
			return negZscore;
		}

		public GeneticVariant getVariant() {
			return variant;
		}

		public void setVariant(GeneticVariant variant) {
			this.variant = variant;
		}

		@Override
		public int hashCode() {
			return variantId.hashCode();
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
			final VariantZscore other = (VariantZscore) obj;
			if (Double.doubleToLongBits(this.negZscore) != Double.doubleToLongBits(other.negZscore)) {
				return false;
			}
			if (!Objects.equals(this.variantId, other.variantId)) {
				return false;
			}
			return true;
		}

		@Override
		public int compareTo(VariantZscore o) {
			int x = Double.compare(this.negZscore, o.negZscore);
			if (x == 0) {
				//for tree map equal compare must be equals is true
				return this.variantId.compareTo(o.variantId);
			} else {
				return x;
			}
		}

	}

}
