/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.downstreamer.runners;

import com.opencsv.CSVWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.stream.IntStream;
import nl.systemsgenetics.downstreamer.DownstreamerOptions;
import nl.systemsgenetics.downstreamer.io.IoUtils;
import org.apache.log4j.Logger;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.util.LdCalculatorException;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variantFilter.VariantIdIncludeFilter;
import umcg.genetica.math.matrix2.DoubleMatrixDatasetRowIterable;
import umcg.genetica.math.stats.ZScores;

/**
 *
 * @author patri
 */
public class PruneToIndependentTopHits {

	private static final Logger LOGGER = Logger.getLogger(PruneToIndependentTopHits.class);

	public static void prune(final DownstreamerOptions options) throws IOException{
		prune(options, null);
	}
	
	public static void prune(final DownstreamerOptions options, final RandomAccessGenotypeData referenceGenotypeData) throws IOException {

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

					//Significant variant
					String variant = allGwasVariants.get(variantI);
					significantVariantsPerTrait.get(p).add(new VariantZscore(variant, variantZscores[p]));
					allSignificantVariants.add(variant);

				}

			}

		}
		

		final HashMap<String, GeneticVariant> genotypeIdVariantMap;
		
		if(referenceGenotypeData == null){
			RandomAccessGenotypeData referenceGenotypeData2 = IoUtils.loadGenotypes(options, allSignificantVariants);
			genotypeIdVariantMap = referenceGenotypeData2.getVariantIdMap();//only significant variants are loaded so make map of all
		} else {
			//only make map of significant variants, we don't need any others.
			genotypeIdVariantMap = referenceGenotypeData.getVariantIdMap(new VariantIdIncludeFilter(allSignificantVariants));
		}
		
		IntStream.range(0, numberPheno).parallel().forEach(p -> {

			TreeSet<VariantZscore> signficantVariants = significantVariantsPerTrait.get(p);

			final ArrayList<VariantZscore> topVariants = topVariantsPerTrait.get(p);

			VariantZscore topVariant;
			while ((topVariant = signficantVariants.pollFirst()) != null) {
				
				final GeneticVariant topVariantGeno = genotypeIdVariantMap.get(topVariant.getVariant());
				
				if(topVariantGeno == null || !topVariantGeno.isBiallelic() ){
					continue;
				}

				LOGGER.debug("Top variant: " + topVariant.getVariant() + " - " + topVariant.getNegZscore());

				topVariants.add(topVariant);

				Iterator<VariantZscore> otherVariantIterator = signficantVariants.iterator();

				while (otherVariantIterator.hasNext()) {
					VariantZscore otherVariant = otherVariantIterator.next();
					GeneticVariant otherVariantGeno = genotypeIdVariantMap.get(otherVariant.getVariant());

					if(otherVariantGeno == null || !otherVariantGeno.isBiallelic()){
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
							LOGGER.fatal(topVariantGeno.getVariantId() + "\t" + topVariantGeno.getMinorAlleleFrequency()+ "\t" + topVariantGeno.getVariantAlleles().toString());
							LOGGER.fatal(otherVariantGeno.getVariantId() + "\t" + otherVariantGeno.getMinorAlleleFrequency()+ "\t" + otherVariantGeno.getVariantAlleles().toString());
							throw new RuntimeException(ex);
						}
					}

				}

				//LOGGER.info("Significant variants left: " + signficantVariants.size());

			}

		});
		

		CSVWriter writer = new CSVWriter(new FileWriter(options.getGwasTopHitsFile()), '\t', '\0', '\0', "\n");

		for (int p = 0; p < numberPheno; ++p) {

			LOGGER.info("Identified: " + topVariantsPerTrait.get(p).size() + " independent loci for: " + phenotypeNames.get(p));
			
			final ArrayList<VariantZscore> topVariants = topVariantsPerTrait.get(p);
			
			String[] outputLine = new String[1 + topVariants.size()];
			int c = 0;
			outputLine[c++] = phenotypeNames.get(p);
			for (VariantZscore topVariant : topVariants) {
				outputLine[c++] = topVariant.getVariant();
			}
			writer.writeNext(outputLine);

		}

		writer.close();

	}

	private static class VariantZscore implements Comparable<VariantZscore> {

		final String variant;
		final double negZscore;

		public VariantZscore(String variant, double zscore) {
			this.variant = variant;
			this.negZscore = zscore <= 0 ? zscore : (zscore * -1);
		}

		public String getVariant() {
			return variant;
		}

		public double getNegZscore() {
			return negZscore;
		}

		@Override
		public int hashCode() {
			return variant.hashCode();
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
			if (!Objects.equals(this.variant, other.variant)) {
				return false;
			}
			return true;
		}

		@Override
		public int compareTo(VariantZscore o) {
			int x = Double.compare(this.negZscore, o.negZscore);
			if(x == 0){
				//for tree map equal compare must be equals is true
				return this.variant.compareTo(o.variant);
			} else {
				return x;
			}
		}


		
		
	}

}
