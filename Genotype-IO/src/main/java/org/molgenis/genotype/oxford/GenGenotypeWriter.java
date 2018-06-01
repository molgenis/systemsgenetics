/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.oxford;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.nio.charset.Charset;
import java.util.HashMap;
import org.apache.log4j.Logger;
import org.molgenis.genotype.Allele;
import org.molgenis.genotype.GenotypeData;
import org.molgenis.genotype.GenotypeWriter;
import org.molgenis.genotype.Sample;
import org.molgenis.genotype.util.Utils;
import org.molgenis.genotype.variant.GeneticVariant;

/**
 *
 * @author Patrick Deelen
 */
public class GenGenotypeWriter implements GenotypeWriter {

	public static final Charset FILE_ENCODING = Charset.forName("UTF-8");
	public static final char LINE_ENDING = '\n';
	private static final char SEPARATOR = ' ';
	private static final Logger LOGGER = Logger.getLogger(GenGenotypeWriter.class);
	private final GenotypeData genotypeData;

	public GenGenotypeWriter(GenotypeData genotypeData) {
		this.genotypeData = genotypeData;
	}

	@Override
	public void write(String basePath) throws IOException {

		if (!genotypeData.isOnlyContaingSaveProbabilityGenotypes()) {
			LOGGER.warn("WARNING!!! writing dosage genotype data to .gen posterior probabilities file. Using heuristic method to convert to probabilities, this is not guaranteed to be accurate. See manual for more details.");
		}

		write(new File(basePath + ".gen"), new File(basePath + ".sample"));
	}

	public void write(File genFile, File sampleFile) throws IOException {
		LOGGER.info("Writing gen file " + genFile.getAbsolutePath() + " and sample file "
				+ sampleFile.getAbsolutePath());

		Utils.createEmptyFile(genFile, "gen");
		Utils.createEmptyFile(sampleFile, "sample");

		HashMap<Sample, Float> sampleMissingness = writeGenFile(genFile);
		OxfordSampleFileWriter.writeSampleFile(sampleFile, genotypeData, sampleMissingness);

	}

	private HashMap<Sample, Float> writeGenFile(File hapsFile) throws IOException {

		BufferedWriter hapsFileWriter = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(hapsFile), FILE_ENCODING));

		float[] sampleMissingCount = new float[genotypeData.getSamples().size()];
		int totalVariants = 0;

		for (GeneticVariant variant : genotypeData) {

			++totalVariants;

			if (variant.getAlleleCount() > 2) {
				LOGGER.warn("Skipping variant: " + variant.getPrimaryVariantId() + " at " + variant.getSequenceName() + ":" + variant.getStartPos() + " with more than 2 alleles: " + variant.getVariantAlleles());
			}

			Allele allele0 = variant.getVariantAlleles().get(0);
			Allele allele1 = variant.getAlleleCount() == 1 ? Allele.ZERO : variant.getVariantAlleles().get(1);

			hapsFileWriter.append(variant.getSequenceName());
			hapsFileWriter.append(SEPARATOR);
			hapsFileWriter.append(variant.getPrimaryVariantId());
			hapsFileWriter.append(SEPARATOR);
			hapsFileWriter.append(String.valueOf(variant.getStartPos()));
			hapsFileWriter.append(SEPARATOR);
			hapsFileWriter.append(allele0.getAlleleAsString());
			hapsFileWriter.append(SEPARATOR);
			hapsFileWriter.append(allele1.getAlleleAsString());

			float[][] probs = variant.getSampleGenotypeProbilities();

			for (int i = 0; i < probs.length; i++) {
				boolean isMissing = true;
				for (float prob : probs[i]) {
					if (prob > 0) {
						isMissing = false;
					}
					hapsFileWriter.append(SEPARATOR);
					if(prob == 0){
						hapsFileWriter.append('0');
					} else if (prob == 1){
						hapsFileWriter.append('1');
					} else {
						hapsFileWriter.append(String.valueOf(prob));
					}
					
				}
				if (isMissing) {
					sampleMissingCount[i]++;
				}
			}

			hapsFileWriter.append(LINE_ENDING);

		}

		hapsFileWriter.close();

		HashMap<Sample, Float> sampleMissingness = new HashMap<Sample, Float>();
		for (int i = 0; i < sampleMissingCount.length; ++i) {
			sampleMissingness.put(genotypeData.getSamples().get(i), sampleMissingCount[i] / (float) totalVariants);
		}
		return sampleMissingness;

	}
}
