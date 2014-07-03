package eqtlmappingpipeline.util;

import com.google.common.collect.Lists;
import eqtlmappingpipeline.ase.AseConfiguration;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import org.apache.commons.lang3.StringUtils;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.RandomAccessGenotypeDataReaderFormats;
import org.molgenis.genotype.util.Ld;
import org.molgenis.genotype.variant.GeneticVariant;

/**
 *
 * @author Patrick Deelen
 */
public class NoLdSnpProbeListCreator {

	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) throws UnsupportedEncodingException, FileNotFoundException, IOException, Exception {

		File probeFile = null;
		String genotypePath = "";
		String genotypeType = "VCF_FOLDER";
		int windowHalfSize = 250000;
		int probeMargin = 2;
		double maxDprime = 0.2;
		double maxR2 = 0.2;
		File outputFile = null;

		RandomAccessGenotypeData genotypeData = RandomAccessGenotypeDataReaderFormats.valueOf(genotypeType).createGenotypeData(genotypePath, 10000);

		BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(probeFile), "UTF-8"));
		final BufferedWriter snpProbeToTestWriter = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outputFile), AseConfiguration.ENCODING));


		String line;
		String[] elements;
		while ((line = reader.readLine()) != null) {
			elements = StringUtils.splitPreserveAllTokens(line, '\t');

			if (elements.length < 4) {
				throw new Exception();
			}

			String chr = elements[0];//TODO remove chr
			int probeStartPos = Integer.parseInt(elements[1]);
			int probeStopPos = Integer.parseInt(elements[2]);
			String probeName = elements[3];


			int windowsStart = probeStartPos - windowHalfSize;
			int windowStop = probeStopPos + windowHalfSize;


			ArrayList<GeneticVariant> probeVariants = Lists.newArrayList(genotypeData.getVariantsByRange(chr, (probeStartPos - probeMargin), (probeStopPos + probeMargin)));

			variants:
			for (GeneticVariant variant : genotypeData.getVariantsByRange(chr, windowsStart, windowStop)) {

				//skip over probe variants
				if (variant.getStartPos() >= (probeStartPos - probeMargin) && variant.getStartPos() <= (probeStopPos + probeMargin)) {
					continue variants;
				}

				for (GeneticVariant probeVariant : probeVariants) {

					Ld ld = variant.calculateLd(probeVariant);

					if (ld.getDPrime() >= maxDprime || ld.getR2() >= maxR2) {
						//Exclude
						continue variants;
					}

				}

				//Probe SNP combination is okay
				snpProbeToTestWriter.append(variant.getPrimaryVariantId());
				snpProbeToTestWriter.append('\t');
				snpProbeToTestWriter.append(probeName);
				snpProbeToTestWriter.append('\n');

			}

			snpProbeToTestWriter.close();

		}


	}
}
