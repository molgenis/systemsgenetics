package eqtlmappingpipeline.ase;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import org.apache.commons.lang3.StringUtils;
import org.molgenis.genotype.snpeff.SnpEffEffect;
import umcg.genetica.collections.ChrPosMap;

/**
 *
 * @author Patrick Deelen
 */
public class AnnotateAseWithSnpEffVcf {

	public static void annotateAseWithSnpEffVcf(String aseFilePath, String snpEffVcfFilePath, String aseOutputPath) throws IOException, Exception {

		System.out.println("In file: " + aseFilePath);
		System.out.println("SnpEff VCF: " + snpEffVcfFilePath);
		System.out.println("Output: " + aseOutputPath);

		ChrPosMap<SnpEffEffect[]> snpEffAnnotations = SnpEffAnnotationMap.loadSnpEffAnnotationMap(snpEffVcfFilePath);
		
		System.out.println("Loading SnpEff VCF completed");

		final BufferedWriter outputWriter = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(aseOutputPath), AseConfiguration.ENCODING));
		final BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(aseFilePath), "UTF-8"));

		outputWriter.append(reader.readLine());
		outputWriter.append("\tSnpEff_EffectTypes\tSnpEff_EffectImpacts\tSnpEff_FunctionalClasses\tSnpEff_StrongestImpact\n");

		String line;
		while ((line = reader.readLine()) != null) {

			outputWriter.append(line);

			String[] aseLineElements = StringUtils.split(line, '\t');

			String chr = aseLineElements[5];
			int pos;
			try {
				pos = Integer.parseInt(aseLineElements[6]);
			} catch(NumberFormatException ex) {
				System.err.println("Error " + aseLineElements[6] + " is not an int for line: " + line );
				return;
			}

			SnpEffEffect[] snpEffAnnotation = snpEffAnnotations.get(chr, pos);

			if (snpEffAnnotation == null || snpEffAnnotation.length == 0) {
				outputWriter.append("\t\t\t\t\n");
			} else {
				StringBuffer effectTypeString = new StringBuffer("\t");
				StringBuffer effectImpactString = new StringBuffer("\t");
				StringBuffer effectFunctionalClassString = new StringBuffer("\t");

				SnpEffEffect.FunctionalClass strongestFunctionalClass = SnpEffEffect.FunctionalClass.NONE;

				boolean notFirst = false;
				for (SnpEffEffect snpEffEffect : snpEffAnnotation) {
					if (notFirst) {
						effectTypeString.append(';');
						effectImpactString.append(';');
						effectFunctionalClassString.append(';');
					}
					notFirst = true;

					effectTypeString.append(snpEffEffect.getEffectType());
					effectImpactString.append(snpEffEffect.getEffectImpact());
					effectFunctionalClassString.append(snpEffEffect.getFunctionalClass());

					if (strongestFunctionalClass.ordinal() < snpEffEffect.getFunctionalClass().ordinal()) {
						strongestFunctionalClass = snpEffEffect.getFunctionalClass();
					}

				}

				outputWriter.append(effectTypeString);
				outputWriter.append(effectImpactString);
				outputWriter.append(effectFunctionalClassString);
				outputWriter.append('\t');
				outputWriter.append(strongestFunctionalClass.toString());
				outputWriter.append('\n');

			}


		}

		outputWriter.close();
		reader.close();


	}
}
