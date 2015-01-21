package eqtlmappingpipeline.ase;

import au.com.bytecode.opencsv.CSVWriter;
import cern.colt.list.tint.IntArrayList;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.regex.Pattern;

/**
 *
 * @author Patrick Deelen
 */
public class AseRecalculate {

	private static final Pattern TAB_PATTERN = Pattern.compile("\\t");

	public static void main(String[] args) throws Exception {

		File oldResultsFile = new File(args[0]);
		File newResultsFile = new File(args[1]);
		int threads = Integer.parseInt(args[2]);

		System.out.println("Old: " + oldResultsFile.getAbsolutePath());
		System.out.println("New: " + newResultsFile.getAbsolutePath());
		System.out.println("Threads: " + threads);

		BufferedReader aseReader = new BufferedReader(new FileReader(oldResultsFile));

		ArrayList<AseVariantBean> inputAse = new ArrayList<AseVariantBean>();

		String line;
		String[] elements;
		//Header
		aseReader.readLine();
		while ((line = aseReader.readLine()) != null) {
			elements = TAB_PATTERN.split(line);
			inputAse.add(new AseVariantBean(elements));
		}
		aseReader.close();

		AseVariantRecalculate[] aseVariants = new AseVariantRecalculate[inputAse.size()];
		{
			int i = 0;
			for (AseVariantBean aseVariant : inputAse) {
				aseVariants[i] = new AseVariantRecalculate(aseVariant);
				++i;
			}
		}

		AseCalculator.startAseCalculators(aseVariants, threads);

		System.out.println("Completed ASE calculations");

		CSVWriter mappingReportWriter = new CSVWriter(new FileWriter(newResultsFile), '\t', CSVWriter.NO_QUOTE_CHARACTER);

		final String[] newResults = new String[17];
		int c = 0;
		newResults[c++] = "OldPvalue";
		newResults[c++] = "OldRatioD";
		newResults[c++] = "OldP";
		newResults[c++] = "NewPvalue";
		newResults[c++] = "NewRatioD";
		newResults[c++] = "NewP";
		newResults[c++] = "NewTheta";
		newResults[c++] = "SampleCount";
		newResults[c++] = "Chr";
		newResults[c++] = "Pos";
		newResults[c++] = "Id";
		newResults[c++] = "Ref_Allele";
		newResults[c++] = "Alt_Allele";
		newResults[c++] = "Id";
		newResults[c++] = "Genes";
		newResults[c++] = "Ref_Counts";
		newResults[c++] = "Alt_Counts";
		mappingReportWriter.writeNext(newResults);

		for (AseVariantRecalculate ase : aseVariants) {

			c = 0;
			newResults[c++] = String.valueOf(ase.getOriginalLikelihoodRatioP());
			newResults[c++] = String.valueOf(ase.getOriginalLikelihoodRatioD());
			newResults[c++] = String.valueOf(ase.getOriginalEffect());
			newResults[c++] = String.valueOf(ase.getLikelihoodRatioP());
			newResults[c++] = String.valueOf(ase.getLikelihoodRatioD());
			newResults[c++] = String.valueOf(ase.getEffect());
			newResults[c++] = String.valueOf(ase.getMle().getMaxLogLikelihoodTheta());
			newResults[c++] = String.valueOf(ase.getSampleCount());
			newResults[c++] = ase.getChr();
			newResults[c++] = String.valueOf(ase.getPos());
			newResults[c++] = ase.getId().getPrimairyId();
			newResults[c++] = ase.getA1().getAlleleAsString();
			newResults[c++] = ase.getA2().getAlleleAsString();
			newResults[c++] = String.valueOf(ase.getPos());
			newResults[c++] = ase.getGenes();
			newResults[c++] = createCountString(ase.getA1Counts());
			newResults[c++] = createCountString(ase.getA2Counts());
			mappingReportWriter.writeNext(newResults);


		}

		mappingReportWriter.close();


	}

	public static String createCountString(IntArrayList counts) {

		StringBuilder builder = new StringBuilder();

		for (int i = 0; i < counts.size(); ++i) {
			if (i > 0) {
				builder.append(',');
			}
			builder.append(String.valueOf(counts.getQuick(i)));
		}
		
		return builder.toString();
	}
}
