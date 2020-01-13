package umcg.genetica.io.trityper;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import org.apache.commons.cli.*;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.bin.BinaryFile;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.util.Primitives;
import umcg.genetica.util.RankArray;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashSet;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Pattern;

/**
 * @author Lude, Marc Jan
 */
public class ConvertDoubleMatrixDataToTriTyper {


	private static final Pattern SPLIT_ON_TAB = Pattern.compile("\t");

	public static void main(String[] args) throws IOException, Exception {

		CommandLineParser parser = new GnuParser();
		Options options = new Options();

		Option datMatrix = OptionBuilder.withArgName("path").hasArg().withDescription("Location of the input file. Needs to be a tab seperated file with samples on the columns and traits on the rows.").withLongOpt("dataMatrix").create("d");
		Option mapFile = OptionBuilder.withArgName("path").hasArg().withDescription("Location of the mapping file describing the chromosomal locations of the traits.").withLongOpt("mappingFile").create("m");
		Option folderOut = OptionBuilder.withArgName("path").hasArg().withDescription("Location and name of the output TriTyper folder.").withLongOpt("OutputFile").create("o");
		Option ranking = OptionBuilder.withArgName("boolean").withDescription("If set first rank the input data, before scaling.").create("r");
		Option removeNan = OptionBuilder.withArgName("boolean").withDescription("If set first remove full non-numeric rows.").create("R");
		options.addOption(folderOut).addOption(datMatrix).addOption(mapFile).addOption(ranking).addOption(removeNan);

		String dataMatrix = null;
		String outputFolder = null;
		String mappingFile = null;
		boolean rank = false;
		boolean removeNanRow = false;
		CommandLine cmd;
		try {
			cmd = parser.parse(options, args);
			HelpFormatter formatter = new HelpFormatter();

			if (cmd.hasOption("OutputFile") || cmd.hasOption("o")) {
				// initialise the member variable
				outputFolder = cmd.getOptionValue("OutputFile");
			} else {
				System.out.println("Missing necesarray information: output loc");
				formatter.printHelp("ant", options);
				System.exit(0);
			}
			if (cmd.hasOption("dataMatrix") || cmd.hasOption("d")) {
				// initialise the member variable
				dataMatrix = cmd.getOptionValue("dataMatrix");
			} else {
				System.out.println("Missing necessary information: data matrix");
				formatter.printHelp("ant", options);
				System.exit(0);
			}

			if (cmd.hasOption("mappingFile") || cmd.hasOption("m")) {
				// initialise the member variable
				mappingFile = cmd.getOptionValue("mappingFile");
			}
			rank = cmd.hasOption("r");
			removeNanRow = cmd.hasOption("R");

		} catch (org.apache.commons.cli.ParseException ex) {
			Logger.getLogger(ConvertDoubleMatrixDataToTriTyper.class.getName()).log(Level.SEVERE, null, ex);
		}


		if (!(new File(outputFolder).exists())) {
			Gpio.createDir(outputFolder);
		} else if (!(new File(outputFolder).isDirectory())) {
			System.out.println("Error file is already there but not a directory!");
			System.exit(0);
		}

		HashSet<String> hashCpGSites = new HashSet<String>();
		try {
			System.out.println("Writing SNPMappings.txt to file:");
			int nrSites = 0;
			java.io.BufferedReader in = new java.io.BufferedReader(new java.io.FileReader(new File(mappingFile)));
			java.io.BufferedWriter outSNPMappings = new java.io.BufferedWriter(new java.io.FileWriter(new File(outputFolder + "/SNPMappings.txt")));
			String str;
			in.readLine();
			while ((str = in.readLine()) != null) {
				String[] data = SPLIT_ON_TAB.split(str);
				hashCpGSites.add(data[1]);
				outSNPMappings.write(data[3] + '\t' + data[4] + '\t' + data[1] + '\n');
				nrSites++;
			}
			System.out.println("Number of sites in mapping file:\t" + nrSites);
			outSNPMappings.close();
		} catch (Exception e) {
			System.out.println("Error:\t" + e.getMessage());
			e.printStackTrace();
			System.exit(0);
		}


		System.out.println("Input: " + dataMatrix);
		TextFile tf = new TextFile(dataMatrix, TextFile.R);
		String[] header = tf.readLineElems(TextFile.tab);
		System.out.println("Writing individuals and phenotype info: " + outputFolder + "/Individuals.txt");
		BufferedWriter outIndNew = new BufferedWriter(new FileWriter(outputFolder + "/Individuals.txt"));
		BufferedWriter outPhenoNew = new BufferedWriter(new FileWriter(outputFolder + "/PhenotypeInformation.txt"));
		for (int i = 0; i < header.length; i++) {
			outIndNew.write(header[i] + '\n');
			outPhenoNew.write(header[i] + "\tcontrol\tinclude\tmale\n");
		}
		outIndNew.close();
		outPhenoNew.close();


		System.out.println("Now processing matrix.");
		TextFile tfsnp = new TextFile(outputFolder + "/SNPs.txt", TextFile.W);

		BinaryFile bfgt = new BinaryFile(outputFolder + "/GenotypeMatrix.dat", BinaryFile.W);
		BinaryFile bfdo = new BinaryFile(outputFolder + "/ImputedDosageMatrix.dat", BinaryFile.W);
		String[] elems = tf.readLineElems(TextFile.tab);
		byte[] alleles = new byte[2];
		alleles[0] = 84;
		alleles[1] = 67;

		int nrRowsRead = 0;
		int nrRowsWritten = 0;
		while (elems != null) {
			String cpg = elems[0];
			if (hashCpGSites == null || hashCpGSites.contains(cpg)) {
				double[] vals = new double[elems.length - 1];
				for (int i = 1; i < elems.length; i++) {
					double d = Double.parseDouble(elems[i]);
					vals[i - 1] = d;
				}

				if (rank) {
					vals = rankRow(vals);
				}
				// rescale
				vals = rescaleValue(vals, 200d);
				int nrSamples = elems.length - 1;

				byte[] allele1 = new byte[nrSamples];
				byte[] allele2 = new byte[nrSamples];
				byte[] dosageValues = new byte[nrSamples];
				for (int ind = 0; ind < nrSamples; ind++) {
					if (vals[ind] > 100) {
						allele1[ind] = alleles[1];
						allele2[ind] = alleles[1];
					} else {
						allele1[ind] = alleles[0];
						allele2[ind] = alleles[0];
					}

					int dosageInt = (int) Math.round(vals[ind]);
					byte value = (byte) (Byte.MIN_VALUE + dosageInt);
					dosageValues[ind] = value;
				}

				bfgt.write(allele1);
				bfgt.write(allele2);
				bfdo.write(dosageValues);

				tfsnp.writeln(cpg);
				nrRowsWritten++;
			}

			nrRowsRead++;

			if (nrRowsRead % 10000 == 0) {
				System.out.println(nrRowsRead + " rows read, " + nrRowsWritten + " written.");
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		tfsnp.close();
		bfgt.close();
		bfdo.close();

		System.out.println("Finished.");
	}


	public static double[] rescaleValue(double[] row, Double multiplier) {
		double[] output = new double[row.length];
		double min = Primitives.min(row);
		double max = Primitives.max(row);
		double denominator = max - min;

		if (multiplier != null) {
			denominator = (max - min) / (1 / multiplier);
		}

		for (int i = 0; i < output.length; i++) {
			output[i] = (row[i] - min) / denominator;
		}
		return output;
	}

	public static DoubleMatrix2D rescaleValue(DoubleMatrix2D matrix, Double multiplier) {
		if (multiplier != null) {
			for (int p = 0; p < matrix.rows(); p++) {
				double min = matrix.viewRow(p).getMinLocation()[0];
				double denominator = (matrix.viewRow(p).getMaxLocation()[0] - min) * (1 / multiplier);
				for (int s = 0; s < matrix.columns(); s++) {
					matrix.setQuick(p, s, ((matrix.getQuick(p, s) - min) / denominator));
				}
			}
		} else {
			for (int p = 0; p < matrix.rows(); p++) {
				double min = matrix.viewRow(p).getMinLocation()[0];
				double denominator = matrix.viewRow(p).getMaxLocation()[0] - min;
				for (int s = 0; s < matrix.columns(); s++) {
					matrix.setQuick(p, s, ((matrix.getQuick(p, s) - min) / denominator));
				}
			}
		}
		return matrix;
	}

	public static DoubleMatrix2D rankRows(DoubleMatrix2D matrix) {

		RankArray rda = new RankArray();
		for (int p = 0; p < matrix.rows(); p++) {
			double[] rankedValues = rda.rank(matrix.viewRow(p).toArray(), true);
			for (int s = 0; s < matrix.columns(); s++) {
				matrix.setQuick(p, s, rankedValues[s]);
			}
		}
		return matrix;
	}

	public static double[] rankRow(double[] row) {

		RankArray rda = new RankArray();
		return rda.rank(row, true);

	}

}
