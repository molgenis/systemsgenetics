/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.metaqtl3;

import cern.colt.matrix.tint.impl.DenseIntMatrix2D;
import cern.jet.math.tint.IntFunctions;
import eqtlmappingpipeline.metaqtl3.graphics.QQPlot;
import gnu.trove.map.hash.TDoubleIntHashMap;
import gnu.trove.set.hash.TDoubleHashSet;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.eQTLTextFile;
import umcg.genetica.text.Strings;
import umcg.genetica.util.RankDoubleArray;

/**
 *
 * @author harmjan
 */
public class FDR {

//    public static String permutationDir = null;
//    public static String outputDir = null;
	public enum FDRMethod {

		PROBELEVEL, GENELEVEL, FULL
	};

	public enum FileFormat {

		LARGE, REDUCED
	};

	/**
	 * calculate the FalseDiscoveryRate for the discovered eQTLS
	 *
	 * @param eQTLTextFileLoc the location where the eQTL text files are stored
	 * @param nrPermutationsFDR number of permutations performed
	 * @param maxNrMostSignificantEQTLs maximum number of eQTLs to output
	 * @param fdrcutoff the FDR cutoff
	 * @param createQQPlot create a QQ plot after performing FDR calculations
	 * @param outputDir set an alternate directory for output
	 * @param permutationDir set an alternate directory for permutation files
	 * @throws IOException
	 */
	public static void calculateFDR(String eQTLTextFileLoc, int nrPermutationsFDR, int maxNrMostSignificantEQTLs, double fdrcutoff, boolean createQQPlot, String outputDir, String permutationDir) throws IOException {

		if (eQTLTextFileLoc == null || eQTLTextFileLoc.length() == 0) {
			throw new IllegalArgumentException("File containing real effects is not specified.");
		}
		if (nrPermutationsFDR < 1) {
			throw new IllegalArgumentException("Need at least one permutation to determine FDR");
		}
		if (maxNrMostSignificantEQTLs < 1) {
			throw new IllegalArgumentException("Need at least a single effect to perform FDR estimation");
		}
		if (fdrcutoff < 0 || fdrcutoff > 1) {
			throw new IllegalArgumentException("FDR threshold should be between 0.0 and 1.0! (Specified: " + fdrcutoff + ")");
		}

		//Load permuted data:
//        // load values for each permutation round:
		if (permutationDir == null) {
			permutationDir = eQTLTextFileLoc;
		}

		if (outputDir == null) {
			outputDir = eQTLTextFileLoc;
		}

		String fileString = permutationDir + "/PermutedEQTLsPermutationRound" + 1 + ".txt.gz";
		TextFile tf = new TextFile(fileString, TextFile.R);
		tf.readLine();
		String[] elems = tf.readLineElems(TextFile.tab);
		int nrColsInPermutedFiles = elems.length;
		tf.close();

		System.out.println(nrColsInPermutedFiles + " columns in permuted QTL file.");
		// new permutationfile format requires different column layout...
		if (nrColsInPermutedFiles > 7) {
			System.out.println("Large permutation files detected.");
			runFDR(eQTLTextFileLoc, nrPermutationsFDR, maxNrMostSignificantEQTLs, fdrcutoff, FileFormat.LARGE, FDRMethod.FULL, outputDir, permutationDir, createQQPlot);
			runFDR(eQTLTextFileLoc, nrPermutationsFDR, maxNrMostSignificantEQTLs, fdrcutoff, FileFormat.LARGE, FDRMethod.PROBELEVEL, outputDir, permutationDir, createQQPlot);
			runFDR(eQTLTextFileLoc, nrPermutationsFDR, maxNrMostSignificantEQTLs, fdrcutoff, FileFormat.LARGE, FDRMethod.GENELEVEL, outputDir, permutationDir, createQQPlot);
		} else {
			runFDR(eQTLTextFileLoc, nrPermutationsFDR, maxNrMostSignificantEQTLs, fdrcutoff, FileFormat.REDUCED, FDRMethod.FULL, outputDir, permutationDir, createQQPlot);
			runFDR(eQTLTextFileLoc, nrPermutationsFDR, maxNrMostSignificantEQTLs, fdrcutoff, FileFormat.REDUCED, FDRMethod.PROBELEVEL, outputDir, permutationDir, createQQPlot);
			if (nrColsInPermutedFiles >= 4) {
				runFDR(eQTLTextFileLoc, nrPermutationsFDR, maxNrMostSignificantEQTLs, fdrcutoff, FileFormat.REDUCED, FDRMethod.GENELEVEL, outputDir, permutationDir, createQQPlot);
			}
		}
	}

	private static void runFDR(String baseDir, int nrPermutationsFDR, int maxNrMostSignificantEQTLs,
			double fdrcutoff, FileFormat f, FDRMethod m, String outputDir, String permutationDir, boolean createQQPlot) throws IOException {
		//Load permuted data:
		// load values for each permutation round:
		System.out.println("");
		if (m == FDRMethod.GENELEVEL) {
			System.out.println("Performing gene level FDR");
		} else if (m == FDRMethod.PROBELEVEL) {
			System.out.println("Performing probe level FDR");
		} else if (m == FDRMethod.FULL) {
			System.out.println("Determining the FDR using all data");
		}


		TDoubleIntHashMap permutedPvalues = new TDoubleIntHashMap(10000, 0.5f);

//        ProgressBar pb = new ProgressBar(nrPermutationsFDR, "Reading permuted data:");
		System.out.println("Reading permuted files");

		for (int permutationRound = 0; permutationRound < nrPermutationsFDR; permutationRound++) {
			String fileString = permutationDir + "/PermutedEQTLsPermutationRound" + (permutationRound + 1) + ".txt.gz";
			System.out.println(fileString);
			// read the permuted eqtl output
			TextFile gz = new TextFile(fileString, TextFile.R);


			String[] header = gz.readLineElems(TextFile.tab);
			int snpcol = -1;
			int pvalcol = -1;
			int probecol = -1;
			int genecol = -1;

			if (f == FileFormat.REDUCED) {



				//PValue  SNP     Probe   Gene 
				for (int col = 0; col < header.length; col++) {
					if (header[col].equals("PValue")) {
						pvalcol = col;
					}
					if (header[col].equals("SNP")) {
						snpcol = col;
					}
					if (header[col].equals("Probe")) {
						probecol = col;
					}
					if (header[col].equals("Gene")) {
						genecol = col;
					}
				}

				//PValue  SNP     Probe   Gene
				if (snpcol == -1 || pvalcol == -1 || probecol == -1 && genecol == -1) {
					System.out.println("Column not found in permutation file: " + fileString);
					System.out.println("PValue: " + pvalcol);
					System.out.println("SNP: " + snpcol);
					System.out.println("Probe: " + probecol);
					System.out.println("Gene: " + genecol);
				}
			}
			String[] data = gz.readLineElemsReturnReference(TextFile.tab);
			int itr = 0;

			HashSet<String> visitedEffects = new HashSet<String>();
			while (data != null) {

				if (data.length != 0) {
					if (itr > maxNrMostSignificantEQTLs - 1) {
						System.out.println("Breaking because: " + itr);
						break;
					} else {
						int filteronColumn;
						String fdrId = null;
						if (f == FileFormat.REDUCED) {
							if (m == FDRMethod.FULL) {
								//fdrId = data[snpcol] + "-" + data[probecol];
								filteronColumn = probecol;
							} else if (m == FDRMethod.GENELEVEL && data.length > 3) {
								fdrId = data[genecol];
								filteronColumn = genecol;
							} else {
								fdrId = data[probecol];
								filteronColumn = probecol;
							}

						} else {
							if (m == FDRMethod.GENELEVEL) {
								fdrId = data[eQTLTextFile.HUGO];
								filteronColumn = eQTLTextFile.HUGO;
							} else if (m == FDRMethod.PROBELEVEL) {
								fdrId = data[4];
								filteronColumn = 4;
							}
						}

						// take top effect per gene / probe
						if (m == FDRMethod.FULL || (!fdrId.equals("-") && !visitedEffects.contains(fdrId))) {

							if (m != FDRMethod.FULL) {
								visitedEffects.add(fdrId);
							}

							double permutedP = Double.parseDouble(data[0]);
							if (permutedPvalues.containsKey(permutedP)) {
								permutedPvalues.increment(permutedP);
							} else {
								permutedPvalues.put(permutedP, 1);
							}

							itr++;
						}

						data = gz.readLineElemsReturnReference(TextFile.tab);
					}
				}
			}
			gz.close();


		}




		double[] uniquePermutedPvalues = permutedPvalues.keys();
		Arrays.sort(uniquePermutedPvalues);

		double[] uniquePermutedPvaluesCounts = new double[uniquePermutedPvalues.length];

		int cummulativeCount = 0;
		double nrPermutationsFDRd = (double) nrPermutationsFDR;
		for (int i = 0; i < uniquePermutedPvalues.length; ++i) {

			cummulativeCount += permutedPvalues.get(uniquePermutedPvalues[i]);
			uniquePermutedPvaluesCounts[i] = cummulativeCount / nrPermutationsFDRd;

		}


		String outFileName = "";
		String outFileNameAll = "";

		if (outputDir == null) {
			outputDir = baseDir;
		}

		if (m == FDRMethod.GENELEVEL) {
			outFileName = outputDir + "/eQTLsFDR" + fdrcutoff + "-GeneLevel.txt";
			outFileNameAll = outputDir + "/eQTLsFDR-GeneLevel.txt.gz";
		} else if (m == FDRMethod.PROBELEVEL) {
			outFileName = outputDir + "/eQTLsFDR" + fdrcutoff + "-ProbeLevel.txt";
			outFileNameAll = outputDir + "/eQTLsFDR-ProbeLevel.txt.gz";
		} else {
			outFileName = outputDir + "/eQTLsFDR" + fdrcutoff + ".txt";
			outFileNameAll = outputDir + "/eQTLsFDR.txt.gz";
		}

		BufferedWriter outputWriterSignificant = new BufferedWriter(new FileWriter(outFileName));
		BufferedWriter outputWriterAll = new BufferedWriter(new FileWriter(outFileNameAll));

		String fileString = baseDir + "/eQTLs.txt.gz";
		if (!Gpio.exists(fileString)) {
			System.out.println("Could not find file: " + fileString + " trying un-GZipped file....");
			fileString = baseDir + "/eQTLs.txt";
		}
		if (!Gpio.exists(fileString)) {
			System.out.println("Could not find file: " + fileString);
			System.exit(0);
		}


		TextFile realEQTLs = new TextFile(fileString, TextFile.R);


		String header = realEQTLs.readLine();

		outputWriterAll.append(header);
		outputWriterAll.append("\tFDR\n");

		outputWriterSignificant.append(header);
		outputWriterSignificant.append("\tFDR\n");

		String str = realEQTLs.readLine();

// REAL DATA PROCESSING
		int itr = 0;
		HashSet<String> visitedEffects = new HashSet<String>();
		double lastEqtlPvalue = 0;

		double currentPvalue = 0;
		ArrayList<String> currentPvalueEqtls = new ArrayList<String>();


		int lastUsedPermutedPvalueIndex = 0;


		int nrSignificantEQTLs = 0;

		while (str != null) {
			if (itr > maxNrMostSignificantEQTLs - 1) {
				break;
			} else {

				String fdrId = null;
				String[] data = Strings.tab.split(str);

				if (m == FDRMethod.GENELEVEL) {
					fdrId = data[eQTLTextFile.HUGO];
				} else if (m == FDRMethod.PROBELEVEL) {
					fdrId = data[4];
				}


				if (m == FDRMethod.FULL || (!fdrId.equals("-") && !visitedEffects.contains(fdrId))) {




					double eQtlPvalue = Double.parseDouble(data[0]);

					if (itr > 0 && lastEqtlPvalue > eQtlPvalue) {
						System.err.println("Sorted P-Value list is not perfectly sorted!!!!");
						System.exit(-1);
					}



					if (eQtlPvalue > currentPvalue) {
						//Process old results for current pvalue


						while (uniquePermutedPvalues[lastUsedPermutedPvalueIndex] <= currentPvalue) {
							++lastUsedPermutedPvalueIndex;
						}
						double fdr = uniquePermutedPvaluesCounts[lastUsedPermutedPvalueIndex] / itr;

						for (String cachedEqtls : currentPvalueEqtls) {
							outputWriterAll.append(cachedEqtls);
							outputWriterAll.append('\t');
							outputWriterAll.append(String.valueOf(fdr));
							outputWriterAll.append('\n');

							if (fdr <= fdrcutoff) {
								outputWriterSignificant.append(cachedEqtls);
								outputWriterSignificant.append('\t');
								outputWriterSignificant.append(String.valueOf(fdr));
								outputWriterSignificant.append('\n');
								++nrSignificantEQTLs;
							}

						}


						//Create new temp list for this pvalue
						currentPvalue = eQtlPvalue;
						currentPvalueEqtls.clear();
						currentPvalueEqtls.add(str);

					} else {
						//add to current pvalue list
						currentPvalueEqtls.add(str);
					}



					lastEqtlPvalue = eQtlPvalue;
					visitedEffects.add(fdrId);
					itr++;
				}

				str = realEQTLs.readLine();
			}

		}

		//Write buffer to files
		while (uniquePermutedPvalues[lastUsedPermutedPvalueIndex] <= currentPvalue) {
			++lastUsedPermutedPvalueIndex;
		}
		double fdr = uniquePermutedPvaluesCounts[lastUsedPermutedPvalueIndex] / itr;

		for (String cachedEqtls : currentPvalueEqtls) {
			outputWriterAll.append(cachedEqtls);
			outputWriterAll.append('\t');
			outputWriterAll.append(String.valueOf(fdr));
			outputWriterAll.append('\n');

			if (fdr <= fdrcutoff) {
				outputWriterSignificant.append(cachedEqtls);
				outputWriterSignificant.append('\t');
				outputWriterSignificant.append(String.valueOf(fdr));
				outputWriterSignificant.append('\n');
				++nrSignificantEQTLs;
			}

		}




		realEQTLs.close();
		outputWriterAll.close();
		outputWriterSignificant.close();

		//Process a certain P-Value, determine how many are significant:
		//System.out.println("\n");
		//System.out.println("Significant detected eQTLs:");


		//System.out.println("");
		String output = "Number of significant eQTLs:\t" + nrSignificantEQTLs;
		System.out.println(output);

		String fileSuffix = "";
		if (m == FDRMethod.GENELEVEL) {
			fileSuffix = "-GeneLevel";
		} else if (m == FDRMethod.PROBELEVEL) {
			fileSuffix = "-ProbeLevel";
		}

		if (createQQPlot) {

			System.err.println("Sorry, QQ plot function is temporarily (or for a very long time) unavailable.");

			//System.out.println("Creating QQ plot. This might take a while...");
			//QQPlot qq = new QQPlot();
			//String fileName = baseDir + "/eQTLsFDR" + fdrcutoff + fileSuffix + "-QQPlot.pdf";
			//qq.draw(fileName, fdrcutoff, nrPermutationsFDR,
			//		maxNrMostSignificantEQTLs, permutedPValues.toArray(), nrRealDataEQTLs, pValues,
			//		pValueSignificant, nrSignificantEQTLs);
		}

		generateESNPsFile(outputDir + "/eQTLsFDR" + fdrcutoff + fileSuffix + ".txt", outputDir + "/eQTLSNPsFDR" + fdrcutoff + fileSuffix + ".txt");
		generateEProbesFile(outputDir + "/eQTLsFDR" + fdrcutoff + fileSuffix + ".txt", outputDir + "/eQTLProbesFDR" + fdrcutoff + fileSuffix + ".txt");

	}

	/**
	 * Generates an eQTL SNPs file
	 *
	 * @param inputFile
	 * @param outputFile
	 */
	private static void generateESNPsFile(String inputFile, String outputFile) throws IOException {

		//Determine the number of unique SNPs and what the most significant eQTL for this SNP is:
		HashMap<String, Double> hashTopSNPsAbsZScore = new HashMap<String, Double>();
		HashMap<String, String> hashTopSNPsAnnotation = new HashMap<String, String>();
		ArrayList<String> vecTopSNPs = new ArrayList<String>();
		String header = "";

		TextFile in = new TextFile(inputFile, TextFile.R);
		String str = in.readLine();
		header = str;
		while ((str = in.readLine()) != null) {
			String[] data = str.split("\t");
			double absZScore = Math.abs(Double.parseDouble(data[10]));
			String rsName = data[1];
			if (hashTopSNPsAbsZScore.get(rsName) != null) {
				double absZScorePrevious = hashTopSNPsAbsZScore.get(rsName);
				if (absZScore > absZScorePrevious) {
					hashTopSNPsAbsZScore.put(rsName, absZScore);
					hashTopSNPsAnnotation.put(rsName, str);
				}
			} else {
				vecTopSNPs.add(rsName);
				hashTopSNPsAbsZScore.put(rsName, absZScore);
				hashTopSNPsAnnotation.put(rsName, str);
			}
		}
		in.close();

		//Per eSNP write the most significant eQTL that has been recorded to file:
		TextFile out = new TextFile(outputFile, TextFile.W);
		out.writeln(header);
		for (int s = 0; s < vecTopSNPs.size(); s++) {
			String rsName = vecTopSNPs.get(s);
			String annotation = hashTopSNPsAnnotation.get(rsName);
			out.writeln(annotation);
		}
		out.close();

		String output = " - Number of unique SNPs, constituting an eQTL:\t" + vecTopSNPs.size();
		System.out.println(output);

	}

	/**
	 * Generates a file containing combinations of the most significant eQTL
	 * with probes
	 *
	 * @param inputFile location of text file containing data on eQTLs
	 * @param outputFile the location where the output will be written
	 */
	private static void generateEProbesFile(String inputFile, String outputFile) throws IOException {

		//Determine the number of unique probes and what the most significant eQTL for each probe is:
		HashMap<String, Double> hashTopProbesAbsZScore = new HashMap<String, Double>();
		HashMap<String, String> hashTopProbesAnnotation = new HashMap<String, String>();
		ArrayList<String> vecTopProbes = new ArrayList<String>();
		String header = "";

		TextFile in = new TextFile(inputFile, TextFile.R);
		String str = in.readLine();
		header = str;
		while ((str = in.readLine()) != null) {
			String[] data = str.split("\t");
			double absZScore = Math.abs(Double.parseDouble(data[10]));
			String probe = data[4];
			if (hashTopProbesAbsZScore.get(probe) != null) {
				double absZScorePrevious = hashTopProbesAbsZScore.get(probe);
				if (absZScore > absZScorePrevious) {
					hashTopProbesAbsZScore.put(probe, absZScore);
					hashTopProbesAnnotation.put(probe, str);
				}
			} else {
				vecTopProbes.add(probe);
				hashTopProbesAbsZScore.put(probe, absZScore);
				hashTopProbesAnnotation.put(probe, str);
			}
		}
		in.close();

		//Per eSNP write the most significant eQTL that has been recorded to file:
		TextFile out = new TextFile(outputFile, TextFile.W);
		out.writeln(header);
		for (int s = 0; s < vecTopProbes.size(); s++) {
			String probe = vecTopProbes.get(s);
			String annotation = hashTopProbesAnnotation.get(probe);
			out.writeln(annotation);

		}

		out.close();

		String output = " - Number of unique probes, constituting an eQTL:\t" + vecTopProbes.size();
		System.out.println(output);
	}
}
