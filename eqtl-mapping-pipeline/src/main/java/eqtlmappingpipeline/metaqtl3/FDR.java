/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.metaqtl3;

import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.tdouble.impl.DenseLargeDoubleMatrix2D;
import eqtlmappingpipeline.metaqtl3.graphics.QQPlot;
import gnu.trove.list.array.TCharArrayList;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.map.hash.TDoubleIntHashMap;
import gnu.trove.set.hash.THashSet;

import java.io.IOException;
import java.lang.management.MemoryUsage;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.concurrent.*;

import umcg.genetica.console.MultiThreadProgressBar;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.QTLTextFile;
import umcg.genetica.text.Strings;

/**
 * @author harmjan, Patrick Deelen, Marc Jan Bonder
 */
public class FDR {

	//    public static String permutationDir = null;
//    public static String outputDir = null;
	public enum FDRMethod {
		PROBELEVEL, SNPLEVEL, GENELEVEL, FULL, ALL
	}

	;

	public enum FileFormat {

		LARGE, REDUCED
	}

	;

	/**
	 * calculate the FalseDiscoveryRate for the discovered eQTLS
	 *
	 * @param eQTLTextFileLoc           the location where the eQTL text files are stored
	 * @param nrPermutationsFDR         number of permutations performed
	 * @param maxNrMostSignificantEQTLs maximum number of eQTLs to output
	 * @param fdrcutoff                 the FDR cutoff
	 * @param createQQPlot              create a QQ plot after performing FDR calculations
	 * @param outputDir                 set an alternate directory for output
	 * @param permutationDir            set an alternate directory for permutation files
	 * @throws IOException
	 */
	public static void calculateFDR(String eQTLTextFileLoc, int nrPermutationsFDR, int maxNrMostSignificantEQTLs, double fdrcutoff, boolean createQQPlot, String outputDir, String permutationDir, FDRMethod fdrType, boolean createLargeFdrFiles) throws IOException {

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
		int nrColsInPermutedFiles = 0;
		if (elems != null) {
			nrColsInPermutedFiles = elems.length;
		}
		tf.close();


		System.out.println(nrColsInPermutedFiles + " columns in permuted QTL file.");
		if (nrColsInPermutedFiles > 0) {
			// new permutationfile format requires different column layout...
			FileFormat format = FileFormat.LARGE;
			if (nrColsInPermutedFiles <= 7) {
				System.out.println("Reduced format permutation files detected.");
				format = FileFormat.REDUCED;
			}

			if (fdrType.equals(FDRMethod.ALL)) {
				runFDR(eQTLTextFileLoc, nrPermutationsFDR, maxNrMostSignificantEQTLs, fdrcutoff, format, FDRMethod.FULL, outputDir, permutationDir, createQQPlot, createLargeFdrFiles);
				runFDR(eQTLTextFileLoc, nrPermutationsFDR, maxNrMostSignificantEQTLs, fdrcutoff, format, FDRMethod.GENELEVEL, outputDir, permutationDir, createQQPlot, createLargeFdrFiles);
				runFDR(eQTLTextFileLoc, nrPermutationsFDR, maxNrMostSignificantEQTLs, fdrcutoff, format, FDRMethod.PROBELEVEL, outputDir, permutationDir, createQQPlot, createLargeFdrFiles);
				runFDR(eQTLTextFileLoc, nrPermutationsFDR, maxNrMostSignificantEQTLs, fdrcutoff, format, FDRMethod.SNPLEVEL, outputDir, permutationDir, createQQPlot, createLargeFdrFiles);
			} else {
				runFDR(eQTLTextFileLoc, nrPermutationsFDR, maxNrMostSignificantEQTLs, fdrcutoff, format, fdrType, outputDir, permutationDir, createQQPlot, createLargeFdrFiles);
			}
		}

	}

	//
	private static void runFDR(String baseDir, int nrPermutationsFDR, int maxNrMostSignificantEQTLs,
							   double fdrcutoff, FileFormat f, FDRMethod m, String outputDir, String permutationDir, boolean createQQPlot, boolean createLargeFdrFiles) throws IOException {
		//Load permuted data:
		// load values for each permutation round:
		System.out.println("");

		if (m == FDRMethod.GENELEVEL) {
			System.out.println("Performing gene level FDR");
		} else if (m == FDRMethod.PROBELEVEL) {
			System.out.println("Performing probe level FDR");
		} else if (m == FDRMethod.SNPLEVEL) {
			System.out.println("Performing SNP level FDR");
		} else if (m == FDRMethod.FULL) {
			System.out.println("Determining the FDR using all data");
		}


//        ProgressBar pb = new ProgressBar(nrPermutationsFDR, "Reading permuted data:");
		System.out.println("Reading permuted files from: " + permutationDir + ".\nAssuming " + maxNrMostSignificantEQTLs + " effects");
		TDoubleIntHashMap permutedPvalues = new TDoubleIntHashMap(10000, 0.5f);
		Executor ex = Executors.newWorkStealingPool();
		ExecutorCompletionService<TDoubleIntHashMap> service = new ExecutorCompletionService<>(ex);
		int submit = 0;
		MultiThreadProgressBar pb = new MultiThreadProgressBar(nrPermutationsFDR);
		for (int permutationRound = 0; permutationRound < nrPermutationsFDR; permutationRound++) {
			// run files in parallel
			// String permutationDir, int permutationRound, FileFormat f, int maxNrMostSignificantEQTLs, FDRMethod m
			ReadPermutationFile t = new ReadPermutationFile(permutationDir, permutationRound, f, maxNrMostSignificantEQTLs, m, pb);
			service.submit(t);
			submit++;
		}

		while (!pb.allCompleted()) {
			pb.display();
			try {
				Thread.sleep(2000);
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}

		System.out.println();
		int returned = 0;
		while (returned < submit) {
			try {
				Future<TDoubleIntHashMap> h = service.take();
				TDoubleIntHashMap result = h.get();
				// update final hash
				if (h != null) {
					for (double d : result.keys()) {
						int val;
						if (!permutedPvalues.containsKey(d)) {
							val = 0;
						} else {
							val = permutedPvalues.get(d);
						}
						permutedPvalues.put(d, val + result.get(d));
					}
					returned++;
				}
			} catch (InterruptedException e) {
				e.printStackTrace();
			} catch (ExecutionException e) {
				e.printStackTrace();
			}
		}

		int nreqtl = 0;
		for (double k : permutedPvalues.keys()) {
			nreqtl += permutedPvalues.get(k);
		}
		System.out.println("After loading: " + submit + " permutation files, found " + permutedPvalues.keys().length + " unique p-values and " + nreqtl + " overall eqtl, or " + (nreqtl / nrPermutationsFDR) + " per file.");
		double[] uniquePermutedPvalues = permutedPvalues.keys();
		Arrays.sort(uniquePermutedPvalues);

		double[] uniquePermutedPvaluesCounts = new double[uniquePermutedPvalues.length];

		long cummulativeCount = 0;
		double nrPermutationsFDRd = (double) nrPermutationsFDR;
		for (int i = 0; i < uniquePermutedPvalues.length; ++i) {
			cummulativeCount += permutedPvalues.get(uniquePermutedPvalues[i]);
			uniquePermutedPvaluesCounts[i] = cummulativeCount / nrPermutationsFDRd;
		}
		permutedPvalues = null;
		System.out.println("Number of unique permutation p-values: " + uniquePermutedPvalues.length);

		if (outputDir == null) {
			outputDir = baseDir;
		}

		String fileSuffix = "";
		if (m == FDRMethod.GENELEVEL) {
			fileSuffix = "-GeneLevel";
		} else if (m == FDRMethod.PROBELEVEL) {
			fileSuffix = "-ProbeLevel";
		} else if (m == FDRMethod.SNPLEVEL) {
			fileSuffix = "-SNPLevel";
		}

		String outFileName = outputDir + "/eQTLsFDR" + fdrcutoff + fileSuffix + ".txt.gz";
		String outFileNameSnps = outputDir + "/eQTLSNPsFDR" + fdrcutoff + fileSuffix + ".txt.gz";
		String outFileNameProbes = outputDir + "/eQTLProbesFDR" + fdrcutoff + fileSuffix + ".txt.gz";
		String outFileNameAll = outputDir + "/eQTLsFDR" + fileSuffix + ".txt.gz";

		TextFile outputWriterSignificant = new TextFile(outFileName, TextFile.W, 10 * 1048576);
		TextFile outputWriterESNPs = new TextFile(outFileNameSnps, TextFile.W, 10 * 1048576);
		TextFile outputWriterEProbes = new TextFile(outFileNameProbes, TextFile.W, 10 * 1048576);
		System.out.println("writing the following files: ");
		System.out.println(outFileName);
		System.out.println(outFileNameSnps);
		System.out.println(outFileNameProbes);
		TextFile outputWriterAll = null;
		if (createLargeFdrFiles) {
			System.out.println("Will also output large FDR file here: " + outFileNameAll);
			outputWriterAll = new TextFile(outFileNameAll, TextFile.W, 10 * 1048576);
		}
		String fileString = baseDir + "/eQTLs.txt.gz";
		if (!Gpio.exists(fileString)) {
			System.out.println("Could not find file: " + fileString + " trying un-GZipped file....");
			fileString = baseDir + "/eQTLs.txt";
		}
		if (!Gpio.exists(fileString)) {
			System.out.println("Could not find file: " + fileString);
			System.exit(0);
		}


		TextFile realEQTLs = new TextFile(fileString, TextFile.R, 10 * 1048576);

		String header = realEQTLs.readLine();
		// CHECK IF THERE IS ALREADY AN FDR column
		String[] headerElems = header.split("\t");
		int fdrcol = -1;
		for (int h = 0; h < headerElems.length; h++) {
			if (headerElems[h].toLowerCase().equals("fdr")) {
				fdrcol = h;
			}
		}

		if (fdrcol > -1) {
			System.out.println("Found FDR at column: " + fdrcol);
			if (createLargeFdrFiles) {
				outputWriterAll.append(header);
				outputWriterAll.append("\n");
			}
			outputWriterEProbes.append(header);
			outputWriterEProbes.append("\n");

			outputWriterESNPs.append(header);
			outputWriterESNPs.append("\n");

			outputWriterSignificant.append(header);
			outputWriterSignificant.append("\n");
		} else {
			if (createLargeFdrFiles) {
				outputWriterAll.append(header);
				outputWriterAll.append("\tFDR\n");
			}
			outputWriterEProbes.append(header);
			outputWriterEProbes.append("\tFDR\n");

			outputWriterESNPs.append(header);
			outputWriterESNPs.append("\tFDR\n");

			outputWriterSignificant.append(header);
			outputWriterSignificant.append("\tFDR\n");
		}


// REAL DATA PROCESSING
		int itr = 0;
		THashSet<String> visitedEffects = new THashSet<>();
		THashSet<String> visitedSnps = new THashSet<>();
		THashSet<String> visitedProbes = new THashSet<>();
		double lastEqtlPvalue = 0;

		double currentPvalue = 0;
		ArrayList<String> currentPvalueEqtls = new ArrayList<String>();
		ArrayList<String> currentPvalueEqtlSnps = new ArrayList<String>();
		ArrayList<String> currentPvalueEqtlProbes = new ArrayList<String>();

		TDoubleArrayList pValueRealData = new TDoubleArrayList();
		ArrayList<Boolean> significantPvalue = new ArrayList<Boolean>();
		int lastUsedPermutedPvalueIndex = 0;

		int nrSignificantEQTLs = 0;

		// start reading
		String str = realEQTLs.readLine();
		int lnctr = 0;
		while (str != null) {
			if (itr > maxNrMostSignificantEQTLs - 1) {
				break;
			} else {

				String fdrId = null;
//				String[] data = Strings.tab.split(str);

				if (m == FDRMethod.GENELEVEL) {
//					fdrId = data[QTLTextFile.HUGO];
					fdrId = Strings.subsplit(str, Strings.tab, QTLTextFile.HUGO, QTLTextFile.HUGO + 1)[0];
					fdrId = new String(fdrId.getBytes("UTF-8")).intern();
				} else if (m == FDRMethod.SNPLEVEL) {
//					fdrId = data[QTLTextFile.SNP];
					fdrId = Strings.subsplit(str, Strings.tab, QTLTextFile.SNP, QTLTextFile.SNP + 1)[0];
					fdrId = new String(fdrId.getBytes("UTF-8")).intern();
				} else if (m == FDRMethod.PROBELEVEL) {
//					fdrId = data[4];
					fdrId = Strings.subsplit(str, Strings.tab, 4, 5)[0];
					fdrId = new String(fdrId.getBytes("UTF-8")).intern();
				}

				double eQtlPvalue = Double.parseDouble(Strings.subsplit(str, Strings.tab, 0, 1)[0]);

				if (itr > 0 && lastEqtlPvalue > eQtlPvalue) {
					System.err.println("Sorted P-Value list is not perfectly sorted!!!!");
					System.out.println("Last p: " + lastEqtlPvalue + "\tCurrent: " + eQtlPvalue);
					System.exit(-1);
				}

				if (eQtlPvalue > currentPvalue) {
					//Process old results for current pvalue

					double fdr = 0;
					if (currentPvalue >= uniquePermutedPvalues[0]) {

						while (uniquePermutedPvalues[lastUsedPermutedPvalueIndex + 1] <= currentPvalue && lastUsedPermutedPvalueIndex < uniquePermutedPvalues.length - 2) {
							++lastUsedPermutedPvalueIndex;
						}
						fdr = uniquePermutedPvaluesCounts[lastUsedPermutedPvalueIndex] / itr;

						if (fdr > 1) {
							fdr = 1;
						}

					}

					for (int i = 0; i < currentPvalueEqtls.size(); ++i) {
						String cachedEqtls = currentPvalueEqtls.get(i);
						String cachedEqtlsProbe = currentPvalueEqtlProbes.get(i);
						String cachedEqtlsSnps = currentPvalueEqtlSnps.get(i);

						StringBuilder currentString = new StringBuilder();

						if (fdrcol > -1) {
							// there is already an FDR column; replace
							String[] elems = Strings.subsplit(cachedEqtls, Strings.tab, 0, headerElems.length);
							elems[fdrcol] = "" + fdr;
							currentString = new StringBuilder(Strings.concat(elems, Strings.tab)).append("\n");
						} else {
							currentString.append(cachedEqtls).append('\t').append(String.valueOf(fdr)).append('\n');
						}

						pValueRealData.add(currentPvalue);

						if (createLargeFdrFiles) {
							outputWriterAll.append(currentString.toString());
						}

						if (fdr <= fdrcutoff) {
							if (!visitedProbes.contains(cachedEqtlsProbe)) {
								outputWriterEProbes.append(currentString.toString());
								visitedProbes.add(cachedEqtlsProbe);
							}
							if (!visitedSnps.contains(cachedEqtlsSnps)) {
								outputWriterESNPs.append(currentString.toString());
								visitedSnps.add(cachedEqtlsSnps);

							}

							significantPvalue.add(true);
							outputWriterSignificant.append(currentString.toString());
							++nrSignificantEQTLs;
						} else {
							significantPvalue.add(false);
						}

					}

					//Create new temp list for this pvalue
					currentPvalue = eQtlPvalue;
					currentPvalueEqtls.clear();
					currentPvalueEqtlProbes.clear();
					currentPvalueEqtlSnps.clear();
					currentPvalueEqtls.add(str);
					String probe = new String(Strings.subsplit(str, Strings.tab, QTLTextFile.PROBE, QTLTextFile.PROBE + 1)[0].getBytes("UTF-8")).intern();
					String snp = new String(Strings.subsplit(str, Strings.tab, QTLTextFile.SNP, QTLTextFile.SNP + 1)[0].getBytes("UTF-8")).intern();
					currentPvalueEqtlProbes.add(probe);
					currentPvalueEqtlSnps.add(snp);

				} else {
					//add to current pvalue list
					currentPvalueEqtls.add(str);
					String probe = new String(Strings.subsplit(str, Strings.tab, QTLTextFile.PROBE, QTLTextFile.PROBE + 1)[0].getBytes("UTF-8")).intern();
					String snp = new String(Strings.subsplit(str, Strings.tab, QTLTextFile.SNP, QTLTextFile.SNP + 1)[0].getBytes("UTF-8")).intern();
					currentPvalueEqtlProbes.add(probe);
					currentPvalueEqtlSnps.add(snp);

				}

				lastEqtlPvalue = eQtlPvalue;

				if (m == FDRMethod.FULL || (!fdrId.equals("-") && !visitedEffects.contains(fdrId))) {
					itr++;
					visitedEffects.add(fdrId);
				}

				str = realEQTLs.readLine();
				lnctr++;
				if (lnctr % 100000 == 0) {
					System.out.print(lnctr + " lines parsed\r");
				}
			}

		}
		System.out.println("");
		System.out.println("Done.");

		//Write buffer to files
		double fdr = 0;
		if (currentPvalue >= uniquePermutedPvalues[0]) {

			while (uniquePermutedPvalues[lastUsedPermutedPvalueIndex + 1] <= currentPvalue && lastUsedPermutedPvalueIndex < uniquePermutedPvalues.length - 2) {
				++lastUsedPermutedPvalueIndex;
			}
			fdr = uniquePermutedPvaluesCounts[lastUsedPermutedPvalueIndex] / itr;

			if (fdr > 1) {
				fdr = 1;
			}
		}

		// write last bit of buffer
		for (int i = 0; i < currentPvalueEqtls.size(); ++i) {
			String cachedEqtls = currentPvalueEqtls.get(i);
			String cachedEqtlsProbe = currentPvalueEqtlProbes.get(i);
			String cachedEqtlsSnps = currentPvalueEqtlSnps.get(i);

			StringBuilder currentString = new StringBuilder();
			if (fdrcol > -1) {
				// there is already an FDR column; replace
				String[] elems = Strings.subsplit(cachedEqtls, Strings.tab, 0, headerElems.length);
				elems[fdrcol] = "" + fdr;
				currentString = new StringBuilder(Strings.concat(elems, Strings.tab)).append("\n");
				;
			} else {
				currentString.append(cachedEqtls).append('\t').append(String.valueOf(fdr)).append('\n');
			}


			pValueRealData.add(currentPvalue);
			if (createLargeFdrFiles) {
				outputWriterAll.append(currentString.toString());
			}

			if (fdr <= fdrcutoff) {
				if (!visitedProbes.contains(cachedEqtlsProbe)) {
					outputWriterEProbes.append(currentString.toString());
					visitedProbes.add(cachedEqtlsProbe);
				}
				if (!visitedSnps.contains(cachedEqtlsSnps)) {
					outputWriterESNPs.append(currentString.toString());
					visitedSnps.add(cachedEqtlsSnps);
				}

				significantPvalue.add(true);
				outputWriterSignificant.append(currentString.toString());
				++nrSignificantEQTLs;
			} else {
				significantPvalue.add(false);
			}
		}

		realEQTLs.close();
		if (createLargeFdrFiles) {
			outputWriterAll.close();
		}
		outputWriterEProbes.close();
		outputWriterESNPs.close();
		outputWriterSignificant.close();

		//System.out.println("");
		System.out.println("Number of significant eQTLs:\t" + nrSignificantEQTLs);
		System.out.println(" - Number of unique SNPs, constituting an eQTL:\t" + visitedSnps.size());
		System.out.println(" - Number of unique probes, constituting an eQTL:\t" + visitedProbes.size());

		if (createQQPlot) {
			System.out.println("Creating QQ plot. This might take a while...");
			String fileName = baseDir + "/eQTLsFDR" + fdrcutoff + fileSuffix + "-QQPlot.pdf";
			if (maxNrMostSignificantEQTLs > pValueRealData.size()) {
				createQQPlots(permutationDir, nrPermutationsFDR, pValueRealData.size(), fdrcutoff, f, m, pValueRealData.toArray(), significantPvalue, nrSignificantEQTLs, fileName);
			} else if (maxNrMostSignificantEQTLs > 100000) {
				System.out.println("Only taking the top 100,000 for QQplot creation.");
				createQQPlots(permutationDir, nrPermutationsFDR, 100000, fdrcutoff, f, m, pValueRealData.toArray(), significantPvalue, nrSignificantEQTLs, fileName);
			} else {
				createQQPlots(permutationDir, nrPermutationsFDR, maxNrMostSignificantEQTLs, fdrcutoff, f, m, pValueRealData.toArray(), significantPvalue, nrSignificantEQTLs, fileName);
			}

		}
	}

	private static void createQQPlots(String permutationDir, int nrPermutationsFDR, int maxNrMostSignificantEQTLs,
									  double fdrcutoff, FileFormat f, FDRMethod m, double[] pValueRealData, ArrayList<Boolean> significantPvalue, int nrSignificantEQTLs, String fileName) throws IOException {
		DoubleMatrix2D permutedPValues;

		if ((nrPermutationsFDR * (long) maxNrMostSignificantEQTLs) < (Integer.MAX_VALUE - 2)) {
			permutedPValues = new DenseDoubleMatrix2D(nrPermutationsFDR, maxNrMostSignificantEQTLs);
		} else {
			permutedPValues = new DenseLargeDoubleMatrix2D(nrPermutationsFDR, maxNrMostSignificantEQTLs);
		}

		int nrEQTLs = -1;
		permutedPValues.assign(1);

		for (int permutationRound = 0; permutationRound < nrPermutationsFDR; permutationRound++) {
			String fileString = permutationDir + "/PermutedEQTLsPermutationRound" + (permutationRound + 1) + ".txt.gz";
			// read the permuted eqtl output
			TextFile gz = new TextFile(fileString, TextFile.R);

			String[] header = gz.readLineElems(TextFile.tab);
			int snpcol = -1;
			int pvalcol = -1;
			int probecol = -1;
			int genecol = -1;

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
			if (f == FileFormat.REDUCED) {
				//PValue  SNP     Probe   Gene
				if (snpcol == -1 || pvalcol == -1 || probecol == -1 && genecol == -1) {
					System.out.println("Column not found in permutation file: " + fileString);
					System.out.println("PValue: " + pvalcol);
					System.out.println("SNP: " + snpcol);
					System.out.println("Probe: " + probecol);
					System.out.println("Gene: " + genecol);
				}
			}
//			String[] data = gz.readLineElemsReturnReference(TextFile.tab);
			int itr = 0;
			String permln = gz.readLine();
			THashSet<String> visitedEffects = new THashSet<String>();
			while (permln != null) {

				if (permln.length() != 0) {
					if (itr > maxNrMostSignificantEQTLs - 1) {
						break;
					} else {
						int filteronColumn;
						String fdrId;
						if (f == FileFormat.REDUCED) {
							if (m == FDRMethod.PROBELEVEL) {
								// fdrId = data[probecol];
								fdrId = Strings.subsplit(permln, Strings.tab, probecol, probecol + 1)[0];
								filteronColumn = probecol;
							} else if (m == FDRMethod.SNPLEVEL) {
//								fdrId = data[snpcol];
								fdrId = Strings.subsplit(permln, Strings.tab, snpcol, snpcol + 1)[0];
								filteronColumn = snpcol;
							} else if (m == FDRMethod.GENELEVEL && header.length > 3) {
//								fdrId = data[genecol];
								fdrId = Strings.subsplit(permln, Strings.tab, genecol, genecol + 1)[0];
								filteronColumn = genecol;
							} else {
								fdrId = Strings.subsplit(permln, Strings.tab, probecol, probecol + 1)[0];
								filteronColumn = probecol;
							}
						} else {
							if (m == FDRMethod.GENELEVEL) {
//								fdrId = data[QTLTextFile.HUGO];
								fdrId = Strings.subsplit(permln, Strings.tab, QTLTextFile.HUGO, QTLTextFile.HUGO + 1)[0];
								filteronColumn = QTLTextFile.HUGO;
							} else if (m == FDRMethod.SNPLEVEL) {
//								fdrId = data[QTLTextFile.SNP];
								fdrId = Strings.subsplit(permln, Strings.tab, QTLTextFile.SNP, QTLTextFile.SNP + 1)[0];
								filteronColumn = QTLTextFile.HUGO;
							} else if (m == FDRMethod.PROBELEVEL) {
//								fdrId = data[4];
								fdrId = Strings.subsplit(permln, Strings.tab, 4, 5)[0];
								filteronColumn = 4;
							} else {
//								fdrId = data[1] + "\t" + data[4];
								fdrId = Strings.subsplit(permln, Strings.tab, QTLTextFile.SNP, QTLTextFile.SNP + 1)[0] + "-" + Strings.subsplit(permln, Strings.tab, 4, 5)[0];
								filteronColumn = 4;
							}
						}


//						if (f == FileFormat.REDUCED) {
//							if (m == FDRMethod.FULL) {
//								fdrId = data[snpcol] + "-" + data[probecol];
//								filteronColumn = probecol;
//							} else if (m == FDRMethod.GENELEVEL && data.length > 3) {
//								fdrId = data[genecol];
//								filteronColumn = genecol;
//							} else {
//								fdrId = data[probecol];
//								filteronColumn = probecol;
//							}
//
//						} else {
//							if (m == FDRMethod.GENELEVEL) {
//								fdrId = data[QTLTextFile.HUGO];
//								filteronColumn = QTLTextFile.HUGO;
//							} else if (m == FDRMethod.PROBELEVEL) {
//								fdrId = data[4];
//								filteronColumn = 4;
//							} else {
//								fdrId = data[1] + "-" + data[4];
//								filteronColumn = 4;
//							}
//						}
//
						// take top effect per gene / probe
						if (Strings.countSeparatorOccurrences(permln, Strings.tab) > filteronColumn) {
							if (!fdrId.equals("-") && !visitedEffects.contains(fdrId)) {
								permutedPValues.setQuick(permutationRound, itr, Double.parseDouble(Strings.subsplit(permln, Strings.tab, 0, 1)[0]));
//                                permutedPValues[permutationRound][itr] = Double.parseDouble(data[0]);
								visitedEffects.add(fdrId);
								if (itr > 0 && permutedPValues.getQuick(permutationRound, (itr - 1)) > permutedPValues.getQuick(permutationRound, itr)) {
									System.err.println("Sorted P-Value list is not perfectly sorted!!!!");
									System.exit(-1);
								}
								itr++;
							}
						} else {
//							System.out.println(Strings.concat(data, Strings.tab));
						}
						permln = gz.readLine();
					}
				}
			}
			gz.close();

			if (nrEQTLs == -1) {
				nrEQTLs = itr;
			}
		}
		boolean[] significant = new boolean[significantPvalue.size()];

		int pos = 0;
		for (Boolean i : significantPvalue) {
			significant[pos] = i;
			pos++;
		}

		QQPlot qq = new QQPlot();
		qq.draw(fileName, fdrcutoff, nrPermutationsFDR, maxNrMostSignificantEQTLs, permutedPValues.toArray(), pValueRealData, significant, nrSignificantEQTLs);
	}


//	public static void calculateFDR(String eQTLTextFileLoc,
//									int nrPermutationsFDR,
//									int maxNrMostSignificantEQTLs,
//									double fdrcutoff,
//									boolean createQQPlot,
//									String outputDir,
//									String permutationDir,
//									FDRMethod fdrType,
//									boolean createLargeFdrFiles) throws IOException {
//		calculateFDRAdvanced(eQTLTextFileLoc, nrPermutationsFDR, maxNrMostSignificantEQTLs, fdrcutoff, createQQPlot, outputDir, permutationDir, fdrType, createLargeFdrFiles, null, null, null);
//	}

	public static void calculateFDRAdvanced(String eQTLTextFileLoc,
											int nrPermutationsFDR,
											int maxNrMostSignificantEQTLs,
											double fdrcutoff,
											boolean createQQPlot,
											String outputDir,
											String permutationDir,
											FDRMethod fdrType,
											boolean createLargeFdrFiles,
											String snpselectionlist,
											String probeselectionlist,
											String snpprobeselectionlist) throws IOException {
		System.out.println("Using advance FDR calculation.\n");
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

		FileFormat format = FileFormat.REDUCED;
		if (nrColsInPermutedFiles > 7) {
			format = FileFormat.LARGE;
		}

		if (fdrType.equals(FDRMethod.ALL)) {
			for (FDRMethod method : FDRMethod.values()) {
				if (!method.equals(FDRMethod.ALL)) {
					runFDRAdvanced(eQTLTextFileLoc, nrPermutationsFDR, maxNrMostSignificantEQTLs, fdrcutoff, format, method, outputDir, permutationDir, createQQPlot, createLargeFdrFiles, snpselectionlist, probeselectionlist, snpprobeselectionlist);
				}
			}
		} else {
			runFDRAdvanced(eQTLTextFileLoc, nrPermutationsFDR, maxNrMostSignificantEQTLs, fdrcutoff, format, fdrType, outputDir, permutationDir, createQQPlot, createLargeFdrFiles, snpselectionlist, probeselectionlist, snpprobeselectionlist);
		}

	}

	private static void runFDRAdvanced(String baseDir, int nrPermutationsFDR, int maxNrMostSignificantEQTLs,
									   double fdrcutoff, FileFormat f, FDRMethod m, String outputDir, String permutationDir, boolean createQQPlot, boolean createLargeFdrFiles, String snpselectionlist, String probeselectionlist, String snpprobeselectionlist) throws IOException {
		//Load permuted data:
		// load values for each permutation round:
		System.out.println("");
		if (m == FDRMethod.GENELEVEL) {
			System.out.println("Performing gene level FDR");
		} else if (m == FDRMethod.SNPLEVEL) {
			System.out.println("Performing SNP level FDR");
		} else if (m == FDRMethod.PROBELEVEL) {
			System.out.println("Performing probe level FDR");
		} else if (m == FDRMethod.FULL) {
			System.out.println("Determining the FDR using all data");
		}

		int selectionCriteria = 0;

		THashSet<String> selectionOfSnps = null;
		if (snpselectionlist != null) {
			System.out.println("Reading: " + snpselectionlist);
			TextFile t = new TextFile(snpselectionlist, TextFile.R);

			selectionOfSnps = new THashSet<String>(100000000, 4f);
			for (String s : t) {
				selectionOfSnps.add(s);
			}

			t.close();
			System.out.println("Size selection list:" + selectionOfSnps.size());
			selectionCriteria++;
		}

		THashSet<String> selectionOfProbes = null;
		if (probeselectionlist != null) {
			System.out.println("Reading: " + probeselectionlist);
			TextFile t = new TextFile(probeselectionlist, TextFile.R);

			selectionOfProbes = new THashSet<String>(100000000, 4f);
			for (String s : t) {
				selectionOfProbes.add(s);
			}

			t.close();
			System.out.println("Size selection list:" + selectionOfProbes.size());
			selectionCriteria++;
		}

		THashSet<String> selectionOfSnpProbes = null;
		if (snpprobeselectionlist != null) {
			System.out.println("Reading: " + snpprobeselectionlist);
			TextFile t = new TextFile(snpprobeselectionlist, TextFile.R);

			selectionOfSnpProbes = new THashSet<String>(100000000, 4f);
			for (String s : t) {
				selectionOfSnpProbes.add(s);
			}
			t.close();
			System.out.println("Size selection list:" + selectionOfSnpProbes.size());
			selectionCriteria++;
		}

		if (selectionCriteria > 1) {
			System.out.println("Error, only one selection criteria at the same time allowed.");
			System.exit(0);
		}

		TDoubleIntHashMap permutedPvalues = new TDoubleIntHashMap(1000000, 2f);

//        ProgressBar pb = new ProgressBar(nrPermutationsFDR, "Reading permuted data:");


		System.out.println("Reading permuted files");

		for (int permutationRound = 0; permutationRound < nrPermutationsFDR; permutationRound++) {
			String fileString = permutationDir + "/PermutedEQTLsPermutationRound" + (permutationRound + 1) + ".txt.gz";
			System.out.println(fileString);
			// read the permuted eqtl output
			TextFile gz = new TextFile(fileString, TextFile.R, 1048576);

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

			THashSet<String> visitedEffects = new THashSet<String>();
			while (data != null) {
				if (data.length != 0) {
					if (itr > maxNrMostSignificantEQTLs - 1) {
						System.out.println("Breaking because: " + itr);
						break;
					} else {
						String fdrId = null;
						if (f == FileFormat.REDUCED) {
							if (m == FDRMethod.PROBELEVEL) {
								fdrId = new String(data[probecol].getBytes("UTF-8")).intern();
							} else if (m == FDRMethod.GENELEVEL && data.length > 3) {
								fdrId = new String(data[genecol].getBytes("UTF-8")).intern();
							} else if (m == FDRMethod.SNPLEVEL) {
								fdrId = new String(data[snpcol].getBytes("UTF-8")).intern();
							}

						} else {
							if (m == FDRMethod.GENELEVEL) {
								fdrId = new String(data[QTLTextFile.HUGO].getBytes("UTF-8")).intern();
							} else if (m == FDRMethod.SNPLEVEL) {
								fdrId = new String(data[QTLTextFile.SNP].getBytes("UTF-8")).intern();
							} else if (m == FDRMethod.PROBELEVEL) {
								fdrId = new String(data[4].getBytes("UTF-8")).intern();
							}
						}

						// take top effect per gene / probe
						if (selectionOfSnpProbes != null && selectionOfSnpProbes.contains(data[snpcol] + "-" + data[probecol])) {
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
						}
						if (selectionOfSnps != null && selectionOfSnps.contains(data[snpcol])) {
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
						}
						if (selectionOfProbes != null && selectionOfProbes.contains(data[probecol])) {
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
						}
						if (selectionOfProbes != null && selectionOfSnps != null && selectionOfSnpProbes != null) {
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
						}

						data = gz.readLineElemsReturnReference(TextFile.tab);
					}

				}
			}
			System.out.println("\tUsed from permutation " + (permutationRound + 1) + " : " + itr + " rows.");
			gz.close();
		}

		double[] uniquePermutedPvalues = permutedPvalues.keys();
		Arrays.sort(uniquePermutedPvalues);

		double[] uniquePermutedPvaluesCounts = new double[uniquePermutedPvalues.length];

		long cummulativeCount = 0;
		double nrPermutationsFDRd = (double) nrPermutationsFDR;
		for (int i = 0; i < uniquePermutedPvalues.length; ++i) {
			cummulativeCount += permutedPvalues.get(uniquePermutedPvalues[i]);
			uniquePermutedPvaluesCounts[i] = cummulativeCount / nrPermutationsFDRd;
		}
		permutedPvalues = null;
		System.out.println("Number of unique permutation p-values: " + uniquePermutedPvalues.length);

		if (outputDir == null) {
			outputDir = baseDir;
		}

		String fileSuffix = "";
		if (m == FDRMethod.GENELEVEL) {
			fileSuffix = "-GeneLevel";
		} else if (m == FDRMethod.SNPLEVEL) {
			fileSuffix = "-SNPLevel";
		} else if (m == FDRMethod.PROBELEVEL) {
			fileSuffix = "-ProbeLevel";
		}

		String outFileName = outputDir + "/eQTLsFDR" + fdrcutoff + fileSuffix + ".txt.gz";
		String outFileNameSnps = outputDir + "/eQTLSNPsFDR" + fdrcutoff + fileSuffix + ".txt.gz";
		String outFileNameProbes = outputDir + "/eQTLProbesFDR" + fdrcutoff + fileSuffix + ".txt.gz";
		String outFileNameAll = outputDir + "/eQTLsFDR" + fileSuffix + ".txt.gz";

		TextFile outputWriterSignificant = new TextFile(outFileName, TextFile.W);
		TextFile outputWriterESNPs = new TextFile(outFileNameSnps, TextFile.W);
		TextFile outputWriterEProbes = new TextFile(outFileNameProbes, TextFile.W);

		TextFile outputWriterAll = null;
		if (createLargeFdrFiles) {
			outputWriterAll = new TextFile(outFileNameAll, TextFile.W);
		}
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

		if (createLargeFdrFiles) {
			outputWriterAll.append(header);
			outputWriterAll.append("\tFDR\n");
		}
		outputWriterEProbes.append(header);
		outputWriterEProbes.append("\tFDR\n");

		outputWriterESNPs.append(header);
		outputWriterESNPs.append("\tFDR\n");

		outputWriterSignificant.append(header);
		outputWriterSignificant.append("\tFDR\n");

		String str = realEQTLs.readLine();


// REAL DATA PROCESSING
		int itr = 0;
		THashSet<String> visitedEffects = new THashSet<String>();
		THashSet<String> visitedSnps = new THashSet<String>();
		THashSet<String> visitedProbes = new THashSet<String>();
		double lastEqtlPvalue = 0;
		boolean foundHigherFDRThanDesiredCutOff = false;

		double currentPvalue = 0;
		ArrayList<String> currentPvalueEqtls = new ArrayList<String>();
		ArrayList<String> currentPvalueEqtlSnps = new ArrayList<String>();
		ArrayList<String> currentPvalueEqtlProbes = new ArrayList<String>();

		TDoubleArrayList pValueRealData = new TDoubleArrayList();
		ArrayList<Boolean> significantPvalue = new ArrayList<Boolean>();
		int lastUsedPermutedPvalueIndex = 0;

		int nrSignificantEQTLs = 0;
		int lnctr = 0;
		System.out.println("Parsing: " + realEQTLs.getFileName());
		while ((str = realEQTLs.readLine()) != null) {
			if (itr > maxNrMostSignificantEQTLs - 1) {
				break;
			} else {

				String fdrId = null;
				String[] data = Strings.tab.split(str);

				if (selectionOfSnps != null && !selectionOfSnps.contains(data[QTLTextFile.SNP])) {
					continue;
				}
				if (selectionOfProbes != null && !selectionOfProbes.contains(data[QTLTextFile.PROBE])) {
					continue;
				}
				if (selectionOfSnpProbes != null && !selectionOfSnpProbes.contains(data[QTLTextFile.SNP] + "-" + data[QTLTextFile.PROBE])) {
					continue;
				}

				if (m == FDRMethod.GENELEVEL) {
					fdrId = data[QTLTextFile.HUGO];
				} else if (m == FDRMethod.PROBELEVEL) {
					fdrId = data[4];
				} else if (m == FDRMethod.SNPLEVEL) {
					fdrId = data[QTLTextFile.SNP];
				}

				double eQtlPvalue = Double.parseDouble(data[0]);

				if (itr > 0 && lastEqtlPvalue > eQtlPvalue) {
					System.err.println("Sorted P-Value list is not perfectly sorted!!!!");
					System.exit(-1);
				}

				if (eQtlPvalue > currentPvalue) {
					//Process old results for current pvalue

					double fdr = 0;
					if (currentPvalue >= uniquePermutedPvalues[0]) {

						while (uniquePermutedPvalues[lastUsedPermutedPvalueIndex + 1] <= currentPvalue && lastUsedPermutedPvalueIndex < uniquePermutedPvalues.length - 2) {
							++lastUsedPermutedPvalueIndex;
						}
						fdr = uniquePermutedPvaluesCounts[lastUsedPermutedPvalueIndex] / itr;

						if (fdr > 1) {
							fdr = 1;
						}

					}

					for (int i = 0; i < currentPvalueEqtls.size(); ++i) {
						String cachedEqtls = currentPvalueEqtls.get(i);
						String cachedEqtlsProbe = currentPvalueEqtlProbes.get(i);
						String cachedEqtlsSnps = currentPvalueEqtlSnps.get(i);

						StringBuilder currentString = new StringBuilder();
						currentString.append(cachedEqtls).append('\t').append(String.valueOf(fdr)).append('\n');
						pValueRealData.add(currentPvalue);

						if (createLargeFdrFiles) {
							outputWriterAll.append(currentString.toString());
						}

						if (fdr <= fdrcutoff) {
							foundHigherFDRThanDesiredCutOff = true;
							if (!visitedProbes.contains(cachedEqtlsProbe)) {
								outputWriterEProbes.append(currentString.toString());
								visitedProbes.add(cachedEqtlsProbe);
							}
							if (!visitedSnps.contains(cachedEqtlsSnps)) {
								outputWriterESNPs.append(currentString.toString());
								visitedSnps.add(cachedEqtlsSnps);

							}

							significantPvalue.add(true);
							outputWriterSignificant.append(currentString.toString());
							++nrSignificantEQTLs;
						} else {
							significantPvalue.add(false);
						}

					}

					//Create new temp list for this pvalue
					currentPvalue = eQtlPvalue;
					currentPvalueEqtls.clear();
					currentPvalueEqtlProbes.clear();
					currentPvalueEqtlSnps.clear();
					currentPvalueEqtls.add(str);
					currentPvalueEqtlProbes.add(data[QTLTextFile.PROBE]);
					currentPvalueEqtlSnps.add(data[QTLTextFile.SNP]);

				} else {
					//add to current pvalue list
					currentPvalueEqtls.add(str);
					currentPvalueEqtlProbes.add(data[QTLTextFile.PROBE]);
					currentPvalueEqtlSnps.add(data[QTLTextFile.SNP]);
				}

				lastEqtlPvalue = eQtlPvalue;

				if (m == FDRMethod.FULL || (!fdrId.equals("-") && !visitedEffects.contains(fdrId))) {
					itr++;
					visitedEffects.add(fdrId);
				}
			}
			if (lnctr % 10000 == 0) {
				System.out.print(lnctr + " lines parsed.\r");
			}
			lnctr++;

		}
		System.out.println();
		System.out.println();
		//Write buffer to files
		double fdr = 0;
		if (currentPvalue >= uniquePermutedPvalues[0]) {

			while (uniquePermutedPvalues[lastUsedPermutedPvalueIndex + 1] <= currentPvalue && lastUsedPermutedPvalueIndex < uniquePermutedPvalues.length - 2) {
				++lastUsedPermutedPvalueIndex;
			}
			fdr = uniquePermutedPvaluesCounts[lastUsedPermutedPvalueIndex] / itr;

			if (fdr > 1) {
				fdr = 1;
			}
		}

		for (int i = 0; i < currentPvalueEqtls.size(); ++i) {
			String cachedEqtls = currentPvalueEqtls.get(i);
			String cachedEqtlsProbe = currentPvalueEqtlProbes.get(i);
			String cachedEqtlsSnps = currentPvalueEqtlSnps.get(i);

			StringBuilder currentString = new StringBuilder();
			currentString.append(cachedEqtls).append('\t').append(String.valueOf(fdr)).append('\n');

			pValueRealData.add(currentPvalue);
			if (createLargeFdrFiles) {
				outputWriterAll.append(currentString.toString());
			}

			if (fdr <= fdrcutoff) {
				foundHigherFDRThanDesiredCutOff = true;
				if (!visitedProbes.contains(cachedEqtlsProbe)) {
					outputWriterEProbes.append(currentString.toString());
					visitedProbes.add(cachedEqtlsProbe);
				}
				if (!visitedSnps.contains(cachedEqtlsSnps)) {
					outputWriterESNPs.append(currentString.toString());
					visitedSnps.add(cachedEqtlsSnps);
				}

				significantPvalue.add(true);
				outputWriterSignificant.append(currentString.toString());
				++nrSignificantEQTLs;
			} else {
				significantPvalue.add(false);
			}
		}

		realEQTLs.close();
		if (createLargeFdrFiles) {
			outputWriterAll.close();
		}
		outputWriterEProbes.close();
		outputWriterESNPs.close();
		outputWriterSignificant.close();

		if (!foundHigherFDRThanDesiredCutOff) {
			System.out.println("Warning: Not enough results stored. Need more results for desired FDR threshold.");
		}
		//System.out.println("");
		System.out.println("Number of significant eQTLs:\t" + nrSignificantEQTLs);
		System.out.println(" - Number of unique SNPs, constituting an eQTL:\t" + visitedSnps.size());
		System.out.println(" - Number of unique probes, constituting an eQTL:\t" + visitedProbes.size());

		if (createQQPlot) {
			System.out.println("Creating QQ plot. This might take a while...");
			String fileName = baseDir + "/eQTLsFDR" + fdrcutoff + fileSuffix + "-QQPlot.pdf";
			if (maxNrMostSignificantEQTLs > pValueRealData.size()) {
				createQQPlots(permutationDir, nrPermutationsFDR, pValueRealData.size(), fdrcutoff, f, m, pValueRealData.toArray(), significantPvalue, nrSignificantEQTLs, fileName);
			} else if (maxNrMostSignificantEQTLs > 100000) {
				System.out.println("Only taking the top 100,000 for QQplot creation.");
				createQQPlots(permutationDir, nrPermutationsFDR, 100000, fdrcutoff, f, m, pValueRealData.toArray(), significantPvalue, nrSignificantEQTLs, fileName);
			} else {
				createQQPlots(permutationDir, nrPermutationsFDR, maxNrMostSignificantEQTLs, fdrcutoff, f, m, pValueRealData.toArray(), significantPvalue, nrSignificantEQTLs, fileName);
			}

		}
	}


}
