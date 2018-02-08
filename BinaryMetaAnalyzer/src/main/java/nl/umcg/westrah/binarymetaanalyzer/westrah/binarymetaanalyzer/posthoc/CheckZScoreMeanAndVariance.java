package nl.umcg.westrah.binarymetaanalyzer.westrah.binarymetaanalyzer.posthoc;

import nl.umcg.westrah.binarymetaanalyzer.BinaryMetaAnalysis;
import nl.umcg.westrah.binarymetaanalyzer.BinaryMetaAnalysisDataset;
import nl.umcg.westrah.binarymetaanalyzer.BinaryMetaAnalysisSettings;
import nl.umcg.westrah.binarymetaanalyzer.MetaQTL4TraitAnnotation;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.util.Primitives;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;

public class CheckZScoreMeanAndVariance {
	
	private MetaQTL4TraitAnnotation probeAnnotation;
	
	private BinaryMetaAnalysisSettings settings;
	
	
	public static void main(String[] args) {
		String in = "C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\eqtl\\met\\";
		String out = "C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\eqtl\\met\\check\\";
		int nrperm = 1;
		CheckZScoreMeanAndVariance c = new CheckZScoreMeanAndVariance();
		try {
			c.checkZScoreTable(in, out, nrperm);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public CheckZScoreMeanAndVariance(String settingsFile, String textToReplace, String replaceTextWith) {
		// initialize settings
		settings = new BinaryMetaAnalysisSettings();
		settings.parse(settingsFile, textToReplace, replaceTextWith);
		
		try {
			runPerDataset();
		} catch (IOException ex) {
			ex.printStackTrace();
			Logger.getLogger(BinaryMetaAnalysis.class.getName()).log(Level.SEVERE, null, ex);
		}
		
	}
	
	public CheckZScoreMeanAndVariance() {
	
	}
	
	public void checkZScoreTable(String inLoc, String outLoc, int nrPerm) throws IOException {
		
		double[][] zMeans = new double[nrPerm][];
		double[][] zVars = new double[nrPerm][];
		double[][] zNs = new double[nrPerm][];
		int maxNrProbes = 0;
		
		ArrayList<ArrayList<String>> probeListPerFile = new ArrayList<ArrayList<String>>();
		for (int permutation = 0; permutation < nrPerm; permutation++) {
//			header += "\tNZ\tMeanZ\tVarZ";
			
			
			// perm0\tperm1\tperm2
			// gene0
			// gene1
			
			String zmatfileloc = inLoc + "ZScoreMatrix-Permutation" + permutation + ".txt.gz";
			if (permutation == 0) {
				zmatfileloc = inLoc + "ZScoreMatrix.txt.gz";
//				zmatfileloc = inLoc + "ZScoreMatrix2.txt";
			}
			System.out.println(zmatfileloc + " parsing");
			TextFile tf = new TextFile(zmatfileloc, TextFile.R);
			String[] headerIn = tf.readLineElems(TextFile.tab);
			String ln = tf.readLine();
			int nrsnps = 0;
			while (ln != null) {
				nrsnps++;
				if (nrsnps % 10 == 0) {
					System.out.print("\r" + nrsnps + " lines read");
				}
				ln = tf.readLine();
			}
			System.out.println();
			tf.close();
			tf.open();
			tf.readLine();
			
			// format: SNP     Alleles AlleleAssessed  66023_ENST00000294984
			int nrprobes = headerIn.length - 3;
			ArrayList<String> fileprobes = new ArrayList<String>();
			for (int g = 3; g < headerIn.length; g++) {
				fileprobes.add(headerIn[g]);
			}
			probeListPerFile.add(fileprobes);
			double[][] zmat = new double[nrprobes][nrsnps];
			System.out.println("size: " + nrprobes + " genes \t" + nrsnps + " snps");
			String[] elems = tf.readLineElems(TextFile.tab);
			int snp = 0;
			while (elems != null) {
				
				double[] probes = new double[elems.length - 3];
				for (int i = 3; i < elems.length; i++) {
					probes[i - 3] = Double.parseDouble(elems[i]);
				}
				
				for (int i = 0; i < probes.length; i++) {
					zmat[i][snp] = probes[i];
				}
				
				snp++;
				if (snp % 100 == 0) {
					System.out.print("\r" + snp + " out of " + nrsnps);
				}
				elems = tf.readLineElems(TextFile.tab);
			}
			tf.close();
			
			System.out.println();
			
			double[] zmean = new double[nrprobes];
			double[] zvar = new double[nrprobes];
			double[] zn = new double[nrprobes];
			int nrWithOutVals = 0;
			for (int p = 0; p < nrprobes; p++) {
				ArrayList<Double> nonan = new ArrayList<Double>();
				for (int s = 0; s < zmat[p].length; s++) {
					if (!Double.isNaN(zmat[p][s])) {
						nonan.add(zmat[p][s]);
					}
				}
				if (nonan.size() > 0) {
					double[] arr = Primitives.toPrimitiveArr(nonan.toArray(new Double[0]));
					zmean[p] = JSci.maths.ArrayMath.mean(arr);
					zvar[p] = JSci.maths.ArrayMath.variance(arr);
					zn[p] = nonan.size();
					
					if (zn.length > maxNrProbes) {
						maxNrProbes = zn.length;
					}
				} else {
					nrWithOutVals++;
					// System.out.println(p + " has length 0 " + zmat[p].length + " vals before filter\t" + nonan.size() + "\tafter filter");
				}
			}
			System.out.println(nrWithOutVals + " probes have no vals out of " + nrprobes);
			
			zMeans[permutation] = zmean;
			zVars[permutation] = zvar;
			zNs[permutation] = zn;
			
		}
		
		// write to disk
		String header = "ProbeRank";
		for (int permutation = 0; permutation <= nrPerm; permutation++) {
			header += "\tGene-Perm\tN-Perm\tMeanZ-Perm" + permutation + "\tVarZ-Perm" + permutation;
		}
		
		TextFile out = new TextFile(outLoc + "ZScoreMeanAndVariancePerProbe.txt", TextFile.W);
		System.out.println("Writing to: " + outLoc + "ZScoreMeanAndVariancePerProbe.txt");
		out.writeln(header);
		for (int i = 0; i < maxNrProbes; i++) {
			String ln = "" + i;
			for (int p = 0; p < zMeans.length; p++) {
				if (i >= zMeans[p].length) {
					ln += "\t-\t0\tNotTested\tNotTested";
				} else {
					ArrayList<String> geneNames = probeListPerFile.get(p);
					String gene = geneNames.get(i);
					ln += "\t" + gene + "\t" + zNs[p][i] + "\t" + zMeans[p][i] + "\t" + zVars[p][i];
				}
			}
			out.writeln(ln);
		}
		out.close();
		
	}
	
	public void runPerDataset() throws IOException {
		
		String outdir = settings.getOutput();
		System.out.println("Placing output here: " + outdir);
		outdir = Gpio.formatAsDirectory(outdir);
		Gpio.createDir(outdir);
		// load probe annotation and index
		// this particular probe annotation can take multiple probes for a single location into account.
		System.out.println("Loading probe annotation from: " + settings.getProbetranslationfile());
		
		HashSet<String> platforms = new HashSet<String>();
		platforms.addAll(settings.getDatasetannotations());
		System.out.println("Defined platforms in settings file: ");
		for (String s : platforms) {
			System.out.println(s);
		}
		probeAnnotation = new MetaQTL4TraitAnnotation(new File(settings.getProbetranslationfile()), platforms);
		
		
		System.out.println("Permutations: " + settings.getStartPermutations() + " until " + settings.getNrPermutations());
		
		
		int nrDataset = settings.getDatasetlocations().size();
		for (int d = 0; d < nrDataset; d++) {
			// map genes to Z and Var
			// double[] z = new double[nrPerm]
			double[][] zMeans = new double[settings.getNrPermutations() + 1][];
			double[][] zVars = new double[settings.getNrPermutations() + 1][];
			double[][] zN = new double[settings.getNrPermutations() + 1][];
			int maxNrProbes = 0;
			for (int permutation = settings.getStartPermutations(); permutation <= settings.getNrPermutations(); permutation++) {
				BinaryMetaAnalysisDataset dataset = new BinaryMetaAnalysisDataset(settings.getDatasetlocations().get(d),
						settings.getDatasetnames().get(d),
						settings.getDatasetPrefix().get(d),
						permutation,
						settings.getDatasetannotations().get(d),
						probeAnnotation);
				
				String[] probeList = dataset.getProbeList();
				String[] snpList = dataset.getSNPs();
				double[][] zmat = new double[probeList.length][snpList.length];
				
				System.out.println("Loading z-scores");
				for (int i = 0; i < snpList.length; i++) {
					float[] z = dataset.getZScores(i);
					for (int p = 0; p < probeList.length; p++) {
						zmat[p][i] = (double) z[p];
					}
					System.out.print("\rPerm " + permutation + "\tsnp " + i + "/" + snpList.length);
				}
				
				System.out.println();
				
				double[] zmean = new double[probeList.length];
				double[] zvar = new double[probeList.length];
				double[] zn = new double[probeList.length];
				
				for (int p = 0; p < probeList.length; p++) {
					// remove nans
					ArrayList<Double> nonan = new ArrayList<Double>();
					for (double z : zmat[p]) {
						if (!Double.isNaN(z)) {
							nonan.add(z);
						}
					}
					double[] arr = Primitives.toPrimitiveArr(nonan.toArray(new Double[0]));
					zmean[p] = JSci.maths.ArrayMath.mean(arr);
					zvar[p] = JSci.maths.ArrayMath.variance(arr);
					zn[p] = nonan.size();
				}
				
				if (probeList.length > maxNrProbes) {
					maxNrProbes = probeList.length;
				}
				zN[permutation] = zn;
				zMeans[permutation] = zmean;
				zVars[permutation] = zvar;
			}
			
			
			System.out.println();
			
			// write to disk
			String header = "ProbeRank";
			for (int permutation = settings.getStartPermutations(); permutation <= settings.getNrPermutations(); permutation++) {
				header += "\tN-Perm\tMeanZ-Perm" + permutation + "\tVarZ-Perm" + permutation;
			}
			
			String datasetname = settings.getDatasetnames().get(d);
			
			TextFile out = new TextFile(settings.getOutput() + datasetname + "-ZScoreMeanAndVariancePerProbe.txt", TextFile.W);
			System.out.println("Writing to: " + settings.getOutput() + datasetname + "-ZScoreMeanAndVariancePerProbe.txt");
			out.writeln(header);
			for (int i = 0; i < maxNrProbes; i++) {
				String ln = "" + i;
				for (int p = 0; p < zMeans.length; p++) {
					if (i >= zMeans[p].length) {
						ln += "\t0\tNotTested\tNotTested";
					} else {
						ln += "\t" + zN[p][i] + "\t" + zMeans[p][i] + "\t" + zVars[p][i];
					}
				}
				out.writeln(ln);
			}
			out.close();
		}
		
	}
	
}
