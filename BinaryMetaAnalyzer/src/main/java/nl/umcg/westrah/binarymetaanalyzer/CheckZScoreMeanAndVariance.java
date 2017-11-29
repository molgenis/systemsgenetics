package nl.umcg.westrah.binarymetaanalyzer;

import gnu.trove.map.hash.TObjectIntHashMap;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.util.BaseAnnot;
import umcg.genetica.math.stats.Descriptives;
import umcg.genetica.math.stats.ZScores;
import umcg.genetica.text.Strings;
import umcg.genetica.util.Primitives;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;

public class CheckZScoreMeanAndVariance {
	
	private MetaQTL4TraitAnnotation probeAnnotation;
	private BinaryMetaAnalysisDataset[] datasets = new BinaryMetaAnalysisDataset[0];
	private int[][] snpIndex;
	private String[] snpList;
	private final BinaryMetaAnalysisSettings settings;
	private String[] snpChr;
	private int[] snpPositions;
	private Integer[][] probeIndex;
	
	private boolean bufferHasOverFlown;
	private double maxSavedPvalue = -Double.MAX_VALUE;
	private boolean sorted;
	private int locationToStoreResult;
	private MetaQTL4MetaTrait[][] snpprobeCombos;
	
	TObjectIntHashMap<MetaQTL4MetaTrait> traitMap = new TObjectIntHashMap<MetaQTL4MetaTrait>();
	MetaQTL4MetaTrait[] traitList = null;
	
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
