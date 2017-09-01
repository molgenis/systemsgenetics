package eqtlmappingpipeline.normalization;

import org.apache.commons.math3.stat.ranking.NaNStrategy;
import org.apache.commons.math3.stat.ranking.NaturalRanking;
import org.apache.commons.math3.stat.ranking.TiesStrategy;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.PCA;
import umcg.genetica.math.matrix.DoubleMatrixDataset;
import umcg.genetica.math.matrix.MatrixHandling;
import umcg.genetica.math.matrix.MatrixTools;
import umcg.genetica.math.stats.Descriptives;
import umcg.genetica.math.stats.Log2Transform;
import umcg.genetica.math.stats.QuantileNormalization;
import umcg.genetica.math.stats.Regression;
import umcg.genetica.math.stats.concurrent.ConcurrentCorrelation;
import umcg.genetica.math.stats.concurrent.ConcurrentCovariation;
import umcg.genetica.methylation.ConvertBetaAndMvalues;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

/**
 * @author harmjan
 */
public class Normalizer {

	//nrIntermediatePCAsOverSamplesToRemoveToOutput = 5
	//nrPCAsOverSamplesToRemove = 100
	public void normalize(String expressionFile, String probeIncludeList, String sampleIncludeList, int nrPCAsOverSamplesToRemove, int nrIntermediatePCAsOverSamplesToRemoveToOutput, String covariatesToRemove, boolean orthogonalizecovariates, String outdir,
						  boolean runQQNorm, boolean runLog2Transform, boolean runMTransform, boolean runCenterScale, boolean runPCA, boolean adjustCovariates, boolean forceMissingValues, boolean forceReplacementOfMissingValues,
						  boolean forceReplacementOfMissingValues2, boolean treatZerosAsNulls, boolean forceNormalDistribution) throws IOException {

		System.out.println("Running normalization.");
		if (outdir != null) {
			outdir = Gpio.formatAsDirectory(outdir);
			Gpio.createDir(outdir);
		} else {
			if (Gpio.getParentDir(expressionFile) == null) {
				//This happens for relative paths in current dir
				outdir = "";
			} else {
				outdir = Gpio.getParentDir(expressionFile) + Gpio.getFileSeparator();
			}

		}

		String parentDir = Gpio.getParentDir(expressionFile);
		String expressionFileName = Gpio.getFileName(expressionFile);
		if (parentDir == null) {
			parentDir = "";
		}

		if (expressionFileName.contains(".txt.gz")) {
			expressionFileName = expressionFileName.replaceAll(".txt.gz", "");
		} else {
			expressionFileName = expressionFileName.replaceAll(".txt", "");
		}

		String outputFileNamePrefix = outdir + expressionFileName;


		Set<String> s = null;
		if (sampleIncludeList != null) {
			TextFile t = new TextFile(sampleIncludeList, TextFile.R);
			s = new HashSet<String>(t.readAsArrayList());
		}
		Set<String> p = null;
		if (probeIncludeList != null) {
			TextFile t = new TextFile(probeIncludeList, TextFile.R);
			p = new HashSet<String>(t.readAsArrayList());
		}
		DoubleMatrixDataset<String, String> dataset = null;

		if (s != null || p != null) {
			dataset = new DoubleMatrixDataset<String, String>(expressionFile, p, s);
			//Check if samples are correclty loaded.
			boolean breakAfterCheck = false;
			if (s != null) {
				outputFileNamePrefix = outputFileNamePrefix + ".SampleSelection";
				HashSet<String> tmpNames = new HashSet<String>();
				tmpNames.addAll(dataset.colObjects);
				tmpNames.addAll(s);
				HashSet<String> missingNames = new HashSet<String>();
				HashSet<String> extraNames = new HashSet<String>();
				for (String colName : tmpNames) {
					if (!s.contains(colName)) {
						extraNames.add(colName);
					}
					if (!dataset.colObjects.contains(colName)) {
						missingNames.add(colName);
					}
				}
				if (!missingNames.isEmpty()) {
					System.err.println("\nMatrix does not contains desired columns, please check filtering list.");
					System.err.println(missingNames.toString() + "\n");
					breakAfterCheck = true;
				} else if (!extraNames.isEmpty()) {
					System.err.println("\nMatrix contains unwanted columns, please check filtering list.");
					System.err.println(extraNames.toString() + "\n");
					breakAfterCheck = true;
				}
			}
			//Check if probes are correclty loaded.
			if (p != null) {
				outputFileNamePrefix = outputFileNamePrefix + ".ProbeSelection";
				HashSet<String> tmpNames = new HashSet<String>();
				tmpNames.addAll(dataset.rowObjects);
				tmpNames.addAll(p);
				HashSet<String> missingNames = new HashSet<String>();
				HashSet<String> extraNames = new HashSet<String>();
				for (String rowName : tmpNames) {
					if (!p.contains(rowName)) {
						extraNames.add(rowName);
					}
					if (!dataset.rowObjects.contains(rowName)) {
						missingNames.add(rowName);
					}
				}
				if (!missingNames.isEmpty()) {
					System.err.println("\nMatrix does not contains desired rows, please check filtering list.");
					System.err.println(missingNames.toString() + "\n");
					breakAfterCheck = true;
				} else if (!extraNames.isEmpty()) {
					System.err.println("\nMatrix contains unwanted rows, please check filtering list.");
					System.err.println(extraNames.toString() + "\n");
					breakAfterCheck = true;
				}
			}

//            if(breakAfterCheck){
//                System.exit(-1);
//            }

			dataset.save(outputFileNamePrefix + ".txt.gz");
		} else {
			dataset = new DoubleMatrixDataset<String, String>(expressionFile);
		}


		// check for probes with zero variance, if there > 3 samples in the dataset
		if (dataset.nrCols > 3) {
			outputFileNamePrefix = removeProbesWithZeroVariance(dataset, outputFileNamePrefix);
		}

		if (runQQNorm) {
			outputFileNamePrefix = quantileNormalize(dataset, outputFileNamePrefix, forceMissingValues, forceReplacementOfMissingValues, forceReplacementOfMissingValues2, treatZerosAsNulls);
		}
		if (runLog2Transform) {
			outputFileNamePrefix = log2transform(dataset, outputFileNamePrefix);
		}
		if (runMTransform) {
			outputFileNamePrefix = mValueTransform(dataset, outputFileNamePrefix);
		}
		if (runCenterScale) {
			outputFileNamePrefix = centerAndScale(dataset, outputFileNamePrefix);
		}

		if (adjustCovariates && covariatesToRemove != null) {
			outputFileNamePrefix = adjustCovariates(dataset, outputFileNamePrefix, covariatesToRemove, orthogonalizecovariates, 1E-10);
		}

		if (runPCA) {
			ConcurrentCorrelation c = new ConcurrentCorrelation(2);
			double[][] correlationMatrix = c.pairwiseCorrelation(dataset.getRawDataTransposed());
			Pair<DoubleMatrixDataset<String, String>, DoubleMatrixDataset<String, String>> PCAResults = calculatePCA(dataset, correlationMatrix, outputFileNamePrefix, null);
			if (nrPCAsOverSamplesToRemove != 0 || nrIntermediatePCAsOverSamplesToRemoveToOutput != 0) {
				correctDataForPCs(dataset, outputFileNamePrefix, nrPCAsOverSamplesToRemove, nrIntermediatePCAsOverSamplesToRemoveToOutput, PCAResults.getLeft(), PCAResults.getRight());
			}
		}

		if (forceNormalDistribution) {
			outputFileNamePrefix = forceNormalDistribution(dataset, outputFileNamePrefix);
		}
	}


	NaturalRanking ranking = new NaturalRanking(NaNStrategy.FAILED, TiesStrategy.AVERAGE);

    
    public Pair<String,String> calculatePcaOnly(String expressionFile) throws IOException{
        String outdir = Gpio.getParentDir(expressionFile) + Gpio.getFileSeparator();
		
		String parentDir = Gpio.getParentDir(expressionFile);
		String expressionFileName = Gpio.getFileName(expressionFile);
		if (parentDir == null) {
			parentDir = "";
		}

		if (expressionFileName.contains(".txt.gz")) {
			expressionFileName = expressionFileName.replaceAll(".txt.gz", "");
		} else {
			expressionFileName = expressionFileName.replaceAll(".txt", "");
		}
        
        DoubleMatrixDataset<String, String> dataset = new DoubleMatrixDataset<String, String>(expressionFile);
        
		String outputFileNamePrefix = outdir + expressionFileName;
        
        ConcurrentCorrelation c = new ConcurrentCorrelation(2);
        double[][] correlationMatrix = c.pairwiseCorrelation(dataset.getRawDataTransposed());
        calculatePCA(dataset, correlationMatrix, outputFileNamePrefix, null);
        return new Pair<String,String>(outputFileNamePrefix + ".PCAOverSamplesEigenvectorsTransposed.txt.gz",outputFileNamePrefix + ".PCAOverSamplesPrincipalComponents.txt.gz");
    }
    
    
    
	public double[] forceNormal(double[] data) {
		double[] rankedValues = ranking.rank(data);
		for (int s = 0; s < data.length; s++) {
			//Convert the rank to a proportion, with range <0, 1>
			double pValue = (0.5d + rankedValues[s] - 1d) / (double) (rankedValues.length);
			//Convert the pValue to a Z-Score:
			data[s] = cern.jet.stat.Probability.normalInverse(pValue);
		}
		return data;
	}


	public String forceNormalDistribution(DoubleMatrixDataset<String, String> dataset, String fileNamePrefix) throws IOException {
		double[][] rawData = dataset.getRawData();
		for (int p = 0; p < dataset.rowObjects.size(); p++) {
			rawData[p] = forceNormal(rawData[p]);
		}

		DoubleMatrixDataset<String, String> datasetNormalized = new DoubleMatrixDataset<String, String>(rawData, dataset.rowObjects, dataset.colObjects);
		fileNamePrefix += ".ForcedNormal";
		datasetNormalized.save(fileNamePrefix + ".txt.gz");
		return fileNamePrefix;


	}

	public String quantileNormalize(DoubleMatrixDataset<String, String> dataset, String fileNamePrefix, boolean forceMissingValues, boolean forceReplacementOfMissingValues, boolean forceReplacementOfMissingValues2, boolean treatZerosAsNulls) throws IOException {
		double[][] rawData = dataset.getRawData();

		boolean dataContainsNulls = MatrixTools.containsNaNs(rawData);

		if (treatZerosAsNulls && dataContainsNulls) {
			System.out.println("Warning: Data already contains nulls before treating zeros as nulls.\n Later on it will not be possible to distinguish between those two!");
		}
		if (treatZerosAsNulls) {
			MatrixHandling.ReplaceZerosToNull(rawData);
			dataContainsNulls = MatrixTools.containsNaNs(rawData);
		}

		if (!dataContainsNulls) {
			QuantileNormalization.quantilenormalize(rawData);
		} else if (forceReplacementOfMissingValues) {
			QuantileNormalization.QuantileNormAdressingNaValuesAfterInitialQN(dataset, false, false, false);
		} else if (forceReplacementOfMissingValues2) {
			QuantileNormalization.QuantileNormAdressingNaValuesAfterInitialQN(dataset, false, true, false);
		} else if (forceMissingValues && treatZerosAsNulls) {
			QuantileNormalization.QuantileNormAdressingNaValuesAfterInitialQN(dataset, true, false, true);
		} else if (forceMissingValues) {
			QuantileNormalization.QuantileNormAdressingNaValuesAfterInitialQN(dataset, true, false, false);
		} else {
			System.out.println("Warning: Your data contains missing values and missing value treatment is not selected.\n"
					+ "If desired please supply additional flag: --forceMissingValues or --forceReplacementOfMissingValues");
			System.exit(0);
		}

		if (treatZerosAsNulls) {
			MatrixHandling.ReplaceNullToZero(rawData);
		}

		DoubleMatrixDataset<String, String> datasetNormalized = new DoubleMatrixDataset<String, String>(rawData, dataset.rowObjects, dataset.colObjects);
		fileNamePrefix += ".QuantileNormalized";
		datasetNormalized.save(fileNamePrefix + ".txt.gz");

		return fileNamePrefix;
	}

	public String log2transform(DoubleMatrixDataset<String, String> dataset, String fileNamePrefix) throws IOException {
		double[][] rawData = dataset.getRawData();
		Log2Transform.log2transform(rawData);
		DoubleMatrixDataset<String, String> datasetNormalized = new DoubleMatrixDataset<String, String>(rawData, dataset.rowObjects, dataset.colObjects);
		fileNamePrefix += ".Log2Transformed";
		datasetNormalized.save(fileNamePrefix + ".txt.gz");
		return fileNamePrefix;
	}

	public String mValueTransform(DoubleMatrixDataset<String, String> dataset, String fileNamePrefix) throws IOException {
		double[][] rawData = dataset.getRawData();
		ConvertBetaAndMvalues.transformToMvalue(rawData);
		DoubleMatrixDataset<String, String> datasetNormalized = new DoubleMatrixDataset<String, String>(rawData, dataset.rowObjects, dataset.colObjects);
		fileNamePrefix += ".MvalueTransformed";
		datasetNormalized.save(fileNamePrefix + ".txt.gz");
		return fileNamePrefix;
	}

	public String centerAndScale(DoubleMatrixDataset<String, String> dataset, String fileNamePrefix) throws IOException {
		double[][] rawData = dataset.getRawData();
		System.out.println("Standardizing probe mean");
		for (int p = 0; p < dataset.rowObjects.size(); p++) {
			double mean = Descriptives.mean(rawData[p]);
			//double stdev = Math.sqrt(Descriptives.variance(rawData[p], mean));
			for (int s = 0; s < dataset.colObjects.size(); s++) {
				rawData[p][s] -= mean;
			}
		}

		dataset.setRawData(rawData);
		fileNamePrefix += ".ProbesCentered";
		dataset.save(fileNamePrefix + ".txt.gz");

		System.out.println("- Standardizing sample mean and standard deviation");
		for (int s = 0; s < dataset.colObjects.size(); s++) {
			double[] vals = new double[dataset.rowObjects.size()];
			for (int p = 0; p < dataset.rowObjects.size(); p++) {
				vals[p] = dataset.getRawData()[p][s];
			}
			double mean = Descriptives.mean(vals);
			for (int p = 0; p < dataset.rowObjects.size(); p++) {
				vals[p] -= mean;
			}
			double var = Descriptives.variance(vals, mean);
			double stdev = Math.sqrt(var);
			for (int p = 0; p < dataset.rowObjects.size(); p++) {
				dataset.getRawData()[p][s] = (vals[p] / stdev);
			}
		}

		DoubleMatrixDataset<String, String> datasetNormalized = new DoubleMatrixDataset<String, String>(rawData, dataset.rowObjects, dataset.colObjects);
		fileNamePrefix += ".SamplesZTransformed";
		datasetNormalized.save(fileNamePrefix + ".txt.gz");
		return fileNamePrefix;
	}

	public String adjustCovariates(DoubleMatrixDataset<String, String> traitData, String fileNamePrefix, String covariatesToRemove, boolean orthogonalizecovariates, double varianceExplainedCutoff) throws IOException {
		// load covariate data, and remove samples for which there is missing covariate data.
		Pair<DoubleMatrixDataset<String, String>, DoubleMatrixDataset<String, String>> covariateData = loadCovariateValues(covariatesToRemove, traitData);
		DoubleMatrixDataset<String, String> covariateDataset = covariateData.getLeft();
		DoubleMatrixDataset<String, String> traitDataUpdated = covariateData.getRight();

		traitData.rawData = traitDataUpdated.rawData;
		traitData.colObjects = traitDataUpdated.colObjects;
		traitData.rowObjects = traitDataUpdated.rowObjects;
		traitData.recalculateHashMaps();

		double[][] covariateValues = null;
		double[] pcaExpVar = null;

		System.out.println("Covariate data has " + covariateDataset.nrRows + " rows and " + covariateDataset.nrCols + " columns.");

		for (int p = 0; p < covariateDataset.rowObjects.size(); p++) {
			double mean = Descriptives.mean(covariateDataset.getRawData()[p]);
			double stdev = Math.sqrt(Descriptives.variance(covariateDataset.getRawData()[p], mean));
			for (int s = 0; s < covariateDataset.colObjects.size(); s++) {
				covariateDataset.getRawData()[p][s] -= mean;
				covariateDataset.getRawData()[p][s] /= stdev;
			}
		}

		//Covariation on a centered and scaled matrix equals the correlation.
		//Covariation is faster to compute.
		ConcurrentCovariation c = new ConcurrentCovariation(2);
		double[][] correlationMatrix = c.pairwiseCovariation(covariateDataset.getRawData());
		covariateDataset.transposeDataset();
		Pair<DoubleMatrixDataset<String, String>, DoubleMatrixDataset<String, String>> PCAResults = calculatePCA(covariateDataset, correlationMatrix, covariatesToRemove, null);

		// replace covariateValues with orthogonal ones...
		covariateDataset = PCAResults.getLeft();


		covariateDataset.transposeDataset();
		covariateValues = covariateDataset.getRawData();

		System.out.println(covariateDataset.nrRows + " covariates finally loaded.");

		// load the eigenvalues
		pcaExpVar = new double[covariateValues.length];
		System.out.println("Loading eigenvalues from: " + covariatesToRemove + ".PCAOverSamplesEigenvalues.txt.gz");
		TextFile tf = new TextFile(covariatesToRemove + ".PCAOverSamplesEigenvalues.txt.gz", TextFile.R); //
		// skip header
		tf.readLine();
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			if (elems.length > 2) {
				int pcanr = Integer.parseInt(elems[0]);
				double expvar = Double.parseDouble(elems[1]);
				pcaExpVar[pcanr - 1] = expvar;
				System.out.println(pcanr + "\t" + expvar);
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
//        } else {
//            // PCA has been performed a-priori. Just check whether the user has supplied proper covariates.
//            if (covariateValues.length > 1) {
//                // check whether the covariates are orthogonal, by calculating the sum of products (inner product)
//                System.out.println("Determining whether covariates are orthogonal, since you defined > 1 covariate:");
//                System.out.println("Covariate1\tCovariate2\tInnerProduct\tCorrelation");
//                double dotproductthreshold = 1E-5;
//                if (covariateValues.length < 100) {
//                    dotproductthreshold = 0.05;
//                }
//                for (int i = 0; i < covariateValues.length; i++) {
//
//                    for (int j = i + 1; j < covariateValues.length; j++) {
//                        double dotproduct = 0;
//                        for (int v = 0; v < covariateValues[i].length; v++) {
//                            dotproduct += covariateValues[i][v] * covariateValues[j][v];
//                        }
//                        double corr = JSci.maths.ArrayMath.correlation(covariateValues[i], covariateValues[j]);
//
//                        if (Math.abs(dotproduct) > dotproductthreshold) {
//                            System.out.println("Innerproduct > 1E-5 for covariates " + covariateDataset.rowObjects.get(i) + " and " + covariateDataset.rowObjects.get(j) + ", InnerProduct: " + Math.abs(dotproduct) + "\tCorrelation: " + corr);
//                            System.out.println("If you want, we can orthogonalize the covariates for you: use --covpca in your command line.");
//                            System.exit(0);
//                        }
//
//                        System.out.println(covariateDataset.rowObjects.get(i) + "\t" + covariateDataset.rowObjects.get(j) + "\t" + dotproduct + "\t" + corr);
//                    }
//
//                }
//
//                System.out.println("Covariates are orthogonal. Now adjusting for covariates.");
//            }
//        }


		double[][] rawdata = traitData.getRawData();
		for (int i = 0; i < covariateValues.length; i++) {
			if (pcaExpVar == null || pcaExpVar[i] > varianceExplainedCutoff) {
				correctForCovariate(rawdata, covariateValues, i);
			} else {
				System.out.println("Not regressing covariate: " + i + " because explained variance < " + varianceExplainedCutoff + ": " + pcaExpVar[i]);
			}
		}

		traitData.rawData = rawdata;
		
		//Why was this done???????
		//DoubleMatrixDataset<String, String> datasetNormalized = new DoubleMatrixDataset<String, String>(rawdata, traitData.rowObjects, traitData.colObjects);
		fileNamePrefix += ".CovariatesRemoved";
		traitData.save(fileNamePrefix + ".txt.gz");

		

		return fileNamePrefix;
	}

	/**
	 * Calculate correlation over columns in DoubleMatrixDataset. WARNING: this
	 * method assumes that SD == 1 and mean == 0 (which makes the covariance
	 * equal to the correlation).
	 *
	 * @param dataset
	 * @return
	 */
	private double[][] correlateSamples(DoubleMatrixDataset<String, String> dataset) {
		double[][] correlationMatrix = new double[dataset.colObjects.size()][dataset.colObjects.size()];
		double probeCountMinusOne = dataset.rowObjects.size() - 1;

		ProgressBar pb = new ProgressBar(dataset.colObjects.size(), "- Calculating correlations: " + dataset.colObjects.size() + " x " + dataset.colObjects.size());

		for (int f = 0; f < dataset.colObjects.size(); f++) {


			for (int g = f; g < dataset.colObjects.size(); g++) {
				double covarianceInterim = 0;
				for (int p = 0; p < dataset.rowObjects.size(); p++) {
					covarianceInterim += dataset.getRawData()[p][f] * dataset.getRawData()[p][g];
				}
				double covariance = covarianceInterim / probeCountMinusOne;
				correlationMatrix[f][g] = covariance;
				correlationMatrix[g][f] = covariance;
//                System.out.println(f + "\t" + g + "\t" + covariance);
			}
			pb.iterate();
		}
		pb.close();
		return correlationMatrix;
	}

	public double[][] correlateProbes(DoubleMatrixDataset<String, String> dataset) {

		double[][] correlationMatrix = new double[dataset.rowObjects.size()][dataset.rowObjects.size()];
		double probeCountMinusOne = dataset.rowObjects.size() - 1;

		ProgressBar pb = new ProgressBar(dataset.rowObjects.size(), "- Calculating correlations: " + dataset.rowObjects.size() + " x " + dataset.rowObjects.size());
		for (int f = 0; f < dataset.rowObjects.size(); f++) {
			for (int g = f; g < dataset.rowObjects.size(); g++) {
				double covarianceInterim = 0;
				for (int p = 0; p < dataset.rowObjects.size(); p++) {
					covarianceInterim += dataset.getRawData()[p][f] * dataset.getRawData()[p][g];
				}
				double covariance = covarianceInterim / probeCountMinusOne;
				correlationMatrix[f][g] = covariance;
				correlationMatrix[g][f] = covariance;
				System.out.println(f + "\t" + g + "\t" + covariance);
			}
			pb.iterate();
		}
		pb.close();
		return correlationMatrix;
	}

	public Pair<DoubleMatrixDataset<String, String>, DoubleMatrixDataset<String, String>> calculatePCA(DoubleMatrixDataset<String, String> dataset, double[][] correlationMatrix, String fileNamePrefix, Integer nrOfPCsToCalculate) throws IOException {
		String expressionFile = fileNamePrefix;
		System.out.println("Calculating PCA over file: " + fileNamePrefix);
		System.out.println("- Performing PCA over correlation matrix of size: " + correlationMatrix.length + "x" + correlationMatrix.length);
		Jama.EigenvalueDecomposition eig = PCA.eigenValueDecomposition(correlationMatrix);

		if (nrOfPCsToCalculate == null || nrOfPCsToCalculate > dataset.colObjects.size()) {
			nrOfPCsToCalculate = dataset.colObjects.size();
		} else if (nrOfPCsToCalculate < 1) {
			throw new IllegalArgumentException("Number of PCs to calculate should be at least 1");
		}

		DoubleMatrixDataset<String, String> datasetEV = new DoubleMatrixDataset<String, String>(dataset.colObjects.size(), nrOfPCsToCalculate);
		datasetEV.rowObjects = dataset.colObjects;
		double[] eigenValues = eig.getRealEigenvalues();
		System.out.println("Eigenvalue results:");

		System.out.println("PCA\tPCANr\tEigenValue\tExplainedVariance\tTotalExplainedVariance");

		TextFile out = new TextFile(expressionFile + ".PCAOverSamplesEigenvalues.txt.gz", TextFile.W);
		double cumExpVarPCA = 0;

		out.writeln("PCA\tPCANr\tEigenValue\tExplainedVariance\tTotalExplainedVariance");

		for (int pca = 0; pca < nrOfPCsToCalculate; pca++) {
			double expVarPCA = PCA.getEigenValueVar(eigenValues, pca);
			double[] pca1ExpEigenVector = PCA.getEigenVector(eig, eigenValues, pca);
			for (int s = 0; s < dataset.colObjects.size(); s++) {
				datasetEV.getRawData()[s][pca] = pca1ExpEigenVector[s];
			}
			int pcaNr = pca + 1;
			cumExpVarPCA += expVarPCA;
			out.write(pcaNr + "\t" + eigenValues[eigenValues.length - 1 - pca] + "\t" + expVarPCA + "\t" + cumExpVarPCA + "\n");
			datasetEV.colObjects.set(pca, "Comp" + String.valueOf(pcaNr));
			System.out.println("PCA:\t" + pcaNr + "\t" + eigenValues[eigenValues.length - 1 - pca] + "\t" + expVarPCA + "\t" + cumExpVarPCA);
		}
		out.close();

		datasetEV.save(expressionFile + ".PCAOverSamplesEigenvectors.txt.gz");

		datasetEV.transposeDataset();

		datasetEV.save(expressionFile + ".PCAOverSamplesEigenvectorsTransposed.txt.gz");

		datasetEV.transposeDataset();
		System.out.println("Calculating PCs");
		System.out.println("Initializing PCA matrix");
		DoubleMatrixDataset<String, String> datasetPCAOverSamplesPCAs = new DoubleMatrixDataset<String, String>(dataset.rowObjects.size(), nrOfPCsToCalculate);
		datasetPCAOverSamplesPCAs.rowObjects = dataset.rowObjects;
		for (int s = 0; s < nrOfPCsToCalculate; s++) {
			datasetPCAOverSamplesPCAs.colObjects.set(s, "Comp" + String.valueOf(s + 1));
		}
		for (int p = 0; p < dataset.rowObjects.size(); p++) {
			for (int t = 0; t < nrOfPCsToCalculate; t++) {
				datasetPCAOverSamplesPCAs.getRawData()[p][t] = 0;
			}
		}


		ProgressBar pb = new ProgressBar(dataset.rowObjects.size(), "Calculating the PCA scores per probe: ");
		for (int probe = 0; probe < dataset.rowObjects.size(); probe++) {
			for (int sample1 = 0; sample1 < nrOfPCsToCalculate; sample1++) {
				for (int sample2 = 0; sample2 < dataset.colObjects.size(); sample2++) {
					double probeCoefficient = datasetEV.getRawData()[sample2][sample1];
					datasetPCAOverSamplesPCAs.getRawData()[probe][sample1] += probeCoefficient * dataset.getRawData()[probe][sample2];
				}
			}
			pb.iterate();
		}
		pb.close();

		String outfilename = expressionFile + ".PCAOverSamplesPrincipalComponents.txt.gz";
		System.out.println("Saving PCA scores: " + outfilename);
		datasetPCAOverSamplesPCAs.save(outfilename);

		return new Pair<DoubleMatrixDataset<String, String>, DoubleMatrixDataset<String, String>>(datasetPCAOverSamplesPCAs, datasetEV);
	}

	public void correctDataForPCs(DoubleMatrixDataset<String, String> dataset, String fileNamePrefix, int nrPCAsOverSamplesToRemove, int nrIntermediatePCAsOverSamplesToRemoveToOutput,
								  DoubleMatrixDataset<String, String> datasetPCAOverSamplesPCAs, DoubleMatrixDataset<String, String> datasetEV) throws IOException {
		String expressionFile = fileNamePrefix;
		System.out.println("\nInitializing residual gene expression matrix");

		if (dataset.colObjects.size() < nrPCAsOverSamplesToRemove) {
			int remainder = dataset.colObjects.size() % nrIntermediatePCAsOverSamplesToRemoveToOutput;
			nrPCAsOverSamplesToRemove = dataset.colObjects.size() - remainder;
		}

		for (int t = 0; t < nrPCAsOverSamplesToRemove; t++) {
			for (int p = 0; p < dataset.rowObjects.size(); p++) {
				for (int s = 0; s < dataset.colObjects.size(); s++) {
					dataset.getRawData()[p][s] -= datasetPCAOverSamplesPCAs.getRawData()[p][t] * datasetEV.getRawData()[s][t];
				}
			}
			int nrPCAs = t + 1;
			if (nrIntermediatePCAsOverSamplesToRemoveToOutput > 0 && nrPCAs % nrIntermediatePCAsOverSamplesToRemoveToOutput == 0) {
				dataset.save(expressionFile + "." + nrPCAs + "PCAsOverSamplesRemoved.txt.gz");
				System.out.println("Removed\t" + nrPCAs + "\tPCs. File:\t" + expressionFile + "." + nrPCAs + "PCAsOverSamplesRemoved.txt.gz");
			}

		}
		dataset.save(expressionFile + "." + nrPCAsOverSamplesToRemove + "PCAsOverSamplesRemoved.txt.gz");

	}

	public void repeatPCAOmitCertainPCAs(HashSet<Integer> pcasNotToRemove, String parentDir, String expressionFile, String nextInExp, String nextInExp2, int nrPCAsOverSamplesToRemove, int nrIntermediatePCAsOverSamplesToRemoveToOutput) throws IOException {
		System.out.println("Will write output to: " + parentDir);
		String[] files = Gpio.getListOfFiles(parentDir);
		String startExpressionFileName = expressionFile;
		File st = new File(startExpressionFileName);

		// strip the parent dir name
		parentDir += Gpio.getFileSeparator();
		String minimalFilename = st.getName();
		String[] expressionFileNameElems = minimalFilename.split("\\.");
		String eigenvectorFile = null;
		String principalComponentsFile = null;

		if (minimalFilename.contains("PCAsOverSamplesRemoved")) {
                    System.out.println("Warning, it seems like this data is already normalized for PCA's over sampels.");
		}

                eigenvectorFile = nextInExp;
                principalComponentsFile = nextInExp2;

		System.out.println("Detected core file name to be: " + minimalFilename);

		DoubleMatrixDataset<String, String> expressionDataset = new DoubleMatrixDataset<String, String>(parentDir + minimalFilename);
		DoubleMatrixDataset<String, String> datasetPCAOverSamplesPCAs = new DoubleMatrixDataset<String, String>(principalComponentsFile);
		DoubleMatrixDataset<String, String> datasetEV = new DoubleMatrixDataset<String, String>(eigenvectorFile);

		if (expressionDataset.colObjects.size() < nrPCAsOverSamplesToRemove) {
			int remainder = expressionDataset.colObjects.size() % nrIntermediatePCAsOverSamplesToRemoveToOutput;
			nrPCAsOverSamplesToRemove = expressionDataset.colObjects.size() - remainder;
		}

//        DoubleMatrixDataset<String, String> datasetResidualExpressionBasedOnPCAOverSamples = new DoubleMatrixDataset<String, String>(expressionDataset.rowObjects.size(), expressionDataset.colObjects.size());
//        datasetResidualExpressionBasedOnPCAOverSamples.rowObjects = expressionDataset.rowObjects;
//        datasetResidualExpressionBasedOnPCAOverSamples.colObjects = expressionDataset.colObjects;
//
//        for (int p = 0; p < expressionDataset.rowObjects.size(); p++) {
//            System.arraycopy(expressionDataset.getRawData()[p], 0, datasetResidualExpressionBasedOnPCAOverSamples.getRawData()[p], 0, expressionDataset.colObjects.size());
//        }

		if (minimalFilename.endsWith(".txt")) {
			minimalFilename = minimalFilename.substring(0, minimalFilename.length() - 4);
		} else if (minimalFilename.endsWith(".txt.gz")) {
			minimalFilename = minimalFilename.substring(0, minimalFilename.length() - 7);
		}

		for (int t = 0; t < nrPCAsOverSamplesToRemove; t++) {
			if (!pcasNotToRemove.contains(t + 1)) {

				for (int p = 0; p < expressionDataset.rowObjects.size(); p++) {
					for (int s = 0; s < expressionDataset.colObjects.size(); s++) {
						//datasetResidualExpressionBasedOnPCAOverSamples.rawData[p][s]-= datasetPCAOverSamplesPCAs.rawData[p][t] * datasetEV.rawData[s][t];
						expressionDataset.getRawData()[p][s] -= datasetPCAOverSamplesPCAs.getRawData()[p][t] * datasetEV.getRawData()[s][t];
					}
				}
			} else {
				System.out.println("Omitting PCA: " + (t + 1) + " since this component is under genetic control");
			}

			int nrPCAs = t + 1;

			if (nrIntermediatePCAsOverSamplesToRemoveToOutput > 0 && nrPCAs % nrIntermediatePCAsOverSamplesToRemoveToOutput == 0) {
				//datasetResidualExpressionBasedOnPCAOverSamples.save(expressionFile + "." + nrPCAs + "PCAsOverSamplesRemoved.txt");
				expressionDataset.save(parentDir + minimalFilename + "." + nrPCAs + "PCAsOverSamplesRemoved-GeneticVectorsNotRemoved.txt.gz");
				System.out.println("Removed\t" + nrPCAs + "\tPCs. File:\t" + minimalFilename + "." + nrPCAs + "PCAsOverSamplesRemoved-GeneticVectorsNotRemoved.txt.gz");
			}

		}
		//datasetResidualExpressionBasedOnPCAOverSamples.save(expressionFile + "." + nrPCAsOverSamplesToRemove + "PCAsOverSamplesRemoved.txt");
		expressionDataset.save(parentDir + minimalFilename + "." + nrPCAsOverSamplesToRemove + "PCAsOverSamplesRemoved-GeneticVectorsNotRemoved.txt.gz");

		System.out.println("Done\n");
	}

	private void correctForCovariate(double[][] rawdata, double[][] covariateValues, int covariateToCorrect) {
		for (int probe = 0; probe < rawdata.length; probe++) {
			double[] y = rawdata[probe];
			double meanY = JSci.maths.ArrayMath.mean(y);
			double varianceY = JSci.maths.ArrayMath.variance(y);
			double[] x = covariateValues[covariateToCorrect];


			double[] rc = Regression.getLinearRegressionCoefficients(x, y);
			double correlation = JSci.maths.ArrayMath.correlation(x, y);
			double propExplainedVarianceTrait = correlation * correlation - 1.0d / (double) y.length;

			if (propExplainedVarianceTrait < 0) {
				propExplainedVarianceTrait = 0;
			}

//            explainedVariancePerEQTLProbe[d][(int) Math.round(propExplainedVarianceTrait * 100d)]++;
			double[] rawDataUpdated = new double[x.length];
			for (int s = 0; s < x.length; s++) {
				double residual = y[s] - x[s] * rc[0];
				rawDataUpdated[s] = residual;
			}

			double meanUpdated = JSci.maths.ArrayMath.mean(rawDataUpdated);
			double stdDevRatio = JSci.maths.ArrayMath.standardDeviation(rawDataUpdated) / Math.sqrt(varianceY);
			for (int s = 0; s < x.length; s++) {
				rawDataUpdated[s] -= meanUpdated;
				rawDataUpdated[s] /= stdDevRatio;
				rawDataUpdated[s] += meanY;
			}
			System.arraycopy(rawDataUpdated, 0, rawdata[probe], 0, x.length);
		}
	}

	// NOTE: this new code switches around columns and rows for the covariate matrix
	private Pair<DoubleMatrixDataset<String, String>, DoubleMatrixDataset<String, String>> loadCovariateValues(String covariatesToRemove, DoubleMatrixDataset<String, String> dataset) throws IOException {
		System.out.println("- Removing covariates as defined in: " + covariatesToRemove);
		TextFile covariates = new TextFile(covariatesToRemove, TextFile.R);
		int numRows = covariates.countLines() - 1; // minus the header :)
		int numCols = covariates.countCols(TextFile.tab) - 1; // minus the header's row identifier (if any)

		if (numRows == 0 || numCols == 0) {
			System.err.println("Covariate file is empty, but no covariates found in file! Is your file format correct?");
			System.err.println("The program is expecting the following: tab separated, one covariate per row, one sample per column, with sample identifiers identical to your --in file.");
			System.exit(0);
		} else {
			System.out.println("Covariate file has " + numRows + " rows and " + numCols + " columns");
		}


		// first hash up which samples are in the dataset
		HashMap<String, Integer> samplesInDatasetIndex = new HashMap<String, Integer>();
		String[] allSamplesInDataset = dataset.colObjects.toArray(new String[0]);
		for (int i = 0; i < allSamplesInDataset.length; i++) {
			samplesInDatasetIndex.put(allSamplesInDataset[i], i);
		}

		// read the column names from the covariate file
		// expect the samples on the columns
		String[] elems = covariates.readLineElemsReturnReference(TextFile.tab); // header

		int ctr = 0;
		boolean[] sampleInDatasetIncludedInCovariates = new boolean[dataset.colObjects.size()];
		ArrayList<String> columnNames = new ArrayList<String>();
		for (int i = 1; i < elems.length; i++) {
			Integer index = samplesInDatasetIndex.get(elems[i]);
			columnNames.add(elems[i]);
			if (index != null) {
				sampleInDatasetIncludedInCovariates[index] = true;
				ctr++;
			}
		}

		// read the covariate names, expect them to be on the rows
		ArrayList<String> rowNames = new ArrayList<String>();
		elems = covariates.readLineElemsReturnReference(TextFile.tab); // first line
		while (elems != null) {
			rowNames.add(elems[0]);
			elems = covariates.readLineElemsReturnReference(TextFile.tab);
		}
		covariates.close();

		boolean isTransposed = false;
		if (ctr == 0) {
			System.err.println("No matching samples detected between covariate file and dataset. Maybe your covariate file needs to be transposed? Will test that for you now:");
			for (String rowName : rowNames) {
				Integer index = samplesInDatasetIndex.get(rowName);
				if (index != null) {
					sampleInDatasetIncludedInCovariates[index] = true;
					ctr++;
				}
			}

			if (ctr == 0) {
				System.err.println("Transposing the data does not seem to resolve the issue. Please check your sample identifiers.");
				System.exit(0);
			} else {
				System.out.println("Transposing the covariate file reveals: " + ctr + " samples present.");
				isTransposed = true;

			}


		}

//        if (dataset.colObjects.size() != numSamples) {
//            System.out.println("Covariates loaded from: " + covariatesToRemove + ", but the number of samples does not correspond! " + numSamples + " in covariates file, " + dataset.colObjects.size() + " in dataset...");
//            System.out.println("Please note that missing samples will be removed from your eventual corrected --in file.");
//        }
        if(!isTransposed && ctr<=numRows+2 || isTransposed && ctr<=numCols+2){
            System.err.println("Less samples present than minimaly required for the normalization, (minimaly covariats+3 samples needed).");
            System.exit(0);
        }
		if (ctr < dataset.colObjects.size()) {
			System.err.println("Covariates loaded from: " + covariatesToRemove + ", but not all samples present in covariates file! " + ctr + " present in covariates file, out of " + dataset.colObjects.size() + " in dataset...");
			System.out.println("Your dataset will be adjusted accordingly.");
		}
		int nrCovariates = numRows;
		if (isTransposed) {
			nrCovariates = numCols;
		}

		// make matrix with equal sample size
		double[][] covariateValues = new double[nrCovariates][dataset.colObjects.size()];
		for (int row = 0; row < covariateValues.length; row++) {
			for (int col = 0; col < covariateValues[row].length; col++) {
				covariateValues[row][col] = Double.NaN;
			}
		}

		int lineCtr = 0;
		covariates.open();
		String[] headerElems = covariates.readLineElemsReturnReference(TextFile.tab); // header
		elems = covariates.readLineElemsReturnReference(TextFile.tab);
		while (elems != null) {
			if (isTransposed) {
				String sampleName = elems[0];
				Integer sampleIdInDataset = samplesInDatasetIndex.get(sampleName);
				if (sampleIdInDataset != null) {
					for (int i = 1; i < elems.length; i++) {
						try {
							covariateValues[i - 1][sampleIdInDataset] = Double.parseDouble(elems[i]);
						} catch (NumberFormatException e) {
//                            System.out.println("WARNING: " + elems[i] + " is not a numeric value! in " + covariatesToRemove + " at line: " + (lineCtr + 1) + ".");
//                            covariateValues[i - 1][sampleIdInDataset] = Double.NaN;
//                            sampleInDatasetIncludedInCovariates[sampleIdInDataset] = false;
						}
					}
				}
			} else {
				for (int i = 1; i < elems.length; i++) {
					String sampleName = headerElems[i];
					Integer sampleIdInDataset = samplesInDatasetIndex.get(sampleName);
					if (sampleIdInDataset != null) {
						try {
							covariateValues[lineCtr][sampleIdInDataset] = Double.parseDouble(elems[i]);
						} catch (NumberFormatException e) {
//                            System.out.println("WARNING: " + elems[i] + " is not a numeric value at line: " + (lineCtr + 1) + "\tcolumn: " + i);
						}
					}
				}
			}
			elems = covariates.readLineElemsReturnReference(TextFile.tab);
			lineCtr++;
		}
		covariates.close();

		// investigate how many covariates there actually is data for.
		int covariateCtr = 0;
		boolean[] includeCovariate = new boolean[covariateValues.length];
		for (int row = 0; row < covariateValues.length; row++) {
			int nrColsFilled = 0;
			for (int col = 0; col < covariateValues[row].length; col++) {
				if (!Double.isNaN(covariateValues[row][col])) {
					nrColsFilled++;
				}
			}

			if (nrColsFilled == 0) {
				// there's no data for this covariate....
				includeCovariate[row] = false;
			} else {
				includeCovariate[row] = true;
				covariateCtr++;
			}
		}

		if (covariateCtr == 0) {
			System.err.println("ERROR: none of your covariates seem to have valid numerical values.. Please check your covariate file.");
			System.exit(0);
		} else {
			System.out.println("After removing covariates without data, your dataset will have " + covariateCtr + " covariates (out of: " + covariateValues.length + ") .");
		}

		ArrayList<String> covariateNames = null;
		if (isTransposed) {
			covariateNames = columnNames;
		} else {
			covariateNames = rowNames;
		}

		if (covariateCtr != covariateValues.length) {
			// remove covariates with missing values
			System.out.println("Removing covariates that have no data at all.");
			double[][] newCovariateData = new double[covariateCtr][dataset.colObjects.size()];
			ArrayList<String> newCovariateNames = new ArrayList<String>();
			int newCovariateCTR = 0;
			for (int row = 0; row < covariateValues.length; row++) {
				if (includeCovariate[row]) {
					newCovariateNames.add(covariateNames.get(row));

					for (int col = 0; col < covariateValues[row].length; col++) {
						newCovariateData[newCovariateCTR][col] = covariateValues[row][col];

						// check whether we should include all samples, but don't remove yet: sync this with the expression/whatever dastaset
						if (Double.isNaN(covariateValues[row][col])) {
							sampleInDatasetIncludedInCovariates[col] = false;
						}
					}
					newCovariateCTR++;
				} else {
					System.out.println(covariateNames.get(row) + " removed.");
				}
			}


			nrCovariates = newCovariateCTR;
			covariateValues = newCovariateData;
			covariateNames = newCovariateNames;
		}
		System.out.println("");
		System.out.println("Remaining covariates: ");
		for (String s : covariateNames) {
			System.out.println(s);
		}
		System.out.println("");
		// investigate how many samples there actually is data for.
		for (int row = 0; row < covariateValues.length; row++) {
			for (int col = 0; col < covariateValues[row].length; col++) {
				if (Double.isNaN(covariateValues[row][col])) {
					sampleInDatasetIncludedInCovariates[col] = false;
				}
			}
		}

		int sampleCtr = 0;
		for (int q = 0; q < sampleInDatasetIncludedInCovariates.length; q++) {
			if (sampleInDatasetIncludedInCovariates[q]) {
				sampleCtr++;
			}
		}

		// remove samples that have a missing value for at least one covariate
//        if (sampleCtr == sampleInDatasetIncludedInCovariates.length) {
//            System.out.println("There were no missing values or samples in your covariate file. Sample size will remain unchanged.");
//            DoubleMatrixDataset<String, String> covariateDataset = new DoubleMatrixDataset<String, String>(covariateValues, dataset.rowObjects, covariateNames);
//            return new Pair<DoubleMatrixDataset<String, String>, DoubleMatrixDataset<String, String>>(covariateDataset, dataset);
//        } else {
		System.out.println("Your covariate corrected dataset will have " + sampleCtr + " samples, after removing samples with missing covariate values.");
		double[][] rawData = dataset.getRawData();
		double[][] newRawData = new double[rawData.length][sampleCtr];
		double[][] finalCovariateData = new double[nrCovariates][sampleCtr];
		ArrayList<String> newColObjects = new ArrayList<String>();

		for (int col = 0; col < dataset.colObjects.size(); col++) {
			if (sampleInDatasetIncludedInCovariates[col]) {
				newColObjects.add(dataset.colObjects.get(col));
			}
		}

		for (int row = 0; row < rawData.length; row++) {
			int includedSampleCtr = 0;
			for (int col = 0; col < dataset.colObjects.size(); col++) {
				if (sampleInDatasetIncludedInCovariates[col]) {
					// include sample
					newRawData[row][includedSampleCtr] = rawData[row][col];
					includedSampleCtr++;
				}
			}
		}

		for (int row = 0; row < covariateValues.length; row++) {
			int includedCovariateSampleCtr = 0;
			for (int col = 0; col < dataset.colObjects.size(); col++) {
				// replace covariate data...
				if (sampleInDatasetIncludedInCovariates[col]) {
					finalCovariateData[row][includedCovariateSampleCtr] = covariateValues[row][col];
					includedCovariateSampleCtr++;
				}
			}
		}

		DoubleMatrixDataset<String, String> covariateDataset = new DoubleMatrixDataset<String, String>(finalCovariateData, covariateNames, newColObjects);
		covariateDataset.save(covariatesToRemove + "-asLoadedByNormalizer.txt");
		DoubleMatrixDataset<String, String> newDataset = new DoubleMatrixDataset<String, String>(newRawData, dataset.rowObjects, newColObjects);
		newDataset.save(dataset.fileName + "-SampleSizeCorrectedForCovariates.txt");
		return new Pair<DoubleMatrixDataset<String, String>, DoubleMatrixDataset<String, String>>(covariateDataset, newDataset);
//        }
	}

	private String removeProbesWithZeroVariance(DoubleMatrixDataset<String, String> dataset, String outputFileNamePrefix) throws IOException {
		boolean[] dataHasZeroVariance = new boolean[dataset.nrRows];
		int nrRowsWithZeroVariance = 0;
		for (int row = 0; row < dataset.nrRows; row++) {
			double[] data = dataset.rawData[row];
			double var = JSci.maths.ArrayMath.variance(data);
			if (var == 0d) {
				System.out.println("Removing probe with zero variance: " + dataset.rowObjects.get(row) + " on line " + (row + 1));
				nrRowsWithZeroVariance++;
				dataHasZeroVariance[row] = true;
			}
		}

		if (nrRowsWithZeroVariance > 0) {
			int newNrRows = dataset.nrRows - nrRowsWithZeroVariance;
			if (newNrRows == 0) {
				System.err.println("ERROR: all probes have zero variance!");
				System.exit(-1);
			}


			double[][] newData = new double[newNrRows][dataset.nrCols];
			int ctr = 0;
			ArrayList<String> newRowHeader = new ArrayList<String>();
			for (int row = 0; row < dataset.nrRows; row++) {
				if (!dataHasZeroVariance[row]) {
					newData[ctr] = dataset.rawData[row];
					newRowHeader.add(dataset.rowObjects.get(row));
					ctr++;
				}
			}

			dataset.rawData = newData;
			dataset.rowObjects = newRowHeader;
			dataset.recalculateHashMaps();
			String outputFileName = outputFileNamePrefix + ".ProbesWithZeroVarianceRemoved";
			dataset.save(outputFileName + ".txt.gz");
			return outputFileName;
		}

		return outputFileNamePrefix;
	}
}
