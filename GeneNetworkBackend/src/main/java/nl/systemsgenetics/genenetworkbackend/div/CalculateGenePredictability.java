/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.genenetworkbackend.div;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import com.opencsv.CSVWriter;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashSet;
import org.apache.commons.math3.stat.descriptive.moment.Kurtosis;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.Skewness;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.stats.MannWhitneyUTest2;

/**
 *
 * @author patri
 */
public class CalculateGenePredictability {

	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) throws IOException {

		File predictionMatrixFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Data31995Genes05-12-2017\\PCA_01_02_2018\\predictions\\reactome_predictions.txt.gz");
		File annotationMatrixFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Data31995Genes05-12-2017\\PCA_01_02_2018\\PathwayMatrix\\Ensembl2Reactome_All_Levels.txt_matrix.txt.gz");
		File significantTermsFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Data31995Genes05-12-2017\\PCA_01_02_2018\\predictions\\reactome_predictions_bonSigTerms_alsoInGoP.txt");
		File outputFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Data31995Genes05-12-2017\\PCA_01_02_2018\\predictions\\reactome_predictions_genePredictability_alsoInGoP.txt");

//		File predictionMatrixFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Data31995Genes05-12-2017\\PCA_01_02_2018\\predictions\\hpo_predictions.txt.gz");
//		File annotationMatrixFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Data31995Genes05-12-2017\\PCA_01_02_2018\\PathwayMatrix\\ALL_SOURCES_ALL_FREQUENCIES_phenotype_to_genes.txt_matrix.txt.gz");
//		File significantTermsFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Data31995Genes05-12-2017\\PCA_01_02_2018\\predictions\\hpo_predictions_bonSigTerms.txt");
//		File outputFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Data31995Genes05-12-2017\\PCA_01_02_2018\\predictions\\hpo_predictions_genePredictability.txt");
//		
//		File predictionMatrixFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Data31995Genes05-12-2017\\PCA_01_02_2018\\predictions\\go_F_predictions.txt.gz");
//		File annotationMatrixFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Data31995Genes05-12-2017\\PCA_01_02_2018\\PathwayMatrix\\goa_human.gaf_F_matrix.txt.gz");
//		File significantTermsFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Data31995Genes05-12-2017\\PCA_01_02_2018\\predictions\\go_F_predictions_bonSigTerms.txt");
//		File outputFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Data31995Genes05-12-2017\\PCA_01_02_2018\\predictions\\go_F_predictions_genePredictability.txt");
//		
//		File predictionMatrixFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Data31995Genes05-12-2017\\PCA_01_02_2018\\predictions\\kegg_predictions.txt.gz");
//		File annotationMatrixFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Data31995Genes05-12-2017\\PCA_01_02_2018\\PathwayMatrix\\c2.cp.kegg.v6.1.entrez.gmt_matrix.txt.gz");
//		File significantTermsFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Data31995Genes05-12-2017\\PCA_01_02_2018\\predictions\\kegg_predictions_bonSigTerms.txt");
//		File outputFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Data31995Genes05-12-2017\\PCA_01_02_2018\\predictions\\kegg_predictions_genePredictability.txt");
//		
//		File predictionMatrixFile = new File("/groups/umcg-wijmenga/tmp04/umcg-svandam/GeneNetwork/Data31995Genes05-12-2017/GeneNetwork_V2_01-02-2018/Covariates/PCA/predictions/go_P_predictions.txt.gz");
//		File annotationMatrixFile = new File("/groups/umcg-wijmenga/tmp04/umcg-svandam/GeneNetwork/Data31995Genes05-12-2017/GeneNetwork_V2_01-02-2018/Covariates/PCA/PathwayMatrix/goa_human.gaf_P_matrix.txt.gz");
//		File significantTermsFile = new File("/groups/umcg-wijmenga/tmp04/umcg-svandam/GeneNetwork/Data31995Genes05-12-2017/GeneNetwork_V2_01-02-2018/Covariates/PCA/predictions/go_P_predictions_bonSigTerms_alsoInReactome.txt");
//		File outputFile = new File("/groups/umcg-wijmenga/tmp04/umcg-svandam/GeneNetwork/Data31995Genes05-12-2017/GeneNetwork_V2_01-02-2018/Covariates/PCA/predictions/go_P_predictions_genePredictability_alsoInReactome.txt");
		
		
		LinkedHashSet<String> significantTerms = loadSignificantTerms(significantTermsFile);

		DoubleMatrixDataset<String, String> predictionMatrix = DoubleMatrixDataset.loadDoubleData(predictionMatrixFile.getAbsolutePath());
		DoubleMatrixDataset<String, String> annotationMatrix = DoubleMatrixDataset.loadDoubleData(annotationMatrixFile.getAbsolutePath());

		DoubleMatrixDataset<String, String> predictionMatrixSignificant = predictionMatrix.viewColSelection(significantTerms);
		DoubleMatrixDataset<String, String> annotationMatrixSignificant = annotationMatrix.viewColSelection(significantTerms);

		if (!predictionMatrixSignificant.getColObjects().equals(annotationMatrixSignificant.getColObjects())) {
			System.err.println("Differnce in terms");
			return;
		}

		if (!predictionMatrixSignificant.getRowObjects().equals(annotationMatrixSignificant.getRowObjects())) {
			System.err.println("Differnce in genes");
			return;
		}

		MannWhitneyUTest2 uTest = new MannWhitneyUTest2();
		Kurtosis kurtosisCalculator = new Kurtosis();
		Skewness skewnessCalculator = new Skewness();
		Mean annotatedMeanCalculator = new Mean();
		Mean notAnnotatedMeanCalculator = new Mean();

		double[] genePredictabilityZscores = new double[predictionMatrixSignificant.rows()];
		int[] pathwayCount = new int[predictionMatrixSignificant.rows()];
		double[] geneKurtosis = new double[predictionMatrixSignificant.rows()];
		double[] geneSkewness = new double[predictionMatrixSignificant.rows()];
		double[] geneAnnotatedMean = new double[predictionMatrixSignificant.rows()];
		double[] geneNotAnnotatedMean = new double[predictionMatrixSignificant.rows()];

		for (int g = 0; g < predictionMatrixSignificant.rows(); g++) {

			kurtosisCalculator.clear();
			skewnessCalculator.clear();
			annotatedMeanCalculator.clear();
			notAnnotatedMeanCalculator.clear();

			DoubleMatrix1D geneAnnotations = annotationMatrixSignificant.getRow(g);

			int geneAnnotationCount = geneAnnotations.cardinality();

			pathwayCount[g] = geneAnnotationCount;

			double[] zScoresAnnotatedPathways = new double[geneAnnotationCount];
			double[] zScoresOtherPathways = new double[annotationMatrixSignificant.columns() - geneAnnotationCount];

			int x = 0;
			int y = 0;

			for (int p = 0; p < geneAnnotations.size(); p++) {
				
				double z = predictionMatrixSignificant.getElementQuick(g, p);
				
				if (geneAnnotations.getQuick(p) != 0) {
					annotatedMeanCalculator.increment(z);
					zScoresAnnotatedPathways[x++] = z;
				} else {
					notAnnotatedMeanCalculator.increment(z);
					zScoresOtherPathways[y++] = z;
				}
				kurtosisCalculator.increment(z);
				skewnessCalculator.increment(z);
			}

			if (geneAnnotationCount >= 10) {
				uTest.setData(zScoresOtherPathways, zScoresAnnotatedPathways);
				genePredictabilityZscores[g] = uTest.getZ();
			} else {
				genePredictabilityZscores[g] = Double.NaN;
			}

			geneKurtosis[g] = kurtosisCalculator.getResult();
			geneSkewness[g] = skewnessCalculator.getResult();
			geneAnnotatedMean[g] = annotatedMeanCalculator.getResult();
			geneNotAnnotatedMean[g] = notAnnotatedMeanCalculator.getResult();

		}

		CSVWriter writer = new CSVWriter(new FileWriter(outputFile), '\t', '\0', '\0', "\n");

		String[] outputLine = new String[7];
		int c = 0;
		outputLine[c++] = "Gene";
		outputLine[c++] = "Z-score";
		outputLine[c++] = "Skewness";
		outputLine[c++] = "Kurtosis";
		outputLine[c++] = "MeanNotAnnotated";
		outputLine[c++] = "MeanAnnotated";
		outputLine[c++] = "Annoted_pathways";
		writer.writeNext(outputLine);

		ArrayList<String> geneNames = predictionMatrixSignificant.getRowObjects();
		for (int g = 0; g < predictionMatrixSignificant.rows(); g++) {
			c = 0;
			outputLine[c++] = geneNames.get(g);
			outputLine[c++] = String.valueOf(genePredictabilityZscores[g]);
			outputLine[c++] = String.valueOf(geneSkewness[g]);
			outputLine[c++] = String.valueOf(geneKurtosis[g]);
			outputLine[c++] = String.valueOf(geneNotAnnotatedMean[g]);
			outputLine[c++] = String.valueOf(geneAnnotatedMean[g]);
			outputLine[c++] = String.valueOf(pathwayCount[g]);
			writer.writeNext(outputLine);
		}

		writer.close();

	}

	public static LinkedHashSet<String> loadSignificantTerms(File significantTermsFile) throws IOException {

		LinkedHashSet<String> significantTerms = new LinkedHashSet<>();

		BufferedReader reader = new BufferedReader(new FileReader(significantTermsFile));

		String line;
		while ((line = reader.readLine()) != null) {
			significantTerms.add(line);
		}

		return significantTerms;

	}

}
