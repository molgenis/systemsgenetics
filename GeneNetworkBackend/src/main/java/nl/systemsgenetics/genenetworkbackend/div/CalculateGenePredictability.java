/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.genenetworkbackend.div;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.LinkedHashSet;
import org.apache.commons.math3.stat.inference.MannWhitneyUTest;
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
		File significantTermsFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Data31995Genes05-12-2017\\PCA_01_02_2018\\predictions\\reactome_predictions_bonSigTerms.txt");


		LinkedHashSet<String> significantTerms = loadSignificantTerms(significantTermsFile);		
		
		DoubleMatrixDataset predictionMatrix = DoubleMatrixDataset.loadDoubleData(predictionMatrixFile.getAbsolutePath());
		DoubleMatrixDataset annotationMatrix = DoubleMatrixDataset.loadDoubleData(annotationMatrixFile.getAbsolutePath());
		
		DoubleMatrixDataset predictionMatrixSignificant = predictionMatrix.viewColSelection(significantTerms);
		DoubleMatrixDataset annotationMatrixSignificant = annotationMatrix.viewColSelection(significantTerms);
		
	
		if(!predictionMatrixSignificant.getColObjects().equals(annotationMatrixSignificant.getColObjects())){
			System.err.println("Differnce in terms");
			return;
		}
		
		if(!predictionMatrixSignificant.getRowObjects().equals(annotationMatrixSignificant.getRowObjects())){
			System.err.println("Differnce in genes");
			return;
		}
		
		for(int g = 0; g < predictionMatrixSignificant.rows() ; g++){
			
			DoubleMatrix1D geneAnnotations = annotationMatrixSignificant.getRow(g);
			
			int geneAnnotationCount = geneAnnotations.cardinality();
			
			if(geneAnnotationCount >= 10){
				//MannWhitneyUTest
				//WilcoxonMannWhitney wmw = new WilcoxonMannWhitney();
				//double pValue = wmw.returnWilcoxonMannWhitneyPValue(vals1, vals2);
				
			}
			
		}
		
		

	}

	private static LinkedHashSet<String> loadSignificantTerms(File significantTermsFile) throws IOException {

		LinkedHashSet<String> significantTerms = new LinkedHashSet<>();

		BufferedReader reader = new BufferedReader(new FileReader(significantTermsFile));

		String line;
		while ((line = reader.readLine()) != null) {
			significantTerms.add(line);
		}

		return significantTerms;

	}

}
