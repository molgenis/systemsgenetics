/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.functionenrichmentoftransqtls;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.jet.random.tdouble.StudentT;
import cern.jet.random.tdouble.engine.DRand;
import cern.jet.random.tdouble.engine.DoubleRandomEngine;
import cern.jet.stat.tdouble.Probability;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.LinkedHashSet;
import org.apache.commons.math3.stat.regression.SimpleRegression;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

/**
 *
 * @author patri
 */
public class CorrelateSumChi2ToPathways {

	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) throws IOException, Exception {

		final File pathwayMatrixFile = new File(args[0]);
		final File significantTermsFile = new File(args[1]);
		final File sumChi2MatrixFile = new File(args[2]);
		final File transQtlEnrichmentsMatrixFile = new File(args[3]);
		
		System.out.println("Pathway file: " + pathwayMatrixFile.getPath());
		System.out.println("Pathway significant terms file: " + significantTermsFile.getPath());
		System.out.println("SumChi2 file: " + sumChi2MatrixFile.getPath());
		System.out.println("Output file: " + transQtlEnrichmentsMatrixFile.getPath());
		
		LinkedHashSet<String> significantTerms = loadSignificantTerms(significantTermsFile);
		
		DoubleMatrixDataset<String, String> pathwayMatrix = DoubleMatrixDataset.loadDoubleData(pathwayMatrixFile.getPath());
		DoubleMatrixDataset<String, String> sumChi2Matrix = DoubleMatrixDataset.loadDoubleData(sumChi2MatrixFile.getPath());

		LinkedHashSet<String> genesInBoth = new LinkedHashSet<String>();

		for (String gene : pathwayMatrix.getHashRows().keySet()) {
			if (sumChi2Matrix.containsRow(gene)) {
				genesInBoth.add(gene);
			}
		}

		pathwayMatrix = pathwayMatrix.viewColSelection(significantTerms);
		
		pathwayMatrix = pathwayMatrix.viewRowSelection(genesInBoth);
		DoubleMatrixDataset<String, String> transQtlEnrichmentsMatrix = new DoubleMatrixDataset<String, String>(pathwayMatrix.getHashCols(), sumChi2Matrix.getHashCols());
		sumChi2Matrix = sumChi2Matrix.viewRowSelection(genesInBoth);

		System.out.println("Genes in both datasets: " + genesInBoth.size());

		System.out.println("Pathways to test: " + pathwayMatrix.columns());
		
		final SimpleRegression regression = new SimpleRegression();
		final DoubleRandomEngine randomEngine = new DRand();
		StudentT tDistColt = new StudentT(sumChi2Matrix.rows() / 2 - 2, randomEngine);
		
		for (String trait : sumChi2Matrix.getColObjects()) {
			
			System.out.println("Trait: " + trait);

			DoubleMatrix1D traitSumChi2 = sumChi2Matrix.getCol(trait);

			for (String pathway : pathwayMatrix.getColObjects()) {

				DoubleMatrix1D pathwayScores = pathwayMatrix.getCol(pathway);

				regression.clear();

				for (int i = 0; i < traitSumChi2.size(); ++i) {

					//System.out.println(traitSumChi2.get(i) + " & " + pathwayScores.get(i));
					
					regression.addData(traitSumChi2.get(i), pathwayScores.get(i));

				}

				double r = regression.getR();

				//System.out.println(trait + " " + pathway + " " + r);
				
				double t = r / (Math.sqrt((1 - r * r) / (double) (traitSumChi2.size() / 2 - 2)));
				double pValue;
				double zScore;
				if (t < 0) {
					pValue = tDistColt.cdf(t);
					if (pValue < 2.0E-323) {
						pValue = 2.0E-323;
					}
					zScore = Probability.normalInverse(pValue);
				} else {
					pValue = tDistColt.cdf(-t);
					if (pValue < 2.0E-323) {
						pValue = 2.0E-323;
					}
					zScore = -Probability.normalInverse(pValue);
				}
				pValue *= 2;

				transQtlEnrichmentsMatrix.setElement(pathway, trait, zScore);

			}

		}
		
		transQtlEnrichmentsMatrix.save(transQtlEnrichmentsMatrixFile);

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
