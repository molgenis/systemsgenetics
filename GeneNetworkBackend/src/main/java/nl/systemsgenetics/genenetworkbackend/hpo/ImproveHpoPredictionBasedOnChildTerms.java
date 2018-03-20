/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.genenetworkbackend.hpo;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import com.opencsv.CSVParser;
import com.opencsv.CSVParserBuilder;
import com.opencsv.CSVReader;
import com.opencsv.CSVReaderBuilder;
import com.opencsv.CSVWriter;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import org.apache.commons.math3.util.FastMath;
import org.biojava.nbio.ontology.Ontology;
import org.biojava.nbio.ontology.Term;
import org.biojava.nbio.ontology.Triple;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.stats.MannWhitneyUTest2;

/**
 *
 * @author patri
 */
public class ImproveHpoPredictionBasedOnChildTerms {

	final HashMap<String, UpdatedPredictionInfo> checkedHpo = new HashMap<>();
	final DoubleMatrixDataset<String, String> predictionMatrix;
	final DoubleMatrixDataset<String, String> annotationMatrix;
	final Ontology hpoOntology;
	final Term is_a;
	final MannWhitneyUTest2 uTest = new MannWhitneyUTest2();
	final double sigCutoff;

	/**
	 * @param args the command line arguments
	 * @throws java.lang.Exception
	 */
	public static void main(String[] args) throws Exception {

		final File predictionMatrixFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Data31995Genes05-12-2017\\PCA_01_02_2018\\predictions\\hpo_predictions.txt.gz");
		final File annotationMatrixFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Data31995Genes05-12-2017\\PCA_01_02_2018\\PathwayMatrix\\ALL_SOURCES_ALL_FREQUENCIES_phenotype_to_genes.txt_matrix.txt.gz");
		final File predictedHpoTermFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Data31995Genes05-12-2017\\PCA_01_02_2018\\predictions\\hpo_predictions_auc_bonferroni.txt");
//		final File predictionMatrixFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Data31995Genes05-12-2017\\PCA_01_02_2018\\predictions\\hpo_predictions_testSet.txt");
//		final File annotationMatrixFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Data31995Genes05-12-2017\\PCA_01_02_2018\\PathwayMatrix\\hpo_annotation_testSet.txt");
//		final File predictedHpoTermFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Data31995Genes05-12-2017\\PCA_01_02_2018\\predictions\\hpo_predictions_testSet_auc_bonferroni.txt");
		final File hpoOboFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\HPO\\135\\hp.obo");
		final File ouputLogFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Data31995Genes05-12-2017\\PCA_01_02_2018\\predictions\\hpo_predictions_improved.log");
		final File updatedPredictionMatrixFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Data31995Genes05-12-2017\\PCA_01_02_2018\\predictions\\hpo_predictions_improved.txt.gz");

		LinkedHashSet<String> predictedHpoTerms = readPredictedHpoTermFile(predictedHpoTermFile);

		DoubleMatrixDataset<String, String> predictionMatrixFull = DoubleMatrixDataset.loadDoubleData(predictionMatrixFile.getAbsolutePath());
		DoubleMatrixDataset<String, String> annotationMatrixFull = DoubleMatrixDataset.loadDoubleData(annotationMatrixFile.getAbsolutePath());

		DoubleMatrixDataset<String, String> predictionMatrixPredicted = predictionMatrixFull.viewColSelection(predictedHpoTerms);
		DoubleMatrixDataset<String, String> annotationMatrixPredicted = annotationMatrixFull.viewColSelection(predictedHpoTerms);

		Ontology hpoOntology = HpoFinder.loadHpoOntology(hpoOboFile);

		ImproveHpoPredictionBasedOnChildTerms improver = new ImproveHpoPredictionBasedOnChildTerms(predictionMatrixPredicted, annotationMatrixPredicted, hpoOntology);

		HashMap<String, UpdatedPredictionInfo> checkedHpoInfo = improver.run();

		System.out.println("Done with improving");

		CSVWriter writer = new CSVWriter(new FileWriter(ouputLogFile), '\t', '\0', '\0', "\n");

		String[] outputLine = new String[11];
		int c = 0;
		outputLine[c++] = "HPO";
		outputLine[c++] = "Gene_count";
		outputLine[c++] = "Origanl_AUC";
		outputLine[c++] = "Orignal_Pvalue";
		outputLine[c++] = "Updated_AUC";
		outputLine[c++] = "Updated_Pvalue";
		outputLine[c++] = "Is_significant";
		outputLine[c++] = "Distance_to_top";
		outputLine[c++] = "Number_of_child_terms";
		outputLine[c++] = "Number_of_child_terms_used";
		outputLine[c++] = "Child_terms_used";
		writer.writeNext(outputLine);

		for (UpdatedPredictionInfo pi : checkedHpoInfo.values()) {
			c = 0;
			outputLine[c++] = pi.getHpo();
			outputLine[c++] = String.valueOf(pi.getGeneCount());
			outputLine[c++] = String.valueOf(pi.getOriginalAuc());
			outputLine[c++] = String.valueOf(pi.getOriginalPvalue());
			outputLine[c++] = String.valueOf(pi.getUpdatedAuc());
			outputLine[c++] = String.valueOf(pi.getUpdatedPvalue());
			outputLine[c++] = String.valueOf(pi.isIsSignificant());
			outputLine[c++] = "-";
			outputLine[c++] = String.valueOf(pi.getChildTermCount());
			outputLine[c++] = String.valueOf(pi.getUsedChildTerms().size());
			outputLine[c++] = String.join(";", pi.getUsedChildTerms());
			writer.writeNext(outputLine);
		}
		writer.close();

		improver.writeUpdatedMatrix(updatedPredictionMatrixFile);
	}

	public ImproveHpoPredictionBasedOnChildTerms(DoubleMatrixDataset<String, String> predictionMatrix, DoubleMatrixDataset<String, String> annotationMatrix, Ontology hpoOntology) {
		this.predictionMatrix = predictionMatrix;
		this.annotationMatrix = annotationMatrix;
		this.hpoOntology = hpoOntology;
		is_a = hpoOntology.getTerm("is_a");
		sigCutoff = 0.05 / (annotationMatrix.columns() * 2);
	}

	private HashMap<String, UpdatedPredictionInfo> run() {

		int nextReport = 100;

		for (String hpo : annotationMatrix.getColObjects()) {
			predict(hpo);
			if (checkedHpo.size() >= nextReport) {
				System.out.println("Processed " + checkedHpo.size() + " HPO terms");
				nextReport += 100;
			}
		}

		return this.checkedHpo;
	}

	private void predict(String hpo) {

		if (checkedHpo.containsKey(hpo)) {
			return;
		}

		final Term hpoTerm = hpoOntology.getTerm(hpo);
		final HashSet<String> hpoChildern = new HashSet<>();

		for (Triple t : hpoOntology.getTriples(null, hpoTerm, is_a)) {
			hpoChildern.add(t.getSubject().getName());
		}

		DoubleMatrix1D hpoAnnotations = annotationMatrix.getCol(hpo);
		DoubleMatrix1D hpoPredictions = predictionMatrix.getCol(hpo);

		int hpoAnnotationCount = hpoAnnotations.cardinality();

		double[] zScoresAnnotatedGenes = new double[hpoAnnotationCount];
		double[] zScoresOtherGenes = new double[annotationMatrix.rows() - hpoAnnotationCount];

		int x = 0;
		int y = 0;

		for (int g = 0; g < hpoAnnotations.size(); g++) {
			if (hpoAnnotations.getQuick(g) != 0) {
				zScoresAnnotatedGenes[x++] = hpoPredictions.getQuick(g);
			} else {
				zScoresOtherGenes[y++] = hpoPredictions.getQuick(g);
			}
		}

		uTest.setData(zScoresAnnotatedGenes, zScoresOtherGenes);

		double originalAuc = uTest.getAuc();
		double originalP = uTest.getP();

		ArrayList<String> usableChildern = new ArrayList<>();
		for (String hpoChild : hpoChildern) {

			if (!predictionMatrix.containsCol(hpoChild)) {
				continue;
			}

			if (!checkedHpo.containsKey(hpoChild)) {
				predict(hpoChild);
			}

			//if (checkedHpo.get(hpoChild).isIsSignificant()) {
			usableChildern.add(hpoChild);
			//}

		}

		double updatedAuc;
		double updatedP;

		final double[] alternativeGeneZscores = new double[annotationMatrix.rows()];
		if (usableChildern.size() >= 1) {

			int sumSquareOfWeigths = 0;

			for (String hpoChild : usableChildern) {

				DoubleMatrix1D childGeneZ = predictionMatrix.getCol(hpoChild);

				int hpoChildGeneCount = checkedHpo.get(hpoChild).geneCount;
				//sumSquareOfWeigths += hpoChildGeneCount;
				//double weight = FastMath.sqrt(hpoChildGeneCount);
				sumSquareOfWeigths += 1;
				double weight = 1;

				for (int i = 0; i < alternativeGeneZscores.length; i++) {
					alternativeGeneZscores[i] += (weight * childGeneZ.getQuick(i));
				}

			}

//			double weight = FastMath.sqrt(hpoAnnotationCount);
//			sumSquareOfWeigths += hpoAnnotationCount;
			double weight = 1;
			sumSquareOfWeigths += 1;

			for (int i = 0; i < alternativeGeneZscores.length; i++) {
				alternativeGeneZscores[i] += (weight * hpoPredictions.getQuick(i));
			}

			double denominator = Math.sqrt(sumSquareOfWeigths);

			for (int i = 0; i < alternativeGeneZscores.length; i++) {
				alternativeGeneZscores[i] = (alternativeGeneZscores[i] / denominator);
				if (Double.isNaN(alternativeGeneZscores[i])) {
					alternativeGeneZscores[i] = 0;
				}
			}

			if (hpo.equals("HP:0001324")) {
				System.out.println("HP:0001324 - RAPSN: " + alternativeGeneZscores[predictionMatrix.getHashRows().get("ENSG00000165917")]);
			}

			zScoresAnnotatedGenes = new double[hpoAnnotationCount];
			zScoresOtherGenes = new double[annotationMatrix.rows() - hpoAnnotationCount];

			x = 0;
			y = 0;

			for (int g = 0; g < hpoAnnotations.size(); g++) {
				if (hpoAnnotations.getQuick(g) != 0) {
					zScoresAnnotatedGenes[x++] = alternativeGeneZscores[g];
				} else {
					zScoresOtherGenes[y++] = alternativeGeneZscores[g];
				}
			}

			uTest.setData(zScoresAnnotatedGenes, zScoresOtherGenes);

			updatedAuc = uTest.getAuc();
			updatedP = uTest.getP();

		} else {
			updatedAuc = Double.NaN;
			updatedP = Double.NaN;
		}

		final boolean isSignificant;
		if (!Double.isNaN(updatedP) && updatedP < originalP) {
			isSignificant = updatedP < sigCutoff;

			for (int g = 0; g < hpoAnnotations.size(); g++) {
				hpoPredictions.setQuick(g, alternativeGeneZscores[g]);
			}

		} else {
			isSignificant = originalP < sigCutoff;
		}

//		System.out.println("");
//		System.out.println("-------------");
//		System.out.println("");
//		System.out.println("Term: " + hpo);
//		System.out.println("U1 " + uTest.getU1());
//		System.out.println("U2 " + uTest.getU2());
//		System.out.println("N1 " + uTest.getN1());
//		System.out.println("N2 " + uTest.getN2());
//		System.out.println("Z " + uTest.getZ());
//		System.out.println("P " + uTest.getP());
//		System.out.println("AUC " + uTest.getAuc());
//
		UpdatedPredictionInfo predictionInfo = new UpdatedPredictionInfo(hpo, hpoAnnotationCount, originalAuc, originalP, updatedAuc, updatedP, usableChildern, usableChildern.size(), isSignificant);

		checkedHpo.put(hpo, predictionInfo);

	}

	private static LinkedHashSet<String> readPredictedHpoTermFile(File predictedHpoTermFile) throws FileNotFoundException, IOException {

		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		final CSVReader reader = new CSVReaderBuilder(new BufferedReader(new FileReader(predictedHpoTermFile))).withSkipLines(1).withCSVParser(parser).build();

		LinkedHashSet<String> hpos = new LinkedHashSet<>();

		String[] nextLine;
		while ((nextLine = reader.readNext()) != null) {

			hpos.add(nextLine[0]);

		}

		reader.close();

		return hpos;

	}

	private void writeUpdatedMatrix(File updatedPredictionMatrixFile) throws IOException {
		predictionMatrix.save(updatedPredictionMatrixFile);
	}

	private class UpdatedPredictionInfo {

		final String hpo;
		final int geneCount;
		final double originalAuc;
		final double originalPvalue;
		final double updatedAuc;
		final double updatedPvalue;
		final ArrayList<String> usedChildTerms;
		final int childTermCount;
		final boolean isSignificant;

		public UpdatedPredictionInfo(String hpo, int geneCount, double originalAuc, double originalPvalue, double updatedAuc, double updatedPvalue, ArrayList<String> usedChildTerms, int childTermCount, boolean isSignificant) {
			this.hpo = hpo;
			this.geneCount = geneCount;
			this.originalAuc = originalAuc;
			this.originalPvalue = originalPvalue;
			this.updatedAuc = updatedAuc;
			this.updatedPvalue = updatedPvalue;
			this.usedChildTerms = usedChildTerms;
			this.childTermCount = childTermCount;
			this.isSignificant = isSignificant;
		}

		public String getHpo() {
			return hpo;
		}

		public double getOriginalAuc() {
			return originalAuc;
		}

		public double getOriginalPvalue() {
			return originalPvalue;
		}

		public double getUpdatedAuc() {
			return updatedAuc;
		}

		public double getUpdatedPvalue() {
			return updatedPvalue;
		}

		public ArrayList<String> getUsedChildTerms() {
			return usedChildTerms;
		}

		public int getChildTermCount() {
			return childTermCount;
		}

		public boolean isIsSignificant() {
			return isSignificant;
		}

		public int getGeneCount() {
			return geneCount;
		}

	}

}
