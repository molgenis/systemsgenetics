/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.genenetworkbackend.hpo;

import umcg.genetica.io.hpo.HpoOntology;
import java.io.File;
import org.biojava.nbio.ontology.Term;

/**
 *
 * @author patri
 */
public class Test {

	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) throws Exception {
//		final File updatedPredictionMatrixFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Data31995Genes05-12-2017\\PCA_01_02_2018\\predictions\\hpo_predictions_improved.txt.gz");
//		
//		DoubleMatrixDataset<String, String> matrix = DoubleMatrixDataset.loadDoubleData(updatedPredictionMatrixFile.getAbsolutePath());
//		
//		System.out.println(matrix.getElement("ENSG00000165917", "HP:0001324"));
//		

		final File hpoOboFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\HPO\\135\\hp.obo");

		HpoOntology hpoOntology = new HpoOntology(hpoOboFile);

		String a = "HP:0001639";
		String b = "HP:0011663";
		
		Term termA = hpoOntology.nameToTerm(a);
		Term termB = hpoOntology.nameToTerm(b);
		
		System.out.println(hpoOntology.isTermABelowTermB(termA, termB));
		System.out.println(hpoOntology.isTermABelowTermB(termB, termA));

	}

}
