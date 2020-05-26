/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.genenetworkbackend.div;

import java.io.File;

/**
 *
 * @author patri
 */
public class CompareGenePredictability {

	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) {

		File goPPredictionMatrixFile = new File("/groups/umcg-wijmenga/tmp04/umcg-svandam/GeneNetwork/Data31995Genes05-12-2017/GeneNetwork_V2_01-02-2018/Covariates/PCA/predictions/go_P_predictions.txt.gz");
		File goPannotationMatrixFile = new File("/groups/umcg-wijmenga/tmp04/umcg-svandam/GeneNetwork/Data31995Genes05-12-2017/GeneNetwork_V2_01-02-2018/Covariates/PCA/PathwayMatrix/goa_human.gaf_P_matrix.txt.gz");
		File goPsignificantTermsFile = new File("/groups/umcg-wijmenga/tmp04/umcg-svandam/GeneNetwork/Data31995Genes05-12-2017/GeneNetwork_V2_01-02-2018/Covariates/PCA/predictions/go_P_predictions_bonSigTerms.txt");
		
		File reactomePredictionMatrixFile = new File("/groups/umcg-wijmenga/tmp04/umcg-svandam/GeneNetwork/Data31995Genes05-12-2017/GeneNetwork_V2_01-02-2018/Covariates/PCA/predictions/reactome_predictions.txt.gz");
		File reactomePannotationMatrixFile = new File("/groups/umcg-wijmenga/tmp04/umcg-svandam/GeneNetwork/Data31995Genes05-12-2017/GeneNetwork_V2_01-02-2018/Covariates/PCA/PathwayMatrix/Ensembl2Reactome_All_Levels.txt_matrix.txt");
		File reactomePsignificantTermsFile = new File("/groups/umcg-wijmenga/tmp04/umcg-svandam/GeneNetwork/Data31995Genes05-12-2017/GeneNetwork_V2_01-02-2018/Covariates/PCA/predictions/reactome_predictions_bonSigTerms.txt");
				
		File outputFile = new File("/groups/umcg-wijmenga/tmp04/umcg-svandam/GeneNetwork/Data31995Genes05-12-2017/GeneNetwork_V2_01-02-2018/Covariates/PCA/predictions/go_P_predictions_genePredictability.txt");

		
		
		
		
	}

}
