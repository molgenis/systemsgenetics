/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.depict2.development;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import java.io.IOException;
import java.util.ArrayList;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.stats.ZScores;

/**
 *
 * @author patri
 */
public class ExploreCoregulationMatrix {

	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) throws IOException {

		String coregFile = "C:\\UMCG\\Genetica\\Projects\\Depict2Pgs\\UKBIO_neale_lab_batch1_B_50\\Coregulation_Enrichment_zscoreExHla";

		DoubleMatrixDataset<String, String> coregData = DoubleMatrixDataset.loadDoubleBinaryData(coregFile);

		final ArrayList<String> traits = coregData.getColObjects();
		final ArrayList<String> genes = coregData.getRowObjects();

		final int geneCount = coregData.rows();

		double bonnP = 0.05 / geneCount;
		double zscoreMin = ZScores.pToZTwoTailed(bonnP);
		double zscorePlus = zscoreMin * -1;

		System.out.println("Bonnferoni p-value:" + bonnP + " zscore: " + zscorePlus);

		System.out.println("Trait\tPosSignifcantGeneCount\tNegSignificantGeneCount");
		for (int phenoI = 0; phenoI < coregData.columns(); ++phenoI) {

			int posSig = 0;
			int negSig = 0;

			for (int geneI = 0; geneI < geneCount; ++geneI) {
				double z = coregData.getElementQuick(geneI, phenoI);
				if (z >= zscorePlus) {
					++posSig;
				} else if (z <= zscoreMin) {
					++negSig;
				}
			}

			System.out.println(traits.get(phenoI) + "\t" + posSig + "\t" + negSig);

		}

		System.out.println("----------");

		for (int geneI = 0; geneI < geneCount; ++geneI) {

			int posSig = 0;
			int negSig = 0;

			for (int phenoI = 0; phenoI < coregData.columns(); ++phenoI) {

				double z = coregData.getElementQuick(geneI, phenoI);
				if (z >= zscorePlus) {
					++posSig;
				} else if (z <= zscoreMin) {
					++negSig;
				}
			}

			System.out.println(genes.get(geneI) + "\t" + posSig + "\t" + negSig);

		}

	}

}
