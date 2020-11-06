/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.downstreamer.development;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import ch.unil.genescore.vegas.Davies;
import ch.unil.genescore.vegas.DaviesBigDecimal;
import java.util.Arrays;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

/**
 *
 * @author patri
 */
public class TestBigDecimal {

	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) throws Exception {

		DoubleMatrixDataset<String, String> eigenMatrix = DoubleMatrixDataset.loadDoubleTextData("C:\\UMCG\\Genetica\\Projects\\Depict2Pgs\\Height_v24_special\\Height_ENSG00000197959_eigenValues.txt", '\t');

		DoubleMatrix1D eigenValues = eigenMatrix.viewCol(0).viewFlip();

		final long eigenValuesLenght = eigenValues.size();

		//Method below if from PASCAL to select relevant eigen values
		double sumPosEigen = 0;
		for (int i = 0; i < eigenValuesLenght; i++) {
			double e = eigenValues.getQuick(i);
			if (e > 0) {
				sumPosEigen += e;
			}
		}

		final double cutoff = sumPosEigen / 10000;//Only use components that explain significant part of variantion

		int eigenValuesToUse = 0;

		for (int i = 0; i < eigenValuesLenght; i++) {
			sumPosEigen -= eigenValues.getQuick(i);
			eigenValuesToUse++;

			if (sumPosEigen < cutoff) {
				break;
			}
		}

		double[] lambdas = eigenValues.viewPart(0, eigenValuesToUse).toArray();

		System.out.println(Arrays.toString(lambdas));

		DaviesBigDecimal f2 = new DaviesBigDecimal(lambdas);

		System.out.println(f2.probQsupx(1400));
		System.out.println("Error: " + f2.getIfault());

		System.out.println("-----------");

		Davies d1 = new Davies(lambdas);

		System.out.println(d1.probQsupx(1400));
		System.out.println("Error: " + d1.getIfault());

	}

}
