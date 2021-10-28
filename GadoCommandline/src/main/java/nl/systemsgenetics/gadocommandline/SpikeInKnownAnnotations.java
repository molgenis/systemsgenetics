/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.gadocommandline;

import java.io.File;
import java.io.IOException;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

/**
 *
 * @author patri
 */
public class SpikeInKnownAnnotations {

	/**
	 * @param options
	 * @throws java.io.IOException
	 */
	public static void spikeIn(GadoOptions options) throws IOException, Exception {
		
		final File hpoPredictionMatrixFile = options.getPredictionMatrixFile();
		final File hpoAnnotationMatrixFile = options.getAnnotationMatrixFile();
		final int minZknownAnnotations = 3;
		
		DoubleMatrixDataset<String, String> hpoPredictionMatrix = DoubleMatrixDataset.loadDoubleBinaryData(hpoPredictionMatrixFile.getAbsolutePath());
		DoubleMatrixDataset<String, String> hpoAnnotationMatrix = DoubleMatrixDataset.loadDoubleData(hpoAnnotationMatrixFile.getAbsolutePath());
		
		if (!hpoPredictionMatrix.getColObjects().equals(hpoAnnotationMatrix.getColObjects())) {
			System.err.println("Differnce in terms");
			return;
		}

		if (!hpoPredictionMatrix.getRowObjects().equals(hpoAnnotationMatrix.getRowObjects())) {
			System.err.println("Differnce in genes");
			return;
		}
		
		for(int i = 0; i < hpoAnnotationMatrix.rows() ; ++i){
			for(int j = 0; j < hpoAnnotationMatrix.columns(); ++j){
				if(hpoAnnotationMatrix.getElementQuick(i, j) > 0 && hpoPredictionMatrix.getElementQuick(i, j) < minZknownAnnotations){
					hpoPredictionMatrix.setElementQuick(i, j, minZknownAnnotations);
				}
			}
		}
		
		hpoPredictionMatrix.saveBinary(options.getOutputBasePath());
		
	}
	
}
