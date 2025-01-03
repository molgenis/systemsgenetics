/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.gadocommandline;

import java.io.IOException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

/**
 *
 * @author patri
 */
public class Converters {
	
	private static final Logger LOGGER = LogManager.getLogger(Converters.class);
	
	public static void convertTxtToBin(GadoOptions options) throws IOException, Exception {

		DoubleMatrixDataset<String, String> matrix = DoubleMatrixDataset.loadDoubleTextData(options.getPredictionMatrixFile().getAbsolutePath(), '\t');
			
		matrix.saveBinary(options.getOutputBasePath());
		
		LOGGER.info("Converted matrix with " + matrix.columns() + " columns and " + matrix.rows() + " rows." );
		
	}

	
}
