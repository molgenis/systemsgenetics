package eqtlmappingpipeline.util;

import umcg.genetica.math.matrix2.DoubleMatrixDataset;

import java.io.IOException;

public class BinaryMatrixConverter {

	public void run(String input, String output) {

		try {
			DoubleMatrixDataset<String, String> ds = DoubleMatrixDataset.loadDoubleBinaryData(input);
			ds.save(output);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}
