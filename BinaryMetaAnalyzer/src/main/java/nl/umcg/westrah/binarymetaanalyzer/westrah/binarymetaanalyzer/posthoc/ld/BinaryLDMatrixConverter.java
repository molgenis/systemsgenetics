package nl.umcg.westrah.binarymetaanalyzer.westrah.binarymetaanalyzer.posthoc.ld;

import umcg.genetica.io.bin.BinaryFile;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

import java.io.IOException;
import java.util.ArrayList;

public class BinaryLDMatrixConverter {
	
	
	public void run(String ldmatrixin, String ldmatrixout) throws IOException {
		
		
		BinaryFile bf = new BinaryFile(ldmatrixin, BinaryFile.R);
		
		int nrsnps = bf.readInt();
		ArrayList<String> snps = new ArrayList<String>();
		for (int i = 0; i < nrsnps; i++) {
			snps.add(bf.readString());
		}
		
		/*
		for (int i = 0; i < snpList.size(); i++) {
						for (int j = i + 1; j < snpList.size(); j++) {
							float f = (float) ldmatrix.getMatrix().getQuick(i, j);
							bf.writeFloat(f);
						}
					}
		 */
		DoubleMatrixDataset<String, String> ds = new DoubleMatrixDataset<String, String>(nrsnps, nrsnps);
		for (int i = 0; i < nrsnps; i++) {
			ds.getMatrix().setQuick(i, i, 1d);
			for (int j = i + 1; j < nrsnps; j++) {
				float f = bf.readFloat();
				ds.getMatrix().setQuick(i, j, f);
				ds.getMatrix().setQuick(j, i, f);
			}
		}
		bf.close();
		
		ds.save(ldmatrixout);
		
		
	}
}
