/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.downstreamer.io;

import com.google.common.io.Files;
import java.io.File;
import java.io.IOException;
import java.util.Random;
import nl.systemsgenetics.downstreamer.DownstreamerOptions;
import org.apache.commons.math3.util.Precision;
import org.apache.log4j.Logger;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

/**
 *
 * @author patri
 */
public class DatToDatg {
	
	private static final Logger LOGGER = Logger.getLogger(DatToDatg.class);
	
	public static void convert(DownstreamerOptions options) throws IOException{
		
		String inputMatrix = options.getGwasZscoreMatrixPath();
		
		if (inputMatrix.endsWith(".dat")) {
			inputMatrix = inputMatrix.substring(0, inputMatrix.length() - 4);
		}
		
		File originalDat = new File(inputMatrix + ".dat");		
		File originalRow = new File(inputMatrix + ".rows.txt");
		File originalCol = new File(inputMatrix + ".cols.txt");
		
		LOGGER.info("Original " + originalDat.getAbsolutePath());
		
		File workdir = new File(inputMatrix).getParentFile();
		File tmpDir = new File(workdir, "tmpdir_" + String.valueOf(Math.abs(new Random().nextInt())));
		tmpDir.mkdir();
		
		LOGGER.info("Tmp dir: " + tmpDir.getAbsolutePath());
		
		
		File originalDatTmp = new File(tmpDir, originalDat.getName());		
		File originalRowTmp = new File(tmpDir, originalRow.getName());
		File originalColTmp = new File(tmpDir,originalCol.getName());
		
		Files.move(originalDat, originalDatTmp);
		Files.move(originalRow, originalRowTmp);
		Files.move(originalCol, originalColTmp);
		
		DoubleMatrixDataset<String, String> data = DoubleMatrixDataset.loadDoubleBinaryData(originalDatTmp.getAbsolutePath());
		
		data.saveBinary(inputMatrix);
		
		DoubleMatrixDataset<String, String> newData = DoubleMatrixDataset.loadDoubleBinaryData(inputMatrix);
		
		compareTwoMatrices(data, newData,0);//This will throw IO exception if not equal
		
		LOGGER.info("New file is identical to original");
		
		originalDatTmp.delete();
		originalRowTmp.delete();
		originalColTmp.delete();
		
		if(tmpDir.listFiles().length == 0){
			tmpDir.delete();
		}
		
	}
	
	public static void compareTwoMatrices(DoubleMatrixDataset<String, String> m1, DoubleMatrixDataset<String, String> m2) throws IOException {

        compareTwoMatrices(m1, m2, 0.00000001);

    }

    public static void compareTwoMatrices(DoubleMatrixDataset<String, String> m1, DoubleMatrixDataset<String, String> m2, double delta) throws IOException {

        if(m1.rows() != m2.rows()){
			throw new IOException("Rows not equal");
		}
        if(m1.columns()!= m2.columns()){
			throw new IOException("Cols not equal");
		}

		if(!m1.getRowObjects().equals(m2.getRowObjects())){
			throw new IOException("Row names not equal");
		}
		
		if(!m1.getColObjects().equals(m2.getColObjects())){
			throw new IOException("Col names not equal");
		}

        for (int r = 0; r < m1.rows(); ++r) {
            for (int c = 0; c < m1.columns(); ++c) {
				if(!Precision.equalsIncludingNaN(m1.getElementQuick(r, c), m2.getElementQuick(r, c), delta)){
					throw new IOException("Difference at r: " + r + " c: " + c);
				}
                
            }
        }

    }

	
}
