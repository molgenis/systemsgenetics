/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.math.matrix2;

import java.io.BufferedInputStream;
import java.io.DataInputStream;
import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.LinkedHashMap;
import java.util.zip.GZIPInputStream;
import org.apache.commons.io.input.RandomAccessFileInputStream;

/**
 *
 * @author patri
 */
public class DoubleMatrixDatasetRowCompressedReader {

	private final LinkedHashMap<String, Integer> rowMap;
	private final LinkedHashMap<String, Integer> colMap;
	private final int numberOfRows;
	private final int numberOfColumns;

	public DoubleMatrixDatasetRowCompressedReader(String path) throws IOException {

		if (path.endsWith(".datg")) {
			path = path.substring(0, path.length() - 5);
		}

		final File matrixFile = new File(path + ".datg");
		final File rowFile = new File(path + ".rows.txt.gz");
		final File colFile = new File(path + ".cols.txt.gz");

		rowMap = DoubleMatrixDataset.loadIdentifiers(rowFile);
		colMap = DoubleMatrixDataset.loadIdentifiers(colFile);

		RandomAccessFile matrixFileReader = new RandomAccessFile(matrixFile, "r");

		matrixFileReader.seek(matrixFileReader.length() - DoubleMatrixDatasetRowCompressedWriter.APPENDIX_BYTE_LENGTH);

		numberOfRows = matrixFileReader.readInt();
		numberOfColumns = matrixFileReader.readInt();

		final long startOfEndBlock = matrixFileReader.readLong();

		if (matrixFileReader.readByte() != DoubleMatrixDatasetRowCompressedWriter.MAGIC_BYTES[0]
				|| matrixFileReader.readByte() != DoubleMatrixDatasetRowCompressedWriter.MAGIC_BYTES[1]
				|| matrixFileReader.readByte() != DoubleMatrixDatasetRowCompressedWriter.MAGIC_BYTES[2]
				|| matrixFileReader.readByte() != DoubleMatrixDatasetRowCompressedWriter.MAGIC_BYTES[3]) {
			throw new IOException("Error parsing datg file.");
		}

		System.out.println(numberOfRows);
		System.out.println(numberOfColumns);
		System.out.println();
		
		matrixFileReader.seek(startOfEndBlock);
		
		DataInputStream test = new DataInputStream(new GZIPInputStream(new BufferedInputStream(new RandomAccessFileInputStream(matrixFileReader, false), 262144)));
		System.out.println(test.readLong());
		System.out.println(test.readLong());
		
		test.close();
		
		matrixFileReader.seek(31);
		
		test = new DataInputStream(new GZIPInputStream(new BufferedInputStream(new RandomAccessFileInputStream(matrixFileReader, false), 262144)));
		System.out.println(test.readDouble());
		System.out.println(test.readDouble());
		System.out.println(test.readDouble());
		
		
		matrixFileReader.seek(0);
		
		test = new DataInputStream(new GZIPInputStream(new BufferedInputStream(new RandomAccessFileInputStream(matrixFileReader, false), 262144)));
		System.out.println(test.readDouble());
		System.out.println(test.readDouble());
		System.out.println(test.readDouble());

		
		

//		DataInputStream test2 = new DataInputStream(new GZIPInputStream(new BufferedInputStream(test)));
//		System.out.println(test2.readDouble());
//		System.out.println(test2.readDouble());
//		System.out.println(test2.readDouble());
//		System.out.println(test2.readDouble());
//		System.out.println(test2.readDouble());
//		System.out.println(test2.readDouble());
//		
	}

}
