package umcg.genetica.math.matrix2;

import umcg.genetica.io.bin.BinaryFile;
import umcg.genetica.io.bin.RandomAccessFile;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.List;

public class DoubleMatrixDatasetAppendableWriter {

	private final TextFile rowIds;
	BinaryFile bf;
	String loc;
	private int rowctr = 0;

	public DoubleMatrixDatasetAppendableWriter(List<String> colids, String loc) throws IOException {
		this(colids, loc, null);
	}

	public DoubleMatrixDatasetAppendableWriter(List<String> colids, String loc, Integer buffersize) throws IOException {
		this.loc = loc;
		TextFile tf = new TextFile(loc + ".cols.txt", TextFile.W);
		for (String s : colids) {
			tf.writeln(s);
		}
		tf.close();

		if (buffersize == null) {
			buffersize = colids.size() * 8;
		}
		bf = new BinaryFile(loc + ".dat", BinaryFile.W, buffersize);

		bf.writeInt(0);
		bf.writeInt(colids.size());
		rowIds = new TextFile(loc + ".rows.txt", TextFile.W);
	}

	synchronized public void append(double[] data, String rowId) throws IOException {
		for (double d : data) {
			bf.writeDouble(d);
		}
		rowctr++;
		rowIds.writeln(rowId);
	}

	public void close() throws IOException {

		bf.close();
		rowIds.close();

		RandomAccessFile rf = new RandomAccessFile(loc + ".dat", "rw");
		rf.seek(0);
		rf.writeInt(rowctr);
		rf.close();
		System.out.println("Final matrix has: " + rowctr + " rows.");


	}


}
