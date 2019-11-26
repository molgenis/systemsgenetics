package umcg.genetica.math.matrix2;

import org.apache.commons.io.input.CountingInputStream;

import java.io.*;
import java.nio.Buffer;
import java.nio.ByteBuffer;
import java.nio.DoubleBuffer;
import java.nio.channels.Channels;
import java.nio.channels.FileChannel;
import java.util.ArrayList;
import java.util.LinkedHashMap;

public class DoubleMatrixDatasetRandomAccessReader<R extends Comparable, C extends Comparable> {

	private final File path;
	private final long headerLen;
	private DataInputStream is;
	private CountingInputStream counter;
	private final FileChannel channel;
	private final long filesize;

	int nrRows;
	int nrCols;
	int bytesPerRow;
	LinkedHashMap<R, Integer> hashRows;
	LinkedHashMap<C, Integer> hashCols;
	long currentPos;

	int buffersize = 32 * 1024;
	private ByteBuffer bytebuffer;


	public DoubleMatrixDatasetRandomAccessReader(String filename) throws IOException {

		this.path = new File(filename + ".dat");
		channel = new RandomAccessFile(path, "r").getChannel();
		filesize = channel.size();
		counter = new CountingInputStream(new BufferedInputStream(Channels.newInputStream(channel), buffersize));
		is = new DataInputStream(counter);

		nrRows = is.readInt();
		nrCols = is.readInt();

		hashRows = (LinkedHashMap<R, Integer>) DoubleMatrixDataset.loadIdentifiers(filename + ".rows.txt");
		hashCols = (LinkedHashMap<C, Integer>) DoubleMatrixDataset.loadIdentifiers(filename + ".cols.txt");

		headerLen = 8;
		bytesPerRow = 8 * nrCols;
		buffersize = bytesPerRow * 10;
		currentPos = headerLen;
		System.out.println(filename + " is " + nrRows + " rows and " + nrCols + " cols; bytes per row: " + bytesPerRow);
	}


	public double[] getNextRow(double[] buffer) throws IOException {
		for (int d = 0; d < nrCols; d++) {
			buffer[d] = is.readDouble();
		}
		currentPos += bytesPerRow;
		return buffer;
	}

	public double[] getNextRow() throws IOException {
		double[] output = new double[nrCols];
		return getNextRow(output);
	}

	public double get(int row, int col) throws IOException {
		return getRow(row)[col];
	}

	public double[] getRow(int row) throws IOException {
		return getRow(row, new double[nrCols]);
	}

	public double[] getRow(int row, double[] buffer) throws IOException {
		long seekLoc = ((long) row * bytesPerRow) + headerLen;

//		System.out.println(row + "\t" + seekLoc + "\t" + bytesPerRow + "\t" + headerLen);
		if (seekLoc > filesize) {
			throw new IllegalArgumentException("Seek location for row: " + row + ", " + seekLoc + " is outside file size: " + filesize);
		}
		// if row is the next row, don't use random access..
		if (seekLoc - currentPos == 0) {
//			System.out.println("Getting using channel");
			return getNextRow(buffer);
		}
//		System.out.println("Getting using RAF");

		// else use random access
		channel.position(seekLoc);
		if (bytebuffer == null) {
//			System.out.println("Creating new bytebuffer: " + bytesPerRow);
			bytebuffer = ByteBuffer.allocate(bytesPerRow);
		} else {
//			System.out.println("Resetting bytebyffer: " + bytesPerRow);
			((Buffer) bytebuffer).clear();
		}

		int nrbytesread = channel.read(bytebuffer);
		if (nrbytesread < bytesPerRow) {
			System.out.println("Error: read " + nrbytesread + " while expected " + bytesPerRow);
			System.exit(-1);
		}
		currentPos = seekLoc + bytesPerRow;


		// TODO: change the position of the datainputstream? this is probably very very slow...
		// question is whether this is needed when the position of the channel underlying a datainputstream is changed?
		// does DataInputStream have a counter of it's own?
		counter = new CountingInputStream(new BufferedInputStream(Channels.newInputStream(channel), buffersize));
		is = new DataInputStream(counter);
		double[] output = new double[nrCols];


//		System.out.println("position: " + bytebuffer.position() + "\tremaining: " + bytebuffer.remaining());
		((Buffer) bytebuffer).flip();
		for (int i = 0; i < nrCols; i++) {
			output[i] = bytebuffer.getDouble();
//			System.out.println(i + "\t" + output[i]);
		}

		return output;

	}

	public void close() throws IOException {
		this.is.close();
		this.counter.close();
		this.channel.close();
	}

	public ArrayList<C> getColObjects() {
		return new ArrayList<C>(hashCols.keySet());
	}

	public ArrayList<R> getRowObjects() {
		return new ArrayList<>(hashRows.keySet());
	}

	public LinkedHashMap<R, Integer> getHashRows() {
		return hashRows;
	}

	public LinkedHashMap<C, Integer> getHashCols() {
		return hashCols;
	}

	public int rows() {
		return nrRows;
	}

	public int cols() {
		return nrCols;
	}
}
