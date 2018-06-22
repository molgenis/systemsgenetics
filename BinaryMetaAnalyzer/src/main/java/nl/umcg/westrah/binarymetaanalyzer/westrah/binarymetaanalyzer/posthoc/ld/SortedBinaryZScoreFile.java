package nl.umcg.westrah.binarymetaanalyzer.westrah.binarymetaanalyzer.posthoc.ld;

import umcg.genetica.io.Gpio;
import umcg.genetica.io.bin.BinaryFile;
import umcg.genetica.io.bin.RandomAccessFile;

import java.io.File;
import java.io.IOException;

public class SortedBinaryZScoreFile extends BinaryFile {
	int size = 0;
	long filesize = 0;
	
	public SortedBinaryZScoreFile(String loc, boolean mode) throws IOException {
		super(loc, mode);
		if (!super.writeable) {
			filesize = is.readLong();
			size = is.readInt();
			File f = new File(loc);
			if (f.length() != filesize) {
				throw new IOException("Expected file size: " + filesize + ",  but found: " + f.length());
			}
		}
	}
	
	public SortedBinaryZScoreFile(String loc, boolean mode, int buffersize) throws IOException {
		super(loc, mode, buffersize);
		if (!super.writeable) {
			filesize = is.readLong();
			size = is.readInt();
			File f = new File(loc);
			if (f.length() != filesize) {
				throw new IOException("Expected file size: " + filesize + ",  but found: " + f.length());
			}
		}
	}
	
	public SortedBinaryZScoreFile(String loc, boolean mode, int buffersize, boolean useHash) throws IOException {
		super(loc, mode, buffersize, useHash);
		if (!super.writeable) {
			filesize = is.readLong();
			size = is.readInt();
			File f = new File(loc);
			if (f.length() != filesize) {
				throw new IOException("Expected file size: " + filesize + ",  but found: " + f.length());
			}
		}
	}
	
	
	public void writeHeader(int size) throws IOException {
		os.writeLong(0l);
		os.writeInt(size);
	}
	
	public void writeZ(float[] outz) throws IOException {
		for (float f : outz) {
			os.writeFloat(f);
		}
	}
	
	public void close() throws IOException {
		super.close();
		File f = new File(this.loc);
		
		long length = f.length();
		System.out.println(length + " bytes (" + Gpio.humanizeFileSize(length) + ") written to " + f.getAbsolutePath());
		RandomAccessFile rf = new RandomAccessFile(f, "rw");
		rf.seek(0);
		rf.writeLong(length);
		rf.close();
	}
}
