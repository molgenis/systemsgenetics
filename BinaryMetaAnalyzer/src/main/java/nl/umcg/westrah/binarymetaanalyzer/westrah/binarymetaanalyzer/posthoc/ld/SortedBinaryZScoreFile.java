package nl.umcg.westrah.binarymetaanalyzer.westrah.binarymetaanalyzer.posthoc.ld;

import umcg.genetica.io.Gpio;
import umcg.genetica.io.bin.BinaryFile;
import umcg.genetica.io.bin.RandomAccessFile;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.File;
import java.io.IOException;

public class SortedBinaryZScoreFile extends BinaryFile {
	int size = 0;
	long filesize = 0;
	private TextFile rowdata;
	
	public SortedBinaryZScoreFile(String loc, boolean mode) throws IOException {
		super(loc, mode);
		if (!super.writeable) {
			filesize = is.readLong();
			size = is.readInt();
			File f = new File(loc);
			if (f.length() != filesize) {
				throw new IOException("Expected file size: " + filesize + ",  but found: " + f.length());
			}
			openrowdata();
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
			openrowdata();
		}
	}
	
	public SortedBinaryZScoreFile(File datasetloc, boolean r) throws IOException {
		this(datasetloc.getAbsolutePath(), r);
	}
	
	private void openrowdata() throws IOException {
		String rowfile = loc.replaceAll("-data.dat", "-combos.txt.gz");
		rowdata = new TextFile(rowfile, TextFile.R);
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
			openrowdata();
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
	
	public SortedBinaryZDataBlock readNextBlock() throws IOException {
		
		String[] elems = rowdata.readLineElems(TextFile.tab);
		if (elems == null) {
			return null;
		} else {
			
			String snp = Strings.cache(elems[1]);
			String gene = Strings.cache(elems[0]);
			String allele = Strings.cache(elems[2]);
			String assessed = Strings.cache(elems[3]);
			
			float[] data = new float[size];
			for (int i = 0; i < data.length; i++) {
				data[i] = is.readFloat();
			}
			
			SortedBinaryZDataBlock b = new SortedBinaryZDataBlock();
			
			b.snp = snp;
			b.gene = gene;
			b.allele = allele;
			b.assessed = assessed;
			b.z = data;
			return b;
		}
		
		
	}
	
	
	public void close() throws IOException {
		if (!writeable) {
			super.close();
			rowdata.close();
		} else {
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
	
	public String getName() {
		File f = new File(this.loc);
		String name = f.getName().replaceAll("-data.dat", "");
		return name;
	}
}
