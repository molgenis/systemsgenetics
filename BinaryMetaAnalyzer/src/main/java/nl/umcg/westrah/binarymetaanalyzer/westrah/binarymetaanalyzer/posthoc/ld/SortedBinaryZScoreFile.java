package nl.umcg.westrah.binarymetaanalyzer.westrah.binarymetaanalyzer.posthoc.ld;

import umcg.genetica.io.Gpio;
import umcg.genetica.io.bin.BinaryFile;
import umcg.genetica.io.bin.RandomAccessFile;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.File;
import java.io.IOException;
import java.util.LinkedHashMap;

public class SortedBinaryZScoreFile extends BinaryFile {
	int nrElemsPerArray = 0;
	long filesize = 0;
	private TextFile rowdata;
	private LinkedHashMap<String, Integer> availableGenes = null;
	private int HEADERLEN = 12;
	
	public SortedBinaryZScoreFile(String loc, boolean mode) throws IOException {
		super(loc, mode, 32 * 1024);
		if (!super.writeable) {
			filesize = is.readLong();
			nrElemsPerArray = is.readInt();
			File f = new File(loc);
			if (f.length() != filesize) {
				throw new IOException("Expected file nrElemsPerArray: " + filesize + ",  but found: " + f.length());
			}
			openrowdata();
		}
		
		
	}
	
	public SortedBinaryZScoreFile(String loc, boolean mode, int buffersize) throws IOException {
		super(loc, mode, buffersize);
		if (!super.writeable) {
			filesize = is.readLong();
			nrElemsPerArray = is.readInt();
			File f = new File(loc);
			if (f.length() != filesize) {
				throw new IOException("Expected file nrElemsPerArray: " + filesize + ",  but found: " + f.length());
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
		
		availableGenes = new LinkedHashMap<>();
		String ln = rowdata.readLine();
		int lnctr = 0;
		while (ln != null) {
			String gene = Strings.subsplit(ln, Strings.tab, 0, 1)[0];
			if (!availableGenes.containsKey(gene)) {
				availableGenes.put(gene, lnctr);
			}
			lnctr++;
			ln = rowdata.readLine();
		}
		
		System.out.println(rowfile + " has " + availableGenes.size() + " genes");
		
		rowdata.close();
		rowdata.open();
		
	}
	
	
	public boolean hasGene(String gene) {
		return availableGenes.containsKey(gene);
	}
	
	public SortedBinaryZScoreFile(String loc, boolean mode, int buffersize, boolean useHash) throws IOException {
		super(loc, mode, buffersize, useHash);
		if (!super.writeable) {
			filesize = is.readLong();
			nrElemsPerArray = is.readInt();
			File f = new File(loc);
			if (f.length() != filesize) {
				throw new IOException("Expected file nrElemsPerArray: " + filesize + ",  but found: " + f.length());
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
	
	int currentposition = 0;
	
	public SortedBinaryZDataBlock readNextBlock() throws IOException {
		
		String[] elems = rowdata.readLineElems(TextFile.tab);
		if (elems == null) {
			return null;
		} else {
			String snp = Strings.cache(elems[1]);
			String gene = Strings.cache(elems[0]);
			String allele = Strings.cache(elems[2]);
			String assessed = Strings.cache(elems[3]);
			
			float[] data = new float[nrElemsPerArray];
			for (int i = 0; i < data.length; i++) {
				data[i] = is.readFloat();
			}
			
			SortedBinaryZDataBlock b = new SortedBinaryZDataBlock();
			
			b.snp = snp;
			b.gene = gene;
			b.allele = allele;
			b.assessed = assessed;
			b.z = data;
			b.n = Integer.parseInt(elems[4]);
			currentposition++;
			return b;
		}
	}
	
	
	public void skipTo(String gene) throws IllegalAccessException, IOException {
		
		Integer geneId = availableGenes.get(gene);
		if (geneId != null) {
			
			long currentPosB = (nrElemsPerArray * 4 * currentposition) + HEADERLEN;
			long lookuppositioninbinaryfile = ((long) nrElemsPerArray * 4 * geneId) + HEADERLEN; // nr floats * geneID + header int + header long
			
			System.out.println("Looking for gene: " + gene + "\t" + geneId + "\tcurrent: " + currentposition + "\tbpos: " + currentPosB + "\tlookup: " + lookuppositioninbinaryfile);
			long difference = lookuppositioninbinaryfile - currentPosB;
			if (difference < 0) {
				throw new IllegalAccessException("Can't skip backwards! " + loc + "\tcurrently at: " + currentPosB + "\tlooking for: " + lookuppositioninbinaryfile + "\tdiff: " + difference);
			}
			
			is.skip(difference);
			
			// skip the text file forward as well...
			// determine the number of lines to read until we hit the geneId
			int nrlinesToread = geneId - currentposition;
			for (int n = 0; n < nrlinesToread; n++) {
				rowdata.readLine();
			}
			
			currentposition = geneId;
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
