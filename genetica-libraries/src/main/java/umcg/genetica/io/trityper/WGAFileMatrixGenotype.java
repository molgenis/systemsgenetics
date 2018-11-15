/*
 * WGAFileMatrixGenotype.java
 *
 * Created on July 9, 2007, 5:08 PM
 *
 */
package umcg.genetica.io.trityper;

import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;

import umcg.genetica.console.ProgressBar;
//import umcg.genetica.io.bin.RandomAccessFile;

/**
 * @author ludefranke
 */
public class WGAFileMatrixGenotype {
	
	public int nrSNPs = 0;
	public int nrInds = 0;
	private RandomAccessFile file = null;
	
	/**
	 * Creates a new instance of WGAFileMatrixGenotype
	 */
	public WGAFileMatrixGenotype(int nrSNPs, int nrInds, File fileName, boolean readOnly) throws IOException {
		this.nrSNPs = nrSNPs;
		this.nrInds = nrInds;
		
		if (readOnly) {
			System.out.println("Opening genotype matrix file: " + fileName);
			file = new RandomAccessFile(fileName, "r");
		} else {
			file = new RandomAccessFile(fileName, "rw");
		}
		
		long fileSize = (long) 2 * nrSNPs * (long) nrInds;
		if (!readOnly) {
			
			if (file.length() != fileSize) {
				
				ProgressBar pb = new ProgressBar(fileSize, "Creating genotype matrix for " + nrSNPs + " SNPS and " + nrInds + " individuals. Eventual size: " + fileSize);
				//Generate file with the size, such that this is appropriate:
				//file.setLength(fileSize);
				file.seek(0);
				int buffersize = 32 * 1024;
				
				byte byteString[] = new byte[buffersize];
				for (int g = 0; g < buffersize; g++) {
					byteString[g] = 0;
				}
				for (long x = 0; x < fileSize - buffersize; x += buffersize) {
					file.write(byteString);
					if (x % 1048576 * 10 == 0) {
						pb.set(x);
					}
				}
				long remainder = fileSize % buffersize;
				if (remainder > 0) {
					byte[] byteSingle = new byte[(int) remainder];
					file.write(byteSingle);
				}
//                byteSingle[0] = 0;
//                for (long x = 0; x < remainder; x++) {
//                    file.write(byteSingle);
//                }
				
				pb.close();
				System.out.println("Size genotype matrix:\t" + fileSize + "\tFile size:\t" + file.length());
			}
			
		}
	}
	
	private long getElement(int snp, int ind) {
		return 2 * (long) snp * (long) nrInds + (long) ind;
	}
	
	public void close() throws IOException {
		file.close();
	}
	
	public byte getAllele1(int snp, int ind) throws IOException {
		file.seek(getElement(snp, ind));
		return file.readByte();
	}
	
	public byte getAllele2(int snp, int ind) throws IOException {
		
		file.seek(getElement(snp, ind) + nrInds);
		return file.readByte();
		
	}
	
	public void setAllele1(int snp, int ind, byte value) throws IOException {
		file.seek(getElement(snp, ind));
		file.write(value);
	}
	
	public void setAllele1(int snp, int ind, byte[] value) throws IOException {
		file.seek(getElement(snp, ind));
		file.write(value);
	}
	
	public void setAllele2(int snp, int ind, byte value) throws IOException {
		
		file.seek(getElement(snp, ind) + nrInds);
		file.write(value);
		
	}
	
	public void setAllele2(int snp, int ind, byte[] value) throws IOException {
		file.seek(getElement(snp, ind) + nrInds);
		file.write(value);
	}
	
	public void setAlleles(int snp, byte[] alleles1, byte[] alleles2) throws IOException {
		// write all alleles at once
		
		// allele1 position
		// 2 * (long) snp * (long) nrInds
		file.seek(2 * (long) snp * (long) nrInds);
		file.write(alleles1);
		
		// allele2 position
		// 2 * (long) snp * (long) nrInds + nrInds
		file.seek((2 * (long) snp * (long) nrInds) + nrInds);
		file.write(alleles2);
		
	}
}
