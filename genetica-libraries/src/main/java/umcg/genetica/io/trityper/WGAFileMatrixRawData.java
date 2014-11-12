/*
 * WGAFileMatrixRawData.java
 *
 * Created on July 10, 2007, 11:15 PM
 *
 */
package umcg.genetica.io.trityper;

import java.io.File;
import java.io.IOException;
import umcg.genetica.io.bin.RandomAccessFile;


/**
 *
 * @author ludefranke
 */
public class WGAFileMatrixRawData {

    public int nrSNPs = 0;
    public int nrInds = 0;
    private RandomAccessFile file = null;
	
	private final byte[] snpCache;
	private int cachedSnp = -1;

    /**
     * Creates a new instance of WGAFileMatrixRawData
     */
    public WGAFileMatrixRawData(int nrSNPs, int nrInds, File fileName, boolean readOnly) throws IOException {
        this.nrSNPs = nrSNPs;
        this.nrInds = nrInds;
		snpCache = new byte[nrInds * 3];
        if (readOnly) {
            file = new RandomAccessFile(fileName, "r");
        } else {
            file = new RandomAccessFile(fileName, "rw");
        }

        long fileSize = (long) 3 * nrSNPs * (long) nrInds;
		
        if (!readOnly) {
            if (file.length() != fileSize) {
                //Generate file with the size, such that this is appropriate:
//                file.setLength(fileSize);
                file.seek(0);
                byte byteString[] = new byte[1000];
                for (int g = 0; g < 1000; g++) {
                    byteString[g] = 0;
                }
                for (long x = 0; x < fileSize - 1000; x += 1000) {
                    file.write(byteString);
                }
                long remainder = fileSize % 1000;
                byte byteSingle[] = new byte[1];
                byteSingle[0] = 0;
                for (long x = 0; x < remainder; x++) {
                    file.write(byteSingle);
                }
                System.out.println("Size matrix:\t" + fileSize + "\tFile size:\t" + file.length());
            }
        }
		
		if(fileSize != file.length()){
			throw new RuntimeException("Raw datafile incorrect size. Expected: " + fileSize + " found: " + file.length());
		}
    }

    private long getElement(int snp, int ind) {
        return 3 * (long) snp * (long) nrInds + (long) ind;
    }
	
	private void readSnpIntoCache(int snp) throws IOException{
		file.seek(3 * (long) snp * (long) nrInds);
		int readBytes = file.read(snpCache);
		if(readBytes != nrInds * 3){
			throw new IOException("RawMatrix read error. expected: " + snpCache.length + " read: " + readBytes + " snp index: " + snp + " file index: " + (3 * (long) snp * (long) nrInds));
		}
		cachedSnp = snp;
	}

    public void close() throws IOException {
        file.close();
    }

    public byte getGCScore(int snp, int ind) throws IOException {
		if(cachedSnp != snp){
			readSnpIntoCache(snp);
		}
        return snpCache[ind];
    }

    public byte getR(int snp, int ind) throws IOException {
        if(cachedSnp != snp){
			readSnpIntoCache(snp);
		}
        return snpCache[ind + nrInds];
    }

    public byte getTheta(int snp, int ind) throws IOException {
        if(cachedSnp != snp){
			readSnpIntoCache(snp);
		}
        return snpCache[ind + nrInds + nrInds];

    }

    public void setGCScore(int snp, int ind, byte[] value) throws IOException {
			cachedSnp = -1;
            file.seek(getElement(snp, ind));
            file.write(value);
        
    }

    public void setR(int snp, int ind, byte[] value) throws IOException {
			cachedSnp = -1;
            file.seek(getElement(snp, ind) + nrInds);
            file.write(value);
    }

    public void setTheta(int snp, int ind, byte[] value) throws IOException {
			cachedSnp = -1;
            file.seek(getElement(snp, ind) + nrInds * 2);
            file.write(value);
        
    }
}
