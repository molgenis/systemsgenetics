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

    /**
     * Creates a new instance of WGAFileMatrixRawData
     */
    public WGAFileMatrixRawData(int nrSNPs, int nrInds, File fileName, boolean readOnly) throws IOException {
        this.nrSNPs = nrSNPs;
        this.nrInds = nrInds;
        if (readOnly) {
            file = new RandomAccessFile(fileName, "r");
        } else {
            file = new RandomAccessFile(fileName, "rw");
        }

        long fileSize = (long) 3 * nrSNPs * (long) nrInds;
		
		if(fileSize != file.length()){
			throw new RuntimeException("Raw datafile incorrect size. Expected: " + fileSize + " found: " + file.length());
		}
		
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
    }

    private long getElement(int snp, int ind) {
        return 3 * (long) snp * (long) nrInds + (long) ind;
    }

    public void close() throws IOException {
        file.close();
    }

    public byte getGCScore(int snp, int ind) throws IOException {

        file.seek(getElement(snp, ind));
        return file.readByte();

    }

    public byte getR(int snp, int ind) throws IOException {
        file.seek(getElement(snp, ind) + nrInds);
        return file.readByte();
    }

    public byte getTheta(int snp, int ind) throws IOException {
        file.seek(getElement(snp, ind) + nrInds * 2);
        return file.readByte();

    }

    public void setGCScore(int snp, int ind, byte[] value) throws IOException {
        
            file.seek(getElement(snp, ind));
            file.write(value);
        
    }

    public void setR(int snp, int ind, byte[] value) throws IOException {
            file.seek(getElement(snp, ind) + nrInds);
            file.write(value);
    }

    public void setTheta(int snp, int ind, byte[] value) throws IOException {
        
            file.seek(getElement(snp, ind) + nrInds * 2);
            file.write(value);
        
    }
}
