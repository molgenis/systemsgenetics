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
//import umcg.genetica.io.bin.RandomAccessFile;

/**
 *
 * @author ludefranke
 */
public class WGAFileMatrixImputedDosage {

    public int nrSNPs = 0;
    public int nrInds = 0;
    private RandomAccessFile file = null;

    /**
     * Creates a new instance of WGAFileMatrixGenotype
     */
    public WGAFileMatrixImputedDosage(int nrSNPs, int nrInds, File fileName, boolean readOnly) throws IOException {
        this.nrSNPs = nrSNPs;
        this.nrInds = nrInds;
        if (readOnly) {
            System.out.println("Opening imputed dosage matrix file: " + fileName);
            file = new RandomAccessFile(fileName, "r");
        } else {
            file = new RandomAccessFile(fileName, "rw");
        }

        long fileSize = (long) 1 * nrSNPs * (long) nrInds;
        if (!readOnly) {
            if (file.length() != fileSize) {
                System.out.println("Creating imputed dosage matrix for " + nrSNPs + " SNPS and " + nrInds + " individuals. Eventual size: " + fileSize);
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
                    if (x % 134217728 == 0) {
                        double perc = Math.floor((double) x / fileSize) * 100;
                        System.out.println("Written " + perc + "%");
                    }
                }
                System.out.println("Size matrix:\t" + fileSize + "\tFile size:\t" + file.length());
            }
        }
    }

    private long getElement(int snp, int ind) {
        return 1 * (long) snp * (long) nrInds + (long) ind;
    }

    public void close() throws IOException {
        file.close();
    }

    public byte getDosage(int snp, int ind) throws IOException {
        file.seek(getElement(snp, ind));
        return file.readByte();

    }

    public void setDosage(int snp, int ind, byte[] value) throws IOException {

        file.seek(getElement(snp, ind));
        file.write(value);

    }

    public void setDosage(int snp, int ind, byte value) throws IOException {

        file.seek(getElement(snp, ind));
        file.write(value);

    }
}
