/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.io.trityper.bin;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.EOFException;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;

/**
 *
 * @author harm-jan
 */
public class BinaryResultSNPSummary {

    private int counter = 0;
    private DataOutputStream out;
    private DataInputStream in;
    private int maxNrSamples;
    public static boolean W = true;
    public static boolean R = false;

    public BinaryResultSNPSummary(String filename, boolean W) throws IOException {
        if (W) {
            out = new DataOutputStream(new FileOutputStream(new File(filename + ".SNPSummary.dat")));
        } else {
            in = new DataInputStream(new FileInputStream(new File(filename)));
        }
    }

    public void close() throws IOException {
        if (out != null) {
            out.flush();
            out.close();
        }

        if (in != null) {
            in.close();
        }
    }

    public void write(String snpname, byte snpchr, Integer snpchrpos, double HWE, double MAF, double CR, byte[] alleles, byte minorallele, byte alleleassessed, Integer numsamples, long gzipIndex) throws IOException {
        out.writeInt(counter);
        out.writeUTF(snpname);
        out.writeByte(snpchr);
        out.writeInt(snpchrpos);
        out.writeDouble(HWE);
        out.writeDouble(MAF);
        out.writeDouble(CR);
        for (int i = 0; i < 2; i++) {
            out.writeByte(alleles[i]);
        }
        out.writeByte(minorallele);
        out.writeByte(alleleassessed);
        out.writeInt(numsamples);
        out.writeLong(gzipIndex);
        counter++;
    }

    public BinaryResultSNP[] readAllSNPs() throws IOException {
        ArrayList<BinaryResultSNP> snps = new ArrayList<BinaryResultSNP>();
        BinaryResultSNP snp = readNextSNP();

        int ct = 0;
        while (snp != null) {
            snps.add(snp);
            snp = readNextSNP();
            ct++;
        }

        BinaryResultSNP[] snplist = new BinaryResultSNP[snps.size()];
        maxNrSamples = 0;
        for (int s = 0; s < snplist.length; s++) {
            snplist[s] = snps.get(s);
            if (snplist[s].getNumsamples() > maxNrSamples) {
                maxNrSamples = snplist[s].getNumsamples();
            }

        }
        return snplist;
    }

    public BinaryResultSNP readNextSNP() throws IOException {
        BinaryResultSNP s = null;
        try {
            byte[] alleles = new byte[2];
            s = new BinaryResultSNP();
            s.setId(in.readInt());
            s.setName(in.readUTF());
            s.setChr(in.readByte());
            s.setChrpos(in.readInt());
            s.setHwe(in.readDouble());
            s.setMaf(in.readDouble());
            s.setCr(in.readDouble());

            alleles[0] = in.readByte();
            alleles[1] = in.readByte();
            s.setAlleles(alleles);
            s.setMinorAllele(in.readByte());
            s.setAssessedAllele(in.readByte());
            s.setNumsamples(in.readInt());
            s.setzScoreIndex(in.readLong());

            return s;
        } catch (EOFException e) {
            return null;
        }

    }

    public int getMaxNrSamples() {
        return maxNrSamples;
    }
}
