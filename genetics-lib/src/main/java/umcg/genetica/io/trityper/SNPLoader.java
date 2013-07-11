/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.io.trityper;

import java.io.IOException;
import java.io.RandomAccessFile;



/**
 *
 * @author harm-jan
 */
public class SNPLoader {

    private RandomAccessFile m_genotypehandle;
    private RandomAccessFile m_dosagehandle;
    private int m_numIndividuals;
    private final Boolean[] m_isIncluded, m_isFemale;

    public SNPLoader(RandomAccessFile genotypehandle, Boolean[] indIsIncluded, Boolean[] isFemale) {
        m_genotypehandle = genotypehandle;
        m_isIncluded = indIsIncluded;
        m_isFemale = isFemale;
    }

    public SNPLoader(RandomAccessFile genotypehandle, RandomAccessFile dosagehandle, Boolean[] indIsIncluded, Boolean[] isFemale) {
        m_genotypehandle = genotypehandle;
        m_dosagehandle = dosagehandle;
        m_isIncluded = indIsIncluded;
        m_isFemale = isFemale;
    }

    public void loadGenotypes(SNP snp) throws IOException {
        byte[] allele1 = new byte[m_numIndividuals];
        byte[] allele2 = new byte[m_numIndividuals];

        int bytesize = m_numIndividuals * 2;

        long seekLoc = (long) snp.getId() * (long) bytesize;

        byte[] alleles = new byte[bytesize];


        if (m_genotypehandle.getFilePointer() != seekLoc) {
            m_genotypehandle.seek(seekLoc);
        }

        m_genotypehandle.read(alleles, 0, bytesize);

        System.arraycopy(alleles, 0, allele1, 0, m_numIndividuals);
        System.arraycopy(alleles, m_numIndividuals, allele2, 0, m_numIndividuals);

        alleles = null;

        snp.setAlleles(allele1, allele2, m_isIncluded, m_isFemale);

    }

    public void loadDosage(SNP snp) throws IOException {
        if (m_dosagehandle != null) {
            byte[] dosageValues = new byte[m_numIndividuals];
            //if (loadedSNP.getGcScores()==null||loadedSNP.getThetaValues()==null||loadedSNP.getRValues()==null) {
            long seekLoc = (long) snp.getId() * (long) m_numIndividuals * 1;

            m_dosagehandle.seek(seekLoc);
            m_dosagehandle.read(dosageValues, 0, m_numIndividuals);

            byte[] genotypes = snp.getGenotypes();

            boolean takeComplement = false;
            for (int ind = 0; ind < dosageValues.length; ind++) {
                double dosagevalue = ((double) (-Byte.MIN_VALUE + dosageValues[ind])) / 100;
                if (genotypes[ind] == 0 && dosagevalue > 1) {
                    takeComplement = true;
                    break;
                }
                if (genotypes[ind] == 2 && dosagevalue < 1) {
                    takeComplement = true;
                    break;
                }
            }
            if (takeComplement) {
                for (int ind = 0; ind < dosageValues.length; ind++) {
                    byte dosageValue = (byte) (200 - (-Byte.MIN_VALUE + dosageValues[ind]) + Byte.MIN_VALUE);
                    dosageValues[ind] = dosageValue;
                }
            }

            snp.setDosage(dosageValues);
        }
    }

    /**
     * @return the numIndividuals
     */
    public int getNumIndividuals() {
        return m_numIndividuals;
    }

    /**
     * @param numIndividuals the numIndividuals to set
     */
    public void setNumIndividuals(int numIndividuals) {
        this.m_numIndividuals = numIndividuals;
    }

    public boolean hasDosageInformation() {
        return (m_dosagehandle != null);
    }

    public double getAverageSNPSize(int numSNPs) throws IOException {
        long size = 0;

        size += m_genotypehandle.length();
        if (m_dosagehandle != null) {
            size += m_dosagehandle.length();
        }


        double avgSNPSize = 0;
        if (size > 0) {
            avgSNPSize = (double) size / numSNPs;
        }

        return avgSNPSize;
    }

    public void close() throws IOException {
        if (m_dosagehandle != null) {
            m_dosagehandle.close();
        }
        m_genotypehandle.close();
    }
}
