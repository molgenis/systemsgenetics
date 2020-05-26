/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.io.trityper.util;

import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class BaseAnnot {

    public static final byte A = 65;
    public static final byte C = 67;
    public static final byte T = 84;
    public static final byte G = 71;
    public static final byte U = 85;
    public static final byte I = 73;
    public static final byte N = 78;

    public static byte[] convertToComplementaryAlleles(byte[] allelesToCompare) {
        byte[] allelesComplementary = new byte[2];
        for (int a = 0; a < 2; a++) {
            allelesComplementary[a] = getComplement(allelesToCompare[a]);
        }
        return allelesComplementary;
    }

    // TODO: AT / GC SNPs??
    public static Boolean flipalleles(String allelesSNP1, String alleleAssessedSNP1, String allelesSNP2, String alleleAssessedSNP2) {

        String[] allelesSNP1Arr = allelesSNP1.split("/");
        byte[] allelesfirst = new byte[2];
        String[] allelesSNP2Arr = allelesSNP2.split("/");
        byte[] allelessecond = new byte[2];
        for (int i = 0; i < 2; i++) {
            allelesfirst[i] = toByte(allelesSNP1Arr[i]);
            allelessecond[i] = toByte(allelesSNP2Arr[i]);
        }

        byte allelefirstassessed = toByte(alleleAssessedSNP1);
        byte allelesecondassessed = toByte(alleleAssessedSNP2);

        int nridenticalalleles = 0;

        for (int i = 0; i < allelesfirst.length; i++) {
            byte allele1 = allelesfirst[i];
            for (int j = 0; j < allelessecond.length; j++) {
                if (allelessecond[j] == allele1) {
                    nridenticalalleles++;
                }
            }
        }

        if (nridenticalalleles == 2) {
            // alleles are identical. check if same allele was assessed...
            if (allelefirstassessed == allelesecondassessed) {
                return false;
            } else {
                return true;
            }
        } else {
            // try complement
            allelessecond = convertToComplementaryAlleles(allelessecond);
            allelesecondassessed = BaseAnnot.getComplement(allelesecondassessed);
            nridenticalalleles = 0;

            for (int i = 0; i < allelesfirst.length; i++) {
                byte allele1 = allelesfirst[i];
                for (int j = 0; j < allelessecond.length; j++) {
                    if (allelessecond[j] == allele1) {
                        nridenticalalleles++;
                    }
                }
            }

            if (nridenticalalleles == 2) {
                // alleles are identical. check if same allele was assessed...
                if (allelefirstassessed == allelesecondassessed) {
                    return false;
                } else {
                    return true;
                }
            }
        }
        return null;
    }

    public static String toString(byte x) {

        if (x == A) {
            return "A";
        }
        if (x == T) {
            return "T";
        }
        if (x == U) {
            return "U";
        }
        if (x == C) {
            return "C";
        }
        if (x == G) {
            return "G";
        }
        if (x == N) {
            return "N";
        }
        if (x == I) {
            return "I";
        }
        return "0";
    }

    public static byte getComplement(byte x) {
        if (x == A) {
            return T;
        }
        if (x == T) {
            return A;
        }
        if (x == U) {
            return A;
        }
        if (x == C) {
            return G;
        }
        if (x == G) {
            return C;
        }
        if(x == I){
            return I;
        }
        if(x == N){
            return N;
        }
        return 0;
    }

    public static String getReverseComplement(String str) {
        String reverse = Strings.reverse(str);
        String reversecomplement = getFullComplement(reverse);
        return reversecomplement;
    }

    public static String getFullComplement(String longstr) {
        StringBuilder buffer = new StringBuilder(longstr.length());
        for (int i = 0; i < longstr.length(); i++) {
            buffer.append(getComplement("" + longstr.charAt(i)));
        }
        return buffer.toString();
    }

    public static String getComplement(String x) {
        if (x.equals("A")) {
            return "T";
        }
        if (x.equals("T")) {
            return "A";
        }
        if (x.equals("U")) {
            return "A";
        }
        if (x.equals("C")) {
            return "G";
        }
        if (x.equals("G")) {
            return "C";
        }
        // inversion genotypes
        if(x.equals("N")){
            return "N";
        }
        if(x.equals("I")){
            return "I";
        }
        return "0";
    }

    public static byte toByte(String x) {
        if (x.equals("A")) {
            return A;
        }
        if (x.equals("C")) {
            return C;
        }
        if (x.equals("T")) {
            return T;
        }
        if (x.equals("G")) {
            return G;
        }
        if (x.equals("N")) {
            return N;
        }
        if (x.equals("U")) {
            return U;
        }
        // inversion genotypes
        if(x.equals("N")){
            return N;
        }
        if(x.equals("I")){
            return I;
        }
        return 0;
    }

    public static String getAllelesDescription(byte[] alleles) {
        return toString(alleles[0]) + "/" + toString(alleles[1]);
    }

    public static byte[] toByteArray(String alleles) {
        byte[] output = new byte[2];
        String[] alleleStr = alleles.split("/");
        for(int i=0;i<alleleStr.length;i++){
            output[i] = toByte(alleleStr[i]);
        }
        return output;
    }
}
