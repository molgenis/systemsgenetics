/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.io.trityper;

import java.util.ArrayList;
import umcg.genetica.io.trityper.util.ChrAnnotation;

/**
 *
 * @author harmjan
 */
public class MinimalEQTL {

    private double pvalue = 1d;
    private String rsName;
    private byte rsChr;
    private int rsChrPos;
    private String probe;
    private byte probeChr;
    private int probeChrPos;
    private double zScore = -9;
    private char assesedAllele = ' ';

    public MinimalEQTL() {
    }
    
    public MinimalEQTL(String[] eQTLrowInformation) {
        this.pvalue = Double.parseDouble(eQTLrowInformation[0]);
        this.rsName = eQTLrowInformation[1];
        this.rsChr = ChrAnnotation.parseChr(eQTLrowInformation[2]);
        this.rsChrPos = Integer.parseInt(eQTLrowInformation[3]);
        this.probe = eQTLrowInformation[4];
        this.probeChr = ChrAnnotation.parseChr(eQTLrowInformation[5]);
        this.probeChrPos = Integer.parseInt(eQTLrowInformation[6]);
        this.zScore = Double.parseDouble(eQTLrowInformation[10]);
        this.assesedAllele = eQTLrowInformation[9].charAt(0);
    }
    
    public MinimalEQTL(double pValue, String rsName, byte rsChr, int rsChrPos,  String probe, byte probeChr, int probeChrPos) {
        this.pvalue = pValue;
        this.rsName = rsName;
        this.rsChr = rsChr;
        this.rsChrPos = rsChrPos;
        this.probe = probe;
        this.probeChr = probeChr;
        this.probeChrPos = probeChrPos;
    }
    
    public MinimalEQTL(EQTL e) {
        this.pvalue = e.getPvalue();
        this.rsName = e.getRsName();
        this.rsChr = e.getRsChr();
        this.rsChrPos = e.getRsChrPos();
        this.probe = e.getProbe();
        this.probeChr = e.getProbeChr();
        this.probeChrPos = e.getProbeChrPos();
    }

    /**
     * @return the pvalue
     */
    public double getPvalue() {
        return pvalue;
    }

    /**
     * @param pvalue the pvalue to set
     */
    public void setPvalue(double pvalue) {
        this.pvalue = pvalue;
    }

    /**
     * @return the rsName
     */
    public String getRsName() {
        return rsName;
    }

    /**
     * @param rsName the rsName to set
     */
    public void setRsName(String rsName) {
        this.rsName = rsName;
    }

    /**
     * @return the rsChr
     */
    public Byte getRsChr() {
        return rsChr;
    }

    /**
     * @param rsChr the rsChr to set
     */
    public void setRsChr(byte rsChr) {
        this.rsChr = rsChr;
    }

    /**
     * @return the rsChrPos
     */
    public Integer getRsChrPos() {
        return rsChrPos;
    }

    /**
     * @param rsChrPos the rsChrPos to set
     */
    public void setRsChrPos(int rsChrPos) {
        this.rsChrPos = rsChrPos;
    }

    /**
     * @return the probe
     */
    public String getProbe() {
        return probe;
    }

    /**
     * @param probe the probe to set
     */
    public void setProbe(String probe) {
        this.probe = probe;
    }

    /**
     * @return the probeChr
     */
    public Byte getProbeChr() {
        return probeChr;
    }

    /**
     * @param probeChr the probeChr to set
     */
    public void setProbeChr(byte probeChr) {
        this.probeChr = probeChr;
    }

    /**
     * @return the probeChrPos
     */
    public Integer getProbeChrPos() {
        return probeChrPos;
    }

    /**
     * @param probeChrPos the probeChrPos to set
     */
    public void setProbeChrPos(int probeChrPos) {
        this.probeChrPos = probeChrPos;
    }

    public Double getPvalueAbs() {
        return Math.abs(pvalue);
    }

    public double getzScore() {
        return zScore;
    }

    public void setzScore(double zScore) {
        this.zScore = zScore;
    }

    public char getAssesedAllele() {
        return assesedAllele;
    }

    public void setAssesedAllele(char assesedAllele) {
        this.assesedAllele = assesedAllele;
    }
    
    

    public static ArrayList<MinimalEQTL> convertArray(ArrayList<EQTL> qtlBuffer) {
        ArrayList<MinimalEQTL> minimalQtlBuffer = new ArrayList<>();

        for (EQTL e : qtlBuffer) {
            minimalQtlBuffer.add(new MinimalEQTL(e));
        }

        return minimalQtlBuffer;
    }

    @Override
    public String toString() {
        char nullstr = '-';
        char tabStr = '\t';
        MinimalEQTL e = this;
        StringBuilder out = new StringBuilder();

        out.append(e.getPvalue());
        out.append(tabStr);

        out.append(e.getRsName());
        out.append(tabStr);
        out.append(e.getRsChr());
        out.append(tabStr);
        out.append(e.getRsChrPos());
        out.append(tabStr);
        out.append(e.getProbe());
        out.append(tabStr);
        out.append(e.getProbeChr());
        out.append(tabStr);
        out.append(e.getProbeChrPos());
        out.append(tabStr);
        out.append(nullstr);
        out.append(tabStr);
        out.append(nullstr);
        out.append(tabStr);
        out.append(nullstr);
        out.append(tabStr);
        out.append(nullstr);
        out.append(tabStr);
        out.append(nullstr);
        out.append(tabStr);
        out.append(nullstr);
        out.append(tabStr);
        out.append(nullstr);
        out.append(tabStr);
        out.append(nullstr);
        out.append(tabStr);
        out.append(nullstr);
        out.append(tabStr);
        out.append(nullstr);
        out.append(tabStr);
        out.append(nullstr);
        out.append(tabStr);
        out.append(nullstr);
        out.append(tabStr);
        out.append(nullstr);
        out.append(tabStr);
        out.append(nullstr);
        out.append(tabStr);
        out.append(nullstr);

        return out.toString();
    }

}
