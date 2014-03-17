/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.io.trityper;

import java.io.IOException;
import java.util.HashMap;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.util.ChrAnnotation;

/**
 *
 * @author harm-jan
 */
public final class ProbeAnnotation {

    private String[] probes;
    private int[] chrStart;
    private int[] chrEnd;
    private short[] chr;
    private String[] annotation;
    private HashMap<String, Integer> probeIndex;

    public ProbeAnnotation() {
    }

    public ProbeAnnotation(String loc) throws IOException {
        load(loc);
    }

    public void load(String loc) throws IOException {
        TextFile in = new TextFile(loc, TextFile.R);
        int probeNum = in.countLines() - 1;

        probes = new String[probeNum];
        annotation = new String[probeNum];
        chr = new short[probeNum];
        chrStart = new int[probeNum];
        chrEnd = new int[probeNum];


        in.readLine();
        probeIndex = new HashMap<String, Integer>();
        String[] elems = in.readLineElemsReturnReference(TextFile.tab);
        int i = 0;
        while (elems != null) {
            if (elems.length >= 6) {
                probes[i] = elems[1];
                annotation[i] = elems[2];
                chr[i] = (short) ChrAnnotation.parseChr(elems[3]);
                if (chr[i] != -1) {
                    chrStart[i] = Integer.parseInt(elems[4]);
                    chrEnd[i] = Integer.parseInt(elems[5]);
                } else {
                    chrStart[i] = -1;
                    chrEnd[i] = -1;
                }
                probeIndex.put(probes[i], i);
                i++;
            }
            elems = in.readLineElemsReturnReference(TextFile.tab);
        }
        in.close();
    }

    /**
     * @return the probes
     */
    public String[] getProbes() {
        return probes;
    }

    /**
     * @return the chrStart
     */
    public int[] getChrStart() {
        return chrStart;
    }

    /**
     * @return the chrEnd
     */
    public int[] getChrEnd() {
        return chrEnd;
    }

    /**
     * @return the chr
     */
    public short[] getChr() {
        return chr;
    }

    /**
     * @return the probeAnnotation
     */
    public String[] getProbeAnnotation() {
        return annotation;
    }

    /**
     * @return the probeToProbeId
     */
    public HashMap<String, Integer> getProbeToProbeId() {
        return probeIndex;
    }

    /**
     * @param probes the probes to set
     */
    public void setProbes(String[] probes) {
        this.probes = probes;
    }

    /**
     * @param chrStart the chrStart to set
     */
    public void setChrStart(int[] chrStart) {
        this.chrStart = chrStart;
    }

    /**
     * @param chrEnd the chrEnd to set
     */
    public void setChrEnd(int[] chrEnd) {
        this.chrEnd = chrEnd;
    }

    /**
     * @param chr the chr to set
     */
    public void setChr(short[] chr) {
        this.chr = chr;
    }

    /**
     * @param probeAnnotation the probeAnnotation to set
     */
    public void setProbeAnnotation(String[] probeAnnotation) {
        this.annotation = probeAnnotation;
    }

    /**
     * @param probeToProbeId the probeToProbeId to set
     */
    public void setProbeToProbeId(HashMap<String, Integer> probeToProbeId) {
        this.probeIndex = probeToProbeId;
    }
}
