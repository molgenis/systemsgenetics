/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.io.ucsc;

import java.io.IOException;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.util.ChrAnnotation;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class WigFile extends TextFile {

    public long size() {
        return file.length();
    }

    public static enum STEPMODE {

        VARIABLESTEP, FIXEDSTEP;
    }
    private STEPMODE stepmode = STEPMODE.VARIABLESTEP;
    private byte currentChr = -1;
    private int currentStart = 0;
    private int currentStep = 1;
    private int currentSpan = 1;
    private long nrLnsReadAfterHeader = 0;

    public WigFile(String name, boolean mode) throws IOException {
        super(name, mode);
        // skip the header line, or actually check the format?
        // readLine();
    }

    public UCSCDataObject parseLn() throws IOException {
        String line = readLine();
        if (line == null) {
//            System.out.println("Line is null");
            return null;
        }
        boolean isHeaderLine = false;
        if (line.toLowerCase().contains("variablestep") || line.toLowerCase().contains("fixedstep")) {
            isHeaderLine = true;
        }

        String[] elems = Strings.whitespace.split(line);
        if (isHeaderLine && elems.length > 1) {
            for (int i = 0; i < elems.length; i++) {
                if (elems[i].toLowerCase().contains("variablestep")) {
                    stepmode = STEPMODE.VARIABLESTEP;
                } else if (elems[i].toLowerCase().contains("fixedstep")) {
                    stepmode = STEPMODE.FIXEDSTEP;
                } else if (elems[i].toLowerCase().contains("chrom")) {
                    currentChr = ChrAnnotation.parseChr(elems[i].substring(9));
                } else if (elems[i].toLowerCase().contains("span")) {
                    currentSpan = Integer.parseInt(elems[i].substring(5));
                } else if (elems[i].toLowerCase().contains("start")) {
                    currentStart = Integer.parseInt(elems[i].substring(6));
                } else if (elems[i].toLowerCase().contains("step")) {
                    currentStep = Integer.parseInt(elems[i].substring(5));
                }
            }
            nrLnsReadAfterHeader = 0;
            // if it is a header line, directly parse the next line and return that as a result.
            // this is however a bit tricky, since it uses recursion.
//            System.out.println(stepmode + "\tchr: " + currentChr + "\tspan: " + currentSpan + "\tstart: " + currentStart + "\tstep: " + currentStep);
            return parseLn();
        } else if (elems.length > 0) {
            int start = -1;
            double value = -1;
            try {
                if (stepmode == STEPMODE.FIXEDSTEP) {
                    value = Double.parseDouble(elems[0]);
                    if (nrLnsReadAfterHeader == 0) {
                        start = currentStart;
                    } else {
                        start = currentStart + currentStep;
                        currentStart += currentStep;
                    }
                } else {
                    start = Integer.parseInt(elems[0]);
                    value = Double.parseDouble(elems[1]);
                }
                UCSCDataObject output = new UCSCDataObject(currentChr, start, start + currentSpan, value, UCSCDataObject.SORTBY.CHRPOS);
                nrLnsReadAfterHeader++;
                return output;
            } catch (NumberFormatException e) {
                System.out.println("Error parsing line! " + stepmode + "\t" + Strings.concat(elems, Strings.tab));
                return parseLn();
            }


        } else {
            return null;
        }
    }
}
