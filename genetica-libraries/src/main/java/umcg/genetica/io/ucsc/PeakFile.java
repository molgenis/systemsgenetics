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
public class PeakFile extends TextFile {

    private byte currentChr = -1;
    private int currentStart = 0;
    private int currentStep = 1;
    private int currentSpan = 1;
    private long nrLnsReadAfterHeader = 0;

    public static enum PEAKFORMAT {

        NARROWPEAK, BROADPEAK, GAPPEDPEAK;
    }
    private PEAKFORMAT peakformat = null;

    public PeakFile(String name, boolean mode) throws IOException {
        super(name, mode);
        // skip the header line, or actually check the format?
        // readLine();
    }

    public PeakFile(String name, boolean mode, PEAKFORMAT format) throws IOException {
        super(name, mode);
        // skip the header line, or actually check the format?
        // readLine();
        peakformat = format;
    }

    public UCSCDataObject parseLn() throws IOException {
        String line = readLine();
        if (line == null) {
//            System.out.println("Line is null");
            return null;
        }
        boolean isHeaderLine = false;
        if (line.toLowerCase().contains("track")
                || line.toLowerCase().contains("narrowpeak")
                || line.toLowerCase().contains("gappedpeak")
                || line.toLowerCase().contains("broadpeak")) {
            if (line.toLowerCase().contains("narrowpeak")) {
                peakformat = PEAKFORMAT.NARROWPEAK;
            } else if (line.toLowerCase().contains("gappedpeak")) {
                peakformat = PEAKFORMAT.GAPPEDPEAK;
            } else if (line.toLowerCase().contains("broadpeak")) {
                peakformat = PEAKFORMAT.BROADPEAK;
            } else {
                System.out.println("Error: unkown peak format!\n" + line);
            }
            System.out.println("File is: " + peakformat);
            isHeaderLine = true;
        }

        if (peakformat == null) {
            throw new IOException("Error: " + file.getAbsolutePath() + " does not adhere to any of the peak formats");
        }

        String[] elems = Strings.whitespace.split(line);
        if (isHeaderLine && elems.length > 1) {
            nrLnsReadAfterHeader = 0;
            return parseLn();
        } else if (elems.length > 0) {
            UCSCDataObject output = null;
            try {

                byte chr = -1;
                int chrStart = -1;
                int chrStop = -1;
                String name = null;
                int score = -1;
                String Strand = null;
                double signal = -1d;
                double pval = -1d;
                double qval = -1d;
                int peak = -1;

                chr = ChrAnnotation.parseChr(elems[0].replace("chr", ""));
                chrStart = Integer.parseInt(elems[1]);
                chrStop = Integer.parseInt(elems[2]);
                name = elems[3];
                score = Integer.parseInt(elems[4]);
                Strand = elems[5];

                if (peakformat == PEAKFORMAT.NARROWPEAK || peakformat == PEAKFORMAT.BROADPEAK) {
                    signal = Double.parseDouble(elems[6]);
                    pval = Double.parseDouble(elems[7]);
                    qval = Double.parseDouble(elems[8]);
                }

                if (peakformat == PEAKFORMAT.NARROWPEAK) {
                    peak = Integer.parseInt(elems[9]);
                }


                output = new UCSCDataObject(chr, chrStart, chrStop, signal, UCSCDataObject.SORTBY.CHRPOS);
                nrLnsReadAfterHeader++;


            } catch (NumberFormatException e) {
//                e.printStackTrace();
//                System.out.println("Error parsing line! " + Strings.concat(elems, Strings.tab));
//                return parseLn();
            }
            if (output == null) {
                return parseLn();
            } else {
                return output;
            }
        } else {
            return null;
        }
    }

    public long size() {
        return file.length();
    }
}
