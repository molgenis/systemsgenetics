/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.io.trityper.probeannotation;

import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.util.ChrAnnotation;

import java.io.IOException;
import java.util.HashMap;

/**
 *
 * @author harmjan
 */
public class ProbeTranslation {

    private String[] probeSymbol;
    private byte[] probeChr;
    private int[] probeChrPos;
    private HashMap<Integer, String> actualMappingPosition;
    private String[] probeName;
    private HashMap<String, Integer> oldToNewProbeAddress;
    private HashMap<String, Integer> toProbeAnnotationId;

    public void load(String probeAnnotation) throws IOException {
        TextFile tf = new TextFile(probeAnnotation, TextFile.R);
        int numSeq = tf.countLines();
        String[] elems = tf.readLineElemsReturnReference(TextFile.tab);
        int numAnnotations = elems.length - 5;

        int num = 0;
        System.out.println("Available annotations: ");
        String[] annotationname = new String[numAnnotations];
        for (int i = 5; i < elems.length; i++) {
            String name = elems[i];
            annotationname[i - 5] = name;
            System.out.println(name);
        }

        probeName = new String[numSeq];
        probeChr = new byte[numSeq];
        probeChrPos = new int[numSeq];
        probeSymbol = new String[numSeq];
        oldToNewProbeAddress = new HashMap<String, Integer>();
        toProbeAnnotationId = new HashMap<String, Integer>();
        elems = tf.readLineElemsReturnObjects(TextFile.tab);

        actualMappingPosition = new HashMap<Integer, String>();

        // array address (platform) --> nieuw gezamenlijk nummer + annotatie
        int probeNum = 0;
        while (elems != null) {
            Integer newProbeNum = null;
            try {
                newProbeNum = Integer.parseInt(elems[0]);
            } catch (NumberFormatException e) {
            }

            if (newProbeNum != null) {
                String seq = elems[1];
                String chr = elems[2];
                String chrpos = elems[3];
                String symbol = elems[4];
                num = 0;

                probeName[probeNum] = elems[0].intern();

                byte bchr = -1;
                try {

                    bchr = ChrAnnotation.parseChr(chr);
                } catch (Exception e) {
                    bchr = -1;
                    System.out.println("Cannot parse chr: " + chr);
                }

                int bchrpos = -1;
                try {
                    bchrpos = Integer.parseInt(chrpos);
                } catch (Exception e) {
                    String[] list = chrpos.split("-");
                    int max = -1;
                    int min = Integer.MAX_VALUE;
                    try {
                        for (String s : list) {
                            String[] list2 = s.split(":");
                            for (String s2 : list2) {
                                int p = Integer.parseInt(s2);
                                if (p > max) {
                                    max = p;
                                }
                                if (p < min) {
                                    min = p;
                                }
                            }
                        }

                        bchrpos = (int) Math.floor((double) (max + min) / 2);

                    } catch (Exception e2) {
                        System.out.println("Could not calculate midpos");
                        bchrpos = -1;
                    }
                }

                actualMappingPosition.put(probeNum, chrpos);
                probeChr[probeNum] = bchr;
                probeChrPos[probeNum] = bchrpos;
                probeSymbol[probeNum] = symbol.intern();

                for (int i = 5; i < elems.length; i++) {

                    String arrayaddress = elems[i];
                    if (!arrayaddress.equals("-")) {
                        try {
                            String[] addresselems = arrayaddress.split(",");
                            for (int q = 0; q < addresselems.length; q++) {
                                String address = addresselems[q].intern();

                                oldToNewProbeAddress.put(annotationname[i - 5] + address.intern(), probeNum);

                            }

                        } catch (Exception e) {
                        }
                    } else {
                        // probe not present
                    }

                    num++;
                }
            }


            elems = tf.readLineElemsReturnObjects(TextFile.tab);
            probeNum++;
        }
        tf.close();
    }

    public HashMap<String, String> getProbeTranslation(String probeTranslationFile, String source, String dest) throws IOException {
        System.out.println("Reading probe annotation table from: " + probeTranslationFile);
        System.out.println("Selecting " + source + " as translation for " + dest);

        HashMap<String, String> output = new HashMap<String, String>();

        TextFile tf = new TextFile(probeTranslationFile, TextFile.R);
        String[] header = tf.readLineElems(TextFile.tab);
        int sourcecol = -1;
        int destcol = -1;
        for (int i = 0; i < header.length; i++) {
            if (header[i].equals(source)) {
                sourcecol = i;
            }
            if (header[i].equals(dest)) {
                destcol = i;
            }
        }


        if (sourcecol < 0) {
            System.err.println("Column: " + source + " not found");
            return null;
        }
        if (destcol < 0) {
            System.err.println("Column: " + dest + " not found");
            return null;
        }

        System.out.println("Source: " + sourcecol + "\tDest: " + destcol);
        if (destcol >= 0 && sourcecol >= 0) {
            String[] elems = tf.readLineElems(TextFile.tab);
            while (elems != null) {
                String sourceStr = elems[sourcecol];
                String destStr = elems[destcol];
//                if (destStr.equals("-")) {
//                    System.out.println("Warning: " + sourceStr + " equals " + destStr + "\t" + destcol + "\t" + sourcecol);
//                }
                output.put(sourceStr, destStr);
                elems = tf.readLineElems(TextFile.tab);
            }

        }
        tf.close();
        return output;
    }

    public byte getProbeChr(int p) {
        return probeChr[p];
    }

    public String getActualMappingPosition(Integer probeId) {
        return actualMappingPosition.get(probeId);
    }

    public int getProbeChrPos(int p) {
        return probeChrPos[p];
    }

    public String getProbeSymbol(int p) {
        return probeSymbol[p];
    }

    public int getNumProbes() {
        return probeName.length;
    }

    public Integer getProbeId(String string) {
        return oldToNewProbeAddress.get(string);
    }

    public HashMap<String, Integer> getProbeTranslationTable() {
        return oldToNewProbeAddress;
    }

    public String[] getProbes() {
        return probeName;
    }
}
