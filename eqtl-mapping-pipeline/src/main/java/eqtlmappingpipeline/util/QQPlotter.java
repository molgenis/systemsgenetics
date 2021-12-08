package eqtlmappingpipeline.util;

import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.tdouble.impl.DenseLargeDoubleMatrix2D;
import eqtlmappingpipeline.metaqtl3.FDR;
import eqtlmappingpipeline.metaqtl3.graphics.QQPlot;
import gnu.trove.set.hash.THashSet;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.QTLTextFile;
import umcg.genetica.text.Strings;
import umcg.genetica.util.Primitives;

import java.io.IOException;
import java.util.ArrayList;

public class QQPlotter {


    public void run(String permutationDir, String outputfilename, double fdrcutoff, int nrPermutationsFDR, int maxNrMostSignificantEQTLs) throws IOException {

        String[] lsof = Gpio.getListOfFiles(permutationDir);

        String eqtlfdrfile = null;
        String prefix = "eQTLsFDR" + fdrcutoff;
        for (String s : lsof) {
            if (s.startsWith(prefix)) {
                eqtlfdrfile = s;
            }
        }
        if (eqtlfdrfile == null) {
            System.out.println("Could not find file starting with " + prefix + " in " + permutationDir);
            System.exit(-1);
        }
        eqtlfdrfile = permutationDir + "/" + eqtlfdrfile;
        FDR.FDRMethod m = FDR.FDRMethod.FULL;
        if (eqtlfdrfile.contains("ProbeLevel")) {
            m = FDR.FDRMethod.PROBELEVEL;
            eqtlfdrfile = permutationDir + "/eQTLsFDR-ProbeLevel.txt.gz";
        } else if (eqtlfdrfile.contains("GeneLevel")) {
            m = FDR.FDRMethod.GENELEVEL;
            eqtlfdrfile = permutationDir + "/eQTLsFDR-GeneLevel.txt.gz";
        } else {
            eqtlfdrfile = permutationDir + "/eQTLsFDR.txt.gz";
        }

        System.out.println("Using " + eqtlfdrfile + " as unpermuted eQTL input.");

        FDR.FileFormat f = FDR.FileFormat.REDUCED;

        DoubleMatrix2D permutedPValues;

        if ((nrPermutationsFDR * (long) maxNrMostSignificantEQTLs) < (Integer.MAX_VALUE - 2)) {
            permutedPValues = new DenseDoubleMatrix2D(nrPermutationsFDR, maxNrMostSignificantEQTLs);
        } else {
            permutedPValues = new DenseLargeDoubleMatrix2D(nrPermutationsFDR, maxNrMostSignificantEQTLs);
        }

        int nrEQTLs = -1;
        permutedPValues.assign(1);

        for (int permutationRound = 0; permutationRound < nrPermutationsFDR; permutationRound++) {
            String fileString = permutationDir + "/PermutedEQTLsPermutationRound" + (permutationRound + 1) + ".txt.gz";
            // read the permuted eqtl output
            System.out.println("Reading " + fileString);
            TextFile gz = new TextFile(fileString, TextFile.R);

            String[] header = gz.readLineElems(TextFile.tab);
            int snpcol = -1;
            int pvalcol = -1;
            int probecol = -1;
            int genecol = -1;

            //PValue  SNP     Probe   Gene
            if (header.length < 10) {
                f = FDR.FileFormat.REDUCED;
            } else {
                f = FDR.FileFormat.LARGE;
            }
            for (int col = 0; col < header.length; col++) {
                if (header[col].equals("PValue")) {
                    pvalcol = col;
                }
                if (header[col].equals("SNP")) {
                    snpcol = col;
                }
                if (header[col].equals("Probe")) {
                    probecol = col;
                }
                if (header[col].equals("Gene")) {
                    genecol = col;
                }
            }
            if (f == FDR.FileFormat.REDUCED) {
                //PValue  SNP     Probe   Gene
                if (snpcol == -1 || pvalcol == -1 || probecol == -1 && genecol == -1) {
                    System.out.println("Column not found in permutation file: " + fileString);
                    System.out.println("PValue: " + pvalcol);
                    System.out.println("SNP: " + snpcol);
                    System.out.println("Probe: " + probecol);
                    System.out.println("Gene: " + genecol);
                }
            }
//			String[] data = gz.readLineElemsReturnReference(TextFile.tab);
            int itr = 0;
            String permln = gz.readLine();
            THashSet<String> visitedEffects = new THashSet<String>();
            while (permln != null) {

                if (permln.length() != 0) {
                    if (itr > maxNrMostSignificantEQTLs - 1) {
                        break;
                    } else {
                        int filteronColumn;
                        String fdrId;
                        if (f == FDR.FileFormat.REDUCED) {
                            if (m == FDR.FDRMethod.PROBELEVEL) {
                                // fdrId = data[probecol];
                                fdrId = Strings.subsplit(permln, Strings.tab, probecol, probecol + 1)[0];
                                filteronColumn = probecol;
                            } else if (m == FDR.FDRMethod.SNPLEVEL) {
//								fdrId = data[snpcol];
                                fdrId = Strings.subsplit(permln, Strings.tab, snpcol, snpcol + 1)[0];
                                filteronColumn = snpcol;
                            } else if (m == FDR.FDRMethod.GENELEVEL && header.length > 3) {
//								fdrId = data[genecol];
                                fdrId = Strings.subsplit(permln, Strings.tab, genecol, genecol + 1)[0];
                                filteronColumn = genecol;
                            } else {
                                fdrId = Strings.subsplit(permln, Strings.tab, probecol, probecol + 1)[0];
                                filteronColumn = probecol;
                            }
                        } else {
                            if (m == FDR.FDRMethod.GENELEVEL) {
//								fdrId = data[QTLTextFile.HUGO];
                                fdrId = Strings.subsplit(permln, Strings.tab, QTLTextFile.HUGO, QTLTextFile.HUGO + 1)[0];
                                filteronColumn = QTLTextFile.HUGO;
                            } else if (m == FDR.FDRMethod.SNPLEVEL) {
//								fdrId = data[QTLTextFile.SNP];
                                fdrId = Strings.subsplit(permln, Strings.tab, QTLTextFile.SNP, QTLTextFile.SNP + 1)[0];
                                filteronColumn = QTLTextFile.HUGO;
                            } else if (m == FDR.FDRMethod.PROBELEVEL) {
//								fdrId = data[4];
                                fdrId = Strings.subsplit(permln, Strings.tab, 4, 5)[0];
                                filteronColumn = 4;
                            } else {
//								fdrId = data[1] + "\t" + data[4];
                                fdrId = Strings.subsplit(permln, Strings.tab, QTLTextFile.SNP, QTLTextFile.SNP + 1)[0] + "-" + Strings.subsplit(permln, Strings.tab, 4, 5)[0];
                                filteronColumn = 4;
                            }
                        }

                        // take top effect per gene / probe
                        if (Strings.countSeparatorOccurrences(permln, Strings.tab) > filteronColumn) {
                            if (!fdrId.equals("-") && !visitedEffects.contains(fdrId)) {
                                permutedPValues.setQuick(permutationRound, itr, Double.parseDouble(Strings.subsplit(permln, Strings.tab, 0, 1)[0]));
//                                permutedPValues[permutationRound][itr] = Double.parseDouble(data[0]);
                                visitedEffects.add(fdrId);
                                if (itr > 0 && permutedPValues.getQuick(permutationRound, (itr - 1)) > permutedPValues.getQuick(permutationRound, itr)) {
                                    System.err.println("Sorted P-Value list is not perfectly sorted!!!!");
                                    System.exit(-1);
                                }
                                itr++;
                            }
                        } else {
//							System.out.println(Strings.concat(data, Strings.tab));
                        }
                        permln = gz.readLine();
                    }
                }
            }
            gz.close();

            if (nrEQTLs == -1) {
                nrEQTLs = itr;
            }
        }


        ArrayList<Boolean> significantPvalueArl = new ArrayList<>();
        ArrayList<Double> pValueRealDataArl = new ArrayList<>();
        TextFile tf = new TextFile(eqtlfdrfile, TextFile.R);
        tf.readLine();
        String[] elems = tf.readLineElems(TextFile.tab);
        int nrSignificantEQTLs = 0;
        int nrTotalEQTLs = 0;
        while (elems != null) {
            if (elems.length > 3) {
                double pval = Double.parseDouble(elems[0]);
                double fdr = Double.parseDouble(elems[elems.length - 1]);
                pValueRealDataArl.add(pval);
                if (fdr < fdrcutoff) {
                    significantPvalueArl.add(true);
                    nrSignificantEQTLs++;
                } else {
                    significantPvalueArl.add(false);
                }
                nrTotalEQTLs++;
            }
            elems = tf.readLineElems(TextFile.tab);
        }

        double[] pValueRealData = Primitives.toPrimitiveArr(pValueRealDataArl);
        boolean[] significant = new boolean[pValueRealDataArl.size()];
        int pos = 0;
        for (Boolean i : significantPvalueArl) {
            significant[pos] = i;
            pos++;
        }

        System.out.println(eqtlfdrfile + " has " + pValueRealData.length + " eqtls, of which " + nrSignificantEQTLs + " are significant");

        QQPlot qq = new QQPlot();
        qq.draw(outputfilename, fdrcutoff, nrPermutationsFDR, nrTotalEQTLs, permutedPValues.toArray(), pValueRealData, significant, nrSignificantEQTLs);

    }

}
