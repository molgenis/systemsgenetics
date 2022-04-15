package mbqtl;

import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.QTLTextFile;
import umcg.genetica.io.trityper.util.ChrAnnotation;
import umcg.genetica.text.Strings;

import java.io.File;
import java.io.IOException;
import java.util.*;

public class QTLFileSorter {

    public static void main(String[] args) {
        QTLFileSorter s = new QTLFileSorter();
        s.test();

    }

    public void test() {
//        QTLObj obj1 = new QTLObj(15, 1E-12, (byte) 1, 1000, "", SORTBY.Z);
//        QTLObj obj2 = new QTLObj(10, 1E-10, (byte) 2, 10000, "", SORTBY.Z);
//        QTLComparator comp = new QTLComparator(SORTBY.P);
//        System.out.println(comp.compare(obj1, obj2));
    }

    public enum SORTBY {
        P, Z, GENEPOS, SNPPOS;
    }

    public void run(String efile, String outfile, SORTBY s) throws IOException {
        run(efile, outfile, 2500000, s);
    }

    public void run(String efile, String outfile, int n, SORTBY s) throws IOException {
        if (!Gpio.exists(efile)) {
            System.err.println("Could not find: " + efile);
        }

        System.out.println("Sorting: " + efile + " by " + s + " in batches of " + n);


        int pvalcol = -1;
        int zcol = -1;
        int geneposcol = -1;
        int genechrcol = -1;
        int snpposcol = -1;
        int snpchrcol = -1;

        TextFile tf = new TextFile(efile, TextFile.R);
        String header = tf.readLine();
        String[] headerElems = header.split("\t");
        for (int i = 0; i < headerElems.length; i++) {
            String h = headerElems[i];
            if (h.equals("MetaP")) {
                pvalcol = i;
                System.out.println("MetaP column: " + i);
            } else if (h.equals("MetaPZ")) {
                zcol = i;
                System.out.println("MetaZ column: " + i);
            } else if (h.equals("GeneChr")) {
                genechrcol = i;
                System.out.println("GeneChr column: " + i);
            } else if (h.equals("GenePos")) {
                geneposcol = i;
                System.out.println("GenePos column: " + i);
            } else if (h.equals("SNPChr")) {
                snpchrcol = i;
                System.out.println("SNPChr column: " + i);
            } else if (h.equals("SNPPos")) {
                snpposcol = i;
                System.out.println("SNPPos column: " + i);
            }
        }

        if (snpchrcol == -1 || snpposcol == -1 || pvalcol == -1 || zcol == -1 || genechrcol == -1 || geneposcol == -1) {
            System.out.println("Could not find all required columns");
            System.out.println("MetaP column: " + zcol);
            System.out.println("MetaZ column: " + pvalcol);
            System.out.println("GeneChr column: " + genechrcol);
            System.out.println("GenePos column: " + geneposcol);
            System.out.println("SNPChr column: " + snpchrcol);
            System.out.println("SNPPos column: " + snpposcol);
            tf.close();
            System.exit(-1);
        }

        int cap = n;
        ArrayList<QTLObj> eQtls = new ArrayList<>(cap);
        int batchctr = 0;
        int total = 0;

        String ln = tf.readLine();
        while (ln != null) {
            if (eQtls.size() == cap) {
                if (eQtls != null) {
                    System.out.println("Sorting: " + eQtls.size() + " to batch " + batchctr);
                    QTLObj[] qtlarr = eQtls.toArray(new QTLObj[0]);
                    Arrays.parallelSort(qtlarr, new QTLComparator(s));
                    System.out.println("Writing: " + eQtls.size() + " to batch " + batchctr + ": " + outfile + "-tmp-" + batchctr + ".txt.gz");
                    TextFile out = new TextFile(outfile + "-tmp-" + batchctr + ".txt.gz", TextFile.W);
                    out.writeln(header);
                    for (QTLObj obj : qtlarr) {
                        out.writeln(obj.ln);
                    }
                    out.close();
                    batchctr++;
                }
                eQtls = new ArrayList<>();
            }


            String[] elems = ln.split("\t");
            double z = Math.abs(Double.parseDouble(elems[zcol]));
            if (!Double.isNaN(z) && Double.isFinite(z)) {
                // Gene    GeneChr GenePos GeneStrand      GeneSymbol      SNP     SNPChr  SNPPos  SNPAlleles      SNPMinorAllele  SNPMinorAlleleFreqOverall       SNPEffectAllele MetaP   MetaPN  MetaPZ  MetaBeta        MetaSE  NrDatasets
                byte genechr = ChrAnnotation.parseChr(elems[genechrcol]);
                Integer genepos = Integer.parseInt(elems[geneposcol]);
                byte snpchr = ChrAnnotation.parseChr(elems[snpchrcol]);
                Integer snppos = Integer.parseInt(elems[snpposcol]);
                double p = Double.parseDouble(elems[pvalcol]);
                QTLObj obj = new QTLObj(z, p, genechr, genepos, snpchr, snppos, ln, s);
                eQtls.add(obj);
                total++;
            } else {
//                System.out.println("Z == NaN: " + z + "\t" + ln);
            }

            ln = tf.readLine();
        }
        tf.close();


        if (eQtls.size() > 0) {
            Collections.sort(eQtls, new QTLComparator(s));
            System.out.println("Writing: " + eQtls.size() + " to batch " + batchctr);
            TextFile out = new TextFile(outfile + "-tmp-" + batchctr + ".txt.gz", TextFile.W);
            out.writeln(header);
            for (QTLObj obj : eQtls) {
                out.writeln(obj.ln);
            }
            out.close();
            batchctr++;
        }
        eQtls = null;

        System.out.println(total + " eqtls over " + batchctr + " batches");

        // merge batch files
        TextFile[] tfs = new TextFile[batchctr];
        String[][] batchlnElems = new String[batchctr][];
//        double[] statArr = new double[batchctr];
        for (int c = 0; c < batchctr; c++) {
            System.out.println("Opening: " + outfile + "-tmp-" + c + ".txt.gz");
            tfs[c] = new TextFile(outfile + "-tmp-" + c + ".txt.gz", QTLTextFile.R);
            tfs[c].readLine();
            batchlnElems[c] = tfs[c].readLineElems(TextFile.tab);
//            if (s.equals(SORTBY.P)) {
//                statArr[c] = Math.abs(Double.parseDouble(lastlnelems[c][0]));
//            } else if (s.equals(SORTBY.Z)) {
//                statArr[c] = Math.abs(Double.parseDouble(lastlnelems[c][10]));
//            } else {
//                statArr[c] = Math.abs(Double.parseDouble(lastlnelems[c][3]));
//            }
        }

        boolean done = false;
        TextFile out = new TextFile(outfile, TextFile.W);
        out.writeln(header);
        QTLComparator comp = new QTLComparator(s);
        int written = 0;
        while (!done) {
            QTLObj maxObj = null;

            // determine best line over all batches
            Integer bestbatch = null;
            for (int c = 0; c < batchctr; c++) {
                if (batchlnElems[c] != null) {
                    String[] elems = batchlnElems[c];
                    double z = Math.abs(Double.parseDouble(elems[zcol]));
                    byte genechr = ChrAnnotation.parseChr(elems[genechrcol]);
                    Integer genepos = Integer.parseInt(elems[geneposcol]);
                    byte snpchr = ChrAnnotation.parseChr(elems[snpchrcol]);
                    Integer snppos = Integer.parseInt(elems[snpposcol]);
                    double p = Double.parseDouble(elems[pvalcol]);

                    QTLObj obj = new QTLObj(z, p, genechr, genepos, snpchr, snppos, null, s);

                    if (maxObj == null) {
                        maxObj = obj;
                        bestbatch = c;
                    } else {
                        int v = comp.compare(obj, maxObj);
                        if (v < 0) {
                            maxObj = obj;
                            bestbatch = c;
                        }

                    }
                }
            }

            // write selected line
            if (bestbatch != null) {
                out.writeln(Strings.concat(batchlnElems[bestbatch], Strings.tab));
                written++;
                batchlnElems[bestbatch] = tfs[bestbatch].readLineElems(TextFile.tab);
            } else {
                done = true;
            }


            if (written % 1000000 == 0) {
                System.out.println(written + " eqtls written.");
            }
        }
        out.close();

        for (int c = 0; c < batchctr; c++) {
            tfs[c].close();
        }

        for (int c = 0; c < batchctr; c++) {
            File f = new File(outfile + "-tmp-" + c + ".txt.gz");
            f.delete();
        }
        System.out.println("Done sorting");
    }


    class QTLObj {
        double z;
        String ln;
        SORTBY s;
        int snppos;
        byte snpchr;
        int genepos;
        byte genechr;
        double p;

        public QTLObj(double z, double p, byte genechr, int genepos, byte snpchr, int snppos, String ln, SORTBY s) {
            this.z = z;
            this.ln = ln;
            this.s = s;
            this.snpchr = snpchr;
            this.snppos = snppos;
            this.genechr = genechr;
            this.genepos = genepos;
            this.p = p;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;
            if (s.equals(SORTBY.P)) {
                QTLObj QTLZObj = (QTLObj) o;
                return Double.compare(QTLZObj.p, p) == 0;
            } else if (s.equals(SORTBY.Z)) {
                QTLObj QTLZObj = (QTLObj) o;
                int comp = Double.compare(Math.abs(QTLZObj.z), Math.abs(z));
                if (comp == 0) {
                    return Double.compare(QTLZObj.p, p) == 0;
                } else {
                    return false;
                }
            } else {
                QTLObj qtlPosObj = (QTLObj) o;
                if (genepos == qtlPosObj.genepos && genechr == qtlPosObj.genechr && snppos == qtlPosObj.snppos && snpchr == qtlPosObj.snpchr)
                    return true;
                return false;
            }
        }

        @Override
        public int hashCode() {
            if (s.equals(SORTBY.P)) {
                return Objects.hash(p);
            } else if (s.equals(SORTBY.Z)) {
                return Objects.hash(z);
            } else {
                return Objects.hash(snppos, snpchr, genepos, genechr);
            }
        }

    }

    class QTLComparator implements Comparator<QTLObj> {
        SORTBY s;

        QTLComparator(SORTBY s) {
            this.s = s;
        }

        @Override
        public int compare(QTLObj o1, QTLObj o2) {
            if (s.equals(SORTBY.P)) {
                if (o1.equals(o2)) {
                    // compare metabeta
                    return 0;
                }
                double p1 = o1.p;
                double p2 = o2.p;
                if (p1 < p2) {
                    return -1;
                } else if (p1 > p2) {
                    return 1;
                } else {
                    return 0;
                }
            } else if (s.equals(SORTBY.Z)) {
                if (o1.equals(o2)) {
                    // compare metabeta
                    return 0;
                }

                double z1 = o1.z;
                double z2 = o2.z;
                if (z1 > z2) {
                    return -1;
                } else if (z1 < z2) {
                    return 1;
                } else {
                    if (o1.p == o2.p) {
                        return 0;
                    } else if (o1.p > o2.p) {
                        return 1;
                    } else {
                        return -1;
                    }
                }
            } else {
                if (o1.equals(o2)) {
                    return 0;
                }

                if (s.equals(SORTBY.GENEPOS)) {
                    int compchr = Integer.compare(o1.genechr, o2.genechr);
                    if (compchr != 0) {
                        return compchr;
                    }
                    int comppos = Integer.compare(o1.genepos, o2.genepos);
                    if (comppos != 0) {
                        return comppos;
                    }
                    int compchr2 = Integer.compare(o1.snpchr, o2.snpchr);
                    if (compchr2 != 0) {
                        return compchr2;
                    }
                    int comppos2 = Integer.compare(o1.snppos, o2.snppos);
                    if (comppos2 != 0) {
                        return comppos2;
                    }
                    return 0;
                } else {
                    int compchr2 = Integer.compare(o1.snpchr, o2.snpchr);
                    if (compchr2 != 0) {
                        return compchr2;
                    }
                    int comppos2 = Integer.compare(o1.snppos, o2.snppos);
                    if (comppos2 != 0) {
                        return comppos2;
                    }
                    int compchr = Integer.compare(o1.genechr, o2.genechr);
                    if (compchr != 0) {
                        return compchr;
                    }
                    int comppos = Integer.compare(o1.genepos, o2.genepos);
                    if (comppos != 0) {
                        return comppos;
                    }
                    return 0;
                }
            }
        }
    }


}
