

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.util;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;

import eqtlmappingpipeline.metaqtl3.containers.QTL;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.EQTL;
import umcg.genetica.io.trityper.QTLTextFile;
import umcg.genetica.text.Strings;


/**
 * @author harmjan
 */
public class QTLFileSorter {

    public void run(String efile, String outfile) throws IOException {
        System.out.println("Loading: " + efile);

        TextFile tf = new TextFile(efile, TextFile.R);
        ArrayList<QTLObj> eQtls = new ArrayList<>(1000000);
        int ctr = 0;
        int total = 0;
        tf.readLine();
        String ln = tf.readLine();
        while (ln != null) {
            if (eQtls.size() == 1000000) {
                if (eQtls != null) {
                    Collections.sort(eQtls, new QTLComparator());
                    System.out.println("Writing: " + eQtls.size() + " to batch " + ctr);
                    TextFile out = new QTLTextFile(outfile + "-tmp-" + ctr + ".txt.gz", QTLTextFile.W);
                    for (QTLObj obj : eQtls) {
                        out.writeln(obj.ln);
                    }
                    out.close();
                    ctr++;
                }
                eQtls = new ArrayList<>();
            }

            String[] elems = ln.split("\t");
            double z = Math.abs(Double.parseDouble(elems[10]));
            if (!Double.isNaN(z) && Double.isFinite(z)) {
                QTLObj obj = new QTLObj(z, ln);
                eQtls.add(obj);
                total++;
            } else {
//                System.out.println("Z == NaN: " + z + "\t" + ln);
            }

            ln = tf.readLine();
        }

        tf.close();


        if (eQtls.size() > 0) {
            Collections.sort(eQtls, new QTLComparator());
            System.out.println("Writing: " + eQtls.size() + " to batch " + ctr);
            TextFile out = new QTLTextFile(outfile + "-tmp-" + ctr + ".txt.gz", QTLTextFile.W);
            for (QTLObj obj : eQtls) {
                out.writeln(obj.ln);
            }
            out.close();
             ctr++;
        }
        eQtls = null;

        System.out.println(total + " eqtls over " + ctr + " batches");

        // merge files
        TextFile[] tfs = new TextFile[ctr];
        String[][] lastlnelems = new String[ctr][];
        double[] zarr = new double[ctr];
        for (int c = 0; c < ctr; c++) {
            System.out.println("Opening: " + outfile + "-tmp-" + c + ".txt.gz");
            tfs[c] = new TextFile(outfile + "-tmp-" + c + ".txt.gz", QTLTextFile.R);
            tfs[c].readLine();
            lastlnelems[c] = tfs[c].readLineElems(TextFile.tab);
            zarr[c] = Math.abs(Double.parseDouble(lastlnelems[c][10]));
        }

        boolean done = false;
        QTLTextFile out = new QTLTextFile(outfile, QTLTextFile.W);

        int written = 0;
        while (!done) {
            double maxZ = 0;
            // determine max ln
            Integer maxc = null;
            for (int c = 0; c < ctr; c++) {
                if (lastlnelems[c] != null) {
                    double z = zarr[c];
                    if (!Double.isNaN(z) && z >= maxZ) {
                        maxc = c;
                        maxZ = z;
                    }
                }
            }

            // write selected line
            if (maxc != null) {
                out.writeln(Strings.concat(lastlnelems[maxc], Strings.tab));
                written++;
                lastlnelems[maxc] = tfs[maxc].readLineElems(TextFile.tab);
                if (lastlnelems[maxc] != null) {
                    zarr[maxc] = Math.abs(Double.parseDouble(lastlnelems[maxc][10]));
                } else {
                    zarr[maxc] = Double.NaN;
                }
            } else {
                done = true;
            }
            if (written % 1000000 == 0) {
                System.out.println(written + " eqtls written.");
            }
        }
        out.close();

        for (int c = 0; c < ctr; c++) {
            tfs[c].close();
        }

        for (int c = 0; c < ctr; c++) {
            File f = new File(outfile + "-tmp-" + c + ".txt.gz");
            f.delete();
        }
        System.out.println("Done sorting");


    }

    class QTLObj {
        double z;
        String ln;

        public QTLObj(double z, String ln) {
            this.z = z;
            this.ln = ln;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;
            QTLObj qtlObj = (QTLObj) o;
            return Double.compare(qtlObj.z, z) == 0;
        }

        @Override
        public int hashCode() {
            return Objects.hash(z);
        }
    }

    class QTLComparator implements Comparator<QTLObj> {

        @Override
        public int compare(QTLObj o1, QTLObj o2) {

            if (o1.equals(o2)) {
                return 0;
            }
            double z1 = o1.z;
            double z2 = o2.z;
            if (z1 > z2) {
                return -1;
            } else if (z1 < z2) {
                return 1;
            } else {
                return 0;
            }
        }
    }
}
