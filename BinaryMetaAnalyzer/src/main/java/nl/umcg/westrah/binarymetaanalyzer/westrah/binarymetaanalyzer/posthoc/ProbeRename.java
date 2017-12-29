/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.westrah.binarymetaanalyzer.westrah.binarymetaanalyzer.posthoc;

import java.io.IOException;
import java.util.Map;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

/**
 *
 * @author Harm-Jan
 */
public class ProbeRename {

    public static void main(String[] args) {
        try {
//            String inf = args[0];
//            String annot = args[1];
//            String innnot = args[2];
//            String outannot = args[3];

            String inf = "C:\\Work\\emptest\\h8v2binarymeta\\eQTLs.txt.gz";
            String annot = "C:\\Work\\GGD\\2013-07-18-ProbeAnnotationFile-H8v2.txt";
            String innnot = "Probe";
            String outannot = "HumanRef-8v2-HT12v3";

            Map<String, String> map;
            TextFile tf = new TextFile(annot, TextFile.R);
            String[] header = tf.readLineElems(TextFile.tab);
            int col1 = -1;
            int col2 = -1;
            for (int col = 0; col < header.length; col++) {
                if (header[col].equals(innnot)) {
                    col1 = col;
                }
                if (header[col].equals(outannot)) {
                    col2 = col;
                }
            }

            if (col1 == -1) {
                System.err.println(innnot + " not found in " + annot);
            }
            if (col2 == -1) {
                System.err.println(outannot + " not found in " + annot);
            }
            System.out.println(col1 + "\t" + col2);
            map = tf.readAsHashMap(col1, col2);
            tf.close();

            TextFile tfin = new TextFile(inf, TextFile.R);
            TextFile out = new TextFile(inf + "-ProbesRewritten-" + outannot + ".txt", TextFile.W);
            out.writeln(tfin.readLine());
            String[] elems = tfin.readLineElems(TextFile.tab);
            while (elems != null) {
                elems[4] = map.get(elems[4]);
                if (elems[4] == null) {
                    // don't write?
                } else {
                    out.writeln(Strings.concat(elems, Strings.tab));
                }
                elems = tfin.readLineElems(TextFile.tab);
            }
            tfin.close();
            out.close();
        } catch (IOException e) {
            e.printStackTrace();
        }

    }
}
