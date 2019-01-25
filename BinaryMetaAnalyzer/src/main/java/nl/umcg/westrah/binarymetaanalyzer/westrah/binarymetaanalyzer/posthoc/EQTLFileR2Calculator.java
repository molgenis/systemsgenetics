package nl.umcg.westrah.binarymetaanalyzer.westrah.binarymetaanalyzer.posthoc;

import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.stats.ZScores;
import umcg.genetica.text.Strings;

import java.io.IOException;

public class EQTLFileR2Calculator {

    public void run(String ein, String eout) throws IOException {


        TextFile tf = new TextFile(ein, TextFile.R);
        TextFile tout = new TextFile(eout, TextFile.W);

        tout.writeln(tf.readLine() + "\tExplainedVariance(R)");
        String[] elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {

            Double z = Double.parseDouble(elems[10]);
            String[] ns = elems[13].split(";");
            int sumN = 0;
            for (String n : ns) {
                try {
                    Integer q = Integer.parseInt(n);
                    sumN += q;
                } catch (NumberFormatException e) {

                }
            }

            double r = ZScores.zToR(z, sumN);

            tout.writeln(Strings.concat(elems, Strings.tab) + "\t" + r);

            elems = tf.readLineElems(TextFile.tab);
        }

        tout.close();
        tf.close();

    }

}
