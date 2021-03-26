package eqtlmappingpipeline.util;

import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.stats.ZScores;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.HashMap;

public class ReplaceBeta {

    // WARNING: this code assumes that N is equal for each eQTL (which may not be the case in a meta-analysis context)
    public void run(String in, String out, String mafqcfile) throws IOException {
        TextFile tf = new TextFile(in, TextFile.R);
        String[] header = tf.readLineElems(TextFile.tab);

        int metazcol = -1;
        int ncol = -1;
        int metabetacol = -1;
        int snpcol = -1;

        for (int i = 0; i < header.length; i++) {
            if (header[i].contains("OverallZScore")) {
                metazcol = i;
            } else if (header[i].contains("Meta-Beta (SE)")) {
                metabetacol = i;
            } else if (header[i].contains("SNPName")) {
                snpcol = i;
            } else if (header[i].contains("DatasetsNrSamples") || header[i].contains("SumNumberOfSamples")) {
                ncol = i;
            }
        }

        if (metazcol < 0 || metabetacol < 0 || ncol < 0 || snpcol < 0) {
            System.out.println("Error: could not find one of the columns. Found the following:");
            System.out.println("MetaZ: " + metazcol);
            System.out.println("SNP: " + snpcol);
            System.out.println("MetaBeta: " + metabetacol);
            System.out.println("Ncol: " + ncol);
            System.exit(-1);
        }

        // load maf per snp
        HashMap<String, Double> mafpersnp = new HashMap<String, Double>();
        TextFile tf1 = new TextFile(mafqcfile, TextFile.R);
        System.out.println("Loading MAF info from: " + mafqcfile);
        tf1.readLine();
        int lnctr2 = 0;
        String[] elems = tf1.readLineElems(TextFile.tab);
        while (elems != null) {
            String snp = elems[0];
            try {
                double maf = Double.parseDouble(elems[elems.length - 1]);
                if (!Double.isNaN(maf)) {
                    mafpersnp.put(snp, maf);
                }
            } catch (NumberFormatException e) {
                // System.out.println(snp + " has unparseable MAF: " + elems[elems.length - 1]);
            }
            elems = tf1.readLineElems(TextFile.tab);
            lnctr2++;
            if (lnctr2 % 1000000 == 0) {
                System.out.println(lnctr2 + " SNPs parsed sofar.");
            }

        }
        tf1.close();
        System.out.println("MAF info found for " + mafpersnp.size() + " variants, out of " + lnctr2 + " total.");

        TextFile tfo = new TextFile(out, TextFile.W);
        tfo.writeln(Strings.concat(header, Strings.tab));
        elems = tf.readLineElems(TextFile.tab);
        int lnctr = 0;
        System.out.println("Processing: " + in);
        while (elems != null) {

            String snp = elems[snpcol];

            Double maf = mafpersnp.get(snp);
            if (maf != null) {
                String nstr = elems[ncol];
                int n = 0;
                String[] nstrelems = nstr.split(";");
                for (String s : nstrelems) {
                    if (!s.equals("-")) {
                        try {
                            n += Integer.parseInt(s);
                        } catch (NumberFormatException e) {
                            try {
                                n += (int) (Double.parseDouble(s));
                            } catch (NumberFormatException e2) {

                            }
                        }
                    }
                }
                double z = Double.parseDouble(elems[metazcol]);

                double[] betaandse = ZScores.zToBeta(z, maf, n);
                elems[metabetacol] = betaandse[0] + " (" + betaandse[1] + ")";
            }
            tfo.writeln(Strings.concat(elems, Strings.tab));
            elems = tf.readLineElems(TextFile.tab);
            lnctr++;
            if (lnctr % 100000 == 0) {
                System.out.println(lnctr + " lines parsed.");
            }
        }

        tf.close();
        tfo.close();
    }

}
