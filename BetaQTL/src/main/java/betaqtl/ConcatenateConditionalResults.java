package betaqtl;

import betaqtl.vcf.DetermineLDVCF;
import betaqtl.vcf.VCFTabix;
import betaqtl.vcf.VCFVariant;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.containers.Pair;
import umcg.genetica.enums.Chromosome;
import umcg.genetica.features.Feature;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;

public class ConcatenateConditionalResults {

    class SNPObj {
        public String content;
        String name;
        int chr;
        int pos;
        int rank;
    }

    public void run(String qtlfile1, String vcfStr, String samplefile, int nrTotalIterations, String output) throws IOException {

        HashMap<String, ArrayList<SNPObj>> geneSnpMap = new HashMap<>();
        HashMap<String, Integer> genechr = new HashMap<>();
        String outheader = null;
        for (int iteration = 1; iteration < nrTotalIterations + 1; iteration++) {
            String actualQTLFile = qtlfile1.replaceAll("ITER", "" + iteration);
            if (Gpio.exists(actualQTLFile)) {
                System.out.println("Found QTL file: " + actualQTLFile);
                TextFile tf = new TextFile(actualQTLFile, TextFile.R);
                int genecol = -1;
                int snpcol = -1;
                int snpchrcol = -1;
                int snpposcol = -1;
                String[] header = tf.readLineElems(TextFile.tab);
                if (outheader == null) {
                    outheader = Strings.concat(header, Strings.tab) + "\tQTLRank\n";
                }
                for (int i = 0; i < header.length; i++) {
                    if (header[i].equals("Gene")) {
                        genecol = i;
                    }
                    if (header[i].equals("SNP")) {
                        snpcol = i;
                    }
                    if (header[i].equals("SNPChr")) {
                        snpchrcol = i;
                    }
                    if (header[i].equals("SNPPos")) {
                        snpposcol = i;
                    }
                }
                String[] elems = tf.readLineElems(TextFile.tab);
                int ectr = 0;
                while (elems != null) {
                    String gene = elems[genecol];
                    String snp = elems[snpcol];
                    int chr = Integer.parseInt(elems[snpchrcol]);
                    int pos = Integer.parseInt(elems[snpposcol]);
                    SNPObj obj = new SNPObj();
                    obj.name = snp;
                    obj.pos = pos;
                    obj.chr = chr;
                    genechr.put(gene, chr);
                    obj.rank = iteration;
                    obj.content = Strings.concat(elems, Strings.tab);
                    ArrayList<SNPObj> geneObjs = geneSnpMap.get(gene);
                    if (geneObjs == null) {
                        geneObjs = new ArrayList<>();
                    }
                    geneObjs.add(obj);
                    geneSnpMap.put(gene, geneObjs);
                    ectr++;
                    elems = tf.readLineElems(TextFile.tab);
                }
                tf.close();
                System.out.println(actualQTLFile + " -- " + geneSnpMap.size() + " unique genes so far, " + ectr + " in this file");
            }

        }

        if (output.endsWith("/")) {
            output = output + "UnlinkedEQTLs.txt.gz";
        }
        TextFile out = new TextFile(output, TextFile.W);
        TextFile outld = new TextFile(output + "-LDLog.txt.gz", TextFile.W);
        out.writeln(outheader);
        outld.writeln("Gene\tSNP1\tSNP1Rank\tSNP2\tSNP2Rank\tRSq");
        // determine linkage per gene
        ProgressBar pb = new ProgressBar(geneSnpMap.keySet().size(), "Checking LD...");
        for (String gene : geneSnpMap.keySet()) {
            ArrayList<SNPObj> snpobjs = geneSnpMap.get(gene);
            if (snpobjs.size() == 1) {
                out.writeln(snpobjs.get(0).content + "\t" + 1);
            } else if (snpobjs.size() > 1) {
                // determine LD between variants
                // open tabix
                String tabixFile = vcfStr.replaceAll("CHR", "" + genechr.get(gene));
                VCFTabix tabix = new VCFTabix(tabixFile);
                boolean[] samplefilter = null;
                if (samplefile != null) {
                    samplefilter = tabix.getSampleFilter(samplefile);
                }
                double[][] ldmatrix = new double[snpobjs.size()][snpobjs.size()];
                for (int i = 0; i < snpobjs.size(); i++) {
                    VCFVariant variant1 = getVariant(tabix, samplefilter, snpobjs.get(i));
                    for (int j = i + 1; j < snpobjs.size(); j++) {
                        VCFVariant variant2 = getVariant(tabix, samplefilter, snpobjs.get(j));
                        if (variant1 != null && variant2 != null) {
                            DetermineLDVCF ld = new DetermineLDVCF();
                            Pair<Double, Double> ldinfo = ld.getLD(variant1, variant2);
                            double rsq = ldinfo.getRight();
                            if (Double.isNaN(rsq) || rsq > 0.8) {
                                ldmatrix[i][j] = 1;
                                ldmatrix[j][i] = 1;
                            }
                            outld.writeln(gene + "\t" + snpobjs.get(i).name + "\t" + snpobjs.get(i).rank + "\t" + snpobjs.get(j).name + "\t" + snpobjs.get(j).rank + "\t" + rsq);
                        } else {
                            System.out.println("Warning variant not found for gene: " + gene + "\t" + variant1 + "\t" + variant2);
                            ldmatrix[i][j] = 1;
                            ldmatrix[j][i] = 1;
                            outld.writeln(gene + "\t" + snpobjs.get(i).name + "\t" + snpobjs.get(i).rank + "\t" + snpobjs.get(j).name + "\t" + snpobjs.get(j).rank + "\t" + Double.NaN);
                        }
                    }
                }
                // always keep the primary variant
                HashSet<Integer> linkedVariants = new HashSet<>();
                for (int i = 1; i < snpobjs.size(); i++) {

                    ArrayList<Integer> allowedIndexes = new ArrayList<>();
                    for (int j = 0; j < i; j++) {
                        if (!linkedVariants.contains(j)) {
                            allowedIndexes.add(j);
                        }
                    }

                    // check LD with previous iterations
                    for (int q = 0; q < allowedIndexes.size(); q++) {
                        int j = allowedIndexes.get(q);
                        if (ldmatrix[i][j] == 1) {
                            linkedVariants.add(i);
                            outld.writeln(gene + "\tExcluding: " + snpobjs.get(i).name + "\t" + snpobjs.get(i).rank);
                        }
                    }
                }

                int newRank = 1;
                for (int i = 0; i < snpobjs.size(); i++) {
                    if (!linkedVariants.contains(i)) {
                        out.writeln(snpobjs.get(i).content + "\t" + newRank);
                        newRank++;
                    }
                }
                tabix.close();
            }
            pb.iterate();
        }
        pb.close();
        out.close();
        outld.close();

    }

    private VCFVariant getVariant(VCFTabix tabix, boolean[] samplefilter, SNPObj snpobj) throws IOException {
        Feature cisRegion = new Feature(Chromosome.parseChr("" + snpobj.chr), snpobj.pos - 1, snpobj.pos + 1);

        Iterator<VCFVariant> snpIterator = tabix.getVariants(cisRegion, samplefilter);
        while (snpIterator.hasNext()) {
            VCFVariant variant = snpIterator.next();
            if (variant != null) {
                String variantId = variant.getId();
                if (variantId.equals(snpobj.name)) {
                    return variant;
                }
            }
        }
        return null;
    }
}
