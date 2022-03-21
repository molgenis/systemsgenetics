package mbqtl;

import mbqtl.vcf.DetermineLDVCF;
import mbqtl.vcf.VCFTabix;
import mbqtl.vcf.VCFVariant;
import umcg.genetica.containers.Pair;
import umcg.genetica.enums.Chromosome;
import umcg.genetica.features.Feature;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.HashMap;
import java.util.Iterator;

public class LDCalculator {

    class SNPObj {
        String name;
        int chr;
        int pos;
    }

    public void run(String qtlfile1, String qtlfile2, String vcfStr, String samplefile, String output) throws IOException {

        HashMap<String, SNPObj> geneSnpMap = new HashMap<>();

        TextFile tf = new TextFile(qtlfile1, TextFile.R);
        int genecol = -1;
        int snpcol = -1;
        int snpchrcol = -1;
        int snpposcol = -1;
        String[] header = tf.readLineElems(TextFile.tab);
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
        while (elems != null) {
            String gene = elems[genecol];
            String snp = elems[snpcol];
            int chr = Integer.parseInt(elems[snpchrcol]);
            int pos = Integer.parseInt(elems[snpposcol]);
            SNPObj obj = new SNPObj();
            obj.name = snp;
            obj.pos = pos;
            obj.chr = chr;
            geneSnpMap.put(gene, obj);
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();
        System.out.println(qtlfile1 + " has " + geneSnpMap.size() + " unique genes.");

        TextFile out = new TextFile(output, TextFile.W);
        out.writeln("Gene\tSNP1\tSNP2\tRsq\tDprime");
        TextFile tf2 = new TextFile(qtlfile2, TextFile.R);
        header = tf2.readLineElems(TextFile.tab);
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
        elems = tf2.readLineElems(TextFile.tab);
        int ctr = 0;
        int ctr2 = 0;
        while (elems != null) {
            String gene = elems[genecol];
            String snp = elems[snpcol];
            int chr = Integer.parseInt(elems[snpchrcol]);
            int pos = Integer.parseInt(elems[snpposcol]);
            SNPObj obj = new SNPObj();
            obj.name = snp;
            obj.pos = pos;
            obj.chr = chr;
            SNPObj other = geneSnpMap.get(gene);
            if (other != null) {
                if (chr == other.chr) {
                    String tabixFile = vcfStr.replaceAll("CHR", "" + chr);
                    VCFTabix tabix = new VCFTabix(tabixFile);
                    boolean[] samplefilter = null;
                    if (samplefile != null) {
                        samplefilter = tabix.getSampleFilter(samplefile);
                    }
                    VCFVariant variant1 = getVariant(tabix, samplefilter, other);
                    VCFVariant variant2 = getVariant(tabix, samplefilter, obj);
                    if (variant1 != null && variant2 != null) {
                        DetermineLDVCF ld = new DetermineLDVCF();
                        Pair<Double, Double> ldinfo = ld.getLD(variant1, variant2);
                        double dpr = ldinfo.getLeft();
                        double rsq = ldinfo.getRight();
                        out.writeln(gene + "\t" + other.name + "\t" + snp + "\t" + dpr + "\t" + rsq);
                        ctr2++;
                    } else {
                        System.out.println("Warning variant not found for gene: " + gene + "\t" + variant1 + "\t" + variant2);
                    }
                    tabix.close();
                }
            }
            ctr++;
            System.out.print(ctr + " lines processed. " + ctr2 + " LD calculations performed.\r");
            elems = tf2.readLineElems(TextFile.tab);
        }
        System.out.print(ctr + " genes processed. " + ctr2 + " genes overlapped.\n");
        tf2.close();
        out.close();
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
