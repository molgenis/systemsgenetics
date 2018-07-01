package nl.umcg.westrah.binarymetaanalyzer.westrah.binarymetaanalyzer.posthoc.ld;

import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.util.BaseAnnot;
import umcg.genetica.text.Strings;
import umcg.genetica.util.Primitives;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.atomic.AtomicInteger;


public class BinaryFileLD {


    public static void main(String[] args) {
        String indir = "D:\\ld\\tmp\\";
        String geneloc = "D:\\ld\\tmp\\sortedGeneSNPCombos.txt.gz-genes.txt";
        int nrthreads = 1;
        boolean writezmat = false;
        int minnrsamples = 0;
        int minnrvals = 0;
        String outdir = "D:\\ld\\tmp\\out\\";
        BinaryFileLD l = new BinaryFileLD();
        try {
            l.run(indir, geneloc, nrthreads, outdir, writezmat, minnrsamples, minnrvals);
        } catch (IOException e) {
            e.printStackTrace();
        }


    }

    public void run(String indir, String geneloc, int nrthreads, String outdir, boolean writezmat, int minnrsamples, int minnrvalues) throws IOException {

        boolean useinnerproduct = false;
        boolean weightforsamplesize = true;

        boolean debug = false;

        File[] datasetlocs = getLocs(indir);
        int nrdatasets = datasetlocs.length;

        // read genes to test
        TextFile tf = new TextFile(geneloc, TextFile.R);
        ArrayList<String> genes = tf.readAsArrayList();
        HashSet<String> geneset = new HashSet<>();
        geneset.addAll(genes);
        tf.close();

        // initialize dataset
        SortedBinaryZScoreFile[] dataset = new SortedBinaryZScoreFile[nrdatasets];
        SortedBinaryZDataBlock[] currentblock = new SortedBinaryZDataBlock[nrdatasets];
        for (int i = 0; i < nrdatasets; i++) {
            dataset[i] = new SortedBinaryZScoreFile(datasetlocs[i], SortedBinaryZScoreFile.R);
            currentblock[i] = dataset[i].readNextBlock();

            // fast forward to next available gene in the list
            if (currentblock[i] != null) {
                boolean fastforward = false;
                if (!geneset.contains(currentblock[i].gene)) {
                    System.out.println("Fast forwarding dataset: " + i + ", current gene: " + currentblock[i].gene);
                    fastforward = true;
                }
                while (currentblock[i] != null && !geneset.contains(currentblock[i].gene)) {
                    currentblock[i] = dataset[i].readNextBlock();
                }
                if (fastforward && currentblock[i] != null) {
                    System.out.println("Fast forwarding dataset: " + i + ", result: " + currentblock[i].gene);
                }

            }

        }

        // iterate genes
        ExecutorService executor = Executors.newFixedThreadPool(nrthreads);
        AtomicInteger nrRunning = new AtomicInteger(0);
        AtomicInteger nrCompleted = new AtomicInteger(0);
        int nrsub = 0;
        for (int g = 0; g < genes.size(); g++) {
            String querygene = genes.get(g);
//            System.out.println("Gene: " + g + " / " + genes.size() + ": " + querygene);
            // get gene data for each dataset
            LinkedHashMap<String, Integer> snpmap = new LinkedHashMap<>();
            ArrayList<ArrayList<SortedBinaryZDataBlock>> geneblocks = new ArrayList<>();
            int ctr = 0;
            for (int d = 0; d < dataset.length; d++) {
                SortedBinaryZDataBlock b = currentblock[d];
                ArrayList<SortedBinaryZDataBlock> blocks = new ArrayList<>();
                if (debug && b != null) {
                    System.out.println("Dataset: " + d + " current gene: " + b.gene);
                }

                // TODO: what happens if a gene is not in the genes "genes" object?
                boolean fastforward = false;
                while (b != null && b.gene.equals(querygene)) {
                    // add
                    blocks.add(b);
                    if (!snpmap.containsKey(b.snp)) {
                        snpmap.put(b.snp, ctr);
                        ctr++;
                    }

                    b = dataset[d].readNextBlock();
                    currentblock[d] = b;
                    if (debug && b != null) {
                        System.out.println("Dataset: " + d + " next gene: " + b.gene);
                    }

                    if (b != null && !geneset.contains(b.gene)) {
                        fastforward = true;
                        System.out.println("Fast forwarding dataset: " + d + ", current gene: " + currentblock[d].gene);
                    }
                    // fast forward to next available gene in the list
                    while (b != null && !geneset.contains(b.gene)) {
                        b = dataset[d].readNextBlock();
                        currentblock[d] = b;
                    }

                    if (fastforward && currentblock[d] != null) {
                        System.out.println("Fast forwarding dataset: " + d + ", result: " + currentblock[d].gene);
                    }

                }
                if (debug) {
                    System.out.println("Done reading ds: " + d);
                }
                geneblocks.add(blocks);
            }


            // prevent clogging of memories.
            int nrrunningi = nrRunning.get();
            while (nrrunningi >= nrthreads) {
                try {
                    System.out.println("Submitted: " + nrsub + "\tCompleted: " + nrCompleted.get() + "\tRunning: " + nrrunningi);
                    Thread.sleep(1000);
                } catch (InterruptedException e) {
                    e.printStackTrace();
                }
                nrrunningi = nrRunning.get();
            }

            if (!snpmap.isEmpty()) {
                SortedZCorrelationTask t = new SortedZCorrelationTask(geneblocks,
                        querygene,
                        snpmap,
                        dataset,
                        outdir,
                        nrRunning,
                        nrCompleted,
                        writezmat,
                        useinnerproduct,
                        weightforsamplesize,
                        g,
                        minnrsamples,
                        minnrvalues
                );

                executor.submit(t);
                nrsub++;
            }
            if (nrsub > 20) {
                break;
            }
        }

        int nrrunningi = nrRunning.get();
        while (nrrunningi > 0) {
            try {
                System.out.println("Submitted: " + nrsub + "\tCompleted: " + nrCompleted.get() + "\tRunning: " + nrrunningi);
                Thread.sleep(10000); // wait 10 secs before checking again
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
            nrrunningi = nrRunning.get();
        }
        executor.shutdown();
        System.out.println("Submitted: " + nrsub + "\tCompleted: " + nrCompleted.get() + "\tRunning: " + nrrunningi);
        System.out.println("Done");
    }

    private File[] getLocs(String indir) {

        File dir = new File(indir);
        File[] list = dir.listFiles();
        ArrayList<File> locs = new ArrayList<>();
        for (File f : list) {
            if (!f.isDirectory() && f.getName().endsWith("-data.dat")) {
                System.out.println("Found: " + f.getName());
                locs.add(f);
            }
        }
        System.out.println(locs.size() + " total datasets.");
        return locs.toArray(new File[0]);
    }

    class SortedZCorrelationTask implements Runnable {

        private final AtomicInteger ctr;
        private final AtomicInteger cmctr;
        private final boolean writezmat;
        private final int id;
        ArrayList<ArrayList<SortedBinaryZDataBlock>> geneblocks;
        String querygene;
        LinkedHashMap<String, Integer> snpmap;
        SortedBinaryZScoreFile[] dataset;
        String outdir;
        int minNrValues = 0;
        int minNrSamples = 0;
        boolean useinnerproduct;
        boolean weightforsamplesize;

        public SortedZCorrelationTask(ArrayList<ArrayList<SortedBinaryZDataBlock>> geneblocks, String querygene,
                                      LinkedHashMap<String, Integer> snpmap, SortedBinaryZScoreFile[] dataset,
                                      String outdir, AtomicInteger ctr, AtomicInteger cmctr, boolean writezmat,
                                      boolean useinnerproduct, boolean weightforsamplesize, int id, int minNrSamples, int minNrValues) {
            this.geneblocks = geneblocks;
            this.querygene = querygene;
            this.snpmap = snpmap;
            this.dataset = dataset;
            this.outdir = outdir;
            this.ctr = ctr;
            this.cmctr = cmctr;
            this.writezmat = writezmat;
            this.id = id;
            this.useinnerproduct = useinnerproduct;
            this.weightforsamplesize = weightforsamplesize;
            this.minNrSamples = minNrSamples;
            this.minNrValues = minNrValues;
        }

        @Override
        public void run() {
            ctr.getAndIncrement();
            // determine which snps were tested in which datasets
            System.out.println("Job" + id + "\tGene " + querygene + "\tBuilding block index: " + snpmap.size() + " x " + dataset.length);
            SortedBinaryZDataBlock[][] blockindex = new SortedBinaryZDataBlock[snpmap.size()][dataset.length];
            for (int d = 0; d < dataset.length; d++) {
                ArrayList<SortedBinaryZDataBlock> blocks = geneblocks.get(d);
                for (SortedBinaryZDataBlock b : blocks) {
                    Integer id = snpmap.get(b.snp);
                    if (id != null) {
                        blockindex[id][d] = b;
                    }
                }
            }

            System.out.println("Job" + id + "\tGene " + querygene + "\tFlipping alleles");
            // flip the z-scores, using a reference dataset (TODO: flip to 1kg reference allele?)
            SortedBinaryZDataBlock[] refBlocks = new SortedBinaryZDataBlock[snpmap.size()];
            for (int s = 0; s < snpmap.size(); s++) {
                String refAlleles = null;
                String refAllelesAssessed = null;
                for (int d = 0; d < dataset.length; d++) {
                    SortedBinaryZDataBlock b = blockindex[s][d];
                    if (b != null) {
                        if (refAlleles == null) {
                            refAlleles = b.allele;
                            refAllelesAssessed = b.assessed;
                            refBlocks[s] = b;
                        } else {

                            Boolean flip = BaseAnnot.flipalleles(refAlleles, refAllelesAssessed, b.allele, b.assessed);
                            if (flip == null) {
                                // something went wrong while flipping the snp
                                blockindex[s][d] = null;
                            } else if (flip) {
                                for (int z = 0; z < b.z.length; z++) {
                                    b.z[z] *= -1;
                                }
                            }
                        }
                    }
                }
            }

            // now correlate
            System.out.println("Job" + id + "\tGene " + querygene + "\thas " + snpmap.size() + " snps");
            try {
                PearsonsCorrelation c = new PearsonsCorrelation();
                TextFile output = new TextFile(outdir + querygene + ".txt.gz", TextFile.W, 32 * 1024);
                TextFile alleleout = new TextFile(outdir + querygene + "-referenceAlleles.txt.gz", TextFile.W, 32 * 1024);
                alleleout.writeln("SNP\tAlleles\tAssessed");

                String header = "SNP1\tSNP2\tNrValues\tNrSamples";
                if (useinnerproduct) {
                    header += "\tInnerProduct\tR\tRSq";
                } else {
                    header += "\tPearsonR\tPearsonRSq\tWeightedR\tWeightedRSq";
                }


                output.writeln(header);
                for (int s1 = 0; s1 < snpmap.size(); s1++) {
                    SortedBinaryZDataBlock[] snp1 = blockindex[s1];
                    SortedBinaryZDataBlock refblock1 = refBlocks[s1];
//                    for (int s2 = (s1 + 1); s2 < snpmap.size(); s2++) {
                    for (int s2 = 0; s2 < snpmap.size(); s2++) {
                        SortedBinaryZDataBlock[] snp2 = blockindex[s2];


                        ReturnObj zs = removeNullsAndConvertToDouble(snp1, snp2);
                        if (zs.x.length >= minNrValues && zs.nx[0] >= minNrSamples) {
                            SortedBinaryZDataBlock refblock2 = refBlocks[s2];
                            String outln = refblock1.snp
                                    + "\t" + refblock2.snp;

                            double wcorr = weightedcorr(zs);
                            double wrsq = wcorr * wcorr;
                            double corr = c.correlation(zs.x, zs.y);
                            double rsq = corr * corr;
                            outln += "\t" + zs.x.length + "\tNA\t" + corr + "\t" + rsq + "\t" + wcorr + "\t" + wrsq;

                            double ip = innerproduct(zs, false);
                            double ipw = innerproduct(zs, true);


                            output.writeln(outln);
                        }


                    }
                    alleleout.writeln(refblock1.snp + "\t" + refblock1.allele + "\t" + refblock1.assessed);
                }
                alleleout.close();
                output.close();

                if (writezmat) {
                    System.out.println("Job" + id + "\tGene " + querygene + "\tInitializing z-mat: " + snpmap.size() + " x " + (dataset.length * 10));
                    float[][] zout = new float[dataset.length * 10][snpmap.size()];
                    for (int s1 = 0; s1 < snpmap.size(); s1++) {
                        SortedBinaryZDataBlock[] snp1 = blockindex[s1];
                        int j = 0;
                        for (int d = 0; d < snp1.length; d++) {
                            if (snp1[d] == null) {
                                for (int q = 0; q < 10; q++) {
                                    zout[j][s1] = Float.NaN;
                                    j++;
                                }

                            } else {
                                float[] zs = snp1[d].z;
                                for (int q = 0; q < 10; q++) {
                                    zout[j][s1] = zs[q];
                                    j++;
                                }

                            }
                        }
                    }

                    System.out.println("Job" + id + "\tGene " + querygene + "\tWriting mat");
                    TextFile outz = new TextFile(outdir + querygene + "-zmat.txt.gz", TextFile.W);
                    header = "-";

                    for (int s1 = 0; s1 < snpmap.size(); s1++) {
                        header += "\t" + refBlocks[s1].snp;
                    }
                    outz.writeln(header);
//					for (int d = 0; d < dataset.length; d++) {
//						for (int p = 1; p < 11; p++) {
//							header += "\t" + dataset[d].getName() + "-perm-" + p;
//						}
//					}
                    int j = 0;
                    for (int d = 0; d < dataset.length; d++) {
                        for (int p = 0; p < 10; p++) {
                            float[] zs = zout[j];
                            String ln = dataset[d].getName() + "-perm-" + (p + 1) + "\t" + Strings.concat(zs, Strings.tab);
                            outz.writeln(ln);
                            j++;
                        }
                    }
                    outz.close();
                    System.out.println("Job" + id + "\tGene " + querygene + "\tDone writing mat");
                }

            } catch (IOException e) {
                e.printStackTrace();
            } catch (Exception e) {
                e.printStackTrace();
            }

            ctr.getAndDecrement();
            cmctr.getAndIncrement();
        }

        private double innerproduct(ReturnObj vals, boolean weightforsamplesize) {
            double sum = 0;
            double wsum = 0;
            for (int i = 0; i < vals.x.length; i++) {
                sum += (vals.x[i] * vals.y[i]);
                wsum += (vals.x[i] * vals.y[i]) * Math.sqrt(vals.nx[i]) * Math.sqrt(vals.ny[i]);
            }

            return sum;
        }

    }


    public class ReturnObj {

        double[] x;
        double[] y;
        double[] nx;
        double[] ny;

        double[] wnx;
        double[] wny;


        public ReturnObj(ArrayList<Double> x, ArrayList<Double> y, ArrayList<Double> n1, ArrayList<Double> n2) {
            this.x = Primitives.toPrimitiveArr(x.toArray(new Double[0]));
            this.y = Primitives.toPrimitiveArr(y.toArray(new Double[0]));
            this.nx = Primitives.toPrimitiveArr(n1.toArray(new Double[0]));
            this.ny = Primitives.toPrimitiveArr(n2.toArray(new Double[0]));

            double sumnwx = 0;
            double sumnwy = 0;

            wnx = new double[nx.length];
            wny = new double[ny.length];

            for (int i = 0; i < nx.length; i++) {
                wnx[i] = Math.sqrt(nx[i]);
                sumnwx += wnx[i];
                wny[i] = Math.sqrt(ny[i]);
                sumnwy += wny[i];
            }

            for (int i = 0; i < nx.length; i++) {
                wnx[i] /= sumnwx;
                wny[i] /= sumnwy;
            }

        }
    }

    private ReturnObj removeNullsAndConvertToDouble(SortedBinaryZDataBlock[] snp1, SortedBinaryZDataBlock[] snp2) {
        ArrayList<Double> x = new ArrayList<>();
        ArrayList<Double> y = new ArrayList<>();
        ArrayList<Double> n1 = new ArrayList<>();
        ArrayList<Double> n2 = new ArrayList<>();

        for (int i = 0; i < snp1.length; i++) {
            SortedBinaryZDataBlock b1 = snp1[i];
            SortedBinaryZDataBlock b2 = snp2[i];
            if (b1 != null && b2 != null) {
                if (b1.n > 0 && b2.n > 0) {
                    for (int q = 0; q < b1.z.length; q++) {
                        x.add((double) b1.z[q]);
                        y.add((double) b2.z[q]);
                        n1.add((double) b1.n);
                        n2.add((double) b2.n);
                    }
                }
            }
        }

        return new ReturnObj(x, y, n1, n2);


    }


    public double weightedcorr(ReturnObj vals) {


        double[] w = new double[vals.x.length];

        double xavg = 0;
        double yavg = 0;

        for (int i = 0; i < vals.x.length; i++) { // for each dataset
            double wi = Math.sqrt((vals.nx[i] + vals.ny[i]) / 2);
            xavg += wi * vals.x[i];
            yavg += wi * vals.y[i];
            w[i] = wi;
        }

        xavg /= vals.x.length;
        yavg /= vals.x.length;

        double lxy = 0;
        double lxx = 0;
        double lyy = 0;

        for (int i = 0; i < vals.x.length; i++) {
            double wi = w[i];
            double xi = vals.x[i];
            double yi = vals.y[i];
            double xminavg = (xi - xavg);
            double yminavg = (yi - yavg);
            lxy += wi * xminavg * yminavg;
            lxx += wi * xminavg * xminavg;
            lyy += wi * yminavg * yminavg;
        }

        double r = lxy / Math.sqrt(lxx * lyy);
        return r;
    }

    public double weightedcorr2(ReturnObj vals) {


        double[] w = new double[vals.x.length];

        double xavg = 0;
        double yavg = 0;
        double sumw = 0;

        for (int i = 0; i < vals.x.length; i++) { // for each dataset
            double wi = Math.sqrt((vals.nx[i] + vals.ny[i]) / 2);
            xavg += wi * vals.x[i];
            yavg += wi * vals.y[i];
            sumw += wi;
            w[i] = wi;
        }

        xavg /= sumw;
        yavg /= sumw;

        double lxy = 0;
        double lxx = 0;
        double lyy = 0;

        for (int i = 0; i < vals.x.length; i++) {
            double wi = w[i];
            double xi = vals.x[i];
            double yi = vals.y[i];
            double xminavg = (xi - xavg);
            double yminavg = (yi - yavg);
            lxy += wi * xminavg * yminavg;
            lxx += wi * xminavg * xminavg;
            lyy += wi * yminavg * yminavg;
        }
        lxy /= sumw;
        lxx /= sumw;
        lyy /= sumw;

        double r = lxy / Math.sqrt(lxx * lyy);
        return r;
    }
}
