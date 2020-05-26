/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.gwas;

import gnu.trove.map.hash.THashMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Vector;
import java.util.logging.Level;
import java.util.logging.Logger;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.SNPLoader;
import umcg.genetica.io.trityper.TriTyperGenotypeData;
import umcg.genetica.io.trityper.util.DetermineLD;
import umcg.genetica.math.matrix.DoubleMatrixDataset;

/**
 *
 * All shall be independent.
 *
 * @author juha
 */
public class Independifier {

    private static final Logger LOGGER = Logger.getLogger(Independifier.class.getName());
    TriTyperGenotypeData genotypeData;
    final SNPLoader snpLoader;
    Map<String, DoubleMatrixDataset<String, String>> permDatasets;

    public Independifier(String datadir) throws IOException {
        genotypeData = new TriTyperGenotypeData();
        genotypeData.load(datadir);
        snpLoader = genotypeData.createSNPLoader();
    }

    public Independifier(TriTyperGenotypeData genotypeData, SNPLoader snpLoader) throws IOException {
        this.genotypeData = genotypeData;
        this.snpLoader = snpLoader;
    }

    public void selectSNPsWithSimilarMAFsAsRealInputSNPs(String[] snps) throws IOException {

        TObjectIntHashMap<String> snpToSNPId = genotypeData.getSnpToSNPId();

        Vector vecMAFs = new Vector();
        for (int s = 0; s < snps.length; s++) {
            if (snpToSNPId.get(snps[s]) != -9) {
                int id = snpToSNPId.get(snps[s]);
                SNP snp = genotypeData.getSNPObject(id);
                synchronized (snpLoader) {
                    snpLoader.loadGenotypes(snp);
                }
                Double maf = snp.getMAF();
                vecMAFs.add(maf);
                System.out.println(snp.getName() + "\t" + maf);
            } else {
                System.out.println("Error! SNP " + snps[s] + " is not present in the genotype data!!!!!!");
            }
        }


        double[] mafs = new double[vecMAFs.size()];
        for (int m = 0; m < mafs.length; m++) {
            mafs[m] = ((Double) vecMAFs.get(m)).doubleValue();
        }
        double median = JSci.maths.ArrayMath.median(mafs);
        System.out.println("Median MAF:\t" + median);

        System.out.println("\n\n\n\n\n\n");
        for (int snpId = 0, len = snpToSNPId.size(); snpId < len; snpId++) {
            SNP snp = genotypeData.getSNPObject(snpId);
            synchronized (snpLoader) {

                snpLoader.loadGenotypes(snp);
            }
            Double maf = snp.getMAF();
            System.out.println(snpId + "\t" + snp.getName() + "\t" + maf);
        }

    }

    /**
     *
     * Independifies given SNPs with given thresholds.
     *
     * @param snps
     * @param r2Threshold
     * @param bpThreshold Base pair distance threshold, ld will be calculated
     * for SNPs
     * @return
     * @throws IOException
     */
    public String[] independify(String[] snps, double r2Threshold, int bpThreshold) throws IOException {

        // limit to SNPs present in genotype data
        TObjectIntHashMap<String> snpToSNPId = genotypeData.getSnpToSNPId();
        Vector vecSNPsPresentInGenotypeData = new Vector();
        for (int s = 0; s < snps.length; s++) {
            if (snpToSNPId.get(snps[s]) == -9) {
                vecSNPsPresentInGenotypeData.add(snps[s]);
            } else {
                System.out.println("Error! SNP " + snps[s] + " is not present in the genotype data!!!!!!");
            }
        }
        snps = new String[vecSNPsPresentInGenotypeData.size()];
        for (int s = 0; s < snps.length; s++) {
            snps[s] = (String) vecSNPsPresentInGenotypeData.get(s);
        }

        // calculate pairwise LD over all SNPs in the same chromosome within the given distance
        DetermineLD ldCalc = new DetermineLD();
        double[][] r2Matrix = new double[snps.length][snps.length];
        int[] snpsChr = new int[snps.length];
        int[] snpsChrPos = new int[snps.length];
        for (int s1 = 0; s1 < snps.length; s1++) {
            Integer get = snpToSNPId.get(snps[s1]);
            SNP snp1 = genotypeData.getSNPObject(get);
            synchronized (snpLoader) {
                snpLoader.loadGenotypes(snp1);
            }
            byte chr1 = snp1.getChr();
            int chrPos1 = snp1.getChrPos();
            snpsChr[s1] = chr1;
            snpsChrPos[s1] = chrPos1;
            for (int s2 = s1 + 1; s2 < snps.length; s2++) {
                SNP snp2 = genotypeData.getSNPObject(snpToSNPId.get(snps[s2]));
                if (s1 != s2 && chr1 == snp2.getChr() && Math.abs(chrPos1 - snp2.getChrPos()) < bpThreshold) {
                    synchronized (snpLoader) {
                        snpLoader.loadGenotypes(snp2);
                    }
                    double r2 = ldCalc.getRSquared(snp1, snp2, genotypeData, DetermineLD.RETURN_R_SQUARED, DetermineLD.INCLUDE_CASES_AND_CONTROLS, false);
                    r2Matrix[s1][s2] = r2;
                    r2Matrix[s2][s1] = r2;
                }
            }
        }

        int nrIndependifiedLeadSNPs = 0;
        int nrSNPsToUse = snps.length;
        System.out.println("Nr of SNPs with available genotype data in imputed dataset:\t" + nrSNPsToUse);

        int[] clusters = new int[snps.length];
        for (int s1 = 0; s1 < snps.length; s1++) {
            clusters[s1] = -1;
        }
        int nrClusters = 0;

        // assign SNPs to clusters
        for (int s1 = 0; s1 < nrSNPsToUse; s1++) {
            boolean SNPInLDWithOtherSNP = false;
            for (int s2 = 0; s2 < nrSNPsToUse; s2++) {
                if (s1 != s2 && snpsChr[s1] == snpsChr[s2] && Math.abs(snpsChrPos[s1] - snpsChrPos[s2]) < bpThreshold) {
                    double r2 = r2Matrix[s1][s2];
                    if (r2 > r2Threshold) {
                        SNPInLDWithOtherSNP = true;
                        if (clusters[s1] == -1 && clusters[s2] == -1) {
                            //Both SNPs have not yet been assigned to any cluster, easy!
                            clusters[s1] = nrClusters;
                            clusters[s2] = nrClusters;
                            nrClusters++;
                        } else {
                            if (clusters[s1] != -1 && clusters[s2] != -1) {
                                //Both SNPS have already been assigned to clusters, merge these clusters:
                                int previousClusterS1 = clusters[s1];
                                int previousClusterS2 = clusters[s2];
                                for (int q = 0; q < nrSNPsToUse; q++) {
                                    if (clusters[q] == previousClusterS1) {
                                        clusters[q] = nrClusters;
                                    }
                                    if (clusters[q] == previousClusterS2) {
                                        clusters[q] = nrClusters;
                                    }
                                }
                                nrClusters++;
                            } else {
                                if (clusters[s1] == -1) {
                                    clusters[s1] = clusters[s2];
                                }
                                if (clusters[s2] == -1) {
                                    clusters[s2] = clusters[s1];
                                }
                            }
                        }
                    }
                }
            }
            if (!SNPInLDWithOtherSNP) {
                clusters[s1] = nrClusters;
                nrClusters++;
            }
        }

        // build string array representation of SNP clusters
        // the first SNP in each cluster will be the most significant one if the SNPs were given in order of significance
        List<String> independifiedSNPs = new ArrayList<String>();
        nrIndependifiedLeadSNPs = 0;
        for (int c = 0; c < nrClusters; c++) {
            boolean clusterContainsSNPs = false;
            for (int s1 = 0; s1 < nrSNPsToUse; s1++) {
                if (clusters[s1] == c) {
                    clusterContainsSNPs = true;
                    break;
                }
            }
            if (clusterContainsSNPs) {
                nrIndependifiedLeadSNPs++;
                String snpsThisCluster = "";
                String delim = "";
                for (int s1 = 0; s1 < nrSNPsToUse; s1++) {
                    if (clusters[s1] == c) {
                        snpsThisCluster += delim + snps[s1];
                        delim = ";";
                    }
                }
                independifiedSNPs.add(snpsThisCluster);
            }
        }

        return independifiedSNPs.toArray(new String[0]);
    }

    public String[] independify(String permutationDataset, int permutation, double pThreshold, double r2Threshold, int bpThreshold, int nrIndependentSNPsWanted) throws IOException {
        return independify(permutationDataset, permutation, pThreshold, r2Threshold, bpThreshold, nrIndependentSNPsWanted, null);
    }

    public void addPermutationDataset(String name, String filename) throws IOException {
        DoubleMatrixDataset<String, String> permDataset = new DoubleMatrixDataset<String, String>(filename);
        if (permDatasets == null) {
            permDatasets = new HashMap<String, DoubleMatrixDataset<String, String>>();
        }
        permDatasets.put(name, permDataset);
        if (filename.endsWith(".txt")) {
            permDataset.save(filename.replaceAll(".txt", ".binary"));
        }
    }

    public String[] independify(String permutationDataset, int permutation, double pThreshold, double r2Threshold, int bpThreshold, int nrIndependentSNPsWanted, String outFile) throws IOException {

        DoubleMatrixDataset<String, String> permDataset = permDatasets.get(permutationDataset);
        if (permDataset == null) {
            throw new IllegalStateException("Dataset '" + permutationDataset + "' not set!");
        }

        TObjectIntHashMap<String> snpToSNPId = genotypeData.getSnpToSNPId();
        Vector vecTopResults = new Vector();
        for (int s = 0; s < permDataset.nrCols; s++) {
            if (permDataset.rawData[permutation][s] < pThreshold) {
                if (snpToSNPId.get(permDataset.colObjects.get(s)) == -9) {
//                    System.err.println("SNP not in data! " + permDataset.colObjects.get(s));
                    continue;
                }
                Byte chr = genotypeData.getChr(s);
                if (chr < 0 || chr > 23) {
                    // bad SNP is bad
                } else {
                    umcg.genetica.containers.StringDoubleObject e = new umcg.genetica.containers.StringDoubleObject(permDataset.colObjects.get(s), permDataset.rawData[permutation][s]);
                    vecTopResults.add(e);
                }
            }
        }
        umcg.genetica.util.StringDoubleObjectSorterSortOnDouble sorter = new umcg.genetica.util.StringDoubleObjectSorterSortOnDouble();
        sorter.sort(vecTopResults);

        String[] snps = new String[vecTopResults.size()];
        for (int q = 0; q < vecTopResults.size(); q++) {
            snps[q] = ((umcg.genetica.containers.StringDoubleObject) vecTopResults.get(q)).stringValue;
        }

        LOGGER.log(Level.INFO, "Permutation {0}:\tUsing {1} SNPs, this will require {2} MB of memory", new Object[]{permutation, snps.length, 4 * snps.length * (snps.length / 1024d / 1024d)});

        DetermineLD ldCalc = new DetermineLD();
        float[][] r2Matrix = new float[snps.length][snps.length];
        int[] snpsChr = new int[snps.length];
        int[] snpsChrPos = new int[snps.length];
        for (int s1 = 0; s1 < snps.length; s1++) {
            Integer get = snpToSNPId.get(snps[s1]);
            SNP snp1 = genotypeData.getSNPObject(get);
            synchronized (snpLoader) {
                snpLoader.loadGenotypes(snp1);
            }
            byte chr1 = snp1.getChr();
            int chrPos1 = snp1.getChrPos();
            snpsChr[s1] = chr1;
            snpsChrPos[s1] = chrPos1;
            for (int s2 = s1 + 1; s2 < snps.length; s2++) {
                SNP snp2 = genotypeData.getSNPObject(snpToSNPId.get(snps[s2]));
                if (s1 != s2 && chr1 == snp2.getChr() && Math.abs(chrPos1 - snp2.getChrPos()) < bpThreshold) {
                    synchronized (snpLoader) {
                        snpLoader.loadGenotypes(snp2);
                    }
                    double r2 = ldCalc.getRSquared(snp1, snp2, genotypeData, DetermineLD.RETURN_R_SQUARED, DetermineLD.INCLUDE_CASES_AND_CONTROLS, false);
                    r2Matrix[s1][s2] = (float) r2;
                    r2Matrix[s2][s1] = (float) r2;
                }
                snp2.clearGenotypes();
            }
            snp1.clearGenotypes();
        }

        int nrIndependifiedLeadSNPs = 0;
        int nrSNPsToUse = snps.length;

        int lower = 0;
        int upper = snps.length;
        List<String> independifiedSNPs = new ArrayList<String>();

        while (nrIndependifiedLeadSNPs != nrIndependentSNPsWanted) {

            nrSNPsToUse = (upper + lower) / 2;

            int[] clusters = new int[snps.length];
            for (int s1 = 0; s1 < snps.length; s1++) {
                clusters[s1] = -1;
            }
            int nrClusters = 0;

            for (int s1 = 0; s1 < nrSNPsToUse; s1++) {
                boolean SNPInLDWithOtherSNP = false;
                for (int s2 = 0; s2 < nrSNPsToUse; s2++) {
                    if (s1 != s2 && snpsChr[s1] == snpsChr[s2] && Math.abs(snpsChrPos[s1] - snpsChrPos[s2]) < bpThreshold) {
                        double r2 = r2Matrix[s1][s2];
                        if (r2 > r2Threshold) {
                            SNPInLDWithOtherSNP = true;
                            if (clusters[s1] == -1 && clusters[s2] == -1) {
                                //Both SNPs have not yet been assigned to any cluster, easy!
                                clusters[s1] = nrClusters;
                                clusters[s2] = nrClusters;
                                nrClusters++;
                            } else {
                                if (clusters[s1] != -1 && clusters[s2] != -1) {
                                    //Both SNPS have already been assigned to clusters, merge these clusters:
                                    int previousClusterS1 = clusters[s1];
                                    int previousClusterS2 = clusters[s2];
                                    for (int q = 0; q < nrSNPsToUse; q++) {
                                        if (clusters[q] == previousClusterS1) {
                                            clusters[q] = nrClusters;
                                        }
                                        if (clusters[q] == previousClusterS2) {
                                            clusters[q] = nrClusters;
                                        }
                                    }
                                    nrClusters++;
                                } else {
                                    if (clusters[s1] == -1) {
                                        clusters[s1] = clusters[s2];
                                    }
                                    if (clusters[s2] == -1) {
                                        clusters[s2] = clusters[s1];
                                    }
                                }
                            }
                        }
                    }
                }
                if (!SNPInLDWithOtherSNP) {
                    clusters[s1] = nrClusters;
                    nrClusters++;
                }
            }

            independifiedSNPs = new ArrayList<String>();
            TextFile out = null;
            if (outFile != null && !outFile.isEmpty()) {
                out = new TextFile(outFile, true);
            }
            nrIndependifiedLeadSNPs = 0;
            for (int c = 0; c < nrClusters; c++) {
                boolean clusterContainsSNPs = false;
                for (int s1 = 0; s1 < nrSNPsToUse; s1++) {
                    if (clusters[s1] == c) {
                        clusterContainsSNPs = true;
                        break;
                    }
                }
                if (clusterContainsSNPs) {
                    nrIndependifiedLeadSNPs++;
                    if (out != null) {
                        out.write(c + "");
                    }
                    String snpsThisCluster = "";
                    String delim = "";
                    for (int s1 = 0; s1 < nrSNPsToUse; s1++) {
                        if (clusters[s1] == c) {
                            snpsThisCluster += delim + snps[s1];
                            delim = ";";
                            if (out != null) {
                                out.write("\t" + snps[s1]);
                            }
                        }
                    }
                    independifiedSNPs.add(snpsThisCluster);
                    if (out != null) {
                        out.writeln();
                    }
                }
            }
            if (out != null) {
                out.close();
            }

//            System.out.println("Number of SNPs to include:\t" + snps.length + "\t" + nrIndependifiedLeadSNPs);
            if (nrIndependifiedLeadSNPs != nrIndependentSNPsWanted) {
                if (nrIndependifiedLeadSNPs < nrIndependentSNPsWanted) {
                    lower = nrSNPsToUse;
                } else {
                    upper = nrSNPsToUse;
                }
            }
            if (upper - lower <= 1) {
                break;
            }

        }

        return independifiedSNPs.toArray(new String[0]);
    }
}
