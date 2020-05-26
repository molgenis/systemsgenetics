/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.util;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import umcg.genetica.containers.Pair;
import umcg.genetica.containers.SortableSNP;
import umcg.genetica.containers.Triple;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.SNPLoader;
import umcg.genetica.io.trityper.TriTyperGenotypeData;
import umcg.genetica.io.trityper.util.DetermineLD;

/**
 *
 * @author harmjan
 */
public class LDCalculator {

    public static void calculatePairwiseLD(String snpfile, String dataset, String output) throws IOException {
        if (snpfile == null) {
            throw new IllegalArgumentException("Error: snp file location should be provided");
        }
        if (dataset == null) {
            throw new IllegalArgumentException("Error: dataset location should be provided");
        }
        if (output == null) {
            throw new IllegalArgumentException("Error: output location should be provided");
        }



        TextFile snpTextFile = new TextFile(snpfile, TextFile.R);

        String[] elems = snpTextFile.readLineElems(TextFile.tab);
        ArrayList<Pair<String, String>> pairs = new ArrayList<Pair<String, String>>();
        while (elems != null) {
            if (elems.length > 1) {
                Pair<String, String> p = new Pair<String, String>(elems[0], elems[1]);
                pairs.add(p);
            }
            elems = snpTextFile.readLineElems(TextFile.tab);
        }

        snpTextFile.close();

        if (pairs.isEmpty()) {
            throw new IllegalArgumentException("Error: " + snpfile + " does not contain any useable SNP pairs. Did you use tab separation?");
        }

        DetermineLD calc = new DetermineLD();
        TriTyperGenotypeData ds = new TriTyperGenotypeData();

        ds.load(dataset);
        SNPLoader loader = ds.createSNPLoader();
        TextFile outputfile = new TextFile(output, TextFile.W);
        String header = "SNP1\tSNP1Chr\tSNP1ChrPos\tSNP1MAF\tSNP1CR\tSNP1HWEP"
                + "\tSNP2\tSNP2Chr\tSNP2ChrPos\tSNP2MAF\tSNP2CR\tSNP2HWEP"
                + "\tDistrance\tR2\tSNPOnSameChr";
        outputfile.writeln(header);
        for (Pair<String, String> p : pairs) {
            Integer snp1Id = ds.getSnpToSNPId().get(p.getLeft());
            Integer snp2Id = ds.getSnpToSNPId().get(p.getRight());

            boolean test = true;
            String teststr = "";
            if (snp1Id == -9) {
                System.out.println(p.getLeft() + "\tnot present in dataset");
                teststr = p.getLeft() + "\tnot present in dataset";
                test = false;
            }

            if (snp2Id == -9) {
                System.out.println(p.getRight() + "\tnot present in dataset");
                if (teststr.length() > 0) {
                    teststr += "\t" + p.getRight() + "\tnot present in dataset";
                } else {
                    teststr = p.getRight() + "\tnot present in dataset";
                }

                test = false;
            }

            if (!test) {
                outputfile.writeln(teststr);
            } else {
                SNP snp1Obj = ds.getSNPObject(snp1Id);
                SNP snp2Obj = ds.getSNPObject(snp2Id);

                boolean snpOnSameChr = true;
                if (snp1Obj.getChr() != snp2Obj.getChr()) {
                    System.out.println(p.getLeft() + "\t" + snp1Obj.getChr() + "\t" + p.getRight() + "\t" + snp2Obj.getChr() + "\tnot on same chr");
//                    outputfile.writeln(p.getLeft() + "\t" + snp1Obj.getChr() + "\t" + p.getRight() + "\t" + snp2Obj.getChr() + "\tnot on same chr");
                    snpOnSameChr = false;
                }
//                else {
                loader.loadGenotypes(snp1Obj);
                loader.loadGenotypes(snp2Obj);

                double r2 = calc.getRSquared(snp1Obj, snp2Obj, ds, DetermineLD.RETURN_R_SQUARED, DetermineLD.INCLUDE_CASES_AND_CONTROLS, false);

                int distance = snp2Obj.getChrPos() - snp1Obj.getChrPos();

                String ln = p.getLeft() + "\t" + snp1Obj.getChr() + "\t" + snp1Obj.getChrPos() + "\t" + snp1Obj.getMAF() + "\t" + snp1Obj.getCR() + "\t" + snp1Obj.getHWEP() + "\t" + p.getRight() + "\t" + snp2Obj.getChr() + "\t" + snp2Obj.getChrPos() + "\t" + snp2Obj.getMAF() + "\t" + snp2Obj.getCR() + "\t" + snp2Obj.getHWEP() + "\t" + distance + "\t" + r2 + "\t" + snpOnSameChr;
                System.out.println(ln);
                outputfile.writeln(ln);
                snp1Obj.clearGenotypes();
                snp2Obj.clearGenotypes();
//                }
            }
        }
        outputfile.close();
    }

    public static void proxyLookUpInReferenceDataset(String reference, String inputList, double mafThreshold, double hwepthreshold, double crthreshold, double r2threshold, String output, int maxdistance) throws IOException {
        TextFile snpList = new TextFile(inputList, TextFile.R);
        String[] snpsToQuery = snpList.readAsArray();
        snpList.close();

        TriTyperGenotypeData refDataset = new TriTyperGenotypeData();
        refDataset.load(reference);
        SNPLoader refLoader = refDataset.createSNPLoader();



        TextFile outputfile = new TextFile(output, TextFile.W);
        // first check which SNPs are actually already in the dataset. We don't have to look for proxies for those SNPs, considering the threshold set as arguments here.
        int ctr = 0;
        for (String snp : snpsToQuery) {
            String outputStr = snp;
            System.out.println(ctr + "/" + snpsToQuery.length);
            Integer refDatasetId = refDataset.getSnpToSNPId().get(snp);
            if (refDatasetId == -9) {
                // this SNP is useless, we cannot find a proxy for a SNP we have no data for.
                outputStr += "\t-\tNo proxy available - Not present in reference";
            } else {
                // possibly there is a proxy in our reference.
                Triple<String, Double, Integer> proxy = proxyLookUpInReferenceDataAlone(refDataset, refLoader, refDatasetId, r2threshold, mafThreshold, crthreshold, hwepthreshold, maxdistance);
                if (proxy == null) {
                    outputStr += "\t-\tNo proxy available - no SNP that meets thresholds";
                } else {
                    outputStr += "\t" + proxy.getLeft() + "\t" + proxy.getMiddle() + "\t" + proxy.getRight();
                }
            }
            outputfile.writeln(outputStr);
            System.out.println(outputStr);
            ctr++;
        }

        outputfile.close();
        refLoader.close();
    }

    public static void proxyLookupInReferenceThatIsAlsoIneQTLDataset(String reference, String dataset, String inputList, double mafThreshold, double hwepthreshold, double crthreshold, double r2threshold, String output, int maxdistance) throws IOException {
        TextFile snpList = new TextFile(inputList, TextFile.R);
        String[] snpsToQuery = snpList.readAsArray();
        snpList.close();


        TriTyperGenotypeData inDataset = new TriTyperGenotypeData();
        inDataset.load(dataset);
        SNPLoader datasetLoader = inDataset.createSNPLoader();



        TriTyperGenotypeData refDataset = new TriTyperGenotypeData();
        refDataset.load(reference);
        SNPLoader refLoader = refDataset.createSNPLoader();



        TextFile outputfile = new TextFile(output, TextFile.W);
        // first check which SNPs are actually already in the dataset. We don't have to look for proxies for those SNPs, considering the threshold set as arguments here.
        for (String snp : snpsToQuery) {
            String outputStr = snp;
            Integer inDatasetId = inDataset.getSnpToSNPId().get(snp);
            Integer refDatasetId = refDataset.getSnpToSNPId().get(snp);
            if (inDatasetId == -9 && refDatasetId == -9) {
                // this SNP is useless, we cannot find a proxy for a SNP we have no data for.
                outputStr += "\t-\tNo proxy available - Not present in reference and dataset";
            } else if (inDatasetId == -9 && refDatasetId != -9) {
                // possibly there is a proxy in our reference.
                Triple<String, Double, Integer> proxy = findProxyThatIsAlsoInDatasetToTest(refDataset, inDataset, refLoader, datasetLoader, refDatasetId, r2threshold, mafThreshold, crthreshold, hwepthreshold, maxdistance);
                if (proxy == null) {
                    outputStr += "\t-\tNo proxy available - no SNP that meets thresholds";
                } else {
                    outputStr += "\t" + proxy.getLeft() + "\t" + proxy.getMiddle() + "\t" + proxy.getRight();
                }
            } else if (inDatasetId != -9 && refDatasetId != -9) {
                // we might need to look for a proxy, if the SNP in the inDataset does not conform to the thresholds set.

                // we need to load the data first
                SNP snpInDataObj = inDataset.getSNPObject(inDatasetId);
                datasetLoader.loadGenotypes(snpInDataObj);
                snpInDataObj.clearGenotypes();

                if (snpInDataObj.getMAF() < mafThreshold || snpInDataObj.getCR() < crthreshold || snpInDataObj.getHWEP() < hwepthreshold) {
                    // we have to look for a proxy
                    Triple<String, Double, Integer> proxy = findProxyThatIsAlsoInDatasetToTest(refDataset, inDataset, refLoader, datasetLoader, refDatasetId, r2threshold, mafThreshold, crthreshold, hwepthreshold, maxdistance);
                    if (proxy == null) {
                        outputStr += "\t-\tNo proxy available - no SNP that meets thresholds";
                    } else {
                        outputStr += "\t" + proxy.getLeft() + "\t" + proxy.getMiddle() + "\t" + proxy.getRight();
                    }
                } else {
                    // its perfectly fine to use this SNP.
                    outputStr += "\t" + snp + "\tNo proxy required";
                }



            } else if (inDatasetId != -9 && refDatasetId == -9) {
                SNP snpInDataObj = inDataset.getSNPObject(inDatasetId);
                datasetLoader.loadGenotypes(snpInDataObj);
                snpInDataObj.clearGenotypes();
                if (snpInDataObj.getMAF() < mafThreshold || snpInDataObj.getCR() < crthreshold || snpInDataObj.getHWEP() < hwepthreshold) {
                    outputStr += "\t-\tNo proxy available - SNP not present in reference";
                } else {
                    outputStr += "\t" + snp + "\tNo proxy required - SNP not present in reference but meets eQTL QC thresholds";
                }

            }
            outputfile.writeln(outputStr);
            System.out.println(outputStr);
        }

        outputfile.close();
        datasetLoader.close();
        refLoader.close();
    }

    // find proxy in reference that is also in the inDataset
    // so we only have to check SNPs that are also in the inDataset
    private static Triple<String, Double, Integer> findProxyThatIsAlsoInDatasetToTest(TriTyperGenotypeData refDataset, TriTyperGenotypeData inDataset, SNPLoader refLoader, SNPLoader datasetLoader, Integer snpIdInReference, double r2threshold,
            double maf, double cr, double hwep, int maxDistance) throws IOException {
        String[] snpsInReference = refDataset.getSNPs();
        byte chromosomeOfStartSNP = refDataset.getChr(snpIdInReference);

        String query = snpsInReference[snpIdInReference];
        System.out.println("");
        System.out.println("Query:\t" + query);
        ArrayList<SortableSNP> snpsToSort = new ArrayList<SortableSNP>();


        for (int s = 0; s < snpsInReference.length; s++) {
            if (refDataset.getChr(s) == chromosomeOfStartSNP && !snpIdInReference.equals(s)) {
                // now check whether the SNP is actually in the in dataset.
                String snpName = snpsInReference[s];
                int chrPos = refDataset.getChrPos(s);
                int distance = Math.abs(refDataset.getChrPos(s) - refDataset.getChrPos(snpIdInReference));
                if (distance < maxDistance) {
                    if (inDataset.getSnpToSNPId().get(snpName) != -9) {
                        snpsToSort.add(new SortableSNP(snpName, s, chromosomeOfStartSNP, refDataset.getChrPos(s), SortableSNP.SORTBY.CHRPOS));
                    }
                }

            }
        }

        System.out.println(snpsToSort.size() + "\t SNPs within 10Mb");


        if (snpsToSort.isEmpty()) {
            return null;
        }
        // sort the SNPs
        Collections.sort(snpsToSort);

        // now calculate the LD with the startSNP.
        // we want the closest SNP, so first find the snp that we started with

        SortableSNP[] sorted = snpsToSort.toArray(new SortableSNP[0]);

        SNP startSNPObj = refDataset.getSNPObject(snpIdInReference);
        refLoader.loadGenotypes(startSNPObj);
        DetermineLD ldcalc = new DetermineLD();
        Integer bestproxy = null;

        double maxR2 = Double.MIN_VALUE;
        int dist = 0;
        for (int s = 0; s < sorted.length; s++) {
            int nextSNP = sorted[s].id;
            SNP nextSNPObj = refDataset.getSNPObject(nextSNP);
            refLoader.loadGenotypes(nextSNPObj);

            double r2 = ldcalc.getRSquared(startSNPObj, nextSNPObj, refDataset, DetermineLD.RETURN_R_SQUARED, DetermineLD.INCLUDE_CASES_AND_CONTROLS, false);
            double dp = ldcalc.getRSquared(startSNPObj, nextSNPObj, refDataset, DetermineLD.RETURN_D_PRIME, DetermineLD.INCLUDE_CASES_AND_CONTROLS, false);
            int distance = nextSNPObj.getChrPos() - startSNPObj.getChrPos();

            if (r2 >= r2threshold) {

                // possible candidate. Now check for MAF, CR and HWEP in the real dataset.
                Integer snpIdInDataset = inDataset.getSnpToSNPId().get(nextSNPObj.getName());
                SNP snpInDatasetObj = inDataset.getSNPObject(snpIdInDataset);
                datasetLoader.loadGenotypes(snpInDatasetObj);
                boolean snpPassesQC = false;
                if (snpInDatasetObj.getMAF() > maf && snpInDatasetObj.getHWEP() > hwep && snpInDatasetObj.getCR() > cr && r2 >= maxR2) {
                    // this is a hit!
                    bestproxy = s;
                    maxR2 = r2;
                    dist = distance;
                    snpPassesQC = true;
//                    System.out.println("Candidate proxy:\t" + nextSNPObj.getName() + "\t" + r2);
                }
//                if ((r2 > 0.1 || dp > 0.1) && snpPassesQC) {
//                    System.out.println("Testing LD\t" + startSNPObj.getName() + "\t" + startSNPObj.getChr() + "\t" + nextSNPObj.getName() + "\t" + nextSNPObj.getChr() + "\t" + distance + "\t" + r2 + "\t" + dp + "\t" + snpPassesQC
//                            + "\t" + snpInDatasetObj.getMAF() + "\t" + snpInDatasetObj.getHWEP() + "\t" + snpInDatasetObj.getCR());
//                }
                snpInDatasetObj.clearGenotypes();
            }
            nextSNPObj.clearGenotypes();
        }

        startSNPObj.clearGenotypes();

        if (bestproxy == null) {
            return null;
        }
        return new Triple<String, Double, Integer>(sorted[bestproxy].name, maxR2, dist); // returns null if nothing found
    }

    private static Triple<String, Double, Integer> proxyLookUpInReferenceDataAlone(TriTyperGenotypeData refDataset, SNPLoader refLoader, Integer refDatasetSNPId, double r2threshold, double maf, double cr, double hwep, int maxdistance) throws IOException {
        String[] snpsInReference = refDataset.getSNPs();
        byte chromosomeOfStartSNP = refDataset.getChr(refDatasetSNPId);

        String query = snpsInReference[refDatasetSNPId];
        System.out.println("");
        System.out.println("Query:\t" + query);
        ArrayList<SortableSNP> snpsToSort = new ArrayList<SortableSNP>();


        for (int s = 0; s < snpsInReference.length; s++) {
            if (refDataset.getChr(s) == chromosomeOfStartSNP && !refDatasetSNPId.equals(s)) {
                // now check whether the SNP is actually in the in dataset.
                String snpName = snpsInReference[s];
                int chrPos = refDataset.getChrPos(s);
                if (chrPos != -1) {
                    int distance = Math.abs(refDataset.getChrPos(s) - refDataset.getChrPos(refDatasetSNPId));
                    if (distance < maxdistance) {
                        snpsToSort.add(new SortableSNP(snpName, s, chromosomeOfStartSNP, refDataset.getChrPos(s), SortableSNP.SORTBY.CHRPOS));
                    }
                }
            }
        }

        System.out.println(snpsToSort.size() + "\t SNPs within 10Mb");


        if (snpsToSort.isEmpty()) {
            return null;
        }
        // sort the SNPs
        Collections.sort(snpsToSort);

        // now calculate the LD with the startSNP.
        // we want the closest SNP, so first find the snp that we started with

        
        SortableSNP[] sorted = snpsToSort.toArray(new SortableSNP[0]);

        SNP startSNPObj = refDataset.getSNPObject(refDatasetSNPId);
        refLoader.loadGenotypes(startSNPObj);
        DetermineLD ldcalc = new DetermineLD();
        Integer bestproxy = null;

        double maxR2 = Double.MIN_VALUE;
        int dist = 0;
        for (int s = 0; s < sorted.length; s++) {
            int nextSNP = sorted[s].id;
            SNP nextSNPObj = refDataset.getSNPObject(nextSNP);
            refLoader.loadGenotypes(nextSNPObj);

            Pair<Double, Double> ld = ldcalc.getLD(startSNPObj, nextSNPObj, refDataset, DetermineLD.INCLUDE_CASES_AND_CONTROLS, false);
            
            double r2 = ld.getLeft();
            double dp = ld.getRight();
            
            int distance = nextSNPObj.getChrPos() - startSNPObj.getChrPos();

//            System.out.println(startSNPObj.getName() + "\t" + nextSNPObj.getName() + "\t" + distance + "\t" + r2 + "\t" + dp);

            if (r2 >= r2threshold) {
                // possible candidate. Now check for MAF, CR and HWEP in the real dataset.
                if (nextSNPObj.getMAF() > maf && nextSNPObj.getHWEP() > hwep && nextSNPObj.getCR() > cr && r2 >= maxR2) {
                    // this is a hit!
                    bestproxy = s;
                    maxR2 = r2;
                    dist = distance;
                }
            }
            nextSNPObj.clearGenotypes();
        }

        startSNPObj.clearGenotypes();

        if (bestproxy == null) {
            return null;
        }
        return new Triple<String, Double, Integer>(sorted[bestproxy].name, maxR2, dist); // returns null if nothing found
    }
}
