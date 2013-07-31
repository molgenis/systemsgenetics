/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.gwas;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.SNPLoader;
import umcg.genetica.io.trityper.TriTyperGenotypeData;
import umcg.genetica.io.trityper.util.DetermineLD;

/**
 *
 * @author harm-jan
 */
public class Dependifier {

    private final SNPLoader loader;
    private final TriTyperGenotypeData genotypeData;
    private final DetermineLD ldCalc = new DetermineLD();
    private final String[] allSNPsInReference;
    private HashSet<String> haystack;

    public Dependifier(String datadir) throws IOException {
        genotypeData = new TriTyperGenotypeData();
        genotypeData.load(datadir);
        loader = genotypeData.createSNPLoader();
        allSNPsInReference = genotypeData.getSNPs();
    }

    public Dependifier(TriTyperGenotypeData dataset, SNPLoader loader) {
        this.loader = loader;
        this.genotypeData = dataset;
        allSNPsInReference = genotypeData.getSNPs();
    }

    public HashMap<String, HashSet<String>> dependifyReturnProxiesPerSNP(String[] inputsnps, double proxyldthreshold, int maxproxydistance) throws IOException {
        HashMap<String, HashSet<String>> proxies = new HashMap<String, HashSet<String>>();
        for (String snp : inputsnps) {
            // include the query itself as well in the final proxy list..
            proxies.put(snp, findProxiesForSNP(snp, proxyldthreshold, maxproxydistance));
        }
        return proxies;
    }

    public HashSet<String> dependify(String[] inputsnps, double proxyldthreshold, int maxproxydistance) throws IOException {
        HashSet<String> proxies = new HashSet<String>();
        for (String snp : inputsnps) {
            // include the query itself as well in the final proxy list..
            proxies.addAll(findProxiesForSNP(snp, proxyldthreshold, maxproxydistance));
        }
        return proxies;
    }

    // find SNPs in needles that are in LD with haystack
    public HashSet<String> dependify(String[] needles, String[] haystack, double proxyldthreshold, int maxproxydistance) throws IOException {
        HashSet<String> proxies = new HashSet<String>();
        this.haystack = new HashSet<String>();
        this.haystack.addAll(Arrays.asList(haystack));
        for (String snp : needles) {
            proxies.addAll(findProxiesForSNP(snp, proxyldthreshold, maxproxydistance));
        }
        this.haystack = null; // clean this up after use..
        return proxies;
    }

    public HashMap<String, HashSet<String>> dependifyReturnProxiesPerSNP(String[] needles, String[] haystack, double proxyldthreshold, int maxproxydistance) throws IOException {
        HashMap<String, HashSet<String>> proxies = new HashMap<String, HashSet<String>>();
        this.haystack = new HashSet<String>();
        this.haystack.addAll(Arrays.asList(haystack));
        for (String snp : needles) {
            proxies.put(snp, findProxiesForSNP(snp, proxyldthreshold, maxproxydistance));
        }
        this.haystack = null; // clean this up after use..
        return proxies;
    }

    public HashSet<String> findProxiesForSNP(String snp, double proxyldthreshold, int maxproxydistance) throws IOException {
        HashSet<String> proxies = new HashSet<String>();
        Integer snpIdInReference = genotypeData.getSnpToSNPId().get(snp);
        if (snpIdInReference != null) {

            SNP snpObj = genotypeData.getSNPObject(snpIdInReference);
            byte chr = snpObj.getChr();
            int chrPos = snpObj.getChrPos();

            // if the SNP has a valid mapping..
            if (chr > 0 && chrPos > 0) {

                HashSet<Integer> allSNPsWithinMaxDistance = new HashSet<Integer>();
                for (int i = 0; i < allSNPsInReference.length; i++) {
                    // check whether we want to limit our search within a preset list of SNPs
                    if (haystack == null || haystack.contains(allSNPsInReference[i])) {
                        byte chr2 = genotypeData.getChr(i);
                        if (chr == chr2) {
                            int chrPos2 = genotypeData.getChrPos(i);
                            if (Math.abs(chrPos2 - chrPos) < maxproxydistance) {
                                allSNPsWithinMaxDistance.add(i);
                            }
                        }
                    }
                }
                // if there are any snps to test within maxdistance
                if (!allSNPsWithinMaxDistance.isEmpty()) {
                    loader.loadGenotypes(snpObj);
                    for (Integer i : allSNPsWithinMaxDistance) {
                        SNP snpObj2 = genotypeData.getSNPObject(i);
                        loader.loadGenotypes(snpObj2);
                        if (snpObj.getMAF() > 0 && snpObj2.getMAF() > 0) {
                            double r2 = ldCalc.getRSquared(snpObj, snpObj2, genotypeData, DetermineLD.RETURN_R_SQUARED, DetermineLD.INCLUDE_CASES_AND_CONTROLS, false);
                            if (!Double.isNaN(r2) && r2 >= proxyldthreshold) {
                                proxies.add(snpObj2.getName());
                            }
                        }
                        snpObj2.clearGenotypes();
                    }
                    snpObj.clearGenotypes();
                }
            }
            snpObj = null;
        }
        return proxies;
    }
}
