package mbqtl;

import mbqtl.vcf.DetermineLDVCF;
import mbqtl.vcf.VCFTabix;
import mbqtl.vcf.VCFVariant;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.containers.Pair;
import umcg.genetica.enums.Chromosome;
import umcg.genetica.features.Feature;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.*;
import java.util.stream.IntStream;

public class GWASLDCalculator {


    private double hwepthreshold = 0.0001;
    private double mafthreshold = 0.01;
    private double ldthreshold = 0.8;
    private int prunedistance = 5000000;
    private double pruneLDthreshold = 0.2;

    public void runPerGWAS3(String eqtlset, String eqtlfile, String gwasset, String gwasAssocFile, String gwasListFile, String vcfPrefixString, String individualSubset, int prunedistance, double prunethreshold, double ldthreshold, String outputfile, boolean skipchr6, boolean matchbyrsid) throws IOException {
        System.out.println("eqtlset: " + eqtlset);
        System.out.println("eqtlfile: " + eqtlfile);
        System.out.println("gwasset: " + gwasset);
        System.out.println("gwasassocfile: " + gwasAssocFile);
        System.out.println("gwaslistfile: " + gwasListFile);
        System.out.println("gwaslistfile: " + vcfPrefixString);
        System.out.println("gte: " + individualSubset);
        System.out.println("prunedistance: " + prunedistance);
        System.out.println("prunethreshold: " + prunethreshold);
        System.out.println("ldthreshold: " + ldthreshold);
        System.out.println("out: " + outputfile);
        System.out.println("skipchr6: " + skipchr6);
        System.out.println("matchbyrsid: " + matchbyrsid);
        System.out.println("mafthreshold: " + mafthreshold);
        System.out.println("hwepthreshold: " + hwepthreshold);

        this.ldthreshold = ldthreshold;
        this.prunedistance = prunedistance;
        this.pruneLDthreshold = prunethreshold;
        HashSet<String> gwasSet = loadSet(gwasset);
        HashSet<String> eqtlSet = loadSet(eqtlset);

        System.out.println("gwasSet: " + gwasset + " has " + gwasSet.size() + " snps");
        System.out.println("eqtlSet: " + eqtlset + " has " + eqtlSet.size() + " snps");

        // load eQTL info
        HashMap<String, ArrayList<eQTLResult>> eqtlsPerSNP = loadEQTLInfo(eqtlfile);

        // load GWAS info
        GWASInfo gwasInfo = loadGWASInfo2(gwasAssocFile, gwasListFile, gwasSet);

        System.out.println("Loading data from: " + vcfPrefixString);

        HashSet<String> allSNPsSet = new HashSet<>();
        allSNPsSet.addAll(gwasSet);
        allSNPsSet.addAll(eqtlSet);
        ArrayList<String> allSNPs = new ArrayList<>();
        allSNPs.addAll(allSNPsSet);
        System.out.println(allSNPs.size() + " SNPs total");
        Collections.sort(allSNPs);

        // determine whether SNPs are in the correct format
        if (allSNPs.size() == 0) {
            System.err.println("No SNPs loaded");
            System.exit(-1);
        }

        System.out.println("Filtering variants for chr:pos:rsid format.");
        ArrayList<String> filteredListOfSNPs = new ArrayList<>();
        for (String s : allSNPs) {
            String[] elems = s.split(":");
            if (elems.length < 3) {
                System.out.println("Skipping " + s + " not in the correct format. Expect: chr:pos:rsid:allele1_allele2");
            }
            try {
                int chr = Integer.parseInt(elems[0]);
                int pos = Integer.parseInt(elems[1]);
                filteredListOfSNPs.add(s);
            } catch (NumberFormatException e) {
                System.out.println("Skipping " + s + " not in the correct format. Expect: chr:pos:rsid:allele1_allele2");
            }
        }
        if (allSNPs.isEmpty()) {
            System.err.println("No snps conform to the chr:pos:rsid format.");
            System.exit(-1);
        }

        allSNPs = filteredListOfSNPs;
        System.out.println(allSNPs.size() + " SNPs total after filter");

        // preload genotypes
        HashMap<String, VCFVariant> genotypes = loadGenotypes(allSNPs, vcfPrefixString, individualSubset, skipchr6, matchbyrsid);
        writeSNPSet(eqtlSet, genotypes, outputfile + "-eQTLSNPsPresent.txt");
        writeSNPSet(gwasSet, genotypes, outputfile + "-GWASSNPsPresent.txt");


        // no need to prune the eqtl snps actually
        Pair<HashMap<String, LDCluster>, ArrayList<LDCluster>> gwasClusterInfo = prune(gwasSet, genotypes, outputfile + "-GWASSNPsPruned.txt");
        HashMap<String, LDCluster> gwasSNPToCluster = gwasClusterInfo.getLeft();
        ArrayList<LDCluster> gwasClusters = gwasClusterInfo.getRight();

        TextFile output = new TextFile(outputfile + "-LD.txt.gz", TextFile.W);
        output.writeln("TraitId\tTrait\tTraitSNP\tTraitP\tEQTLSNP\tEQTLP\tEQTLGene\tEQTLGeneSymbol\tLD(Rsq)\tGwasClusterID\tGwasClusterSize\tGwasClusterSNPs");


        ProgressBar pb = new ProgressBar(gwasClusters.size(), "Calculating LD.");
        // IntStream.range(0, allGWASIDs.size()).forEach(g -> {
        DetermineLDVCF ld = new DetermineLDVCF();
        for (int g = 0; g < gwasClusters.size(); g++) {
            LDCluster gwasCluster = gwasClusters.get(g);
            String proxyStr = gwasCluster.getProxyStr();
            Set<String> gwasProxies = gwasCluster.ldpairs.keySet();
            for (String gwasSNPId : gwasProxies) {
                VCFVariant gwasSNPObj = genotypes.get(gwasSNPId);

                if (gwasSNPObj != null) { // if there are genotypes
                    String chr = gwasSNPObj.getChr();
                    int pos = gwasSNPObj.getPos();
                    for (String eqtlSNP : eqtlSet) {
                        VCFVariant eQTLSNPObj = genotypes.get(eqtlSNP);
                        if (eQTLSNPObj != null) { // if there are genotypes
                            String chr2 = eQTLSNPObj.getChr();
                            int pos2 = eQTLSNPObj.getPos();
                            if (chr.equals(chr2) && Math.abs(pos - pos2) < prunedistance) {

                                Pair<Double, Double> ldvals = ld.getLD(gwasSNPObj, eQTLSNPObj);
                                double rsq = ldvals.getRight();
                                double dpr = ldvals.getLeft();
                                if (rsq >= ldthreshold) {
                                    HashMap<String, GWASResult> traits = gwasInfo.snpIDToResultPerTrait.get(gwasSNPId);

                                    // output all combinations between traits and eqtls for both snps
                                    for (String traitId : traits.keySet()) {
                                        GWASResult gwasresult = traits.get(traitId);
                                        ArrayList<eQTLResult> eQTLResults = eqtlsPerSNP.get(eqtlSNP);
                                        for (eQTLResult eQTLResult : eQTLResults) {
                                            String outln = traitId + "\t" + gwasInfo.traitIdToTrait.get(traitId)
                                                    + "\t" + gwasresult.snp + "\t" + gwasresult.p  // + "\t" + gwasSNPObj.getMAF() + "\t" + gwasSNPObj.getHwep()
                                                    + "\t" + eQTLSNPObj.getId() // + "\t" + eQTLSNPObj.getMAF() + "\t" + eQTLSNPObj.getHwep()
                                                    + "\t" + eQTLResult.p + "\t" + eQTLResult.gene + "\t" + eQTLResult.hugo
                                                    + "\t" + rsq // + "\t" + dpr + "\t" +
                                                    + "\t" + gwasCluster.id + "\t" + gwasCluster.size() + "\t" + proxyStr;
                                            output.writeln(outln);

                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
//            HashMap<String, GWASResult> gwasSNPsPValues = gwasSNPsPerGWAS.get(gwasID);
            pb.set(g);
        }
        // );
        pb.close();
        output.close();

    }

    class GWASInfo {
        HashMap<String, String> traitIdToTrait = new HashMap<>();
        HashMap<String, HashMap<String, GWASResult>> snpIDToResultPerTrait = new HashMap<>();
        HashMap<String, HashMap<String, GWASResult>> traitIdToResultPerSNP = new HashMap<>();
    }

    private GWASInfo loadGWASInfo2(String gwasAssocFile, String gwasListFile, HashSet<String> gwasSet) throws IOException {
        GWASInfo info = new GWASInfo();
        TextFile tf = new TextFile(gwasListFile, TextFile.R);
        String[] header = tf.readLine().split("\t");

        int idcol = -1;
        int traitcol = -1;
        int maxcol = -1;
        for (int c = 0; c < header.length; c++) {
            String v = header[c].toLowerCase();
            if (v.equals("id")) {
                idcol = c;
                if (c > maxcol) {
                    maxcol = c;
                }
            } else if (v.equals("trait")) {
                traitcol = c;
                if (c > maxcol) {
                    maxcol = c;
                }
            }
        }

        System.out.println("Trait col: " + traitcol);
        System.out.println("Id col: " + idcol);
        if (idcol == -1 || traitcol == -1) {
            System.err.println("Error parsing " + gwasListFile);
            System.exit(-1);
        }

        info.traitIdToTrait = new HashMap<>();
        String[] elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {
            if (elems.length >= maxcol) {
                String trait = elems[traitcol];
                String id = elems[idcol];
                info.traitIdToTrait.put(id, trait);
            }
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();
        System.out.println(info.traitIdToTrait.size() + " pairs of IDs and traits loaded.");

        System.out.println("Loading GWAS assoc: " + gwasAssocFile);
        tf = new TextFile(gwasAssocFile, TextFile.R);
        header = tf.readLine().split("\t");
        idcol = -1;
        int otherallelecol = -1;
        int effectallelecol = -1;
        int betacol = -1;
        int secol = -1;
        int snpcol = -1;
        int pvalcol = -1;

        maxcol = -1;
        for (int c = 0; c < header.length; c++) {
            String v = header[c].toLowerCase();
            // ID      RsID    OtherAllele     EffectAllele    EffectAlleleFreq        Beta    SE      Pvalue
            if (v.equals("id")) {
                idcol = c;
                if (c > maxcol) {
                    maxcol = c;
                }
            } else if (v.equals("otherallele")) {
                otherallelecol = c;
                if (c > maxcol) {
                    maxcol = c;
                }
            } else if (v.equals("effectallele")) {
                effectallelecol = c;
                if (c > maxcol) {
                    maxcol = c;
                }
            } else if (v.equals("rsid")) {
                snpcol = c;
                if (c > maxcol) {
                    maxcol = c;
                }
            } else if (v.equals("beta")) {
                betacol = c;
                if (c > maxcol) {
                    maxcol = c;
                }
            } else if (v.equals("se")) {
                secol = c;
                if (c > maxcol) {
                    maxcol = c;
                }
            } else if (v.equals("pvalue")) {
                pvalcol = c;
                if (c > maxcol) {
                    maxcol = c;
                }
            }
        }

        System.out.println("ID col:" + idcol);
        System.out.println("OtherAllele col:" + otherallelecol);
        System.out.println("EffectAllele col:" + effectallelecol);
        System.out.println("Beta col:" + betacol);
        System.out.println("SE col:" + secol);
        System.out.println("RSID col:" + snpcol);
        System.out.println("Pvalue col:" + pvalcol);

        if (otherallelecol == -1 || effectallelecol == -1 || betacol == -1 || secol == -1 || snpcol == -1 || pvalcol == -1) {
            System.err.println("Could not find all required columns");
            System.exit(-1);
        }
        info.traitIdToResultPerSNP = new HashMap<>();
        elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {
            if (elems.length >= maxcol) {
                // ID      RsID    OtherAllele     EffectAllele    EffectAlleleFreq        Beta    SE      Pvalue
                String snp = elems[snpcol];
                if (gwasSet.contains(snp)) {
                    String traitId = elems[idcol];
                    HashMap<String, GWASResult> snpsForTrait = info.traitIdToResultPerSNP.get(traitId);
                    if (snpsForTrait == null) {
                        snpsForTrait = new HashMap<>();
                    }
                    HashMap<String, GWASResult> traitsForSNP = info.snpIDToResultPerTrait.get(snp);
                    if (traitsForSNP == null) {
                        traitsForSNP = new HashMap<>();
                    }

                    String allele = elems[otherallelecol] + "/" + elems[effectallelecol];
                    String beta = elems[betacol];
                    String se = elems[secol];

                    double p = 1;
                    try {
                        p = Double.parseDouble(elems[pvalcol]);
                    } catch (NumberFormatException e) {

                    }

                    GWASResult gwasresult = new GWASResult();
                    gwasresult.traitId = traitId;
                    gwasresult.alleles = allele;
                    gwasresult.beta = beta;
                    gwasresult.se = se;
                    gwasresult.snp = snp;
                    gwasresult.p = p;

                    traitsForSNP.put(traitId, gwasresult);
                    snpsForTrait.put(snp, gwasresult);
                    info.traitIdToResultPerSNP.put(traitId, snpsForTrait);
                    info.snpIDToResultPerTrait.put(snp, traitsForSNP);
                }
            }
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();

        return info;
    }

    class LDCluster {
        int id;
        String startVariant;
        HashMap<String, Double> ldpairs = new HashMap<>();

        public String getProxyStr() {
            ArrayList<String> varlist = new ArrayList<>();
            varlist.addAll(ldpairs.keySet());
            Collections.sort(varlist);
            return Strings.concat(varlist, Strings.semicolon);
        }

        public int size() {
            return ldpairs.size();
        }
    }

    private Pair<HashMap<String, LDCluster>, ArrayList<LDCluster>> prune(HashSet<String> snpset, HashMap<String, VCFVariant> genotypes, String outputfile) throws IOException {
        // this set of pruned eQTLs is not actually used, but here for convenience for other downstream analyses
        System.out.println("Pruning eQTL SNPs");
        // get a list of pruned eQTL SNPs
        ArrayList<String> snpStrList = new ArrayList<>();
        snpStrList.addAll(snpset);
        Collections.sort(snpStrList);

        TextFile outPrune = new TextFile(outputfile + "-PrunedSNPsPerEQTL.txt", TextFile.W);
        outPrune.writeln("ClusterID\tStartSNP\tAliases\tRsquares");
        HashSet<String> visitedEQTLSNPs = new HashSet<>();

        HashMap<String, LDCluster> clusterMap = new HashMap<>();
        ArrayList<LDCluster> clusterList = new ArrayList<>();

        for (int i = 0; i < snpStrList.size(); i++) {
            String snp1 = snpStrList.get(i);
            VCFVariant snpobj1 = genotypes.get(snp1);
            if (snpobj1 != null) {
                if (!visitedEQTLSNPs.contains(snp1)) {
                    DetermineLDVCF ld = new DetermineLDVCF();
                    ArrayList<String> proxies = new ArrayList<>();
                    proxies.add(snp1);

                    ArrayList<String> ldvals = new ArrayList<>();
                    ldvals.add("1.0");
                    LDCluster cluster = new LDCluster();
                    cluster.id = clusterList.size();
                    clusterList.add(cluster);
                    clusterMap.put(snp1, cluster);
                    cluster.startVariant = snp1;
                    cluster.ldpairs.put(snp1, 1d);

                    visitedEQTLSNPs.add(snp1);
                    for (int j = i + 1; j < snpStrList.size(); j++) {
                        String snp2 = snpStrList.get(j);
                        VCFVariant snpobj2 = genotypes.get(snp2);
                        if (snpobj2 != null) {
                            if (!visitedEQTLSNPs.contains(snp2)) {
                                if (snpobj1.getChr().equals(snpobj2.getChr()) && Math.abs(snpobj1.getPos() - snpobj2.getPos()) < prunedistance) {
                                    // check LD
                                    Pair<Double, Double> ldinfo = ld.getLD(snpobj1, snpobj2);
                                    double rsq = ldinfo.getRight();
                                    if (rsq > pruneLDthreshold) {
                                        proxies.add(snp2);
                                        visitedEQTLSNPs.add(snp2);
                                        DecimalFormat df = new DecimalFormat("#.###", new DecimalFormatSymbols(Locale.US));
                                        ldvals.add(df.format(rsq));
                                        cluster.ldpairs.put(snp2, rsq);
                                        clusterMap.put(snp2, cluster);
                                    }
                                }
                            }
                        }
                    }


                    outPrune.writeln(cluster.id + "\t" + cluster.startVariant + "\t" + Strings.concat(proxies, Strings.semicolon) + "\t" + Strings.concat(ldvals, Strings.semicolon));

                }
            }
        }
        outPrune.close();
        return new Pair<>(clusterMap, clusterList);
    }

    private void writeSNPSet(HashSet<String> set, HashMap<String, VCFVariant> subset, String outputfile) throws IOException {
        TextFile tfq = new TextFile(outputfile, TextFile.W);
        int setctr = 0;
        for (String s : set) {
            if (subset.containsKey(s)) {
                tfq.writeln(s);
                setctr++;
            }
        }
        tfq.close();
        System.out.println(setctr + " SNPs found to overlap, written to: " + outputfile);
    }

    public void runPerGWAS2(String eqtlset, String eqtlfile, String gwasset, String gwasAssocFile, String gwasListFile, String vcfPrefixString, String individualSubset, int prunedistance, double prunethreshold, double ldthreshold, String outputfile, boolean skipchr6, boolean matchbyrsid) throws IOException {

        System.out.println("eqtlset: " + eqtlset);
        System.out.println("eqtlfile: " + eqtlfile);
        System.out.println("gwasset: " + gwasset);
        System.out.println("gwasassocfile: " + gwasAssocFile);
        System.out.println("gwaslistfile: " + gwasListFile);
        System.out.println("gwaslistfile: " + vcfPrefixString);
        System.out.println("gte: " + individualSubset);
        System.out.println("prunedistance: " + prunedistance);
        System.out.println("prunethreshold: " + prunethreshold);
        System.out.println("ldthreshold: " + ldthreshold);
        System.out.println("out: " + outputfile);
        System.out.println("skipchr6: " + skipchr6);
        System.out.println("matchbyrsid: " + matchbyrsid);
        System.out.println("mafthreshold: " + mafthreshold);
        System.out.println("hwepthreshold: " + hwepthreshold);

        HashSet<String> gwasSet = loadSet(gwasset);
        HashSet<String> eqtlSet = loadSet(eqtlset);

        System.out.println("gwasSet: " + gwasset + " has " + gwasSet.size() + " snps");
        System.out.println("eqtlSet: " + eqtlset + " has " + eqtlSet.size() + " snps");

        // load eQTL info
        HashMap<String, ArrayList<eQTLResult>> eqtlsPerSNP = loadEQTLInfo(eqtlfile);

        // load GWAS info
        Pair<HashMap<String, String>, HashMap<String, HashMap<String, GWASResult>>> gwasinfo = loadGWASInfo(gwasAssocFile, gwasListFile, gwasSet);
        HashMap<String, String> gwasIdToTrait = gwasinfo.getLeft();
        HashMap<String, HashMap<String, GWASResult>> gwasSNPsPerGWAS = gwasinfo.getRight();

        System.out.println("Loading data from: " + vcfPrefixString);

        HashSet<String> allSNPsSet = new HashSet<>();
        allSNPsSet.addAll(gwasSet);
        allSNPsSet.addAll(eqtlSet);
        ArrayList<String> allSNPs = new ArrayList<>();
        allSNPs.addAll(allSNPsSet);
        System.out.println(allSNPs.size() + " SNPs total");
        Collections.sort(allSNPs);

        // determine whether SNPs are in the correct format
        if (allSNPs.size() == 0) {
            System.err.println("No SNPs loaded");
            System.exit(-1);
        }

        System.out.println("Filtering variants for chr:pos:rsid format.");
        ArrayList<String> filteredListOfSNPs = new ArrayList<>();
        for (String s : allSNPs) {
            String[] elems = s.split(":");
            if (elems.length < 3) {
                System.out.println("Skipping " + s + " not in the correct format. Expect: chr:pos:rsid:allele1_allele2");
            }
            try {
                int chr = Integer.parseInt(elems[0]);
                int pos = Integer.parseInt(elems[1]);
                filteredListOfSNPs.add(s);
            } catch (NumberFormatException e) {
                System.out.println("Skipping " + s + " not in the correct format. Expect: chr:pos:rsid:allele1_allele2");
            }
        }
        if (allSNPs.isEmpty()) {
            System.err.println("No snps conform to the chr:pos:rsid format.");
            System.exit(-1);
        }

        allSNPs = filteredListOfSNPs;
        System.out.println(allSNPs.size() + " SNPs total after filter");

        // preload genotypes
        HashMap<String, VCFVariant> genotypes = loadGenotypes(allSNPs, vcfPrefixString, individualSubset, skipchr6, matchbyrsid);

        // write an output file listing available eQTL SNPs
        writeSNPSet(eqtlSet, genotypes, outputfile + "-eQTLSNPsPresent.txt");
        writeSNPSet(gwasSet, genotypes, outputfile + "-GWASSNPsPresent.txt");

        // get a list of pruned eQTL SNPs
        ArrayList<String> alleqtl = new ArrayList<>();
        alleqtl.addAll(eqtlSet);

        // this set of pruned eQTLs is not actually used, but here for convenience for other downstream analyses
        System.out.println("Pruning eQTL SNPs");
        TextFile outPruneEQTL = new TextFile(outputfile + "-PrunedSNPsPerEQTL.txt", TextFile.W);
        outPruneEQTL.writeln("SNP\tAliases");
        HashSet<String> visitedEQTLSNPs = new HashSet<>();
        for (int i = 0; i < alleqtl.size(); i++) {
            String snp1 = alleqtl.get(i);
            VCFVariant snpobj1 = genotypes.get(snp1);
            if (snpobj1 != null) {
                if (!visitedEQTLSNPs.contains(snp1)) {

                    DetermineLDVCF ld = new DetermineLDVCF();
                    ArrayList<String> aliases = new ArrayList<>();
                    visitedEQTLSNPs.add(snp1);
                    for (int j = i + 1; j < alleqtl.size(); j++) {
                        String snp2 = alleqtl.get(j);
                        VCFVariant snpobj2 = genotypes.get(snp2);
                        if (snpobj2 != null) {
                            if (!visitedEQTLSNPs.contains(snp2)) {
                                if (snpobj1.getChr().equals(snpobj2.getChr()) && Math.abs(snpobj1.getPos() - snpobj2.getPos()) < prunedistance) {
                                    // check LD
                                    Pair<Double, Double> ldinfo = ld.getLD(snpobj1, snpobj2);
                                    double rsq = ldinfo.getRight();
                                    if (rsq > ldthreshold) {
                                        aliases.add(snp2);
                                        visitedEQTLSNPs.add(snp2);
                                    }
                                }
                            }
                        }
                    }
                    if (aliases.isEmpty()) {
                        outPruneEQTL.writeln(snp1 + "\t" + "-");
                    } else {
                        outPruneEQTL.writeln(snp1 + "\t" + Strings.concat(aliases, Strings.semicolon));
                    }
                }
            }
        }
        outPruneEQTL.close();

        // for each GWAS, determine loaded SNPs, and their associated p-values
        ArrayList<String> allGWASIDs = new ArrayList<>();
        allGWASIDs.addAll(gwasSNPsPerGWAS.keySet());
        System.out.println(allGWASIDs.size() + " GWASes to test");

        TextFile outAll = new TextFile(outputfile + "-all.txt.gz", TextFile.W);
        String header = "GWASId\tTrait\tTraitP\tTraitSNP1\tMAF\tHWEP\tSNP2\tMAF\tHWEP\tRSQ\tDPR\tENSG\tGene";
        outAll.writeln(header);

        TextFile outSum = new TextFile(outputfile + "-summary.txt", TextFile.W);
        String header2 = "GWASID\tTrait\tNrGWASSNPsTotal\tNrPrunedGWASSNPsLinked\tNrPrunedGWASSNPs\tPercGWASSNPsLinked\tnrEQTLSNPsLinked\tTotalEQTLSNPs\tPercEQTLSNPsLinked";
        outSum.writeln(header2);

        TextFile outPrune = new TextFile(outputfile + "-PrunedSNPsPerTrait.txt", TextFile.W);
        String pruneheader = "GWASID\tTrait\tIndexVariant\tIndexVariantP\tIndexVariantAlleles\tIndexVariantEffect\tLinkedEQTLSNP\tLD(rsq)\tLinkedEQTLGenes\tLinkedEQTLHgncIDs\tLinkedEQTLAlleles\tLinkedEQTLZScores\tLinkedEQTLP\tGWASClusterSize\tSNPsInGWASCluster";
        outPrune.writeln(pruneheader);

        TextFile clusteroutput = new TextFile(outputfile + "-GWASClusters.txt", TextFile.W);
        clusteroutput.writeln("Index\tNrSNPs\tSNPsInCluster");
        // count eQTL SNPs with genotypes
        int nrEQTLSNPs = 0;
        for (String snp : eqtlSet) {
            if (genotypes.containsKey(snp)) {
                nrEQTLSNPs++;
            }
        }
        int finalNrEQTLSNPs = nrEQTLSNPs;


        ProgressBar pb = new ProgressBar(allGWASIDs.size(), "Calculating LD.");
        IntStream.range(0, allGWASIDs.size()).forEach(g -> {
            String gwasID = allGWASIDs.get(g);
            HashMap<String, GWASResult> gwasSNPsPValues = gwasSNPsPerGWAS.get(gwasID);

            // find gwasSNPsPValues that are within 1mb of each other
            ArrayList<String> allGWASSNPs = new ArrayList<>();
            allGWASSNPs.addAll(gwasSNPsPValues.keySet());
            Collections.sort(allGWASSNPs);

            // determine linked GWAS variants linked to eQTL variants
            HashMap<String, HashSet<String>> gwasSNPToLinkedEQTLSNPs = new HashMap<>();
            HashMap<String, HashMap<String, Double>> gwasSNPToLinkedEQTLSNPsLD = new HashMap<>();
            // for each GWAS SNP
            // first determine LD between all available GWAS and eQTL SNPs, without pruning them
            IntStream.range(0, allGWASSNPs.size()).forEach(i -> {
                DetermineLDVCF ld = new DetermineLDVCF();
                String gwasSNPId = allGWASSNPs.get(i);
                VCFVariant gwasSNPObj = genotypes.get(gwasSNPId); // get genotypes

                HashSet<String> linkedSNPs = new HashSet<>();
                if (gwasSNPObj != null) { // if there are genotypes
                    String chr = gwasSNPObj.getChr();
                    int pos = gwasSNPObj.getPos();
                    for (String eqtlSNP : eqtlSet) { // compare to each eQTL SNP
                        VCFVariant eQTLSNPObj = genotypes.get(eqtlSNP);
                        if (eQTLSNPObj != null) { // if there are genotypes
                            String chr2 = eQTLSNPObj.getChr();
                            int pos2 = eQTLSNPObj.getPos();
                            if (chr.equals(chr2) && Math.abs(pos - pos2) < prunedistance) {
                                Pair<Double, Double> ldvals = ld.getLD(gwasSNPObj, eQTLSNPObj);
                                double rsq = ldvals.getRight();
                                double dpr = ldvals.getLeft();

                                synchronized (gwasSNPToLinkedEQTLSNPsLD) { // this is some logic for multithreading. isn't used anymore
                                    HashMap<String, Double> ldsnps = gwasSNPToLinkedEQTLSNPsLD.get(gwasSNPId);
                                    if (ldsnps == null) {
                                        ldsnps = new HashMap<>();
                                    }
                                    ldsnps.put(eqtlSNP, rsq);
                                    gwasSNPToLinkedEQTLSNPsLD.put(gwasSNPId, ldsnps);
                                }

                                if (!Double.isNaN(rsq)) {
                                    if (rsq > ldthreshold) {
                                        linkedSNPs.add(eqtlSNP);
                                    }
                                    ArrayList<eQTLResult> genedata = eqtlsPerSNP.get(eqtlSNP);
                                    ArrayList<String> genes = new ArrayList<>();
                                    ArrayList<String> hugos = new ArrayList<>();
                                    String ensg = "-";
                                    String hugo = "-";
                                    if (genedata != null) {
                                        for (eQTLResult p : genedata) {
                                            genes.add(p.gene);
                                            hugos.add(p.hugo);
                                        }
                                        ensg = Strings.concat(genes, Strings.semicolon);
                                        hugo = Strings.concat(hugos, Strings.semicolon);
                                    }

                                    try {
                                        outAll.writelnsynced(gwasID + "\t" + gwasIdToTrait.get(gwasID) + "\t" + gwasSNPsPValues.get(gwasSNPObj.getId()).p + "\t" + gwasSNPObj.getId() + "\t" + gwasSNPObj.getMAF() + "\t" + gwasSNPObj.getHwep() + "\t" + eqtlSNP + "\t" + eQTLSNPObj.getMAF() + "\t" + eQTLSNPObj.getHwep() + "\t" + rsq + "\t" + dpr + "\t" + ensg + "\t" + hugo);
                                    } catch (IOException e) {
                                        e.printStackTrace();
                                    }

                                }
                            }
                        }
                    }
                }
                synchronized (gwasSNPToLinkedEQTLSNPs) {
                    HashSet<String> linked = gwasSNPToLinkedEQTLSNPs.get(gwasSNPId);
                    if (linked == null) {
                        linked = new HashSet<>();
                    }
                    linked.addAll(linkedSNPs);
                    gwasSNPToLinkedEQTLSNPs.put(gwasSNPId, linked);
                }
            });

            // now prune / clump GWAS set
            HashSet<String> visited = new HashSet<String>();
            int nrGWASClustersLinked = 0;
            int nrGWASClusters = 0;
            HashSet<String> eqtlSNPsLinkedToGWASSNPClusters = new HashSet<String>();

            // iterate GWAS SNPs
            for (int s = 0; s < allGWASSNPs.size(); s++) {
                String gwasSNP1 = allGWASSNPs.get(s);
                VCFVariant gwasSNP1Obj = genotypes.get(gwasSNP1);
                String indexVariant = gwasSNP1;

                if (gwasSNP1Obj != null && !visited.contains(gwasSNP1)) { // if we haven't seen this variant before and there are genotypes
                    String chr = gwasSNP1Obj.getChr();
                    Integer pos = gwasSNP1Obj.getPos(); // Integer.parseInt(elems[1]);

                    // check whether there are any SNPs linked to the current one
                    // determine which gwasSNPsPValues are within prunedistance of each other that are in LD
                    // if so, gather all linked eQTLs
                    ArrayList<Pair<String, Double>> gwasCluster = new ArrayList<>(); // start a new 'cluster' of GWAS SNPs
                    HashSet<String> eqtlSNPsLinkedToCluster = new HashSet<>(); // determine which eQTLs should be linked to this cluster

                    // add SNPs linked to this variant
                    HashSet<String> linked = gwasSNPToLinkedEQTLSNPs.get(gwasSNP1);
                    if (linked != null) {
                        eqtlSNPsLinkedToCluster.addAll(linked);
                    }

                    GWASResult gwasResult = gwasSNPsPValues.get(gwasSNP1);
                    Pair<String, Double> pair = new Pair<String, Double>(gwasSNP1, gwasResult.p, Pair.SORTBY.RIGHT); // keep this GWAS variants' GWAS sumstats somewhere
                    gwasCluster.add(pair);
                    visited.add(gwasSNP1);
                    // iterate all other GWAS SNPs
                    for (int s2 = s + 1; s2 < allGWASSNPs.size(); s2++) {
                        String gwasSNP2 = allGWASSNPs.get(s2);
                        VCFVariant gwasSNP2Obj = genotypes.get(gwasSNP2);
                        if (gwasSNP2Obj != null && !visited.contains(gwasSNP2)) { // if this second variant has genotypes and hasn't been seen before
                            String chr2 = gwasSNP2Obj.getChr();
                            Integer pos2 = gwasSNP2Obj.getPos(); // Integer.parseInt(elems2[1]);
                            if (chr.equals(chr2)) { // if SNP1 and 2 are on the same chromosome
                                int distance = Math.abs(pos - pos2);
                                if (distance < prunedistance) { // and the distance between them is < prunedistance
                                    // determine LD
                                    DetermineLDVCF ld = new DetermineLDVCF();
                                    VCFVariant snpjobj = genotypes.get(gwasSNP2);
                                    Pair<Double, Double> ldvals = ld.getLD(gwasSNP1Obj, snpjobj);
                                    double rsq = ldvals.getRight();
                                    if (rsq > prunethreshold) { // add SNP2 to cluster if LD is above the prune threshold
                                        GWASResult gwasResult2 = gwasSNPsPValues.get(gwasSNP2);
                                        Pair<String, Double> pair2 = new Pair<String, Double>(gwasSNP2, gwasResult2.p, Pair.SORTBY.RIGHT);
                                        gwasCluster.add(pair2);
                                        visited.add(gwasSNP2);

                                        // get all linked eQTLs
                                        linked = gwasSNPToLinkedEQTLSNPs.get(gwasSNP2);

                                        // add all linked eQTLs to 'gwasCluster' set
                                        if (linked != null) {
                                            eqtlSNPsLinkedToCluster.addAll(linked);
                                        }
                                    }
                                }
                            }
                        }
                    }

                    // prune gwasCluster
                    ArrayList<VCFVariant> remainingSNPs = new ArrayList<>();
                    if (gwasCluster.size() == 1) { // this GWAS SNPs is not in LD with other GWAS variants
                        String snpi = gwasCluster.get(0).getLeft();
                        VCFVariant snpiobj = genotypes.get(snpi);
                        remainingSNPs.add(snpiobj);
                    } else if (gwasCluster.size() > 1) {
                        // determine index variant in this cluster using summary statistics
                        // sort gwasCluster on P
                        Collections.sort(gwasCluster); // sort by GWAS P-value
                        System.out.println();
                        System.out.println("Cluster of " + gwasCluster.size() + " GWAS variants found for " + gwasSNP1);
                        int min = Math.min(5, gwasCluster.size());
                        for (int i = 0; i < gwasCluster.size(); i++) {
                            System.out.println(gwasCluster.get(i).getLeft() + "\t" + gwasCluster.get(i).getRight());
                        }
                        System.out.println();
                        indexVariant = gwasCluster.get(0).getLeft(); // take the top variant as index variant for this cluster
//                        System.out.println(gwasCluster.size() + " SNPs in gwasCluster " + clusterctr.get() + ", " + remainingSNPs.size() + " remain after pruning for trait: " + gwasID + "\t" + gwasIdToTrait.get(gwasID));
//                        waitForEnter("Check this out!");
                    }

                    // eqtlSNPsLinkedToCluster has all linked SNPs for this gwasCluster.
                    if (!eqtlSNPsLinkedToCluster.isEmpty()) {
                        nrGWASClustersLinked++;
                        eqtlSNPsLinkedToGWASSNPClusters.addAll(eqtlSNPsLinkedToCluster);
                    }
                    nrGWASClusters++;

                    // concatenate GWAS SNP IDs
                    String clusterStr = "";
                    for (Pair<String, Double> clusterpair : gwasCluster) {
                        String snp = clusterpair.getLeft();
                        if (!snp.equals(gwasSNP1)) {
                            if (clusterStr.length() == 0) {
                                clusterStr += snp;
                            } else {
                                clusterStr += ";" + snp;
                            }
                        }
                    }

                    // write data to disk
                    for (String eqtl : eqtlSNPsLinkedToCluster) {
                        ArrayList<eQTLResult> eqtlinfo = eqtlsPerSNP.get(eqtl);
                        if (eqtlinfo != null) {
                            for (eQTLResult r : eqtlinfo) {
                                // write index variant, clustered SNPs, and linked eQTL SNPs
                                // get list of genes as well?
                                // 	GWASID	Trait	IndexVariant	IndexVariantP	IndexVariantAlleles	IndexVariantEffect	nrEQTLSNPsLinked
                                // 	LinkedEQTLSNPs	LinkedEQTLGenes	LinkedEQTLHgncIDs	LinkedEQTLAlleles	LinkedEQTLZScores	GWASClusterSize	SNPsInGWASCluster
                                HashMap<String, Double> ldsnps = gwasSNPToLinkedEQTLSNPsLD.get(indexVariant);
                                String actualIndexVariant = indexVariant;
                                if (ldsnps == null) {
                                    System.err.println("Error: no LD snps for " + indexVariant);
                                }
                                Double rsq = ldsnps.get(r.snp);

                                if (rsq == null) { // this can happen if the index variant of the cluster is > 1mb away from the eqtl variant
                                    if (gwasCluster.size() > 1) {
                                        System.out.println("No matching LD found for " + r.snp + " and gwas cluster index: " + indexVariant + " trying next best variant(s)");
                                        int currentindex = 1; // skip 0, because that one is clearly not working
                                        while (rsq == null && currentindex < gwasCluster.size()) {
                                            String tmpIndexVariant = gwasCluster.get(currentindex).getLeft();
                                            ldsnps = gwasSNPToLinkedEQTLSNPsLD.get(tmpIndexVariant);
                                            Double tmprsq = ldsnps.get(r.snp);
                                            if (tmprsq == null) {
                                                currentindex++;
                                                System.out.println(currentindex + "\t" + tmpIndexVariant + "\t" + gwasCluster.get(currentindex).getRight() + " has no LD info for " + r.snp);
                                            } else {
                                                actualIndexVariant = tmpIndexVariant;
                                                System.out.println(currentindex + "\t" + tmpIndexVariant + "\t" + gwasCluster.get(currentindex).getRight() + " has LD info for " + r.snp + " - " + tmprsq);
                                                System.out.println("Using this variant as the index for this cluster.");
                                                rsq = tmprsq;
                                            }
                                        }
                                    }

//                                    System.err.println("Error: no LD info for " + r.snp + " in eqtl: " + eqtl);
//                                    System.out.println("Index: " + indexVariant);
//                                    for (String key : ldsnps.keySet()) {
//                                        System.out.println(key + "\t" + ldsnps.get(key));
//                                    }
//                                    System.out.println("");
//                                    System.out.println("GWAS cluster: ");
//                                    for (Pair<String, Double> clusterpair : gwasCluster) {
//                                        String snp = clusterpair.getLeft();
//                                        System.out.println(clusterpair.getLeft() + "\t" + clusterpair.getRight());
//                                    }
//
//
//                                    System.exit(-1);
                                }

                                if (rsq != null) {
                                    // check direction?
                                    String[] allelesGWAS = gwasSNPsPValues.get(actualIndexVariant).alleles.split("/");
                                    String assessedGWAS = "";
                                    String gwasAlleleStr = gwasSNPsPValues.get(actualIndexVariant).alleles;
//                                direction dir = direction.UNKNOWN;
//                                if (allelesGWAS.length > 1) {
//                                    assessedGWAS = allelesGWAS[1];
//                                    Boolean alleleflip = BaseAnnot.flipalleles(gwasSNPsPValues.get(gwasSNP1).alleles, assessedGWAS, r.alleles, r.assessed);
//                                    if (alleleflip != null) {
//                                        try {
//                                            double beta = Double.parseDouble(gwasSNPsPValues.get(gwasSNP1).beta);
//                                            double z = Double.parseDouble(r.z);
//                                            if (alleleflip) {
//                                                z *= -1;
//                                            }
//                                            if ((beta >= 0 && z >= 0) || (beta < 0 && z <= 0)) {
//                                                dir = direction.SAME;
//                                            } else {
//                                                dir = direction.OPPOSITE;
//                                            }
//                                        } catch (NumberFormatException e) {
//                                        }
//                                    }
//                                }

                                    String outln = gwasID + "\t" + gwasIdToTrait.get(gwasID) + "\t" + actualIndexVariant + "\t" + gwasSNPsPValues.get(actualIndexVariant).p + "\t" + gwasAlleleStr + "\t" + gwasSNPsPValues.get(actualIndexVariant).beta + "\t" + r.snp + "\t" + rsq + "\t" + r.gene + "\t" + r.hugo + "\t" + r.alleles + "\t" + r.assessed + "\t" + r.z + "\t" + r.p + "\t" + gwasCluster.size() + "\t" + clusterStr;
                                    try {
                                        outPrune.writelnsynced(outln);
                                    } catch (IOException e) {
                                        e.printStackTrace();
                                    }
                                }
                            }
                        }
                    }


                } // end: if snpobj!=null and notvisited(snp)
            } // end loop: iterate GWAS SNPs


            // write a summary
            double perc1 = (double) nrGWASClustersLinked / nrGWASClusters;
            if (Double.isNaN(perc1)) {
                perc1 = 0;
            }
            double perc2 = (double) eqtlSNPsLinkedToGWASSNPClusters.size() / finalNrEQTLSNPs;
            if (Double.isNaN(perc2)) {
                perc2 = 0;
            }

            String outln = gwasID + "\t" + gwasIdToTrait.get(gwasID) + "\t" + gwasSNPsPValues.size() + "\t" + nrGWASClustersLinked + "\t" + nrGWASClusters + "\t" + perc1 + "\t" + eqtlSNPsLinkedToGWASSNPClusters.size() + "\t" + finalNrEQTLSNPs + "\t" + perc2;

            try {
                outSum.writelnsynced(outln);
            } catch (IOException e) {
                e.printStackTrace();
            }
            pb.iterateSynched();
        });
        pb.close();
        outAll.close();
        outSum.close();
        outPrune.close();
        clusteroutput.close();

    }

    public HashSet<String> loadSet(String setfile) throws IOException {
        System.out.println("Reading set: " + setfile);
        TextFile tf = new TextFile(setfile, TextFile.R);
        String ln = tf.readLine();
        HashSet<String> set = new HashSet<>();
        while (ln != null) {
            set.add(ln);
            ln = tf.readLine();
        }
        tf.close();
        System.out.println(set.size() + " items in " + setfile);
        return set;
    }

    // expects a beta-qtl style QTL output file
    private HashMap<String, ArrayList<eQTLResult>> loadEQTLInfo(String eqtlfile) throws IOException {
        System.out.println("Parsing: " + eqtlfile);
        TextFile tf = new TextFile(eqtlfile, TextFile.R);
        HashMap<String, ArrayList<eQTLResult>> output = new HashMap<>();
        String[] header = tf.readLine().split("\t");

        int snpcol = -1;
        int genecol = -1;
        int zcol = -1;
        int pcol = -1;
        int allelescol = -1;
        int assessedcol = -1;
        int hugocol = -1;

        int maxc = -1;
        for (int c = 0; c < header.length; c++) {
            String v = header[c].toLowerCase();
            if (v.equals("SNP".toLowerCase())) {
                snpcol = c;
                if (c > maxc) {
                    maxc = c;
                }
            } else if (v.equals("Gene".toLowerCase())) {
                genecol = c;
                if (c > maxc) {
                    maxc = c;
                }
            } else if (v.equals("MetaPZ".toLowerCase())) {
                zcol = c;
                if (c > maxc) {
                    maxc = c;
                }
            } else if (v.equals("MetaP".toLowerCase())) {
                pcol = c;
                if (c > maxc) {
                    maxc = c;
                }
            } else if (v.equals("SNPAlleles".toLowerCase())) {
                allelescol = c;
                if (c > maxc) {
                    maxc = c;
                }
            } else if (v.equals("SNPEffectAllele".toLowerCase())) {
                assessedcol = c;
                if (c > maxc) {
                    maxc = c;
                }
            } else if (v.equals("GeneSymbol".toLowerCase())) {
                hugocol = c;
                if (c > maxc) {
                    maxc = c;
                }
            }
        }

        System.out.println("SNP col: " + snpcol);
        System.out.println("Gene col: " + genecol);
        System.out.println("MetaP col: " + pcol);
        System.out.println("MetaPZ col: " + zcol);
        System.out.println("SNPAlleles col: " + allelescol);
        System.out.println("SNPEffectAllele col: " + assessedcol);
        System.out.println("GeneSymbol col: " + hugocol);

        if (snpcol == -1 || genecol == -1 || pcol == -1 || zcol == -1 || allelescol == -1 || assessedcol == -1 || hugocol == -1) {
            System.err.println("Error: did not find all expected columns in qtlfile: " + eqtlfile);
            System.exit(-1);
        }

        String[] elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {
            if (elems.length >= maxc) {
                String snp = elems[snpcol];
                String gene = elems[genecol];
                String z = elems[zcol];
                String alleles = elems[allelescol];
                String assessed = elems[assessedcol];
                String hugo = elems[hugocol];

                eQTLResult r = new eQTLResult();
                r.p = elems[pcol];
                r.alleles = alleles;
                r.assessed = assessed;
                r.z = z;
                r.snp = snp;
                r.gene = gene;
                r.hugo = hugo;
                ArrayList<eQTLResult> pairs = output.get(snp);
                if (pairs == null) {
                    pairs = new ArrayList<>();
                }
                pairs.add(r);
                output.put(snp, pairs);
            }
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();
        return output;

    }

    class eQTLResult {
        String p;
        String snp;
        String gene;
        String hugo;
        String z;
        String alleles;
        String assessed;
    }

    // read gwas info
    // requires a list of gwas associations (gwasAssocFile) with ID      RsID    OtherAllele     EffectAllele    EffectAlleleFreq        Beta    SE      Pvalue columns, tab separated
    // requires a list of gwases with ID and trait columns, IDs matching the associations in the gwas associations file
    private Pair<HashMap<String, String>, HashMap<String, HashMap<String, GWASResult>>> loadGWASInfo(String gwasAssocFile, String gwasListFile, HashSet<String> gwasSet) throws IOException {
        System.out.println("Loading GWAS list: " + gwasListFile);
        TextFile tf = new TextFile(gwasListFile, TextFile.R);
        String[] header = tf.readLine().split("\t");

        int idcol = -1;
        int traitcol = -1;
        int maxcol = -1;
        for (int c = 0; c < header.length; c++) {
            String v = header[c].toLowerCase();
            if (v.equals("id")) {
                idcol = c;
                if (c > maxcol) {
                    maxcol = c;
                }
            } else if (v.equals("trait")) {
                traitcol = c;
                if (c > maxcol) {
                    maxcol = c;
                }
            }
        }

        System.out.println("Trait col: " + traitcol);
        System.out.println("Id col: " + idcol);
        if (idcol == -1 || traitcol == -1) {
            System.err.println("Error parsing " + gwasListFile);
            System.exit(-1);
        }

        HashMap<String, String> idToTrait = new HashMap<>();
        String[] elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {
            if (elems.length >= maxcol) {
                String trait = elems[traitcol];
                String id = elems[idcol];
                idToTrait.put(id, trait);
            }
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();
        System.out.println(idToTrait.size() + " pairs of IDs and traits loaded.");

        System.out.println("Loading GWAS assoc: " + gwasAssocFile);
        tf = new TextFile(gwasAssocFile, TextFile.R);
        header = tf.readLine().split("\t");
        idcol = -1;
        int otherallelecol = -1;
        int effectallelecol = -1;
        int betacol = -1;
        int secol = -1;
        int snpcol = -1;
        int pvalcol = -1;

        maxcol = -1;
        for (int c = 0; c < header.length; c++) {
            String v = header[c].toLowerCase();
            // ID      RsID    OtherAllele     EffectAllele    EffectAlleleFreq        Beta    SE      Pvalue
            if (v.equals("id")) {
                idcol = c;
                if (c > maxcol) {
                    maxcol = c;
                }
            } else if (v.equals("otherallele")) {
                otherallelecol = c;
                if (c > maxcol) {
                    maxcol = c;
                }
            } else if (v.equals("effectallele")) {
                effectallelecol = c;
                if (c > maxcol) {
                    maxcol = c;
                }
            } else if (v.equals("rsid")) {
                snpcol = c;
                if (c > maxcol) {
                    maxcol = c;
                }
            } else if (v.equals("beta")) {
                betacol = c;
                if (c > maxcol) {
                    maxcol = c;
                }
            } else if (v.equals("se")) {
                secol = c;
                if (c > maxcol) {
                    maxcol = c;
                }
            } else if (v.equals("pvalue")) {
                pvalcol = c;
                if (c > maxcol) {
                    maxcol = c;
                }
            }
        }

        System.out.println("ID col:" + idcol);
        System.out.println("OtherAllele col:" + otherallelecol);
        System.out.println("EffectAllele col:" + effectallelecol);
        System.out.println("Beta col:" + betacol);
        System.out.println("SE col:" + secol);
        System.out.println("RSID col:" + snpcol);
        System.out.println("Pvalue col:" + pvalcol);

        if (otherallelecol == -1 || effectallelecol == -1 || betacol == -1 || secol == -1 || snpcol == -1 || pvalcol == -1) {
            System.err.println("Could not find all required columns");
            System.exit(-1);
        }
        HashMap<String, HashMap<String, GWASResult>> gwasSNPsPerTrait = new HashMap<>();

        elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {
            if (elems.length >= maxcol) {
                // ID      RsID    OtherAllele     EffectAllele    EffectAlleleFreq        Beta    SE      Pvalue
                String snp = elems[snpcol];
                if (gwasSet.contains(snp)) {
                    String id = elems[idcol];
                    HashMap<String, GWASResult> snps = gwasSNPsPerTrait.get(id);
                    if (snps == null) {
                        snps = new HashMap<>();
                    }

                    String allele = elems[otherallelecol] + "/" + elems[effectallelecol];
                    String beta = elems[betacol];
                    String se = elems[secol];

                    double p = 1;
                    try {
                        p = Double.parseDouble(elems[pvalcol]);
                    } catch (NumberFormatException e) {

                    }

                    GWASResult r = new GWASResult();
                    r.alleles = allele;
                    r.beta = beta;
                    r.se = se;
                    r.snp = snp;
                    r.p = p;

                    snps.put(snp, r);
                    gwasSNPsPerTrait.put(id, snps);
                }
            }
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();
        return new Pair<>(idToTrait, gwasSNPsPerTrait);
    }

    class GWASResult {
        public String traitId;
        double p;
        String snp;
        String alleles;
        String beta;
        String se;
    }

    // this assumes the variant ids to be in the format: chr:pos:rsid:allele1_allele2
    public HashMap<String, VCFVariant> loadGenotypes(ArrayList<String> allSNPs, String vcfPrefix, String sampleLimit, boolean skipchr6, boolean matchByRsId) throws IOException {
        HashMap<String, VCFVariant> output = new HashMap<>();

        // sort variants
        Collections.sort(allSNPs);

        VCFTabix currentVCF = null;
        String currentChr = null;
        boolean[] currentSampleLimit = null;
        int notfound = 0;
        int chr6 = 0;
        for (int v = 0; v < allSNPs.size(); v++) {
            String variantid = allSNPs.get(v);
            String[] variantidelems = variantid.split(":");
            String variantChr = variantidelems[0];
            SNPObj obj = new SNPObj();
            String variantName = variantid;
            if (matchByRsId) {
                variantName = variantidelems[2];
            }
            obj.name = variantName;
            obj.chr = Chromosome.parseChr(variantChr).getNumber();
            obj.pos = Integer.parseInt(variantidelems[1]);

            if (!skipchr6 || obj.chr != 6) { // if we don't skip chr6, or chromosome is not 6
                // initialize tabix filehandle
                if (currentVCF == null || !currentChr.equals(variantChr)) {
                    String actualVCF = vcfPrefix.replaceAll("CHR", "" + variantChr);
                    System.out.println("Initializing tabix: " + actualVCF);
                    currentVCF = new VCFTabix(actualVCF);
                    if (sampleLimit != null) {
                        currentSampleLimit = currentVCF.getSampleFilter(sampleLimit);
                    }
                    currentChr = variantidelems[0];
                }

                VCFVariant variant = getVariant(currentVCF, currentSampleLimit, obj);
                if (variant != null && variant.getMAF() > mafthreshold && variant.getHwep() > hwepthreshold) {
//                    if () {
                    output.put(variantid, variant);
//                    } else {
////                    System.out.println(variantid + " not found in " + currentVCF);
//                        notfound++;
//                    }
                } else {
                    notfound++;
                }
            } else {
//                System.out.println("Skipping variant: " + variantid + " because on chr6");
                chr6++;
            }
            if (output.size() % 1000 == 0) {
                System.out.println(output.size() + " variants loaded sofar. " + chr6 + " skipped on chr6 and " + notfound + " not found.");
            }
        }
        System.out.println(output.size() + " variants loaded total. " + chr6 + " skipped on chr6 and " + notfound + " not found.");
        return output;
    }

    class SNPObj {
        String name;
        int chr;
        int pos;
    }

    enum direction {
        SAME, OPPOSITE, UNKNOWN
    }

    private VCFVariant getVariant(VCFTabix tabix, boolean[] samplefilter, SNPObj snpobj) throws IOException {
        Feature cisRegion = new Feature(Chromosome.parseChr("" + snpobj.chr), snpobj.pos - 1, snpobj.pos + 1); // you never know where to start looking in a VCF...

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
