package eqtlmappingpipeline.util;

import eqtlmappingpipeline.metaqtl3.containers.Settings;
import org.apache.commons.configuration.ConfigurationException;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.io.bin.BinaryFile;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.*;
import umcg.genetica.io.trityper.util.BaseAnnot;
import umcg.genetica.io.trityper.util.ChrAnnotation;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.*;
import java.util.stream.IntStream;

public class TriTyperToDosageMatrix {

    public void runTriTyper(String settingsfile, String settingstexttoreplace, String settingstexttoreplacewith, boolean sortbyId) throws IOException, ConfigurationException {
        System.out.println("TriTyper to dosagematrix converter.");
        Settings settings = new Settings();
        settings.settingsTextToReplace = settingstexttoreplace;
        settings.settingsTextReplaceWith = settingstexttoreplacewith;
        settings.load(settingsfile);

        // initialize datasets
        ArrayList<TriTyperGeneticalGenomicsDatasetSettings> datasetsettings = settings.datasetSettings;
        TriTyperGenotypeData[] m_gg = new TriTyperGenotypeData[datasetsettings.size()];
        SNPLoader[] loaders = new SNPLoader[datasetsettings.size()];
        System.out.println("Initializing " + datasetsettings.size() + " datasets in parallel.");
        IntStream.range(0, m_gg.length).parallel().forEach(i -> {
            m_gg[i] = new TriTyperGenotypeData();
            try {
                m_gg[i].load(datasetsettings.get(i).genotypeLocation, datasetsettings.get(i).snpmapFileLocation, datasetsettings.get(i).snpFileLocation);
                loaders[i] = m_gg[i].createSNPLoader(1);
            } catch (IOException e) {
                e.printStackTrace();
            }
        });

        System.out.println();
        System.out.println("Selecting samples per dataset");
        ArrayList<HashSet<String>> allowedSamplesPerDs = new ArrayList<HashSet<String>>();
        for (int d = 0; d < m_gg.length; d++) {
            TriTyperGeneticalGenomicsDatasetSettings dataset = datasetsettings.get(d);
            HashSet<String> allowedSamples = new HashSet<String>();
            if (dataset.genotypeToExpressionCoupling == null) {
                allowedSamples.addAll(Arrays.asList(m_gg[d].getIndividuals()));
            } else {
                TextFile tf = new TextFile(dataset.genotypeToExpressionCoupling, TextFile.R);
                String[] elems = tf.readLineElems(TextFile.tab);
                while (elems != null) {
                    allowedSamples.add(elems[0]);
                    elems = tf.readLineElems(TextFile.tab);
                }
                tf.close();
            }
            System.out.println(datasetsettings.get(d).name + " has " + allowedSamples.size() + " samples selected.");
            allowedSamplesPerDs.add(allowedSamples);
        }

        // create sample index
        System.out.println("Creating sample index..");
        HashMap<String, Integer> sampleIdMap = new HashMap<String, Integer>();
        int sampleCounter = 0;
        for (int i = 0; i < allowedSamplesPerDs.size(); i++) {
            HashSet<String> samplesInDs = allowedSamplesPerDs.get(i);
            for (String s : samplesInDs) {
                if (sampleIdMap.containsKey(s)) {
                    System.out.println("Warning: duplicate sample ID: " + s);
                    sampleIdMap.put(s, sampleCounter);
                    sampleCounter++;
                } else {
                    sampleIdMap.put(s, sampleCounter);
                    sampleCounter++;
                }
            }
        }
        System.out.println(sampleIdMap.size() + " unique sample IDs.");

        // write data to disk
        String[] header = new String[sampleIdMap.size()];
        for (Map.Entry<String, Integer> samplePair : sampleIdMap.entrySet()) {
            header[samplePair.getValue()] = samplePair.getKey();
        }
        String outdir = settings.outputReportsDir;

        TextFile indout = new TextFile(outdir + "Individuals.txt", TextFile.W);
        TextFile phenoout = new TextFile(outdir + "PhenotypeInformation.txt", TextFile.W);
        for (int i = 0; i < header.length; i++) {
            indout.writeln(header[i]);
            phenoout.writeln(header[i] + "\tunknown\tinclude\tunknown");
        }
        indout.close();
        phenoout.close();

        System.out.println("Inventorizing snps");
        // sorting SNPs
        Integer selectChr = null;
        if (settings.confineToSNPsThatMapToChromosome != null) {
            selectChr = (int) settings.confineToSNPsThatMapToChromosome;
            System.out.println("Selecting variants from chr " + selectChr);
        }

        HashSet<String> allSNPsHash = new HashSet<String>();
        for (int d = 0; d < m_gg.length; d++) {
            String[] datasetSNPs = m_gg[d].getSNPs();
            if (selectChr != null) {
                for (int s = 0; s < datasetSNPs.length; s++) {
                    int chr = m_gg[d].getChr(s);
                    if (selectChr == chr) {
                        allSNPsHash.add(datasetSNPs[s]);
                    }
                }
            } else {
                allSNPsHash.addAll(Arrays.asList(datasetSNPs));
            }
        }

        System.out.println(allSNPsHash.size() + " unique SNPs");

        System.out.println("Sorting SNPs");
        ArrayList<String> allSNPs = new ArrayList<>();
        allSNPs.addAll(allSNPsHash);

        ArrayList<String> snpsToQuery = null;
        if (settings.tsSNPsConfine != null) {
            snpsToQuery = new ArrayList<>();
            HashSet<String> confineset = settings.tsSNPsConfine;
            for (String s : confineset) {
                if (allSNPsHash.contains(s)) {
                    snpsToQuery.add(s);
                }
            }
        } else {
            snpsToQuery = allSNPs;
        }
        System.out.println(snpsToQuery.size() + " SNPs remain after filtering for SNP confinement list");


        snpsToQuery = sortSNPs(snpsToQuery, sortbyId, m_gg);
        System.out.println(snpsToQuery.size() + " SNPs to output");

        double[] dosageData = new double[sampleIdMap.size()];
        double[] genotypeData = new double[sampleIdMap.size()];

        ProgressBar pb = new ProgressBar(snpsToQuery.size(), "Querying SNPs...");
        BinaryFile bfg = new BinaryFile(outdir + "GenotypeMatrix.dat", BinaryFile.W);
        BinaryFile bfd = new BinaryFile(outdir + "ImputedDosageMatrix.dat", BinaryFile.W);
        TextFile snpo = new TextFile(outdir + "SNPs.txt.gz", TextFile.W);
        TextFile snpm = new TextFile(outdir + "SNPMappings.txt.gz", TextFile.W);
        TextFile log = new TextFile(outdir + "MergeLog.txt.gz", TextFile.W);
        int written = 0;
        for (int snpctr = 0; snpctr < snpsToQuery.size(); snpctr++) {
            String snp = snpsToQuery.get(snpctr);
            SNP[] snpObjs = new SNP[m_gg.length];
            // load genotype data per dataset
            String refSNPAlleles = null;
            String refSNPRefAllele = null;
            String refSNPAlleleAssessed = null;
            SNP refSNP = null;
            Boolean[] flip = new Boolean[m_gg.length];
            for (int d = 0; d < m_gg.length; d++) {
                Integer id = m_gg[d].getSnpToSNPId().get(snp);
                if (id >= 0) {
                    SNP snpobj = m_gg[d].getSNPObject(id);
                    loaders[d].loadGenotypes(snpobj);
//                    if (snpobj.getMAF() > settings.snpQCMAFThreshold && snpobj.getHWEP() > settings.snpQCHWEThreshold && snpobj.getCR() > settings.snpQCCallRateThreshold) {
                    snpObjs[d] = snpobj;
                    if (loaders[d].hasDosageInformation()) {
                        loaders[d].loadDosage(snpObjs[d]);
                    }
                    if (refSNPAlleles == null) {
                        refSNP = snpobj;
                        refSNPAlleles = BaseAnnot.getAllelesDescription(snpobj.getAlleles());
                        refSNPRefAllele = BaseAnnot.toString(snpobj.getAlleles()[0]);
                        refSNPAlleleAssessed = BaseAnnot.toString(snpobj.getAlleles()[1]);
                        flip[d] = false;
                    } else {
                        String SNPAlleles = BaseAnnot.getAllelesDescription(snpobj.getAlleles());
                        String SNPAlleleAssessed = BaseAnnot.toString(snpobj.getAlleles()[1]);
                        flip[d] = BaseAnnot.flipalleles(refSNPAlleles, refSNPAlleleAssessed, SNPAlleles, SNPAlleleAssessed);
                        if (flip[d] == null) {
                            log.writeln(snpobj.getName() + " has incompatible alleles. Expected: " + refSNPAlleles + " found " + SNPAlleles);
                            snpobj.clearGenotypes();
                        }
                    }
//                    } else {
//                        snpobj.clearGenotypes();
//                    }
                }
            }

            // clear array
            Arrays.fill(dosageData, -1);
            Arrays.fill(genotypeData, -1);

            // get data per dataset
            for (int d = 0; d < m_gg.length; d++) {
                if (flip[d] != null) {
                    boolean flipalleles = flip[d];
                    String[] individuals = m_gg[d].getIndividuals();
                    HashSet<String> allowedInds = allowedSamplesPerDs.get(d);

                    byte[] genotypes = null;
                    double[] dosages = null;
                    genotypes = snpObjs[d].getGenotypes();
                    if (snpObjs[d].hasDosageInformation()) {
                        dosages = snpObjs[d].getDosageValues();
                    }


                    for (int i = 0; i < individuals.length; i++) {
                        String ind = individuals[i];
                        if (allowedInds.contains(ind)) {
                            Integer newId = sampleIdMap.get(ind);
                            if (newId != null) {
                                if (genotypes[i] == -1) {
                                    genotypeData[newId] = -1;
                                    dosageData[newId] = -1;
                                } else {
                                    if (flipalleles) {
                                        genotypeData[newId] = 2 - genotypes[i];
                                        if (dosages != null) {
                                            dosageData[newId] = 2 - dosages[i];
                                        } else {
                                            dosageData[newId] = 2 - genotypes[i];
                                        }
                                    } else {
                                        genotypeData[newId] = genotypes[i];
                                        if (dosages != null) {
                                            dosageData[newId] = dosages[i];
                                        } else {
                                            dosageData[newId] = genotypes[i];
                                        }
                                    }
                                }
                            } else {
                                System.out.println("Sample " + ind + " has null id in index, but is present in GTE for dataset?");
                                System.exit(-1);
                            }
                        }
                    }
                } else if (snpObjs[d] != null) {
                    log.writeln("Excluding\t" + snp + "\tin dataset\t" + datasetsettings.get(d).name + "\tsince it has incompatible alleles: " + BaseAnnot.getAllelesDescription(snpObjs[d].getAlleles()));
                }
            }

            // count missing values
            int missing = 0;
            for (int d = 0; d < genotypeData.length; d++) {
                if (genotypeData[d] == -1) {
                    missing++;
                }
            }

            // skip empty variants
            if (missing == genotypeData.length) {
                System.out.println();

                log.writeln("Excluding\t" + snp + "\tsince it has no values");
            } else {

                // write SNP ID and mappings

                // write snp bytes
                byte[] a1 = new byte[genotypeData.length];
                byte[] a2 = new byte[genotypeData.length];
                byte[] dos = new byte[genotypeData.length];

                byte ba1 = BaseAnnot.toByte(refSNPRefAllele);
                byte ba2 = BaseAnnot.toByte(refSNPAlleleAssessed);

                for (int d = 0; d < genotypeData.length; d++) {
                    double vg = genotypeData[d];
                    double vd = dosageData[d];
                    int dosageInt = (int) Math.round(vd * 100d);
                    byte dosageByte = (byte) (Byte.MIN_VALUE + dosageInt);
                    if (vg < 0) {
                        a1[d] = 0;
                        a2[d] = 0;
                        dos[d] = 0;
                    } else {
                        dos[d] = dos[d] = dosageByte;
                        if (vg == 0) {
                            a1[d] = ba1;
                            a2[d] = ba1;
                        } else if (vg == 2) {
                            a1[d] = ba2;
                            a2[d] = ba2;
                        } else {
                            a1[d] = ba1;
                            a2[d] = ba2;
                        }
                    }
                }


                bfg.write(a1);
                bfg.write(a2);
                bfd.write(dos);
                snpm.writeln(refSNP.getChr() + "\t" + refSNP.getChrPos() + "\t" + refSNP.getName());
                snpo.writeln(refSNP.getName());
                written++;
            }

            for (int d = 0; d < snpObjs.length; d++) {
                if (snpObjs[d] != null) {
                    snpObjs[d].clearGenotypes();
                }
            }
            pb.iterate();
            if (snpctr % 1000000 == 0) {
                System.out.println("Written: " + written + " out of " + snpctr);
            }
        }
        pb.close();
        System.out.println("Written: " + written + " out of " + snpsToQuery.size());

        bfg.close();
        bfd.close();
        snpm.close();
        snpo.close();
        log.close();

    }

    private ArrayList<String> sortSNPs(ArrayList<String> snpsToQuery, boolean byId, TriTyperGenotypeData[] m_gg) {
        if (byId) {
            ArrayList<SNPObj> objs = new ArrayList<>();
            for (String s : snpsToQuery) {
                String[] elems = s.split(":");
                if (elems.length < 2) {
                    objs.add(new SNPObj(1, 1, s));
                } else if (elems.length >= 3) {
                    int chr = ChrAnnotation.parseChr(elems[0]);

                    int pos = Integer.parseInt(elems[1]);
                    objs.add(new SNPObj(chr, pos, s));

                }
            }

            Collections.sort(objs);
            ArrayList<String> output = new ArrayList<>();
            for (SNPObj s : objs) {
                output.add(s.id);
            }
            return output;
        } else {
            ArrayList<SNPObj> objs = new ArrayList<>();
            for (String s : snpsToQuery) {
                for (int d = 0; d < m_gg.length; d++) {
                    int id = m_gg[d].getSnpToSNPId().get(s);
                    if (id >= 0) {
                        int chr = m_gg[d].getChr(id);

                        int pos = m_gg[d].getChrPos(id);
                        objs.add(new SNPObj(chr, pos, s));
                        break;

                    }
                }
            }

            Collections.sort(objs);
            ArrayList<String> output = new ArrayList<>();
            for (SNPObj s : objs) {
                output.add(s.id);
            }
            return output;

        }
    }

    class SNPObj implements Comparable<SNPObj> {
        int chr;
        int pos;
        String id;

        public SNPObj(int chr, int pos, String s) {
            this.chr = chr;
            this.pos = pos;
            this.id = s;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;
            SNPObj snpObj = (SNPObj) o;
            return chr == snpObj.chr &&
                    pos == snpObj.pos &&
                    Objects.equals(id, snpObj.id);
        }

        @Override
        public int hashCode() {
            return Objects.hash(chr, pos, id);
        }

        @Override
        public int compareTo(SNPObj snpObj) {
            if (this.chr > snpObj.chr) {
                return 1;
            } else if (this.chr < snpObj.chr) {
                return -1;
            } else {
                if (this.pos > snpObj.pos) {
                    return 1;
                } else if (this.pos < snpObj.pos) {
                    return -1;
                } else {
                    return this.id.compareTo(snpObj.id);
                }
            }

        }
    }

    public void run(String settingsfile, String settingstexttoreplace, String settingstexttoreplacewith, boolean vcf, boolean sortbyId) throws IOException, ConfigurationException {
        System.out.println("TriTyper to dosagematrix converter.");
        Settings settings = new Settings();
        settings.settingsTextToReplace = settingstexttoreplace;
        settings.settingsTextReplaceWith = settingstexttoreplacewith;
        settings.load(settingsfile);


        // initialize datasets
        ArrayList<TriTyperGeneticalGenomicsDatasetSettings> datasetsettings = settings.datasetSettings;
        TriTyperGenotypeData[] m_gg = new TriTyperGenotypeData[datasetsettings.size()];
        SNPLoader[] loaders = new SNPLoader[datasetsettings.size()];
        System.out.println("Initializing " + datasetsettings.size() + " datasets in parallel.");
        IntStream.range(0, m_gg.length).parallel().forEach(i -> {
            m_gg[i] = new TriTyperGenotypeData();
            try {
                m_gg[i].load(datasetsettings.get(i).genotypeLocation, datasetsettings.get(i).snpmapFileLocation, datasetsettings.get(i).snpFileLocation);
                loaders[i] = m_gg[i].createSNPLoader(1);
            } catch (IOException e) {
                e.printStackTrace();
            }
        });


        System.out.println("Selecting samples per dataset");
        ArrayList<HashSet<String>> allowedSamplesPerDs = new ArrayList<HashSet<String>>();
        for (int d = 0; d < m_gg.length; d++) {
            TriTyperGeneticalGenomicsDatasetSettings dataset = datasetsettings.get(d);
            HashSet<String> allowedSamples = new HashSet<String>();
            if (dataset.genotypeToExpressionCoupling == null || dataset.genotypeToExpressionCoupling.length() == 0) {
                allowedSamples.addAll(Arrays.asList(m_gg[d].getIndividuals()));
            } else {
                TextFile tf = new TextFile(dataset.genotypeToExpressionCoupling, TextFile.R);
                String[] elems = tf.readLineElems(TextFile.tab);
                while (elems != null) {
                    allowedSamples.add(elems[0]);
                    elems = tf.readLineElems(TextFile.tab);
                }
                tf.close();
            }
            System.out.println(datasetsettings.get(d).name + " has " + allowedSamples.size() + " samples selected.");
            allowedSamplesPerDs.add(allowedSamples);
        }

        // create sample index
        System.out.println("Creating sample index..");
        HashMap<String, Integer> sampleIdMap = new HashMap<String, Integer>();
        int sampleCounter = 0;
        for (int i = 0; i < allowedSamplesPerDs.size(); i++) {
            HashSet<String> samplesInDs = allowedSamplesPerDs.get(i);
            for (String s : samplesInDs) {
                sampleIdMap.put(s, sampleCounter);
                sampleCounter++;
            }
        }
        System.out.println(sampleIdMap.size() + " unique sample IDs.");

        // sorting SNPs
        Integer selectChr = null;
        if (settings.confineToSNPsThatMapToChromosome != null) {
            selectChr = (int) settings.confineToSNPsThatMapToChromosome;
            System.out.println("Selecting variants from chr " + selectChr);
        }

        System.out.println("Inventorizing snps");
        HashSet<String> allSNPsHash = new HashSet<String>();
        for (int d = 0; d < m_gg.length; d++) {
            String[] datasetSNPs = m_gg[d].getSNPs();
            if (selectChr != null) {
                for (int s = 0; s < datasetSNPs.length; s++) {
                    int chr = m_gg[d].getChr(s);
                    if (selectChr == chr) {
                        allSNPsHash.add(datasetSNPs[s]);
                    }
                }
            } else {
                allSNPsHash.addAll(Arrays.asList(datasetSNPs));
            }
        }

        System.out.println(allSNPsHash.size() + " unique SNPs");

        System.out.println("Sorting SNPs");
        ArrayList<String> allSNPs = new ArrayList<>();
        allSNPs.addAll(allSNPsHash);

        ArrayList<String> snpsToQuery = null;
        if (settings.tsSNPsConfine != null) {
            snpsToQuery = new ArrayList<>();
            HashSet<String> confineset = settings.tsSNPsConfine;
            for (String s : confineset) {
                if (allSNPsHash.contains(s)) {
                    snpsToQuery.add(s);
                }
            }
        } else {
            snpsToQuery = allSNPs;
        }
        System.out.println(snpsToQuery.size() + " SNPs remain after filtering for SNP confinement list");


        snpsToQuery = sortSNPs(snpsToQuery, sortbyId, m_gg);
        System.out.println(snpsToQuery.size() + " SNPs to output");

        // write data to disk
        String[] header = new String[sampleIdMap.size()];
        for (Map.Entry<String, Integer> samplePair : sampleIdMap.entrySet()) {
            header[samplePair.getValue()] = samplePair.getKey();
        }
        TextFile tf = null;
        TextFile log = new TextFile(settings.outputReportsDir + "MergeLog.txt.gz", TextFile.W);
        if (vcf) {
            tf = new TextFile(settings.outputReportsDir + "GenotypeData.vcf.gz", TextFile.W);
            tf.writeln("##fileformat=VCFv4.3");
            tf.writeln("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Non phased Genotype\">");
            tf.writeln("##FORMAT=<ID=DS,Number=1,Type=String,Description=\"Genotype dosage\">");
            // #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT
            tf.writeln("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + Strings.concat(header, Strings.tab));
        } else {
            tf = new TextFile(settings.outputReportsDir + "GenotypeData.txt.gz", TextFile.W);
            tf.writeln("SNP\tAlleles\tMinorAllele\t" + Strings.concat(header, Strings.tab));
        }


        double[] data = new double[sampleIdMap.size()];
        log.write(snpsToQuery.size() + " SNPs to output");
        ProgressBar pb = new ProgressBar(snpsToQuery.size(), "Querying SNPs...");
        for (int snpctr = 0; snpctr < snpsToQuery.size(); snpctr++) {
            String snp = snpsToQuery.get(snpctr);
            SNP[] snpObjs = new SNP[m_gg.length];
            // load genotype data per dataset
            String refSNPAlleles = null;
            String refSNPRefAllele = null;
            String refSNPAlleleAssessed = null;
            SNP refSNP = null;
            Boolean[] flip = new Boolean[m_gg.length];
            for (int d = 0; d < m_gg.length; d++) {
                Integer id = m_gg[d].getSnpToSNPId().get(snp);
                if (id >= 0) {
                    SNP snpobj = m_gg[d].getSNPObject(id);
                    loaders[d].loadGenotypes(snpobj);
                    if (snpobj.getMAF() > settings.snpQCMAFThreshold && snpobj.getHWEP() > settings.snpQCHWEThreshold && snpobj.getCR() > settings.snpQCCallRateThreshold) {
                        snpObjs[d] = snpobj;
                        if (loaders[d].hasDosageInformation()) {
                            loaders[d].loadDosage(snpObjs[d]);
                        }
                        if (refSNPAlleles == null) {
                            refSNP = snpobj;
                            refSNPAlleles = BaseAnnot.getAllelesDescription(snpobj.getAlleles());
                            refSNPRefAllele = BaseAnnot.toString(snpobj.getAlleles()[0]);
                            refSNPAlleleAssessed = BaseAnnot.toString(snpobj.getAlleles()[1]);
                            flip[d] = false;
                        } else {
                            String SNPAlleles = BaseAnnot.getAllelesDescription(snpobj.getAlleles());
                            String SNPAlleleAssessed = BaseAnnot.toString(snpobj.getAlleles()[1]);
                            flip[d] = BaseAnnot.flipalleles(refSNPAlleles, refSNPAlleleAssessed, SNPAlleles, SNPAlleleAssessed);
                            if (flip[d] == null) {
                                snpobj.clearGenotypes();
                            }
                        }
                    } else {
                        snpobj.clearGenotypes();
                    }
                }
            }

            // flip alleles

            for (int d = 0; d < data.length; d++) {
                data[d] = -1;
            }

            // get data per dataset
            for (int d = 0; d < m_gg.length; d++) {
                if (flip[d] != null) {
                    boolean flipalleles = flip[d];
                    String[] individuals = m_gg[d].getIndividuals();
                    HashSet<String> allowedInds = allowedSamplesPerDs.get(d);

                    double[] dosages = null;
                    if (snpObjs[d].hasDosageInformation()) {
                        dosages = snpObjs[d].getDosageValues();
                    } else {
                        byte[] genotypes = snpObjs[d].getGenotypes();
                        dosages = new double[genotypes.length];
                        for (int g = 0; g < dosages.length; g++) {
                            dosages[g] = genotypes[g];
                        }
                    }


                    for (int i = 0; i < individuals.length; i++) {
                        String ind = individuals[i];
                        if (allowedInds.contains(ind)) {
                            Integer newId = sampleIdMap.get(ind);
                            if (newId != null) {
                                if (dosages[i] == -1) {
                                    data[newId] = -1;
                                } else {
                                    if (flipalleles) {
                                        data[newId] = 2 - dosages[i];
                                    } else {
                                        data[newId] = dosages[i];
                                    }
                                }
                            } else {
                                System.out.println("Sample " + ind + " has null id in index, but is present in GTE for dataset?");
                                System.exit(-1);
                            }
                        }
                    }
                } else if (snpObjs[d] != null) {
                    log.writeln("Excluding\t" + snp + "\tin dataset\t" + datasetsettings.get(d).name + "\tsince it has incompatible alleles: " + BaseAnnot.getAllelesDescription(snpObjs[d].getAlleles()));
                }
            }

            int missing = 0;
            for (int d = 0; d < data.length; d++) {
                if (data[d] == -1) {
                    missing++;
                }
            }
            if (missing == data.length) {
                log.write(snp + " has no data");
            } else if (missing != data.length) {
                // System.out.println();
                // System.out.println("Excluding\t" + snp + "\tsince it has no values");
                if (vcf) {
                    // #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT


                    // ./.:0,0:0
                    String[] gts = new String[data.length];
                    for (int d = 0; d < data.length; d++) {
                        double v = data[d];
                        if (v < 0) {
                            gts[d] = "./.:0";
                        } else {
                            if (v < 0.5) {
                                gts[d] = "0/0:" + data[d];
                            } else if (v > 1.5) {
                                gts[d] = "1/1:" + data[d];
                            } else {
                                gts[d] = "0/1:" + data[d];
                            }
                        }
                    }
                    tf.writeln(refSNP.getChr()
                            + "\t" + refSNP.getChrPos()
                            + "\t" + refSNP.getName()
                            + "\t" + refSNPRefAllele
                            + "\t" + refSNPAlleleAssessed
                            + "\t.\t.\t.\tGT:DS\t"
                            + Strings.concat(gts, Strings.tab));

                } else {
                    tf.writeln(snp + "\t" + refSNPAlleles + "\t" + refSNPAlleleAssessed + "\t" + Strings.concat(data, Strings.tab));
                }
            }
            pb.iterate();
        }
        pb.close();

        tf.close();
    }
}
