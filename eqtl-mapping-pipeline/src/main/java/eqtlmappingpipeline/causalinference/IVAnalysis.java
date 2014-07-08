/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.causalinference;

import eqtlmappingpipeline.metaqtl3.containers.Settings;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import umcg.genetica.console.ConsoleGUIElems;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.containers.Triple;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.SNPLoader;
import umcg.genetica.io.trityper.TriTyperGeneticalGenomicsDataset;
import umcg.genetica.io.trityper.TriTyperGeneticalGenomicsDatasetSettings;
import umcg.genetica.math.stats.Correlation;
import umcg.genetica.math.stats.Regression;
import umcg.genetica.math.stats.TwoStepLeastSquares;

/**
 *
 * @author harmjan
 */
public class IVAnalysis {

    protected HashSet<Triple<String, String, String>> snpProbeCombos = new HashSet<Triple<String, String, String>>();
    protected Settings m_settings;
    // snpProbeCombinationList (SNP cis trans)
    protected TriTyperGeneticalGenomicsDataset[] m_gg;
    protected final String outDir;

    public IVAnalysis(String xmlSettingsFile,
            String ingt, String inexp, String inexpplatform, String inexpannot,
            String gte, String out, int perm, String snpProbeCombinationList) throws IOException, Exception {

        if (xmlSettingsFile == null && (ingt == null || inexp == null)) {
            throw new IllegalArgumentException("Supply settingsfile for IV Analysis");
        }

        if (snpProbeCombinationList == null) {
            throw new IllegalArgumentException("Supply SNP probe combination file for IV Analysis");
        }

        loadSNPProbeCombos(snpProbeCombinationList);
        if (snpProbeCombos.isEmpty()) {
            throw new IllegalArgumentException("SNP Probe combination file is empty!");
        }

        initializeDatasets(xmlSettingsFile, ingt, inexp, inexpplatform, inexpannot, gte, out, perm);

        if (!out.endsWith("/")) {
            out += "/";
        }

        Gpio.createDir(out);
        this.outDir = out;


    }

    private void loadSNPProbeCombos(String snpProbeCombinationList) throws IOException {
        System.out.println("Loading SNP-Cis-Trans combo's from file: " + snpProbeCombinationList);
        TextFile tf = new TextFile(snpProbeCombinationList, TextFile.R);
        String[] elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {
            if (elems.length >= 3) {
                Triple<String, String, String> t = new Triple<String, String, String>(elems[0], elems[1], elems[2]);
                snpProbeCombos.add(t);
            }
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();
        System.out.println("Loaded " + snpProbeCombos.size() + " combinations.");
    }

    private void initializeDatasets(String xmlSettingsFile,
            String ingt, String inexp, String inexpplatform, String inexpannot,
            String gte, String out, int perm) throws IOException, Exception {


        if (m_settings == null && xmlSettingsFile == null && ingt != null) {

            // check the settings
            boolean settingsOk = true;
            if (inexp == null || inexp.trim().length() == 0) {
                System.err.println("ERROR: you did not specify a gene expression file.");
                settingsOk = false;
            }

            if (inexpannot != null && inexpannot.trim().length() != 0) {
                if (inexpplatform == null || inexpplatform.trim().length() == 0) {
                    System.err.println("ERROR: you specified " + inexpannot + " but you did not specify the platform (using --inexpplatform)!");
                    settingsOk = false;
                }
            }

            if (out == null || out.trim().length() == 0) {
                System.err.println("ERROR: you did not specify an output directory.");
                settingsOk = false;
            }

            if (!settingsOk) {
                System.out.println();
                System.exit(0);
            }

            m_settings = new Settings();
            TriTyperGeneticalGenomicsDatasetSettings s = new TriTyperGeneticalGenomicsDatasetSettings();

            s.name = "Dataset";
            s.expressionLocation = inexp;
            s.expressionplatform = inexpplatform;
            s.probeannotation = inexpannot;
            s.genotypeLocation = ingt;
            s.genotypeToExpressionCoupling = gte;
            s.cisAnalysis = true;
            s.transAnalysis = true;

            m_settings.cisAnalysis = true;
            m_settings.transAnalysis = true;

            boolean cistrans = false;
            if (m_settings.cisAnalysis && m_settings.transAnalysis) {
                m_settings.confineProbesToProbesMappingToAnyChromosome = true;
            }

            m_settings.datasetSettings = new ArrayList<TriTyperGeneticalGenomicsDatasetSettings>();

            m_settings.regressOutEQTLEffectFileName = null;
            m_settings.datasetSettings.add(s);
            m_settings.nrThreads = 1;
            m_settings.cisAnalysis = true;
            m_settings.transAnalysis = true;
            m_settings.nrPermutationsFDR = perm;
            if (!out.endsWith("/")) {
                out += "/";
            }
            if (!Gpio.exists(out)) {
                Gpio.createDir(out);
            }


            m_settings.outputReportsDir = out;
            m_settings.createTEXTOutputFiles = true;
            m_settings.createBinaryOutputFiles = false;

        } else if (m_settings == null && xmlSettingsFile != null) {
            // parse settings
            m_settings = new Settings();
            m_settings.load(xmlSettingsFile);
        } else if (m_settings == null) {
            System.out.println("ERROR: No input specified");
            System.exit(0);
        }

        // initialize dataset
        if (!m_settings.cisAnalysis && !m_settings.transAnalysis) {
            System.err.println("! WARNING: defaulting to CIS analysis (override with --trans or --trans and --cis))");
            m_settings.cisAnalysis = true;
        }

        m_settings.writeSettingsToDisk();


        int numDatasets = m_settings.datasetSettings.size();
        m_gg = new TriTyperGeneticalGenomicsDataset[numDatasets];

        int nrOfDatasetsWithGeneExpressionData = 0;
        for (int i = 0; i < numDatasets; i++) {

            System.out.println("- Loading dataset: " + m_settings.datasetSettings.get(i).name + "");
            m_settings.datasetSettings.get(i).confineProbesToProbesMappingToAnyChromosome = m_settings.confineProbesToProbesMappingToAnyChromosome;
            System.out.println(ConsoleGUIElems.LINE);
            m_gg[i] = new TriTyperGeneticalGenomicsDataset(m_settings.datasetSettings.get(i));

            if (m_gg[i].isExpressionDataLoadedCorrectly()) {
                nrOfDatasetsWithGeneExpressionData++;
            }

        }

        if (nrOfDatasetsWithGeneExpressionData == 0) {
            System.out.println("Error: none of your datasets contain any gene expression data for the settings you have specified");
            System.exit(0);
        }


    }

    public void run() throws IOException {
        for (int d = 0; d < m_gg.length; d++) {
            // now test all triples
            SNPLoader snpLoader = m_gg[d].getGenotypeData().createSNPLoader();
            int[] indWGA = m_gg[d].getExpressionToGenotypeIdArray();

            for (int perm = 0; perm < m_settings.nrPermutationsFDR + 1; perm++) {
                String outfile = null;
                if (perm == 0) {
                    outfile = outDir + m_gg[d].getSettings().name + "_IVAnalysis-RealData.txt";
                } else {
                    outfile = outDir + m_gg[d].getSettings().name + "_IVAnalysis-PermutationRound-" + perm + ".txt";
                    m_gg[d].permuteSampleLables(m_settings.r);
                }
                TextFile out = new TextFile(outfile, TextFile.W);
                Iterator<Triple<String, String, String>> it = snpProbeCombos.iterator();
                Triple<String, String, String> next = it.next();
                ProgressBar pb = new ProgressBar(snpProbeCombos.size(), "Running IV Analysis - Permutation " + perm);
//                out.writeln("SNP\tSNPId\tCisProbe\tCisProbeId\tCisGeneName\tTransProbe\tTransProbeId\tTransGeneName\tSNPTrans-Beta\tSNPTrans-Alpha\tSNPTrans-SE\tCisTrans-Beta\tCisTrans-Alpha\tCisTrans-SE\tIV-Beta\tIV-Alpha\tIV-SE\tcorcistrans\tr2cistrans\tcorResCisResTrans\tr2rescisrestrans\tcorCisResTrans\tcorSNPTrans\tcorSNPResTrans\tcorrSNPResCisResTrans");
                out.writeln("SNP\t"
                        + "CisArrayAddress\t"
                        + "CisGeneName"
                        + "TransArrayAddress\t"
                        + "TransGeneName"
                        + "cisEQTL\t"
                        + "r2cisEQTL\t"
                        + "transEQTL\t"
                        + "r2TransEQTL\t"
                        + "corCisTrans\t"
                        + "r2CisTrans\t"
                        + "corResCisResTrans\t"
                        + "r2ResCisResTrans\t"
                        + "IV-Beta\t"
                        + "IV-SE");
                while (next != null) {
                    String snp = next.getLeft();
                    String cisprobe = next.getMiddle();
                    String transprobe = next.getRight();

                    Integer snpId = m_gg[d].getGenotypeData().getSnpToSNPId().get(snp);
                    Integer cisProbeId = m_gg[d].getExpressionData().getProbeToId().get(cisprobe);
                    Integer transProbeId = m_gg[d].getExpressionData().getProbeToId().get(transprobe);

                    if (snpId == -9 || cisProbeId == null || transProbeId == null) {
//                        out.writeln(snp + "\t" + snpId + "\t" + cisprobe + "\t" + cisProbeId + "\t" + null + "\t" + transprobe + "\t" + transProbeId + "\t" + null + "\t" + null + "\t" + null + "\t" + null + "\t" + null + "\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA");
                    } else {

                        SNP snpObj = m_gg[d].getGenotypeData().getSNPObject(snpId);
                        snpLoader.loadGenotypes(snpObj);
                        snpLoader.loadDosage(snpObj);



                        double[] origCisVals = m_gg[d].getExpressionData().getMatrix()[cisProbeId];
                        double[] origTransVals = m_gg[d].getExpressionData().getMatrix()[transProbeId];

                        int nrequal = 0;
                        int calledGenotypes = 0;
                        for (int i = 0; i < m_gg[d].getExpressionData().getIndividuals().length; i++) {
                            int genotypeId = indWGA[i];
                            short gt = snpObj.getGenotypes()[genotypeId];
//                            if(m_gg[d].getExpressionData().getIndividuals()[i].equals(m_gg[d].getGenotypeData().getIndividuals()[genotypeId])){
//                                nrequal++;
//                            }
//                            System.out.println(m_gg[d].getExpressionData().getIndividuals()[i] + "\t" + m_gg[d].getGenotypeData().getIndividuals()[genotypeId]);
                            if (genotypeId > -1 && gt > -1) {
                                calledGenotypes++;

                            }
                        }
//                        System.out.println(nrequal+" inds equal..");


                        double[] genotypes = new double[calledGenotypes];
                        double[] cisvals = new double[calledGenotypes];
                        double[] transvals = new double[calledGenotypes];

                        calledGenotypes = 0;
                        for (int i = 0; i < m_gg[d].getExpressionData().getIndividuals().length; i++) {
                            int genotypeId = indWGA[i];
                            short gt = snpObj.getGenotypes()[genotypeId];
                            if (genotypeId > -1 && gt > -1) {
                                genotypes[calledGenotypes] = snpObj.getDosageValues()[genotypeId];
                                cisvals[calledGenotypes] = origCisVals[i];
                                transvals[calledGenotypes] = origTransVals[i];
                                calledGenotypes++;
                            }
                        }




                        double corrCisTrans = JSci.maths.ArrayMath.correlation(cisvals, transvals);

                        double[] snpCisRCs = Regression.getLinearRegressionCoefficients(genotypes, cisvals);
                        double[] snpTransRCs = Regression.getLinearRegressionCoefficients(genotypes, transvals);

                        double[] resCis = new double[cisvals.length];
                        double[] resTrans = new double[cisvals.length];
                        for (int i = 0; i < resCis.length; i++) {
                            resCis[i] = cisvals[i] - snpCisRCs[0] * genotypes[i];
                            resTrans[i] = transvals[i] - snpTransRCs[0] * genotypes[i];
                        }

                        double corrResCisResTrans = JSci.maths.ArrayMath.correlation(resCis, resTrans);



                        double[] cisTransRCs = Regression.getLinearRegressionCoefficients(cisvals, transvals);
                        double[] resCisTransRCs = Regression.getLinearRegressionCoefficients(resCis, transvals);
                        double[] resCisTrans = new double[cisvals.length];
                        double[] resResCisTrans = new double[cisvals.length];
                        for (int i = 0; i < resCisTrans.length; i++) {
                            resCisTrans[i] = transvals[i] - (cisTransRCs[0] * cisvals[i]);
                            resResCisTrans[i] = transvals[i] - (resCisTransRCs[0] * cisvals[i]);
                        }

                        double corrCisResTrans = Correlation.correlate(cisvals, resCisTrans);

                        double transEQTL = Correlation.correlate(genotypes, transvals);
                        double corrSNPResTrans = Correlation.correlate(genotypes, resCisTrans);
                        double corrSNPResCisResTrans = Correlation.correlate(genotypes, resResCisTrans);


//                        if (snp.equals("rs12718597") && transprobe.equals("4900309")) {
////                            for (int i = 0; i < genotypes.length; i++) {
////                                System.out.println(genotypes[i] + "\t" + cisvals[i] + "\t" + transvals[i]);
////                            }
//                            double[] result = TwoStepLeastSquares.tsls(transvals, cisvals, genotypes, true);
//                            System.exit(0);
//                        }

//                        double m1 = JSci.maths.ArrayMath.mean(genotypes);
//                        double m2 = JSci.maths.ArrayMath.mean(transvals);
//                        double m3 = JSci.maths.ArrayMath.mean(cisvals);
//
//                        double s1 = JSci.maths.ArrayMath.standardDeviation(genotypes);
//                        double s2 = JSci.maths.ArrayMath.standardDeviation(transvals);
//                        double s3 = JSci.maths.ArrayMath.standardDeviation(cisvals);
//
//                        for (int i = 0; i < genotypes.length; i++) {
//                            genotypes[i] -= m1;
//                            genotypes[i] /= s1;
//
//                            transvals[i] -= m2;
//                            transvals[i] /= s2;
//
//                            cisvals[i] -= m3;
//                            cisvals[i] /= s3;
//                        }
//
//                        RankDoubleArray rda = new RankDoubleArray();
//                        genotypes = rda.rank(genotypes);
//                        cisvals = rda.rank(cisvals);
//                        transvals = rda.rank(transvals);

                        // we now have all the values we need.. perform two-step OLS

                        // perform regression of snp against trans probe
//                        double[] stats = Regression.getLinearRegressionCoefficients(genotypes, transvals);
//                        double[] stats2 = Regression.getLinearRegressionCoefficients(cisvals, transvals);

                        double[] result = TwoStepLeastSquares.tsls(transvals, cisvals, genotypes);

                        double cisEQTL = Correlation.correlate(genotypes, cisvals);

//                        out.writeln(snp + "\t" + snpId + "\t" + cisprobe + "\t" + cisProbeId + "\t" + m_gg[d].getExpressionData().getAnnotation()[cisProbeId] + "\t" + transprobe + "\t" + transProbeId + "\t" + m_gg[d].getExpressionData().getAnnotation()[transProbeId] + "\t" + stats[0] + "\t" + stats[1] + "\t" + stats[2] + "\t" + stats2[0] + "\t" + stats2[1] + "\t" + stats2[2] + "\t" + result[0] + "\t" + result[1] + "\t" + result[2] + "\t" + corrCisTrans + "\t" + r2CisTrans + "\t" + corrResCisResTrans + "\t" + r2ResCisResTrans + "\t" + corrCisResTrans + "\t" + corSNPTrans + "\t" + corrSNPResTrans+"\t"+corrSNPResCisResTrans);

                        out.write(snp
                                + "\t" + cisprobe
                                + "\t" + m_gg[d].getExpressionData().getAnnotation()[cisProbeId]
                                + "\t" + transprobe
                                + "\t" + m_gg[d].getExpressionData().getAnnotation()[transProbeId]
                                + "\t" + cisEQTL
                                + "\t" + (cisEQTL * cisEQTL)
                                + "\t" + transEQTL
                                + "\t" + (transEQTL * transEQTL)
                                + "\t" + corrCisTrans
                                + "\t" + (corrCisTrans * corrCisTrans)
                                + "\t" + corrResCisResTrans
                                + "\t" + (corrResCisResTrans * corrResCisResTrans)
                                + "\t" + result[0]
                                + "\t" + result[2]);
                        snpObj.clearGenotypes();
                    }


                    if (it.hasNext()) {
                        next = it.next();
                    } else {
                        next = null;
                    }
                    pb.iterate();
                }
                pb.close();
                out.close();
            }

            snpLoader.close();


        }

    }
}
