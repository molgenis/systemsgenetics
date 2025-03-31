package eqtlmappingpipeline.interactionMeta;

import com.opencsv.CSVParser;
import com.opencsv.CSVParserBuilder;
import com.opencsv.CSVReader;
import com.opencsv.CSVReaderBuilder;
import me.tongfei.progressbar.ProgressBarStyle;
import me.tongfei.progressbar.ProgressBar;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.matrix2.DoubleMatrixDatasetRowCompressedReader;

import java.io.*;
import java.util.*;
import java.util.stream.IntStream;
import java.util.zip.GZIPInputStream;

public class InteractionMeta {

    public static void runMeta() throws IOException {

        File resultFolder = new File("/groups/umcg-fg/tmp04/projects/eqtlgen-phase2/interactions/replicationMarcJan5");
        String[] cohorts = new String[]{"LL", "LLS_660Q", "LLS_OmniExpr", "NTR_Affy6", "NTR_GoNL", "RS"};//
        ArrayList<File> cohortFolders = new ArrayList<File>();

        for (String cohort : cohorts) {
            cohortFolders.add(new File(resultFolder, cohort));
        }

        File allTestEffectsFile = new File("/groups/umcg-fg/tmp04/projects/eqtlgen-phase2/interactions/replicationMarcJan3/co_eQTL_significantTests_inc_all_sceQTLGen_eGenes2.txt.gz");
        //File allTestEffectsFile = new File("/groups/umcg-fg/tmp04/projects/eqtlgen-phase2/interactions/biosInteractionsPcCor50/nod2Qtls.txt.gz");

        File metaZFile = new File(resultFolder, "metaZTest.txt");
        File metaSampleCountFile = new File(resultFolder, "metaSampleCountTest.txt");

        runMeta(cohortFolders, allTestEffectsFile, metaZFile, metaSampleCountFile, "/interactionZscoreTest.datg");

        metaZFile = new File(resultFolder, "metaZPermutation.txt");
        metaSampleCountFile = new File(resultFolder, "metaSampleCountPermutation.txt");

        //this does not work yet. In permutation row names contain round number
//        runMeta(cohortFolders, allTestEffectsFile, metaZFile, metaSampleCountFile, "/interactionZscorePermutation.datg");

    }


    private static void runMeta(List<File> cohortResultFolders, File allTestEffectsFile, File metaZFile, File metaSampleCountFile, String zscoreFile) throws IOException {


        int cores = Runtime.getRuntime().availableProcessors();
        System.out.println("Number of available cores: " + cores);

        //make sure there are other treads that can run if one is waiting for IO
        System.setProperty("java.util.concurrent.ForkJoinPool.common.parallelism", String.valueOf(cores * 4));

        List<Eqtl> testedEqtls = loadTestedEqtls(allTestEffectsFile);

        ArrayList<String> testedEqtlsNames = new ArrayList<>(testedEqtls.size());
        HashSet<String> eqtlGenes = new LinkedHashSet<>();
        for (Eqtl eqtl : testedEqtls) {
            testedEqtlsNames.add(eqtl.getSnpGene());
            eqtlGenes.add(eqtl.getGene());
        }
        System.out.println("Total eqtls: " + testedEqtlsNames.size());
        System.out.println("Unique eqtl genes: " + eqtlGenes.size());


        DoubleMatrixDatasetRowCompressedReader[] cohortsZscores = new DoubleMatrixDatasetRowCompressedReader[cohortResultFolders.size()];
        int[][] cohortsColumnIndices = new int[cohortResultFolders.size()][];
        String[] cohortNames = new String[cohortResultFolders.size()];

        HashSet<String> foundCovariates = new LinkedHashSet<>();
        for (int cohortI = 0; cohortI < cohortResultFolders.size(); cohortI++) {

            DoubleMatrixDatasetRowCompressedReader cohortZscores = new DoubleMatrixDatasetRowCompressedReader(cohortResultFolders.get(cohortI).getPath() + zscoreFile);
            cohortsZscores[cohortI] = cohortZscores;
            cohortNames[cohortI] = cohortResultFolders.get(cohortI).getName();

            foundCovariates.addAll(cohortZscores.getColumnIdentifiers());

        }

//        System.out.println("DEBUG MODE!!!!!!!!!!!!!!!");
//
//        foundCovariates.clear();
//        foundCovariates.add("ENSG00000166900");

        final String[] covariates = foundCovariates.toArray(new String[0]);

        System.out.println("Covariates: " + covariates.length);

        DoubleMatrixDataset<String, String> eqtlgeneNonMissingCount = new DoubleMatrixDataset<>(eqtlGenes, Arrays.asList(cohortNames));

        for (int cohortI = 0; cohortI < cohortResultFolders.size(); cohortI++) {
            System.out.println(cohortResultFolders.get(cohortI).getName());

            loadCohortEqtlgeneCounts(new File(cohortResultFolders.get(cohortI).getPath(), "/feature_metadata.txt.gz"), cohortResultFolders.get(cohortI).getName(), eqtlgeneNonMissingCount);

            int[] cohortColumnIndices = new int[covariates.length];

            DoubleMatrixDatasetRowCompressedReader cohortZscores = cohortsZscores[cohortI];

            Map<String, Integer> cohortColumnMap = cohortZscores.getColumnMap();

            for (int covariateI = 0; covariateI < covariates.length; covariateI++) {

                Integer colIndex = cohortColumnMap.get(covariates[covariateI]);

                if (colIndex == null) {
                    cohortColumnIndices[covariateI] = -1;
                } else {
                    cohortColumnIndices[covariateI] = colIndex;
                }

            }

            cohortsColumnIndices[cohortI] = cohortColumnIndices;

        }


        DoubleMatrixDataset<String, String> metaZscores = new DoubleMatrixDataset<>(testedEqtlsNames, foundCovariates);
        DoubleMatrixDataset<String, String> totalSamples = new DoubleMatrixDataset<>(testedEqtlsNames, foundCovariates);

        ProgressBar pb = new ProgressBar("Meta analysis", testedEqtlsNames.size(), ProgressBarStyle.COLORFUL_UNICODE_BLOCK);


        // int cohortWithAnyData;
        IntStream.range(0, testedEqtls.size()).parallel().forEach(eqtlI -> {


            DoubleMatrixDataset<String, String> intermediates = new DoubleMatrixDataset<>(Arrays.asList("sumWZ", "sumW2"), foundCovariates);

            final Eqtl eqtl = testedEqtls.get(eqtlI);

 //           try {

//                BufferedWriter writer = new BufferedWriter(new FileWriter(new File(metaZFile.getParent(), "debug_" + eqtl.snp + "-" + eqtl.gene)));
//
//
//                writer.write(eqtl.getSnpGene());
//                writer.newLine();

                //cohortWithAnyData = 0;

                //reset intermediates
                intermediates.getMatrix().assign(0);

                for (int cohortI = 0; cohortI < cohortResultFolders.size(); cohortI++) {

//                    writer.write(cohortNames[cohortI]);
//                    writer.newLine();

                    final DoubleMatrixDatasetRowCompressedReader cohortZscores = cohortsZscores[cohortI];

//                boolean anyHasZ = false;

                    //System.out.println(eqtl.getSnpGene());

                    if (cohortZscores.hasRow(eqtl.getSnpGene())) {


                        final String cohort = cohortNames[cohortI];


                        final double[] cohortZscoresArray;
                        try {
                            cohortZscoresArray = cohortZscores.loadSingleRow(eqtl.getSnpGene());
                        } catch (IOException e) {
                            throw new RuntimeException(e);
                        }

//                    final double[] cohortZscoresArray = new double[cohortZscores.getNumberOfColumns()];


                        final double sampleSize = eqtlgeneNonMissingCount.getElement(eqtl.getGene(), cohort);
                        final double sqrtSampleSize = Math.sqrt(sampleSize);

                        for (int covariateI = 0; covariateI < covariates.length; covariateI++) {


                            // System.out.println(covariate);

                            final int covariantIndexInCohort = cohortsColumnIndices[cohortI][covariateI];
                            if (covariantIndexInCohort >= 0) {

                                final double zScore = cohortZscoresArray[covariantIndexInCohort];
//
//                                writer.write("Zscore: " + zScore);
//                                writer.newLine();
//
//                                writer.write("Sample size: " + sampleSize);
//                                writer.newLine();


////                            System.out.println(zScore);
//
                                if (!Double.isNaN(zScore)) {

                                    intermediates.setElementQuick(0, covariateI, intermediates.getElementQuick(0, covariateI) + (sqrtSampleSize * zScore));
                                    intermediates.setElementQuick(1, covariateI, intermediates.getElementQuick(1, covariateI) + (sampleSize));
                                    totalSamples.setElementQuick(eqtlI, covariateI, totalSamples.getElementQuick(eqtlI, covariateI) + sampleSize);

//
                                }

                            }

                        }

                    }
//                    else {
//                        writer.write("No data");
//                    }

                }

                //Here we have looped through all cohorts

                //System.out.println(eqtl.getSnpGene() + "\t" + cohortWithAnyData);

                //intermediates.viewDice().save("/groups/umcg-fg/tmp04/projects/eqtlgen-phase2/interactions/biosInteractions/intermediates.txt");


                for (int covariateI = 0; covariateI < foundCovariates.size(); covariateI++) {
                    double sumWZ = intermediates.getElementQuick(0, covariateI);
                    double sumW2 = intermediates.getElementQuick(1, covariateI);

//                    writer.write("SumWZ: " + sumWZ);
//                    writer.newLine();
//                    writer.write("SumW2: " + sumW2);
//                    writer.newLine();
//                    writer.write("Z: " + sumWZ / Math.sqrt(sumW2));
//                    writer.newLine();

                    if (sumW2 > 0) {

                        metaZscores.setElementQuick(eqtlI, covariateI, sumWZ / Math.sqrt(sumW2));

                    }

                }

//                writer.close();

//            } catch (IOException e) {
//                throw new RuntimeException(e);
//            }

            pb.step();

        });


        pb.close();

        metaZscores.save(metaZFile);
        totalSamples.save(metaSampleCountFile);

    }

    private static List<Eqtl> loadTestedEqtls(File allTestEffectsFile) throws IOException {

        ArrayList<Eqtl> testedEqtls = new ArrayList<>();

        final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
        final CSVReader reader;
        if (allTestEffectsFile.getName().endsWith(".gz")) {
            reader = new CSVReaderBuilder((new BufferedReader(new InputStreamReader((new GZIPInputStream(new FileInputStream(allTestEffectsFile))))))).withCSVParser(parser).build();
        } else {
            reader = new CSVReaderBuilder((new BufferedReader(new InputStreamReader((new FileInputStream(allTestEffectsFile)))))).withCSVParser(parser).build();
        }

        //First nextLine contains the header
        String[] nextLine = reader.readNext();

        while ((nextLine = reader.readNext()) != null) {
            testedEqtls.add(new Eqtl(nextLine[1], nextLine[0]));
        }

        return testedEqtls;

    }

    private static void loadCohortEqtlgeneCounts(File featureMetadataFile, String cohort, DoubleMatrixDataset<String, String> covariateNonMissingCount) throws IOException {

        final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
        final CSVReader reader;
        if (featureMetadataFile.getName().endsWith(".gz")) {
            reader = new CSVReaderBuilder((new BufferedReader(new InputStreamReader((new GZIPInputStream(new FileInputStream(featureMetadataFile))))))).withCSVParser(parser).build();
        } else {
            reader = new CSVReaderBuilder((new BufferedReader(new InputStreamReader((new FileInputStream(featureMetadataFile)))))).withCSVParser(parser).build();
        }

        //First nextLine contains the header
        String[] nextLine = reader.readNext();

        while ((nextLine = reader.readNext()) != null) {

            if (covariateNonMissingCount.containsRow(nextLine[0])) {
                covariateNonMissingCount.setElement(nextLine[0], cohort, Double.parseDouble(nextLine[6]));
            }

        }


    }


    private static class Eqtl {

        private final String snp;
        private final String gene;

        public Eqtl(String snp, String gene) {
            this.snp = snp;
            this.gene = gene;
        }

        public String getSnp() {
            return snp;
        }

        public String getGene() {
            return gene;
        }

        public String getSnpGene() {
            return snp + "-" + gene;
        }
    }

}