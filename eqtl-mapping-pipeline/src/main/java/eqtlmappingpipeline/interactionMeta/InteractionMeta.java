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
import java.util.zip.GZIPInputStream;

public class InteractionMeta {

    public static void runMeta() throws IOException {

        File resultFolder = new File("/groups/umcg-fg/tmp04/projects/eqtlgen-phase2/interactions/biosInteractionsPcCor50");
        String[] cohorts = new String[]{"LL", "LLS_660Q", "LLS_OmniExpr", "NTR_Affy6", "NTR_GoNL", "RS"};
        ArrayList<File> cohortFolders = new ArrayList<File>();

        for (String cohort : cohorts) {
            cohortFolders.add(new File(resultFolder, cohort));
        }

        File allTestEffectsFile = new File("/groups/umcg-fg/tmp04/projects/eqtlgen-phase2/interactions/eQtlgenLeadVariants.txt.gz");

        File metaZFile = new File(resultFolder, "metaZ.txt");
        File metaSampleCountFile = new File(resultFolder, "metaSampleCount.txt");

        runMeta(cohortFolders, allTestEffectsFile, metaZFile, metaSampleCountFile);

    }


    private static void runMeta(List<File> cohortResultFolders, File allTestEffectsFile, File metaZFile, File metaSampleCountFile) throws IOException {

        List<Eqtl> testedEqtls = loadTestedEqtls(allTestEffectsFile);

        ArrayList<String> testedEqtlsNames = new ArrayList<>(testedEqtls.size());
        HashSet<String> eqtlGenes = new LinkedHashSet<>();
        for(Eqtl eqtl : testedEqtls) {
            testedEqtlsNames.add(eqtl.getSnpGene());
            eqtlGenes.add(eqtl.getGene());
        }
        System.out.println("Total eqtls: " + testedEqtlsNames.size());
        System.out.println("Unique eqtl genes: " + eqtlGenes.size());


        HashMap<String, DoubleMatrixDatasetRowCompressedReader> cohortZscoresMap = new HashMap<>(cohortResultFolders.size());
        HashMap<String, Map<String, Integer>> cohortColumnMaps = new HashMap<>(cohortResultFolders.size());


        HashSet<String> foundCovariates = new LinkedHashSet<>();

        for (File cohortResultFolder : cohortResultFolders) {

            DoubleMatrixDatasetRowCompressedReader cohortZscores = new DoubleMatrixDatasetRowCompressedReader(cohortResultFolder.getPath() + "/interactionZscoreTest.datg");
            cohortZscoresMap.put(cohortResultFolder.getName(), cohortZscores);
            cohortColumnMaps.put(cohortResultFolder.getName(), cohortZscores.getColumnMap());

            foundCovariates.addAll(cohortZscores.getColumnIdentifiers());

        }

        DoubleMatrixDataset<String, String> eqtlgeneNonMissingCount = new DoubleMatrixDataset<>(eqtlGenes, cohortZscoresMap.keySet());

        for (File cohortResultFolder : cohortResultFolders) {
            System.out.println(cohortResultFolder.getName());
            loadCohortEqtlgeneCounts(new File(cohortResultFolder.getPath(), "/feature_metadata.txt.gz"), cohortResultFolder.getName(), eqtlgeneNonMissingCount);

        }



        DoubleMatrixDataset<String, String> metaZscores = new DoubleMatrixDataset<>(testedEqtlsNames,foundCovariates);
        DoubleMatrixDataset<String, String> totalSamples = new DoubleMatrixDataset<>(testedEqtlsNames,foundCovariates);

        DoubleMatrixDataset<String, String> intermediates = new DoubleMatrixDataset<>(Arrays.asList("sumWZ", "sumW2"),foundCovariates);

        ProgressBar pb = new ProgressBar("Meta analysis", testedEqtlsNames.size(), ProgressBarStyle.COLORFUL_UNICODE_BLOCK);

       // int cohortWithAnyData;
        for(Eqtl eqtl : testedEqtls) {

            //cohortWithAnyData = 0;

            //reset intermediates
            intermediates.getMatrix().assign(0);

            for(Map.Entry<String, DoubleMatrixDatasetRowCompressedReader> cohortZscoresEntry : cohortZscoresMap.entrySet()){

                DoubleMatrixDatasetRowCompressedReader cohortZscores = cohortZscoresEntry.getValue();

//                boolean anyHasZ = false;

                if(cohortZscores.hasRow(eqtl.getSnpGene())) {

                    String cohort = cohortZscoresEntry.getKey();
                    Map<String, Integer> cohortColumnMap = cohortColumnMaps.get(cohort);

                    double[] cohortZscoresArray = cohortZscores.loadSingleRow(eqtl.getSnpGene());

                    double sampleSize = eqtlgeneNonMissingCount.getElement(eqtl.getGene(), cohort);

                    for(String covariate : foundCovariates){

                        if(cohortColumnMap.containsKey(covariate)){

                            double zScore = cohortZscoresArray[cohortColumnMap.get(covariate)];


                            if(!Double.isNaN(zScore)){

//                                if(!anyHasZ){
//                                    //cohortWithAnyData++;
//                                    anyHasZ=true;
//                                }




                                int intermediateColumn = intermediates.getColIndex(covariate);

                                intermediates.setElementQuick(0, intermediateColumn, intermediates.getElementQuick(0, intermediateColumn) + (sampleSize * zScore));
                                intermediates.setElementQuick(1, intermediateColumn, intermediates.getElementQuick(1, intermediateColumn) + (sampleSize * sampleSize));

                                totalSamples.setElement(eqtl.getSnpGene(), covariate, totalSamples.getElement(eqtl.getSnpGene(), covariate) + sampleSize);
                            }

                        }

                    }

                }

            }

            //Here we have looped through all cohorts

            //System.out.println(eqtl.getSnpGene() + "\t" + cohortWithAnyData);

            intermediates.viewDice().save("/groups/umcg-fg/tmp04/projects/eqtlgen-phase2/interactions/biosInteractions/intermediates.txt");



            int eqtlI = metaZscores.getRowIndex(eqtl.getSnpGene());
            for(int covariateI = 0; covariateI < foundCovariates.size(); covariateI++){
                double sumWZ = intermediates.getElementQuick(0,covariateI);
                double sumW2 = intermediates.getElementQuick(1,covariateI);

                if(sumW2 > 0){

                    metaZscores.setElementQuick(eqtlI, covariateI, sumWZ / Math.sqrt(sumW2));

                }

            }

            pb.step();

        }

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

            covariateNonMissingCount.setElement(nextLine[0], cohort, Double.parseDouble(nextLine[6]));


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