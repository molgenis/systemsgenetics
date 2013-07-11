package umcg.genetica.methylation;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Pattern;
import org.apache.commons.collections.primitives.ArrayDoubleList;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix.DoubleMatrixDataset;
import umcg.genetica.math.stats.Correlation;
import umcg.genetica.math.stats.TTest;
import umcg.genetica.math.stats.ZScores;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author MarcJan & Juha
 */
public class AssociatingPcasWithAnnotation {

    private static Pattern SPLIT_ON_TAB = Pattern.compile("\\t");

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws IOException, ClassNotFoundException {

//        String fileWithAnnotation = "/Data/MJ/Annotation/GPL8490_family_annotation_mesh_2013_2.txt";
        String fileWithAnnotation = "/Data/Sasha/GPL96GPL570AgeSamplesWithRangesAveragedInfantsLeftOut.txt";
//        String eigenVectorFile = "/Data/MJ/PCA_GPL8490_19102012/eigenvectors_Filtered.txt";
//        String eigenVectorFile = "/Data/MJ/PCA_GPL8490_SexChrs-Filtered/eigenvectors_Filtered.txt";
//        String eigenVectorFile = "/Data/MJ/PCA_GPL8490_19102012/GPL8490_family_all.quantilenormalized-missingvaluesreplaced.txt";
//        String eigenVectorFile = "/Data/MJ/GPL8490_family_SexProbesRemoved.quantilenormalized.missingvaluesreplaced-transposed.binary";
        String datafile = "/Data/GeneExpressionFinal/PCA/GPL570/GPL570ExpressiondataQNOnlyHumanSamplesOnlyENSGsCollapsed.binary";
//        String datafile = "/Data/GeneExpressionFinal/PCA/GPL96/GPL96ExpressiondataQNOnlyHumanSamplesOnlyENSGsCollapsed.binary";

        System.out.print("Read annotation file .... ");
        HashMap<String, SoftfileAnnotation> sampleAnnotation = readAnnotationFile(fileWithAnnotation);
        System.out.println("done");

//        TextFile tf = new TextFile("/Data/MJ/Top500AgeProbes.txt", TextFile.R);
//        Set<String> probes = new HashSet<String>(tf.readAsArrayList());
//        System.out.println(probes.size() + " probes read");

        TextFile tf = new TextFile("/Data/GeneExpressionFinal/SampleAnnotation/GPL570/GPL570CellLineSamplesAsPerTextMiningAndCorrelationWithCellLineProfile.txt", TextFile.R);
//        TextFile tf = new TextFile("/Data/GeneExpressionFinal/SampleAnnotation/GPL96/GPL96CellLineSamplesAsPerTextMiningAndCorrelationWithCellLineProfile.txt", TextFile.R);
        ArrayList<String> kickOutSamples = new ArrayList<String>(tf.readAsArrayList());
        System.out.println(kickOutSamples.size() + " samples will be kicked out");

        System.out.print("Read data file .... ");
//        DoubleMatrixDataset<String, String> data = readDoubleMatrixFile(datafile);
        DoubleMatrixDataset<String, String> data = readDoubleMatrixFileWithOutGivenColumns(datafile, kickOutSamples);
//        eigenVectors = eigenVectors.getTransposedDataset();
//        eigenVectors.save("/Data/MJ/GPL8490_family_SexProbesRemoved.quantilenormalized.missingvaluesreplaced-transposed.binary");
//        data.save("/Data/GeneExpressionFinal/PCA/GPL96/GPL96ExpressiondataQNOnlyHumanSamplesOnlyENSGsCollapsed.binary");
        System.out.println("done");

//        String infoKey = "Gender";
        String infoKey = "Age";
        ArrayList<String> entries = new ArrayList<String>();

        //entries.addAll(Arrays.asList("Male", "Female"));

        HashMap<String, HashMap<String, String>> interestSets;
        interestSets = selectSamplesWithInformationOfInterest(sampleAnnotation, infoKey, data, false);
        //interestSets = selectSamplesWithSeriesInformation(sampleAnnotation, eigenVectors);

        System.out.println("Number of interest sets: " + interestSets.size());

        //associateScoreAndItemOfInterest(eigenVectors, interestSets, entries);
        correlateScoreAndItemOfInterest(data, interestSets, "/Data/Sasha/GenesCorrelatedWithAgeGPL570CellLinesExcluded.txt", false);
    }

    /**
     * Read annotation file Tab separated file containing sample annotation
     * 
     * @param fileWithAnnotation
     * @return Sample annotation
     */
    private static HashMap<String, SoftfileAnnotation> readAnnotationFile(String fileWithAnnotation) throws IOException {

        TextFile tf = new TextFile("/Data/GeneExpressionFinal/SampleAnnotation/GSMToGenericGSEName-GSE2109SplitPerTissue.txt", false);
        Map<String, String> gsm2gse = tf.readAsHashMap(0, 1);

        HashMap<String, SoftfileAnnotation> sampleInfo = new HashMap<String, SoftfileAnnotation>();

        try {
            TextFile in = new TextFile(fileWithAnnotation, TextFile.R);

            String str = in.readLine();

            String[] headers = SPLIT_ON_TAB.split(str);


            int meshInfoIndex = -1;

            for (int i = 1; i < headers.length; ++i) {
                if (headers[i].toLowerCase().contains("mesh")) {
                    meshInfoIndex = i;
                    break;
                }
            }

            while ((str = in.readLine()) != null) {
                String[] entries = SPLIT_ON_TAB.split(str);
                String gse = gsm2gse.get(entries[0]);
                if (gse == null) {
                    System.out.println("problem");
                }
                entries[2] = gse;
                SoftfileAnnotation tmp = new SoftfileAnnotation();

                tmp.setAccession(entries[0]);

                if (!(meshInfoIndex < 0)) {
                    tmp.setMeshTerms(entries[meshInfoIndex]);
                }

                for (int i = 1; i < entries.length; ++i) {

                    tmp.putAnnotationInformation(headers[i], entries[i]);
                }

                sampleInfo.put(entries[0], tmp);
            }
            in.close();
        } catch (IOException e) {
            System.out.println(e.getMessage());
            System.exit(-1);
        }

        return (sampleInfo);
    }

    /**
     * Read double matrix file
     * Eigenvector file / pc file / probe matrix
     * @param eigenVectorFile
     * @return 
     */
    private static DoubleMatrixDataset<String, String> readDoubleMatrixFile(String eigenVectorFile) {

        return readDoubleMatrixFile(eigenVectorFile, null);
    }

    /**
     * Read double matrix file restricting to given rows
     * Eigenvector file / pc file / probe matrix
     * @param eigenVectorFile
     * @return 
     */
    private static DoubleMatrixDataset<String, String> readDoubleMatrixFile(String eigenVectorFile, Set<String> rowsToInclude) {

        DoubleMatrixDataset<String, String> tmp = new DoubleMatrixDataset<String, String>();
        try {
            if (rowsToInclude == null) {
                tmp = new DoubleMatrixDataset<String, String>(eigenVectorFile);//, "\t");                
            } else {
                tmp = new DoubleMatrixDataset<String, String>(eigenVectorFile, null, rowsToInclude);//, "\t");
            }
        } catch (IOException ex) {
            Logger.getLogger(AssociatingPcasWithAnnotation.class.getName()).log(Level.SEVERE, null, ex);
        }

        return (tmp);
    }

    /**
     * Read double matrix file not including given columns
     * Eigenvector file / pc file / probe matrix
     * @param eigenVectorFile
     * @return 
     */
    private static DoubleMatrixDataset<String, String> readDoubleMatrixFileWithOutGivenColumns(String eigenVectorFile, ArrayList<String> columnsToExclude) throws IOException, ClassNotFoundException {

        List<Object> columnObjectsOnly = DoubleMatrixDataset.getColumnObjectsOnly(eigenVectorFile);
        Set colsToRetain = new HashSet(columnObjectsOnly);
        colsToRetain.removeAll(columnsToExclude);
        DoubleMatrixDataset<String, String> tmp = new DoubleMatrixDataset<String, String>(eigenVectorFile, null, colsToRetain);

        return (tmp);
    }

    /**
     * Filter out interest GSE sets.
     * Sets need to have at least 2 different values for the infoKey of interest.
     * Automagicaly it checks if samples are in the double matrix dataset
     * 
     * @param sampleAnnotation Annotation information
     * @param infoKey Key for annotation of interest
     * @param doubleMatrix double matrix dataset
     * @return 
     */
    private static HashMap<String, HashMap<String, String>> selectSamplesWithInformationOfInterest(HashMap<String, SoftfileAnnotation> sampleAnnotation, String infoKey, DoubleMatrixDataset<String, String> eigenVectors, boolean samplesOnRows) {
        HashMap<String, HashMap<String, String>> gseSets = new HashMap<String, HashMap<String, String>>();

        ArrayList<String> removeSamples = new ArrayList<String>();

        for (Entry<String, SoftfileAnnotation> tmp : sampleAnnotation.entrySet()) {
            if (!(tmp.getValue().getAnnotationInformation().containsKey(infoKey))) {
                System.out.print("No " + infoKey + " information");
                System.exit(0);
            }
            break;
        }

        for (Entry<String, SoftfileAnnotation> sample : sampleAnnotation.entrySet()) {
            if (!(sample.getValue().getAnnotationInformation().get(infoKey).isEmpty()) || !(sample.getValue().getAnnotationInformation().get(infoKey).equals(""))) {
                boolean contains;
                if (samplesOnRows) {
                    contains = eigenVectors.rowObjects.contains(sample.getKey());
                } else {
                    contains = eigenVectors.colObjects.contains(sample.getKey());
                }
                if (contains) {
                    String seriesId = sample.getValue().getAnnotationInformation().get("series id");
                    if (gseSets.containsKey(seriesId)) {
                        gseSets.get(seriesId).put(sample.getKey(), sample.getValue().getAnnotationInformation().get(infoKey));
                    } else {
                        HashMap<String, String> tmp = new HashMap<String, String>();
                        tmp.put(sample.getKey(), sample.getValue().getAnnotationInformation().get(infoKey));
                        gseSets.put(seriesId, tmp);
                    }
                } else {
                    removeSamples.add(sample.getKey());
                }
            } else {
                removeSamples.add(sample.getKey());
            }
        }

        ArrayList<String> removeGseSets = new ArrayList<String>();

        int numberOfInterestSets = 0;
        int numberOfInterestSamples = 0;
        if (gseSets.size() > 0) {
            for (Entry<String, HashMap<String, String>> gse : gseSets.entrySet()) {

                ArrayList<String> uniqueValues = new ArrayList<String>();
                for (Entry<String, String> sample : gse.getValue().entrySet()) {
                    if (!uniqueValues.contains(sample.getValue())) {
                        uniqueValues.add(sample.getValue());
                    }
                }

                if (uniqueValues.size() >= 2 && gse.getValue().size() >= 10) {
                    numberOfInterestSets++;
                    // System.out.println(gse.getKey());
                    for (Entry<String, String> sample : gse.getValue().entrySet()) {
                        numberOfInterestSamples++;
                        //System.out.println("\t" + sample.getKey() + "\t" + sample.getValue());
                    }
                } else {
                    removeGseSets.add(gse.getKey());
                    for (Entry<String, String> sample : gse.getValue().entrySet()) {
                        removeSamples.add(sample.getKey());
                    }
                }
            }

        } else {
            System.out.println("Unforeseen error check Key and code");
            System.exit(0);
        }



        System.out.println("Number of sets: " + numberOfInterestSets);
        System.out.println("Total samples of interest: " + numberOfInterestSamples);

        for (String removeEntry : removeGseSets) {
            gseSets.remove(removeEntry);
        }

        for (String removeEntry : removeSamples) {
            sampleAnnotation.remove(removeEntry);
        }

        return (gseSets);

    }

    /**
     * Filter out interest GSE sets. Samples need to have a GSE id. Automagicaly
     * it checks if samples are in the double matrix dataset
     *
     * @param sampleAnnotation Annotation information
     * @param infoKey Key for annotation of interest
     * @param doubleMatrix double matrix dataset
     * @return
     */
    private static HashMap<String, HashMap<String, String>> selectSamplesWithSeriesInformation(HashMap<String, SoftfileAnnotation> sampleAnnotation, DoubleMatrixDataset<String, String> eigenVectors) {
        HashMap<String, HashMap<String, String>> gseSets = new HashMap<String, HashMap<String, String>>();

        ArrayList<String> removeSamples = new ArrayList<String>();

        for (Entry<String, SoftfileAnnotation> sample : sampleAnnotation.entrySet()) {
            if (!(sample.getValue().getAnnotationInformation().get("series id").isEmpty()) || !(sample.getValue().getAnnotationInformation().get("series id").equals(""))) {
                if (eigenVectors.rowObjects.contains(sample.getKey())) {
                    String seriesId = sample.getValue().getAnnotationInformation().get("series id");
                    if (gseSets.containsKey(seriesId)) {
                        gseSets.get(seriesId).put(sample.getKey(), sample.getValue().getAnnotationInformation().get("series id"));
                    } else {
                        HashMap<String, String> tmp = new HashMap<String, String>();
                        tmp.put(sample.getKey(), sample.getValue().getAnnotationInformation().get("series id"));
                        gseSets.put(seriesId, tmp);
                    }
                }
            } else {
                removeSamples.add(sample.getKey());
            }

        }

        for (String removeEntry : removeSamples) {
            sampleAnnotation.remove(removeEntry);
        }

        return (gseSets);

    }

    /**
     * Test for difference between 2 groups with T-test.
     * 
     * @param doubleMatrix 
     * @param interestSets
     * @param entries names of the two groups
     */
    public static void associateScoreAndItemOfInterest(DoubleMatrixDataset<String, String> doubleMatrix, HashMap<String, HashMap<String, String>> interestSets, ArrayList<String> entries) {
        //System.out.println(eigenVectors.nrCols);
        //System.out.println(eigenVectors.nrRows);
        HashMap<String, Double> scorePerGse = new HashMap<String, Double>();
        HashMap<String, Integer> indeces = new HashMap<String, Integer>();

        for (Entry<String, HashMap<String, String>> set : interestSets.entrySet()) {
            for (Entry<String, String> sample : set.getValue().entrySet()) {
                if (doubleMatrix.rowObjects.contains(sample.getKey())) {
                    int index = doubleMatrix.rowObjects.indexOf(sample.getKey());
                    indeces.put(sample.getKey(), index);
                } else {
                    System.out.println("Potential mismatch between annotation and samples");
                    System.out.println(sample.getKey() + " is not in value matrix");
                    System.out.println("\n However :" + indeces.size() + " are in the matrix");
                    System.exit(0);
                }
            }
        }

        //System.out.println(indeces.size());

        for (int i = 0; i < doubleMatrix.nrCols; ++i) {
            //System.out.println(doubleMatrix.colObjects.get(i));
            for (Entry<String, HashMap<String, String>> set : interestSets.entrySet()) {

                ArrayDoubleList valueSet1 = new ArrayDoubleList();
                ArrayDoubleList valueSet2 = new ArrayDoubleList();

                //ArrayList<String> keySet1 = new ArrayList<String>();
                //ArrayList<String> keySet2 = new ArrayList<String>();

                for (Entry<String, String> sample : set.getValue().entrySet()) {
                    if (sample.getValue().equals(entries.get(0))) {
                        //keySet1.add(sample.getKey());
                        valueSet1.add(doubleMatrix.rawData[indeces.get(sample.getKey())][i]);
                    } else if (sample.getValue().equals(entries.get(1))) {
                        //keySet2.add(sample.getKey());
                        valueSet2.add(doubleMatrix.rawData[indeces.get(sample.getKey())][i]);
                    }
                }
                double[] set1 = valueSet1.toArray(new double[0]);
                double[] set2 = valueSet2.toArray(new double[0]);

                if (set1.length > 2 && set2.length > 2) {
                    double zScore = TTest.testZscore(set1, set2);
                    //System.out.println(doubleMatrix.colObjects.get(i)+"_"+set.getKey()+"\t"+pValue);
                    scorePerGse.put(doubleMatrix.colObjects.get(i) + "_" + set.getKey(), zScore);
                }
            }
        }
    }

    public static String[] readGeneNamesForProbes(String filename, List<String> probes) throws IOException {
        String[] geneNames = new String[probes.size()];
        TextFile in = new TextFile(filename, TextFile.R);
        String line = in.readLine();
        while ((line = in.readLine()) != null) {
            String[] split = SPLIT_ON_TAB.split(line);
            int indexOf = probes.indexOf(split[0]);
            if (indexOf >= 0) {
                geneNames[indexOf] = split[1];
            }
        }
        in.close();
        return geneNames;
    }

    /**
     * Correlate values of interest to age
     * 
     * @param doubleMatrix 
     * @param interestSets
     */
    public static void correlateScoreAndItemOfInterest(DoubleMatrixDataset<String, String> doubleMatrix, HashMap<String, HashMap<String, String>> interestSets, String outfile, boolean samplesOnRows) throws IOException {
        HashMap<String, Integer> indeces = new HashMap<String, Integer>();
        int largestSet = 0;
        for (Entry<String, HashMap<String, String>> set : interestSets.entrySet()) {
            for (Entry<String, String> sample : set.getValue().entrySet()) {
                int index;
                if (samplesOnRows) {
                    index = doubleMatrix.rowObjects.indexOf(sample.getKey());
                } else {
                    index = doubleMatrix.colObjects.indexOf(sample.getKey());
                }
                if (index > -1) {
                    indeces.put(sample.getKey(), index);
                } else {
                    System.out.println("Potential mismatch between annotation and samples");
                    System.out.println(sample.getKey() + " is not in value matrix");
                    System.out.println("\n However :" + indeces.size() + " are in the matrix");
                    System.exit(0);
                }
            }
            if (largestSet < set.getValue().size()) {
                largestSet = set.getValue().size();
            }
        }
        Correlation.correlationToZScore(largestSet);
//        String[] geneNames = readGeneNamesForProbes("/Data/MJ/gpl_8490.chip", doubleMatrix.colObjects);
        TextFile plos = new TextFile("/Data/MJ/Epigenome-Wide_Scans.txt", TextFile.R);
        Map<String, String> plosPValues = plos.readAsHashMap(0, 5);

        int nrProbes;
        if (samplesOnRows) {
            nrProbes = doubleMatrix.nrCols;
        } else {
            nrProbes = doubleMatrix.nrRows;
        }
        double[] metaZ = new double[nrProbes];
        double[][] probeGSEZ = new double[nrProbes][interestSets.size()];
        double[][] leaveOneOutProbeGSEZ = new double[nrProbes][interestSets.size() + 1];
        SpearmansCorrelation sc = new SpearmansCorrelation();
        TextFile out = new TextFile(outfile, TextFile.W);
        String[] setNames = new String[interestSets.size()];
        int[] setSizes = new int[interestSets.size()];
        for (int i = 0; i < nrProbes; ++i) {
            double[] zScores = new double[interestSets.size()];
            int index = 0;
            for (Entry<String, HashMap<String, String>> set : interestSets.entrySet()) {
                int sizeOfGseSet = set.getValue().size();
                setSizes[index] = sizeOfGseSet;
                setNames[index] = set.getKey();
                ArrayDoubleList valueSet = new ArrayDoubleList();
                ArrayDoubleList ageSet = new ArrayDoubleList();

                for (Entry<String, String> sample : set.getValue().entrySet()) {
                    if (samplesOnRows) {
                        valueSet.add(doubleMatrix.rawData[indeces.get(sample.getKey())][i]);
                    } else {
                        valueSet.add(doubleMatrix.rawData[i][indeces.get(sample.getKey())]);
                    }
                    try {
                        ageSet.add(Double.parseDouble(sample.getValue()));
                    } catch (NumberFormatException ex) { // data not numerical, assume gender here as a quick fix
                        ageSet.add("male".equals(sample.getValue().toLowerCase()) ? 1 : 2);
                    }
                }
                double[] setValues = valueSet.toArray(new double[0]);
                double[] setAges = ageSet.toArray(new double[0]);

                double spearman = sc.correlation(setValues, setAges);
//                    double correlation = JSci.maths.ArrayMath.correlation(setValues, setAges);
//                    double zScore = Correlation.convertCorrelationToZScore(sizeOfGseSet, correlation);
                double zScore = Correlation.convertCorrelationToZScore(sizeOfGseSet, spearman);
                zScores[index] = zScore;
                index++;
            }
            probeGSEZ[i] = zScores;

            // leave-one-out z weighting
            for (int leave = 0; leave < zScores.length; leave++) { // leave leave'th z score out
                double[] zScoresLeft = new double[zScores.length - 1];
                int[] setSizesLeft = new int[zScores.length - 1];
                int zi = 0;
                for (int j = 0; j < zScores.length; j++) {
                    if (j != leave) {
                        zScoresLeft[zi] = zScores[j];
                        setSizesLeft[zi] = setSizes[j];
                        zi++;
                    }
                }
                double leftZ = ZScores.getWeightedZ(zScoresLeft, setSizesLeft);
//                double p = ZScores.zToP(leftZ);
                leaveOneOutProbeGSEZ[i][leave + 1] = leftZ;
//                out.writeln(doubleMatrix.colObjects.get(i) + "\t" + geneNames[i] + "\t" + setNames[leave] + "\t" + leftZ + "\t" + p);
            }

            metaZ[i] = ZScores.getWeightedZ(zScores, setSizes);
            leaveOneOutProbeGSEZ[i][0] = metaZ[i];
            double p = ZScores.zToP(metaZ[i]);
//            String plosP = plosPValues.get(doubleMatrix.colObjects.get(i));
//            out.writeln(doubleMatrix.colObjects.get(i) + "\t" + geneNames[i] + "\t-\t" + metaZ[i] + "\t" + p);
            if (samplesOnRows) {
                out.writeln(doubleMatrix.colObjects.get(i) + "\t-\t" + metaZ[i] + "\t" + p);
            } else {
                out.writeln(doubleMatrix.rowObjects.get(i) + "\t-\t" + metaZ[i] + "\t" + p);
            }
        }
        out.close();

        for (int i = 0; i < setNames.length; i++) {
            System.out.println(setNames[i] + "\t" + setSizes[i]);
        }

//        DoubleMatrixDataset<String, String> probeGSEDataset = new DoubleMatrixDataset<String, String>(probeGSEZ);
//        probeGSEDataset.rowObjects = doubleMatrix.colObjects;
//        probeGSEDataset.colObjects = Arrays.asList(setNames);
//        probeGSEDataset.removeColumnsWithNaNs();
//        probeGSEDataset.save("/Data/MJ/ProbeGSEAgeCorrelationZScores.txt");
//
//        DoubleMatrixDataset<String, String> leaveOneOutDataset = new DoubleMatrixDataset<String, String>(leaveOneOutProbeGSEZ);
//        leaveOneOutDataset.rowObjects = doubleMatrix.colObjects;
//        List<String> setNamesAbstractList = Arrays.asList(setNames);
//        List<String> setNamesList = new ArrayList<String>(setNamesAbstractList);
//        setNamesList.add(0, "meta");
//        leaveOneOutDataset.colObjects = setNamesList;
//        leaveOneOutDataset.removeColumnsWithNaNs();
//        leaveOneOutDataset.save("/Data/MJ/ProbeGSEAgeCorrelationLeaveOneOutZScores.txt");

//        for (Entry<String, HashMap<String, String>> set : interestSets.entrySet()) {
//            for (Entry<String, String> e : set.getValue().entrySet()) {
//                Integer sampleIndex = indeces.get(e.getKey());
//                double correlation = ArrayMath.correlation(doubleMatrix.rawData[sampleIndex], metaZ);
//                System.out.println(e.getKey() + "\t" + e.getValue() + "\t" + correlation);
//            }
//        }
    }

    /**
     * Test for difference between all GSE groups with Anova.
     *
     * @param doubleMatrix
     * @param interestSets
     * @param entries names of the two groups
     */
    public static void associateAnovaScoreAndItemOfInterest(DoubleMatrixDataset<String, String> doubleMatrix, HashMap<String, HashMap<String, String>> interestSets) {
        //System.out.println(eigenVectors.nrCols);
        //System.out.println(eigenVectors.nrRows);
        HashMap<String, Double> scorePerGse = new HashMap<String, Double>();
        HashMap<String, Integer> indeces = new HashMap<String, Integer>();

        for (Entry<String, HashMap<String, String>> set : interestSets.entrySet()) {
            for (Entry<String, String> sample : set.getValue().entrySet()) {
                if (doubleMatrix.rowObjects.contains(sample.getKey())) {
                    int index = doubleMatrix.rowObjects.indexOf(sample.getKey());
                    indeces.put(sample.getKey(), index);
                } else {
                    System.out.println("Potential mismatch between annotation and samples");
                    System.out.println(sample.getKey() + " is not in value matrix");
                    System.out.println("\n However :" + indeces.size() + " are in the matrix");
                    System.exit(0);
                }
            }
        }
    }
}
