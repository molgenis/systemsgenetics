/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.bondermj.pcoa;

import cern.colt.matrix.tdouble.algo.DenseDoubleAlgebra;
import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix2D;
//import cern.colt.matrix.tdouble.algo.decomposition.FastSymetricDenseDoubleEigenvalueDecomposition;
import cern.colt.matrix.tdouble.algo.decomposition.DenseDoubleEigenvalueDecomposition;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix1D;
import cern.colt.matrix.tdouble.impl.DenseLargeDoubleMatrix2D;
import cern.jet.math.tdouble.DoubleFunctions;
import edu.emory.mathcs.utils.ConcurrencyUtils;
import java.io.File;
import java.io.IOException;
import java.text.DateFormat;
import java.text.NumberFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.LinkedHashMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import nl.umcg.bondermj.pcoa.PcoaParamaters.MatrixType;
import nl.umcg.bondermj.pcoa.PcoaParamaters.PackageToUse;
import org.apache.commons.cli.ParseException;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.matrix2.MatrixTools;
import umcg.genetica.math.stats.concurrent.*;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.EVD;
import no.uib.cipr.matrix.NotConvergedException;
import org.apache.commons.math3.stat.ranking.NaNStrategy;
import org.apache.commons.math3.stat.ranking.NaturalRanking;
import org.apache.commons.math3.stat.ranking.RankingAlgorithm;
import org.apache.commons.math3.stat.ranking.TiesStrategy;

/**
 *
 * @author MarcJan
 */
public class DoPcoa {

    private static final RankingAlgorithm RANKER_TIE = new NaturalRanking(NaNStrategy.FAILED, TiesStrategy.SEQUENTIAL);

    private static final String HEADER
            = "  /---------------------------------------\\\n"
            + "  |              Parallel PCoA            |\n"
            + "  |                                       |\n"
            + "  |             Marc Jan Bonder           |\n"
            + "  |            bondermj@gmail.com         |\n"
            + "  |                                       |\n"
            + "  |          Genetics department          |\n"
            + "  |  University Medical Center Groningen  |\n"
            + "  \\---------------------------------------/";
    private static final DateFormat DATE_TIME_FORMAT = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
    public static final NumberFormat DEFAULT_NUMBER_FORMATTER = NumberFormat.getInstance();

    /**
     * @param args the command line arguments
     */
    public static void main(String... args) throws IOException, Exception {

        System.out.println(HEADER);
        System.out.println();

        Date currentDataTime = new Date();
        System.out.println("Current date and time: " + DATE_TIME_FORMAT.format(currentDataTime));
        System.out.println();

        System.out.flush(); //flush to make sure header is before errors
        Thread.sleep(25); //Allows flush to complete

        if (args.length == 0) {
            PcoaParamaters.printHelp();
            System.exit(1);
        }

        final PcoaParamaters paramaters;
        try {
            paramaters = new PcoaParamaters(args);
        } catch (ParseException ex) {
            System.err.println("Invalid command line arguments: ");
            System.err.println(ex.getMessage());
            System.err.println();
            PcoaParamaters.printHelp();
            System.exit(1);
            return;
        }

        MatrixType matrixType = paramaters.getMatrixType();
        PackageToUse packageToUse = paramaters.getPackageToUse();

        boolean overRows = paramaters.isOverRows();
        String nrComp = paramaters.getNumberComponentsToCalc();
        String fileIn = paramaters.getFileIn();
        File fileInput = new File(fileIn);
        if (!fileInput.exists()) {
            System.out.println("Warning: input for input file not found.");
            System.exit(-1);
        }
        String prefix = fileIn.replace(".txt", "");
        prefix = prefix.replace(".gz", "");

//To add later?
//        String componentsToSkip = null;
//        boolean regresOutPcs = false;
//        int maxNrToRegres = -1;
//        int stepsize = -1;
//        int startSize=0;
//
        int usable = paramaters.getThreads();
        if (usable > Runtime.getRuntime().availableProcessors()) {
            int nrProc = Runtime.getRuntime().availableProcessors();
            usable = (int) Math.floor(nrProc - (nrProc / 10.0d));
        }
        ConcurrencyUtils.setNumberOfThreads(usable);

        DoubleMatrixDataset<String, String> inputData;

        System.out.println("Reading the input data ....\n");
        try {
            inputData = DoubleMatrixDataset.loadDoubleData(fileIn);
        } catch (IOException ex) {
            Logger.getLogger(DoPcoa.class.getName()).log(Level.SEVERE, null, ex);
            throw new IOException("Problem reading data.");
        } catch (Exception ex) {
            Logger.getLogger(DoPcoa.class.getName()).log(Level.SEVERE, null, ex);
            throw new IOException("Problem reading data.");
        }
      

        if (overRows) {
            inputData = inputData.viewDice();
            System.out.println("\nConducting Analysis over rows\n");
        } else {
            System.out.println("\nConducting Analysis over columns\n");
        }

//        if ((inputData.columns() * (long) inputData.columns()) > (Integer.MAX_VALUE - 2)) {
//            throw new Exception("Matrix has to many columns for PCA");
//        }
        
        if(paramaters.isRank()){
            MatrixTools.rankColumns(inputData.getMatrix());
        }
        
        if (paramaters.isCenter() && paramaters.isScale()) {
            MatrixTools.centerAndScaleColum(inputData.getMatrix());
        } else if (paramaters.isScale()) {
            MatrixTools.scaleColum(inputData.getMatrix());
        } else if (paramaters.isCenter()) {
            MatrixTools.centerColum(inputData.getMatrix());
        }

        DoubleMatrix2D matrix = null;

        try {
            switch (matrixType) {
                case COVARIATION: {
                    DoubleMatrixDataset<String, String> t;
                    File f = new File(prefix + ".CovariationMatrix.txt.gz");
                    if (f.exists()) {
                        System.out.println("Warning found an old distance file, taking this for the distances. ");
                        t = DoubleMatrixDataset.loadDoubleData(prefix + ".CorrelationMatrix.txt.gz");
                        if (inputData.getHashCols().size() == t.getHashCols().size()) {
                            int i = 0;
                            for (String s : inputData.getColObjects()) {
                                if (!s.equals(t.getColObjects().get(i))) {
                                    break;
                                }
                                i++;
                            }
                            if (i == (inputData.getHashCols().size() - 1)) {
                                matrix = t.getMatrix();
                                break;
                            }
                        }
                    }
                    ConcurrentCovariation c = new ConcurrentCovariation(usable);
                    matrix = c.pairwiseCovariationDoubleMatrix(inputData.getMatrix().viewDice().toArray());
                    t = new DoubleMatrixDataset<>(matrix, inputData.getHashCols(), inputData.getHashCols());
                    t.save(prefix + ".CovariationMatrix.txt.gz");

                    break;
                }
                case CITYBLOCK: {
                    DoubleMatrixDataset<String, String> t;
                    File f = new File(prefix + ".CityBlockMatrix.txt.gz");
                    if (f.exists()) {
                        System.out.println("Warning found an old distance file, taking this for the distances. ");
                        t = DoubleMatrixDataset.loadDoubleData(prefix + ".CityBlockMatrix.txt.gz");
                        if (inputData.getHashCols().size() == t.getHashCols().size()) {
                            int i = 0;
                            for (String s : inputData.getColObjects()) {
                                if (!s.equals(t.getColObjects().get(i))) {
                                    break;
                                }
                                i++;
                            }
                            if (i == (inputData.getHashCols().size() - 1)) {
                                matrix = t.getMatrix();
                                break;
                            }
                        }
                    }
                    ConcurrentCityBlock c = new ConcurrentCityBlock(usable);
                    matrix = c.pairwiseCityBlockDoubleMatrix(inputData.getMatrix().viewDice().toArray());
                    t = new DoubleMatrixDataset<>(matrix, inputData.getHashCols(), inputData.getHashCols());
                    t.save(prefix + ".CityBlockMatrix.txt.gz");
                    break;
                }
                case BRAYCURTIS: {
                    DoubleMatrixDataset<String, String> t;
                    File f = new File(prefix + ".BrayCurtisMatrix.txt.gz");
                    if (f.exists()) {
                        System.out.println("Warning found an old distance file, taking this for the distances. ");
                        t = DoubleMatrixDataset.loadDoubleData(prefix + ".BrayCurtisMatrix.txt.gz");
                        if (inputData.getHashCols().size() == t.getHashCols().size()) {
                            int i = 0;
                            for (String s : inputData.getColObjects()) {
                                if (!s.equals(t.getColObjects().get(i))) {
                                    break;
                                }
                                i++;
                            }
                            if (i == (inputData.getHashCols().size() - 1)) {
                                matrix = t.getMatrix();
                                break;
                            }
                        }
                    }
                    ConcurrentBrayCurtis c = new ConcurrentBrayCurtis(usable);
                    matrix = c.pairwiseBrayCurtisDoubleMatrix(inputData.getMatrix().viewDice().toArray());
                    t = new DoubleMatrixDataset<>(matrix, inputData.getHashCols(), inputData.getHashCols());
                    t.save(prefix + ".BrayCurtisMatrix.txt.gz");
                    break;
                }
                default: {
                    DoubleMatrixDataset<String, String> t;
                    File f = new File(prefix + ".CorrelationMatrix.txt.gz");
                    if (f.exists()) {
                        System.out.println("Warning found an old distance file, taking this for the distances. ");
                        t = DoubleMatrixDataset.loadDoubleData(prefix + ".CorrelationMatrix.txt.gz");
                        if (inputData.getHashCols().size() == t.getHashCols().size()) {
                            int i = 0;
                            for (String s : inputData.getColObjects()) {
                                if (!s.equals(t.getColObjects().get(i))) {
                                    break;
                                }
                                i++;
                            }
                            if (i == (inputData.getHashCols().size() - 1)) {
                                matrix = t.getMatrix();
                                break;
                            }
                        }
                    }
                    ConcurrentCorrelation c = new ConcurrentCorrelation(usable);
                    matrix = c.pairwiseCorrelationDoubleMatrix(inputData.getMatrix().viewDice().toArray());
                    t = new DoubleMatrixDataset<>(matrix, inputData.getHashCols(), inputData.getHashCols());
                    t.save(prefix + ".CorrelationMatrix.txt.gz");
                    break;
                }
            }

            System.out.println("done");

        } catch (IOException ex) {
            System.out.println("Failed to write distance matrix");
        }
        if (matrix == null) {
            System.out.println("Failed to read or make distance matrix");
            System.exit(0);
        }
        
        Integer nrComponents;
        if (isInteger(nrComp)) {
            nrComponents = Integer.parseInt(nrComp);
        } else if (nrComp.equals("all")) {
            nrComponents = inputData.columns();
        } else {
            nrComponents = inputData.columns();
        }

        if (nrComponents < 1) {
            System.out.println("Warning: Number of PCs to calculate should be at least 1");
            System.out.println("Reset nrComp to all");
            nrComponents = inputData.columns();
        } else if (nrComponents > inputData.columns()) {
            System.out.println("Warning: Number of PCs to calculate is larger dan dim of matrix");
            System.out.println("Reset nrComp to all");
            nrComponents = inputData.columns();
        }

        try {
            if (packageToUse.equals(PackageToUse.COLT)) {
                calculateColtPCA(inputData, matrix, prefix, nrComponents);
            } else {
                calculateMtjPCA(inputData, matrix, prefix, nrComponents);
            }

        } catch (IOException ex) {
            Logger.getLogger(DoPcoa.class.getName()).log(Level.SEVERE, null, ex);
        }
        matrix = null;

//        if(regresOutPcs){
//            if(nrComponents != inputData.columns()){
//                if(componentsToSkip!=null){
//                    throw new UnsupportedOperationException("Not supported yet.");
//                }
//                try {
//                    regressOutPCs(inputData, prefix, prefix + ".PCAOverSamplesPrincipalComponents.txt.gz" , prefix+ ".PCAOverSamplesEigenvectors.txt.gz", maxNrToRegres, stepsize);
//                } catch (IOException ex) {
//                    Logger.getLogger(DoPcoa.class.getName()).log(Level.SEVERE, null, ex);
//                }
//                
//            } else {
//                inputData=null;
//                if(componentsToSkip!=null){
//                    TextFile cmpSkipListFile = new TextFile(componentsToSkip, TextFile.R);
//                    ArrayList<String> compSkipList = cmpSkipListFile.readAsArrayList();
//                    try {
//                        generatePcCorrectedDataWithSkipping(prefix, prefix + ".PCAOverSamplesPrincipalComponents.txt.gz" , prefix+ ".PCAOverSamplesEigenvectorsTransposed.txt.gz", startSize, maxNrToRegres, stepsize, compSkipList);
//                    } catch (IOException ex) {
//                        Logger.getLogger(DoPcoa.class.getName()).log(Level.SEVERE, null, ex);
//                    }
//                } else {
//                    try {
//                        generatePcCorrectedData(prefix, prefix + ".PCAOverSamplesPrincipalComponents.txt.gz" , prefix+ ".PCAOverSamplesEigenvectorsTransposed.txt.gz", startSize, maxNrToRegres, stepsize);
//                    } catch (IOException ex) {
//                        Logger.getLogger(DoPcoa.class.getName()).log(Level.SEVERE, null, ex);
//                    }
//                }
//            } 
//        }
    }

    /**
     * Calculate Colt PCA scores Matrix contains the correlation, covariation,
     * bray curtis distance or city block distance.
     *
     * @param dataset
     * @param CorMatrix
     * @param fileNamePrefix
     * @param nrOfPCsToCalculate
     * @return
     * @throws IOException
     */
    public static void calculateColtPCA(DoubleMatrixDataset<String, String> dataset, DoubleMatrix2D CorMatrix, String fileNamePrefix, Integer nrOfPCsToCalculate) throws IOException {
        String expressionFile = fileNamePrefix;
        System.out.println("- Performing PCA over matrix of size: " + CorMatrix.columns() + "x" + CorMatrix.rows());

        DenseDoubleEigenvalueDecomposition eig = new DenseDoubleEigenvalueDecomposition(CorMatrix);

        //System.out.println(dataset.columns());
        System.out.println("- Number of components to be written: " + nrOfPCsToCalculate);

        DoubleMatrix1D eigenValues = eig.getRealEigenvalues();

        TextFile out = new TextFile(expressionFile + ".PCAOverSamplesEigenvalues_pc.txt.gz", TextFile.W);
        out.writeln("PCA\tEigenValue\tExplained variance\tTotal variance");

        double cumExpVarPCA = 0;

        double eigenValueSum = eigenValues.zSum();

        LinkedHashMap<String, Integer> tmpNameBuffer = new LinkedHashMap<>();
        for (int pca = 0; pca < nrOfPCsToCalculate; pca++) {
            double expVarPCA = eigenValues.getQuick((int) eigenValues.size() - 1 - pca) / eigenValueSum;
            int pcaNr = pca + 1;
            cumExpVarPCA += expVarPCA;
            out.write(pcaNr + "\t" + eigenValues.getQuick((int) eigenValues.size() - 1 - pca) + "\t" + expVarPCA + "\t" + cumExpVarPCA + "\n");
            tmpNameBuffer.put("Comp" + String.valueOf(pcaNr), pca);
        }
        out.close();

        DoubleMatrix2D eigenValueMatrix = eig.getV();
        DoubleMatrixDataset<String, String> datasetEV = new DoubleMatrixDataset<>((DenseDoubleMatrix2D) eigenValueMatrix.viewColumnFlip().viewPart(0, 0, dataset.columns(), nrOfPCsToCalculate), dataset.getHashCols(), tmpNameBuffer);
        eig = null;

        datasetEV.save(expressionFile + ".PCAOverSamplesEigenvectors_pc.txt.gz");
//        datasetEV.saveDice(expressionFile + ".PCAOverSamplesEigenvectorsTransposed_pc.txt.gz");

        System.out.println("Calculating PCs");
        System.out.println("Initializing PCA matrix");

        DoubleMatrixDataset<String, String> datasetPCAOverSamplesPCAs;
        if ((dataset.columns() * (long) dataset.rows()) < (Integer.MAX_VALUE - 2)) {
            datasetPCAOverSamplesPCAs = new DoubleMatrixDataset<String, String>((DenseDoubleMatrix2D) dataset.getMatrix().zMult(datasetEV.getMatrix(), null), dataset.getHashRows(), tmpNameBuffer);
        } else {
            datasetPCAOverSamplesPCAs = new DoubleMatrixDataset<String, String>((DenseLargeDoubleMatrix2D) dataset.getMatrix().zMult(datasetEV.getMatrix(), null), dataset.getHashRows(), tmpNameBuffer);
        }

        System.out.println("Saving PCA scores: " + expressionFile + ".PCAOverSamplesPrincipalComponents.txt.gz");
        datasetPCAOverSamplesPCAs.save(expressionFile + ".PCAOverSamplesPrincipalComponents.txt.gz");
    }

    /**
     * Calculate MTJ PCA scores Matrix contains the correlation, covariation,
     * bray curtis distance or city block distance.
     *
     * @param dataset
     * @param CorMatrix
     * @param fileNamePrefix
     * @param nrOfPCsToCalculate
     * @return
     * @throws IOException
     */
    public static void calculateMtjPCA(DoubleMatrixDataset<String, String> dataset, DoubleMatrix2D CorMatrix, String fileNamePrefix, Integer nrOfPCsToCalculate) throws IOException {
        String expressionFile = fileNamePrefix;
        System.out.println("- Performing PCA over matrix of size: " + CorMatrix.columns() + "x" + CorMatrix.rows());

        EVD evd = new EVD(CorMatrix.columns(), false, true);
        try {
            evd.factor(new DenseMatrix(CorMatrix.toArray()));
        } catch (NotConvergedException ex) {
            Logger.getLogger(DoPcoa.class.getName()).log(Level.SEVERE, null, ex);
        }

        //System.out.println(dataset.columns());
        System.out.println("- Number of components to be written: " + nrOfPCsToCalculate);

        DoubleMatrix1D eigenValues = new DenseDoubleMatrix1D(evd.getRealEigenvalues());

        TextFile out = new TextFile(expressionFile + ".PCAOverSamplesEigenvalues.txt.gz", TextFile.W);
        out.writeln("PCA\tEigenValue\tExplained variance\tTotal variance");

        double cumExpVarPCA = 0;

        DenseMatrix eigenVectorsMatrix = evd.getRightEigenvectors();

        //Fix sorting of eigenValues
        double[] eigenValuesRanked = RANKER_TIE.rank(eigenValues.toArray());

        DenseDoubleMatrix1D rankedeigenValues = new DenseDoubleMatrix1D((int) eigenValues.size());
        DenseMatrix rankedeigenVectorsMatrix = new DenseMatrix(eigenVectorsMatrix.numRows(), eigenVectorsMatrix.numColumns());

        //Rows are the samples, cols are the eigenvectors
        for (int i = 0; i < eigenValuesRanked.length; i++) {
            int correctColumn = eigenValuesRanked.length - ((int) eigenValuesRanked[i]);
            for (int r = 0; r < eigenVectorsMatrix.numRows(); r++) {
                rankedeigenVectorsMatrix.set(r, correctColumn, eigenVectorsMatrix.get(r, i));
            }
            rankedeigenValues.setQuick(correctColumn, eigenValues.get(i));
        }

        eigenValues = rankedeigenValues;
        eigenVectorsMatrix = rankedeigenVectorsMatrix;
        // End fix sort.

        double eigenValueSum = eigenValues.zSum();
        LinkedHashMap<String, Integer> tmpNameBuffer = new LinkedHashMap<>();
        for (int pca = 0; pca < nrOfPCsToCalculate; pca++) {
            double expVarPCA = eigenValues.getQuick(pca) / eigenValueSum;
            int pcaNr = pca + 1;
            cumExpVarPCA += expVarPCA;
            out.write(pcaNr + "\t" + eigenValues.getQuick(pca) + "\t" + expVarPCA + "\t" + cumExpVarPCA + "\n");
            tmpNameBuffer.put("Comp" + String.valueOf(pcaNr), pca);
        }
        out.close();

        DoubleMatrixDataset<String, String> datasetEV = new DoubleMatrixDataset<>(MatrixTools.toDenseDoubleMatrix(eigenVectorsMatrix).viewPart(0, 0, dataset.columns(), nrOfPCsToCalculate), dataset.getHashCols(), tmpNameBuffer);
        evd = null;

        datasetEV.save(expressionFile + ".PCAOverSamplesEigenvectors.txt.gz");

//        datasetEV.saveDice(expressionFile + ".PCAOverSamplesEigenvectorsTransposed.txt.gz");
        datasetEV = null;
        System.out.println("Calculating PCs");
        System.out.println("Initializing PCA matrix");

        DenseMatrix scoreMatrix = new DenseMatrix(dataset.rows(), dataset.columns());

        new DenseMatrix((dataset.getMatrix().toArray())).mult(eigenVectorsMatrix.transpose(), scoreMatrix);

        DoubleMatrixDataset<String, String> datasetPCAOverSamplesPCAs = new DoubleMatrixDataset<String, String>(MatrixTools.toDenseDoubleMatrix(scoreMatrix).viewPart(0, 0, dataset.rows(), nrOfPCsToCalculate), dataset.getHashRows(), tmpNameBuffer);
        System.out.println("Saving PCA scores: " + expressionFile + ".PCAOverSamplesPrincipalComponents.txt.gz");
        datasetPCAOverSamplesPCAs.save(expressionFile + ".PCAOverSamplesPrincipalComponents.txt.gz");
    }

    public static void regressOutPCs(DoubleMatrixDataset<String, String> dataset, String fileNamePrefix, String PcaFile, String eigenVectorFile, int nrPCAsOverSamplesToRemove, int nrIntermediatePCAsOverSamplesToRemoveToOutput) throws IOException {
        DenseDoubleAlgebra Alg = new DenseDoubleAlgebra();

        DoubleMatrixDataset<String, String> datasetPCAOverSamplesPCAs = DoubleMatrixDataset.loadDoubleData(PcaFile);
        DoubleMatrix2D datasetPCAOverSamplesPCAsSmal = Alg.subMatrix(datasetPCAOverSamplesPCAs.getMatrix(), 0, datasetPCAOverSamplesPCAs.rows() - 1, 0, nrPCAsOverSamplesToRemove - 1);
        datasetPCAOverSamplesPCAs = null;
        System.gc();
        System.gc();

        DoubleMatrixDataset<String, String> datasetEV = DoubleMatrixDataset.loadDoubleData(eigenVectorFile);
        DoubleMatrix2D datasetEVSmal = Alg.subMatrix(datasetEV.getMatrix(), 0, datasetEV.rows() - 1, 0, nrPCAsOverSamplesToRemove - 1);
        datasetEV = null;
        System.gc();
        System.gc();
        System.gc();

        System.out.println("\nInitializing residual gene expression matrix");

        for (int t = 0; t < nrPCAsOverSamplesToRemove; t++) {

            dataset.getMatrix().assign(Alg.multOuter(datasetPCAOverSamplesPCAsSmal.viewColumn(t), datasetEVSmal.viewColumn(t), null), DoubleFunctions.minus);

            int nrPCAs = t + 1;
            if (nrPCAs % nrIntermediatePCAsOverSamplesToRemoveToOutput == 0) {
                dataset.save(fileNamePrefix + "." + nrPCAs + "PCAsOverSamplesRemoved.txt.gz");
                System.out.println("Removed\t" + nrPCAs + "\tPCs. File:\t" + fileNamePrefix + "." + nrPCAs + "PCAsOverSamplesRemoved.txt.gz");
            }
        }

    }

    public static void generatePcCorrectedDataWithSkipping(String fileNamePrefix, String PcaFile, String eigenVectorFile, int startSize, int nrPCAsOverSamplesToRemove, int nrIntermediatePCAsOverSamplesToRemoveToOutput, ArrayList<String> compSkipList) throws IOException {

        DoubleMatrixDataset<String, String> datasetPCAOverSamplesPCAs = DoubleMatrixDataset.loadDoubleData(PcaFile);

        DoubleMatrixDataset<String, String> datasetEV = DoubleMatrixDataset.loadDoubleData(eigenVectorFile);

        DoubleMatrixDataset dataset = new DoubleMatrixDataset(null, datasetPCAOverSamplesPCAs.getHashRows(), datasetEV.getHashCols());

        String expressionFile = fileNamePrefix;
        System.out.println("\nInitializing residual gene expression matrix");

        DenseDoubleAlgebra Alg = new DenseDoubleAlgebra();

        //Clean and optimize for calculation.
        DoubleMatrix2D datasetPCAOverSamplesPCAs_Mat = datasetPCAOverSamplesPCAs.getMatrix();
        DoubleMatrix2D datasetEV_Mat = datasetEV.getMatrix();

        datasetPCAOverSamplesPCAs = null;
        datasetEV = null;

        for (int r = 0; r < (startSize - 1); r++) {
            if (!compSkipList.contains("Comp" + (r + 1))) {
                datasetEV_Mat.viewRow(r).assign(0);
            }
        }

        boolean dataChanged = false;
        for (int r = (startSize - 1); r < nrPCAsOverSamplesToRemove; r++) {
            int nrPCAs = r + 1;

            if (!compSkipList.contains("Comp" + nrPCAs)) {
                datasetEV_Mat.viewRow(r).assign(0);
                dataChanged = true;
            }

            if (nrPCAs % nrIntermediatePCAsOverSamplesToRemoveToOutput == 0 && dataChanged == true) {
                dataset.setMatrix(Alg.mult(datasetPCAOverSamplesPCAs_Mat, datasetEV_Mat));
                dataset.save(expressionFile + "." + nrPCAs + "PCAsOverSamplesRemoved.txt.gz");
                System.out.println("Removed\t" + nrPCAs + "\tPCs. File:\t" + expressionFile + "." + nrPCAs + "PCAsOverSamplesRemoved.txt.gz");
                dataChanged = false;
            }

        }
    }

    public static void generatePcCorrectedData(String fileNamePrefix, String PcaFile, String eigenVectorFile, int startSize, int nrPCAsOverSamplesToRemove, int nrIntermediatePCAsOverSamplesToRemoveToOutput) throws IOException {

        DoubleMatrixDataset<String, String> datasetPCAOverSamplesPCAs = DoubleMatrixDataset.loadDoubleData(PcaFile);

        DoubleMatrixDataset<String, String> datasetEV = DoubleMatrixDataset.loadDoubleData(eigenVectorFile);

        DoubleMatrixDataset dataset = new DoubleMatrixDataset(null, datasetPCAOverSamplesPCAs.getHashRows(), datasetEV.getHashCols());

        String expressionFile = fileNamePrefix;
        System.out.println("\nInitializing residual gene expression matrix");

        DenseDoubleAlgebra Alg = new DenseDoubleAlgebra();

        //Clean and optimize for calculation.
        DoubleMatrix2D datasetPCAOverSamplesPCAs_Mat = datasetPCAOverSamplesPCAs.getMatrix();
        DoubleMatrix2D datasetEV_Mat = datasetEV.getMatrix();

        datasetPCAOverSamplesPCAs = null;
        datasetEV = null;

        for (int r = (startSize - 1); r < nrPCAsOverSamplesToRemove; r++) {

//            datasetEV_Mat.viewRow(r).assign(0);
            int nrPCAs = r + 1;
            if (nrPCAs % nrIntermediatePCAsOverSamplesToRemoveToOutput == 0) {
                dataset.setMatrix(Alg.mult(datasetPCAOverSamplesPCAs_Mat.viewPart(0, nrPCAs, datasetPCAOverSamplesPCAs_Mat.rows(), datasetPCAOverSamplesPCAs_Mat.columns() - nrPCAs), datasetEV_Mat.viewPart(nrPCAs, 0, datasetEV_Mat.rows() - nrPCAs, datasetEV_Mat.columns())));
                dataset.save(expressionFile + "." + nrPCAs + "PCAsOverSamplesRemoved.txt.gz");
                System.out.println("Removed\t" + nrPCAs + "\tPCs. File:\t" + expressionFile + "." + nrPCAs + "PCAsOverSamplesRemoved.txt.gz");
            }

        }
    }

    public static boolean isInteger(String string) {
        try {
            Integer.valueOf(string);
            return true;
        } catch (NumberFormatException e) {
            return false;
        }
    }

}
