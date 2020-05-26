/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.legacy;

//import java.io.*;
//import java.util.*;
//import java.util.concurrent.CompletionService;
//import java.util.concurrent.ExecutorCompletionService;
//import java.util.concurrent.ExecutorService;
//import java.util.concurrent.Executors;
//import umcg.genetica.containers.Triple;
//import umcg.genetica.io.concurrent.DoubleParseTask;
//import umcg.genetica.io.text.TextFile;
//import umcg.genetica.util.RunTimer;

/**
 *
 * @author lude
 */
public class ExpressionDataset {

    // DEAR COLLABS:
    // THIS DOES NOT WORK BECAUSE EXPRESSIONDATASET IS A PIECE OF OLD STICKY
    // SHIT.. PLEASE LET ME FIX THIS FOR IT WILL BE REPLACED BY THE NEW
    // AWESOME DOUBLEMATRIXDATASET<U,V>
    // WITH KIND REGARDS,
    // HARM-JAN
    
    
    
//    private double[][] rawData = null;
//    private int nrSamples = 0;
//    private int nrProbes = 0;
//    private String[] probeNames = null;
//    private String[] sampleNames = null;
//    private HashMap<String, Integer> hashSamples = new HashMap<String, Integer>();
//    private HashMap<String, Integer> hashProbes = new HashMap<String, Integer>();
//    private HashSet<String> hashProbesToInclude = null;
//    private HashSet<String> hashSamplesToInclude = null;
//    private String fileName = null;
//
//    public ExpressionDataset(){
//        
//    }
//    
//    public ExpressionDataset(String fileName) throws IOException {
//        System.out.println("Loading: " + fileName);
//        if (fileName.endsWith(".binary")) {
//            loadExpressionDataInBinaryFormat(fileName);
//        } else {
//            loadExpressionData(fileName, "\t");
//        }
//    }
//
//    public ExpressionDataset(String fileName, String delimiter) throws IOException {
//        if (fileName.endsWith(".binary")) {
//            loadExpressionDataInBinaryFormat(fileName);
//        } else {
//            loadExpressionData(fileName, delimiter);
//        }
//    }
//
//    public ExpressionDataset(String fileName, String delimiter, HashSet<String> hashProbesToInclude) throws IOException {
//        this.hashProbesToInclude = hashProbesToInclude;
//        if (fileName.endsWith(".binary")) {
//            loadExpressionDataInBinaryFormat(fileName);
//        } else {
//            loadExpressionData(fileName, delimiter);
//        }
//    }
//
//    public ExpressionDataset(String fileName, String delimiter, HashSet<String> hashProbesToInclude, HashSet<String> hashSamplesToInclude) throws IOException {
//        this.hashProbesToInclude = hashProbesToInclude;
//        this.hashSamplesToInclude = hashSamplesToInclude;
//        if (fileName.endsWith(".binary")) {
//            loadExpressionDataInBinaryFormat(fileName);
//        } else {
//            loadExpressionData(fileName, delimiter);
//        }
//    }
//
//    public ExpressionDataset(int nrProbes, int nrSamples) {
//        this.nrProbes = nrProbes;
//        this.nrSamples = nrSamples;
//        sampleNames = new String[nrSamples];
//        for (int s = 0; s < nrSamples; s++) {
//            sampleNames[s] = "Sample_" + String.valueOf(s);
//        }
//        probeNames = new String[nrProbes];
//        for (int p = 0; p < nrProbes; p++) {
//            probeNames[p] = "Probe_" + String.valueOf(p);
//        }
//        rawData = new double[nrProbes][nrSamples];
//    }
//
//    public void loadExpressionDataInBinaryFormat(String fileName) throws IOException {
//        //First load the raw binary data:
//        this.fileName = fileName;
//        File fileBinary = new File(fileName + ".dat");
//        BufferedInputStream in = null;
//        int nrProbesThisBinaryFile = -1;
//        int nrSamplesThisBinaryFile = -1;
//        try {
//            in = new BufferedInputStream(new FileInputStream(fileBinary));
//            byte[] bytes = new byte[4];
//            in.read(bytes, 0, 4);
//            nrProbesThisBinaryFile = byteArrayToInt(bytes);
//            in.read(bytes, 0, 4);
//            nrSamplesThisBinaryFile = byteArrayToInt(bytes);
//        } catch (FileNotFoundException e) {
//            System.err.println("File " + fileBinary.getName() + " not found: " + e.getMessage());
//            System.exit(-1);
//        }
//
//        if (hashProbesToInclude == null && hashSamplesToInclude == null) {
//
//            //We want to load all the data:
//            nrProbes = nrProbesThisBinaryFile;
//            nrSamples = nrSamplesThisBinaryFile;
//            rawData = new double[nrProbes][nrSamples];
//
//            //Now load the row identifiers from file:
//            probeNames = new String[nrProbes];
//            File fileProbes = new File(fileName + ".rows.txt");
//            try {
//                java.io.BufferedReader inProbes = new java.io.BufferedReader(new java.io.FileReader(fileProbes));
//                for (int p = 0; p < nrProbes; p++) {
//                    probeNames[p] = inProbes.readLine();
//                }
//                inProbes.close();
//            } catch (FileNotFoundException e) {
//                System.err.println("File " + fileProbes.getName() + " not found: " + e.getMessage());
//                System.exit(-1);
//            }
//
//            //Now load the column identifiers from file:
//            sampleNames = new String[nrSamples];
//            File fileSamples = new File(fileName + ".columns.txt");
//            try {
//                java.io.BufferedReader inColumns = new java.io.BufferedReader(new java.io.FileReader(fileSamples));
//                for (int s = 0; s < nrSamples; s++) {
//                    sampleNames[s] = inColumns.readLine();
//                }
//                inColumns.close();
//            } catch (FileNotFoundException e) {
//                System.err.println("File " + fileSamples.getName() + " not found: " + e.getMessage());
//                System.exit(-1);
//            }
//
//            byte[] buffer = new byte[nrSamples * 8];
//            long bits = 0;
//
//            for (int row = 0; row < nrProbes; row++) {
//                in.read(buffer, 0, nrSamples * 8);
//                int bufferLoc = 0;
//                for (int col = 0; col < nrSamples; col++) {
//                    bits = (long) (0xff & buffer[bufferLoc + 7])
//                            | (long) (0xff & buffer[bufferLoc + 6]) << 8
//                            | (long) (0xff & buffer[bufferLoc + 5]) << 16
//                            | (long) (0xff & buffer[bufferLoc + 4]) << 24
//                            | (long) (0xff & buffer[bufferLoc + 3]) << 32
//                            | (long) (0xff & buffer[bufferLoc + 2]) << 40
//                            | (long) (0xff & buffer[bufferLoc + 1]) << 48
//                            | (long) (buffer[bufferLoc]) << 56;
//
//                    rawData[row][col] = Double.longBitsToDouble(bits);
//                    bufferLoc += 8;
//                }
//            }
//            in.close();
//
//
//        } else {
//
//            //We want to confine the set of probes and samples to a subset. Deal with this in a different way.
//
//            //Now load the row identifiers from file:
//            File fileProbes = new File(fileName + ".rows.txt");
//            HashMap<String, Integer> hashProbesPresentAndRequested = new HashMap<String, Integer>();
//            int[] probeSubsetIndex = new int[nrProbesThisBinaryFile];
//            for (int p = 0; p < probeSubsetIndex.length; p++) {
//                probeSubsetIndex[p] = -1;
//            }
//            try {
//                java.io.BufferedReader inProbes = new java.io.BufferedReader(new java.io.FileReader(fileProbes));
//                String str = "";
//                int probeIndex = 0;
//                while ((str = inProbes.readLine()) != null) {
//                    if (hashProbesToInclude == null || hashProbesToInclude.contains(str)) {
//                        probeSubsetIndex[probeIndex] = hashProbesPresentAndRequested.size();
//                        hashProbesPresentAndRequested.put(str, probeIndex);
//                    }
//                    probeIndex++;
//                }
//                inProbes.close();
//                if (nrProbesThisBinaryFile != probeIndex) {
//                    System.out.println("Number of probes in binary data file does not correspond to the number of probes in the .rows.txt file!!!");
//                    System.exit(-1);
//                }
//            } catch (FileNotFoundException e) {
//                System.err.println("File " + fileProbes.getName() + " not found: " + e.getMessage());
//                System.exit(-1);
//            }
//
//            nrProbes = hashProbesPresentAndRequested.size();
//            probeNames = new String[nrProbes];
//            try {
//                java.io.BufferedReader inProbes = new java.io.BufferedReader(new java.io.FileReader(fileProbes));
//                String str = "";
//                int probeCounter = 0;
//                while ((str = inProbes.readLine()) != null) {
//                    if (hashProbesToInclude == null || hashProbesToInclude.contains(str)) {
//                        probeNames[probeCounter] = str;
//                        probeCounter++;
//                    }
//                }
//                inProbes.close();
//            } catch (FileNotFoundException e) {
//                System.err.println("File " + fileProbes.getName() + " not found: " + e.getMessage());
//                System.exit(-1);
//            }
//
//            //Now load the column identifiers from file:
//            sampleNames = new String[nrSamples];
//            File fileSamples = new File(fileName + ".columns.txt");
//            HashMap<String, Integer> hashSamplesPresentAndRequested = new HashMap<String, Integer>();
//            int[] sampleSubsetIndex = new int[nrSamplesThisBinaryFile];
//            for (int s = 0; s < sampleSubsetIndex.length; s++) {
//                sampleSubsetIndex[s] = -1;
//            }
//            try {
//                java.io.BufferedReader inColumns = new java.io.BufferedReader(new java.io.FileReader(fileSamples));
//                String str = "";
//                int sampleIndex = 0;
//                while ((str = inColumns.readLine()) != null) {
//                    if (hashSamplesToInclude == null || hashSamplesToInclude.contains(str)) {
//                        sampleSubsetIndex[sampleIndex] = hashSamplesPresentAndRequested.size();
//                        hashSamplesPresentAndRequested.put(str, sampleIndex);
//                    }
//                    sampleIndex++;
//                }
//                inColumns.close();
//                if (nrSamplesThisBinaryFile != sampleIndex) {
//                    System.out.println("Number of samples in binary data file does not correspond to the number of samples in the .columns.txt file!!!");
//                    System.exit(-1);
//                }
//
//            } catch (FileNotFoundException e) {
//                System.err.println("File " + fileSamples.getName() + " not found: " + e.getMessage());
//                System.exit(-1);
//            }
//            nrSamples = hashSamplesPresentAndRequested.size();
//            sampleNames = new String[nrSamples];
//            try {
//                java.io.BufferedReader inColumns = new java.io.BufferedReader(new java.io.FileReader(fileSamples));
//                String str = "";
//                int sampleCounter = 0;
//                while ((str = inColumns.readLine()) != null) {
//                    if (hashSamplesToInclude == null || hashSamplesToInclude.contains(str)) {
//                        sampleNames[sampleCounter] = str;
//                        sampleCounter++;
//                    }
//                }
//                inColumns.close();
//            } catch (FileNotFoundException e) {
//                System.err.println("File " + fileSamples.getName() + " not found: " + e.getMessage());
//                System.exit(-1);
//            }
//
//            //Now load the binary data:
//            rawData = new double[nrProbes][nrSamples];
//
//            byte[] buffer = new byte[nrSamplesThisBinaryFile * 8];
//            long bits = 0;
//
//            for (int row = 0; row < nrProbesThisBinaryFile; row++) {
//                in.read(buffer, 0, nrSamplesThisBinaryFile * 8);
//                int bufferLoc = 0;
//                for (int col = 0; col < nrSamplesThisBinaryFile; col++) {
//                    bits = (long) (0xff & buffer[bufferLoc + 7])
//                            | (long) (0xff & buffer[bufferLoc + 6]) << 8
//                            | (long) (0xff & buffer[bufferLoc + 5]) << 16
//                            | (long) (0xff & buffer[bufferLoc + 4]) << 24
//                            | (long) (0xff & buffer[bufferLoc + 3]) << 32
//                            | (long) (0xff & buffer[bufferLoc + 2]) << 40
//                            | (long) (0xff & buffer[bufferLoc + 1]) << 48
//                            | (long) (buffer[bufferLoc]) << 56;
//
//                    int rowIndex = probeSubsetIndex[row];
//                    int colIndex = sampleSubsetIndex[col];
//                    if (rowIndex != -1 && colIndex != -1) {
//                        rawData[rowIndex][colIndex] = Double.longBitsToDouble(bits);
//                    }
//                    bufferLoc += 8;
//                }
//            }
//            in.close();
//        }
//
//        recalculateHashMaps();
//        System.out.println("Binary file:\t" + fileName + "\thas been loaded, nrProbes:\t" + nrProbes + "\tnrSamples:\t" + nrSamples);
//    }
//
//    public void loadExpressionData(String fileName, String delimiter) throws IOException {
//        TextFile in = new TextFile(fileName, TextFile.R);
//
//        this.fileName = fileName;
//        boolean dataIsInTriTyperFormat = false;
//
//        int sampleOffset = 1;
//        int[] sampleIndex = null;
//
//        String str = in.readLine();
//        String[] data = str.split(delimiter);
//        if (data.length > 2 && data[1].length() > 0 && data[1].equals("MultipleHits")) {
//            dataIsInTriTyperFormat = true;
//            sampleOffset = 9;
//
//        }
//
//        if (hashSamplesToInclude == null) {
//            nrSamples = data.length - sampleOffset;
//            sampleNames = new String[nrSamples];
//            sampleIndex = new int[nrSamples];
//            for (int s = 0; s < nrSamples; s++) {
//                String samplename = data[s + sampleOffset];
//                sampleNames[s] = samplename;
//                if (hashSamples.containsKey(samplename)) {
//                    System.err.println("WARNING: duplicate column name detected in file: " + samplename);
//                }
//                hashSamples.put(samplename, s);
//                sampleIndex[s] = s;
//            }
//        } else {
//            ArrayList<Integer> vecSampleIndex = new ArrayList<Integer>();
//            ArrayList<String> vecSampleName = new ArrayList<String>();
//            for (int s = 0; s < data.length - sampleOffset; s++) {
//                String sample = data[s + sampleOffset];
//                if (hashSamplesToInclude.contains(sample)) {
//                    vecSampleIndex.add(s);
//                    vecSampleName.add(sample);
//                }
//            }
//            nrSamples = vecSampleIndex.size();
//            sampleNames = new String[nrSamples];
//            sampleIndex = new int[nrSamples];
//            for (int s = 0; s < nrSamples; s++) {
//                String sample = vecSampleName.get(s);
//                sampleNames[s] = sample;
//                hashSamples.put(sample, s);
//                sampleIndex[s] = vecSampleIndex.get(s);
//            }
//        }
//        nrProbes = 0;
//
//        while ((str = in.readLine()) != null) {
//            if (hashProbesToInclude == null) {
//                nrProbes++;
////                if (nrProbes % 1000 == 0) {
////                    System.out.println(nrProbes);
////                }
//            } else {
//                data = str.split(delimiter);
//                if (hashProbesToInclude.contains(data[0])) {
//                    nrProbes++;
//                    //if (nrProbes%1000==0) System.out.println(nrProbes);
//                }
//            }
//        }
//        in.close();
//
//        rawData = new double[nrProbes][nrSamples];
//        probeNames = new String[nrProbes];
//
//        in.open();
//
//        int nrprocs = Runtime.getRuntime().availableProcessors();
//        System.out.println("Using " + nrprocs + " threads for parsing the file");
//        ExecutorService threadPool = Executors.newFixedThreadPool(nrprocs);
//        CompletionService<Triple<Integer, String, double[]>> pool = new ExecutorCompletionService<Triple<Integer, String, double[]>>(threadPool);
//
//        str = in.readLine(); // skip header
//        nrProbes = 0;
//        int lnctr = 0;
//        int tasksSubmitted = 0;
//        int returnedResults = 0;
//
//        RunTimer timer = new RunTimer();
//        timer.start();
//        while ((str = in.readLine()) != null) {
//
//            DoubleParseTask task = new DoubleParseTask(str, sampleOffset, lnctr, sampleIndex, hashProbesToInclude);
//            pool.submit(task);
//
//            tasksSubmitted++;
//
//            if (lnctr % (nrprocs * 2) == 0) {
//                while (returnedResults < tasksSubmitted) {
//                    try {
//                        Triple<Integer, String, double[]> result = pool.take().get();
//                        if (result != null) {
//                            int rownr = result.getLeft(); //  < 0 when row is not to be included because of hashProbesToInclude.
//                            if (rownr >= 0) {
//                                probeNames[rownr] = result.getMiddle();
//                                double[] doubles = result.getRight();
//                                rawData[rownr] = doubles;
//                            }
//                            result = null;
//                            returnedResults++;
//                        }
//                    } catch (Exception e) {
//                        e.printStackTrace();
//                    }
//                }
//            }
//
//            lnctr++;
//
//            nrProbes++;
////            if (nrProbes % 100 == 0) {
////                System.out.println(nrProbes);
////            }
////            }
//        }
//        in.close();
//
//        while (returnedResults < tasksSubmitted) {
//            try {
//                Triple<Integer, String, double[]> result = pool.take().get();
//                if (result != null) {
//                    int rownr = result.getLeft(); //  < 0 when row is not to be included because of hashProbesToInclude.
//                    if (rownr >= 0) {
//                        probeNames[rownr] = result.getMiddle();
//                        double[] doubles = result.getRight();
//                        rawData[rownr] = doubles;
//
////                            System.out.println(rownr + "\t " + returnedResults);
//                    }
//                    result = null;
//                    returnedResults++;
//                }
//            } catch (Exception e) {
//                e.printStackTrace();
//            }
//        }
//
//        timer.stop();
//
//
//        threadPool.shutdown();
//
//        System.out.println(fileName + "\thas been loaded, nrProbes:\t" + nrProbes + "\tnrSamples:\t" + nrSamples);
//        System.out.println("Time needed: " + timer.getTimeDesc());
//    }
//
//    public void recalculateHashMaps() {
//
//        hashProbes.clear();
//        for (int probeItr = 0; probeItr < nrProbes; probeItr++) {
//            hashProbes.put(probeNames[probeItr], probeItr);
//        }
//
//        hashSamples.clear();
//        for (int sampleItr = 0; sampleItr < nrSamples; sampleItr++) {
//            hashSamples.put(sampleNames[sampleItr], sampleItr);
//        }
//
//    }
//
//    public double[][] getRawDataTransposed() {
//        double[][] rawDataTransposed = new double[nrSamples][nrProbes];
//        for (int s = 0; s < nrSamples; s++) {
//            for (int p = 0; p < nrProbes; p++) {
//                rawDataTransposed[s][p] = rawData[p][s];
//            }
//        }
//        return rawDataTransposed;
//    }
//
//    public void standardNormalizeData() {
//
//        /*
//         * System.out.println("\nNormalizing data:"); //Calculate the average
//         * expression, when per sample all raw expression levels have been
//         * ordered: double[] rankedMean = new double[nrProbes]; for (int
//         * probeID=0; probeID<nrProbes; probeID++) { double quantile = ((double)
//         * probeID + 1.0d) / ((double) nrProbes + 1d); rankedMean[probeID] =
//         * cern.jet.stat.Probability.normalInverse(quantile); }
//         *
//         * //Iterate through each sample: hgea.math.RankDoubleArray
//         * rankDoubleArray = new hgea.math.RankDoubleArray(); for (int s=0;
//         * s<nrSamples; s++) { double[] probes = new double[nrProbes]; for (int
//         * p=0; p<nrProbes; p++) probes[p]=rawData[p][s]; double[] probesRanked
//         * = rankDoubleArray.rank(probes); double[] probesQuantileNormalized =
//         * new double[nrProbes]; for (int p=0; p<nrProbes; p++) {
//         * probesQuantileNormalized[p] = rankedMean[(int) probesRanked[p]]; }
//         * for (int p=0; p<nrProbes; p++) rawData[p][s] = (float)
//         * probesQuantileNormalized[p]; }
//         */
//
//
//        System.out.println("Setting probe mean to zero and stdev to one for every probe:");
//        for (int probeID = 0; probeID < nrProbes; probeID++) {
//            double vals[] = new double[nrSamples];
//            for (int s = 0; s < nrSamples; s++) {
//                vals[s] = rawData[probeID][s];
//            }
//            double mean = JSci.maths.ArrayMath.mean(vals);
//            for (int s = 0; s < nrSamples; s++) {
//                vals[s] -= (double) mean;
//            }
//            double standardDeviation = JSci.maths.ArrayMath.standardDeviation(vals);
//            for (int s = 0; s < nrSamples; s++) {
//                rawData[probeID][s] = (float) (vals[s] / standardDeviation);
//            }
//        }
//
//    }
//
//    public void save(String fileName, boolean gz) {
//        if (fileName.endsWith(".binary")) {
//
//            //Create binary file:
//            BufferedOutputStream out = null;
//            File fileBinary = new File(fileName + ".dat");
//            try {
//                out = new BufferedOutputStream(new FileOutputStream(fileBinary));
//                out.write(intToByteArray(nrProbes));
//                out.write(intToByteArray(nrSamples));
//            } catch (IOException e) {
//                System.err.println("Can't write to " + fileBinary.getName() + ": " + e.getMessage());
//                System.exit(1);
//            }
//            byte[] buffer = new byte[rawData[0].length * 8];
//            for (int row = 0; row < rawData.length; row++) { // rows
//                int bufferLoc = 0;
//                for (int col = 0; col < rawData[0].length; col++) { // columns
//                    long bits = Double.doubleToLongBits(rawData[row][col]);
//                    buffer[bufferLoc] = (byte) (bits >> 56);
//                    buffer[bufferLoc + 1] = (byte) (bits >> 48 & 0xff);
//                    buffer[bufferLoc + 2] = (byte) (bits >> 40 & 0xff);
//                    buffer[bufferLoc + 3] = (byte) (bits >> 32 & 0xff);
//                    buffer[bufferLoc + 4] = (byte) (bits >> 24 & 0xff);
//                    buffer[bufferLoc + 5] = (byte) (bits >> 16 & 0xff);
//                    buffer[bufferLoc + 6] = (byte) (bits >> 8 & 0xff);
//                    buffer[bufferLoc + 7] = (byte) (bits & 0xff);
//                    bufferLoc += 8;
//                }
//                try {
//                    out.write(buffer);
//                } catch (IOException e) {
//                    System.err.println("Can't write to " + fileBinary.getName() + ": " + e.getMessage());
//                    System.exit(1);
//                }
//            }
//            try {
//                out.close();
//            } catch (IOException e) {
//                e.printStackTrace();
//            }
//            File fileProbes = new File(fileName + ".rows.txt");
//            try {
//                java.io.BufferedWriter outProbes = new java.io.BufferedWriter(new java.io.FileWriter(fileProbes));
//                for (int p = 0; p < nrProbes; p++) {
//                    outProbes.write(probeNames[p] + "\n");
//                }
//                outProbes.flush();
//                outProbes.close();
//            } catch (Exception e) {
//                System.out.println("Error:\t" + e.getMessage());
//                e.printStackTrace();
//            }
//            File fileSamples = new File(fileName + ".columns.txt");
//            try {
//                java.io.BufferedWriter outSamples = new java.io.BufferedWriter(new java.io.FileWriter(fileSamples));
//                for (int s = 0; s < nrSamples; s++) {
//                    outSamples.write(sampleNames[s] + "\n");
//                }
//                outSamples.flush();
//                outSamples.close();
//            } catch (Exception e) {
//                System.out.println("Error:\t" + e.getMessage());
//                e.printStackTrace();
//            }
//
//
//        } else {
//            try {
//                TextFile out = new TextFile(fileName, TextFile.W);
//
//                // java.io.BufferedWriter out = new java.io.BufferedWriter(new java.io.FileWriter(fileOut));
//
//                StringBuffer sb = new StringBuffer("-");
//                for (int s = 0; s < nrSamples; s++) {
//                    sb.append("\t").append(sampleNames[s]);
//                }
//                out.write(sb.toString() + "\n");
//                for (int p = 0; p < nrProbes; p++) {
//                    sb = new StringBuffer(probeNames[p]);
//                    for (int s = 0; s < nrSamples; s++) {
//                        sb.append("\t").append(rawData[p][s]);
//                    }
//                    out.write(sb.toString() + "\n");
//                }
//                out.close();
//            } catch (Exception e) {
//                System.out.println("Error:\t" + e.getMessage());
//                e.printStackTrace();
//            }
//        }
//    }
//
//    public void transposeDataset() {
//        rawData = getRawDataTransposed();
//
//        int nrProbesTemp = nrProbes;
//        String[] probeNamesTemp = probeNames;
//        HashMap<String, Integer> hashProbesTemp = hashProbes;
//
//        nrProbes = nrSamples;
//        probeNames = sampleNames;
//        hashProbes = hashSamples;
//
//        nrSamples = nrProbesTemp;
//        sampleNames = probeNamesTemp;
//        hashSamples = hashProbesTemp;
//
//    }
//
//    private byte[] intToByteArray(int value) {
//        return new byte[]{(byte) (value >>> 24),
//                    (byte) (value >>> 16),
//                    (byte) (value >>> 8),
//                    (byte) value};
//    }
//
//    private int byteArrayToInt(byte[] b) {
//        return (b[0] << 24)
//                + ((b[1] & 0xff) << 16)
//                + ((b[2] & 0xff) << 8)
//                + (b[3] & 0xff);
//    }
//
//    /**
//     * @return the rawData
//     */
//    public double[][] getRawData() {
//        return rawData;
//    }
//
//    /**
//     * @param rawData the rawData to set
//     */
//    public void setRawData(double[][] rawData) {
//        this.rawData = rawData;
//        if (rawData != null) {
//            this.nrProbes = rawData.length;
//            if (this.nrProbes > 0) {
//                this.nrSamples = rawData[0].length;
//            }
//            System.out.println("New matrix is "+nrProbes+"x"+nrSamples);
//        }
//    }
//
//    /**
//     * @return the nrSamples
//     */
//    public int getNrSamples() {
//        return nrSamples;
//    }
//
//    /**
//     * @param nrSamples the nrSamples to set
//     */
//    public void setNrSamples(int nrSamples) {
//        this.nrSamples = nrSamples;
//    }
//
//    /**
//     * @return the nrProbes
//     */
//    public int getNrProbes() {
//        return nrProbes;
//    }
//
//    /**
//     * @param nrProbes the nrProbes to set
//     */
//    public void setNrProbes(int nrProbes) {
//        this.nrProbes = nrProbes;
//    }
//
//    /**
//     * @return the probeNames
//     */
//    public String[] getProbeNames() {
//        return probeNames;
//    }
//
//    /**
//     * @param probeNames the probeNames to set
//     */
//    public void setProbeNames(String[] probeNames) {
//        this.probeNames = probeNames;
//    }
//
//    /**
//     * @return the sampleNames
//     */
//    public String[] getSampleNames() {
//        return sampleNames;
//    }
//
//    /**
//     * @param sampleNames the sampleNames to set
//     */
//    public void setSampleNames(String[] sampleNames) {
//        this.sampleNames = sampleNames;
//    }
//
//    /**
//     * @return the hashSamples
//     */
//    public HashMap getHashSamples() {
//        return hashSamples;
//    }
//
//    /**
//     * @param hashSamples the hashSamples to set
//     */
//    public void setHashSamples(HashMap<String, Integer> hashSamples) {
//        this.hashSamples = hashSamples;
//    }
//
//    /**
//     * @return the hashProbes
//     */
//    public HashMap getHashProbes() {
//        return hashProbes;
//    }
//
//    /**
//     * @param hashProbes the hashProbes to set
//     */
//    public void setHashProbes(HashMap<String, Integer> hashProbes) {
//        this.hashProbes = hashProbes;
//    }
//
//    /**
//     * @return the hashProbesToInclude
//     */
//    public HashSet<String> getHashProbesToInclude() {
//        return hashProbesToInclude;
//    }
//
//    /**
//     * @param hashProbesToInclude the hashProbesToInclude to set
//     */
//    public void setHashProbesToInclude(HashSet<String> hashProbesToInclude) {
//        this.hashProbesToInclude = hashProbesToInclude;
//    }
//
//    /**
//     * @return the hashSamplesToInclude
//     */
//    public HashSet<String> getHashSamplesToInclude() {
//        return hashSamplesToInclude;
//    }
//
//    /**
//     * @param hashSamplesToInclude the hashSamplesToInclude to set
//     */
//    public void setHashSamplesToInclude(HashSet<String> hashSamplesToInclude) {
//        this.hashSamplesToInclude = hashSamplesToInclude;
//    }
//
//    /**
//     * @return the fileName
//     */
//    public String getFileName() {
//        return fileName;
//    }
//
//    /**
//     * @param fileName the fileName to set
//     */
//    public void setFileName(String fileName) {
//        this.fileName = fileName;
//    }
}
