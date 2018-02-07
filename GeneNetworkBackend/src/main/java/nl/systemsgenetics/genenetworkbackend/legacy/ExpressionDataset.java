package nl.systemsgenetics.genenetworkbackend.legacy;

import java.io.*;
import java.net.*;
import java.util.*;
import java.awt.image.BufferedImage;
import java.awt.image.*;
import java.awt.*;
import java.awt.geom.*;
import java.lang.Math;
import javax.imageio.*;

/**
 *
 * @author lude
 */
public class ExpressionDataset {

    public double[][] rawData = null;
    public int nrSamples = 0;
    public int nrProbes = 0;
    public String[] probeNames = null;
    public String[] sampleNames = null;
    public HashMap hashSamples = new HashMap();
    public HashMap hashProbes = new HashMap();
    private HashMap hashProbesToInclude = null;
    private HashMap hashSamplesToInclude = null;
    public String fileName = null;

    public ExpressionDataset(String fileName) {
        if (fileName.endsWith(".binary")) {
            loadExpressionDataInBinaryFormat(fileName);
        } else {
            loadExpressionData(fileName, "\t");
        }
    }

    public ExpressionDataset(String fileName, String delimiter) {
        if (fileName.endsWith(".binary")) {
            loadExpressionDataInBinaryFormat(fileName);
        } else {
            loadExpressionData(fileName, delimiter);
        }
    }

    public ExpressionDataset(String fileName, String delimiter, HashMap hashProbesToInclude) {
        this.hashProbesToInclude = hashProbesToInclude;
        if (fileName.endsWith(".binary")) {
            loadExpressionDataInBinaryFormat(fileName);
        } else {
            loadExpressionData(fileName, delimiter);
        }
    }

    public ExpressionDataset(String fileName, String delimiter, HashMap hashProbesToInclude, HashMap hashSamplesToInclude) {
        this.hashProbesToInclude = hashProbesToInclude;
        this.hashSamplesToInclude = hashSamplesToInclude;
        if (fileName.endsWith(".binary")) {
            loadExpressionDataInBinaryFormat(fileName);
        } else {
            loadExpressionData(fileName, delimiter);
        }
    }

    public ExpressionDataset(int nrProbes, int nrSamples) {
        this.nrProbes = nrProbes;
        this.nrSamples = nrSamples;
        sampleNames = new String[nrSamples];
        for (int s=0; s<nrSamples; s++) sampleNames[s] = "Sample_" + String.valueOf(s);
        probeNames = new String[nrProbes];
        for (int p=0; p<nrProbes; p++) probeNames[p] = "Probe_" + String.valueOf(p);
        rawData = new double[nrProbes][nrSamples];
    }

    public void loadExpressionDataInBinaryFormat(String fileName) {
        //First load the raw binary data:
        this.fileName = fileName;
        File fileBinary = new File(fileName + ".dat");
        BufferedInputStream in = null;
        int nrProbesThisBinaryFile = -1;
        int nrSamplesThisBinaryFile = -1;
        try {
            in = new BufferedInputStream(new FileInputStream(fileBinary));
            byte[] bytes = new byte[4];
            in.read(bytes, 0, 4);
            nrProbesThisBinaryFile = byteArrayToInt(bytes);
            in.read(bytes, 0, 4);
            nrSamplesThisBinaryFile = byteArrayToInt(bytes);
        } catch (FileNotFoundException e) {
            System.err.println("File " + fileBinary.getName() + " not found: " + e.getMessage());
            System.exit(-1);
        } catch (IOException e) {
            System.err.println("Can't read " + fileBinary.getName() + ": " + e.getMessage());
            System.exit(-1);
        }

        if (hashProbesToInclude == null && hashSamplesToInclude == null) {

            //We want to load all the data:
            nrProbes = nrProbesThisBinaryFile;
            nrSamples = nrSamplesThisBinaryFile;
            rawData = new double[nrProbes][nrSamples];

            //Now load the row identifiers from file:
            probeNames = new String[nrProbes];
            File fileProbes = new File(fileName + ".rows.txt");
            try {
                java.io.BufferedReader inProbes = new java.io.BufferedReader(new java.io.FileReader(fileProbes));
                for (int p=0; p<nrProbes; p++) {
                    probeNames[p] = inProbes.readLine();
                }
                inProbes.close();
            } catch (FileNotFoundException e) {
                System.err.println("File " + fileProbes.getName() + " not found: " + e.getMessage());
                System.exit(-1);
            } catch (IOException e) {
                System.err.println("Can't read " + fileProbes.getName() + ": " + e.getMessage());
                System.exit(-1);
            }

            //Now load the column identifiers from file:
            sampleNames = new String[nrSamples];
            File fileSamples = new File(fileName + ".columns.txt");
            try {
                java.io.BufferedReader inColumns = new java.io.BufferedReader(new java.io.FileReader(fileSamples));
                for (int s=0; s<nrSamples; s++) {
                    sampleNames[s] = inColumns.readLine();
                }
                inColumns.close();
            } catch (FileNotFoundException e) {
                System.err.println("File " + fileSamples.getName() + " not found: " + e.getMessage());
                System.exit(-1);
            } catch (IOException e) {
                System.err.println("Can't read " + fileSamples.getName() + ": " + e.getMessage());
                System.exit(-1);
            }

            byte[] buffer = new byte[nrSamples * 8];
            long bits = 0;
            try {
                for (int row = 0; row < nrProbes; row++) {
                    in.read(buffer, 0, nrSamples * 8);
                    int bufferLoc = 0;
                    for (int col = 0; col < nrSamples; col++) {
                        bits = (long) (0xff & buffer[bufferLoc + 7])
                                | (long) (0xff & buffer[bufferLoc + 6]) << 8
                                | (long) (0xff & buffer[bufferLoc + 5]) << 16
                                | (long) (0xff & buffer[bufferLoc + 4]) << 24
                                | (long) (0xff & buffer[bufferLoc + 3]) << 32
                                | (long) (0xff & buffer[bufferLoc + 2]) << 40
                                | (long) (0xff & buffer[bufferLoc + 1]) << 48
                                | (long) (buffer[bufferLoc]) << 56;

                        rawData[row][col] = Double.longBitsToDouble(bits);
                        bufferLoc += 8;
                    }
                }
                in.close();
            } catch (IOException e) {
                System.err.println("Can't read from binary file: " + fileBinary.getName() + ": " + e.getMessage());
                System.exit(1);
            }

        } else {

            //We want to confine the set of probes and samples to a subset. Deal with this in a different way.

            //Now load the row identifiers from file:
            File fileProbes = new File(fileName + ".rows.txt");
            HashMap hashProbesPresentAndRequested = new HashMap();
            int[] probeSubsetIndex = new int[nrProbesThisBinaryFile];
            for (int p=0; p<probeSubsetIndex.length; p++) probeSubsetIndex[p] = -1;
            try {
                java.io.BufferedReader inProbes = new java.io.BufferedReader(new java.io.FileReader(fileProbes));
                String str = "";
                int probeIndex = 0;
                while ((str = inProbes.readLine()) != null) {
                    if (hashProbesToInclude==null || hashProbesToInclude.containsKey(str)) {
                        probeSubsetIndex[probeIndex]  = hashProbesPresentAndRequested.size();
                        hashProbesPresentAndRequested.put(str, probeIndex);
                    }
                    probeIndex++;
                }
                inProbes.close();
                if (nrProbesThisBinaryFile != probeIndex) {
                    System.out.println("Number of probes in binary data file does not correspond to the number of probes in the .rows.txt file!!!");
                    System.exit(-1);
                }
            } catch (FileNotFoundException e) {
                System.err.println("File " + fileProbes.getName() + " not found: " + e.getMessage());
                System.exit(-1);
            } catch (IOException e) {
                System.err.println("Can't read " + fileProbes.getName() + ": " + e.getMessage());
                System.exit(-1);
            }

            nrProbes = hashProbesPresentAndRequested.size();
            probeNames = new String[nrProbes];
            try {
                java.io.BufferedReader inProbes = new java.io.BufferedReader(new java.io.FileReader(fileProbes));
                String str = "";
                int probeCounter = 0;
                while ((str = inProbes.readLine()) != null) {
                    if (hashProbesToInclude==null || hashProbesToInclude.containsKey(str)) {
                        probeNames[probeCounter] = str;
                        probeCounter++;
                    }
                }
                inProbes.close();
            } catch (FileNotFoundException e) {
                System.err.println("File " + fileProbes.getName() + " not found: " + e.getMessage());
                System.exit(-1);
            } catch (IOException e) {
                System.err.println("Can't read " + fileProbes.getName() + ": " + e.getMessage());
                System.exit(-1);
            }

            //Now load the column identifiers from file:
            sampleNames = new String[nrSamples];
            File fileSamples = new File(fileName + ".columns.txt");
            HashMap hashSamplesPresentAndRequested = new HashMap();
            int[] sampleSubsetIndex = new int[nrSamplesThisBinaryFile];
            for (int s=0; s<sampleSubsetIndex.length; s++) sampleSubsetIndex[s] = -1;
            try {
                java.io.BufferedReader inColumns = new java.io.BufferedReader(new java.io.FileReader(fileSamples));
                String str = "";
                int sampleIndex = 0;
                while ((str = inColumns.readLine()) != null) {
                    if (hashSamplesToInclude==null || hashSamplesToInclude.containsKey(str)) {
                        sampleSubsetIndex[sampleIndex] = hashSamplesPresentAndRequested.size();
                        hashSamplesPresentAndRequested.put(str, sampleIndex);
                    }
                    sampleIndex++;
                }
                inColumns.close();
                if (nrSamplesThisBinaryFile != sampleIndex) {
                    System.out.println("Number of samples in binary data file does not correspond to the number of samples in the .columns.txt file!!!");
                    System.exit(-1);
                }

            } catch (FileNotFoundException e) {
                System.err.println("File " + fileSamples.getName() + " not found: " + e.getMessage());
                System.exit(-1);
            } catch (IOException e) {
                System.err.println("Can't read " + fileSamples.getName() + ": " + e.getMessage());
                System.exit(-1);
            }
            nrSamples = hashSamplesPresentAndRequested.size();
            sampleNames = new String[nrSamples];
            try {
                java.io.BufferedReader inColumns = new java.io.BufferedReader(new java.io.FileReader(fileSamples));
                String str = "";
                int sampleCounter = 0;
                while ((str = inColumns.readLine()) != null) {
                    if (hashSamplesToInclude==null || hashSamplesToInclude.containsKey(str)) {
                        sampleNames[sampleCounter] = str;
                        sampleCounter++;
                    }
                }
                inColumns.close();
            } catch (FileNotFoundException e) {
                System.err.println("File " + fileSamples.getName() + " not found: " + e.getMessage());
                System.exit(-1);
            } catch (IOException e) {
                System.err.println("Can't read " + fileSamples.getName() + ": " + e.getMessage());
                System.exit(-1);
            }

            //Now load the binary data:
            rawData = new double[nrProbes][nrSamples];

            byte[] buffer = new byte[nrSamplesThisBinaryFile * 8];
            long bits = 0;
            try {
                for (int row = 0; row < nrProbesThisBinaryFile; row++) {
                    in.read(buffer, 0, nrSamplesThisBinaryFile * 8);
                    int bufferLoc = 0;
                    for (int col = 0; col < nrSamplesThisBinaryFile; col++) {
                        bits = (long) (0xff & buffer[bufferLoc + 7])
                                | (long) (0xff & buffer[bufferLoc + 6]) << 8
                                | (long) (0xff & buffer[bufferLoc + 5]) << 16
                                | (long) (0xff & buffer[bufferLoc + 4]) << 24
                                | (long) (0xff & buffer[bufferLoc + 3]) << 32
                                | (long) (0xff & buffer[bufferLoc + 2]) << 40
                                | (long) (0xff & buffer[bufferLoc + 1]) << 48
                                | (long) (buffer[bufferLoc]) << 56;

                        int rowIndex = probeSubsetIndex[row];
                        int colIndex = sampleSubsetIndex[col];
                        if (rowIndex!=-1 && colIndex!=-1) {
                            rawData[rowIndex][colIndex] = Double.longBitsToDouble(bits);
                        }
                        bufferLoc += 8;
                    }
                }
                in.close();
            } catch (IOException e) {
                System.err.println("Can't read from binary file: " + fileBinary.getName() + ": " + e.getMessage());
                System.exit(1);
            }
            
            

        }

        recalculateHashMaps();
        System.out.println("Binary file:\t" + fileName + "\thas been loaded, nrProbes:\t" + nrProbes + "\tnrSamples:\t" + nrSamples);
    }

    public void loadExpressionData(String fileName, String delimiter) {
        this.fileName = fileName;
        boolean dataIsInTriTyperFormat = false;
        File file = new File(fileName);
        if (!file.canRead()) {
            System.out.println("Error! Cannot open file:\t" + fileName);
            System.exit(0);
        }
        int sampleOffset = 1;
        int[] sampleIndex = null;
        try {
            java.io.BufferedReader in = new java.io.BufferedReader(new java.io.FileReader(file));
            String str = in.readLine();
            String[] data = str.split(delimiter);
            if (data.length>2 && data[1].length() > 0 && data[1].equals("MultipleHits")) {
                dataIsInTriTyperFormat = true;
                sampleOffset = 9;

            }

            if (hashSamplesToInclude==null) {
                nrSamples = data.length - sampleOffset;
                sampleNames = new String[nrSamples];
                sampleIndex = new int[nrSamples];
                for (int s=0; s<nrSamples; s++) {
                    sampleNames[s] = data[s + sampleOffset];
                    hashSamples.put(sampleNames[s], s);
                    sampleIndex[s] = s;
                }
            } else {
                Vector vecSampleIndex = new Vector();
                Vector vecSampleName = new Vector();
                for (int s=0; s<data.length - sampleOffset; s++) {
                    String sample = data[s + sampleOffset];
                    if (hashSamplesToInclude.containsKey(sample)) {
                        vecSampleIndex.add(s);
                        vecSampleName.add(sample);
                    }
                }
                nrSamples = vecSampleIndex.size();
                sampleNames = new String[nrSamples];
                sampleIndex = new int[nrSamples];
                for (int s=0; s<nrSamples; s++) {
                    String sample = (String) vecSampleName.get(s);
                    sampleNames[s] = sample;
                    hashSamples.put(sampleNames[s], s);
                    sampleIndex[s] = ((Integer) vecSampleIndex.get(s)).intValue();
                }
            }
            nrProbes = 0;
            while ((str = in.readLine()) != null) {
                if (hashProbesToInclude==null) {
                    nrProbes++;
                    //if (nrProbes%1000==0) System.out.println(nrProbes);
                } else {
                    data = str.split(delimiter);
                    if (hashProbesToInclude.containsKey(data[0])) {
                        nrProbes++;
                        //if (nrProbes%1000==0) System.out.println(nrProbes);
                    }
                }
            }
            in.close();
        } catch (Exception e) {
            System.out.println("Error:\t" + e.getMessage());
            e.printStackTrace();
        }
        rawData = new double[nrProbes][nrSamples];
        probeNames = new String[nrProbes];
        try {
            java.io.BufferedReader in = new java.io.BufferedReader(new java.io.FileReader(file));
            String str = in.readLine();
            nrProbes = 0;
            while ((str = in.readLine()) != null) {
                String[] data = str.split(delimiter);
                if (hashProbesToInclude==null || hashProbesToInclude.containsKey(data[0])) {
                    probeNames[nrProbes] = new String(data[0].getBytes());
                    hashProbes.put(probeNames[nrProbes], nrProbes);
                    for (int s=0; s<nrSamples; s++) {
                        rawData[nrProbes][s] = Double.parseDouble(data[sampleIndex[s] + sampleOffset]);
                    }
                    nrProbes++;
                    //if (nrProbes%100==0) System.out.println(nrProbes);
                }
            }
            in.close();
        } catch (Exception e) {
            System.out.println("Error:\t" + e.getMessage());
            e.printStackTrace();
        }
        System.out.println(fileName + "\thas been loaded, nrProbes:\t" + nrProbes + "\tnrSamples:\t" + nrSamples);
    }

    public void recalculateHashMaps() {

        hashProbes.clear();
        for (int probeItr = 0; probeItr < nrProbes; probeItr++) {
            hashProbes.put(probeNames[probeItr], probeItr);
        }

        hashSamples.clear();
        for (int sampleItr = 0; sampleItr < nrSamples; sampleItr++) {
            hashSamples.put(sampleNames[sampleItr], sampleItr);
        }

    }

    public double[][] getRawDataTransposed() {
        double[][] rawDataTransposed = new double[nrSamples][nrProbes];
        for (int s=0; s<nrSamples; s++) {
            for (int p=0; p < nrProbes; p++) {
                rawDataTransposed[s][p] = rawData[p][s];
            }
        }
        return rawDataTransposed;
    }


    public void standardNormalizeData() {

        System.out.println("Setting probe mean to zero and stdev to one for every probe:");
        for (int probeID=0; probeID<nrProbes; probeID++) {
            double vals[] = new double[nrSamples];
            for (int s=0; s<nrSamples; s++) {
                vals[s] = rawData[probeID][s];
            }
            double mean = JSci.maths.ArrayMath.mean(vals);
            for (int s=0; s<nrSamples; s++) {
                vals[s]-=(double) mean;
            }
            double standardDeviation = JSci.maths.ArrayMath.standardDeviation(vals);
            for (int s=0; s<nrSamples; s++) {
                rawData[probeID][s] = (float) (vals[s] / standardDeviation);
            }
        }

    }

    public void save (String fileName) {
        if (fileName.endsWith(".binary")) {

            //Create binary file:
            BufferedOutputStream out = null;
            File fileBinary = new File(fileName + ".dat");
            try {
                out = new BufferedOutputStream(new FileOutputStream(fileBinary));
                out.write(intToByteArray(nrProbes));
                out.write(intToByteArray(nrSamples));
            } catch (IOException e) {
                System.err.println("Can't write to " + fileBinary.getName() + ": " + e.getMessage());
                System.exit(1);
            }
            byte[] buffer = new byte[rawData[0].length * 8];
            for (int row = 0; row < rawData.length; row++) { // rows
                int bufferLoc = 0;
                for (int col = 0; col < rawData[0].length; col++) { // columns
                    long bits = Double.doubleToLongBits(rawData[row][col]);
                    buffer[bufferLoc] = (byte) (bits >> 56);
                    buffer[bufferLoc + 1] = (byte) (bits >> 48 & 0xff);
                    buffer[bufferLoc + 2] = (byte) (bits >> 40 & 0xff);
                    buffer[bufferLoc + 3] = (byte) (bits >> 32 & 0xff);
                    buffer[bufferLoc + 4] = (byte) (bits >> 24 & 0xff);
                    buffer[bufferLoc + 5] = (byte) (bits >> 16 & 0xff);
                    buffer[bufferLoc + 6] = (byte) (bits >> 8 & 0xff);
                    buffer[bufferLoc + 7] = (byte) (bits & 0xff);
                    bufferLoc += 8;
                }
                try {
                    out.write(buffer);
                } catch (IOException e) {
                    System.err.println("Can't write to " + fileBinary.getName() + ": " + e.getMessage());
                    System.exit(1);
                }
            }
            try {
                out.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
            File fileProbes = new File(fileName + ".rows.txt");
            try {
                java.io.BufferedWriter outProbes = new java.io.BufferedWriter(new java.io.FileWriter(fileProbes));
                for (int p=0; p<nrProbes; p++) {
                    outProbes.write(probeNames[p] + "\n");
                }
                outProbes.flush();
                outProbes.close();
            } catch (Exception e) {
                System.out.println("Error:\t" + e.getMessage());
                e.printStackTrace();
            }
            File fileSamples = new File(fileName + ".columns.txt");
            try {
                java.io.BufferedWriter outSamples = new java.io.BufferedWriter(new java.io.FileWriter(fileSamples));
                for (int s=0; s<nrSamples; s++) {
                    outSamples.write(sampleNames[s] + "\n");
                }
                outSamples.flush();
                outSamples.close();
            } catch (Exception e) {
                System.out.println("Error:\t" + e.getMessage());
                e.printStackTrace();
            }


        } else {
            try {
                File fileOut = new File(fileName);
                java.io.BufferedWriter out = new java.io.BufferedWriter(new java.io.FileWriter(fileOut));
                StringBuffer sb = new StringBuffer("-");
                for (int s=0; s<nrSamples; s++) {
                    sb.append("\t" + sampleNames[s]);
                }
                out.write(sb.toString() + "\n");
                for (int p=0; p<nrProbes; p++) {
                    sb = new StringBuffer(probeNames[p]);
                    for (int s=0; s<nrSamples; s++) {
                        sb.append("\t" + rawData[p][s]);
                    }
                    out.write(sb.toString() + "\n");
                }
                out.flush();
                out.close();
            } catch (Exception e) {
                System.out.println("Error:\t" + e.getMessage());
                e.printStackTrace();
            }
        }
    }




    public void transposeDataset() {
        rawData = getRawDataTransposed();

        int nrProbesTemp = nrProbes;
        String[] probeNamesTemp = probeNames;
        HashMap hashProbesTemp = hashProbes;

        nrProbes = nrSamples;
        probeNames = sampleNames;
        hashProbes = hashSamples;

        nrSamples = nrProbesTemp;
        sampleNames = probeNamesTemp;
        hashSamples = hashProbesTemp;

    }

    private byte[] intToByteArray(int value) {
        return new byte[]{(byte) (value >>> 24),
                    (byte) (value >>> 16),
                    (byte) (value >>> 8),
                    (byte) value};
    }

    private int byteArrayToInt(byte[] b) {
        return (b[0] << 24)
                + ((b[1] & 0xff) << 16)
                + ((b[2] & 0xff) << 8)
                + (b[3] & 0xff);
    }


}
