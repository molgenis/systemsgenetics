/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.westrah.binarymetaanalyzer;

import umcg.genetica.io.Gpio;
import umcg.genetica.io.bin.BinaryFile;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;
import umcg.genetica.util.Primitives;

import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.ByteBuffer;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * @author Harm-Jan
 */
public class BinaryMetaAnalysisDataset {

    private final int permutation;
    private final AtomicInteger snpctr;
    private boolean isCisDataset = false;
    private final String datasetLoc;
    private MetaQTL4MetaTrait[][] snpCisProbeMap;

    private long[] snpBytes;
    private String[] alleles;
    private String[] allelesAssessed;
    private String[] minorAlleles;
    private int[] n;
    private float[] callrates;
    private float[] hwes;
    private float[] mafs;
    private String[] probeList;

    private String[] snps;
    private final MetaQTL4TraitAnnotation probeAnnotation;
    private final int platformId;
    private RandomAccessFile raf;
    private HashMap<String, Double> featureOccuranceScaleMap = null;

    private String name = null;
    private String platform = null;

    public BinaryMetaAnalysisDataset(String dir, String name, String prefix, int permutation, String platform,
                                     MetaQTL4TraitAnnotation probeAnnotation, String featureOccuranceScaleMapFile,
                                     boolean loadsnpstats, AtomicInteger c) throws IOException {
        dir = Gpio.formatAsDirectory(dir);
        String matrix = dir;
        this.snpctr = c;
        String probeFile = dir;
        String snpFile = dir;
        this.platform = platform;
        this.name = name;
        this.probeAnnotation = probeAnnotation;
        this.platformId = probeAnnotation.getPlatformId(platform);
        String pref = "Dataset";
        if (prefix != null) {
            pref = prefix;
        }
        if (permutation > 0) {
            matrix += pref + "-PermutationRound-" + permutation + ".dat";
            probeFile += pref + "-PermutationRound-" + permutation + "-ColNames.txt.gz";
            snpFile += pref + "-PermutationRound-" + permutation + "-RowNames.txt.gz";
        } else {
            matrix += pref + ".dat";
            probeFile += pref + "-ColNames.txt.gz";
            snpFile += pref + "-RowNames.txt.gz";
        }

        this.datasetLoc = dir;
        // check presence of files
        if (!Gpio.exists(matrix)) {
            throw new IOException("Could not find file: " + matrix);
        }
        if (!Gpio.exists(probeFile)) {
            throw new IOException("Could not find file: " + probeFile);
        }
        if (!Gpio.exists(snpFile)) {
            throw new IOException("Could not find file: " + snpFile);
        }

        // TODO: if dataset is < 2gb, then just buffer the whole thing in memory.

        this.permutation = permutation;
        BinaryFile f = new BinaryFile(matrix, BinaryFile.R);
        int firstInt = f.readInt();
        f.close();
        isCisDataset = (firstInt == 1);
        System.out.println(name+"\tMatrix: " + matrix);
        System.out.println(name+"\tSNPFile: " + snpFile);
        System.out.println(name+"\tProbeFile: " + probeFile);
        if (isCisDataset) {
            System.out.println(name+"\tThis dataset is a CIS dataset.");
        } else {
            System.out.println(name+"\tThis dataset is a full size dataset.");
        }
        loadProbes(probeFile);

        System.out.println(name+"\tPermutation: " + permutation + "\t" + probeList.length + " probes loaded");

        loadSNPs(snpFile, loadsnpstats);
        if (snpctr == null) {
            System.out.println(name+"\tPermutation: " + permutation + "\t" + snps.length + " SNPs loaded");
        }
        if (featureOccuranceScaleMapFile != null) {
            TextFile readScaleInfo = new TextFile(featureOccuranceScaleMapFile, TextFile.R);
            HashMap<String, String> temporaryMap = new HashMap<>(readScaleInfo.readAsHashMap(0, 1));
            readScaleInfo.close();
            featureOccuranceScaleMap = new HashMap<>();
            for (String k : temporaryMap.keySet()) {
                featureOccuranceScaleMap.put(k, Double.parseDouble(temporaryMap.get(k)));
            }
        }

        raf = new RandomAccessFile(matrix, "r");
        long nrBytesPerBuffer = 1048576 * 2024; // 2Gb
        if (nrBytesPerBuffer > raf.length()) {
            nrBytesPerBuffer = raf.length();
        }

        mappedzscorehandle = raf.getChannel().map(FileChannel.MapMode.READ_ONLY, 0, nrBytesPerBuffer);
        mappedzscorehandle.load();
        bZs = new byte[(int) nrBytesPerBuffer];
        mappedzscorehandle.get(bZs);

        currentmapstart = 0;
        currentmapend = currentmapstart + nrBytesPerBuffer;
        System.out.println(name+"\tPermutation: " + permutation + "\t" + "File size: " + raf.length() + "\tNew buffer: " + bZs.length + "\tsta: " + currentmapstart + "\tsto: " + currentmapend +
                "\tlen: " + (currentmapend - currentmapstart) + "\tnrbytes: " + nrBytesPerBuffer);
    }

    private void loadSNPs(String snpFile, boolean loadstats) throws IOException {
        // get nr of lines
        ArrayList<String> snpsAl = new ArrayList<>(10000000);
        ArrayList<Long> snpbytesAl = new ArrayList<>(10000000);
        ArrayList<String> allelesAl = new ArrayList<>(10000000);
        ArrayList<String> allelesAssessedAl = new ArrayList<>(10000000);
        ArrayList<String> minorAllelesAl = new ArrayList<>(10000000);

        ArrayList<MetaQTL4MetaTrait[]> snpCisProbeMapAl = new ArrayList<>(10000000);

//
//		tf.readLine(); // skip header
//		int nrSNPs = tf.countLines();
//		tf.close();
        ArrayList<Float> mafsAl = null;
        ArrayList<Integer> nAl = null;
        ArrayList<Float> rCallRatesAl = null;
        ArrayList<Float> hwesAl = null;
        if (loadstats) {
            mafsAl = new ArrayList<>(10000000);
            nAl = new ArrayList<>(10000000);
            rCallRatesAl = new ArrayList<>(10000000);
            hwesAl = new ArrayList<>(10000000);
        }
        TextFile tf = new TextFile(snpFile, TextFile.R, 1048576);
        tf.readLine(); // skip header
        String[] elems = tf.readLineElems(TextFile.tab);
        int ln = 0;

        Map<String, MetaQTL4MetaTrait> traithash = probeAnnotation.getTraitHashForPlatform(platformId);

        snpbytesAl.add(4l);
        while (elems != null) {

            String snp = Strings.cache(elems[0]);//new String(elems[0].getBytes("UTF-8")).intern();
            String allelesStr = Strings.cache(elems[1]); //new String(elems[1].getBytes("UTF-8")).intern();
            String minorAlleleStr = Strings.cache(elems[2]);//new String(elems[2].getBytes("UTF-8")).intern();
            String alleleAssessedStr = Strings.cache(elems[3]); //new String(elems[3].getBytes("UTF-8")).intern();

            snpsAl.add(snp);
            allelesAl.add(allelesStr);
            allelesAssessedAl.add(alleleAssessedStr);
            minorAllelesAl.add(minorAlleleStr);
//			snps[ln] = snp;
//			alleles[ln] = allelesStr;
//			allelesAssessed[ln] = alleleAssessedStr;
//			minorAlleles[ln] = minorAlleleStr;

            int nrCalled = 0;
            float maf = 0;
            float cr = 0;
            float hwe = 0;
            int nrZScores = 0;

//			if (loadstats) {
            try {
                nrCalled = Integer.parseInt(elems[4]);
            } catch (NumberFormatException e) {
                System.err.println(name+"\tERROR: nrCalled is not an int (input: " + elems[4] + ") for dataset: " + datasetLoc + " on line: " + ln);
            }
            try {
                maf = Float.parseFloat(elems[5]);
            } catch (NumberFormatException e) {
                System.err.println(name+"\tERROR: maf is not a double (" + elems[5] + ") for dataset: " + datasetLoc + " on line: " + ln);
            }
            try {
                cr = Float.parseFloat(elems[6]);
            } catch (NumberFormatException e) {
                System.err.println(name+"\tERROR: cr is not a double (" + elems[6] + ") for dataset: " + datasetLoc + " on line: " + ln);
            }
            try {
                hwe = Float.parseFloat(elems[7]);
            } catch (NumberFormatException e) {
                System.err.println(name+"\tERROR: hwe is not a double (" + elems[7] + ") for dataset: " + datasetLoc + " on line: " + ln);
            }
            hwesAl.add(hwe);
            mafsAl.add(maf);
            rCallRatesAl.add(cr);
            nAl.add(nrCalled);
//			}
            try {
                nrZScores = Integer.parseInt(elems[8]);
            } catch (NumberFormatException e) {
                System.err.println(name+"\tERROR: nrZScores is not an int (input: " + elems[8] + ") for dataset: " + datasetLoc + " on line: " + ln);
            }


            snpbytesAl.add(snpbytesAl.get(ln) + (nrZScores * 4));

//			n[ln] = nrCalled;
//			callrates[ln] = cr;
//			hwes[ln] = hwe;
//			mafs[ln] = maf;

//			if (ln + 1 < nrSNPs) {
//
//				snpBytes[ln + 1] = snpBytes[ln] + (nrZScores * 4);
//			}

            if (isCisDataset) {
                MetaQTL4MetaTrait[] snpProbeList = new MetaQTL4MetaTrait[(elems.length - 9)];
                for (int e = 9; e < elems.length; e++) {
                    // get the list of probes for this particular SNP.
                    String probe = elems[e]; // Strings.cache(elems[e]); // new String(elems[e].getBytes("UTF-8")).intern();
                    MetaQTL4MetaTrait t = traithash.get(probe);
                    // System.out.println(snp+"\t"+elems[e]);
                    snpProbeList[e - 9] = t;
                }

                snpCisProbeMapAl.add(snpProbeList);
            }
            elems = tf.readLineElems(TextFile.tab);

            ln++;
            if (ln % 1000000 == 0) {
                System.out.print(name+"\tPermutation: " + permutation + "\t" + ln + " lines parsed sofar.\r");
            }
        }
        tf.close();
        if (snpctr != null) {
            snpctr.getAndSet(ln);
        }
        System.out.println();
        // copy data.
        snps = snpsAl.toArray(new String[0]);
        snpBytes = Primitives.toPrimitiveArr(snpbytesAl.toArray(new Long[0]));
        alleles = allelesAl.toArray(new String[0]);
        allelesAssessed = allelesAssessedAl.toArray(new String[0]);
        minorAlleles = minorAllelesAl.toArray(new String[0]);

//		if (loadstats) {
        n = Primitives.toPrimitiveArr(nAl.toArray(new Integer[0]));
        callrates = Primitives.toPrimitiveArr(rCallRatesAl.toArray(new Float[0]));
        hwes = Primitives.toPrimitiveArr(hwesAl.toArray(new Float[0]));
        mafs = Primitives.toPrimitiveArr(mafsAl.toArray(new Float[0]));
//		}
        if (isCisDataset) {
            // jagged array, hurrah
            snpCisProbeMap = new MetaQTL4MetaTrait[snps.length][0];
            for (int s = 0; s < snps.length; s++) {
                snpCisProbeMap[s] = snpCisProbeMapAl.get(s);
            }
//
        }

    }

    private void loadProbes(String columns) throws IOException {
        TextFile tf = new TextFile(columns, TextFile.R);
//        nrProbes = tf.countLines();
//        tf.close();
//
//        tf.open();
        ArrayList<String> allProbes = tf.readAsArrayList();
        tf.close();
        probeList = allProbes.toArray(new String[0]);
    }

    public  MetaQTL4MetaTrait[] getCisProbes(int snp) {
        return snpCisProbeMap[snp];
    }


    private MappedByteBuffer mappedzscorehandle = null;
    long currentmapend;
    byte[] bZs = null;
    long currentmapstart;

    public synchronized float[] getZScores(int snp) throws IOException {
        long snpBytePos = snpBytes[snp];

        long snpByteNextPos = 0;
        if (snp == snpBytes.length - 1) {
            snpByteNextPos = raf.length();
        } else {
            snpByteNextPos = snpBytes[snp + 1];
        }

        int readlen = (int) (snpByteNextPos - snpBytePos);

//		System.out.println("snp pos: " + snpBytePos);

        if (mappedzscorehandle == null || snpBytePos > currentmapend || snpByteNextPos > currentmapend || snpBytePos < currentmapstart) {

            long nrBytesPerBuffer = 1048576 * 10; // 256 mb

            if (snpBytePos + nrBytesPerBuffer > raf.length()) {
                nrBytesPerBuffer = raf.length() - snpBytePos;
            }
            mappedzscorehandle = raf.getChannel().map(FileChannel.MapMode.READ_ONLY, snpBytePos, nrBytesPerBuffer);
            mappedzscorehandle.load();
            bZs = new byte[(int) nrBytesPerBuffer];
            mappedzscorehandle.get(bZs);

            currentmapstart = snpBytePos;
            currentmapend = currentmapstart + nrBytesPerBuffer;
            System.out.println(name+"\tNew buffer: " + bZs.length + "\tsta: " + currentmapstart + "\tsto: " + currentmapend +
                    "\tlen: " + (currentmapend - currentmapstart) + "\tnrbytes: " + nrBytesPerBuffer);
        }

        int offset = (int) (snpBytePos - currentmapstart);
//		System.out.println("byte pos: " + snpBytePos + "\toffset " + offset + "\tlen " + readlen);
        byte[] bytesToRead = new byte[readlen];

        System.arraycopy(bZs, offset, bytesToRead, 0, readlen);
//
////
////        if (snpCisProbeMap != null) {
////            nrZ = snpCisProbeMap[snp].length;
////        }
////        System.out.println(snp + "\t" + snpBytePos + "\t" + snpByteNextPos);
//		raf.seek(snpBytePos);
//

//		raf.read(bytesToRead);
        ByteBuffer bytebuffer = ByteBuffer.wrap(bytesToRead);
        float[] output = new float[readlen / 4];
        for (int i = 0; i < output.length; i++) {
            output[i] = bytebuffer.getFloat();
        }
        
        
        return output;
    }

    public String[] getSNPs() {
        return snps;
    }

    public String[] getProbeList() {
        return probeList;
    }

    public int getSampleSize(int datasetSNPId) {
        return n[datasetSNPId];
    }


    public String getAlleles(int datasetSNPId) {
        return alleles[datasetSNPId];
    }

    public String getAlleleAssessed(int datasetSNPId) {
        return allelesAssessed[datasetSNPId];
    }

    public boolean getIsCisDataset() {
        return isCisDataset;
    }

    public void close() throws IOException {
        raf.close();

    }

    public String getName() {
        return name;
    }

    public String getPlatform() {
        return platform;
    }

    public HashMap<String, Double> getFeatureOccuranceScaleMap() {
        return featureOccuranceScaleMap;
    }

    public int getPermutation() {

        return permutation;
    }
}
