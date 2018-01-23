/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.westrah.binarymetaanalyzer;

import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.ByteBuffer;
import java.util.ArrayList;

import umcg.genetica.io.Gpio;
import umcg.genetica.io.bin.BinaryFile;
import umcg.genetica.io.text.TextFile;

/**
 * @author Harm-Jan
 */
public class BinaryMetaAnalysisDataset {
	
	private boolean isCisDataset = false;
	private final String datasetLoc;
	private MetaQTL4MetaTrait[][] snpCisProbeMap;
	
	private long[] snpBytes;
	private String[] alleles;
	private String[] allelesAssessed;
	private String[] minorAlleles;
	private int[] n;
	private double[] callrates;
	private double[] hwes;
	private double[] mafs;
	private String[] probeList;
	
	private String[] snps;
	private final MetaQTL4TraitAnnotation probeAnnotation;
	private final int platformId;
	private RandomAccessFile raf;
	
	private String name = null;
	private String platform = null;
	
	public BinaryMetaAnalysisDataset(String dir, String name, String prefix, int permutation, String platform, MetaQTL4TraitAnnotation probeAnnotation) throws IOException {
		dir = Gpio.formatAsDirectory(dir);
		String matrix = dir;
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
		
		
		BinaryFile f = new BinaryFile(matrix, BinaryFile.R);
		int firstInt = f.readInt();
		f.close();
		isCisDataset = (firstInt == 1);
		System.out.println("Matrix: " + matrix);
		System.out.println("SNPFile: " + snpFile);
		System.out.println("ProbeFile: " + probeFile);
		if (isCisDataset) {
			System.out.println("This dataset is a cis- dataset.");
		} else {
			System.out.println("This dataset is a full size dataset.");
		}
		loadSNPs(snpFile);
		System.out.println(snps.length + " SNPs loaded");
		loadProbes(probeFile);
		System.out.println(probeList.length + " probes loaded");
		
		raf = new RandomAccessFile(matrix, "r");
	}
	
	private void loadSNPs(String snpFile) throws IOException {
		// get nr of lines
		TextFile tf = new TextFile(snpFile, TextFile.R);
		tf.readLine(); // skip header
		int nrSNPs = tf.countLines();
		tf.close();
		
		System.out.println(snpFile + "\t has " + nrSNPs + " SNPs");
		
		snpBytes = new long[nrSNPs];
		alleles = new String[nrSNPs];
		allelesAssessed = new String[nrSNPs];
		minorAlleles = new String[nrSNPs];
		n = new int[nrSNPs];
		callrates = new double[nrSNPs];
		hwes = new double[nrSNPs];
		mafs = new double[nrSNPs];
		
		if (isCisDataset) {
			
			// jagged array, hurrah
			snpCisProbeMap = new MetaQTL4MetaTrait[nrSNPs][0];
		}
		
		tf.open();
		tf.readLine(); // skip header
		String[] elems = tf.readLineElems(TextFile.tab);
		int ln = 0;
		
		snps = new String[nrSNPs];
		snpBytes[0] = 4; // account for magic number.
		while (elems != null) {
			String snp = new String(elems[0].getBytes("UTF-8")).intern();
			String allelesStr = new String(elems[1].getBytes("UTF-8")).intern();
			String minorAlleleStr = new String(elems[2].getBytes("UTF-8")).intern();
			String alleleAssessedStr = new String(elems[3].getBytes("UTF-8")).intern();
			
			snps[ln] = snp;
			alleles[ln] = allelesStr;
			allelesAssessed[ln] = alleleAssessedStr;
			minorAlleles[ln] = minorAlleleStr;
			
			int nrCalled = 0;
			double maf = 0;
			double cr = 0;
			double hwe = 0;
			int nrZScores = 0;
			
			try {
				nrCalled = Integer.parseInt(elems[4]);
			} catch (NumberFormatException e) {
				System.err.println("ERROR: nrCalled is not an int (input: " + elems[4] + ") for dataset: " + datasetLoc + " on line: " + ln);
			}
			try {
				maf = Double.parseDouble(elems[5]);
			} catch (NumberFormatException e) {
				System.err.println("ERROR: maf is not a double (" + elems[5] + ") for dataset: " + datasetLoc + " on line: " + ln);
			}
			try {
				cr = Double.parseDouble(elems[6]);
			} catch (NumberFormatException e) {
				System.err.println("ERROR: cr is not a double (" + elems[6] + ") for dataset: " + datasetLoc + " on line: " + ln);
			}
			try {
				hwe = Double.parseDouble(elems[7]);
			} catch (NumberFormatException e) {
				System.err.println("ERROR: hwe is not a double (" + elems[7] + ") for dataset: " + datasetLoc + " on line: " + ln);
			}
			try {
				nrZScores = Integer.parseInt(elems[8]);
			} catch (NumberFormatException e) {
				System.err.println("ERROR: nrZScores is not an int (input: " + elems[8] + ") for dataset: " + datasetLoc + " on line: " + ln);
			}
			
			n[ln] = nrCalled;
			callrates[ln] = cr;
			hwes[ln] = hwe;
			mafs[ln] = maf;
			
			if (ln + 1 < nrSNPs) {
				snpBytes[ln + 1] = snpBytes[ln] + (nrZScores * 4);
			}
			
			if (isCisDataset) {
				MetaQTL4MetaTrait[] snpProbeList = new MetaQTL4MetaTrait[(elems.length - 9)];
				for (int e = 9; e < elems.length; e++) {
					// get the list of probes for this particular SNP.
					
					String probe = elems[e];
					MetaQTL4MetaTrait t = probeAnnotation.getTraitForPlatformId(platformId, probe);
					// System.out.println(snp+"\t"+elems[e]);
					snpProbeList[e - 9] = t;
				}
				
				snpCisProbeMap[ln] = snpProbeList;
			}
			elems = tf.readLineElems(TextFile.tab);
			ln++;
		}
		tf.close();
		
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
	
	public MetaQTL4MetaTrait[] getCisProbes(int snp) {
		return snpCisProbeMap[snp];
	}
	
	public synchronized float[] getZScores(int snp) throws IOException {
		long snpBytePos = snpBytes[snp];
		
		long snpByteNextPos = 0;
		if (snp == snpBytes.length - 1) {
			snpByteNextPos = raf.length();
		} else {
			snpByteNextPos = snpBytes[snp + 1];
		}

//        
//        if (snpCisProbeMap != null) {
//            nrZ = snpCisProbeMap[snp].length;
//        }
//        System.out.println(snp + "\t" + snpBytePos + "\t" + snpByteNextPos);
		raf.seek(snpBytePos);
		int readlen = (int) (snpByteNextPos - snpBytePos);
		byte[] bytesToRead = new byte[readlen];
		raf.read(bytesToRead);
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
	
}
