/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.io.trityper.bin;


import umcg.genetica.console.ConsoleGUIElems;
import umcg.genetica.console.ProgressBar;

import java.io.IOException;
import java.util.HashMap;
import java.util.zip.DataFormatException;

/**
 * @author harmjan
 */
public class BinaryResultDataset {

	private String m_name;
	private String m_location;
	private BinaryResultSNP[] snps;
	private HashMap<String, BinaryResultSNP> stringToSNP = new HashMap<String, BinaryResultSNP>();
	private BinaryResultProbe[] probes;
	private HashMap<String, BinaryResultProbe> stringToProbe = new HashMap<String, BinaryResultProbe>();
	private BinaryGZipFloatMatrix bgfm;
	private int maxNrSamples;
	//    private long[] filepointers;
	private float maxfloat = Float.MIN_VALUE;
	private float minfloat = Float.MIN_VALUE;
	private int numprobes;

	public BinaryResultDataset(String location, String name, int permutation) throws IOException {
		m_location = location;
		m_name = name;
		System.out.println("Loading " + name + " from " + location);
		if (permutation == 0) {
			load(m_location + m_name + ".ProbeSummary.dat", m_location + m_name + ".SNPSummary.dat", m_location + m_name + ".ZScoreMatrix.dat");
		} else {
			load(m_location + m_name + ".ProbeSummary.dat", m_location + m_name + "-PermutationRound-" + permutation + ".SNPSummary.dat", m_location + m_name + "-PermutationRound-" + permutation + ".ZScoreMatrix.dat");
		}

	}

	private void load(String probesummaryloc, String snpsummaryloc, String zscoreloc) throws IOException {
		System.out.println("Loading files: \n - " + probesummaryloc + "\n - " + snpsummaryloc + "\n - " + zscoreloc);
		BinaryResultProbeSummary ps = new BinaryResultProbeSummary(probesummaryloc, BinaryResultProbeSummary.R);
		BinaryResultSNPSummary ss = new BinaryResultSNPSummary(snpsummaryloc, BinaryResultSNPSummary.R);

		snps = ss.readAllSNPs();

		probes = ps.readAllProbes();

		for (BinaryResultSNP s : snps) {
			stringToSNP.put(s.getName().intern(), s);
		}
		for (BinaryResultProbe p : probes) {
			stringToProbe.put(p.getName().intern(), p);
		}
		System.out.print("Dataset\t" + m_name + "\n" + ConsoleGUIElems.LINE);
		System.out.println(snps.length + "\t\tSNPs read.");
		System.out.println(probes.length + "\t\tProbes read.");
		System.out.println(ss.getMaxNrSamples() + " samples.");

		this.maxNrSamples = ss.getMaxNrSamples();


		ps.close();
		ss.close();


		bgfm = new BinaryGZipFloatMatrix(zscoreloc, BinaryGZipFloatMatrix.R);
//        checkMatrix();
		numprobes = probes.length;
		System.out.println(ConsoleGUIElems.LINE);
	}

	public void closeMatrix() throws IOException {
		if (bgfm != null) {
			bgfm.close();
			bgfm = null;
		}
	}


	public void openMatrix(int permutation) throws IOException {
		closeMatrix();
	}

	/**
	 * @return the m_name
	 */
	public String getM_name() {
		return m_name;
	}

	/**
	 * @param m_name the m_name to set
	 */
	public void setM_name(String m_name) {
		this.m_name = m_name;
	}

	/**
	 * @return the m_location
	 */
	public String getM_location() {
		return m_location;
	}

	/**
	 * @param m_location the m_location to set
	 */
	public void setM_location(String m_location) {
		this.m_location = m_location;
	}

	/**
	 * @return the snps
	 */
	public BinaryResultSNP[] getSnps() {
		return snps;
	}

	/**
	 * @param snps the snps to set
	 */
	public void setSnps(BinaryResultSNP[] snps) {
		this.snps = snps;
	}

	/**
	 * @return the stringToSNP
	 */
	public HashMap<String, BinaryResultSNP> getStringToSNP() {
		return stringToSNP;
	}

	/**
	 * @param stringToSNP the stringToSNP to set
	 */
	public void setStringToSNP(HashMap<String, BinaryResultSNP> stringToSNP) {
		this.stringToSNP = stringToSNP;
	}

	/**
	 * @return the probes
	 */
	public BinaryResultProbe[] getProbes() {
		return probes;
	}

	/**
	 * @param probes the probes to set
	 */
	public void setProbes(BinaryResultProbe[] probes) {
		this.probes = probes;
	}

	/**
	 * @return the stringToProbe
	 */
	public HashMap<String, BinaryResultProbe> getStringToProbe() {
		return stringToProbe;
	}

	/**
	 * @param stringToProbe the stringToProbe to set
	 */
	public void setStringToProbe(HashMap<String, BinaryResultProbe> stringToProbe) {
		this.stringToProbe = stringToProbe;
	}

	private void checkMatrix() throws IOException {
		System.out.println("Detecting whether binary matrix corresponds to SNP and Probe definition.");
		long expectedsize = (long) snps.length * probes.length;
		System.out.println("Expected matrix size:\t" + expectedsize + " Z-scores");
		System.out.println("Checking matrix: ");
		ProgressBar pb = new ProgressBar(snps.length);
		long count = 0;
		for (int i = 0; i < snps.length; i++) {
			long index = snps[i].getzScoreIndex();
			long next = -1;
			if (i + 1 < snps.length) {
				next = snps[i + 1].getzScoreIndex();
			}

			try {
				bgfm.read(index, next, probes.length);
//                for(int f=0; f<floats.length; f++){
//                    Float fs = floats[f];
//                    if(!Float.isNaN(fs)){
//                        if(fs > maxfloat){
//                            maxfloat = fs;
//                        }
//                        if(fs < minfloat){
//                            minfloat = fs;
//                        }
//                    }
//                }
			} catch (DataFormatException e) {
				System.out.println("");
				e.printStackTrace();
				System.exit(-1);
			}
			pb.iterate();
		}
		pb.close();
		System.out.println("");
		System.out.println("All Probes are present");
		System.out.println("Matrix is OK");
	}

	public synchronized BinaryGZipFloatMatrix getMatrix() {
		return bgfm;
	}

	public int getNumProbes() {
		return numprobes;
	}

	public int getMaxNrSamples() {
		return maxNrSamples;
	}

//    public long getFilePointer(int id) {
//        return filepointers[id];
//    }

	public Float[] readSNPZScores(BinaryResultSNP snp) throws IOException, DataFormatException {
		if (snp == null) {
			return null;
		} else {
			long index = snp.getzScoreIndex();
			long next = -1;
			if (snp.getId() + 1 < snps.length) {
				next = snps[snp.getId() + 1].getzScoreIndex();
			}

			Float[] output = bgfm.read(index, next, probes.length);

			return output;
		}

	}

	/**
	 * @return the maxfloat
	 */
	public float getMaxfloat() {
		return maxfloat;
	}

	/**
	 * @param maxfloat the maxfloat to set
	 */
	public void setMaxfloat(float maxfloat) {
		this.maxfloat = maxfloat;
	}

	/**
	 * @return the minfloat
	 */
	public float getMinfloat() {
		return minfloat;
	}

	/**
	 * @param minfloat the minfloat to set
	 */
	public void setMinfloat(float minfloat) {
		this.minfloat = minfloat;
	}

	public void clearProbeObjects() {
		for (BinaryResultProbe p : probes) {
			p = null;
		}
		probes = null;
	}

	public void close() throws IOException {
		bgfm.close();
	}


}
