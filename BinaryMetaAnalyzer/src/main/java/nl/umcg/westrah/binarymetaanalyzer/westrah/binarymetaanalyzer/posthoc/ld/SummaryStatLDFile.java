package nl.umcg.westrah.binarymetaanalyzer.westrah.binarymetaanalyzer.posthoc.ld;

import nl.umcg.westrah.binarymetaanalyzer.BinaryMetaAnalysisDataset;
import umcg.genetica.containers.Triple;
import umcg.genetica.io.bin.BinaryFile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

public class SummaryStatLDFile extends BinaryFile {
	
	HashMap<String, Integer> genemap = null;
	ArrayList<String> genelist = new ArrayList<>();
	HashMap<String, Integer> snpmap = null;
	ArrayList<String> snplist = new ArrayList<>();
	
	HashMap<String, Integer> datasetmap = null;
	
	int genectr = 0;
	int snpctr = 0;
	
	public SummaryStatLDFile(String loc, boolean mode) throws IOException {
		super(loc, mode);
	}
	
	public SummaryStatLDFile(String loc, boolean mode, int buffersize) throws IOException {
		super(loc, mode, buffersize);
	}
	
	public SummaryStatLDFile(String loc, boolean mode, int buffersize, boolean useHash) throws IOException {
		super(loc, mode, buffersize, useHash);
	}
	
	public void writeHeader(BinaryMetaAnalysisDataset[] datasets) throws IOException {
		os.writeInt(datasets.length);
		for (BinaryMetaAnalysisDataset d : datasets) {
			os.writeChars(d.getName());
		}
		genemap = new HashMap<>();
		snpmap = new HashMap<>();
		genelist = new ArrayList<>();
		snplist = new ArrayList<>();
	}
	
	public void write(String snp, String gene, float[] data) throws IOException {
		
		Integer snpid = snpmap.get(snp);
		Integer geneid = genemap.get(gene);
		if (snpid == null) {
			snpmap.put(snp, snpctr);
			snplist.add(snp);
			snpid = snpctr;
			snpctr++;
		}
		if (geneid == null) {
			genemap.put(gene, genectr);
			genelist.add(gene);
			geneid = genectr;
			genectr++;
		}
		
		os.writeInt(snpid);
		os.writeInt(geneid);
		os.writeInt(data.length);
		for (Float f : data) {
			os.writeFloat(f);
		}
	}
	
	public String[] readHeader() throws IOException {
		int nrdatasets = is.readInt();
		String[] output = new String[nrdatasets];
		for (int i = 0; i < nrdatasets; i++) {
			output[i] = is.readUTF();
		}
		
		// initialize hashmaps..
		
		return output;
	}
	
	public Triple<Integer, Integer, Float[]> read() throws IOException {
		int snp = is.readInt();
		int gene = is.readInt();
		int nrfloats = is.readInt();
		Float[] data = new Float[nrfloats];
		for (int i = 0; i < data.length; i++) {
			data[i] = is.readFloat();
		}
		return new Triple<Integer, Integer, Float[]>(snp, gene, data);
	}
	
	public void writeHashes(String basefilename) throws IOException {
		BinaryFile bf1 = new BinaryFile(basefilename + "-GeneList.dat", BinaryFile.W, 32 * 1024, false);
		for (String gene : genelist) {
			bf1.writeString(gene);
		}
		bf1.close();
		
		BinaryFile bf2 = new BinaryFile(basefilename + "-SNPList.dat", BinaryFile.W, 32 * 1024, false);
		for (String snp : snplist) {
			bf2.writeString(snp);
		}
		bf2.close();
	}
	
	
	public synchronized void writesync(String snpname, String genename, float[] datasetZ) throws IOException {
		this.write(snpname, genename, datasetZ);
	}
}
