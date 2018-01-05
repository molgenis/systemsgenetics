package nl.umcg.westrah.binarymetaanalyzer;

import gnu.trove.map.hash.THashMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import org.apache.commons.collections.primitives.ArrayIntList;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.io.bin.BinaryFile;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.stats.ZScores;

import javax.xml.bind.annotation.adapters.HexBinaryAdapter;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.regex.Pattern;

public class InternalMetaAnalysisTask implements Runnable {
	
	
	private final MultiThreadProgressBar mtPb;
	private final int taskid;
	private int[] snpIndex;
	private InternalMetaAnalysisDataset dataset;
	private final InternalMetaAnalysisSettings settings;
	private final THashMap<String, ArrayList<String>> traitMap;
	private final ArrayList<String> outputEntries;
	
	private final LinkedHashMap<String, Integer> traitLocationMap;
	private String[] snpList;
	private TextFile zScoreRowNamesFile;
	private BinaryFile zScoreBinaryFile;
	int permutation;
	String outdir;
	
	public InternalMetaAnalysisTask(
			InternalMetaAnalysisSettings settings,
			THashMap<String, ArrayList<String>> traitMap,
			ArrayList<String> outputEntries,
			LinkedHashMap<String, Integer> traitLocationMap,
			int permutation,
			String outdir,
			MultiThreadProgressBar mtPb,
			int taskid
	) {
		this.settings = settings;
		this.traitMap = traitMap;
		this.outputEntries = outputEntries;
		this.traitLocationMap = traitLocationMap;
		this.permutation = permutation;
		this.outdir = outdir;
		this.mtPb = mtPb;
		this.taskid = taskid;
	}
	
	
	@Override
	public void run() {
		
		try {
			
			boolean runningPermutation = true;
			if (permutation == 0) {
				runningPermutation = false;
			}
			
			//Initialize the new binaryOutput
			//Original initialize binary matrix
			String fileName = settings.getOutput() + settings.getDatasetname();
			
			if (runningPermutation) {
				fileName += "-PermutationRound-" + permutation;
			}
			
			boolean useHash = false;
			zScoreBinaryFile = new BinaryFile(fileName + ".dat", BinaryFile.W, 1048576 * 50, useHash); // 50mb output buffer
			
			System.out.println("Loading dataset");
			dataset = new InternalMetaAnalysisDataset(settings.getDatasetlocation(),
					settings.getDatasetname(),
					settings.getDatasetPrefix(),
					permutation);
			
			System.out.println("Loaded");
			
			// create meta-analysis SNP index. have to recreate this every permutation,
			// since the order of SNPs is generated at random.
			System.out.println("Creating SNP index");
			createSNPIndex(outdir);
			System.out.println("Total of " + snpIndex.length + " SNPs");
			
			// write magic number
			if (dataset.getIsCisDataset()) {
				zScoreBinaryFile.writeInt(1);
			} else {
				zScoreBinaryFile.writeInt(0);
			}
			
			zScoreRowNamesFile = new TextFile(fileName + "-RowNames.txt.gz", TextFile.W);
			zScoreRowNamesFile.writeln("SNP\tAlleles\tMinorAllele\tAlleleAssessed\tNrCalled\tMaf\tHWE\tCallRate");
			
			TextFile tf = new TextFile(fileName + "-ColNames.txt.gz", TextFile.W);
			tf.writeList(outputEntries);
			tf.close();
			
			// run analysis
			mtPb.setSubtasks(taskid, snpList.length);
			for (int snp = 0; snp < snpList.length; snp++) {
				metaAnalyze(snp);
				mtPb.iterate(taskid);
			}
			
			zScoreBinaryFile.close();
			zScoreRowNamesFile.close();
			
			fileName = "check";
			if (runningPermutation) {
				fileName += "-PermutationRound-" + permutation;
			}
			if (useHash) {
				//Close binaryOutput & Create binary check file.
				fileName += ".md5";
				HexBinaryAdapter md5Parser = new HexBinaryAdapter();
				BufferedWriter md5writer = new BufferedWriter(new FileWriter(settings.getOutput() + fileName));
				md5writer.write(md5Parser.marshal(zScoreBinaryFile.getWrittenHash()) + "  " + fileName + '\n');
				md5writer.close();
			} else {
				// lame attempt to check whether process is completed.
				fileName += ".txt";
				TextFile check = new TextFile(settings.getOutput() + fileName, TextFile.W);
				check.writeln("Number of snps analyzed: " + snpList.length);
				check.close();
			}
			
			mtPb.complete(taskid);
			
		} catch (IOException e) {
			
			e.printStackTrace();
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		
	}
	
	
	private void metaAnalyze(int snp) throws IOException {
		// Here now. Need to check what probes to meta-analyze, do the meta-analysis and write.
		// do we need to check if alleles are different in the same dataset?
		// get ZScores for this SNP, no matter what.
		// get list of probes to test
		int datasetSNPId = snpIndex[snp];
		
		if (datasetSNPId != -9) { // -9 means: snp not available
			float[] datasetZScores = null;
			try {
				datasetZScores = dataset.getZScores(datasetSNPId);
			} catch (IOException e) {
				e.printStackTrace();
			}
			THashMap<String, DescriptiveStatistics> remappedEntries = new THashMap<String, DescriptiveStatistics>();
			THashMap<String, SummaryStatistics> remappedEntriesAbs = new THashMap<String, SummaryStatistics>();
//                    System.out.println(dataset.getSNPs()[datasetSNPId]);
//                    System.out.println(datasetZScores.length);
			
			// count NaNs
			
			
			if (datasetZScores != null && datasetZScores.length > 0) {
				int nrNaN = 0;
				for (int c = 0; c < datasetZScores.length; c++) {
					if (Float.isNaN(datasetZScores[c])) {
						nrNaN++;
					}
				}
				// no need to run analysis on a snp that has no data associated with it.
				if (nrNaN != datasetZScores.length) {
					if (dataset.getIsCisDataset()) {
						// this requires us to retrieve the z-scores differently
						// we need to figure out which probes match up, but their orders might be different
						// and the number of probes tested in each dataset might differ as well
						
						// get the probes tested against the SNP
						String[] datasetCisProbes = dataset.getCisProbes(datasetSNPId);
//                        System.out.println(datasetCisProbes.length);
						for (int i = 0; i < datasetCisProbes.length; i++) {
							String p = datasetCisProbes[i];
							if (traitMap.containsKey(p)) {
								for (String feature : traitMap.get(p)) {
									if (!remappedEntries.containsKey(feature)) {
										remappedEntries.put(feature, new DescriptiveStatistics());
										remappedEntriesAbs.put(feature, new SummaryStatistics());
									}
									remappedEntries.get(feature).addValue(datasetZScores[i]);
									remappedEntriesAbs.get(feature).addValue(Math.abs(datasetZScores[i]));
//                                System.out.print(p+"\t"+datasetZScores[i]);
								}
							}
						}
						if (!remappedEntries.isEmpty()) {
							double[] zScoresOut = new double[remappedEntries.size()];
							// initialize with NaN
							for (int q = 0; q < zScoresOut.length; q++) {
								zScoresOut[q] = Double.NaN;
							}
							
							String[] keysOut = new String[remappedEntries.size()];
							int counter = 0;
							for (Map.Entry<String, DescriptiveStatistics> e : remappedEntries.entrySet()) {
								keysOut[counter] = e.getKey();
								if (e.getValue().getN() == 0) {
									zScoresOut[counter] = Double.NaN;
								} else if (e.getValue().getN() == 1) {
									zScoresOut[counter] = e.getValue().getValues()[0];
								} else {
									//Now we need to meta-analyze it.
									if (settings.getzScoreMergeOption().equals("weightedzscore")) {
										//Need to get sample size
										ArrayIntList sampleSizes = new ArrayIntList();
										for (int sC = 0; sC < e.getValue().getN(); sC++) {
											sampleSizes.add(dataset.getSampleSize(datasetSNPId));
										}
										zScoresOut[counter] = ZScores.getWeightedZ(e.getValue().getValues(), sampleSizes.toArray());
									} else if (settings.getzScoreMergeOption().equals("mean")) {
										zScoresOut[counter] = e.getValue().getMean();
									} else if (settings.getzScoreMergeOption().equals("median")) {
										zScoresOut[counter] = e.getValue().getPercentile(50);
									} else if (settings.getzScoreMergeOption().equals("min")) {
										zScoresOut[counter] = remappedEntriesAbs.get(e.getKey()).getMin();
									} else if (settings.getzScoreMergeOption().equals("max")) {
										zScoresOut[counter] = remappedEntriesAbs.get(e.getKey()).getMax();
									} else {
										System.out.println("Not supported merging.");
										System.exit(0);
									}
									
								}
								counter++;
							}
							//                        writeBinaryResult(String snpname, double hwe, double cr, double maf, int numberCalled, String alleles, String minorAllele, String alleleassessed, double[] datasetZScores, String[] probeNames, BinaryFile outfile, TextFile snpfile)
							
							writeBinaryResult(dataset.getSNPs()[datasetSNPId],
									dataset.getHwes()[datasetSNPId],
									dataset.getCallrates()[datasetSNPId],
									dataset.getMafs()[datasetSNPId],
									dataset.getN()[datasetSNPId],
									dataset.getAlleles()[datasetSNPId],
									dataset.getMinorAlleles()[datasetSNPId],
									dataset.getAlleleAssessed(datasetSNPId),
									zScoresOut,
									keysOut);
							
						}
					} else { // this is not a cis dataset
						// use the full probe index
						
						// probeIndex[t.getMetaTraitId()][d] = p;
						for (int i = 0; i < dataset.getProbeList().length; i++) {
							String p = dataset.getProbeList()[i];
							if (traitMap.containsKey(p)) {
								for (String feature : traitMap.get(p)) {
									if (!remappedEntries.containsKey(feature)) {
										remappedEntries.put(feature, new DescriptiveStatistics());
										remappedEntriesAbs.put(feature, new SummaryStatistics());
									}
									remappedEntries.get(feature).addValue(datasetZScores[i]);
									remappedEntriesAbs.get(feature).addValue(datasetZScores[i]);
								}
							}
						}
						if (!remappedEntries.isEmpty()) {
							//This is almost idenitcal to the cis-data.
							//But here we need to make sure the remapped data is stored in the correct relative to the file description. i.e. outputEntries
							
							//Not observed zScores are kept at 0.
							double[] zScoresOut = new double[outputEntries.size()];
							for (int q = 0; q < zScoresOut.length; q++) {
								zScoresOut[q] = Double.NaN;
							}
							
							for (Map.Entry<String, DescriptiveStatistics> e : remappedEntries.entrySet()) {
								
								int arrayLoc = traitLocationMap.get(e.getKey());
								
								if (e.getValue().getN() == 0) {
									zScoresOut[arrayLoc] = Double.NaN;
								} else if (e.getValue().getN() == 1) {
									zScoresOut[arrayLoc] = e.getValue().getValues()[0];
								} else {
									//Now we need to meta-analyze it.
									if (settings.getzScoreMergeOption().equals("weightedzscore")) {
										//Need to get sample size
										ArrayIntList sampleSizes = new ArrayIntList();
										for (int sC = 0; sC < e.getValue().getN(); sC++) {
											sampleSizes.add(dataset.getSampleSize(datasetSNPId));
										}
										zScoresOut[arrayLoc] = ZScores.getWeightedZ(e.getValue().getValues(), sampleSizes.toArray());
									} else if (settings.getzScoreMergeOption().equals("mean")) {
										zScoresOut[arrayLoc] = e.getValue().getMean();
									} else if (settings.getzScoreMergeOption().equals("median")) {
										zScoresOut[arrayLoc] = e.getValue().getPercentile(50);
									} else if (settings.getzScoreMergeOption().equals("min")) {
										zScoresOut[arrayLoc] = remappedEntriesAbs.get(e.getKey()).getMin();
									} else if (settings.getzScoreMergeOption().equals("max")) {
										zScoresOut[arrayLoc] = remappedEntriesAbs.get(e.getKey()).getMax();
									} else {
										System.out.println("Not supported merging.");
										System.exit(0);
									}
									
								}
							}
							//                        writeBinaryResult(String snpname, double hwe, double cr, double maf, int numberCalled, String alleles, String minorAllele, String alleleassessed, double[] datasetZScores, String[] probeNames, BinaryFile outfile, TextFile snpfile)
							
							writeBinaryResult(dataset.getSNPs()[datasetSNPId],
									dataset.getHwes()[datasetSNPId],
									dataset.getCallrates()[datasetSNPId],
									dataset.getMafs()[datasetSNPId],
									dataset.getN()[datasetSNPId],
									dataset.getAlleles()[datasetSNPId],
									dataset.getMinorAlleles()[datasetSNPId],
									dataset.getAlleleAssessed(datasetSNPId),
									zScoresOut,
									null
							);
							
						}
					}
				}
				
			}
			
		}
	}
	
	
	private void createSNPIndex(String outdir) throws IOException {
		// create a list of all available SNPs
		ArrayList<String> allSNPs = new ArrayList<String>();
		String[] snps = dataset.getSNPs();
		{
			HashSet<String> visitedSNPs = new HashSet<>();
			for (String snp : snps) {
				if (!visitedSNPs.contains(snp)) {
					allSNPs.add(snp);
				}
			}
		}
		
		// create a temporary map that maps each SNP to a meta-analysis position
		int ctr = 0;
		TObjectIntHashMap<String> snpMap = new TObjectIntHashMap<String>(allSNPs.size(), 0.85f, -9);
		snpList = new String[allSNPs.size()];
		for (String s : allSNPs) {
			snpMap.put(s, ctr);
			snpList[ctr] = s;
			ctr++;
		}
		
		// fill index
		snpIndex = new int[allSNPs.size()];
		for (int s = 0; s < allSNPs.size(); s++) {
			snpIndex[s] = -9;
		}
		
		for (int s = 0; s < snps.length; s++) {
			String snp = snps[s];
			int id = snpMap.get(snp);
			if (id != -9) {
				snpIndex[id] = s;
			}
		}

//        TextFile tf = new TextFile(outdir + "snpindex.txt", TextFile.W);
//        String header = "metaID";
//        header += "\t" + dataset.getName() + "-sid";
//        tf.writeln(header);
//
//        for (int s = 0; s < snpList.length; s++) {
//            String ln = snpList[s];
//            ln += "\t" + snpIndex[s];
//            tf.writeln(ln);
//        }
//        tf.close();
	}
	
	private void writeBinaryResult(String snpname,
								   double hwe,
								   double cr,
								   double maf,
								   int numberCalled,
								   String alleles,
								   String minorAllele,
								   String alleleassessed,
								   double[] datasetZScores,
								   String[] probeNames
	
	) throws IOException {
		StringBuilder sb = null;
		for (int p = 0; p < datasetZScores.length; p++) {
			float z = (float) datasetZScores[p];
			// System.out.println(p + "\t" + alleleassessed + "\t" + m_probeList[p] + "\t" + z + "\t" + currentWP.getFlipSNPAlleles()[d]);
			if (probeNames != null) {
				// take into account that not all probes have been tested..
				String probeName = probeNames[p];
				zScoreBinaryFile.writeFloat(z);
				if (sb == null) {
					sb = new StringBuilder();
				} else {
					sb.append("\t");
				}
				sb.append(probeName);
			} else {
				zScoreBinaryFile.writeFloat(z);
			}
		}
		
		StringBuilder buffer = new StringBuilder();
		buffer.append(snpname)
				.append("\t")
				.append(alleles)
				.append("\t")
				.append(minorAllele)
				.append("\t")
				.append(alleleassessed)
				.append("\t")
				.append(numberCalled)
				.append("\t")
				.append(maf)
				.append("\t")
				.append(hwe)
				.append("\t")
				.append(cr)
				.append("\t")
				.append(datasetZScores.length)
				.append("\t");
		if (sb != null) {
			buffer.append(sb.toString());
		} else {
			buffer.append("-");
		}
		
		zScoreRowNamesFile.writeln(buffer.toString());
	}
	
}
