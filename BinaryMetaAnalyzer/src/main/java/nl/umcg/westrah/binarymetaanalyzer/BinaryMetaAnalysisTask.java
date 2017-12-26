package nl.umcg.westrah.binarymetaanalyzer;

import gnu.trove.map.hash.TObjectIntHashMap;
import umcg.genetica.containers.Triple;
import umcg.genetica.io.trityper.util.BaseAnnot;
import umcg.genetica.math.stats.Descriptives;
import umcg.genetica.math.stats.ZScores;
import umcg.genetica.text.Strings;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Set;
import java.util.concurrent.Callable;

public class BinaryMetaAnalysisTask implements Callable<Triple<ArrayList<QTL>, String, String>> {
	
	private final BinaryMetaAnalysisSettings settings;
	private final MetaQTL4TraitAnnotation probeAnnotation;
	private final BinaryMetaAnalysisDataset[] datasets;
	private final int[][] snpIndex;
	private final String[] snpList;
	private final String[] snpChr;
	private final int[] snpPositions;
	private final Integer[][] probeIndex;
	private final MetaQTL4MetaTrait[][] snpprobeCombos;
	private final TObjectIntHashMap<MetaQTL4MetaTrait> traitMap;
	private final MetaQTL4MetaTrait[] traitList;
	private final int snp;
	
	public BinaryMetaAnalysisTask(BinaryMetaAnalysisSettings settings,
								  MetaQTL4TraitAnnotation probeAnnotation,
								  BinaryMetaAnalysisDataset[] datasets,
								  int[][] snpIndex,
								  String[] snpList,
								  String[] snpChr,
								  int[] snpPositions,
								  Integer[][] probeIndex,
								  MetaQTL4MetaTrait[][] snpprobeCombos,
	
								  TObjectIntHashMap<MetaQTL4MetaTrait> traitMap,
	
								  MetaQTL4MetaTrait[] traitList,
								  int snp) {
		this.settings = settings;
		this.probeAnnotation = probeAnnotation;
		this.datasets = datasets;
		this.snpIndex = snpIndex;
		this.snpList = snpList;
		this.snpChr = snpChr;
		this.snpPositions = snpPositions;
		this.probeIndex = probeIndex;
		this.snpprobeCombos = snpprobeCombos;
		this.traitMap = traitMap;
		this.traitList = traitList;
		this.snp = snp;
	}
	
	@Override
	public Triple<ArrayList<QTL>, String, String> call() throws Exception {
		double[] zscoretableoutput = null;
		int[] zscorenrsamplestableoutput = null;
		if (settings.isMakezscoretable()) {
			
			if (zscoretableoutput != null) {
				for (int i = 0; i < zscoretableoutput.length; i++) {
					zscoretableoutput[i] = Double.NaN;
					zscorenrsamplestableoutput[i] = 0;
				}
			} else {
				zscoretableoutput = new double[probeIndex.length];
				zscorenrsamplestableoutput = new int[probeIndex.length];
			}
		}
		
		boolean printed = false;
		int[] sampleSizes = new int[datasets.length];
		Boolean[] flipZScores = new Boolean[datasets.length];
		String alleles = null;
		String alleleAssessed = null;
		
		ArrayList<QTL> qtlOutput = new ArrayList<>();
		Triple<ArrayList<QTL>, String, String> output = new Triple<ArrayList<QTL>, String, String>(qtlOutput, null, null);
		// determine whether to flip the alleles for a certain dataset
		for (int d = 0; d < datasets.length; d++) {
			int datasetSNPId = snpIndex[snp][d];
			if (datasetSNPId != -9) {
				sampleSizes[d] = datasets[d].getSampleSize(datasetSNPId);
				if (alleles == null) {
					flipZScores[d] = false;
					alleles = datasets[d].getAlleles(datasetSNPId);
					alleleAssessed = datasets[d].getAlleleAssessed(datasetSNPId);
				} else {
					String alleles2 = datasets[d].getAlleles(datasetSNPId);
					String alleleAssessed2 = datasets[d].getAlleleAssessed(datasetSNPId);
					flipZScores[d] = BaseAnnot.flipalleles(alleles, alleleAssessed, alleles2, alleleAssessed2);
				}
			}
		}
		
		// get ZScores for this SNP
		// get list of probes to test
		if (settings.getAnalysisType().equals(BinaryMetaAnalysisSettings.Analysis.CIS) || snpprobeCombos != null) {
			// do cis stuff, or stuff to specific sets of probes...
			
			// get all the possible traits near the SNP
//                    Set<MetaQTL4MetaTrait> cisProbesForSNP = probeAnnotation.getMetatraits().getTraitInWindow(snpChr[snp], snpPositions[snp], settings.getCisdistance());
//                    MetaQTL4MetaTrait[] cisProbeArray = cisProbesForSNP.toArray(new MetaQTL4MetaTrait[0]);
			HashMap<MetaQTL4MetaTrait, Integer> cisProbeMap = new HashMap<MetaQTL4MetaTrait, Integer>();
			MetaQTL4MetaTrait[] cisProbeArray = null;
			if (snpprobeCombos != null) {
				cisProbeArray = snpprobeCombos[snp];
			} else {
				Set<MetaQTL4MetaTrait> cisProbesForSNP = probeAnnotation.getMetatraits().getTraitInWindow(snpChr[snp], snpPositions[snp], settings.getCisdistance());
				cisProbeArray = cisProbesForSNP.toArray(new MetaQTL4MetaTrait[0]);
			}
			
			if (cisProbeArray == null || cisProbeArray.length == 0) {
				// nothing to do. skip variant //
			} else {
				int ctr = 0;
				for (MetaQTL4MetaTrait cisProbe : cisProbeArray) {
					cisProbeMap.put(cisProbe, ctr);
					ctr++;
				}
				
				
				float[][] finalZScores = new float[cisProbeMap.size()][datasets.length];
				
				// initialize with NaN
				for (int q = 0; q < probeIndex.length; q++) {
					for (int r = 0; r < datasets.length; r++) {
						finalZScores[q][r] = Float.NaN;
					}
				}
				
				// get list of probes to test for each dataset
				for (int d = 0; d < datasets.length; d++) {
					if (flipZScores[d] == null) {
						// the allele could not be flipped. set the Z to NaN
						for (float[] zScore : finalZScores) {
							zScore[d] = Float.NaN;
						}
					} else {
						//initialize z-score
						for (int p = 0; p < cisProbeMap.size(); p++) {
							finalZScores[p][d] = Float.NaN; // this is not very nice, but does prevent the metaZ method from going nuts
						}
						
						// load the z-scores for the dataset
						int datasetSNPId = snpIndex[snp][d];
						
						if (datasetSNPId != -9) { // -9 means: snp not available
							
							// TODO: for faster disk access, we should wrap this into a buffer of some sort
							float[] datasetZScores = datasets[d].getZScores(datasetSNPId);
							
							if (datasets[d].getIsCisDataset()) {
								// this requires us to retrieve the z-scores differently:
								// a cis dataset only stores z-scores for the tested probes/traits
								// position 0 in the datasetZScores array may therefore point to meta-trait 1000 in our annotation
								// we need to figure out which probes match up.. their orders might be different
								// and the number of probes tested in each dataset might differ as well
								
								// get the probes tested against the SNP
								MetaQTL4MetaTrait[] datasetCisProbes = datasets[d].getCisProbes(datasetSNPId);
								
								for (int i = 0; i < datasetCisProbes.length; i++) {
									MetaQTL4MetaTrait p = datasetCisProbes[i];
									if (p != null) {
										Integer index = cisProbeMap.get(p);
										if (index != null) {
											float datasetZ = datasetZScores[i];
											finalZScores[index][d] = datasetZ;
											if (flipZScores[d]) {
												finalZScores[index][d] *= -1;
											}
										}
									}
								}
							} else { // this is not a cis dataset
								// use the full probe index
								for (int probe = 0; probe < cisProbeArray.length; probe++) {
									MetaQTL4MetaTrait cisProbe = cisProbeArray[probe];
									Integer metaProbeIndex = traitMap.get(cisProbe);
									Integer datasetProbeId = probeIndex[metaProbeIndex][d];
									if (datasetProbeId != null) {
										finalZScores[probe][d] = datasetZScores[datasetProbeId];
										if (flipZScores[d]) {
											finalZScores[probe][d] *= -1;
										}
									}
								}
							}
							
							
						}
					}
				}
				
				// meta-analyze!
				for (int probe = 0; probe < finalZScores.length; probe++) {
					MetaQTL4MetaTrait t = cisProbeArray[probe];
					
					double metaZ = ZScores.getWeightedZ(finalZScores[probe], sampleSizes);
					int metaZNrSamples = 0;
					for (int s = 0; s < sampleSizes.length; s++) {
						if (!Float.isNaN(finalZScores[probe][s])) {
							metaZNrSamples += sampleSizes[s];
						}
					}
					double p = Descriptives.convertZscoreToPvalue(metaZ);
					
					if (settings.isMakezscoretable()) {
						// get the correct index for trait t
						int metaid = t.getCurrentMetaId();
						zscoretableoutput[metaid] = metaZ;
						zscorenrsamplestableoutput[metaid] = metaZNrSamples;
					}
					
					if (!Double.isNaN(p) && !Double.isNaN(metaZ)) {
						// create output object
						QTL q = new QTL(p, t, snp, BaseAnnot.toByte(alleleAssessed), metaZ, BaseAnnot.toByteArray(alleles), finalZScores[probe], sampleSizes); // sort buffer if needed.
//                            System.out.println(q.getSNPId()+"\t"+q.getMetaTrait().getMetaTraitName()+"\t"+q.toString());
						qtlOutput.add(q);
					} else {
//                            if (!printed) {
//                                printed = true;
//                                System.out.println("SNP creates NaN pval: " + snpList[snp] + "\n");//z " + Strings.concat(zScores, Strings.semicolon) + "\tn " + Strings.concat(sampleSizes, Strings.semicolon));
//                            }
					}
				}
			}
		} else {
			// analysis is not cis, but may be cis/trans
			Set<MetaQTL4MetaTrait> cisProbes = null;
			if (!settings.getAnalysisType().equals(BinaryMetaAnalysisSettings.Analysis.CISTRANS)) {
				// do not test the cis probes if not cistrans
				cisProbes = probeAnnotation.getMetatraits().getTraitInWindow(snpChr[snp], snpPositions[snp], settings.getTransdistance());
			}
			
			
			// iterate over the probe index
			float[][] finalZScores = new float[probeIndex.length][datasets.length];
			
			// initialize with NaN
			for (int q = 0; q < probeIndex.length; q++) {
				for (int r = 0; r < datasets.length; r++) {
					finalZScores[q][r] = Float.NaN;
				}
			}
			
			for (int d = 0; d < datasets.length; d++) {
				if (datasets[d].getIsCisDataset()) {
					System.err.println("ERROR: cannot run trans analysis on a cis dataset: " + settings.getDatasetlocations().get(d));
					System.exit(-1);
				}
				
				if (flipZScores[d] == null) {
					for (float[] zScore : finalZScores) {
						zScore[d] = Float.NaN;
					}
				} else {
					int datasetSNPId = snpIndex[snp][d];
					float[] datasetZScores = datasets[d].getZScores(datasetSNPId);

// probeIndex[t.getMetaTraitId()][d] = p;
					for (int p = 0; p < traitList.length; p++) {
						MetaQTL4MetaTrait t = traitList[p];
						if (cisProbes != null && cisProbes.contains(t)) {
							finalZScores[p][d] = Float.NaN;
						} else {
							Integer datasetProbeId = probeIndex[p][d];
							if (datasetProbeId != null) {
								finalZScores[p][d] = datasetZScores[datasetProbeId];
								if (flipZScores[d]) {
									finalZScores[p][d] *= -1;
								}
							} else {
								finalZScores[p][d] = Float.NaN;
							}
						}
					}
				}
			}
			
			// meta-analyze!
			if (settings.isMakezscoretable()) {
				// initialize with NaN
				for (int i = 0; i < zscoretableoutput.length; i++) {
					zscoretableoutput[i] = Double.NaN;
					zscorenrsamplestableoutput[i] = 0;
				}
			}
			for (int probe = 0; probe < traitList.length; probe++) {
				MetaQTL4MetaTrait t = traitList[probe];
				double metaAnalysisZ = ZScores.getWeightedZ(finalZScores[probe], sampleSizes);
				int metaAnalysisZNrSamples = 0;
				for (int s = 0; s < sampleSizes.length; s++) {
					if (!Float.isNaN(finalZScores[probe][s])) {
						metaAnalysisZNrSamples += sampleSizes[s];
					}
				}
				double metaAnalysisP = Descriptives.convertZscoreToPvalue(metaAnalysisZ);
				
				if (settings.isMakezscoretable()) {
					zscoretableoutput[probe] = metaAnalysisZ;
					zscorenrsamplestableoutput[probe] = metaAnalysisZNrSamples;
				}
				
				// create output object
				if (!Double.isNaN(metaAnalysisP) && !Double.isNaN(metaAnalysisZ)) {
					QTL q = new QTL(metaAnalysisP, t, snp, BaseAnnot.toByte(alleleAssessed), metaAnalysisZ, BaseAnnot.toByteArray(alleles), finalZScores[probe], sampleSizes); // sort buffer if needed.
//                            System.out.println(q.getSNPId()+"\t"+q.getMetaTrait().getMetaTraitName()+"\t"+q.toString());
					qtlOutput.add(q);
				}
			}
		}

//		if (qtlOutput.size() > 1) {
//			Collections.sort(qtlOutput);
//		}
		
		// write z-score output
		if (settings.isMakezscoretable()) {
			String snpName = snpList[snp];
			// get alleles
			String zscoreTableTf = snpName + "\t" + alleles + "\t" + alleleAssessed + "\t" + Strings.concat(zscoretableoutput, Strings.tab);
			String zscoreTableTfNrSamples = snpName + "\t" + alleles + "\t" + alleleAssessed + "\t" + Strings.concat(zscorenrsamplestableoutput, Strings.tab);
			output = new Triple<>(output.getLeft(), zscoreTableTf, zscoreTableTfNrSamples);
		}
		
		
		return output;
		
	}
}
