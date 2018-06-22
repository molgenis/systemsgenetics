package nl.umcg.westrah.binarymetaanalyzer.westrah.binarymetaanalyzer.posthoc;

import nl.umcg.westrah.binarymetaanalyzer.BinaryMetaAnalysis;
import nl.umcg.westrah.binarymetaanalyzer.BinaryMetaAnalysisSettings;
import nl.umcg.westrah.binarymetaanalyzer.MetaQTL4MetaTrait;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.trityper.util.BaseAnnot;

import java.io.IOException;
import java.util.HashMap;
import java.util.Set;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.atomic.AtomicInteger;

public class ConvertSummaryStats extends BinaryMetaAnalysis {
	
	
	public ConvertSummaryStats(String settingsFile, String textToReplace, String replaceTextWith, boolean usetmp) {
		super(settingsFile, textToReplace, replaceTextWith, usetmp);
		
	}
	
	public void run() throws IOException {
		super.initialize();
		
		String outdir = settings.getOutput();
		if (usetmp) {
			outdir = tempDir;
		}
		
		System.out.println("Placing output here: " + outdir);
		outdir = Gpio.formatAsDirectory(outdir);
		Gpio.createDir(outdir);
		
		
		System.out.println("Permutations: " + settings.getStartPermutations() + " until " + settings.getNrPermutations());
		
		int availableProcessors = Runtime.getRuntime().availableProcessors();
		int cores = settings.getNrThreads();
		if (cores < 1) {
			cores = 1;
		} else if (cores > availableProcessors) {
			cores = availableProcessors;
		}
		
		System.out.println("Will try to make use of " + cores + " CPU cores");
		System.out.println();
		
		
		// get data for region
		for (int permutation = settings.getStartPermutations(); permutation <= settings.getNrPermutations(); permutation++) {
			// load probe annotation and index
			// this particular probe annotation can take multiple probes for a single location into account.
			System.out.println("Permutation: " + permutation);
			
			// don't think we need this
			// Descriptives.initializeZScoreToPValue();
			
			// re-intialize for each permutation, just to be sure
			if (permutation > settings.getStartPermutations()) {
				initialize();
				System.out.println("Loading probe annotation from: " + settings.getProbetranslationfile());
				if (traitList.length == 0) {
					System.err.println("Error: no annotation loaded.");
					System.exit(-1);
				}
			}
			//			clearResultsBuffer();
			
			
			initdataset(permutation, outdir, false);
			
			System.out.println("Writing output here: " + settings.getOutput() + "Permutation-" + permutation + "-ZScores.dat");
			SummaryStatLDFile lfo = new SummaryStatLDFile(settings.getOutput() + "Permutation-" + permutation + "-ZScores.dat", SummaryStatLDFile.W, 1048576, false);
			lfo.writeHeader(datasets);
			
			boolean debug = false;
			ExecutorService executor = Executors.newFixedThreadPool(cores);
			System.out.println("Booting threadpool with " + cores + " threads");
			AtomicInteger ctr = new AtomicInteger();
			for (int snp = 0; snp < snpList.length; snp++) {
				ConversionTask t = new ConversionTask(snp, lfo, ctr);
				executor.submit(t);
				
			}
			
			int returned = ctr.get();
			ProgressBar pb = new ProgressBar(snpList.length, "Converting SNPs");
			while (returned < snpList.length) {
				pb.set(returned);
				try {
					Thread.sleep(1000);
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
				returned = ctr.get();
			}
			pb.close();
			System.out.println("Writing hashes");
			lfo.writeHashes(settings.getOutput() + "Permutation-" + permutation);
			lfo.close();
			
			System.out.println("Done");
			
			System.out.println();
//			}
		
		} // end iterate permutations
		
	}
	
	public class ConversionTask implements Runnable {
		private final AtomicInteger jobctr;
		boolean debug = false;
		int snp;
		SummaryStatLDFile lfo;
		
		public ConversionTask(int snp, SummaryStatLDFile lfo, AtomicInteger ctr) {
			this.lfo = lfo;
			this.snp = snp;
			this.jobctr = ctr;
		}
		
		@Override
		public void run() {
			try {
				//				int[] sampleSizes = new int[datasets.length];
				
				Boolean[] flipZScores = new Boolean[datasets.length];
				String alleles = null;
				String alleleAssessed = null;
				
				// determine whether to flip the alleles for a certain dataset
				for (int d = 0; d < datasets.length; d++) {
					int datasetSNPId = snpIndex[snp][d];
					if (datasetSNPId != -9) {
//						sampleSizes[d] = datasets[d].getSampleSize(datasetSNPId);
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
				if (settings.getAnalysisType().equals(BinaryMetaAnalysisSettings.Analysis.CIS)) {
					// do cis stuff, or stuff to specific sets of probes...
					if (debug) {
						System.out.println("Entering snpprobe/cis mode");
					}
					// get all the possible traits near the SNP
//                    Set<MetaQTL4MetaTrait> cisProbesForSNP = probeAnnotation.getMetatraits().getTraitInWindow(snpChr[snp], snpPositions[snp], settings.getCisdistance());
//                    MetaQTL4MetaTrait[] cisProbeArray = cisProbesForSNP.toArray(new MetaQTL4MetaTrait[0]);
					HashMap<MetaQTL4MetaTrait, Integer> cisProbeMap = new HashMap<MetaQTL4MetaTrait, Integer>();
					MetaQTL4MetaTrait[] cisProbeArray = null;
					
					Set<MetaQTL4MetaTrait> cisProbesForSNP = probeAnnotation.getMetatraits().getTraitInWindow(snpChr[snp], snpPositions[snp], settings.getCisdistance());
					cisProbeArray = cisProbesForSNP.toArray(new MetaQTL4MetaTrait[0]);
					
					
					if (cisProbeArray == null || cisProbeArray.length == 0) {
						// nothing to do. skip variant //
						if (debug) {
							System.out.println(snpList[snp] + "\thas no probes");
						}
						
						
					} else {
						if (debug) {
							System.out.println(snpList[snp] + "\thas " + cisProbeArray.length + " probes");
						}
						int ctr = 0;
						
						for (MetaQTL4MetaTrait cisProbe : cisProbeArray) {
							cisProbeMap.put(cisProbe, ctr);
							ctr++;
						}
						
						float[][] finalZScores = new float[cisProbeMap.size()][datasets.length];
						
						// writeHeader with NaN
						for (int q = 0; q < cisProbeMap.size(); q++) {
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
								//writeHeader z-score
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
												if (debug) {
													System.out.println("Dataset " + d + "\tProbe map: " + cisProbe.getMetaTraitName().toString() + "\tProbeId: " + datasetProbeId + "\tZBeforeFlip: " + datasetZScores[datasetProbeId]);
												}
												if (flipZScores[d]) {
													finalZScores[probe][d] *= -1;
												}
												if (debug) {
													System.out.println("Dataset " + d + "\tProbe map: " + cisProbe.getMetaTraitName().toString() + "\tProbeId: " + datasetProbeId + "\tZAfterFlip: " + finalZScores[probe][d]);
												}
											}
										}
									}
									
									// now we've loaded the z-scores, we should be able to save them
									String snpname = snpList[snp];
									for (int p = 0; p < finalZScores.length; p++) {
										float[] datasetZ = finalZScores[p];
										
										MetaQTL4MetaTrait cisProbe = cisProbeArray[p];
										// write to disk
										String genename = cisProbe.getMetaTraitName();
										lfo.writesync(snpname, genename, datasetZ);
									}
								}
							}
						}
						
					}
				}
				
			} catch (IOException e) {
				e.printStackTrace();
			}
			jobctr.getAndIncrement();
		}
	}
}
