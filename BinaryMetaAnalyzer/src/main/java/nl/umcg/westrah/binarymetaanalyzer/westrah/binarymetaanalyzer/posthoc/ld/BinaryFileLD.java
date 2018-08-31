package nl.umcg.westrah.binarymetaanalyzer.westrah.binarymetaanalyzer.posthoc.ld;

import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.util.BaseAnnot;
import umcg.genetica.text.Strings;
import umcg.genetica.util.Primitives;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;


public class BinaryFileLD {
	
	
	public static void main(String[] args) {
		String indir = "D:\\ld\\tmp\\";
		String geneloc = "D:\\ld\\tmp\\sortedGeneSNPCombos.txt.gz-genes.txt";
		int nrthreads = 1;
		boolean writezmat = false;
		int minnrsamples = 0;
		int minnrvals = 0;
		String outdir = "D:\\ld\\tmp\\out\\";
		String mode = "pearson";
		BinaryFileLD l = new BinaryFileLD();
		try {
			l.run(indir, geneloc, nrthreads, outdir, writezmat, minnrsamples, minnrvals, mode);
		} catch (IOException e) {
			e.printStackTrace();
		} catch (IllegalAccessException e) {
			e.printStackTrace();
		}
		
		
	}
	
	public enum MODE {
		PEARSON,
		WEIGHTEDPEARSON,
		INNERPRODUCT,
		WEIGHTEDINNERPRODUCT,
		NOLD
	}
	
	public MODE method;
	
	public void run(String indir, String geneloc, int nrthreads, String outdir, boolean writezmat, int minnrsamples, int minnrvalues, String mode) throws IOException, IllegalAccessException {
		System.out.println("In: " + indir);
		System.out.println("Genes: " + geneloc);
		System.out.println("Threads: " + nrthreads);
		System.out.println("Out: " + outdir);
		System.out.println("Zmat: " + writezmat);
		System.out.println("MinSamples: " + minnrsamples);
		System.out.println("MinValues: " + minnrvalues);
		
		Gpio.createDir(outdir);
		
		if (mode == null) {
			method = MODE.PEARSON;
		} else {
			if (mode.equals("pearson")) {
				method = MODE.PEARSON;
			} else if (mode.equals("weightedpearson")) {
				method = MODE.WEIGHTEDPEARSON;
			} else if (mode.equals("innerproduct")) {
				method = MODE.INNERPRODUCT;
			} else if (mode.equals("weightedinnerproduct")) {
				method = MODE.WEIGHTEDINNERPRODUCT;
			} else if (mode.equals("nold")) {
				method = MODE.NOLD;
			} else {
				method = MODE.PEARSON;
			}
		}
		
		boolean debug = false;
		
		File[] datasetlocs = getLocs(indir);
		int nrdatasets = datasetlocs.length;
		
		// read genes to test
		TextFile tf = new TextFile(geneloc, TextFile.R);
		ArrayList<String> genes = tf.readAsArrayList();
		HashSet<String> geneset = new HashSet<>();
		geneset.addAll(genes);
		tf.close();
		
		
		List<File> datasetlist = Arrays.asList(datasetlocs);
		
		
		List<SortedBinaryZScoreFile> datasets = datasetlist.parallelStream()
				.map(v -> loadDataset(v))
				.filter(Objects::nonNull)
				.collect(Collectors.toList());
		SortedBinaryZScoreFile[] dataset = datasets.toArray(new SortedBinaryZScoreFile[0]);
		
		// iterate genes
		ExecutorService executor = Executors.newFixedThreadPool(nrthreads);
		AtomicInteger nrRunning = new AtomicInteger(0);
		AtomicInteger nrCompleted = new AtomicInteger(0);
		int nrsub = 0;
		SortedBinaryZDataBlock[] currentblock = new SortedBinaryZDataBlock[nrdatasets];
		int nrdatasetswithgene = 0;
		for (int g = 0; g < genes.size(); g++) {
			String querygene = genes.get(g);
//            System.out.println("Gene: " + g + " / " + genes.nrElemsPerArray() + ": " + querygene);
			// get gene data for each dataset
			LinkedHashMap<String, Integer> snpmap = new LinkedHashMap<>();
			ArrayList<ArrayList<SortedBinaryZDataBlock>> geneblocks = new ArrayList<>();
			int ctr = 0;
			for (int d = 0; d < dataset.length; d++) {
				if (dataset[d].hasGene(querygene)) {
					nrdatasetswithgene++;
					// get last read block
					SortedBinaryZDataBlock b = currentblock[d];
					if (b == null) {
						// first block seen .. initialize
						b = dataset[d].readNextBlock();
					}
					
					// dataset has gene, but apparently it's not on the current line. continue reading until we hit it, or we're at end of file
					// (which should not be possible)
					if (b != null && !b.gene.equals(querygene)) {
						
						// fast forward
						dataset[d].skipTo(querygene);
						
						
						b = dataset[d].readNextBlock();
						
						if (!b.gene.equals(querygene)) {
							System.out.println("Dataset is supposed to have " + querygene + " but I skipped to the wrong position?");
							System.exit(-1);
						}
						
					}
					
					ArrayList<SortedBinaryZDataBlock> blocks = new ArrayList<>();
					
					while (b != null && b.gene.equals(querygene)) {
						// add
						blocks.add(b);
						if (!snpmap.containsKey(b.snp)) {
							snpmap.put(b.snp, ctr);
							ctr++;
						}
						
						b = dataset[d].readNextBlock();
						if (debug) {
							if (b != null && currentblock[d] != null) {
								System.out.println("Dataset: " + d + " current gene:" + currentblock[d].gene + ". Next gene: " + b.gene);
							} else if (b != null) {
								System.out.println("Dataset: " + d + " Next gene: " + b.gene);
							} else if (currentblock[d] != null) {
								System.out.println("Dataset: " + d + " current gene:" + currentblock[d].gene);
								
							}
						}
						currentblock[d] = b;
						
						
					}
					if (debug) {
						System.out.println("Done reading ds: " + d);
					}
					geneblocks.add(blocks);
				} else {
					geneblocks.add(new ArrayList<>());
				}
			}
			
			if (nrdatasetswithgene > 0 && snpmap.size() > 1) {
				SortedZCorrelationTask t = new SortedZCorrelationTask(geneblocks,
						querygene,
						snpmap,
						dataset,
						outdir,
						nrRunning,
						nrCompleted,
						writezmat,
						method,
						g,
						minnrsamples,
						minnrvalues
				);
				
				executor.submit(t);
				nrsub++;
				
				if (nrsub > 0 && nrsub % (nrthreads * 2) == 0) {
					// prevent clogging of memories.
					int nrrunningi = nrRunning.get();
					while (nrrunningi >= nrthreads) {
						try {
							System.out.println("Submitted: " + nrsub + "\tCompleted: " + nrCompleted.get() + "\tRunning: " + nrrunningi);
							Thread.sleep(10000);
						} catch (InterruptedException e) {
							e.printStackTrace();
						}
						nrrunningi = nrRunning.get();
					}
				}
			}
		}
		
		int nrrunningi = nrRunning.get();
		while (nrrunningi > 0) {
			try {
				System.out.println("Submitted: " + nrsub + "\tCompleted: " + nrCompleted.get() + "\tRunning: " + nrrunningi);
				Thread.sleep(10000); // wait 10 secs before checking again
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
			nrrunningi = nrRunning.get();
		}
		
		executor.shutdown();
		
		System.out.println("Submitted: " + nrsub + "\tCompleted: " + nrCompleted.get() + "\tRunning: " + nrrunningi);
		System.out.println("Done");
	}
	
	private SortedBinaryZScoreFile loadDataset(File dataset) {
		try {
			return new SortedBinaryZScoreFile(dataset, SortedBinaryZScoreFile.R);
		} catch (IOException e) {
			e.printStackTrace();
		}
		return null;
	}
	
	private File[] getLocs(String indir) {
		
		File dir = new File(indir);
		File[] list = dir.listFiles();
		ArrayList<File> locs = new ArrayList<>();
		for (File f : list) {
			if (!f.isDirectory() && f.getName().endsWith("-data.dat")) {
				System.out.println("Found: " + f.getName());
				locs.add(f);
			}
		}
		System.out.println(locs.size() + " total datasets.");
		return locs.toArray(new File[0]);
	}
	
	class SortedZCorrelationTask implements Runnable {
		
		private final AtomicInteger ctr;
		private final AtomicInteger cmctr;
		private final boolean writezmat;
		private final int id;
		private final MODE method;
		ArrayList<ArrayList<SortedBinaryZDataBlock>> geneblocks;
		String querygene;
		LinkedHashMap<String, Integer> snpmap;
		SortedBinaryZScoreFile[] dataset;
		String outdir;
		int minNrValues = 0;
		int minNrSamples = 0;
		
		
		public SortedZCorrelationTask(ArrayList<ArrayList<SortedBinaryZDataBlock>> geneblocks, String querygene,
									  LinkedHashMap<String, Integer> snpmap, SortedBinaryZScoreFile[] dataset,
									  String outdir, AtomicInteger ctr, AtomicInteger cmctr, boolean writezmat,
									  MODE method, int id, int minNrSamples, int minNrValues) {
			this.geneblocks = geneblocks;
			this.querygene = querygene;
			this.snpmap = snpmap;
			this.dataset = dataset;
			this.outdir = outdir;
			this.ctr = ctr;
			this.cmctr = cmctr;
			this.writezmat = writezmat;
			this.id = id;
			this.method = method;
			this.minNrSamples = minNrSamples;
			this.minNrValues = minNrValues;
		}
		
		@Override
		public void run() {
			try {
				
				ctr.getAndIncrement();
				
				// determine which snps were tested in which datasets
				System.out.println("Job" + id + "\tGene " + querygene + "\tBuilding block index: " + snpmap.size() + " x " + dataset.length);
				SortedBinaryZDataBlock[][] blockindex = new SortedBinaryZDataBlock[snpmap.size()][dataset.length];
				
				for (int d = 0; d < dataset.length; d++) {
					ArrayList<SortedBinaryZDataBlock> blocks = geneblocks.get(d);
					System.out.println("Job" + id + "\t" + blocks.size() + " blocks ");
					for (SortedBinaryZDataBlock b : blocks) {
						Integer id = snpmap.get(b.snp);
						if (id != null) {
							blockindex[id][d] = b;
						}
					}
				}
				
				System.out.println("Job" + id + "\tGene " + querygene + "\tFlipping alleles");
				
				// flip the z-scores, using a reference dataset (TODO: flip to 1kg reference allele?)
				SortedBinaryZDataBlock[] refBlocks = new SortedBinaryZDataBlock[snpmap.size()];
				for (int s = 0; s < snpmap.size(); s++) {
					String refAlleles = null;
					String refAllelesAssessed = null;
					for (int d = 0; d < dataset.length; d++) {
						SortedBinaryZDataBlock b = blockindex[s][d];
						if (b != null) {
							if (refAlleles == null) {
								refAlleles = b.allele;
								refAllelesAssessed = b.assessed;
								refBlocks[s] = b;
							} else {
								
								Boolean flip = BaseAnnot.flipalleles(refAlleles, refAllelesAssessed, b.allele, b.assessed);
								if (flip == null) {
									// something went wrong while flipping the snp
									blockindex[s][d] = null;
								} else if (flip) {
									for (int z = 0; z < b.z.length; z++) {
										b.z[z] *= -1;
									}
								}
							}
						}
					}
				}
				
				
				if (!method.equals(MODE.NOLD)) {
					PearsonsCorrelation c = new PearsonsCorrelation();
					TextFile output = new TextFile(outdir + querygene + ".txt.gz", TextFile.W, 32 * 1024);
					
					String header = "SNP1\tSNP2\tNrValues\tNrSamples\tR(" + method + ")\tRsq";
					output.writeln(header);
					for (int s1 = 0; s1 < snpmap.size(); s1++) {
						SortedBinaryZDataBlock[] snp1 = blockindex[s1];
						SortedBinaryZDataBlock refblock1 = refBlocks[s1];

//					int s2 = s1;
						for (int s2 = (s1 + 1); s2 < snpmap.size(); s2++) {
//					for (int s2 = 0; s2 < snpmap.nrElemsPerArray(); s2++) {
							SortedBinaryZDataBlock[] snp2 = blockindex[s2];
							
							ReturnObj zs = removeNullsAndConvertToDouble(snp1, snp2);
							if (zs != null && zs.x.length >= minNrValues && zs.nx[0] >= minNrSamples) {
								SortedBinaryZDataBlock refblock2 = refBlocks[s2];
								String outln = refblock1.snp
										+ "\t" + refblock2.snp;
								
								double r = 0;
								double rsq = 0;
								if (method.equals(MODE.WEIGHTEDPEARSON)) {
									r = weightedcorr(zs);
								} else if (method.equals(MODE.INNERPRODUCT)) {
									r = innerproduct(zs, false);
								} else if (method.equals(MODE.WEIGHTEDINNERPRODUCT)) {
									r = innerproduct(zs, true);
								} else {
									r = c.correlation(zs.x, zs.y);
								}
								
								rsq = r * r;
								
								outln += "\t"
										+ zs.x.length + "\t"
										+ zs.weight + "\t"
										+ r + "\t"
										+ rsq;
								
								output.writeln(outln);
							}
							
							
						}
					}
					output.close();
					
				}
				
				// write alleles
				System.out.println("Job" + id + "\tGene " + querygene + "\thas " + snpmap.size() + " snps");
				TextFile alleleout = new TextFile(outdir + querygene + "-referenceAlleles.txt.gz", TextFile.W, 32 * 1024);
				alleleout.writeln("SNP\tAlleles\tAssessed");
				
				for (int s1 = 0; s1 < snpmap.size(); s1++) {
					SortedBinaryZDataBlock[] snp1 = blockindex[s1];
					SortedBinaryZDataBlock refblock1 = refBlocks[s1];
					alleleout.writeln(refblock1.snp + "\t" + refblock1.allele + "\t" + refblock1.assessed);
				}
				alleleout.close();
				
				
				// write zscore matrices.
				if (writezmat) {
					System.out.println("Job" + id + "\tGene " + querygene + "\tInitializing z-mat: " + snpmap.size() + " x " + (dataset.length * 10));
					float[][] zout = new float[dataset.length * 10][snpmap.size()];
					float[][] nout = new float[dataset.length][snpmap.size()];
					
					for (int s1 = 0; s1 < snpmap.size(); s1++) {
						SortedBinaryZDataBlock[] snp1 = blockindex[s1];
						int j = 0;
						for (int d = 0; d < snp1.length; d++) {
							
							if (snp1[d] == null) {
								for (int q = 0; q < 10; q++) {
									zout[j][s1] = Float.NaN;
									
									j++;
								}
								
							} else {
								nout[d][s1] = snp1[d].n;
								float[] zs = snp1[d].z;
								for (int q = 0; q < 10; q++) {
									zout[j][s1] = zs[q];
									j++;
								}
								
							}
						}
					}
					
					System.out.println("Job" + id + "\tGene " + querygene + "\tWriting mat");
					TextFile outz = new TextFile(outdir + querygene + "-zmat.txt.gz", TextFile.W);
					TextFile outn = new TextFile(outdir + querygene + "-nmat.txt.gz", TextFile.W);
					String header = "-";
					
					for (int s1 = 0; s1 < snpmap.size(); s1++) {
						header += "\t" + refBlocks[s1].snp;
					}
					outz.writeln(header);
					outn.writeln(header);
//					for (int d = 0; d < dataset.length; d++) {
//						for (int p = 1; p < 11; p++) {
//							header += "\t" + dataset[d].getName() + "-perm-" + p;
//						}
//					}
					int j = 0;
					for (int d = 0; d < dataset.length; d++) {
						String lnn = dataset[d].getName() + "\t" + Strings.concat(nout[d], Strings.tab);
						for (int p = 0; p < 10; p++) {
							float[] zs = zout[j];
							String ln = dataset[d].getName() + "-perm-" + (p + 1) + "\t" + Strings.concat(zs, Strings.tab);
							outz.writeln(ln);
							j++;
						}
						outn.writeln(lnn);
					}
					outz.close();
					outn.close();
					System.out.println("Job" + id + "\tGene " + querygene + "\tDone writing mat");
				}
				
			} catch (IOException e) {
				e.printStackTrace();
			} catch (Exception e) {
				e.printStackTrace();
			}
			
			ctr.getAndDecrement();
			cmctr.getAndIncrement();
		}
		
		private double innerproduct(ReturnObj vals, boolean weightforsamplesize) {
			double productsum = 0;
			
			double samplesum = 0;
			for (int i = 0; i < vals.x.length; i++) {
				if (weightforsamplesize) {
					double weight = Math.sqrt((vals.nx[i] + vals.ny[i]) / 2);
					productsum += (vals.x[i] * vals.y[i]) * weight;
					samplesum += weight;
				} else {
					productsum += (vals.x[i] * vals.y[i]);
				}
			}
			
			double ip = 0;
			
			
			if (weightforsamplesize) {
//				return sum / Math.sqrt(samplesum);
				ip = productsum / samplesum;
			} else {
				ip = productsum / vals.x.length;
			}
			
			if (ip > 1) {
				ip = 1;
			} else if (ip < -1) {
				ip = -1;
			}
			
			return ip;
		}
		
	}
	
	
	public class ReturnObj {
		
		public double weight;
		double[] x;
		double[] y;
		double[] nx;
		double[] ny;
		
		double mx = 0;
		double my = 0;
		double sx = 0;
		double sy = 0;
		double mnx = 0;
		double mny = 0;
		
		public ReturnObj(ArrayList<Double> x, ArrayList<Double> y, ArrayList<Double> n1, ArrayList<Double> n2) {
			this.x = Primitives.toPrimitiveArr(x.toArray(new Double[0]));
			this.y = Primitives.toPrimitiveArr(y.toArray(new Double[0]));
			this.nx = Primitives.toPrimitiveArr(n1.toArray(new Double[0]));
			this.ny = Primitives.toPrimitiveArr(n2.toArray(new Double[0]));
			this.weight = average(nx, ny);

//			mx = JSci.maths.ArrayMath.mean(this.x);
//			my = JSci.maths.ArrayMath.mean(this.y);
//			mnx = JSci.maths.ArrayMath.mean(this.nx);
//			mny = JSci.maths.ArrayMath.mean(this.ny);
//			sx = JSci.maths.ArrayMath.standardDeviation(this.x);
//			sy = JSci.maths.ArrayMath.standardDeviation(this.y);
		}
		
		private double average(double[] nx, double[] ny) {
			double sum = 0;
			for (int i = 0; i < nx.length; i++) {
				sum += nx[i] + ny[i];
			}
			sum /= 2;
			return sum;
		}
		
		
	}
	
	private ReturnObj removeNullsAndConvertToDouble(SortedBinaryZDataBlock[] snp1, SortedBinaryZDataBlock[] snp2) {
		ArrayList<Double> x = new ArrayList<>();
		ArrayList<Double> y = new ArrayList<>();
		ArrayList<Double> n1 = new ArrayList<>();
		ArrayList<Double> n2 = new ArrayList<>();
		
		for (int i = 0; i < snp1.length; i++) {
			SortedBinaryZDataBlock b1 = snp1[i];
			SortedBinaryZDataBlock b2 = snp2[i];
			if (b1 != null && b2 != null) {
				if (b1.n > 0 && b2.n > 0) {
					for (int q = 0; q < b1.z.length; q++) {
						x.add((double) b1.z[q]);
						y.add((double) b2.z[q]);
						n1.add((double) b1.n);
						n2.add((double) b2.n);
					}
				}
			}
		}
		
		if (x.size() > 0) {
			return new ReturnObj(x, y, n1, n2);
		} else {
			return null;
		}
		
		
	}
	
	
	public double weightedcorr(ReturnObj vals) {
		
		
		double[] w = new double[vals.x.length];
		
		double xavg = 0;
		double yavg = 0;
		
		for (int i = 0; i < vals.x.length; i++) { // for each dataset
			double wi = Math.sqrt((vals.nx[i] + vals.ny[i]) / 2);
			xavg += wi * vals.x[i];
			yavg += wi * vals.y[i];
			w[i] = wi;
		}
		
		xavg /= vals.x.length;
		yavg /= vals.x.length;
		
		double lxy = 0;
		double lxx = 0;
		double lyy = 0;
		
		for (int i = 0; i < vals.x.length; i++) {
			double wi = w[i];
			double xi = vals.x[i];
			double yi = vals.y[i];
			double xminavg = (xi - xavg);
			double yminavg = (yi - yavg);
			lxy += wi * xminavg * yminavg;
			lxx += wi * xminavg * xminavg;
			lyy += wi * yminavg * yminavg;
		}
		
		double r = lxy / Math.sqrt(lxx * lyy);
		return r;
	}
	
	public double weightedcorr2(ReturnObj vals) {
		
		
		double[] w = new double[vals.x.length];
		
		double xavg = 0;
		double yavg = 0;
		double sumw = 0;
		
		for (int i = 0; i < vals.x.length; i++) { // for each dataset
			double wi = Math.sqrt((vals.nx[i] + vals.ny[i]) / 2);
			xavg += wi * vals.x[i];
			yavg += wi * vals.y[i];
			sumw += wi;
			w[i] = wi;
		}
		
		xavg /= sumw;
		yavg /= sumw;
		
		double lxy = 0;
		double lxx = 0;
		double lyy = 0;
		
		for (int i = 0; i < vals.x.length; i++) {
			double wi = w[i];
			double xi = vals.x[i];
			double yi = vals.y[i];
			double xminavg = (xi - xavg);
			double yminavg = (yi - yavg);
			lxy += wi * xminavg * yminavg;
			lxx += wi * xminavg * xminavg;
			lyy += wi * yminavg * yminavg;
		}
		lxy /= sumw;
		lxx /= sumw;
		lyy /= sumw;
		
		double r = lxy / Math.sqrt(lxx * lyy);
		return r;
	}
}
