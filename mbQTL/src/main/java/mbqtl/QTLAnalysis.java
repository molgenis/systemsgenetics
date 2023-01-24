package mbqtl;

import mbqtl.datastructures.Dataset;
import mbqtl.datastructures.GeneAnnotation;
import mbqtl.datastructures.GeneExpressionData;
import mbqtl.datastructures.SNPAnnotation;
import umcg.genetica.containers.Triple;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.stats.Correlation;
import umcg.genetica.math.stats.HWE;

import java.io.IOException;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.*;

public class QTLAnalysis {

	protected final String origvcfFile;
	protected String vcfFile;

	protected final GeneExpressionData expressionData;
	protected Dataset[] datasets;
	protected final GeneAnnotation geneAnnotation;
	protected boolean[] genotypeSamplesToInclude;
	protected final String outputPrefix;
	protected final int chromosome;
	protected Dataset combinedDataset;
	private final String linkfile;
	protected HashMap<String, HashSet<String>> snpGeneLimitSet;
	protected Set<String> snpLimitSet;
	protected HashMap<String, String> RNAToDNA = new HashMap<>();
	protected HashMap<String, String> DNAtoRNA = new HashMap<>();
	protected HashMap<String, String> RNAToDataset = new HashMap<>();
	protected HashMap<String, String> RNAsPerDataset;

	protected HashMap<String, HashSet<String>> geneGroups = null;

	protected double mafthreshold = 0.01;
	protected double callratethreshold = 0.95;
	protected double hwepthreshold = 0.0001;
	protected int minNumberOfDatasets = 2;

	public void setMinObservations(int minObservations) {
		this.minObservations = minObservations;
	}

	protected int minObservations = 10;

	protected SNPAnnotation snpAnnotation = null;

	public void loadSNPAnnotation(String snpannotation) throws IOException {
		snpAnnotation = new SNPAnnotation(snpannotation);
	}

	protected void updateDatasets(String vcfFile) throws IOException {
		// TODO: check if sample size is equal to previously loaded VCF, other wise crash

		ArrayList<String> genotypeSamples = getGenotypeSamples(vcfFile);
		HashMap<String, Integer> genotypeSampleHash = Util.hash(genotypeSamples);

		// determine genotype samples with links to an RNA-seq sample
		// load sample links
		loadSampleLinks(linkfile, genotypeSampleHash.keySet());

		// determine RNA samples to load
		HashSet<String> RNASamplesMatchedToDNA = new HashSet<String>();
		for (int i = 0; i < genotypeSamples.size(); i++) {
			String RNASample = DNAtoRNA.get(genotypeSamples.get(i));
			if (RNASample != null) {
				RNASamplesMatchedToDNA.add(RNASample);
			}
		}

		// determine which genotype samples to load using the available RNA samples
		HashSet<String> availableExpressionSamples = new HashSet<>();
		for (String s : expressionData.samples) {
			availableExpressionSamples.add(s);
		}

		// make a list of columns to load from the VCF file
		genotypeSamplesToInclude = new boolean[genotypeSamples.size()];
		ArrayList<String> loadedGenotypeSamples = new ArrayList<>();
		for (int i = 0; i < genotypeSamples.size(); i++) {
			String DNA = genotypeSamples.get(i);
			String RNA = DNAtoRNA.get(DNA);
			if (availableExpressionSamples.contains(RNA)) {
				genotypeSamplesToInclude[i] = true;
				loadedGenotypeSamples.add(DNA);
			}
		}

		// replace genotype sample list, and rehash
		genotypeSamples = loadedGenotypeSamples;
		genotypeSampleHash = Util.hash(genotypeSamples);
		System.out.println(genotypeSamples.size() + " DNA samples match loaded RNA samples");

		// divide samples up into datasets
		datasets = createDatasets(genotypeSampleHash, expressionData.sampleMap);

		combinedDataset = createCombinedDataset(genotypeSampleHash, expressionData.sampleMap);
		System.out.println("Combined dataset has: " + combinedDataset.expressionIds.length + " samples.");

	}

	public QTLAnalysis(String vcfFile,
					   int chromosome,
					   String linkfile,
					   String snpLimitFile,
					   String geneLimitFile,
					   String snpGeneLimitFile,
					   String geneExpressionDataFile,
					   String geneAnnotationFile,
					   String outputPrefix) throws IOException {
		System.out.println("-------------------------");
		System.out.println("Initializing QTL analysis");
		System.out.println("-------------------------");

		DateTimeFormatter dtf = DateTimeFormatter.ofPattern("yyyy/MM/dd HH:mm:ss");
		LocalDateTime now = LocalDateTime.now();
		System.out.println("Date:\t" + dtf.format(now));
		System.out.println("VCF:\t" + vcfFile);
		System.out.println("Chrom:\t" + chromosome);
		System.out.println("Expression: " + geneExpressionDataFile);
		System.out.println("Linkfile: " + linkfile);
		System.out.println("SNP limit: " + snpLimitFile);
		System.out.println("Gene limit: " + geneLimitFile);
		System.out.println("SNP/Gene limit: " + snpGeneLimitFile);
		System.out.println("Annotation: " + geneAnnotationFile);
		System.out.println("Outdir: " + outputPrefix);
		System.out.println();

		this.linkfile = linkfile;
		this.vcfFile = vcfFile;
		this.origvcfFile = vcfFile;
		this.chromosome = chromosome;

		this.outputPrefix = outputPrefix;
		// load genotyped samples

		if (!Gpio.exists(vcfFile)) {
			// try a replacement
			System.out.println("Could not find: " + vcfFile);
			System.out.println("Trying to replace CHR with chromosome number.");
			if (chromosome > 0) {
				vcfFile = vcfFile.replaceAll("CHR", "" + chromosome);
				if (!Gpio.exists(vcfFile)) {
					System.out.println("Could not find VCF file for chromosome: " + chromosome);
					System.exit(0);
				}
				this.vcfFile = vcfFile;
			} else {
				System.out.println("Could not find " + vcfFile);
				for (int i = 1; i < 23; i++) {
					vcfFile = vcfFile.replaceAll("CHR", "" + i);
					System.out.println("Trying: " + vcfFile);
					if (Gpio.exists(vcfFile)) {
						break;
					}
				}
			}

			if (!Gpio.exists(vcfFile)) {
				System.out.println("Could not find VCF file or chromosome replacement.");
				System.exit(-1);
			} else {
				this.vcfFile = vcfFile;
			}
		}

		ArrayList<String> genotypeSamples = getGenotypeSamples(vcfFile);
		HashMap<String, Integer> genotypeSampleHash = Util.hash(genotypeSamples);

		// determine genotype samples with links to an RNA-seq sample
		// load sample links
		loadSampleLinks(linkfile, genotypeSampleHash.keySet());

		// determine RNA samples to load
		HashSet<String> RNASamplesMatchedToDNA = new HashSet<String>();
		for (int i = 0; i < genotypeSamples.size(); i++) {
			String RNASample = DNAtoRNA.get(genotypeSamples.get(i));
			if (RNASample != null) {
				RNASamplesMatchedToDNA.add(RNASample);
			}
		}

		// load set of genes to limit to
		Set<String> geneLimitSet = null;
		if (geneLimitFile != null) {
			System.out.println("Limiting to genes found in : " + geneLimitFile);
			TextFile tf2 = new TextFile(geneLimitFile, TextFile.R);
			geneLimitSet = tf2.readAsSet(0, TextFile.tab);
			tf2.close();
			System.out.println("File contains " + geneLimitSet.size() + " genes.");
		}

		if (snpLimitFile != null) {
			System.out.println("Limiting to snps found in : " + snpLimitFile);
			TextFile tf2 = new TextFile(snpLimitFile, TextFile.R);
			snpLimitSet = tf2.readAsSet(0, TextFile.tab);
			tf2.close();
			System.out.println("File contains " + snpLimitSet.size() + " SNPs.");
		}

		if (snpGeneLimitFile != null) {
			System.out.println("Limiting to snp/gene pairs found in : " + snpGeneLimitFile);
			snpGeneLimitSet = new HashMap<String, HashSet<String>>();
			TextFile tf2 = new TextFile(snpGeneLimitFile, TextFile.R);
			String[] elems = tf2.readLineElems(TextFile.tab);
			Set<String> snpset = new HashSet<>();
			while (elems != null) {
				String snp = elems[0];
				String gene = elems[1];
				HashSet<String> set = snpGeneLimitSet.get(gene);
				if (set == null) {
					set = new HashSet<>();
				}
				set.add(snp);
				snpset.add(snp);
				snpGeneLimitSet.put(gene, set);
				elems = tf2.readLineElems(TextFile.tab);
			}
			snpLimitSet = snpset;
			tf2.close();
			System.out.println("File contains " + snpLimitSet.size() + " SNPs.");
			System.out.println("File contains " + snpGeneLimitSet.size() + " genes.");
			if (geneLimitSet == null) {
				geneLimitSet = snpGeneLimitSet.keySet();
				System.out.println("Limiting to " + geneLimitSet.size() + " genes");
			}
//            else {
//                HashSet<String> tmpgenelimit = new HashSet<>();
//                for (String gene : geneLimitSet) {
//                    if (snpGeneLimitSet.containsKey(gene)) {
//                        tmpgenelimit.add(gene);
//                    }
//                }
////				for (String s : snpGeneLimitSet.keySet()) {
////					if (!geneLimitSet.contains(s)) {
////						System.out.println("Could not find " + s);
////					}
////				}
//                geneLimitSet = tmpgenelimit;
//                System.out.println("File contains " + geneLimitSet.size() + " genes after matching with " + geneLimitFile);
//            }

		}

		// load gene annotation
		geneAnnotation = new GeneAnnotation(geneAnnotationFile, chromosome, geneLimitSet);

		// load expression data
		expressionData = new GeneExpressionData(geneExpressionDataFile, geneAnnotation.getAllGenes(), RNASamplesMatchedToDNA);

		// determine which genotype samples to load using the available RNA samples
		HashSet<String> availableExpressionSamples = new HashSet<>();
		for (String s : expressionData.samples) {
			availableExpressionSamples.add(s);
		}

		// make a list of columns to load from the VCF file
		genotypeSamplesToInclude = new boolean[genotypeSamples.size()];
		ArrayList<String> loadedGenotypeSamples = new ArrayList<>();
		for (int i = 0; i < genotypeSamples.size(); i++) {
			String DNA = genotypeSamples.get(i);
			String RNA = DNAtoRNA.get(DNA);
			if (availableExpressionSamples.contains(RNA)) {
				genotypeSamplesToInclude[i] = true;
				loadedGenotypeSamples.add(DNA);
			}
		}

		// replace genotype sample list, and rehash
		genotypeSamples = loadedGenotypeSamples;
		genotypeSampleHash = Util.hash(genotypeSamples);
		System.out.println(genotypeSamples.size() + " DNA samples match loaded RNA samples");

		// divide samples up into datasets
		datasets = createDatasets(genotypeSampleHash, expressionData.sampleMap);

		combinedDataset = createCombinedDataset(genotypeSampleHash, expressionData.sampleMap);
		System.out.println("Combined dataset has: " + combinedDataset.expressionIds.length + " samples.");

		// initialize Z-score lookup table
		Correlation.correlationToZScore(combinedDataset.expressionIds.length);

		if (datasets.length < minNumberOfDatasets) {
			System.out.println("Warning: min number of datasets " + minNumberOfDatasets + " smaller than actual number of datasets " + datasets.length);
			System.out.println("Setting variable to match.");
			minNumberOfDatasets = datasets.length;
		}
		System.out.println("Initialization done.");
		System.out.println();
	}

	protected Dataset[] createDatasets(HashMap<String, Integer> DNASampleMap, HashMap<String, Integer> RNASampleMap) {
		System.out.println("Creating datasets.");
		HashMap<String, Dataset> datasetMap = new HashMap<>();
		ArrayList<Dataset> datasets = new ArrayList<>();
		for (String RNA : RNASampleMap.keySet()) {
			String DNA = RNAToDNA.get(RNA);

			if (DNA != null) {
				String datasetName = RNAToDataset.get(RNA);
				Dataset dataset = datasetMap.get(datasetName);
				if (dataset == null) {
					dataset = new Dataset();
					dataset.name = datasetName;
					datasets.add(dataset);
					datasetMap.put(datasetName, dataset);
				}
				Integer RNAId = RNASampleMap.get(RNA);
				Integer DNAId = DNASampleMap.get(DNA);
				if (RNAId != null && DNAId != null) {
					dataset.append(DNAId, RNAId);
				}
			}
		}

		int nrDatasetsWithData = 0;
		for (int d = 0; d < datasets.size(); d++) {
			if (!datasets.get(d).RNAIds.isEmpty()) {
				nrDatasetsWithData++;
			} else {
				System.out.println(datasets.get(d).name + " has no data!");
			}
		}

		Dataset[] output = new Dataset[nrDatasetsWithData];
		int dctr = 0;
		for (int d = 0; d < datasets.size(); d++) {
			if (!datasets.get(d).RNAIds.isEmpty()) {
				output[dctr] = datasets.get(d);
				dctr++;
			}
		}
		Arrays.sort(output);
		System.out.println();
		System.out.println("Available datasets: " + datasets.size());
		for (Dataset d : output) {
			System.out.println(d.name + "\t" + d.RNAIds.size() + " samples.");
			d.toStr();
		}
		return output;
	}

	protected void loadSampleLinks(String linkfile, Set<String> includedDNAs) throws IOException {
		System.out.println("Reading: " + linkfile);
		TextFile tf = new TextFile(linkfile, TextFile.R);
		String[] elems = tf.readLineElems(TextFile.tab);
		RNAsPerDataset = new HashMap<>();
		while (elems != null) {
			String gt = elems[0];
			if (includedDNAs.contains(gt)) {
				String rna = elems[1];
				String dataset = elems[2];
				RNAToDNA.put(rna, gt);
				DNAtoRNA.put(gt, rna);
				RNAToDataset.put(rna, dataset);
			}
			elems = tf.readLineElems(TextFile.tab);
		}

		tf.close();
		System.out.println(RNAToDataset.size() + " links between RNA and datasets.");
		System.out.println(RNAToDNA.size() + " links between RNA and DNA");
		System.out.println(DNAtoRNA.size() + " links between DNA and RNA");
	}

	protected ArrayList<String> getGenotypeSamples(String vcfFile) throws IOException {
		TextFile tf = new TextFile(vcfFile, TextFile.R);
		String ln = tf.readLine();
		ArrayList<String> samples = new ArrayList<>();
		while (ln != null) {
			if (ln.startsWith("#")) {
				if (ln.startsWith("#CHROM")) {
					String[] elems = ln.split("\t");
					for (int i = 9; i < elems.length; i++) {
						samples.add(elems[i]);
					}
				}
			} else {
				break;
			}
			ln = tf.readLine();
		}
		tf.close();
		System.out.println(samples.size() + " samples in " + vcfFile);
		return samples;
	}

	protected Dataset createCombinedDataset(HashMap<String, Integer> DNASampleMap, HashMap<String, Integer> RNASampleMap) {
		Dataset dataset = new Dataset();
		dataset.name = "Combined";
		for (String RNA : RNASampleMap.keySet()) {
			String DNA = RNAToDNA.get(RNA);

			if (DNA != null) {
				Integer RNAId = RNASampleMap.get(RNA);
				Integer DNAId = DNASampleMap.get(DNA);
				if (RNAId != null && DNAId != null) {
					dataset.append(DNAId, RNAId);
				}
			}
		}
		dataset.toStr();
		return dataset;
	}

	protected VariantQCObj checkVariant(double[] gt) {
		int obsAA = 0;
		int obsAB = 0;
		int obsBB = 0;

		double freqA = 0;
		double called = 0;
		for (int i = 0; i < gt.length; i++) {
			double gti = gt[i];
			if (gti != -1) {
				called++;
				if (gt[i] == 0) {
					freqA += 2;
					obsAA++;
				} else if (gt[i] == 1) {
					freqA += 1;
					obsAB++;
				} else {
					obsBB++;
				}
			}
		}
		double maf = 0;

		double hwep = HWE.calculateExactHWEPValue(obsAB, obsAA, obsBB);

		if (called == 0) {
			maf = 0;
		} else {
			maf = freqA / (called * 2);
			called /= gt.length;
			if (maf > 0.5) {
				maf = 1 - maf;
			}
		}

		VariantQCObj obj = new VariantQCObj();
		obj.maf = maf;
		obj.cr = called;
		obj.hwep = hwep;
		obj.passqc = ((obsAA > 0 && obsAB > 0) || (obsAB > 0 && obsBB > 0) || (obsAA > 0 && obsBB > 0))
				&& (maf >= mafthreshold && called >= callratethreshold && hwep >= hwepthreshold);
		return obj;
	}

	protected Triple<double[], double[], double[]> pruneMissingValues(double[] datasetGenotypeData,
																	  double[] datasetGenotypeDosages,
																	  double[] datasetExpressionData) {
		int nrmissing = 0;
		for (int i = 0; i < datasetGenotypeData.length; i++) {
			if (datasetGenotypeData[i] == -1 || Double.isNaN(datasetExpressionData[i])) {
				nrmissing++;
			}
		}
		double[] genotypes = new double[datasetGenotypeData.length - nrmissing];
		double[] dosages = new double[datasetGenotypeData.length - nrmissing];
		double[] expression = new double[datasetGenotypeData.length - nrmissing];
		int ctr = 0;
		for (int i = 0; i < datasetGenotypeData.length; i++) {
			if (datasetGenotypeData[i] != -1 && !Double.isNaN(datasetExpressionData[i])) {
				expression[ctr] = datasetExpressionData[i];
				genotypes[ctr] = datasetGenotypeData[i];
				dosages[ctr] = datasetGenotypeDosages[i];
//                dosages[ctr] = datasetGenotypeData[i];
				ctr++;
			}
		}
		return new Triple<>(genotypes, dosages, expression);
	}

	public void setMafthreshold(double mafthreshold) {
		this.mafthreshold = mafthreshold;
	}

	public void setCallratethreshold(double callratethreshold) {
		this.callratethreshold = callratethreshold;
	}

	public void setHwepthreshold(double hwepthreshold) {
		this.hwepthreshold = hwepthreshold;
	}

	public void setMinNumberOfDatasets(int minNumberOfDatasets) {
		this.minNumberOfDatasets = minNumberOfDatasets;
	}

	protected double[] getDosage(double[][] dosage) {
		double[] output = new double[dosage.length];
		for (int d = 0; d < output.length; d++) {
			output[d] = dosage[d][0];
		}
		return output;
	}

	protected double[] getGenotype(byte[] genotypesAsByteVector) {
		double[] genotypes = new double[genotypesAsByteVector.length];
		for (int d = 0; d < genotypes.length; d++) {
			genotypes[d] = genotypesAsByteVector[d];
		}
		return genotypes;
	}

// TODO: merge group annotation with gene annotation object.
	public void loadExpGroupAnnotation(String expgroups) throws IOException {
		System.out.println("Loading gene groupings from: " + expgroups);
		geneGroups = new HashMap<>();
		TextFile tf = new TextFile(expgroups, TextFile.R);
		String[] elems = tf.readLineElems(TextFile.tab);
		HashSet<String> geneSetTmp = new HashSet<>();

		while (elems != null) {
			if (elems.length > 1) {
				String geneID = elems[0];
				String groupID = elems[1];
				if (expressionData.geneMap.containsKey(geneID)) {
					HashSet<String> groupSet = geneGroups.get(groupID);
					if (groupSet == null) {
						groupSet = new HashSet<>();
					}
					groupSet.add(geneID);
					geneSetTmp.add(geneID);
					geneGroups.put(groupID, groupSet);
				}
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		System.out.println(geneGroups.size() + " groups loaded, encompassing " + geneSetTmp.size() + " genes that overlap with loaded expression data.");
		if (geneGroups.isEmpty()) {
			System.out.println("List of gene groups is empty. This can happen if none of the genes in groups is on the specified chromosome, or there is no ID overlap with the expression data. ");
			System.exit(-1);
		}
	}
}
