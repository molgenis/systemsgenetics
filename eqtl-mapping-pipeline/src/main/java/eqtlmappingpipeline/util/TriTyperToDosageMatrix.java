package eqtlmappingpipeline.util;

import eqtlmappingpipeline.metaqtl3.containers.Settings;
import org.apache.commons.configuration.ConfigurationException;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.*;
import umcg.genetica.io.trityper.util.BaseAnnot;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.*;
import java.util.stream.IntStream;

public class TriTyperToDosageMatrix {


	public void run(String settingsfile) throws IOException, ConfigurationException {
		System.out.println("TriTyper to dosagematrix converter.");
		Settings settings = new Settings();
		settings.load(settingsfile);

		// initialize datasets
		ArrayList<TriTyperGeneticalGenomicsDatasetSettings> datasetsettings = settings.datasetSettings;
		TriTyperGenotypeData[] m_gg = new TriTyperGenotypeData[datasetsettings.size()];
		SNPLoader[] loaders = new SNPLoader[datasetsettings.size()];
		System.out.println("Initializing " + datasetsettings.size() + " datasets in parallel.");
		IntStream.range(0, m_gg.length).parallel().forEach(i -> {
			m_gg[i] = new TriTyperGenotypeData();
			try {
				m_gg[i].load(datasetsettings.get(i).genotypeLocation, datasetsettings.get(i).snpmapFileLocation, datasetsettings.get(i).snpFileLocation);
				loaders[i] = m_gg[i].createSNPLoader(1);
			} catch (IOException e) {
				e.printStackTrace();
			}
		});


		System.out.println("Selecting samples per dataset");
		ArrayList<HashSet<String>> allowedSamplesPerDs = new ArrayList<HashSet<String>>();
		for (int d = 0; d < m_gg.length; d++) {
			TriTyperGeneticalGenomicsDatasetSettings dataset = datasetsettings.get(d);
			HashSet<String> allowedSamples = new HashSet<String>();
			if (dataset.genotypeToExpressionCoupling == null) {
				allowedSamples.addAll(Arrays.asList(m_gg[d].getIndividuals()));
			} else {
				TextFile tf = new TextFile(dataset.genotypeToExpressionCoupling, TextFile.R);
				String[] elems = tf.readLineElems(TextFile.tab);
				while (elems != null) {
					allowedSamples.add(elems[0]);
					elems = tf.readLineElems(TextFile.tab);
				}
				tf.close();
			}
			System.out.println(datasetsettings.get(d).name + " has " + allowedSamples.size() + " samples selected.");
			allowedSamplesPerDs.add(allowedSamples);
		}

		// create sample index
		System.out.println("Creating sample index..");
		HashMap<String, Integer> sampleIdMap = new HashMap<String, Integer>();
		int sampleCounter = 0;
		for (int i = 0; i < allowedSamplesPerDs.size(); i++) {
			HashSet<String> samplesInDs = allowedSamplesPerDs.get(i);
			for (String s : samplesInDs) {
				sampleIdMap.put(s, sampleCounter);
				sampleCounter++;
			}
		}
		System.out.println(sampleIdMap.size() + " unique sample IDs.");

		// write data to disk
		String[] header = new String[sampleIdMap.size()];
		for (Map.Entry<String, Integer> samplePair : sampleIdMap.entrySet()) {
			header[samplePair.getValue()] = samplePair.getKey();
		}
		TextFile tf = new TextFile(settings.outputReportsDir + "GenotypeData.txt.gz", TextFile.W);
		tf.writeln("SNP\tAlleles\tMinorAllele\t" + Strings.concat(header, Strings.tab));

		ArrayList<String> snpsToQuery = new ArrayList<>();
		snpsToQuery.addAll(settings.tsSNPsConfine);
		Collections.sort(snpsToQuery);
		double[] data = new double[sampleIdMap.size()];
		ProgressBar pb = new ProgressBar(snpsToQuery.size(), "Querying SNPs...");
		for (int snpctr = 0; snpctr < snpsToQuery.size(); snpctr++) {
			String snp = snpsToQuery.get(snpctr);
			SNP[] snpObjs = new SNP[m_gg.length];
			// load genotype data per dataset
			String refSNPAlleles = null;
			String refSNPAlleleAssessed = null;
			Boolean[] flip = new Boolean[m_gg.length];
			for (int d = 0; d < m_gg.length; d++) {
				Integer id = m_gg[d].getSnpToSNPId().get(snp);
				if (id >= 0) {
					SNP snpobj = m_gg[d].getSNPObject(id);
					loaders[d].loadGenotypes(snpobj);
					if (snpobj.getMAF() > settings.snpQCMAFThreshold && snpobj.getHWEP() > settings.snpQCHWEThreshold && snpobj.getCR() > settings.snpQCCallRateThreshold) {
						snpObjs[d] = snpobj;
						if (loaders[d].hasDosageInformation()) {
							loaders[d].loadDosage(snpObjs[d]);
						}
						if (refSNPAlleles == null) {
							refSNPAlleles = BaseAnnot.getAllelesDescription(snpobj.getAlleles());
							refSNPAlleleAssessed = BaseAnnot.toString(snpobj.getAlleles()[1]);
							flip[d] = false;
						} else {
							String SNPAlleles = BaseAnnot.getAllelesDescription(snpobj.getAlleles());
							String SNPAlleleAssessed = BaseAnnot.toString(snpobj.getAlleles()[1]);
							flip[d] = BaseAnnot.flipalleles(refSNPAlleles, refSNPAlleleAssessed, SNPAlleles, SNPAlleleAssessed);
							if (flip[d] == null) {
								snpobj.clearGenotypes();
							}
						}
					} else {
						snpobj.clearGenotypes();
					}
				}
			}

			// flip alleles

			for (int d = 0; d < data.length; d++) {
				data[d] = -1;
			}

			// get data per dataset
			for (int d = 0; d < m_gg.length; d++) {
				if (flip[d] != null) {
					boolean flipalleles = flip[d];
					String[] individuals = m_gg[d].getIndividuals();
					HashSet<String> allowedInds = allowedSamplesPerDs.get(d);

					double[] dosages = null;
					if (snpObjs[d].hasDosageInformation()) {
						dosages = snpObjs[d].getDosageValues();
					} else {
						byte[] genotypes = snpObjs[d].getGenotypes();
						dosages = new double[genotypes.length];
						for (int g = 0; g < dosages.length; g++) {
							dosages[g] = genotypes[g];
						}
					}


					for (int i = 0; i < individuals.length; i++) {
						String ind = individuals[i];
						if (allowedInds.contains(ind)) {
							Integer newId = sampleIdMap.get(ind);
							if (newId != null) {
								if (dosages[i] == -1) {
                                    data[newId] = -1;
								} else {
									if (flipalleles) {
										data[newId] = 2 - dosages[i];
									} else {
										data[newId] = dosages[i];
									}
								}
							} else {
								System.out.println("Sample " + ind + " has null id in index, but is present in GTE for dataset?");
								System.exit(-1);
							}
						}
					}
				} else if (snpObjs[d] != null) {
					System.out.println("Excluding\t" + snp + "\tin dataset\t" + datasetsettings.get(d).name + "\tsince it has incompatible alleles: " + BaseAnnot.getAllelesDescription(snpObjs[d].getAlleles()));
				}
			}

			int missing = 0;
			for (int d = 0; d < data.length; d++) {
				if (data[d] == -1) {
					missing++;
				}
			}
			if (missing == data.length) {
				System.out.println();
				System.out.println("Excluding\t" + snp + "\tsince it has no values");
			} else {
				tf.writeln(snp + "\t" + refSNPAlleles + "\t" + refSNPAlleleAssessed + "\t" + Strings.concat(data, Strings.tab));
			}
			pb.set(snpctr);
		}
		pb.close();

		tf.close();
	}
}
