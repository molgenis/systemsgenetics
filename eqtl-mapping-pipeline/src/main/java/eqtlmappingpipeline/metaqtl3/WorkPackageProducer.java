/* * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.metaqtl3;


import cern.colt.matrix.tint.IntMatrix2D;
import eqtlmappingpipeline.metaqtl3.containers.Settings;
import eqtlmappingpipeline.metaqtl3.containers.WorkPackage;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.SNPLoader;
import umcg.genetica.io.trityper.TriTyperGeneticalGenomicsDataset;
import umcg.genetica.io.trityper.util.BaseAnnot;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.stream.IntStream;

/**
 * @author harmjan
 */
class WorkPackageProducer extends Thread {

	private String[] m_snpList;
	//    private String[] m_probeList;
//    private IntMatrix2D m_probeTranslationTable;
//    private IntMatrix2D m_snpTranslationTable;
	private TriTyperGeneticalGenomicsDataset[] m_gg;
	private LinkedBlockingQueue<WorkPackage> m_queue;
	private WorkPackage[] m_workPackages;
	private SNPLoader[] m_SNPLoaders;
	private Settings m_settings;
	//    private int m_name, nrThreads;
	private double m_mafthreshold, m_hwethreshold, m_callratethreshold;
	public boolean done, semaphore;
	private boolean m_permuting;
	private final String m_outputdir;

	WorkPackageProducer(LinkedBlockingQueue<WorkPackage> packageQueue, WorkPackage[] workPackages, String[] snpList, String[] probeList, IntMatrix2D probeTranslationTable,
						IntMatrix2D snpTranslationTable, TriTyperGeneticalGenomicsDataset[] gg, SNPLoader[] snploaders, Settings settings, boolean permuting) {
		this.m_workPackages = workPackages;
		this.m_queue = packageQueue;
		this.m_snpList = snpList;
//        this.m_probeList = probeList;
//        this.m_probeTranslationTable = probeTranslationTable;
//        this.m_snpTranslationTable = snpTranslationTable;
		this.m_gg = gg;
		this.m_SNPLoaders = snploaders;
		this.m_settings = settings;
		this.m_mafthreshold = settings.snpQCMAFThreshold;
		this.m_hwethreshold = settings.snpQCHWEThreshold;
		this.m_callratethreshold = settings.snpQCCallRateThreshold;
		this.m_permuting = permuting;
		this.m_outputdir = settings.outputReportsDir;
	}

	@Override
	public void run() {
//        System.out.println("Starting production facility.");
		done = false;
		semaphore = false;

		// create workpackage objects: determine to which probe each snp should be mapped

		int workPackageBufferSize = 1;

		double sumaveragesnpsize = 0;
		try {
			for (int d = 0; d < m_SNPLoaders.length; d++) {
				sumaveragesnpsize += m_SNPLoaders[d].getAverageSNPSize(m_gg[d].getGenotypeData().getSNPs().length);
			}
		} catch (IOException e) {
			e.printStackTrace();
		}

		// make sure not to load more than 500mb worth of SNP data per block
		// SNPLoader has it's own buffer of a set number of variants
		sumaveragesnpsize /= m_SNPLoaders.length;
		workPackageBufferSize = 1; // 256; // (int) Math.floor((double) (16 * 1048576) / sumaveragesnpsize);


		if (m_workPackages.length < workPackageBufferSize) {
			workPackageBufferSize = m_workPackages.length;
		}

//		System.out.println("Loading " + workPackageBufferSize + " SNPs per buffer.");
//		System.out.println();
		int workPackagesPassingQC = 0;
		int numProcessed = 0;

//		workPackageBufferSize = 1;

		TextFile snplog = null;
		try {
			if (!m_permuting && m_settings.writeSNPQCLog) {
				snplog = new TextFile(m_outputdir + "SNPQCLog.txt.gz", TextFile.W);

				String ln = "-";
				// SNPId, MAF, HWE, CR, passesQC
				for (int d = 0; d < m_gg.length; d++) {
					ln += "\t" + m_gg[d].getSettings().name + "\t\t\t\t\t\t\t";
				}
				snplog.writeln(ln);

				ln = "SNP";
				for (int d = 0; d < m_gg.length; d++) {
					ln += "\tSNPId\tAlleles\tFreqAA\tFreqAB\tFreqBB\tCR\tMAF\tHWE\tPassesQC";
				}
				snplog.writeln(ln);
			}


			while (numProcessed < m_workPackages.length) {

				if (numProcessed + workPackageBufferSize > m_workPackages.length) {
					workPackageBufferSize = m_workPackages.length - numProcessed;
				}

				WorkPackage[] workPackageBuffer = new WorkPackage[workPackageBufferSize];
				StringBuilder[][] qcBuffer = null;

				if (!m_permuting && m_settings.writeSNPQCLog) {
					qcBuffer = new StringBuilder[workPackageBufferSize][m_gg.length];
				}

				// load a set of workpackages in the buffer
				int numInBuffer = 0;
				while ((numInBuffer < workPackageBufferSize) && (numProcessed < m_workPackages.length)) {
					WorkPackage toAdd = m_workPackages[numProcessed];
					if (toAdd != null) {
						workPackageBuffer[numInBuffer] = toAdd;
						numInBuffer++;
					}
					numProcessed++;
				}

				// load the SNPs for each dataset
				int finalWorkPackageBufferSize = workPackageBufferSize;
				StringBuilder[][] finalQcBuffer = qcBuffer;
				IntStream.range(0, m_gg.length).parallel().forEach(d -> {
					SNPLoader loader = m_SNPLoaders[d];

					boolean dosageAvailable = loader.hasDosageInformation();

					for (int i = 0; i < finalWorkPackageBufferSize; i++) {

						WorkPackage wp = workPackageBuffer[i];

						if (!m_permuting && m_settings.writeSNPQCLog && finalQcBuffer[i][d] == null) {
							finalQcBuffer[i][d] = new StringBuilder();
						}
						// update sorting dataset
//                        if (m_gg.length > 1) {
//                            wp.setDatasetToSortSNPs(d + 1);
//                        }

						SNP[] snps = wp.getSnps();
						SNP dSNP = snps[d];

//						if (i % 5000 == 0) {
//							System.out.println("d " + d + "\ti " + i + "\tnp: " + numProcessed + "/" + m_workPackages.length + "\t" + workPackageBuffer.length);
//						}
						if (dSNP != null) {
							try {
								loader.loadGenotypes(dSNP);
							} catch (IOException e) {
								e.printStackTrace();
							}
							if (!dSNP.passesQC() || dSNP.getCR() < m_callratethreshold || dSNP.getMAF() < m_mafthreshold || dSNP.getHWEP() < m_hwethreshold || dSNP.getAlleleItr() > 2) {
								snps[d].setPassesQC(false);
							} else {
								wp.incrementDatasetsPassingQC();
//								short dsPassingQC = wp.getDatasetsPassingQC();
//								dsPassingQC++;
//								wp.setDatasetsPassingQC(dsPassingQC);
							}

							if (!m_permuting && m_settings.writeSNPQCLog) {
								Integer snpid = m_gg[d].getGenotypeData().getSnpToSNPId().get(dSNP.getName());
								String allele1;
								String allele2;
								String alleleDesc;
								if (dSNP.hasAlleleEncoding()) {
									String[] alleleEncoding = dSNP.getAlleleEncoding();
									allele1 = alleleEncoding[0];
									allele2 = alleleEncoding[1];
									alleleDesc = Strings.concat(alleleEncoding, Strings.forwardslash);
								} else {
									allele1 = BaseAnnot.toString(dSNP.getAlleles()[0]);
									allele2 = BaseAnnot.toString(dSNP.getAlleles()[1]);
									alleleDesc = BaseAnnot.getAllelesDescription(dSNP.getAlleles());
								}
								finalQcBuffer[i][d].append("\t").
										append(snpid).append("\t").append(alleleDesc).append("\t").
										append(dSNP.getGenotypeFreq()[0]).append(" (").append(allele1).append(allele1).append(")").append("\t").
										append(dSNP.getGenotypeFreq()[1]).append(" (").append(allele1).append(allele2).append(")").append("\t").
										append(dSNP.getGenotypeFreq()[2]).append(" (").append(allele2).append(allele2).append(")").append("\t").
										append(dSNP.getCR()).append("\t").append(dSNP.getMAF()).append("\t").append(dSNP.getHWEP()).append("\t").append(dSNP.passesQC());
							}

							if (!dSNP.passesQC()) {
								snps[d].clearGenotypes();
								snps[d] = null;
							}
						} else {
							if (!m_permuting && m_settings.writeSNPQCLog) {
								finalQcBuffer[i][d].append("\tNA\t-\t-\t-\t-\t-\t-\t-\t-");
							}
						}
					}

					if (dosageAvailable) {
						for (int i = 0; i < finalWorkPackageBufferSize; i++) {
							WorkPackage wp = workPackageBuffer[i];
							SNP[] snps = wp.getSnps();
							SNP dSNP = snps[d];
							if (wp.getDatasetsPassingQC() > 0 && dSNP != null && dSNP.passesQC()) {
								try {
									loader.loadDosage(dSNP);
								} catch (IOException e) {
									e.printStackTrace();
								}
							}
						}
					}
				});

				qcBuffer = finalQcBuffer;

				// done QC-ing and parsing SNPs
				for (int i = 0; i < workPackageBufferSize; i++) {
					WorkPackage wp = workPackageBuffer[i];
					if (!m_permuting && m_settings.writeSNPQCLog) {
						String snpName = m_snpList[wp.getMetaSNPId()];
						SNP[] snps = wp.getSnps();
						boolean notequal = false;
						for (int d = 0; d < snps.length; d++) {
							SNP snp = snps[d];
							if (snp != null) {
								if (!snpName.equals(snp.getName())) {
									System.out.println("ERROR! SNP names not equal: " + d + "-" + snp.getName());
									notequal = true;
								}
							}
						}
						if (notequal) {
							System.exit(0);
						}
						StringBuilder finalQCString = new StringBuilder().append(m_snpList[wp.getMetaSNPId()]);
						for (int d = 0; d < m_gg.length; d++) {
							finalQCString.append(qcBuffer[i][d].toString());
						}
						snplog.writeln(finalQCString.toString());

					}

					int dsPassingQC = wp.getDatasetsPassingQC();
					if ((!m_settings.confineSNPsToSNPsPresentInAllDatasets && dsPassingQC > 0)
							|| (m_settings.confineSNPsToSNPsPresentInAllDatasets && dsPassingQC == m_gg.length)
							|| (m_settings.confineSNPsToSNPsPresentInAllDatasets && m_permuting && dsPassingQC > 0)
							|| (wp.getDatasetsPassingQC() >= m_settings.requireAtLeastNumberOfDatasets)) {
						// check whether alleles should be flipped.
						boolean allelesOk = detmermineAlleleFlips(wp, snplog);


						if (allelesOk) {
							// put the fully loaded WP in the queue for further processing...
							try {
								m_queue.put(wp);
								workPackagesPassingQC++;
							} catch (InterruptedException ex) {
								ex.printStackTrace();
							}
						}
					} else {
						wp = null;
					}
				}
			}

			if (!m_permuting && snplog != null) {
				snplog.close();
			}
		} catch (IOException e) {
			e.printStackTrace();
		}


		int remaining = m_queue.remainingCapacity();


//        System.out.println("Done loading work packages.");

		int prevSize = Integer.MAX_VALUE;
		int realsize = m_queue.size() - remaining;
		while (realsize > 0) {
			realsize = m_queue.size() - remaining;
			if (realsize < prevSize && realsize > 0) {
				if (m_settings.transAnalysis) {
//                    if(realsize % 100 == 0){
//                        System.out.println(realsize +" workpackages left in queue");
//                    }
				}
				prevSize = realsize;
			}
		}

		for (int i = 0; i < m_settings.nrThreads * 2; i++) {
			WorkPackage bomb = new WorkPackage();
			bomb.setIsKillPackage(true);
			m_queue.offer(bomb);
		}

	}

	private byte[] convertToComplementaryAlleles(byte[] allelesToCompare) {
		byte[] allelesComplementary = new byte[2];
		for (int a = 0; a < 2; a++) {
			allelesComplementary[a] = BaseAnnot.getComplement(allelesToCompare[a]);
		}
		return allelesComplementary;
	}

	private boolean detmermineAlleleFlips(WorkPackage wp, TextFile snplog) throws IOException {
		SNP[] snps = wp.getSnps();

		// check whether there are any indels
		int nrDsWithVars = 0;
		int nrDsWithIndels = 0;
		int firstds = -1;
		Boolean[] flipAlleles = new Boolean[m_gg.length];
		for (int d = 0; d < snps.length; d++) {
			SNP s = snps[d];
			if (s != null) {
				if (firstds == -1) {
					firstds = d;
				}
				if (s.isIndel()) {
					nrDsWithIndels++;
				}
				nrDsWithVars++;
			}
		}

		if (nrDsWithIndels > 0) {
			if (nrDsWithIndels != nrDsWithVars) {
				if (!m_permuting && m_settings.writeSNPQCLog) {
					String dsWithIndel = "";
					String dsWithOutIndel = "";

					for (int d = 0; d < snps.length; d++) {
						SNP s = snps[d];
						if (s != null) {
							if (s.isIndel()) {
								if (dsWithIndel.length() > 0) {
									dsWithIndel += ";" + m_gg[d].getSettings().name;
								} else {
									dsWithIndel += m_gg[d].getSettings().name;
								}
							} else {
								if (dsWithOutIndel.length() > 0) {
									dsWithOutIndel += ";" + m_gg[d].getSettings().name;
								} else {
									dsWithOutIndel += m_gg[d].getSettings().name;
								}

							}
						}
					}

					String output = "Variant: " + snps[firstds].getChr() + ":" + snps[firstds].getChrPos() + ":" + snps[firstds].getName() + " is an indel in datasets: " + dsWithIndel + " but not " + dsWithOutIndel;
					System.err.println(output);
					snplog.writeln(output);
				}
				return false;
			} else {
				// do indel comparison
				// no fancy pancy alignment here. let's just check whether the strings are equal
				String[] reference = null;
				byte referenceEffectAllele = 0;
				for (int d = 0; d < snps.length; d++) {
					SNP s = snps[d];
					if (s != null) {
						if (s.isIndel()) {
							if (reference == null) {
								reference = s.getAlleleEncoding();
								// determine if allele 1 is minor
								referenceEffectAllele = s.getMinorAllele();
								// if not, flip alleles
								if (referenceEffectAllele == s.getAlleles()[0]) {
									flipAlleles[d] = true;
								} else {
									flipAlleles[d] = false;
								}

							} else {
								// compare against reference
								int nrequal = 0;
								String[] enc = s.getAlleleEncoding();

								for (int b = 0; b < enc.length; b++) {
									for (int q = 0; q < reference.length; q++) {
										if (enc[b].equals(reference[q])) {
											nrequal++;
											break;
										}
									}
								}
								if (nrequal < 2) {
									String output = "Variant: " + s.getChr() + ":" + s.getChrPos() + ":" + s.getName()
											+ " is an indel and has different alleles in dataset: " + m_gg[d].getSettings().name + " (" + Strings.concat(enc, Strings.comma) + ") " +
											"compared to (" + Strings.concat(reference, Strings.comma) + ")";
//									System.err.println(output);
									snplog.writeln(output);
									return false;
								} else {

								}
							}
						}
					}
				}
			}
		} else {

			// SNP comparison
			int firstDatasetToPassQC = -1;
			byte[] firstDatasetPassinQCAlleles = null;

			byte firstminor = -1;

			for (int d = 0; d < m_gg.length; d++) {
				SNP dSNP = snps[d];

				if (dSNP != null) {
					// check if the alleles are identical with previously loaded SNP...
					if (firstDatasetToPassQC == -1) {
						firstDatasetToPassQC = d;
						firstDatasetPassinQCAlleles = dSNP.getAlleles();
						byte minor = dSNP.getMinorAllele();
						if (firstDatasetPassinQCAlleles[1] == minor) {
							flipAlleles[d] = false;
						} else {
							flipAlleles[d] = true;
						}

						firstminor = minor;

						if (dSNP.hasAlleleEncoding()) {
							// replace byte codes
							String[] enc = dSNP.getAlleleEncoding();
							for (int a = 0; a < firstDatasetPassinQCAlleles.length; a++) {
								byte b = firstDatasetPassinQCAlleles[a];
								firstDatasetPassinQCAlleles[a] = BaseAnnot.toByte(enc[b - 100]);
							}
							firstminor = BaseAnnot.toByte(enc[minor - 100]);
						}
					} else {

						byte[] allelesToCompare = dSNP.getAlleles();
						byte minor = dSNP.getMinorAllele();
						if (dSNP.hasAlleleEncoding()) {
							// replace byte codes
							String[] enc = dSNP.getAlleleEncoding();
							for (int a = 0; a < allelesToCompare.length; a++) {
								byte b = allelesToCompare[a];
								allelesToCompare[a] = BaseAnnot.toByte(enc[b - 100]);
							}
							minor = BaseAnnot.toByte(enc[minor - 100]);
						}

						int nrAllelesIdentical = 0;

						boolean flipalleles = false;
						int minorAlleleNum = 0;
						if (allelesToCompare[0] != minor) {
							minorAlleleNum = 1;
						}

						for (int a = 0; a < 2; a++) {
							for (int b = 0; b < 2; b++) {
								if (firstDatasetPassinQCAlleles[a] == allelesToCompare[b]) {
									nrAllelesIdentical++;
								}
							}
						}

						if (nrAllelesIdentical != 2) {
							//Alleles are different, take complimentary:
							allelesToCompare = convertToComplementaryAlleles(allelesToCompare);
							minor = BaseAnnot.getComplement(minor);
						}

						nrAllelesIdentical = 0;

						for (int a = 0; a < 2; a++) {
							for (int b = 0; b < 2; b++) {
								if (firstDatasetPassinQCAlleles[a] == allelesToCompare[b]) {
									nrAllelesIdentical++;
								}
							}
						}

						if (nrAllelesIdentical != 2) {
							if (!m_permuting) {
								String snp1Alleles = BaseAnnot.toString(firstDatasetPassinQCAlleles[0]) + "/" + BaseAnnot.toString(firstDatasetPassinQCAlleles[1]);
								String snp2Alleles = BaseAnnot.toString(allelesToCompare[0]) + "/" + BaseAnnot.toString(allelesToCompare[1]);
								String output = "SNP alleles are not identical between datasets for SNP: " + wp.getSnps()[d].getName()
										+ "\tSNP1 (" + m_gg[firstDatasetToPassQC].getSettings().name + "): " + snp1Alleles
										+ "\tSNP2 (" + m_gg[d].getSettings().name + "): " + snp2Alleles;
								System.err.println(output);
								if (m_settings.writeSNPQCLog) {
									snplog.writeln(output);
								}
							}
							return false;
						} else {
							if (minor != firstminor) {
								// error or warning or whatever
								//                        System.out.println("WARNING: minor allele is different for identical SNP: "+dSNP.getName() + ", probably due to high MAF.\nWill conform to allelic direction of dataset: "+m_gg[firstDatasetToPassQC].getSettings().name);
								//                        double[] allelefreq = dSNP.getAlleleFreq();
								//                        byte[] origAlleles = snps[firstDatasetToPassQC].getAlleles();
								//                        double[] origAlleleFreq = snps[firstDatasetToPassQC].getAlleleFreq();
								//                        System.out.println("Reference MAF:"+snps[firstDatasetToPassQC].getMAF()+"\tAssessed MAF:"+dSNP.getMAF());
								//                        for(int i=0; i<2; i++){
								//                            System.out.println("ref ds: "+m_gg[firstDatasetToPassQC].getSettings().name+"\t"+BaseAnnot.toString(origAlleles[i])+"\t("+origAlleleFreq[i]+")\tAssessed: "+m_gg[d].getSettings().name+"\t"+BaseAnnot.toString(allelesToCompare[i])+"\t("+allelefreq[i]+")");
								//                        }
								//                        System.out.println("");
								// take the orientation of the first dataset..., which is dataset
								if (minorAlleleNum == 1) {
									flipalleles = true;
								} else {
									flipalleles = false;
								}

							} else {
								if (allelesToCompare[0] == minor) {
									flipalleles = true;
								} else {
									flipalleles = false;
								}
							}
							flipAlleles[d] = flipalleles;
						}
					}
				}
			}
		}
		wp.setFlipSNPAlleles(flipAlleles);

		return true;
	}

	void kill() {
		done = true;
	}

	WorkPackage[] getWorkPackages() {
		return m_workPackages;
	}
}
