/* * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.metaqtl3;


import eqtlmappingpipeline.metaqtl3.containers.Settings;
import cern.colt.matrix.tint.IntMatrix2D;
import eqtlmappingpipeline.metaqtl3.containers.WorkPackage;

import java.io.IOException;
import java.util.concurrent.LinkedBlockingQueue;

import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.SNPLoader;
import umcg.genetica.io.trityper.TriTyperGeneticalGenomicsDataset;
import umcg.genetica.io.trityper.util.BaseAnnot;

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
		
		sumaveragesnpsize /= m_SNPLoaders.length;
		workPackageBufferSize = (int) Math.floor((double) (50 * 1048576) / sumaveragesnpsize);
		
		
		if (m_workPackages.length < workPackageBufferSize) {
			workPackageBufferSize = m_workPackages.length;
		}
		
		int workPackagesPassingQC = 0;
		int numProcessed = 0;
		
		TextFile snplog = null;
		try {
			if (!m_permuting) {
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
				StringBuilder[] qcBuffer = null;
				
				if (!m_permuting) {
					qcBuffer = new StringBuilder[workPackageBufferSize];
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
				for (int d = 0; d < m_gg.length; d++) {
					
					SNPLoader loader = m_SNPLoaders[d];
					
					boolean dosageAvailable = loader.hasDosageInformation();
					
					// too bad, but we need sorting by SNP id for optimal sequential access of the RandomAccessFile handlers in the SNPLoader objects
//                    if (d > 0) {
//                    	for(WorkPackage p : workPackageBuffer){
//                    		p.setDatasetToSortSNPs(d);
//						}
//                        java.util.Arrays.sort(workPackageBuffer);
//                    }
					
					
					for (int i = 0; i < workPackageBufferSize; i++) {
						
						WorkPackage wp = workPackageBuffer[i];
						
						if (!m_permuting && qcBuffer[i] == null) {
							qcBuffer[i] = new StringBuilder();
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
							loader.loadGenotypes(dSNP);
							if (!dSNP.passesQC() || dSNP.getCR() < m_callratethreshold || dSNP.getMAF() < m_mafthreshold || dSNP.getHWEP() < m_hwethreshold || dSNP.getAlleleItr() > 2) {
								snps[d].setPassesQC(false);
							} else {
								short dsPassingQC = wp.getDatasetsPassingQC();
								dsPassingQC++;
								wp.setDatasetsPassingQC(dsPassingQC);
							}
							
							if (!m_permuting) {
								Integer snpid = m_gg[d].getGenotypeData().getSnpToSNPId().get(dSNP.getName());
								qcBuffer[i].append("\t").
										append(snpid).append("\t").append(BaseAnnot.getAllelesDescription(dSNP.getAlleles())).append("\t").
										append(dSNP.getGenotypeFreq()[0]).append(" (").append(BaseAnnot.toString(dSNP.getAlleles()[0])).append(BaseAnnot.toString(dSNP.getAlleles()[0])).append(")").append("\t").
										append(dSNP.getGenotypeFreq()[1]).append(" (").append(BaseAnnot.toString(dSNP.getAlleles()[0])).append(BaseAnnot.toString(dSNP.getAlleles()[1])).append(")").append("\t").
										append(dSNP.getGenotypeFreq()[2]).append(" (").append(BaseAnnot.toString(dSNP.getAlleles()[1])).append(BaseAnnot.toString(dSNP.getAlleles()[1])).append(")").append("\t").
										append(dSNP.getCR()).append("\t").append(dSNP.getMAF()).append("\t").append(dSNP.getHWEP()).append("\t").append(dSNP.passesQC());
							}
							
							if (!dSNP.passesQC()) {
								snps[d].clearGenotypes();
								snps[d] = null;
							}
						} else {
							if (!m_permuting) {
								qcBuffer[i].append("\tNA\t-\t-\t-\t-");
							}
						}
					}
					
					if (dosageAvailable) {
						for (int i = 0; i < workPackageBufferSize; i++) {
							WorkPackage wp = workPackageBuffer[i];
							SNP[] snps = wp.getSnps();
							SNP dSNP = snps[d];
							if (wp.getDatasetsPassingQC() > 0 && dSNP != null && dSNP.passesQC()) {
								loader.loadDosage(dSNP);
							}
						}
					}
					
					
				}
				
				
				// done QC-ing and parsing SNPs
				for (int i = 0; i < workPackageBufferSize; i++) {
					WorkPackage wp = workPackageBuffer[i];
					if (!m_permuting) {
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
						StringBuilder finalQCString = new StringBuilder().append(m_snpList[wp.getMetaSNPId()]).append(qcBuffer[i].toString());
						snplog.writeln(finalQCString.toString());
						
					}
					
					if ((!m_settings.confineSNPsToSNPsPresentInAllDatasets && wp.getDatasetsPassingQC() > 0) || (m_settings.confineSNPsToSNPsPresentInAllDatasets && wp.getDatasetsPassingQC() == m_gg.length) || (m_settings.confineSNPsToSNPsPresentInAllDatasets && m_permuting && wp.getDatasetsPassingQC() > 0)) {
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
		int firstDatasetToPassQC = -1;
		byte[] firstDatasetPassinQCAlleles = null;
		Boolean[] flipAlleles = new Boolean[m_gg.length];
		byte firstminor = -1;
		
		for (int d = 0; d < m_gg.length; d++) {
			SNP dSNP = snps[d];
			
			if (dSNP != null) {
				// check if the alleles are identical with previuously loaded SNP...
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
					
				} else {
					
					byte[] allelesToCompare = dSNP.getAlleles();
					int nrAllelesIdentical = 0;
					byte minor = dSNP.getMinorAllele();
					
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
							snplog.writeln(output);
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
