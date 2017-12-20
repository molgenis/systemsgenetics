/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.textmeta;

import eqtlmappingpipeline.util.QTLFileSorter;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.concurrent.CompletionService;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.logging.Level;
import java.util.logging.Logger;

import umcg.genetica.console.ProgressBar;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.EQTL;
import umcg.genetica.io.trityper.QTLTextFile;
import umcg.genetica.text.Strings;

/**
 * @author harmjan
 */
public class FixedEffectMetaAnalysis {
	
	public void run(String filesDir, String output, Integer minimalNrDatasets, Integer minimalNrSamples) throws IOException {
		if (filesDir == null || output == null) {
			throw new IllegalArgumentException("Both input and output directory should be set!");
		}
		
		if (minimalNrDatasets == null) {
			minimalNrDatasets = 0;
		}
		if (minimalNrSamples == null) {
			minimalNrSamples = 0;
		}
		
		filesDir = Gpio.formatAsDirectory(filesDir);
		output = Gpio.formatAsDirectory(output);
		
		// do we need to convert the probe identifiers?
		// no. do this before running this tool
		// read each eQTL file, assess which eQTLs are in which files..
		// load: z-score, alleles (+ assessed), sample size
		if (!Gpio.exists(filesDir)) {
			throw new IllegalArgumentException("Input directory does not seem to exist!");
		}
		
		String[] filesInDir = Gpio.getListOfFiles(filesDir, "txt");
		
		if (filesInDir.length == 0) {
			filesInDir = Gpio.getListOfFiles(filesDir, "gz");
			if (filesInDir.length == 0) {
				System.err.println("No parseable files found in directory: " + filesDir);
				System.exit(0);
			}
		}
		
		for (String s : filesInDir) {
			System.out.println("Found text file:\t" + s);
		}
		
		Gpio.createDir(output);
		
		
		// obviously, we want to iterate through permuted files as well at some point,
		// although for the current analysis this is left out.
		EQTL[][] allEQTLs = new EQTL[filesInDir.length][0];
		HashSet<Pair<String, String>> uniqueSNPProbeCombos = new HashSet<Pair<String, String>>();
		HashSet<Pair<String, String>> sharedSNPProbeCombos = new HashSet<Pair<String, String>>();
		for (int f = 0; f < filesInDir.length; f++) {
			String fileName = filesInDir[f];
			TextFile tf = new TextFile(fileName, TextFile.R);
			String header = tf.readLine();
			String[] elems = tf.readLineElemsReturnObjects(TextFile.tab);
			ArrayList<EQTL> eqtls = new ArrayList<EQTL>();
			int ctr = 0;
			
			while (elems != null) {
				
				if (elems.length > 10) {
					
					EQTL e = new EQTL();
					String snp = elems[1];
					String probe = elems[4];
					Pair<String, String> p = new Pair<String, String>(snp, probe);
					if (uniqueSNPProbeCombos.contains(p)) {
						sharedSNPProbeCombos.add(p);
					} else {
						uniqueSNPProbeCombos.add(p);
					}
					
					
					// file may contain multiple datasets
					String[] datasetname = Strings.semicolon.split(elems[QTLTextFile.DATASETNAMES]);
					Double overallZscore = Double.parseDouble(elems[QTLTextFile.METAZ]);
					Integer[] samplesizes = new Integer[datasetname.length];
					
					Double[] datasetZ = new Double[datasetname.length];
					
					String[] samplesizestr = Strings.semicolon.split(elems[QTLTextFile.DATASETSIZE]);
					String[] datasetZstr = Strings.semicolon.split(elems[QTLTextFile.DATASETZSCORE]);
					for (int i = 0; i < datasetname.length; i++) {
						samplesizes[i] = Integer.parseInt(samplesizestr[i]);
						datasetZ[i] = Double.parseDouble(datasetZstr[i]);
					}
					
					String alleleAssessed = elems[QTLTextFile.ASESSEDALLELE];
					String alleles = elems[QTLTextFile.ASESSEDALLELE - 1];
					
					e.setAlleles(alleles);
					e.setAlleleAssessed(alleleAssessed);
					e.setZscore(overallZscore);
					e.setDatasetZScores(datasetZ);
					e.setDatasets(datasetname);
					e.setDatasetsSamples(samplesizes);
					e.setProbe(probe);
					e.setRsName(snp);
					
					eqtls.add(e);
					ctr++;
				}
				elems = tf.readLineElemsReturnObjects(TextFile.tab);
			}
			tf.close();
			allEQTLs[f] = eqtls.toArray(new EQTL[0]);
			System.out.println(eqtls.size() + " QTLs loaded from file: " + fileName);
		}
		
		System.out.println(uniqueSNPProbeCombos.size() + " unique SNP-probe combinations");
		System.out.println(sharedSNPProbeCombos.size() + " SNP-probe combinations shared with > 1 dataset");
		
		// iterate through all eQTLs
		// and just assume that snp-probe combinations are unique for each dataset.
		TextFile outfile = new TextFile(output + "eQTLs.txt", TextFile.W);
		outfile.writeln(QTLTextFile.header);
		int eqctr = 0;
		
		int nrprocs = Runtime.getRuntime().availableProcessors();
		ExecutorService threadPool = Executors.newFixedThreadPool(nrprocs);
		CompletionService<String> pool = new ExecutorCompletionService<String>(threadPool);
		
		HashMap<String, Integer> eqtlindex = new HashMap<String, Integer>();
		for (int i = 0; i < allEQTLs.length; i++) {
			for (int j = 0; j < allEQTLs[i].length; j++) {
				eqtlindex.put(i + "-" + allEQTLs[i][j].getRsName() + "-" + allEQTLs[i][j].getProbe(), j);
			}
		}
		
		int submitted = 0;
		for (Pair<String, String> eqtl : uniqueSNPProbeCombos) {
			FixedEffectMetaAnalysisTask t = new FixedEffectMetaAnalysisTask(eqtlindex, eqtl, filesInDir, allEQTLs, minimalNrDatasets, minimalNrSamples);
			pool.submit(t);
//            ArrayList<EQTL> eqtls = new ArrayList<EQTL>();
//            String snp = eqtl.getLeft();
//            String probe = eqtl.getRight();
//            for (int f = 0; f < filesInDir.length; f++) {
//                for (int e = 0; e < allEQTLs[f].length; e++) {
//                    EQTL eq = allEQTLs[f][e];
//                    if (eq.getRsName().equals(snp) && eq.getProbe().equals(probe)) {
//                        eqtls.add(eq);
//                    }
//                }
//            }
//
//            System.out.println(eqctr + "\t" + snp + "\t" + probe + "\t" + eqtls.size());
//            eqctr++;
//            // meta-analyze the collected EQTLs
//            double[] zscores = new double[eqtls.size()];
//            int[] samplesize = new int[eqtls.size()];
//            int nrSamples = 0;
//            String[] datsets = new String[eqtls.size()];
//            EQTL first = null;
//            if (eqtls.size() >= 3) {
//                for (int q = 0; q < eqtls.size(); q++) {
//                    // if this is not the first eQTL
//                    // check whether we should flip the allele...
//                    EQTL e = eqtls.get(q);
//                    Boolean flipZ = false;
//                    if (q > 0) {
//                        flipZ = BaseAnnot.flipalleles(first.getAlleles(), first.getAlleleAssessed(), e.getAlleles(), e.getAlleleAssessed());
//                        if (flipZ == null) {
//                            System.err.println("ERROR: alleles not compatible! " + e.getRsName() + "\t" + first.getAlleles() + "\t" + e.getAlleles());
//                        }
//                    } else {
//                        first = e;
//                    }
//
//                    // flip the allele if required
//                    if (flipZ) {
//                        zscores[q] = -e.getZscore();
//                    } else {
//                        zscores[q] = e.getZscore();
//                    }
//
//                    samplesize[q] = e.getDatasetsSamples()[0];
//                    nrSamples += samplesize[q];
//                }
//
//                // calculate meta statistics
//                double metaZ = ZScores.getWeightedZ(zscores, samplesize);
//                double pvalue = ZScores.zToP(metaZ);
//
//                // format
//                // PValue  SNPName SNPChr  SNPChrPos       ProbeName       ProbeChr        ProbeCenterChrPos       CisTrans        SNPType AlleleAssessed  OverallZScore   DatasetsWhereSNPProbePairIsAvailableAndPassesQC DatasetsZScores DatasetsNrSamples       IncludedDatasetsMeanProbeExpression     IncludedDatasetsProbeExpressionVariance HGNCName        IncludedDatasetsCorrelationCoefficient
//
//                String outStr =
//                        pvalue + "\t"
//                        + snp + "\t-\t-\t"
//                        + probe + "\t-\t-\ttrans\t"
//                        + eqtls.get(0).getAlleles() + "\t"
//                        + eqtls.get(0).getAlleleAssessed() + "\t"
//                        + metaZ + "\t"
//                        + Strings.concat(datsets, Strings.comma) + "\t"
//                        + Strings.concat(zscores, Strings.comma) + "\t" + Strings.concat(samplesize, Strings.comma) + "\t-\t-";
//
//                outfile.writeln(outStr);
//            }
			submitted++;
		}
		System.out.println(submitted + " eQTLs meta-analyzing");
		
		int received = 0;
		ProgressBar pb = new ProgressBar(submitted, "Running meta-analysis.");
		while (received < submitted) {
			try {
				Future<String> future = pool.take();
				String outputStr = future.get();
				if (outputStr != null) {
					outfile.writeln(outputStr);
				}
				received++;
				pb.set(received);
//                System.out.println(received + " / " + submitted);
			} catch (InterruptedException e) {
				e.printStackTrace();
			} catch (ExecutionException ex) {
				Logger.getLogger(FixedEffectMetaAnalysis.class.getName()).log(Level.SEVERE, null, ex);
			}
			
		}
		
		pb.close();
		outfile.close();
		
		System.out.println("Done. Now sorting results");
		QTLFileSorter sorter = new QTLFileSorter();
		sorter.run(output + "eQTLs.txt", output + "eQTLs_sorted.txt");
		if (Gpio.exists(output + "eQTLs_sorted.txt")) {
			Gpio.moveFile(output + "eQTLs_sorted.txt", output + "eQTLs.txt");
		}
		
	}
}
