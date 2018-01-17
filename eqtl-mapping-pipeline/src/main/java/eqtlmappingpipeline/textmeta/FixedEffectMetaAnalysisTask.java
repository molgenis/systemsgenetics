/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.textmeta;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.concurrent.Callable;

import umcg.genetica.containers.Pair;
import umcg.genetica.io.trityper.EQTL;
import umcg.genetica.io.trityper.util.BaseAnnot;
import umcg.genetica.math.stats.ZScores;
import umcg.genetica.text.Strings;
import umcg.genetica.util.Primitives;

/**
 * @author harmjan
 */
public class FixedEffectMetaAnalysisTask implements Callable<String> {
	
	Pair<String, String> eqtl;
	String[] filesInDir;
	EQTL[][] allEQTLs;
	private final HashMap<String, Integer> index;
	private final int minimalNrDatasets;
	private final int minimalNrSamples;
	
	public FixedEffectMetaAnalysisTask(HashMap<String, Integer> index, Pair<String, String> eqtl, String[] files, EQTL[][] allEQTLs, int minimalNrDatasets, int minimalNrSamples) {
		this.eqtl = eqtl;
		this.index = index;
		this.filesInDir = files;
		this.allEQTLs = allEQTLs;
		this.minimalNrDatasets = minimalNrDatasets;
		this.minimalNrSamples = minimalNrSamples;
	}
	
	@Override
	public String call() throws Exception {
		ArrayList<EQTL> eqtls = new ArrayList<EQTL>();
		String snp = eqtl.getLeft();
		String probe = eqtl.getRight();
		for (int f = 0; f < filesInDir.length; f++) {
			Integer id = index.get(f + "-" + snp + "-" + probe);
			if (id != null) {
				EQTL eq = allEQTLs[f][id];
				eqtls.add(eq);
			}
			
		}
		
		// meta-analyze the collected EQTLs
//        double[] zscores = new double[eqtls.size()];
//        int[] samplesize = new int[eqtls.size()];
		ArrayList<Double> zScores = new ArrayList<>();
		ArrayList<Integer> sampleSizes = new ArrayList<>();
		ArrayList<String> datasetNames = new ArrayList<>();
		
		int nrSamples = 0;
		
		EQTL first = null;
		if (eqtls.size() >= minimalNrDatasets) {
			for (int q = 0; q < eqtls.size(); q++) {
				// if this is not the first eQTL
				// check whether we should flip the allele...
				EQTL e = eqtls.get(q);
				Boolean flipZ = false;
				if (q > 0) {
					flipZ = BaseAnnot.flipalleles(first.getAlleles(), first.getAlleleAssessed(), e.getAlleles(), e.getAlleleAssessed());
					if (flipZ == null) {
						System.err.println("ERROR: alleles not compatible! " + e.getRsName() + "\t" + first.getAlleles() + "\t" + e.getAlleles());
					}
				} else {
					first = e;
				}
				
				// flip the allele if required
				Double[] dsZscores = e.getDatasetZScores();
				if (flipZ != null) {
					for (int d = 0; d < dsZscores.length; d++) {
						Double dz = dsZscores[d];
						if (dz != null) {
							if (flipZ) {
								zScores.add(-dz);
							} else {
								zScores.add(dz);
							}
							sampleSizes.add(e.getDatasetsSamples()[d]);
							nrSamples += e.getDatasetsSamples()[d];
							datasetNames.add(e.getDatasets()[d]);
						}
						
					}
				}
			}
			
			if (nrSamples >= minimalNrSamples) {
				// calculate meta statistics
				double[] zscores = Primitives.toPrimitiveArr(zScores.toArray(new Double[0]));
				int[] samplesizes = Primitives.toPrimitiveArr(sampleSizes.toArray(new Integer[0]));
				double metaZ = ZScores.getWeightedZ(zscores, samplesizes);
				double pvalue = ZScores.zToP(metaZ);
				
				// format
				// PValue  SNPName SNPChr  SNPChrPos       ProbeName       ProbeChr        ProbeCenterChrPos       CisTrans        SNPType AlleleAssessed  OverallZScore   DatasetsWhereSNPProbePairIsAvailableAndPassesQC DatasetsZScores DatasetsNrSamples       IncludedDatasetsMeanProbeExpression     IncludedDatasetsProbeExpressionVariance HGNCName        IncludedDatasetsCorrelationCoefficient
				
				String outStr =
						pvalue + "\t"
								+ snp + "\t-\t-\t"
								+ probe + "\t-\t-\ttrans\t"
								+ eqtls.get(0).getAlleles() + "\t"
								+ eqtls.get(0).getAlleleAssessed() + "\t"
								+ metaZ + "\t"
								+ Strings.concat(datasetNames, Strings.semicolon) + "\t"
								+ Strings.concat(zscores, Strings.semicolon) + "\t" + Strings.concat(samplesizes, Strings.semicolon) + "\t-\t-\t-\t-\t-\t-\t-";
				
				return outStr;
			} else {
				return null;
			}
		}
		return null;
	}
}
