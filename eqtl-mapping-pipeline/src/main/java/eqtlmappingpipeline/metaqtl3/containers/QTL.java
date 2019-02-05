/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.metaqtl3.containers;

import cern.colt.matrix.tint.IntMatrix2D;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.TriTyperGeneticalGenomicsDataset;
import umcg.genetica.io.trityper.util.BaseAnnot;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Locale;

/**
 * @author harmjan
 */
public class QTL implements Comparable<QTL> {
	
	private double pvalue = Double.MAX_VALUE;
	private int pid = -1;
	private int sid = -1;
	private byte alleleAssessed;
	private double zscore = 0;
	private byte[] alleles;
	private double[] datasetZScores;
	private int[] datasetsSamples;
	private double[] correlations;
	private double[] datasetfc;
	private double[] datasetbeta;
	private double[] datasetbetase;
	private double finalbeta;
	private double finalbetase;
	
	public QTL(int datasets) {
		alleles = null;
		datasetZScores = null;
		datasetsSamples = null;
		correlations = null;
	}
	
	public QTL() {
	}
	
	public QTL(double pval, int pid, int sid, byte assessedAllele, double zscore, byte[] alleles, double[] zscores, int[] numSamples, double[] correlations, double[] fc, double[] beta, double[] betase, double finalbeta, double finalbetase) {
		this.pvalue = pval;
		this.pid = pid;
		this.sid = sid;
		this.alleleAssessed = assessedAllele;
		this.zscore = zscore;
		this.alleles = alleles;
		this.datasetZScores = zscores;
		this.datasetsSamples = numSamples;
		this.correlations = correlations;
		this.datasetfc = fc;
		this.datasetbeta = beta;
		this.datasetbetase = betase;
		this.finalbeta = finalbeta;
		this.finalbetase = finalbetase;
	}
	
	@Override
	public int compareTo(QTL o) {
		if (pvalue == o.pvalue) {
			if (Math.abs(zscore) == Math.abs(o.zscore)) {
				if (sid == o.sid) {
					if (pid == o.pid) {
						return 0;
					} else if (pid < o.pid) {
						return 1;
					} else {
						return -1;
					}
				} else if (sid < o.sid) {
					return 1;
				} else {
					return -1;
				}
			} else if (Math.abs(zscore) < Math.abs(o.zscore)) {
				return 1;
			} else {
				return -1;
			}
		} else if (pvalue > o.pvalue) {
			return 1;
		} else {
			return -1;
		}
		
	}
	
	public boolean equals(QTL o) {
		if (pvalue == o.pvalue) {
			if (Math.abs(zscore) == Math.abs(o.zscore)) {
				if (sid == o.sid) {
					if (pid == o.pid) {
						return true;
					} else {
						return false;
					}
				} else {
					return false;
				}
			} else {
				return false;
			}
		} else {
			return false;
		}
	}
	
	public void cleanUp() {
		
		alleles = null;
		if (datasetZScores != null) {
			for (int i = 0; i < datasetZScores.length; i++) {
				datasetZScores[i] = Double.NaN;
			}
			datasetZScores = null;
		}
		if (datasetsSamples != null) {
			for (int i = 0; i < datasetsSamples.length; i++) {
				datasetsSamples[i] = -9;
			}
			datasetsSamples = null;
		}
		
		if (correlations != null) {
			for (int i = 0; i < correlations.length; i++) {
				correlations[i] = Double.NaN;
			}
			correlations = null;
		}
	}
	
	public double getPvalue() {
		return pvalue;
	}
	
	public double getZscore() {
		return zscore;
	}
	
	public double[] getCorrelations() {
		return correlations;
	}
	
	private final static DecimalFormat format = new DecimalFormat("###.#######", new DecimalFormatSymbols(Locale.US));
	private final static DecimalFormat smallFormat = new DecimalFormat("0.#####E0", new DecimalFormatSymbols(Locale.US));
	
	public String getDescription(WorkPackage[] workPackages, IntMatrix2D probeTranslation, TriTyperGeneticalGenomicsDataset[] gg, int maxCisDistance) {
		if (sid == -1 && pid == -1) {
			return null;
		}
		
		String sepStr = ";";
		String nullstr = "-";
		char tabStr = '\t';
		
		StringBuilder out = new StringBuilder();
		
		if (pvalue < 1E-4) {
			out.append(smallFormat.format(pvalue));
		} else {
			out.append(format.format(pvalue));
		}
		
		out.append(tabStr);
		
		String rsName = null;
		String rsChr = nullstr;
		String rsChrPos = nullstr;
		
		WorkPackage currentWP = workPackages[sid];
		SNP[] snps = workPackages[sid].getSnps();
		
		for (int d = 0; d < snps.length; d++) {
			if (snps[d] != null) {
				rsName = snps[d].getName();
				rsChr = String.valueOf(snps[d].getChr());
				rsChrPos = String.valueOf(snps[d].getChrPos());
				break;
			}
		}
		
		
		String probe = nullstr;
		String probeChr = nullstr;
		String probeChrPos = nullstr;
		
		for (int d = 0; d < snps.length; d++) {
			if (probeTranslation.get(d, pid) != -9) {
				int probeId = probeTranslation.get(d, pid);
				probe = gg[d].getExpressionData().getProbes()[probeId];
				probeChr = String.valueOf(gg[d].getExpressionData().getChr()[probeId]);
				probeChrPos = String.valueOf((gg[d].getExpressionData().getChrStart()[probeId] + gg[d].getExpressionData().getChrStop()[probeId]) / 2);
				break;
			}
		}
		
		out.append(rsName);
		out.append(tabStr);
		out.append(rsChr);
		out.append(tabStr);
		out.append(rsChrPos);
		out.append(tabStr);
		out.append(probe);
		out.append(tabStr);
		out.append(probeChr);
		out.append(tabStr);
		out.append(probeChrPos);
		out.append(tabStr);
		
		String eQTLType = "trans";
		if (!rsChr.equals(nullstr) && !probeChr.equals(nullstr) && !probeChrPos.equals(nullstr) && !rsChrPos.equals(nullstr)) {
			if (rsChr.equals(probeChr) && Math.abs(Integer.parseInt(probeChrPos) - Integer.parseInt(rsChrPos)) < maxCisDistance) {
				eQTLType = "cis";
			}
		}
		
		out.append(eQTLType);
		out.append(tabStr);
		
		if (alleles == null) {
			System.err.println(rsName + " has null alleles..?'\n" + out.toString());
			return null;
		}
		
		out.append(BaseAnnot.toString(alleles[0])).append("/").append(BaseAnnot.toString(alleles[1]));
		out.append(tabStr);
		out.append(BaseAnnot.toString(alleleAssessed));
		out.append(tabStr);
		if (Math.abs(zscore) < 1E-4) {
			out.append(smallFormat.format(zscore));
		} else {
			out.append(format.format(zscore));
		}
		
		out.append(tabStr);
		
		String[] ds = new String[gg.length];
		Double[] probevars = new Double[gg.length];
		Double[] probemeans = new Double[gg.length];
		
		String hugo = nullstr;
		for (int d = 0; d < gg.length; d++) {
			if (!Double.isNaN(correlations[d])) {
				ds[d] = gg[d].getSettings().name;
				if (probeTranslation.get(d, pid) != -9) {
					int probeId = probeTranslation.get(d, pid);
					probevars[d] = gg[d].getExpressionData().getOriginalProbeVariance()[probeId];
					probemeans[d] = gg[d].getExpressionData().getOriginalProbeMean()[probeId];
					hugo = gg[d].getExpressionData().getAnnotation()[probeId];
				} else {
					System.out.println("ERROR!!!");
				}
			} else {
				ds[d] = null;
				probevars[d] = null;
				probemeans[d] = null;
			}
		}
		
		StringBuilder outcorrs = new StringBuilder();
		StringBuilder outzscores = new StringBuilder();
		StringBuilder outsamples = new StringBuilder();
		StringBuilder outmeans = new StringBuilder();
		StringBuilder outvars = new StringBuilder();
		StringBuilder outfc = new StringBuilder();
		StringBuilder outbeta = new StringBuilder();
		
		if (ds == null) {
			out.append(nullstr);
			out.append(tabStr);
			out.append(nullstr);
			out.append(tabStr);
			out.append(nullstr);
			out.append(tabStr);
			out.append(nullstr);
			out.append(tabStr);
			out.append(nullstr);
			out.append(tabStr);
			out.append(nullstr);
			out.append(tabStr);
			out.append(nullstr);
		} else {
			for (int d = 0; d < ds.length; d++) {
				if (d == 0) {
					sepStr = "";
				} else {
					sepStr = ";";
				}
				
				if (ds[d] == null) {
					out.append(sepStr).append(nullstr);
				} else {
					out.append(sepStr).append(ds[d]);
				}
				
				if (correlations == null || Double.isNaN(correlations[d])) {
					outcorrs.append(sepStr).append(nullstr);
				} else {
					DecimalFormat f = format;
					if (Math.abs(correlations[d]) < 1E-4) {
						f = smallFormat;
					}
					if (currentWP.getFlipSNPAlleles()[d]) {
						outcorrs.append(sepStr).append(f.format(-correlations[d]));
					} else {
						outcorrs.append(sepStr).append(f.format(correlations[d]));
					}
					
				}
				
				if (datasetZScores == null || Double.isNaN(datasetZScores[d])) {
					outzscores.append(sepStr).append(nullstr);
				} else {
					DecimalFormat f = format;
					if (Math.abs(datasetZScores[d]) < 1E-4) {
						f = smallFormat;
					}
					if (currentWP.getFlipSNPAlleles()[d]) {
						outzscores.append(sepStr).append(f.format(-datasetZScores[d]));
					} else {
						outzscores.append(sepStr).append(f.format(datasetZScores[d]));
					}
					
				}
				
				if (datasetsSamples == null || datasetsSamples[d] == -9) {
					outsamples.append(sepStr).append(nullstr);
				} else {
					outsamples.append(sepStr).append(datasetsSamples[d]);
				}
				
				if (probemeans == null || probemeans[d] == null) {
					outmeans.append(sepStr).append(nullstr);
				} else {
					outmeans.append(sepStr).append(format.format(probemeans[d]));
				}
				
				if (probevars == null || probevars[d] == null) {
					outvars.append(sepStr).append(nullstr);
				} else {
					outvars.append(sepStr).append(format.format(probevars[d]));
				}
				
				if (datasetfc == null || Double.isNaN(datasetfc[d])) {
					outfc.append(sepStr).append(nullstr);
				} else {
					outfc.append(sepStr).append(format.format(datasetfc[d]));
				}
				
				if (datasetbeta == null || Double.isNaN(datasetbeta[d])) {
					outbeta.append(sepStr).append(nullstr);
				} else {
					
					DecimalFormat f = format;
					if (Math.abs(datasetbeta[d]) < 1E-4) {
						f = smallFormat;
					}
					DecimalFormat f2 = format;
					if (Math.abs(datasetbetase[d]) < 1E-4) {
						f2 = smallFormat;
					}
					
					if (currentWP.getFlipSNPAlleles()[d]) {
						outbeta.append(sepStr).append((f.format(-datasetbeta[d]))).append(" (").append(f2.format(datasetbetase[d])).append(")");
					} else {
						outbeta.append(sepStr).append((f.format(datasetbeta[d]))).append(" (").append(f2.format(datasetbetase[d])).append(")");
					}
				}
				
			}
			
			out.append(tabStr);
			out.append(outzscores.toString());
			out.append(tabStr);
			out.append(outsamples.toString());
			out.append(tabStr);
			out.append(outmeans.toString());
			out.append(tabStr);
			out.append(outvars.toString());
			out.append(tabStr);
			out.append(hugo);
			out.append(tabStr);
			out.append(outcorrs.toString());
			out.append(tabStr);
			DecimalFormat f = format;
			if (Math.abs(finalbeta) < 1E-4) {
				f = smallFormat;
			}
			DecimalFormat f2 = format;
			if (Math.abs(finalbetase) < 1E-4) {
				f2 = smallFormat;
			}
			out.append(f.format(finalbeta)).append(" (").append(f2.format(finalbetase)).append(")");
			out.append(tabStr);
			out.append(outbeta.toString());
			out.append(tabStr);
			out.append(outfc.toString());
			out.append(tabStr);
			out.append(Double.NaN);
		}
		return out.toString();
		
	}
	
	public String getPermutationDescription(WorkPackage[] workPackages, IntMatrix2D probeTranslation, TriTyperGeneticalGenomicsDataset[] gg, int maxCisDistance) {
		if (sid == -1 && pid == -1) {
			return null;
		}
		
		String nullstr = "-";
		char tabStr = '\t';
		
		StringBuilder out = new StringBuilder();
		
		if (pvalue < 1E-4) {
			out.append(smallFormat.format(pvalue));
		} else {
			out.append(format.format(pvalue));
		}
//        out.append(pvalue);
		out.append(tabStr);
		
		String rsName = null;
		
		SNP[] snps = workPackages[sid].getSnps();
		
		for (int d = 0; d < snps.length; d++) {
			if (snps[d] != null) {
				rsName = snps[d].getName();
				break;
			}
		}
		
		String probe = nullstr;
		for (int d = 0; d < snps.length; d++) {
			if (probeTranslation.get(d, pid) != -9) {
				int probeId = probeTranslation.get(d, pid);
				probe = gg[d].getExpressionData().getProbes()[probeId];
				break;
			}
		}
		
		out.append(rsName);
		out.append(tabStr);
		
		out.append(probe);
		out.append(tabStr);
		
		String hugo = nullstr;
		for (int d = 0; d < gg.length; d++) {
			if (!Double.isNaN(correlations[d])) {
				if (probeTranslation.get(d, pid) != -9) {
					int probeId = probeTranslation.get(d, pid);
					hugo = gg[d].getExpressionData().getAnnotation()[probeId];
					if (hugo == null) {
						hugo = nullstr;
					}
				}
			}
		}
		
		out.append(hugo);
		out.append(tabStr);
		
		if (alleles == null) {
			System.err.println(rsName + " has null alleles..?'\n" + out.toString());
			return null;
		}
		
		out.append(BaseAnnot.toString(alleles[0])).append("/").append(BaseAnnot.toString(alleles[1]));
		out.append(tabStr);
		
		out.append(BaseAnnot.toString(alleleAssessed));
		out.append(tabStr);
		
		if (zscore < 1E-4) {
			out.append(smallFormat.format(zscore));
		} else {
			out.append(format.format(zscore));
		}
		
		return out.toString();
	}
}
