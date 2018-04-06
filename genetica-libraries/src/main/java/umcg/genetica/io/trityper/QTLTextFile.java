/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.io.trityper;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.regex.Pattern;
import umcg.genetica.collections.ChrPosTreeMap;
import umcg.genetica.io.text.TextFile;
import static umcg.genetica.io.text.TextFile.tab;
import umcg.genetica.io.trityper.util.ChrAnnotation;

/**
 *
 * @author harmjan
 */
public class QTLTextFile extends TextFile {

	public static int PVAL = 0;
	public static int SNP = 1;
	public static int SNPCHR = 2;
	public static int SNPLOC = 3;
	public static int PROBE = 4;
	public static int PROBECHR = 5;
	public static int PROBELOC = 6;
	public static int CISTRANS = 7;
	public static int SNPTYPE = 8;
	public static int ASESSEDALLELE = 9;
	public static int METAZ = 10;
	public static int DATASETNAMES = 11;
	public static int DATASETZSCORE = 12;
	public static int DATASETSIZE = 13;
	public static int HUGO = 16;
	public static int DATASECORR = 17;
	public static int METAB = 18;
	public static int DATASETB = 19;
	/*
     // PValue
     * SNPName
     * SNPChr
     * SNPChrPos
     * ProbeName
     * ProbeChr
     * ProbeCenterChrPos
     * CisTrans
     * SNPType
     * AlleleAssessed
     * OverallZScore
     * DatasetsWhereSNPProbePairIsAvailableAndPassesQC
     * DatasetsZScores
     * DatasetsNrSamples
     * IncludedDatasetsMeanProbeExpression
     * IncludedDatasetsProbeExpressionVariance
     * HGNCName
     * IncludedDatasetsCorrelationCoefficient
     * Meta-Beta (SE)
     * Beta (SE)
     * FoldChange
     * FDR
     * */
	private static String sepStr = ";";
	//private static String tabStr = "\t";
	private static String nullStr = "-";
	private static Pattern separator = Pattern.compile(sepStr);
	public static String header = "PValue\t"
			+ "SNPName\t"
			+ "SNPChr\t"
			+ "SNPChrPos\t"
			+ "ProbeName\t"
			+ "ProbeChr\t"
			+ "ProbeCenterChrPos\t"
			+ "CisTrans\t"
			+ "SNPType\t"
			+ "AlleleAssessed\t"
			+ "OverallZScore\t"
			+ "DatasetsWhereSNPProbePairIsAvailableAndPassesQC\t"
			+ "DatasetsZScores\t"
			+ "DatasetsNrSamples\t"
			+ "IncludedDatasetsMeanProbeExpression\t"
			+ "IncludedDatasetsProbeExpressionVariance\t"
			+ "HGNCName\t"
			+ "IncludedDatasetsCorrelationCoefficient\t"
			+ "Meta-Beta (SE)\t"
			+ "Beta (SE)\t"
			+ "FoldChange\t"
			+ "FDR";

	public QTLTextFile(String loc, boolean W) throws IOException {
		super(loc, W);
		if (W) {
			write(header + '\n');
		}
	}

	public QTLTextFile(String loc, boolean W, boolean gz) throws IOException {
		super(loc, W);
		if (W) {
			write(header + '\n');
		}
	}

	public void write(EQTL[] eqtllist) throws IOException {
		for (EQTL e : eqtllist) {
			write(e.toString() + '\n');
		}
	}

	public void write(ArrayList<EQTL> eqtllist) throws IOException {
		for (EQTL e : eqtllist) {
			write(e.toString() + '\n');
		}
	}

	public void writeMinimal(ArrayList<MinimalEQTL> eqtllist) throws IOException {
		for (MinimalEQTL e : eqtllist) {
			write(e.toString() + '\n');
		}
	}

	public EQTL[] read() throws IOException {

		return readExpectedSize(1000);

	}

	public ArrayList<EQTL> readList() throws IOException {

		ArrayList<EQTL> alEQTLS = new ArrayList<EQTL>();

		for (Iterator<EQTL> it = getEQtlIterator(); it.hasNext();) {
			alEQTLS.add(it.next());
		}

		return alEQTLS;
	}

	public EQTL[] readExpectedSize(int expSize) throws IOException {

		ArrayList<EQTL> alEQTLS = new ArrayList<EQTL>(expSize);

		for (Iterator<EQTL> it = getEQtlIterator(); it.hasNext();) {
			alEQTLS.add(it.next());
		}

		EQTL[] eqtls = new EQTL[alEQTLS.size()];
		return alEQTLS.toArray(eqtls);
	}

	public Iterator<EQTL> getEQtlIterator() throws IOException {
		open();
		return new EQtlIterator();
	}

	public ChrPosTreeMap<ArrayList<EQTL>> readQtlsAsTreeMap() throws IOException {
		
		ChrPosTreeMap<ArrayList<EQTL>> qtlTreeMap = new ChrPosTreeMap<>();
		
		for (Iterator<EQTL> eqtlIterator = this.getEQtlIterator(); eqtlIterator.hasNext();) {
			EQTL qtl = eqtlIterator.next();
			
			ArrayList<EQTL> thisPosQtls = qtlTreeMap.get(qtl.getRsChr().toString(), qtl.getRsChrPos());
			if(thisPosQtls == null){
				thisPosQtls = new ArrayList<>(1);
				qtlTreeMap.put(qtl.getRsChr().toString(), qtl.getRsChrPos(), thisPosQtls);
			}
			thisPosQtls.add(qtl);
			
		}
		
		return qtlTreeMap;
		
	}

	private class EQtlIterator implements Iterator<EQTL> {

		private String[] elems;
		private final boolean fdrpresent;

		public EQtlIterator() throws IOException {

			/*
             * 0 - pval 1 - rs 2 - rs chr 3 - rs chr pos 4 - probe 5 - probe chr 6 - probe chr center pos 7 - cis 8 - alleles 9 - allele assessed 10 - Z-score 11 -
             * dataset 12 - Z-score 13 - nr samples assessed 14 - probe mean 15 - probe variance 16 - probeHugo 17 - correlation 18 - meta-beta 19 - beta 20 - fc 21
             * - fdr
			 */
			elems = readLineElemsReturnReference(tab); // skip headerline

			if (elems.length > 21) {
				fdrpresent = true;
			} else {
				fdrpresent = false;
			}

			//read first line
			elems = readLineElemsReturnReference(tab);

		}

		@Override
		public boolean hasNext() {
			return elems != null;
		}

		@Override
		public EQTL next() {

			String[] subelems;

			EQTL e = new EQTL();
			if (!elems[0].equals(nullStr)) {
				e.setPvalue(Double.parseDouble(elems[0]));
			}

			if (!elems[1].equals(nullStr)) {
				e.setRsName(elems[1]);
			}

			if (!elems[2].equals(nullStr)) {
				e.setRsChr(ChrAnnotation.parseChr(elems[2]));
			}

			if (!elems[3].equals(nullStr)) {
				e.setRsChrPos(Integer.parseInt(elems[3]));
			}

			if (!elems[4].equals(nullStr)) {
				e.setProbe(elems[4]);
			}

			if (!elems[5].equals(nullStr)) {
				e.setProbeChr(ChrAnnotation.parseChr(elems[5]));
			}

			if (!elems[6].equals(nullStr)) {
				e.setProbeChrPos(Integer.parseInt(elems[6]));
			}

			if (!elems[7].equals(nullStr)) {
				e.setType(elems[7]);
			}

			if (!elems[8].equals(nullStr)) {
				e.setAlleles(elems[8]);
			}

			if (!elems[9].equals(nullStr)) {
				e.setAlleleAssessed(elems[9]);
			}

			if (!elems[10].equals(nullStr)) {
				e.setZscore(Double.parseDouble(elems[10]));
			}

			if (!elems[11].equals(nullStr)) {
				e.setDatasets(separator.split(elems[11]));
			}

			if (!elems[12].equals(nullStr)) {
				subelems = separator.split(elems[12]);
				Double[] dsZScores = new Double[subelems.length];
				for (int i = 0; i < subelems.length; i++) {
					try {
						dsZScores[i] = Double.parseDouble(subelems[i]);
					} catch (NumberFormatException ex) {
						dsZScores[i] = Double.NaN;
					}
				}
				e.setDatasetZScores(dsZScores);
			}

			if (!elems[13].equals(nullStr)) {
				String[] samples = elems[13].split(";");
				Integer[] sampleS = new Integer[samples.length];
				for (int i = 0; i < samples.length; i++) {
					try {
						sampleS[i] = Integer.parseInt(samples[i]);
					} catch (NumberFormatException ex) {
						sampleS[i] = null;
					}
				}

				e.setDatasetsSamples(sampleS);

			}

			if (!elems[14].equals(nullStr)) {
				subelems = separator.split(elems[14]);
				Double[] dsPmeans = new Double[subelems.length];
				for (int i = 0; i < subelems.length; i++) {
					try {
						dsPmeans[i] = Double.parseDouble(subelems[i]);
					} catch (NumberFormatException ex) {
						dsPmeans[i] = Double.NaN;
					}
				}
				e.setProbeMeans(dsPmeans);
			}

			if (!elems[15].equals(nullStr)) {
				subelems = separator.split(elems[15]);
				Double[] dsPvars = new Double[subelems.length];
				for (int i = 0; i < subelems.length; i++) {
					try {
						dsPvars[i] = Double.parseDouble(subelems[i]);
					} catch (NumberFormatException ex) {
						dsPvars[i] = Double.NaN;
					}
				}
				e.setProbeVariance(dsPvars);
			}

			if (!elems[16].equals(nullStr)) {
				e.setProbeHUGO(elems[16]);
			}

			if (!elems[17].equals(nullStr)) {
				subelems = separator.split(elems[17]);
				Double[] dsCorrs = new Double[subelems.length];
				for (int i = 0; i < subelems.length; i++) {
					try {
						dsCorrs[i] = Double.parseDouble(subelems[i]);
					} catch (NumberFormatException ex) {
						dsCorrs[i] = Double.NaN;
					}
				}
				e.setCorrelations(dsCorrs);
			}

			if (!elems[18].equals(nullStr)) {
				e.setMetaBeta(elems[18]);
			}

			if (!elems[19].equals(nullStr)) {
				e.setBeta(elems[19]);
			}

			if (elems.length > 20) {
				if (!elems[20].equals(nullStr)) {
					e.setFC(elems[20]);
				}
			}

			if (fdrpresent && elems.length > 21 && !elems[21].equals(nullStr)) {
				try {
					e.setFDR(Double.parseDouble(elems[21]));
				} catch (java.lang.NumberFormatException ex) {
					//do nothing
				}
			}
			try {
				elems = readLineElemsReturnReference(tab);
			} catch (IOException ex) {
				throw new RuntimeException(ex);
			}
			return e;

		}

		@Override
		public void remove() {
			throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
		}
	}
}
