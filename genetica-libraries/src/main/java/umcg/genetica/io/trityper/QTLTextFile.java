/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.io.trityper;

import umcg.genetica.collections.ChrPosTreeMap;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.regex.Pattern;

/**
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
		
		for (Iterator<EQTL> it = getEQtlIterator(); it.hasNext(); ) {
			alEQTLS.add(it.next());
		}
		
		return alEQTLS;
	}
	
	public EQTL[] readExpectedSize(int expSize) throws IOException {
		
		ArrayList<EQTL> alEQTLS = new ArrayList<EQTL>(expSize);
		
		for (Iterator<EQTL> it = getEQtlIterator(); it.hasNext(); ) {
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
		
		for (Iterator<EQTL> eqtlIterator = this.getEQtlIterator(); eqtlIterator.hasNext(); ) {
			EQTL qtl = eqtlIterator.next();
			
			ArrayList<EQTL> thisPosQtls = qtlTreeMap.get(qtl.getRsChr().toString(), qtl.getRsChrPos());
			if (thisPosQtls == null) {
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
			
			EQTL e = EQTL.fromString(elems, nullStr, separator);
			
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
