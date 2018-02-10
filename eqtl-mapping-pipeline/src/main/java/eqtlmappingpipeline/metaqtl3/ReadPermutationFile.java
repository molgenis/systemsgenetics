package eqtlmappingpipeline.metaqtl3;

import gnu.trove.map.hash.TDoubleIntHashMap;
import gnu.trove.set.hash.THashSet;
import umcg.genetica.console.MultiThreadProgressBar;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.QTLTextFile;
import umcg.genetica.text.Strings;

import java.util.concurrent.Callable;

public class ReadPermutationFile implements Callable<TDoubleIntHashMap> {
	
	private final String permutationDir;
	private final int permutationRound;
	private final FDR.FileFormat f;
	private final int maxNrMostSignificantEQTLs;
	private final FDR.FDRMethod m;
	private final MultiThreadProgressBar mpb;
	
	public ReadPermutationFile(String permutationDir, int permutationRound, FDR.FileFormat f, int maxNrMostSignificantEQTLs, FDR.FDRMethod m, MultiThreadProgressBar pb) {
		this.permutationDir = permutationDir;
		this.permutationRound = permutationRound;
		this.f = f;
		this.maxNrMostSignificantEQTLs = maxNrMostSignificantEQTLs;
		this.m = m;
		this.mpb = pb;
	}
	
	@Override
	public TDoubleIntHashMap call() throws Exception {
		TDoubleIntHashMap permutedPvalues = new TDoubleIntHashMap(10000, 0.5f);
		String fileString = permutationDir + "/PermutedEQTLsPermutationRound" + (permutationRound + 1) + ".txt.gz";
//		System.out.println(fileString);
		// read the permuted eqtl output
		if (mpb != null) {
			mpb.setSubtasks(permutationRound, maxNrMostSignificantEQTLs);
		}
		TextFile gz = new TextFile(fileString, TextFile.R, 10 * 1048576);
		
		String[] header = gz.readLineElems(TextFile.tab);
		int snpcol = -1;
		int pvalcol = -1;
		int probecol = -1;
		int genecol = -1;
		
		
		if (f == FDR.FileFormat.REDUCED) {
			
			//PValue  SNP     Probe   Gene
			for (int col = 0; col < header.length; col++) {
				if (header[col].equals("PValue")) {
					pvalcol = col;
				}
				if (header[col].equals("SNP")) {
					snpcol = col;
				}
				if (header[col].equals("Probe")) {
					probecol = col;
				}
				if (header[col].equals("Gene")) {
					genecol = col;
				}
			}
			
			//PValue  SNP     Probe   Gene
			if (snpcol == -1 || pvalcol == -1 || probecol == -1 && genecol == -1) {
				System.out.println("Column not found in permutation file: " + fileString);
				System.out.println("PValue: " + pvalcol);
				System.out.println("SNP: " + snpcol);
				System.out.println("Probe: " + probecol);
				System.out.println("Gene: " + genecol);
			}
		}

//            String[] data = gz.readLineElemsReturnReference(TextFile.tab);
		int itr = 0;
		
		THashSet<String> visitedEffects = new THashSet<String>();
		String permln = gz.readLine();
		int lnctr = 0;
		while (permln != null) {
			if (permln.length() != 0) {
				if (itr > maxNrMostSignificantEQTLs - 1) {
					System.out.println("Breaking because: " + itr);
					break;
				} else {
					String fdrId = null;
					if (f == FDR.FileFormat.REDUCED) {
						if (m == FDR.FDRMethod.PROBELEVEL) {
							// fdrId = data[probecol];
							fdrId = Strings.subsplit(permln, Strings.tab, probecol, probecol + 1)[0];
						} else if (m == FDR.FDRMethod.SNPLEVEL) {
//								fdrId = data[snpcol];
							fdrId = Strings.subsplit(permln, Strings.tab, snpcol, snpcol + 1)[0];
						} else if (m == FDR.FDRMethod.GENELEVEL && header.length > 3) {
//								fdrId = data[genecol];
							fdrId = Strings.subsplit(permln, Strings.tab, genecol, genecol + 1)[0];
						}
						if (fdrId != null) {
							fdrId = new String(fdrId.getBytes("UTF-8")).intern();
						}
					} else {
						if (m == FDR.FDRMethod.GENELEVEL) {
//								fdrId = data[QTLTextFile.HUGO];
							fdrId = Strings.subsplit(permln, Strings.tab, QTLTextFile.HUGO, QTLTextFile.HUGO + 1)[0];
						} else if (m == FDR.FDRMethod.SNPLEVEL) {
//								fdrId = data[QTLTextFile.SNP];
							fdrId = Strings.subsplit(permln, Strings.tab, QTLTextFile.SNP, QTLTextFile.SNP + 1)[0];
						} else if (m == FDR.FDRMethod.PROBELEVEL) {
//								fdrId = data[4];
							fdrId = Strings.subsplit(permln, Strings.tab, 4, 5)[0];
						}
						if (fdrId != null) {
							fdrId = new String(fdrId.getBytes("UTF-8")).intern();
						}
					}
					
					// take top effect per gene / probe
					if (m == FDR.FDRMethod.FULL || (!fdrId.equals("-") && !visitedEffects.contains(fdrId))) {
						
						if (m != FDR.FDRMethod.FULL) {
							visitedEffects.add(fdrId);
						}
						
						double permutedP = Double.parseDouble(Strings.subsplit(permln, Strings.tab, 0, 1)[0]);
						if (permutedPvalues.containsKey(permutedP)) {
							permutedPvalues.increment(permutedP);
						} else {
							permutedPvalues.put(permutedP, 1);
						}
						
						itr++;
						if (mpb != null) {
							mpb.set(permutationRound, itr);
						}
					}
					lnctr++;
//					if (lnctr % 1000000 == 0) {
//						System.out.println(fileString + "\t" + lnctr + " lines parsed.");
//					}
					
					permln = gz.readLine();
				}
			}
		}
		gz.close();
		System.out.println("Done reading: " + fileString + "\t" + permutedPvalues.keys().length + " unique p-values. " + visitedEffects.size() + " IDs seen.");
		if(m.equals(FDR.FDRMethod.GENELEVEL)){
			System.out.println("check this out");
		}
		if (mpb != null) {
			mpb.complete(permutationRound);
		}
		return permutedPvalues;
	}
}