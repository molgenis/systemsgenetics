package eqtlmappingpipeline.util;

import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.util.BaseAnnot;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;

public class GetMAFFromSNPQCLog {

    public static void main(String[] args) {
        String file = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-08-11-MAFCompare\\SNPQCLog.txt.gz";
        String file2 = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-08-11-MAFCompare\\SNPQCLog-test.txt.gz";
        GetMAFFromSNPQCLog q = new GetMAFFromSNPQCLog();
        try {
            q.run(file, file2);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void run(String qclog, String out) throws IOException {
        TextFile tf = new TextFile(qclog, TextFile.R);
        String[] dsnameElems = tf.readLine().split("\t"); // header 1
        ArrayList<String> dsNames = new ArrayList<String>();
        for (int i = 1; i < dsnameElems.length; i++) {
            if (dsnameElems[i].length() > 0) {
                dsNames.add(dsnameElems[i]);
            }
        }

        tf.readLine(); // header 2
        TextFile tfo = new TextFile(out, TextFile.W);
        TextFile tfo2 = new TextFile(out + "-AFPerDS.txt.gz", TextFile.W);
        String outheader = "SNP\tAlleles\tAlleleA\tNumA\tNumB\tAF(A)\tAF(B)\tMinorAllele\tMAF";
        String dsheader = "SNP\tAlleles\tAF-A\t" + Strings.concat(dsNames, Strings.tab);
        tfo2.writeln(dsheader);
        tfo.writeln(outheader);
        String[] elems = tf.readLineElems(TextFile.tab);
        int numelemsperds = 9;
//		while (elems != null) {
//			String snp = elems[0];
//
//			ArrayList<HashMap<Character, Integer>> ctPerAlleleArr = new ArrayList<HashMap<Character, Integer>>();
//
//			for (int i = 1; i < elems.length; i += numelemsperds) {
//				String alleles = elems[i + 1];
//				HashMap<Character, Integer> ctPerAllele = new HashMap<Character, Integer>();
//				for (int q = i + 2; q < i + 5; q++) {
//					String[] freqelems = elems[q].split(" ");
//					if (!freqelems[q].equals("-") && freqelems.length > 1) {
//						Integer ct = Integer.parseInt(freqelems[0]);
//						String allele = freqelems[1];
//						allele = allele.replace("(", "");
//						allele = allele.replace(")", "");
//
//						char a0 = allele.charAt(0);
//						Integer cta = ctPerAllele.get(a0);
//						if (cta == null) {
//							cta = ct;
//						} else {
//							cta += ct;
//						}
//						ctPerAllele.put(a0, cta);
//					}
//					ctPerAlleleArr.add(ctPerAllele);
//				}
//			}
//
//
//
//			elems = tf.readLineElems(TextFile.tab);
//		}

        int nrElems = elems.length;
        int nrheader = nrElems - 1;
        int numds = nrheader / 9;
        int lnctr = 0;
        int[] freqdist = new int[51];
        int[][] freqdistperds = new int[51][numds];
        while (elems != null) {
            try {
                String snp = elems[0];

                double[] afAds = new double[numds];
                Boolean[] flipalleles = new Boolean[numds];
                Boolean[] passqc = new Boolean[numds];
                String refAlleles = null;
                String refAlleleA = null;
                String refAlleleB = null;
                int ds = 0;
                for (int i = 1; i < elems.length; i += numelemsperds) {
                    String freqstr = elems[i + 2];
                    String passqcstr = elems[i + 8];
                    if (!freqstr.equals("-")) {
                        Boolean b = Boolean.parseBoolean(passqcstr);
                        passqc[ds] = b;
                        String alleles = elems[i + 1];
                        String[] allelesSplit = alleles.split("/");
                        String alleleA = allelesSplit[0];
                        String alleleB = allelesSplit[1];
                        if (b != null && b && refAlleles == null) {
                            refAlleles = alleles;
                            refAlleleA = alleleA;
                            refAlleleB = alleleB;
                            flipalleles[ds] = false;
                        } else if (b != null && b) {
                            flipalleles[ds] = BaseAnnot.flipalleles(refAlleles, refAlleleA, alleles, alleleA);
                        }
                    }
                    ds++;
                }

                int numA = 0;
                int numB = 0;
                int total = 0;

                ds = 0;
                // loop over datasets
                for (int i = 1; i < elems.length; i += numelemsperds) {


                    if (flipalleles[ds] != null && passqc[ds] != null && passqc[ds]) {

                        int actr = 0;
                        int numAds = 0;
                        int numBds = 0;
                        int totalds = 0;

                        // loop over alleles
                        for (int q = i + 2; q < i + 5; q++) {
                            String[] freqelems = Strings.whitespace.split(elems[q]);
                            if (!freqelems[0].equals("-") && freqelems.length > 1) {
                                Integer ct = Integer.parseInt(freqelems[0]);
                                total += (2 * ct);
                                totalds += (2 * ct);
                                switch (actr) {
                                    case 0: // AA
                                        numAds += (2 * ct);
                                        if (flipalleles[ds]) {
                                            numB += (2 * ct);
                                        } else {
                                            numA += (2 * ct);
                                        }
                                        break;
                                    case 1: // AB
                                        numAds += (ct);
                                        numBds += (ct);
                                        numA += ct;
                                        numB += ct;
                                        break;
                                    case 2: // BB
                                        numBds += (2 * ct);
                                        if (flipalleles[ds]) {
                                            numA += (2 * ct);
                                        } else {
                                            numB += (2 * ct);
                                        }
                                        break;
                                }
                            }
                            actr++;
                        }

                        double percAds = (double) numAds / totalds;
                        double percBds = (double) numBds / totalds;
                        double mafds = percAds;
                        if (percAds > percBds) {
                            mafds = percBds;
                        }
//						System.out.println(ds + "\t" + mafds);
                        int bin = (int) Math.floor(100 * mafds);
                        if (bin > 50) {
                            System.out.println("Error 1: " + mafds);
                        }
                        freqdistperds[bin][ds]++;
                        if(flipalleles[ds]){
                            afAds[ds] = percBds;
                        } else {
                            afAds[ds] = percAds;
                        }
                    } else {
                        afAds[ds] = Double.NaN;
                    }

                    ds++;
                }
//				System.exit(-1);

                double afA = ((double) numA / total);
                double afB = ((double) numB / total);
                String minorAllele = refAlleleA;
                double maf = afA;
                if (numA > numB) {
                    minorAllele = refAlleleB;
                    maf = afB;

                }

                String outln = snp + "\t" + refAlleles + "\t" + refAlleleA + "\t" + numA + "\t" + numB + "\t" + afA + "\t" + afB + "\t" + minorAllele + "\t" + maf;
                String outln2 = snp + "\t" + refAlleles + "\t" + afA + "\t" + Strings.concat(afAds, Strings.tab);
                tfo.writeln(outln);
                tfo2.writeln(outln2);
                int bin = (int) Math.floor(100 * maf);
                if (bin > 50) {
                    System.out.println("Error 2: " + maf);
                }
                freqdist[bin]++;
                // System.out.println(lnctr + "\t" + outln);
            } catch (ArrayIndexOutOfBoundsException e) {
                System.out.println("Problem with line " + lnctr + "\t" + e.getMessage());
                System.out.println(Strings.concat(elems, Strings.tab));
//				System.exit(-1);
            }
            elems = tf.readLineElems(TextFile.tab);
            lnctr++;
            if (lnctr % 100000 == 0) {
                System.out.print(lnctr + " lines parsed\r");
            }
        }
        System.out.println();

        tfo.close();
        tfo2.close();
        tf.close();
        System.out.println();
        for (int i = 0; i < freqdist.length; i++) {
            String ln = i + "\t" + freqdist[i];
            for (int d = 0; d < numds; d++) {
                ln += "\t" + freqdistperds[i][d];
            }
            System.out.println(ln);
        }
    }
}
