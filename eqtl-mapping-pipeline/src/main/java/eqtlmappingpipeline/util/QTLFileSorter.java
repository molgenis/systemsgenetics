

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.util;

import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.QTLTextFile;
import umcg.genetica.io.trityper.util.ChrAnnotation;
import umcg.genetica.text.Strings;

import java.io.File;
import java.io.IOException;
import java.util.*;


/**
 * @author harmjan
 */
public class QTLFileSorter {

	public enum SORTBY {
		P, Z, POS;
	}

	public void run(String efile, String outfile, SORTBY s) throws IOException {
		run(efile, outfile, 2500000, s);
	}

	public void run(String efile, String outfile, int n, SORTBY s) throws IOException {
		System.out.println("Sorting: " + efile + " by " + s + " in batches of " + n);

		TextFile tf = new TextFile(efile, TextFile.R);
		int cap = n;
		ArrayList<QTLObj> eQtls = new ArrayList<>(cap);
		int ctr = 0;
		int total = 0;
		tf.readLine();
		String ln = tf.readLine();
		while (ln != null) {
			if (eQtls.size() == cap) {
				if (eQtls != null) {
					System.out.println("Sorting: " + eQtls.size() + " to batch " + ctr);
					QTLObj[] qtlarr = eQtls.toArray(new QTLObj[0]);
					Arrays.parallelSort(qtlarr, new QTLComparator(s));
					System.out.println("Writing: " + eQtls.size() + " to batch " + ctr + ": " + outfile + "-tmp-" + ctr + ".txt.gz");
					TextFile out = new QTLTextFile(outfile + "-tmp-" + ctr + ".txt.gz", QTLTextFile.W);
					for (QTLObj obj : qtlarr) {
						out.writeln(obj.ln);
					}
					out.close();
					ctr++;
				}
				eQtls = new ArrayList<>();
			}


			String[] elems = ln.split("\t");
			double z = Math.abs(Double.parseDouble(elems[10]));
			if (!Double.isNaN(z) && Double.isFinite(z)) {

				byte chr = ChrAnnotation.parseChr(elems[2]);
				Integer pos = Integer.parseInt(elems[3]);
				double p = Double.parseDouble(elems[0]);
				QTLObj obj = new QTLObj(z, p, chr, pos, ln, s);
				eQtls.add(obj);
				total++;
			} else {
//                System.out.println("Z == NaN: " + z + "\t" + ln);
			}

			ln = tf.readLine();
		}

		tf.close();


		if (eQtls.size() > 0) {
			Collections.sort(eQtls, new QTLComparator(s));
			System.out.println("Writing: " + eQtls.size() + " to batch " + ctr);
			TextFile out = new QTLTextFile(outfile + "-tmp-" + ctr + ".txt.gz", QTLTextFile.W);
			for (QTLObj obj : eQtls) {
				out.writeln(obj.ln);
			}
			out.close();
			ctr++;
		}
		eQtls = null;

		System.out.println(total + " eqtls over " + ctr + " batches");

		// merge batch files
		TextFile[] tfs = new TextFile[ctr];
		String[][] lastlnelems = new String[ctr][];
		double[] statArr = new double[ctr];
		for (int c = 0; c < ctr; c++) {
			System.out.println("Opening: " + outfile + "-tmp-" + c + ".txt.gz");
			tfs[c] = new TextFile(outfile + "-tmp-" + c + ".txt.gz", QTLTextFile.R);
			tfs[c].readLine();
			lastlnelems[c] = tfs[c].readLineElems(TextFile.tab);
			if (s.equals(SORTBY.P)) {
				statArr[c] = Math.abs(Double.parseDouble(lastlnelems[c][0]));
			} else if (s.equals(SORTBY.Z)) {
				statArr[c] = Math.abs(Double.parseDouble(lastlnelems[c][10]));
			} else {
				statArr[c] = Math.abs(Double.parseDouble(lastlnelems[c][3]));
			}
		}

		boolean done = false;
		QTLTextFile out = new QTLTextFile(outfile, QTLTextFile.W);

		int written = 0;
		while (!done) {

			if (s.equals(SORTBY.P)) {
				double maxP = 1;
				// determine max ln
				Integer maxc = null;
				for (int c = 0; c < ctr; c++) {
					if (lastlnelems[c] != null) {
						double p = statArr[c];
						if (!Double.isNaN(p) && p <= maxP) {
							maxc = c;
							maxP = p;
						}
					}
				}

				// write selected line
				if (maxc != null) {
					out.writeln(Strings.concat(lastlnelems[maxc], Strings.tab));
					written++;
					lastlnelems[maxc] = tfs[maxc].readLineElems(TextFile.tab);
					if (lastlnelems[maxc] != null) {
						statArr[maxc] = Math.abs(Double.parseDouble(lastlnelems[maxc][0]));
					} else {
						statArr[maxc] = Double.NaN;
					}
				} else {
					done = true;
				}
				if (written % 1000000 == 0) {
					System.out.println(written + " eqtls written.");
				}
			} else if (s.equals(SORTBY.Z)) {

				double maxZ = 0;
				// determine max ln
				Integer maxc = null;
				for (int c = 0; c < ctr; c++) {
					if (lastlnelems[c] != null) {
						double z = statArr[c];
						if (!Double.isNaN(z) && z >= maxZ) {
							maxc = c;
							maxZ = z;
						}
					}
				}

				// write selected line
				if (maxc != null) {
					out.writeln(Strings.concat(lastlnelems[maxc], Strings.tab));
					written++;
					lastlnelems[maxc] = tfs[maxc].readLineElems(TextFile.tab);
					if (lastlnelems[maxc] != null) {
						statArr[maxc] = Math.abs(Double.parseDouble(lastlnelems[maxc][10]));
					} else {
						statArr[maxc] = Double.NaN;
					}
				} else {
					done = true;
				}
				if (written % 1000000 == 0) {
					System.out.println(written + " eqtls written.");
				}
			} else {
				double maxZ = Double.MAX_VALUE;

				// determine max ln
				Integer maxc = null;
				for (int c = 0; c < ctr; c++) {
					if (lastlnelems[c] != null) {
						double z = statArr[c];
						if (!Double.isNaN(z) && z <= maxZ) {
							maxc = c;
							maxZ = z;
						}
					}
				}

				// write selected line
				if (maxc != null) {
					out.writeln(Strings.concat(lastlnelems[maxc], Strings.tab));
					written++;
					lastlnelems[maxc] = tfs[maxc].readLineElems(TextFile.tab);
					if (lastlnelems[maxc] != null) {
						statArr[maxc] = Math.abs(Double.parseDouble(lastlnelems[maxc][3]));
					} else {
						statArr[maxc] = Double.NaN;
					}
				} else {
					done = true;
				}
				if (written % 1000000 == 0) {
					System.out.println(written + " eqtls written.");
				}
			}
		}
		out.close();

		for (int c = 0; c < ctr; c++) {
			tfs[c].close();
		}

		for (int c = 0; c < ctr; c++) {
			File f = new File(outfile + "-tmp-" + c + ".txt.gz");
			f.delete();
		}
		System.out.println("Done sorting");
	}


	class QTLObj {
		double z;
		String ln;
		SORTBY s;
		int pos;
		byte chr;
		double p;

		public QTLObj(double z, double p, byte chr, int pos, String ln, SORTBY s) {
			this.z = z;
			this.ln = ln;
			this.s = s;
			this.chr = chr;
			this.pos = pos;
			this.p = p;
		}

		@Override
		public boolean equals(Object o) {

			if (s.equals(SORTBY.P)) {
				if (this == o) return true;
				if (o == null || getClass() != o.getClass()) return false;
				QTLObj QTLZObj = (QTLObj) o;
				return Double.compare(QTLZObj.p, p) == 0;
			} else if (s.equals(SORTBY.Z)) {

				if (this == o) return true;
				if (o == null || getClass() != o.getClass()) return false;
				QTLObj QTLZObj = (QTLObj) o;
				int comp = Double.compare(QTLZObj.z, z);
				if (comp == 0) {
					return Double.compare(QTLZObj.p, p) == 0;
				} else {
					return false;
				}
			} else {
				if (this == o) return true;
				if (o == null || getClass() != o.getClass()) return false;

				QTLObj qtlPosObj = (QTLObj) o;

				if (pos == qtlPosObj.pos && chr == qtlPosObj.chr) return true;
				return false;
			}
		}

		@Override
		public int hashCode() {
			if (s.equals(SORTBY.P)) {
				return Objects.hash(p);
			} else if (s.equals(SORTBY.Z)) {
				return Objects.hash(z);
			} else {
				int result = pos;
				result = 31 * result + chr;
				return result;
			}
		}
	}

	class QTLComparator implements Comparator<QTLObj> {
		SORTBY s;

		QTLComparator(SORTBY s) {
			this.s = s;
		}

		@Override
		public int compare(QTLObj o1, QTLObj o2) {
			if (s.equals(SORTBY.P)) {
				if (o1.equals(o2)) {
					// compare metabeta
					return 0;
				}
				double p1 = o1.p;
				double p2 = o2.p;
				if (p1 < p2) {
					return -1;
				} else if (p1 > p2) {
					return 1;
				} else {
					return 0;
				}
			} else if (s.equals(SORTBY.Z)) {
				if (o1.equals(o2)) {
					// compare metabeta
					return 0;
				}

				double z1 = o1.z;
				double z2 = o2.z;
				if (z1 > z2) {
					return -1;
				} else if (z1 < z2) {
					return 1;
				} else {
					if (o1.p == o2.p) {
						return 0;
					} else if (o1.p > o2.p) {
						return 1;
					} else {
						return -1;
					}
				}
			} else {
				if (o1.equals(o2)) {
					return 0;
				}

				if (o1.chr > o2.chr) {
					return -1;
				} else if (o1.chr < o2.chr) {
					return 1;
				} else {
					if (o1.pos > o2.pos) {
						return 1;
					} else if (o1.pos < o2.pos) {
						return -1;
					}
				}
			}
			return 0;
		}
	}


}
