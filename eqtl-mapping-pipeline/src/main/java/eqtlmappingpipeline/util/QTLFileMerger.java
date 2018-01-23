package eqtlmappingpipeline.util;

import eqtlmappingpipeline.metaqtl3.containers.TextResult;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import javax.xml.soap.Text;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;

public class QTLFileMerger {
	
	public void mergeChr(String loc, String out, int nrqtl) throws IOException {
		ArrayList<String> locs = new ArrayList<>();
		for (int chr = 1; chr < 23; chr++) {
			String f = loc.replaceAll("CHR", "" + chr);
			if (Gpio.exists(f)) {
				locs.add(f);
				System.out.println("Found: " + f);
			}
		}
		merge(locs.toArray(new String[0]), out, nrqtl);
	}
	
	public void merge(String[] listOfQTLFiles, String outloc, int nrQTL) throws IOException {
		
		QTLObject[] r = null;
		
		QTLObject[] rtmp = new QTLObject[nrQTL];
		
		String outputHeader = null;
		int qtlctr = 0;
		double highestP = 2;
		
		for (String s : listOfQTLFiles) {
			System.out.println("Parsing: " + s);
			TextFile tf = new TextFile(s, TextFile.R);
			String fileHeader = tf.readLine();
			if (outputHeader == null) {
				outputHeader = fileHeader;
			}
			String[] origFileHeaderElems = outputHeader.split("\t");
			String[] fileHeaderElems = fileHeader.split("\t");
			int zscoreCol = -1;
			if (origFileHeaderElems.length != fileHeaderElems.length) {
				System.err.println("Error: header of file doesn't match header of first file: " + s);
				System.out.println("This is a very naive file merger ;)");
				System.exit(-1);
			} else {
				for (int i = 0; i < fileHeaderElems.length; i++) {
					if (fileHeaderElems[i].toLowerCase().equals("zscore") || fileHeaderElems[i].toLowerCase().equals("overallzscore")) {
						zscoreCol = i;
					}
				}
			}
			
			String ln = tf.readLine();
			int nrinfile = 0;
			int nrincluded = 0;
			while (ln != null) {
				
				String[] elems = Strings.tab.split(ln);
				
				double p = Double.parseDouble(elems[0]);
				if (p <= highestP) {
					double z = Double.parseDouble(elems[zscoreCol]);
					QTLObject obj = new QTLObject(p, z, ln);
					rtmp[qtlctr] = obj;
					nrincluded++;
					qtlctr++;
					if (qtlctr == rtmp.length) {
						if (r == null) {
							r = rtmp;
							System.out.println();
							System.out.println("Merging results..");
							rtmp = new QTLObject[nrQTL];
							highestP = r[r.length - 1].p;
							System.out.println("Setting initial highest P: " + highestP);
							qtlctr = 0;
						} else {
							System.out.println();
							System.out.println("Merging results..");
							QTLObject[] rtmp2 = new QTLObject[r.length + qtlctr];
							System.arraycopy(r, 0, rtmp2, 0, r.length);
							System.arraycopy(rtmp, 0, rtmp2, r.length, qtlctr);
							Arrays.parallelSort(rtmp2);
							System.arraycopy(rtmp2, 0, r, 0, r.length);
							highestP = r[r.length - 1].p;
							System.out.println("Setting highest P: " + highestP);
							qtlctr = 0;
						}
					}
					
				}
				nrinfile++;
				if (nrinfile % 10000 == 0) {
					System.out.print("\r" + nrinfile + " parsed. " + nrincluded + " included sofar. Current buffer position: " + qtlctr);
				}
				ln = tf.readLine();
			}
			tf.close();
			System.out.println("Done");
			System.out.println(nrinfile + " QTL in file. " + nrincluded + " made it into the list (for now).\tHighest P: " + highestP + "\tRemaining in buffer: " + qtlctr);
		}
		System.out.println();
		System.out.println("Done reading files. " + qtlctr + " QTL remain in buffer.");
		// solve last kliekjes
		if (qtlctr > 0) {
			if (r == null) {
				r = new QTLObject[qtlctr];
				System.arraycopy(rtmp, 0, r, 0, qtlctr);
				Arrays.parallelSort(r);
			} else {
				QTLObject[] rtmp2 = new QTLObject[r.length + qtlctr];
				System.arraycopy(r, 0, rtmp2, 0, r.length);
				System.arraycopy(rtmp, 0, rtmp2, r.length, qtlctr);
				Arrays.parallelSort(rtmp2);
				System.arraycopy(rtmp2, 0, r, 0, r.length);
			}
		}
		
		System.out.println("Writing results.");
		// write all the things.
		TextFile out = new TextFile(outloc, TextFile.W);
		out.writeln(outputHeader);
		int nrWritten = 0;
		for (QTLObject o : r) {
			if (o != null) {
				out.writeln(o.s);
				nrWritten++;
				if (nrWritten % 10000 == 0) {
					System.out.print("\r" + nrWritten + " written.");
				}
			}
		}
		out.close();
		
		System.out.println(nrWritten + " written in total.");
		System.out.println("Done. Result is here: " + outloc);
		
		
	}
	
	private class QTLObject implements Comparable<QTLObject> {
		
		
		private final double p;
		private final double z;
		private final String s;
		
		public QTLObject(double p, double z, String s) {
			this.p = p;
			this.z = z;
			this.s = s;
		}
		
		@Override
		public boolean equals(Object o) {
			if (this == o) return true;
			if (o == null || getClass() != o.getClass()) return false;
			
			QTLObject qtlObject = (QTLObject) o;
			
			if (Double.compare(qtlObject.p, p) != 0) return false;
			if (Double.compare(qtlObject.z, z) != 0) return false;
			return s != null ? s.equals(qtlObject.s) : qtlObject.s == null;
		}
		
		@Override
		public int hashCode() {
			int result;
			long temp;
			temp = Double.doubleToLongBits(p);
			result = (int) (temp ^ (temp >>> 32));
			temp = Double.doubleToLongBits(z);
			result = 31 * result + (int) (temp ^ (temp >>> 32));
			result = 31 * result + (s != null ? s.hashCode() : 0);
			return result;
		}
		
		@Override
		public int compareTo(QTLObject o) {
			if (this.equals(o)) {
				return 0;
			} else if (o == null) {
				return -1;
			} else if (this.p < o.p) {
				return -1;
			} else if (this.p > o.p) {
				return 1;
			} else {
				double thisz = Math.abs(this.z);
				double oz = Math.abs(o.z);
				return Double.compare(thisz, oz);
			}
		}
	}
	
}
