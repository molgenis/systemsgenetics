package eqtlmappingpipeline.util;

import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;

public class QTLFileMerger {

	public void mergeChr(String loc, String out) throws IOException {

		// try comma sep list first;
		String[] locelems = loc.split(",");
		if (locelems.length > 1) {
			merge(locelems, out);
		} else {
			ArrayList<String> locs = new ArrayList<>();
			for (int chr = 1; chr < 23; chr++) {
				String f = loc.replaceAll("CHR", "" + chr);
				if (Gpio.exists(f)) {
					locs.add(f);
					System.out.println("Found: " + f);
				}
			}
			if (locs.size() == 1) {
				System.out.println("Only one QTL file specified. Use comma separated list, or CHR as template.");
			} else {
				merge(locs.toArray(new String[0]), out);
			}

		}


	}

	public void merge(String[] listOfQTLFiles, String outloc) throws IOException {

		String outputHeader = null;
		int qtlctr = 0;
		double highestP = 2;

		TextFile outf = new TextFile(outloc, TextFile.W);


		for (String s : listOfQTLFiles) {
			System.out.println("Parsing: " + s);
			TextFile tf = new TextFile(s, TextFile.R);
			String fileHeader = tf.readLine();
			if (outputHeader == null) {
				outputHeader = fileHeader;
				outf.writeln(outputHeader);
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
			while (ln != null) {
				outf.writeln(ln);
				nrinfile++;
				if (nrinfile % 10000 == 0) {
					System.out.print("\r" + nrinfile + " parsed. Total QTL: " + qtlctr);
				}
				qtlctr++;
				ln = tf.readLine();
			}
			tf.close();
			System.out.println("Done");
			System.out.println(nrinfile + " QTL in file. Total QTL: " + qtlctr);
		}
		System.out.println();
		outf.close();


		System.out.println(qtlctr + " written in total.");
		System.out.println("Now sorting results");
		QTLFileSorter sorter = new QTLFileSorter();
		String tmpfile = outloc + "_tmp-sort.txt.gz";
		sorter.run(outloc, tmpfile, QTLFileSorter.SORTBY.Z);
		Gpio.moveFile(tmpfile, outloc);

		System.out.println("Done. Result is here: " + outloc);


	}


}
