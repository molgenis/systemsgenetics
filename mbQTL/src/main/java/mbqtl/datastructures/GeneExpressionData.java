//
// Source code recreated from a .class file by IntelliJ IDEA
// (powered by FernFlower decompiler)
//

package mbqtl.datastructures;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Set;
import mbqtl.Util;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

public class GeneExpressionData {
	public double[][] data;
	public String[] genes;
	public String[] samples;
	public HashMap<String, Integer> sampleMap;
	public HashMap<String, Integer> geneMap;

	public GeneExpressionData(String geneExpressionDataFile, Set<String> geneSelection, Set<String> requestedSamples) throws IOException {
		System.out.println("Loading expression data from: " + geneExpressionDataFile);
		if (requestedSamples != null) {
			System.out.println("Max number of samples: " + requestedSamples.size());
		}

		if (geneSelection != null) {
			System.out.println("Max number of genes: " + geneSelection.size());
		}

		int genecol = 0;
		TextFile tf = new TextFile(geneExpressionDataFile, false);
		String ln = tf.readLine();
		if (ln.startsWith("#Chr")) {
			genecol = 3;
			System.out.println("FastQTL style expression file");
		}

		String[] header = Strings.whitespace.split(ln);
		boolean[] includeColumn = new boolean[header.length];
		ArrayList<String> sampleTmp = new ArrayList();

		for(int i = 1; i < header.length; ++i) {
			String sample = header[i];
			if (requestedSamples == null || requestedSamples.contains(sample)) {
				includeColumn[i] = true;
				sampleTmp.add(sample);
			}
		}

		if (sampleTmp.isEmpty()) {
			System.out.println("No matching samples found in expression data!");
			System.exit(-1);
		}

		System.out.println(sampleTmp.size() + " samples found.");
		this.samples = (String[])sampleTmp.toArray(new String[0]);
		this.sampleMap = Util.hash(sampleTmp);
		String[] elems = tf.readLineElems(Strings.whitespace);
		ArrayList<double[]> dataList = new ArrayList();
		ArrayList<String> genetmp = new ArrayList();
		int lctr = 0;

		boolean printWarning;
		for(printWarning = false; elems != null; elems = tf.readLineElems(Strings.whitespace)) {
			String gene = Strings.cache(elems[genecol]);
			if (geneSelection == null || geneSelection.contains(gene)) {
				double[] dataln = new double[this.samples.length];
				Arrays.fill(dataln, Double.NaN);
				genetmp.add(gene);
				int sctr = 0;

				for(int i = 1; i < elems.length; ++i) {
					if (elems.length < dataln.length) {
						printWarning = true;
					}

					if (includeColumn[i]) {
						try {
							double d = Double.parseDouble(elems[i]);
							dataln[sctr] = d;
						} catch (NumberFormatException var21) {
							dataln[sctr] = Double.NaN;
						}

						++sctr;
					}
				}

				dataList.add(dataln);
			}

			++lctr;
			if (geneSelection != null && dataList.size() == geneSelection.size()) {
				System.out.print(lctr + " lines parsed, " + dataList.size() + " genes loaded.\r");
				break;
			}

			if (lctr % 2000 == 0) {
				System.out.print(lctr + " lines parsed, " + dataList.size() + " genes loaded.\r");
			}
		}

		tf.close();
		System.out.println(lctr + " lines parsed, " + dataList.size() + " genes loaded.");
		if (printWarning) {
			System.err.println("WARNING: some lines in the file " + geneExpressionDataFile + " have fewer elements than specified in the header. Please check your input for broken lines. Missing elements have been set as missing.");
		}

		this.genes = (String[])genetmp.toArray(new String[0]);
		this.data = new double[this.genes.length][0];

		for(int g = 0; g < this.genes.length; ++g) {
			this.data[g] = (double[])dataList.get(g);
		}

		this.geneMap = Util.hash(genetmp);
	}

	public void save(String s) throws IOException {
		TextFile tf = new TextFile(s, true);
		tf.writeln("-\t" + Strings.concat(this.samples, Strings.tab));

		for(int q = 0; q < this.data.length; ++q) {
			tf.writeln(this.genes[q] + "\t" + Strings.concat(this.data[q], Strings.tab));
		}

		tf.close();
	}
}
