package mbqtl.datastructures;

import umcg.genetica.io.text.TextFile;

import org.apache.commons.lang3.StringUtils;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;

public class GeneExpressionData {
	public HashMap<String, Integer> sampleMap;
	public String[] samples;
	public HashMap<String, Integer> geneMap;
	public String[] genes;
	public double[][] data;

	public GeneExpressionData(String geneExpressionDataFile, ArrayList<String> allGenes, HashSet<String> rnaSamplesMatchedToDNA) throws IOException {

		load(geneExpressionDataFile, allGenes, rnaSamplesMatchedToDNA);


	}

	private void load(String geneExpressionDataFile, ArrayList<String> allGenes, HashSet<String> rnaSamplesMatchedToDNA) throws IOException {
		System.out.println("Loading phenotype matrix: " + geneExpressionDataFile);
		HashSet<String> geneSet = null;
		if (allGenes != null) {
			geneSet = new HashSet<>();
			geneSet.addAll(allGenes);
		}

		TextFile tf = new TextFile(geneExpressionDataFile, TextFile.R);
		String[] header = tf.readLineElems(TextFile.tab);

		boolean[] includeColumn = new boolean[header.length];
		ArrayList<String> samplesArr = new ArrayList<>();


		int nrColsToInclude = 0;
		for (int v = 0; v < header.length; v++) {
			String sample = header[v];
			if (rnaSamplesMatchedToDNA == null || rnaSamplesMatchedToDNA.contains(sample)) {
				if (sampleMap.containsKey(sample)) {
					System.out.println("Duplicate column: " + v + " in header with ID:" + sample);
				} else {
					samplesArr.add(header[v]);
					sampleMap.put(header[v], nrColsToInclude);
					includeColumn[v] = true;
					nrColsToInclude++;
				}
			}
		}

		System.out.println(samplesArr.size() + " selected columns (samples) in header.");
		samples = samplesArr.toArray(new String[0]);

		String line = tf.readLine();
		int lineCounter = 1;
		boolean warningprinted = false;
		ArrayList<double[]> dataTMP = new ArrayList<>();
		ArrayList<String> geneTMP = new ArrayList<>();
		int geneCounter = 0;

		while (line != null) {
			String[] elems = StringUtils.split(line, "\t", 3); // only fully split the line if we actually need the data
			if (elems.length < 2) {
				System.out.println("Skipping line " + lineCounter + ": fewer than 2 tab separated columns.");
			} else {
				String gene = elems[0];
				if (geneSet == null || geneSet.contains(gene)) {
					if (geneMap.containsKey(gene)) {
						System.out.println("Duplicate phenotype found on line " + lineCounter + " with ID: " + gene);
					} else {
						elems = StringUtils.split(line, "\t");
						int columnOffset = 0;

						// check whether there is as much information as in the header for this line
						if (elems.length != header.length) {
							if (elems.length == header.length + 1) {
								// this is fine; probably means that the header didn't have an ID for the first column.
								// this is often the case if the file comes from R for instance
								columnOffset = 1;
								if (!warningprinted) {
									System.out.println();
									System.out.println("Warning: the header of " + geneExpressionDataFile + " has a different number of tab separated elements compared to line " + lineCounter + ". " + "The header is 1 element shorter. It's probably best to fix your file. Nevertheless, the program will try to continue...");
									System.out.println();
									warningprinted = true;
								}
							} else {
								throw new RuntimeException("Error parsing: " + geneExpressionDataFile + ": number of elements in header (" + header.length + ")does not match length of line: " + lineCounter + " (" + elems.length + ")");
							}
						}

						double[] rowdata = new double[nrColsToInclude];
						Arrays.fill(rowdata, Double.NaN); // initialize with NaN
						int ctr = 0;
						for (int v = 1; v < elems.length; v++) {
							if (includeColumn[v - columnOffset]) {
								try {
									rowdata[ctr] = Double.parseDouble(elems[v]);
								} catch (NumberFormatException nfe) {
									// if non-parse able, keep as NaN
								}
								ctr++;
							}
						}

						dataTMP.add(rowdata);
						geneMap.put(gene, geneCounter);
						geneTMP.add(gene);
						geneCounter++;
					}
				}
			}
			line = tf.readLine();
			if (lineCounter % 10000 == 0) {
				System.out.print(lineCounter + " lines read, " + geneTMP.size() + " phenotypes loaded.\r");
			}
			lineCounter++;
		}
		tf.close();
		System.out.print(lineCounter + " lines read, " + geneTMP.size() + " phenotypes loaded. Done reading.\n");

		// copy the data into an array of arrays (might be faster than an arraylist)
		data = new double[dataTMP.size()][];
		for (int i = 0; i < dataTMP.size(); i++) {
			data[i] = dataTMP.get(i);
		}
		genes = geneTMP.toArray(new String[0]);

	}

	public void save(String s) {

	}
}
