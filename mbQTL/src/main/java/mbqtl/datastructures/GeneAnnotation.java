package mbqtl.datastructures;

import umcg.genetica.enums.Chromosome;
import umcg.genetica.enums.Strand;
import umcg.genetica.features.Gene;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Set;

public class GeneAnnotation {

	ArrayList<String> geneIdList = new ArrayList<>();
	HashMap<String, Gene> geneIdToObjMap = new HashMap<>();
	HashMap<String, Integer> geneIndex = new HashMap<>();

	public GeneAnnotation(String geneAnnotationFile, int chromosome, Set<String> geneLimitSet) throws IOException {

		load(geneAnnotationFile, chromosome, geneLimitSet);

	}

	private void load(String geneAnnotationFile, int limitChromosome, Set<String> geneLimitSet) throws IOException {
		System.out.println("Parsing " + geneAnnotationFile);
		TextFile tf = new TextFile(geneAnnotationFile, TextFile.R);
		String[] header = tf.readLineElems(TextFile.tab);
		int idcol = -1;
		int chrcol = -1;
		int stacol = -1;
		int stocol = -1;
		int strandcol = -1;
		int symbolcol = -1;

		int maxCol = 0;
		for (int i = 0; i < header.length; i++) {
			String h = header[i].toLowerCase();
			if (h.equals("ensembl") || h.equals("gene") || h.equals("id") || h.equals("gene_id") || h.equals("arrayaddress")) {
				idcol = i;
				if (i > maxCol) {
					maxCol = i;
				}
			} else if (h.equals("chr") || h.equals("chromosome")) {
				chrcol = i;
				if (i > maxCol) {
					maxCol = i;
				}
			} else if (h.equals("chrstart") || h.equals("start") || h.equals("start_position") || h.equals("startposition")) {
				stacol = i;
				if (i > maxCol) {
					maxCol = i;
				}
			} else if (h.equals("chrend") || h.equals("stop") || h.equals("stop_position") || h.equals("stopposition") || h.equals("end") || h.equals("end_position") || h.equals("endposition")) {
				stocol = i;
				if (i > maxCol) {
					maxCol = i;
				}
			} else if (h.equals("symbol") || h.equals("genesymbol") || h.equals("gene_symbol") || h.equals("hgnc_symbol") || h.equals("hgnc")) {
				symbolcol = i;
				if (i > maxCol) {
					maxCol = i;
				}
			} else if (h.equals("strand")) {
				strandcol = i;
				if (i > maxCol) {
					maxCol = i;
				}
			}
		}

		boolean ok = true;
		if (idcol == -1) {
			System.out.println("ID col not found (ensembl|gene|id|gene_id|arrayaddress) in " + geneAnnotationFile);
			ok = false;
		} else {
			System.out.println("ID column: " + idcol);
		}
		if (chrcol == -1) {
			System.out.println("Chromosome col not found (chr|chromosome) in " + geneAnnotationFile);
			ok = false;
		} else {
			System.out.println("Chromosome column: " + chrcol);
		}
		if (stacol == -1) {
			System.out.println("Chromosome start position col not found (chrstart|start|start_position|startposition) in " + geneAnnotationFile);
			ok = false;
		} else {
			System.out.println("Chromosome start position column: " + stacol);
		}
		if (stocol == -1) {
			System.out.println("Chromosome start position col not found (chrend|stop|stop_position|stopposition) in " + geneAnnotationFile);
			ok = false;
		} else {
			System.out.println("Chromosome stop position column: " + stocol);
		}
		System.out.println("Strand column: " + (strandcol == -1 ? "Not found" : strandcol));
		System.out.println("Gene symbol column: " + (symbolcol == -1 ? "Not found" : symbolcol));
		System.out.println("");

		if (!ok) {
			throw new RuntimeException("Required columns not found in: " + geneAnnotationFile);
		}

		String[] elems = tf.readLineElems(TextFile.tab);
		int gctr = 0;
		int lctr = 1;
		while (elems != null) {
			if (elems.length > maxCol) {


				Chromosome chr = Chromosome.parseChr(elems[chrcol]);
				boolean include = false;
				if (limitChromosome == -1 || chr.getNumber() == limitChromosome) {
					include = true;
				}

				String id = elems[idcol];
				if (geneLimitSet != null && !geneLimitSet.contains(id)) {
					include = false;
				}

				if (include) {
					int start = Integer.parseInt(elems[stacol]);
					int stop = Integer.parseInt(elems[stocol]);
					String symbol = Strings.cache("NA");
					if (symbolcol > -1) {
						symbol = elems[symbolcol];
					}
					Strand strand = Strand.NA;
					if (strandcol > -1) {
						strand = Strand.parseStr(elems[strandcol]);
					}


					Gene g = new Gene(id, chr, strand, start, stop);
					g.setGeneSymbol(symbol);
					Gene g2 = geneIdToObjMap.get(id);
					if (g2 == null) {
						geneIdToObjMap.put(id, g);
						geneIndex.put(id, gctr);
						geneIdList.add(id);
						gctr++;
					} else {
						boolean equal = g.equals(g2);
						System.out.println(geneAnnotationFile + " has multiple annotations for " + id + "; current (line: " + lctr + ") is equal to previous occurence: " + equal);
						System.out.println("Note: loaded annotation for " + id + " will reflect first occurrence in " + geneAnnotationFile);
					}
				}


			}
			elems = tf.readLineElems(TextFile.tab);
			lctr++;
		}
		tf.close();
		System.out.println("Annotations loaded for " + geneIdList.size() + " genes/phenotypes");
	}


	public ArrayList<String> getAllGenes() {
		return geneIdList;
	}

	public Gene getGene(String id) {
		return geneIdToObjMap.get(id);
	}
}
