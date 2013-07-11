/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.io.text.converters;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author juha
 */
public class NetworkConverter {

    private static final Logger LOGGER = Logger.getLogger(NetworkConverter.class.getName());

    // prevent instantiation, only static factory methods
    private NetworkConverter() {
    }

    public static void convertPlainTextDoubleMatrixToPlainTextEdgeList(String matrixfile, String outfile, String infileDelimiter, String outfileDelimiter) throws IOException {

        TextFile in = new TextFile(matrixfile, TextFile.R);
        String line = in.readLine();
        String[] split = line.split(infileDelimiter);
        String[] nodes = Arrays.copyOfRange(split, 1, split.length);
        LOGGER.log(Level.FINE, "{0} nodes in {1}", new Object[]{nodes.length, matrixfile});

        TextFile out = new TextFile(outfile, TextFile.W);
        int lineNr = 1;
        int nrEdges = 0;
        while ((line = in.readLine()) != null) {
            split = line.split(infileDelimiter);
            if (split.length != nodes.length + 1) {
                throw new IllegalArgumentException("The data in file '" + matrixfile + "' are not a matrix. Check line " + (lineNr + 1) + ".");
            }
            String node = split[0];
            for (int i = 1; i < split.length; i++) {
                if (lineNr == i) {
                    continue; // skip diagonal
                }
                try {
                    double weight = Double.parseDouble(split[i]);
                    if (weight > 0) {
                        out.writeln(node + outfileDelimiter + nodes[i - 1] + outfileDelimiter + weight);
                        nrEdges++;
                    }
                } catch (NumberFormatException ex) {
                    throw new IllegalArgumentException("The data in file '" + matrixfile + "' are not numerical (cast to double failed). Check line " + (lineNr + 1) + ", column " + (i + 1) + ".");
                }
            }
            lineNr++;
        }
        in.close();
        out.close();
        LOGGER.log(Level.FINE, "{0} edges from {1} written to ''{2}''", new Object[]{nrEdges, matrixfile, outfile});
    }

    public static void writeGMTFileBasedOnGeneSetFileAndMappingFileRemovingDuplicateGeneSets(String genesetfile, String mappingfile, String mappingdelimiter, String gmtfile) throws IOException {
        TextFile in = new TextFile(mappingfile, TextFile.R);
        String line = in.readLine();
        Map<String, String> code2name = new HashMap<String, String>();
        while ((line = in.readLine()) != null) {
            String[] split = line.split(mappingdelimiter);
            String name = split[2].trim().replace("\"", "") + " (" + split[6].trim() + ")";
            code2name.put(split[0].trim().replace("\"", ""), name);
        }
        in.close();
        System.out.println(code2name.size() + " gene set annotations read");

        in = new TextFile(genesetfile, TextFile.R);
        TextFile out = new TextFile(gmtfile, TextFile.W);
        Map<String, Integer> usedNames = new HashMap<String, Integer>();
        Set<String> usedGeneSets = new HashSet<String>();
        while ((line = in.readLine()) != null) {
            String[] split = line.split("\t");
            if (split.length == 0) {
                continue;
            }
            String code = split[0].trim();
            String name = code2name.get(code);
            if (name == null) {
                LOGGER.log(Level.WARNING, "No annotation for gene set ''{0}''", code);
            } else {
                String genes = Arrays.asList(Arrays.copyOfRange(split, 1, split.length)).toString();
                System.out.println(genes);
                if (!usedGeneSets.contains(genes)) {
                    Integer oldNrItems = usedNames.get(name);
                    if (oldNrItems == null) {
                        out.write(code + "\t" + name);
                        for (int i = 1; i < split.length; i++) {
                            out.write("\t" + split[i]);
                        }
                        out.writeln();
                        usedNames.put(name, split.length - 1);
                        usedGeneSets.add(genes);
                    } else {
//                    if (oldNrItems != split.length - 1) {
                        System.out.println("Sets with different numbers of genes for " + name);
//                    }
                    }
                }
            }
        }
        in.close();
        out.close();
    }

    public static void convertGMTFileToPlainTextBinaryNetwork(String infile, String outfile, boolean hasIds, boolean isActuallyGMT) throws IOException {
        TextFile in = new TextFile(infile, TextFile.R);
        String line = null;
        Map<String, Integer> hashItems = new HashMap<String, Integer>();
        Map<String, Set<Integer>> hashSetIndices = new HashMap<String, Set<Integer>>();
        List<String> sets = new ArrayList<String>();
        int nextIndex = 0;
        int nrItemSets = 0;
        while ((line = in.readLine()) != null) {
            String[] split = line.split("\t");
            if (split.length == 0) {
                continue;
            }
            int firstItemIndex = hasIds ? (isActuallyGMT ? 2 : 1) : 0;
            if (hasIds) {
                sets.add(split[0]);
            }
            for (int i = firstItemIndex; i < split.length; i++) {
                Integer itemIndex = hashItems.get(split[i]);
                if (itemIndex == null) {
                    hashItems.put(split[i], nextIndex);
                    nextIndex++;
                }
                Set<Integer> setIndicesThisItem = hashSetIndices.get(split[i]);
                if (setIndicesThisItem == null) {
                    setIndicesThisItem = new HashSet<Integer>();
                    hashSetIndices.put(split[i], setIndicesThisItem);
                }
                setIndicesThisItem.add(nrItemSets);
            }
            nrItemSets++;
        }
        in.close();

        TextFile out = new TextFile(outfile, TextFile.W);
        for (int i = 0; i < nrItemSets; i++) {
            if (hasIds) {
                out.write("\t" + sets.get(i));
            } else {
                out.write("\tComplex" + (i + 1));
            }
        }
        out.writeln();
        for (String item : hashSetIndices.keySet()) {
            if (!"".equals(item)) {
                out.write(item);
                Set<Integer> setIndicesThisItem = hashSetIndices.get(item);
                for (int set = 0; set < nrItemSets; set++) {
                    if (setIndicesThisItem.contains(set)) {
                        out.write("\t1");
                    } else {
                        out.write("\t0");
                    }
                }
                out.writeln();
            }
        }
        out.close();
    }
}
