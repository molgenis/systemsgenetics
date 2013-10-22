/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.util;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import umcg.genetica.io.gmt.GMTFile;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix.DoubleMatrixDataset;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class ExpressionFileMerger {

    public void merge(String file1, String file2, String outfile) {
//        try {
//            DoubleMatrixDataset<String, String> dataset1 = new DoubleMatrixDataset<String, String>(file1);
//            DoubleMatrixDataset<String, String> dataset2 = new DoubleMatrixDataset<String, String>(file2);
//
//            String[] probes1 = dataset1.rowObjects.toArray(new String[0]);
//            String[] probes2 = dataset2.rowObjects.toArray(new String[0]);
//            int numSamples = dataset1.colObjects.size() + dataset2.colObjects.size();
//            HashMap<String, Integer> probeToProbeId = new HashMap<String, Integer>();
//            HashSet<String> probesInDataset2 = new HashSet<String>();
//            probesInDataset2.addAll(Arrays.asList(probes2));
//
//            String[] newSampleNames = new String[numSamples];
//
//            String[] dataset1samples = dataset1.colObjects.toArray(new String[0]);
//            String[] dataset2samples = dataset2.colObjects.toArray(new String[0]);
//            System.arraycopy(dataset1samples, 0, newSampleNames, 0, dataset1samples.length);
//            System.arraycopy(dataset2samples, 0, newSampleNames, dataset1samples.length, dataset2samples.length);
//
//            int ctr = 0;
//            int pos = 0;
//            for (String probe : probes1) {
//                if (probesInDataset2.contains(probe)) {
//                    ctr++;
//                    probeToProbeId.put(probe, pos);
//                }
//                pos++;
//            }
//
////            String[] sharedProbes = new String[ctr];
////            double[][] newmatrix = new double[ctr][numSamples];
//
//        } catch (Exception e) {
//        }
    }

    public void collapseOnGeneName(String in, String probeToGeneNameFile, String outfilename) throws IOException {

        HashMap<String, String> probeToGeneMap = null;
        TextFile pm = new TextFile(probeToGeneNameFile, TextFile.R);
        probeToGeneMap = (HashMap<String, String>) pm.readAsHashMap(0, 1);
        pm.close();

        DoubleMatrixDataset<String, String> d = new DoubleMatrixDataset<String, String>(in);
        String[] probeNames = d.rowObjects.toArray(new String[0]);
        String[] sampleNames = d.colObjects.toArray(new String[0]);

        HashSet<String> genesAvailable = new HashSet<String>();
        for (int i = 0; i < probeNames.length; i++) {
            String gene = probeToGeneMap.get(probeNames[i]);
            if (gene == null) {
                gene = "-";
                probeToGeneMap.put(probeNames[i], gene);
            }
            genesAvailable.add(gene);
        }
        System.out.println("Total genes in expression file: " + genesAvailable.size());

        HashMap<String, Integer> geneToColMap = new HashMap<String, Integer>();
        int geneCtr = 0;
        for (String gene : genesAvailable) {
            geneToColMap.put(gene, geneCtr);
            geneCtr++;
        }

        HashMap<String, Integer> nrProbesPerGene = new HashMap<String, Integer>();
        for (int i = 0; i < probeNames.length; i++) {
            String gene = probeToGeneMap.get(probeNames[i]);
            Integer id = geneToColMap.get(gene);
            Integer nrProbes = nrProbesPerGene.get(gene);
            if (nrProbes == null) {
                nrProbes = 0;
            }

            nrProbes++;
            nrProbesPerGene.put(gene, nrProbes);
        }

        double[][] collapsedData = new double[geneCtr][sampleNames.length];
        double[][] rawData = d.rawData;

        for (int sample = 0; sample < sampleNames.length; sample++) {
            for (int probe = 0; probe < probeNames.length; probe++) {
                String probeName = probeNames[probe];
                String gene = probeToGeneMap.get(probeName);
                Integer geneId = geneToColMap.get(gene);
                collapsedData[geneId][sample] += rawData[probe][sample];
            }
        }

        String[] newProbeNames = new String[genesAvailable.size()];
        for (String gene : genesAvailable) {
            Integer id = geneToColMap.get(gene);
            Integer nrProbes = nrProbesPerGene.get(gene);
            for (int i = 0; i < sampleNames.length; i++) {
                collapsedData[id][i] /= nrProbes;
            }
            newProbeNames[id] = gene;
        }



        DoubleMatrixDataset<String, String> dout = new DoubleMatrixDataset<String, String>();
        dout.colObjects = Arrays.asList(sampleNames);
        dout.rowObjects = Arrays.asList(newProbeNames);
        dout.rawData = collapsedData;
        dout.save(outfilename + ".gz");

        TextFile out = new TextFile(outfilename + "-NrProbesPerGene.txt", TextFile.W);
        out.writeln("Gene\tNrProbes");
        for (String gene : genesAvailable) {
            Integer nrProbes = nrProbesPerGene.get(gene);
            out.writeln(gene + "\t" + nrProbes);
        }
        out.close();


    }

    public void collapseProbesBasedOnPathwayAnnotation(String infile, String ensemblannot, String probeannot, String pathwayfile, int annotcol, String pathwayname, boolean standardnormalize) throws IOException {

        TextFile ptf = new TextFile(probeannot, TextFile.R);
        String[] header = ptf.readLineElems(TextFile.tab); // skip header

        System.out.println("Assuming annotation " + header[annotcol]);
        String[] elems = ptf.readLineElems(TextFile.tab);
        HashMap<String, String> probeNrToHT12 = new HashMap<String, String>();
        while (elems != null) {

            if (elems.length > annotcol) {
                String ht12probe = elems[annotcol];

                String probeId = elems[0];
                if (ht12probe.length() > 1) {
                    probeNrToHT12.put(probeId, ht12probe);
                }
            } else {
                System.err.println("WARNING: probe annotation file does not contain all elements expected for line: ");
                System.err.println(Strings.concat(elems, Strings.tab));
            }
            elems = ptf.readLineElems(TextFile.tab);
        }
        ptf.close();

        TextFile etf = new TextFile(ensemblannot, TextFile.R);

        elems = etf.readLineElems(TextFile.tab);
        HashMap<String, String> probeToEns = new HashMap<String, String>();
        while (elems != null) {

            if (elems.length >= 5) {
                String probe = elems[0].trim();
                if (probeNrToHT12.get(probe) != null) {
                    String ens = elems[4].trim();

                    probeToEns.put(probeNrToHT12.get(probe), ens);
                }

            }
            elems = etf.readLineElems(TextFile.tab);
        }
        etf.close();

        // "/Volumes/iSnackHD/PathwayData/2012-06-01-PathwayGMTFiles/reactome/Reactome+.gmt"
        GMTFile gmt = new GMTFile(pathwayfile);

        ArrayList<String> pathways = (ArrayList<String>) gmt.getPathways();
        String query = "";
        HashSet<String> genesInPathway = new HashSet<String>();


        DoubleMatrixDataset<String, String> ds = new DoubleMatrixDataset<String, String>(infile);
        String[] probeNames = ds.rowObjects.toArray(new String[0]);
        ArrayList<String> selectedPathways = new ArrayList<String>();
        for (int pw = 0; pw < pathways.size(); pw++) {

            query = pathways.get(pw);
            if (pathways.contains(query)) {
                genesInPathway = (HashSet<String>) gmt.getGenesForPathway(query);
            }

            int nonannotated = 0;
            int nrcollapsed = 0;

            HashMap<String, Integer> probeToProbeId = new HashMap<String, Integer>();
            ArrayList<Integer> probesToMerge = new ArrayList<Integer>();
            int probectr = 1;
            for (int i = 0; i < probeNames.length; i++) {
                String probe = probeNames[i];
                String ensannot = probeToEns.get(probe);
                if (ensannot == null) {
                    // kick out probe
                    nonannotated++;
                    probeToProbeId.put(probe, null);
                } else {
                    if (genesInPathway.contains(ensannot)) {
                        nrcollapsed++;
                        probesToMerge.add(i);
                        probeToProbeId.put(probe, 0);
                    } else {
                        probeToProbeId.put(probe, probectr);
                        probectr++;
                    }
                }
            }
            if (nrcollapsed > 3) {
                selectedPathways.add(query);
            }
        }

        if (standardnormalize) {
            ds.standardNormalizeData();
        }
        String[] newProbeNames = selectedPathways.toArray(new String[0]);
        double[][] newData = new double[newProbeNames.length][ds.colObjects.size()];

        TextFile nrGenesUsed = new TextFile(infile + "-NrGenesUsedForCollapse-" + pathwayname + ".txt", TextFile.W);
        for (int pw = 0; pw < selectedPathways.size(); pw++) {

            query = selectedPathways.get(pw);
            if (pathways.contains(query)) {
                genesInPathway = (HashSet<String>) gmt.getGenesForPathway(query);
            }

            int nonannotated = 0;
            int nrcollapsed = 0;

            HashMap<String, Integer> probeToProbeId = new HashMap<String, Integer>();
            ArrayList<Integer> probesToMerge = new ArrayList<Integer>();
            int probectr = 1;
            for (int i = 0; i < probeNames.length; i++) {
                String probe = probeNames[i];
                String ensannot = probeToEns.get(probe);
                if (ensannot == null) {
                    // kick out probe
                    nonannotated++;
                    probeToProbeId.put(probe, null);
                } else {
                    if (genesInPathway.contains(ensannot)) {
                        nrcollapsed++;
                        probesToMerge.add(i);
                        probeToProbeId.put(probe, 0);
                    } else {
                        probeToProbeId.put(probe, probectr);
                        probectr++;
                    }
                }
            }

            System.out.println(nonannotated + " probes don't have an ensembl annotation, " + nrcollapsed + " will be collapsed");
            int finalNrOfProbes = probectr;
            System.out.println("Final size: " + probectr);

            nrGenesUsed.writeln(query + "\t" + nonannotated + " probes don't have an ensembl annotation, " + nrcollapsed + " will be collapsed. Final size: " + probectr);
            double[][] originalData = ds.rawData;
            double[] mergedProbes = new double[ds.colObjects.size()];



            for (int col = 0; col < ds.colObjects.size(); col++) {
                double probeSum = 0;
                // double[] vals = new double[probesToMerge.size()];
                int q = 0;
                for (Integer i : probesToMerge) {
                    double v = originalData[i][col];
                    probeSum += v;
                    //  vals[q] = v;
                    q++;
                }
                // double sd = JSci.maths.ArrayMath.standardDeviation(vals);
                mergedProbes[col] = probeSum / probesToMerge.size();
                //mergedProbes[col] /= sd;
            }
            newData[pw] = mergedProbes;
        }
        nrGenesUsed.close();

        String outfile = infile + "-PathWayCollapsed-" + pathwayname;
        if (standardnormalize) {
            outfile += "-SDNorm";
        }

        outfile += ".txt.gz";

        System.out.println("Outfile: " + outfile);
        ds.rawData = (newData);
        ds.rowObjects = Arrays.asList(newProbeNames);
        ds.save(outfile);

    }
}
