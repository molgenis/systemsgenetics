/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.util;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map.Entry;
import umcg.genetica.io.text.TextFile;
import cern.jet.stat.tdouble.Probability;

/**
 *
 * @author MarcJan
 */
public class CheckAssocationMultipleGeneProbes {
    public static void main(String[] args) {
        
        HashMap<String, ArrayList<String>> hashGenes = new HashMap<String, ArrayList<String>>();
        LinkedHashMap<String, Double> hashTValues = new LinkedHashMap<String, Double>();
        HashMap<String, Integer> hashCpGSites = new HashMap<String, Integer>();

        try {
            TextFile in = new TextFile("C:\\Users\\MarcJan\\Desktop\\pcas_450K.txt", TextFile.R);
            String str = in.readLine();

            int i=-1;
            while ((str = in.readLine()) != null) {
                String[] data = str.split("\t");
                String CpG = data[0];
//                if (data.length>7 && !data[7].equals("")) {
                if (data.length>6) {
                    ++i;
                    String gene = data[6];
//                    double t = Double.parseDouble(data[1]);
//                    double t = Double.parseDouble(data[2]);
//                    double t = Double.parseDouble(data[3]);
//                    double t = Double.parseDouble(data[4]);
                    double t = Double.parseDouble(data[5]);
                    hashTValues.put(CpG, t);
                    hashCpGSites.put(CpG, i);
                    
                    if (!hashGenes.containsKey(gene)) {
                        ArrayList<String> CpGs = new ArrayList<String>();
                        CpGs.add(CpG);
                       
                        hashGenes.put(gene, CpGs);
                    } else {
                        hashGenes.get(gene).add(CpG);
                    }
                }
               
            }
            in.close();
//            System.out.println(hashGenes.size());
//            System.out.println(hashTValues.size());
        } catch (Exception e) {
            System.out.println("Error:\t" + e.getMessage());
            e.printStackTrace();
        }
        
        double[] pValues = new double[hashTValues.size()];
//        System.out.println(pValues.length);
        int c=0;
        for (Entry<String, Double> ent : hashTValues.entrySet()) {
            double t = ent.getValue();
//            if (t<0) {
//                pValues[c] = t;
//            } else {
//                pValues[c] = 0;
//            }
//            pValues[c] = Math.abs(t);
            pValues[c] = t;
            c++;
            
        }
        
        RankArray r = new RankArray();
        pValues = r.rank(pValues, false);
        
//        double[] rank = r.rank(pValues, false);
//        for (int i=0; i<pValues.length; i++) {
//            System.out.println(rank[i]);
//        }

        for (int i=0; i<pValues.length; i++) {
            pValues[i] = (pValues[i] - 0.5d) / (double) pValues.length;
        }

        JSci.maths.statistics.NormalDistribution normDist = new JSci.maths.statistics.NormalDistribution();
        for (Entry<String, ArrayList<String>> ent : hashGenes.entrySet()) {
            String gene = ent.getKey();
            ArrayList<String> sites = ent.getValue();
            double[] pSites = new double[sites.size()];
            double chiSquare = 0;
            double[] zSites = new double[sites.size()];
            for (int p=0; p<sites.size(); p++) {
                int cpgIndex = hashCpGSites.get(sites.get(p));
                pSites[p] = pValues[cpgIndex];
                chiSquare+=-2*Math.log(pSites[p]);
                zSites[p] = Probability.normalInverse(pSites[p]);
            }
            //JSci.maths.statistics.ChiSqrDistribution dist = new JSci.maths.statistics.ChiSqrDistribution(2*sites.size());
            //double chiSquareP = 1 - dist.cumulative(chiSquare);
            double zOverall = JSci.maths.ArrayMath.mass(zSites) / Math.sqrt(sites.size());
            
            
            double minPValue = JSci.maths.ArrayMath.min(pSites);
            //double minPValueCorrected = 1d - Math.pow(1d - minPValue, (double) sites.size());
//            System.out.println(gene + "\t" + sites.size() + "\t" + minPValue + "\t" + minPValueCorrected + "\t" + chiSquare + "\t" + chiSquareP);
            System.out.println(gene + "\t" + sites.size() + "\t" + minPValue + "\t" + zOverall);
            
        }
    }

}
