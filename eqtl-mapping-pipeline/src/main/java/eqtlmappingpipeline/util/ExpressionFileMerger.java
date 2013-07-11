/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.util;

import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import umcg.genetica.math.matrix.DoubleMatrixDataset;

/**
 *
 * @author harmjan
 */
public class ExpressionFileMerger {
    public void merge(String file1, String file2, String outfile){
        try{
            DoubleMatrixDataset<String,String> dataset1 = new DoubleMatrixDataset<String,String>(file1);
            DoubleMatrixDataset<String,String> dataset2 = new DoubleMatrixDataset<String,String>(file2);
            
            String[] probes1 = dataset1.rowObjects.toArray(new String[0]);
            String[] probes2 = dataset2.rowObjects.toArray(new String[0]);
            int numSamples = dataset1.colObjects.size()+dataset2.colObjects.size();
            HashMap<String, Integer> probeToProbeId = new HashMap<String, Integer>();
            HashSet<String> probesInDataset2 = new HashSet<String>();
            probesInDataset2.addAll(Arrays.asList(probes2));
            
            String[] newSampleNames = new String[numSamples];
            
            String[] dataset1samples = dataset1.colObjects.toArray(new String[0]);
            String[] dataset2samples = dataset2.colObjects.toArray(new String[0]);
            System.arraycopy(dataset1samples, 0, newSampleNames, 0, dataset1samples.length);
            System.arraycopy(dataset2samples, 0, newSampleNames, dataset1samples.length, dataset2samples.length);
            
            int ctr = 0;
            int pos = 0;
            for(String probe: probes1){
                if(probesInDataset2.contains(probe)){
                    ctr++;
                    probeToProbeId.put(probe, pos);
                }
                pos++;
            }
            
            String[] sharedProbes = new String[ctr];
            double[][] newmatrix  = new double[ctr][numSamples];
            
        } catch (Exception e){
            
        }
    }
}
