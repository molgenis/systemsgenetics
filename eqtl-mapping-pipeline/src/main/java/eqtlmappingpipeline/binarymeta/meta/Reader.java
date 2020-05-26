/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.binarymeta.meta;

import umcg.genetica.io.trityper.bin.BinaryResultDataset;
import umcg.genetica.io.trityper.bin.BinaryResultProbe;
import umcg.genetica.io.trityper.bin.BinaryResultSNP;
import umcg.genetica.math.stats.Descriptives;
import umcg.genetica.io.trityper.util.BaseAnnot;

/**
 *
 * @author harmjan
 */
public class Reader {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        Reader r =new Reader();
        r.run("Dataset","/Data/eQTLTest/Meta3-bin-GRNGDataOnly-SS/");
    }

    public void run(String datasetname, String location) {
        try{
            BinaryResultDataset ds = new BinaryResultDataset(location, datasetname,0);
            BinaryResultProbe[] dsProbes = ds.getProbes();
            BinaryResultSNP[] dsSNPs = ds.getSnps();
            int nrTotalSamples = ds.getMaxNrSamples();
            Descriptives.lookupSqrt(nrTotalSamples);
            int snpId = 0;

            BinaryResultSNP snpObject = ds.getSnps()[snpId];
            long pointer = snpObject.getzScoreIndex();
            long nextpointer = -1;

            if (snpId+1 < ds.getSnps().length) {
                BinaryResultSNP snpObject2 = ds.getSnps()[snpId+1];
                nextpointer = snpObject2.getzScoreIndex();
            }


            Float[] zscores = ds.getMatrix().read(pointer, nextpointer, ds.getNumProbes());
            System.out.println("Assessed allele: "+ BaseAnnot.toString(snpObject.getAssessedAllele()));
            for(int p=0; p<dsProbes.length; p++){
                int nrSamples           = snpObject.getNumsamples();
                double weight           = Descriptives.getSqrt(nrSamples);
                double zscore           = zscores[p];
                double zSum = (zscore * weight);
                double finalzScore = zSum / weight;

                System.out.println(p+"\t"+dsProbes[p].getName()+"\t"+zscores[p]+"\t"+finalzScore);
            }


        } catch (Exception e){
            e.printStackTrace();
        }

    }
}
