/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.binarymeta.meta;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.zip.DataFormatException;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.io.trityper.bin.BinaryResultDataset;
import umcg.genetica.io.trityper.bin.BinaryResultSNP;

/**
 *
 * @author harmjan
 */
public class MetaAnalysisLoaderThread extends Thread {

    private final LinkedBlockingQueue<MetaAnalysisWorkPackage> m_queue_output;
    private ArrayList<String> snps;
    protected Integer[][] snpTranslation;
    protected BinaryResultDataset[] ds;

    public MetaAnalysisLoaderThread(LinkedBlockingQueue<MetaAnalysisWorkPackage> output,
            Integer[][] snpTranslation,
            ArrayList<String> snps,
            BinaryResultDataset[] ds) {

        this.snps = snps;
        this.ds = ds;
        this.snpTranslation = snpTranslation;
        m_queue_output = output;
    }

    @Override
    public void run() {

        int buffersize = 8;
        int snpsprocessed = 0;

        ProgressBar pb = new ProgressBar(snps.size());
        while (snpsprocessed < snps.size()) {
//	while (snpsprocessed < 1000) {
            if (snpsprocessed + buffersize > snps.size()) {
                buffersize = snps.size() - snpsprocessed;
            }

            MetaAnalysisWorkPackage[] buffer = new MetaAnalysisWorkPackage[buffersize];

            // load the next buffersize SNPs.
            int bufferpos = 0;
            for (int s = snpsprocessed; s < snpsprocessed + buffersize; s++) {
                buffer[bufferpos] = new MetaAnalysisWorkPackage(s, ds.length);
                for (int d = 0; d < ds.length; d++) {
                    Integer snpid = snpTranslation[d][s];
                    buffer[bufferpos].setSNPId(d, snpid);
                }

                bufferpos++;
                // Float[] zscores = ds[d].getMatrix().read(pointer, nextpointer, ds[d].getNumProbes());
            }

            // initialized wp, now sort and load
            bufferpos = 0;
            for (int d = 0; d < ds.length; d++) {
                for (int b = 0; b < buffersize; b++) {
                    buffer[b].setSortByDataset(d);
                }

                Arrays.sort(buffer);

                for (int b = 0; b < buffersize; b++) {
                    Integer snpidtoload = buffer[b].getSNPId(d);
                    if (snpidtoload != null) {
                        BinaryResultSNP snpObject = ds[d].getSnps()[snpidtoload];
                        buffer[b].setSNPObject(d, snpObject);
                        long pointer = snpObject.getzScoreIndex();
                        long nextpointer = -1;

                        if (snpidtoload + 1 < ds[d].getSnps().length) {
                            BinaryResultSNP snpObject2 = ds[d].getSnps()[snpidtoload + 1];
                            nextpointer = snpObject2.getzScoreIndex();
                        }
                        try {
                            byte[] datasetdata = ds[d].getMatrix().readDeflated(pointer, nextpointer, ds[d].getNumProbes());
                            buffer[b].setData(d, datasetdata);
                        } catch (IOException e) {
                            e.printStackTrace();
                        }
                    }
                }
            }

            for (int b = 0; b < buffersize; b++) {
                // add to queue
                try {

                    m_queue_output.put(buffer[b]);
                } catch (InterruptedException ex) {
                    ex.printStackTrace();
                }
                pb.set(snpsprocessed + b);

            }

            snpsprocessed += buffersize;

        }

        pb.close();
        System.out.println("done loading... ");



//	SNP snpObject = pack.getSNPObject(d); // ds[d].getSnps()[snpId];

//		long pointer = snpObject.getzScoreIndex();
//		long nextpointer = -1;
//
//		if (snpId + 1 < ds[d].getSnps().length) {
//		    SNP snpObject2 = ds[d].getSnps()[snpId + 1];
//		    nextpointer = snpObject2.getzScoreIndex();
//		}
    }
}
