/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.math.stats.concurrent;

import java.util.concurrent.CompletionService;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.containers.Pair;
import umcg.genetica.math.stats.Descriptives;

/**
 *
 * @author harmjan
 */
public class ConcurrentCorrelation {

    private int nrThreads = Runtime.getRuntime().availableProcessors();

    public ConcurrentCorrelation() {
    }

    public ConcurrentCorrelation(int nrProcs) {
        nrThreads = nrProcs;
    }

    public double[][] pairwiseCorrelation(double[][] in) {
        ExecutorService threadPool = Executors.newFixedThreadPool(nrThreads);
        CompletionService<Pair<Integer, double[]>> pool = new ExecutorCompletionService<Pair<Integer, double[]>>(threadPool);
        double meanOfSamples[] = new double[in.length];
        
        for(int i=0; i<meanOfSamples.length; ++i){
            meanOfSamples[i] = Descriptives.mean(in[i]);
        }
        
        for (int i = 0; i < in.length; i++) {
            ConcurrentCorrelationTask task = new ConcurrentCorrelationTask(in, meanOfSamples, i);
            pool.submit(task);
        }

        int returned = 0;

        double[][] correlationMatrix = new double[in.length][0];
        ProgressBar pb = new ProgressBar(in.length, "Calculation of correlation matrix: " + in.length + " x " + in.length);
        while (returned < in.length) {
            try {
                Pair<Integer, double[]> result = pool.take().get();
                if (result != null) {
                    int rownr = result.getLeft(); //  < 0 when row is not to be included because of hashProbesToInclude.
                    if (rownr >= 0) {
                        double[] doubles = result.getRight();
                        correlationMatrix[rownr] = doubles;
                    }
                    result = null;
                    returned++;
                    pb.iterate();
                }
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
        for(int r=1;r<correlationMatrix.length; r++){
            for(int c=0; c<r; c++){
                correlationMatrix[r][c] = correlationMatrix[c][r];
            }
        }
        threadPool.shutdown();
        pb.close();
        return correlationMatrix;
    }
}
