/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.math.stats.concurrent;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.tdouble.impl.DenseLargeDoubleMatrix2D;

import java.util.concurrent.CompletionService;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.stream.IntStream;

import umcg.genetica.console.ProgressBar;
import umcg.genetica.containers.Pair;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.stats.Correlation;
import umcg.genetica.math.stats.Descriptives;

/**
 * @author harmjan
 */
public class ConcurrentCorrelation {

    private int nrThreads = Runtime.getRuntime().availableProcessors();

    public ConcurrentCorrelation() {
    }

    public ConcurrentCorrelation(int nrProcs) {
        nrThreads = nrProcs;
    }

    public DoubleMatrixDataset<String, String> pairwiseCorrelation(DoubleMatrixDataset<String, String> in) throws Exception {

        DoubleMatrixDataset<String, String> output = new DoubleMatrixDataset<>(in.rows(), in.rows());
        output.setRowObjects(in.getRowObjects());
        output.setColObjects(in.getRowObjects());


        double[][] data = in.getMatrix().toArray();
        ProgressBar pb = new ProgressBar(in.rows(), "Calculating correlation matrix");
        IntStream.range(0, in.rows()).parallel().forEach(row -> {
            double[] xvals = data[row];
            for (int i = row + 1; i < in.rows(); i++) {
                double[] yvals = data[i];
                double r = Correlation.correlate(xvals, yvals);
                output.setElementQuick(row, i, r);
                output.setElementQuick(i, row, r);
            }
            output.setElementQuick(row, row, 1d);
            pb.iterateSynched();
        });
        pb.close();

        return output;
    }

    public DoubleMatrix2D pairwiseCorrelationDoubleMatrix(double[][] in) {
        ExecutorService threadPool = Executors.newFixedThreadPool(nrThreads);
        CompletionService<Pair<Integer, double[]>> pool = new ExecutorCompletionService<Pair<Integer, double[]>>(threadPool);
        double meanOfSamples[] = new double[in.length];

        for (int i = 0; i < meanOfSamples.length; ++i) {
            meanOfSamples[i] = Descriptives.mean(in[i]);
        }

        for (int i = 0; i < in.length; i++) {
            ConcurrentCorrelationTask task = new ConcurrentCorrelationTask(in, meanOfSamples, i);
            pool.submit(task);
        }

        int returned = 0;

        DoubleMatrix2D correlationMatrix;
        if ((in.length * (long) in.length) > (Integer.MAX_VALUE - 2)) {
            correlationMatrix = new DenseLargeDoubleMatrix2D(in.length, in.length);
        } else {
            correlationMatrix = new DenseDoubleMatrix2D(in.length, in.length);
        }


        ProgressBar pb = new ProgressBar(in.length, "Calculation of correlation matrix: " + in.length + " x " + in.length);
        while (returned < in.length) {
            try {
                Pair<Integer, double[]> result = pool.take().get();
                if (result != null) {
                    int rownr = result.getLeft(); //  < 0 when row is not to be included because of hashProbesToInclude.
                    if (rownr >= 0) {
                        double[] doubles = result.getRight();
                        for (int i = 0; i < doubles.length; ++i) {
                            correlationMatrix.setQuick(rownr, i, doubles[i]);
                        }
                    }
                    result = null;
                    returned++;
                    pb.iterate();
                }
            } catch (Exception e) {
                e.printStackTrace();
            }
        }

        for (int r = 1; r < correlationMatrix.rows(); r++) {
            for (int c = 0; c < r; c++) {
                correlationMatrix.setQuick(r, c, correlationMatrix.getQuick(c, r));
            }
        }

        threadPool.shutdown();
        pb.close();
        return correlationMatrix;
    }
}
