/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.binarymeta.meta.cis;

import eqtlmappingpipeline.binarymeta.meta.MetaAnalysisCalculationThread;
import java.util.HashMap;
import java.util.concurrent.ArrayBlockingQueue;
import org.apache.logging.log4j.LogManager;
import umcg.genetica.containers.Pair;

/**
 *
 * @author harmjan
 */
public class CalcThread extends Thread {

    private final ArrayBlockingQueue<BinaryUnzipTask> input;
    private final ArrayBlockingQueue<Pair<Integer, HashMap<Integer, Float>>> output;
    private boolean poisoned;

    public CalcThread(ArrayBlockingQueue<BinaryUnzipTask> input, ArrayBlockingQueue<Pair<Integer, HashMap<Integer, Float>>> output) {
        this.input = input;
        this.output = output;
    }

    @Override
    public void run() {

        while (!poisoned) {
            try {
                BinaryUnzipTask task = input.take();
                if (task.isPoison()) {
                    poisoned = true;
                } else {
                    Pair<Integer, HashMap<Integer, Float>> result = task.call();
                    output.offer(result);
                }
            } catch (InterruptedException ex) {
                LogManager.getLogger(MetaAnalysisCalculationThread.class.getName()).log(org.apache.logging.log4j.Level.FATAL, ex);
            } catch (Exception ex) {
                LogManager.getLogger(MetaAnalysisCalculationThread.class.getName()).log(org.apache.logging.log4j.Level.FATAL, ex);
            }

        }
    }
}
