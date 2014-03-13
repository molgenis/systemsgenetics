/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.binarymeta.meta.cis;

import java.io.IOException;
import java.util.HashMap;
import java.util.concurrent.ArrayBlockingQueue;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.bin.BinaryFile;

/**
 *
 * @author harm-jan
 */
public class ResultProcessor extends Thread {

    private final ArrayBlockingQueue<Pair<Integer, HashMap<Integer, Float>>> pool;
    private final int totalNrOfSNPs;
    private final String outputfile;
    private boolean printRemaining;

    public ResultProcessor(ArrayBlockingQueue<Pair<Integer, HashMap<Integer, Float>>> pool, int totalNrOfSNPs, String outputfile) {
	this.pool = pool;
	this.totalNrOfSNPs = totalNrOfSNPs;
	this.outputfile = outputfile;
	System.out.println("Result thread for " + totalNrOfSNPs + " snps calls for duty :)");
    }

    @Override
    public void run() {
	try {
	    BinaryFile bf = new BinaryFile(outputfile, true);
	    int returnedResults = 0;
	    ProgressBar pb = new ProgressBar(totalNrOfSNPs, "Now extracting");
	    while (returnedResults < totalNrOfSNPs) {
		try {
		    Pair<Integer, HashMap<Integer, Float>> result = pool.take();
		    HashMap<Integer, Float> data = result.getRight();
		    Integer snpid = result.getLeft();
		    if (snpid.equals(-1)) {
			printRemaining = true;
			System.out.println("Remaining: " + (totalNrOfSNPs - returnedResults));

		    } else {
			if (printRemaining) {
			    System.out.println("Remaining: " + (totalNrOfSNPs - returnedResults));
			}
			for (int d = 0; d < data.size(); d++) {
			    Float datapoint = data.get(d);
			    if (datapoint != null) {
				bf.writeInt(snpid);
				bf.writeInt(d);
				bf.writeFloat(datapoint);
				datapoint = null;
				data.put(d, null);
			    }
			}
			data = null;
			returnedResults++;
//			if (returnedResults % 100000 == 0) {
//			    System.gc();
//			}
		    }
		    pb.iterate();
		} catch (Exception e) {
		    e.printStackTrace();
		}
	    }
	    pb.close();
	    bf.close();
	    System.out.println("Thread died!?");
	} catch (IOException e) {
	    e.printStackTrace();
	    System.exit(0);
	}
    }
}
