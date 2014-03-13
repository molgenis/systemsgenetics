/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.binarymeta.meta.cis;

import eqtlmappingpipeline.binarymeta.meta.MetaSettings;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.zip.DataFormatException;
import java.util.zip.Inflater;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.bin.BinaryFile;
import umcg.genetica.io.trityper.bin.BinaryResultDataset;
import umcg.genetica.io.trityper.bin.BinaryResultProbe;
import umcg.genetica.io.trityper.bin.BinaryResultSNP;

/**
 *
 * @author harm-jan
 */
public class CisBinaryConverter {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
	// TODO code application logic here
//	if (args.length < 1) {
//	    System.out.println("Usage: converter settings.xml");
//	    System.exit(0);
//	}
        String settingsfile = "/Volumes/Data2/MetaAnalysisSoftware/settings/2011-12-06-CIS-40PCs-4GWAS-GeneticVectorsNotRemoved.xml";
	CisBinaryConverter c = new CisBinaryConverter(settingsfile);
	try {
	    c.run();
	} catch (IOException ex) {
	    Logger.getLogger(CisBinaryConverter.class.getName()).log(Level.SEVERE, null, ex);
	}
    }
    private final MetaSettings m_settings;

    public CisBinaryConverter(String settings) {
	m_settings = new MetaSettings();
	m_settings.parse(settings, null, null);

//         m_settings =  null;
    }

    public void run() throws IOException {

	ArrayList<String> datasetLocations = m_settings.getDatasetlocations();
	ArrayList<String> datasetNames = m_settings.getDatasetnames();
	ArrayList<String> datasetAnnot = m_settings.getDatasetannotations();
//        ArrayList<String> datasetLocations = new ArrayList<String>();
//        datasetLocations.add("d:\\BinData\\");
//        ArrayList<String> datasetNames = new ArrayList<String>();
//        datasetNames.add("Dataset");
	int nrPerm = m_settings.getNrPermutations();
//        int nrPerm = 1;
	int nrprocs = Runtime.getRuntime().availableProcessors() / 2;
//	ExecutorService threadPool = Executors.newFixedThreadPool(nrprocs);

	System.out.println(nrprocs + " procs.");



//	ArrayBlockingQueue<BinaryUnzipTask> inputqueue = new ArrayBlockingQueue<BinaryUnzipTask>(nrprocs);
//	ArrayBlockingQueue<Pair<Integer, HashMap<Integer, Float>>> outputqueue = new ArrayBlockingQueue<Pair<Integer, HashMap<Integer, Float>>>(nrprocs);
////	CompletionService<Pair<Integer, HashMap<Integer, Float>>> pool = new ExecutorCompletionService<Pair<Integer, HashMap<Integer, Float>>>(threadPool, queue);

	Gpio.createDir(m_settings.getOutput());
	System.out.println("Found " + datasetLocations.size() + " datasets..");
	for (int d = 0; d < datasetLocations.size(); d++) {



	    System.out.println("Starting for " + datasetNames.get(d));
	    for (int perm = 0; perm < nrPerm + 1; perm++) {

		String datasetName = "Dataset";

		BinaryResultDataset dataset = new BinaryResultDataset(datasetLocations.get(d), datasetName, perm);
		BinaryResultSNP[] snps = dataset.getSnps();
		BinaryResultProbe[] probes = dataset.getProbes();

		int nrProbes = dataset.getNumProbes();

		String outfile = null;
		if (perm == 0) {
		    outfile = m_settings.getOutput() + datasetNames.get(d) + "-eQTLs.dat";
		} else {
		    outfile = m_settings.getOutput() + datasetNames.get(d) + "-PermutedDataPermutationRound-" + perm + ".dat";
		}

		BinaryFile bf = new BinaryFile(outfile, true);

//		CalcThread[] threads = new CalcThread[nrprocs];
//		for(int i=0; i<threads.length; i++){
//		    threads[i] = new CalcThread(inputqueue, outputqueue);
//		    threads[i].start();
//		}

		System.out.println("Outfile: " + outfile);
//		Thread completor = new ResultProcessor(outputqueue, snps.length, outfile);
////		Thread completor = new ResultProcessor(pool, 1000, outfile);
//		completor.start();
//
//		completor.setName("Completion-service");

//                for (int snp = 0; snp < snps.length; snp+=nrprocs) {
		int snp = 0;
		dataset.clearProbeObjects();
		ProgressBar pb = new ProgressBar(snps.length);
		for (snp = 0; snp < snps.length; snp++) {
//		for (snp = 0; snp < 100000; snp++) {
		    // load the zscores

		    BinaryResultSNP snpObject = snps[snp];
		    long pointer = snpObject.getzScoreIndex();
		    long nextpointer = -1;

		    if (snp + 1 < snps.length) {
			BinaryResultSNP snpObject2 = snps[snp + 1];
			nextpointer = snpObject2.getzScoreIndex();
		    }

		    byte[] bindata = dataset.getMatrix().readDeflated(pointer, nextpointer, dataset.getNumProbes());
		    try {


			Float[] data = inflate(bindata, dataset.getNumProbes());
			bindata = null;
			Integer snpid = snp;


			for (int dp = 0; dp < data.length; dp++) {
			    Float datapoint = data[dp];
			    if (datapoint != null) {
				bf.writeInt(snpid);
				bf.writeInt(dp);
				bf.writeFloat(datapoint);
				datapoint = null;
			    }
			}
			data = null;
//			if (returnedResults % 100000 == 0) {
//			    System.gc();
//			}
			//		    BinaryUnzipTask task = new BinaryUnzipTask(snp, dataset, nrProbes);
			//UnzipTask task = new UnzipTask(snpId, ds[d], d, 0, probespresentarray, probes, probeTranslationLookupTable);


			//		    pool.submit(task);
			//		    inputqueue.offer(task);
			//		    pb.iterate();
			//                    snp++;
		    } catch (DataFormatException ex) {
			Logger.getLogger(CisBinaryConverter.class.getName()).log(Level.SEVERE, null, ex);
		    }
		    pb.iterate();
		}
		pb.close();
		bf.close();
//		pool.submit(new BinaryUnzipTask(-1, null, -1));

//		try {
//		    System.out.println("Waiting for result thread...");
//		    completor.join();
//		} catch (InterruptedException ex) {
//		    Logger.getLogger(CisBinaryConverter.class.getName()).log(Level.SEVERE, null, ex);
//		}
		dataset.close();
		dataset = null;
	    }
	}

	System.out.println("Done.");
    }
    private Inflater inflater = new Inflater();

    protected Float[] inflate(byte[] buffer, int numElems) throws DataFormatException {
	inflater.setInput(buffer);
	inflater.finished();
	byte[] decompressed = new byte[numElems * 4];
	inflater.inflate(decompressed);

	long actuallydecompressed = inflater.getBytesWritten();
	if (actuallydecompressed != numElems * 4) {
	    throw new DataFormatException("IO Error: uncompressed data does not correspond to the size requested\t" + actuallydecompressed + "\t" + numElems * 4);
	}

	inflater.reset();

	ByteBuffer bytebuffer = ByteBuffer.wrap(decompressed);
	Float[] output = new Float[numElems];
	int ctr = 0;
	for (int i = 0; i < numElems; i++) {
	    Float f = bytebuffer.getFloat();
	    if (f.isNaN()) {
		f = null;
	    } else {
		ctr++;
	    }
	    output[i] = f;
	}

	decompressed = null;

	if (ctr == 0) {
	    return null;
	} else {
	    return output;
	}
    }
}
