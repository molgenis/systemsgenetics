/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.binarymeta.meta.cis;

import java.nio.ByteBuffer;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.concurrent.Callable;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.zip.DataFormatException;
import java.util.zip.Inflater;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.trityper.bin.BinaryResultDataset;
import umcg.genetica.io.trityper.bin.BinaryResultSNP;

/**
 *
 * @author harmjan
 */
public class UnzipTask implements Callable<Pair<BinaryResultSNP, HashMap<String, Float>>> {

    private final int snp;
    private final Inflater inflater = new Inflater();
    private final BinaryResultDataset ds;
    private final int d;
    private final boolean[] probePresent;
    private final ArrayList<String> probes;
    private final Integer[][] probeTranslationLookupTable;
    private final BinaryResultSNP[] binsnps;

    public UnzipTask(int snpInDataset, BinaryResultDataset ds, int d, Integer metaSNPId, boolean[] probePresent, ArrayList<String> probes, Integer[][] probeTranslationLookupTable) {
	this.snp = snpInDataset;
	this.ds = ds;
	this.d = d;
	this.probePresent = probePresent;
	this.probes = probes;
	this.probeTranslationLookupTable = probeTranslationLookupTable;
	this.binsnps = ds.getSnps();
    }

    @Override
    public Pair<BinaryResultSNP, HashMap<String, Float>> call() throws Exception {



	HashMap<String, Float> outputData = new HashMap<String, Float>();
//	for (int snp = 0; snp < binsnps.length; snp++) {
//	if (snpsForWhichThereAreEffects[snp]) {
	BinaryResultSNP snpObject = binsnps[snp];
	Integer binSNPId = snpObject.getId();
//	    Integer metaSNPId = reverseSNPLookupTable[d][binSNPId];

	// load the zscores
	long pointer = snpObject.getzScoreIndex();
	long nextpointer = -1;

	if (binSNPId + 1 < ds.getSnps().length) {
	    BinaryResultSNP snpObject2 = ds.getSnps()[binSNPId + 1];
	    nextpointer = snpObject2.getzScoreIndex();
	}

	Float[] zscores = null;
	try {
	    byte[] data = ds.getMatrix().readDeflated(pointer, nextpointer, ds.getNumProbes());
	    zscores = inflate(data, ds.getNumProbes());
	} catch (DataFormatException ex) {
	    Logger.getLogger(CisAnalysis.class.getName()).log(Level.SEVERE, null, ex);
	}

	for (int probe = 0; probe < probes.size(); probe++) {
	    if (probePresent[probe]) {
		Integer dsProbeId = probeTranslationLookupTable[d][probe];
		if (dsProbeId != null) {
		    outputData.put(probes.get(probe), zscores[dsProbeId]);
		}
	    }
	}
	return new Pair<BinaryResultSNP, HashMap<String, Float>>(snpObject, outputData);
    }

    private Float[] inflate(byte[] buffer, int numElems) throws DataFormatException {
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
