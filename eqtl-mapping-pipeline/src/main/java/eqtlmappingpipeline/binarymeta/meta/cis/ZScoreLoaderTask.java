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
public class ZScoreLoaderTask implements Callable<Pair<Integer, Float[]>> {

    private final int d;
    private final boolean[] snpsForWhichThereAreEffects;
    private final Integer[][] reverseSNPLookupTable;
    private final BinaryResultDataset[] ds;
    private final ArrayList<String> probes;
    private final Integer[][] probeTranslationLookupTable;
    private final HashMap<Pair<Integer, Integer>, Integer> SNPProbeToEffectMap;
    private final Float[] output;
    private final Inflater inflater = new Inflater();
    private final boolean[] probePresent;

    public ZScoreLoaderTask(BinaryResultDataset[] ds, int d, boolean[] snpsForWhichThereAreEffects, Integer[][] reverseSNPLookupTable, ArrayList<String> probes, Integer[][] probeTranslationLookupTable, HashMap<Pair<Integer, Integer>, Integer> SNPProbeToEffectMap, Float[] output, boolean[] probePresent) {
	this.ds = ds;
	this.d = d;
	this.snpsForWhichThereAreEffects = snpsForWhichThereAreEffects;
	this.reverseSNPLookupTable = reverseSNPLookupTable;
	this.probes = probes;
	this.probeTranslationLookupTable = probeTranslationLookupTable;
	this.SNPProbeToEffectMap = SNPProbeToEffectMap;
	this.output = output;
	this.probePresent = probePresent;
    }

    @Override
    public Pair<Integer, Float[]> call() throws Exception {
	BinaryResultSNP[] binsnps = ds[d].getSnps();
	for (int snp = 0; snp < binsnps.length; snp++) {
	    if (snpsForWhichThereAreEffects[snp]) {
		BinaryResultSNP snpObject = binsnps[snp];
		Integer binSNPId = snpObject.getId();
		Integer metaSNPId = reverseSNPLookupTable[d][binSNPId];

		// load the zscores
		long pointer = snpObject.getzScoreIndex();
		long nextpointer = -1;

		if (binSNPId + 1 < ds[d].getSnps().length) {
		    BinaryResultSNP snpObject2 = ds[d].getSnps()[binSNPId + 1];
		    nextpointer = snpObject2.getzScoreIndex();
		}
		Float[] zscores = null;
		try {
		    byte[] data = ds[d].getMatrix().readDeflated(pointer, nextpointer, ds[d].getNumProbes());
		    zscores = inflate(data, ds[d].getNumProbes());
		} catch (DataFormatException ex) {
		    Logger.getLogger(CisAnalysis.class.getName()).log(Level.SEVERE, null, ex);
		}


		for (int probe = 0; probe < probes.size(); probe++) {
		    if (probePresent[probe]) {
			Integer dsProbeId = probeTranslationLookupTable[d][probe];
			if (dsProbeId != null) {
			    Pair<Integer, Integer> snpProbePair = new Pair<Integer, Integer>(metaSNPId, probe);
			    Integer pairId = SNPProbeToEffectMap.get(snpProbePair);
			    if (pairId != null) {
				output[pairId] = zscores[dsProbeId];
			    }
			}
		    }
		}
	    }
	    if (snp > 0 && snp % 1000 == 0) {
		System.out.println("Dataset\t" + d + "\t" + snp + "/" + binsnps.length + "\t"+((double) snp/binsnps.length));
	    }
	}

	return new Pair<Integer, Float[]>(d, output);
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
