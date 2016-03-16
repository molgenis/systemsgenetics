/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.binarymeta.meta.cis;

import umcg.genetica.containers.Pair;
import umcg.genetica.io.trityper.bin.BinaryResultDataset;
import umcg.genetica.io.trityper.bin.BinaryResultSNP;

import java.nio.ByteBuffer;
import java.util.HashMap;
import java.util.concurrent.Callable;
import java.util.zip.DataFormatException;
import java.util.zip.Inflater;

/**
 * @author harm-jan
 */
public class BinaryUnzipTask implements Callable<Pair<Integer, HashMap<Integer, Float>>> {

	private final int snp;
	private BinaryResultDataset data;
	private final int numprobes;
	private final Inflater inflater = new Inflater();
	private boolean poison;

	public BinaryUnzipTask(int snp, BinaryResultDataset data, int numprobes) {
		this.snp = snp;
		this.data = data;
		this.numprobes = numprobes;
	}

	BinaryUnzipTask(int snp, int nrProbes, BinaryResultDataset dataset, BinaryResultSNP[] snps) {
		throw new UnsupportedOperationException("Not yet implemented");
	}

	@Override
	public Pair<Integer, HashMap<Integer, Float>> call() throws Exception {
		if (snp < 0) {
			return new Pair<Integer, HashMap<Integer, Float>>(-1, null);
		}
		BinaryResultSNP[] snps = data.getSnps();
		BinaryResultSNP snpObject = snps[snp];
		long pointer = snpObject.getzScoreIndex();
		long nextpointer = -1;

		if (snp + 1 < snps.length) {
			BinaryResultSNP snpObject2 = snps[snp + 1];
			nextpointer = snpObject2.getzScoreIndex();
		}

		byte[] bindata = data.getMatrix().readDeflated(pointer, nextpointer, data.getNumProbes());
		HashMap<Integer, Float> dataUnzipped = inflate(bindata, data.getNumProbes());
		bindata = null;

		return new Pair<Integer, HashMap<Integer, Float>>(snp, dataUnzipped);
	}

	private HashMap<Integer, Float> inflate(byte[] buffer, int numElems) throws DataFormatException {
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
		HashMap<Integer, Float> results = new HashMap<Integer, Float>();
		for (int i = 0; i < numElems; i++) {
			Float f = bytebuffer.getFloat();
			if (f.isNaN()) {
				f = null;
			} else {
				ctr++;
				results.put(i, f);
			}

//	    output[i] = f;
		}

		decompressed = null;
		buffer = null;

		if (ctr == 0) {
			return null;
		} else {
			return results;
		}
	}

	void setIsPoison() {
		poison = true;
	}

	boolean isPoison() {
		return poison;
	}
}
