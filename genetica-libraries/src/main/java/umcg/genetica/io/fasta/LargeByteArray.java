package umcg.genetica.io.fasta;

import cern.colt.matrix.tdouble.impl.DenseLargeDoubleMatrix2D;

/**
 *
 * @author Patrick Deelen
 */
public class LargeByteArray {
	
	private final byte[][] data;
	private long size;
	private static final int BLOCK_SIZE = Integer.MAX_VALUE - 8;

	public LargeByteArray(long size) {
		
		this.size = size;
		
		long blocks = size / BLOCK_SIZE;
		long left = size % BLOCK_SIZE;
		if(left > 0){
			++blocks;
		}
		
		//Also prevent of creating more blocks than allowed by array.
		if(blocks >= BLOCK_SIZE){
			throw new IllegalArgumentException("Array to large");
		}
		
		//Cast is save because of check
		data = new byte[(int) blocks][BLOCK_SIZE];
		
	}
	
	/**
	 * Set value at index. No check is done if within range
	 * 
	 * @param index
	 * @param value 
	 */
	public void setQuick(long index, byte value){
		
		int block = (int) (index / BLOCK_SIZE);
		int blockIndex = (int) (index % BLOCK_SIZE);
		data[block][blockIndex] = value;
		
	}
	
	/**
	 * Get value at index. No check is done if within range
	 * 
	 * @param index
	 * @param value 
	 */
	public byte getQuick(long index){
		
		int block = (int) (index / BLOCK_SIZE);
		int blockIndex = (int) (index % BLOCK_SIZE);
		return data[block][blockIndex];
		
	}
	

}
