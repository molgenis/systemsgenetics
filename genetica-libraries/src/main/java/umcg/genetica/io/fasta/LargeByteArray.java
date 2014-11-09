package umcg.genetica.io.fasta;

/**
 *
 * @author Patrick Deelen
 */
public class LargeByteArray {

	private final byte[][] data;
	private long size;
	private int blockSize;
	private static final int MAX_BLOCK_SIZE = Integer.MAX_VALUE - 8;

	public LargeByteArray(long size) {

		this.size = size;

		long blocks = size / MAX_BLOCK_SIZE;
		if (size % MAX_BLOCK_SIZE > 0) {
			++blocks;
		}

		//Also prevent of creating more blocks than allowed by array.
		if (blocks >= MAX_BLOCK_SIZE) {
			throw new IllegalArgumentException("Array to large");
		}

		blockSize = (int) (size / blocks);
		if (size % blocks > 0) {
			++blockSize;
		}

//		System.out.println("blocks: " + (int) blocks);
//		System.out.println("blockSize: " + blockSize);

		//Cast is save because of check
		data = new byte[(int) blocks][blockSize];

	}

	/**
	 * Set value at index. No check is done if within range
	 *
	 * @param index
	 * @param value
	 */
	public void setQuick(long index, byte value) {

		int block = (int) (index / blockSize);
		int blockIndex = (int) (index % blockSize);
//		try{
		data[block][blockIndex] = value;
//		} catch (ArrayIndexOutOfBoundsException ex){
//			System.err.println("ERROR out of bound");
//			System.out.println("size " + size);
//			System.out.println("blockSize " + blockSize);
//			System.out.println("blocks " + data.length);
//			System.out.println("query index " + index);
//			System.out.println("block " + block);
//			System.out.println("block index " + blockIndex);
//			throw ex;
//		}

	}

	/**
	 * Get value at index. No check is done if within range
	 *
	 * @param index
	 * @param value
	 */
	public byte getQuick(long index) {

		int block = (int) (index / blockSize);
		int blockIndex = (int) (index % blockSize);
		return data[block][blockIndex];

	}

	public long getSize() {
		return size;
	}
}
