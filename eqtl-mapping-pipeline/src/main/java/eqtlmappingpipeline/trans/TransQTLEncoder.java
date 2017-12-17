package eqtlmappingpipeline.trans;

import java.util.BitSet;

public class TransQTLEncoder {
	
	private static byte[] convertShortMatrixToByteArrat(short[][] matrix) {
		int len = matrix.length * matrix[0].length * 2; // 2 bytes per short
		byte[] concate = new byte[len];
		int pos = 0;
		for (int i = 0; i < matrix.length; i++) {
			byte[] conv = convertShortArr(matrix[i]);
			System.arraycopy(conv, 0, concate, pos, conv.length);
			pos += conv.length;
		}
		return concate;
	}
	
	public static byte[] convertShortToByteArray(short x) {
		byte[] ret = new byte[2];
		ret[0] = (byte) x;
		ret[1] = (byte) (x >> 8);
		return ret;
	}
	
	public static byte[] convertShortArr(short[] x) {
		byte[] ret = new byte[x.length * 2];
		for (int i = 0; i < x.length; i++) {
			short s = x[i];
			int idx = i * 2;
			ret[idx] = (byte) (s >> 8);
			ret[idx + 1] = (byte) s;
		}
		return ret;
	}
	
	public BitSet convertLonhToBitset(long value) {
		BitSet bits = new BitSet();
		int index = 0;
		while (value != 0L) {
			if (value % 2L != 0) {
				bits.set(index);
			}
			++index;
			value = value >>> 1;
		}
		return bits;
	}
	
	public long convertBitsetToLong(BitSet bits) {
		long value = 0L;
		for (int i = 0; i < bits.length(); ++i) {
			value += bits.get(i) ? (1L << i) : 0L;
		}
		return value;
	}
	
}
