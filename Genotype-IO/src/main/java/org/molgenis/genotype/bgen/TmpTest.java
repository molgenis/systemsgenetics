package org.molgenis.genotype.bgen;

import org.molgenis.genotype.GenotypeDataException;

/**
 *
 * @author patri
 */
public class TmpTest {

	/**
	 * Index is the number of last bits used from the first byte
	 */
	static final int[] FIRST_BYTE_MASK = {255, 127, 63, 31, 15, 7, 3, 1};
	/**
	 * Index is the number of first bits used from the last byte
	 */
	static final int[] LAST_BYTE_MASK = {0, 128, 192, 224, 240, 248, 252, 254, 255};

	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) {

		//byte[] test = {(byte) 0b0011_0000, (byte) 0b1111_0011, (byte) 0b1011_0011, (byte) 0b0000_1101, (byte) 0b1111_1111};
		//byte[] test = {(byte) 0b0000_1000, (byte) 0b1111_1111, (byte) 0b1111_1111, (byte) 0b11111_1111, (byte) 0b1111_1111};
		byte[] test = {(byte) 0b0011_0000, (byte) 0b1111_0011, (byte) 0b1011_0011, (byte) 0b0000_1101, (byte) 0b1111_1111};
		
		System.out.println(Integer.toBinaryString(test[0] & (255)));
		System.out.println(Integer.toBinaryString(test[1] & (255)));
		System.out.println(Integer.toBinaryString(test[2] & (255)));
		System.out.println(Integer.toBinaryString(test[3] & (255)));
		System.out.println(Integer.toBinaryString(test[4] & (255)));

		int bits = 23;
		long factor = (1L << bits) - 1;

		System.out.println("Prob divide factor: " + factor);

		readProb(test, 1, 3, bits, factor);

		bits = 1;
		factor = (1L << bits) - 1;

		System.out.println("Prob divide factor: " + factor);

		readProb(test, 0, 4, bits, factor);
		
		System.out.println("---------");
		
		bits = 32;
		factor = (1L << bits) - 1;

		System.out.println("Prob divide factor: " + factor);

		readProb(test, 1, 0, bits, factor);

	}

	/**
	 *
	 * @param bytes
	 * @param firstByteIndex zero based index of first byte in byte array
	 * @param indexBitInFirstByte zero based index of first used bit in first
	 * byte
	 * @param totalBits
	 * @param conversionFactor
	 * @return
	 */
	private static double readProb(byte[] bytes, int firstByteIndex, int indexBitInFirstByte, int totalBits, long conversionFactor) {

		int totalBytesMin1 = (totalBits + indexBitInFirstByte - 1) / 8;
		
		
		//int bitsFromLastByte2 = totalBits - ((totalBytesMin1 - 1) * 8) - (8 - indexBitInFirstByte);
		//Below is simplification of real formula above
		int bitsFromLastByte = totalBits - (totalBytesMin1 * 8) + indexBitInFirstByte;
		
		//Below is old way to calculate bitsFromLastByte but I think above is faster
//		int bitsFromLastByte = (totalBits + indexBitInFirstByte) % 8;
//		if (bitsFromLastByte == 0) {
//			bitsFromLastByte = 8;
//		}
		
		

		int bitshiftAfterFirstByte = 8 - indexBitInFirstByte;

		System.out.println("Total bytes min 1: " + totalBytesMin1);
		System.out.println("First byte: " + firstByteIndex);
		System.out.println("Index of first bit first byte: " + indexBitInFirstByte);
		System.out.println("Total bits: " + totalBits);
		System.out.println("Bits from last byte: " + bitsFromLastByte);
		System.out.println("First byte mask: " + Integer.toBinaryString(FIRST_BYTE_MASK[indexBitInFirstByte]));
		System.out.println("Last byte mask: " + Integer.toBinaryString(LAST_BYTE_MASK[bitsFromLastByte]));
		System.out.println("Bit shift after first byte: " + bitshiftAfterFirstByte);

		long encodedProb; // long because values are stored as unsigned int

		//Switch because this is very fast compared to loops or if statements
		switch (totalBytesMin1) {
			case 0:
				encodedProb
						= bytes[firstByteIndex] >> (8 - indexBitInFirstByte - totalBits) & FIRST_BYTE_MASK[8 - totalBits];
				break;
			case 1:
				//Last byte parsing is different then longer encoding
				encodedProb
						= (bytes[firstByteIndex + 1] & LAST_BYTE_MASK[bitsFromLastByte]) >> (8 - bitsFromLastByte) << bitshiftAfterFirstByte
						| (bytes[firstByteIndex] & FIRST_BYTE_MASK[indexBitInFirstByte]);
				break;
			case 2:
				encodedProb
						= (bytes[firstByteIndex + 2] & LAST_BYTE_MASK[bitsFromLastByte]) << (bitshiftAfterFirstByte + bitsFromLastByte)
						| (bytes[firstByteIndex + 1] & 255) << (indexBitInFirstByte)
						| (bytes[firstByteIndex] & FIRST_BYTE_MASK[indexBitInFirstByte]);
				break;
			case 3:
				encodedProb
						= ((long) bytes[firstByteIndex + 3] & LAST_BYTE_MASK[bitsFromLastByte]) << (8 + bitshiftAfterFirstByte + bitsFromLastByte)
						| (bytes[firstByteIndex + 2] & 255) << (8 + bitshiftAfterFirstByte)
						| (bytes[firstByteIndex + 1] & 255) << (bitshiftAfterFirstByte)
						| (bytes[firstByteIndex] & FIRST_BYTE_MASK[indexBitInFirstByte]);
				break;
			case 4:
				encodedProb
						= ((long) bytes[firstByteIndex + 4] & LAST_BYTE_MASK[bitsFromLastByte]) << (16 + bitshiftAfterFirstByte + bitsFromLastByte)
						| (bytes[firstByteIndex + 3] & 255) << (16 + indexBitInFirstByte)
						| (bytes[firstByteIndex + 2] & 255) << (8 + indexBitInFirstByte)
						| (bytes[firstByteIndex + 1] & 255) << (indexBitInFirstByte)
						| (bytes[firstByteIndex] & FIRST_BYTE_MASK[indexBitInFirstByte]);
				break;

			default:
				throw new GenotypeDataException("Error parsing bgen file. Debug info: totalBits=" + totalBits + " indexBitInFirstByte=" + indexBitInFirstByte + " totalBits=" + totalBits);
		}

		System.out.println("Result: " + encodedProb);
		System.out.println("Result in bin: " + Long.toBinaryString(encodedProb));

		return ((double) encodedProb) / conversionFactor;

	}

}

