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
	static final int[] LAST_BYTE_MASK = {1, 3, 7, 15, 31, 63, 127, 255};
	/**
	 * Index is the number of first bits used from the last byte
	 */
	static final int[] FIRST_BYTE_MASK = {0, 128, 192, 224, 240, 248, 252, 254, 255};

	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) {

		//byte[] test = {(byte) 0b0011_0000, (byte) 0b1111_0011, (byte) 0b1011_0011, (byte) 0b0000_1101, (byte) 0b1111_1111};
		//byte[] test = {(byte) 0b0000_1000, (byte) 0b1111_1111, (byte) 0b1111_1111, (byte) 0b11111_1111, (byte) 0b1111_1111};
		byte[] test = {(byte) 0b0011_0000, (byte) 0b1111_0011, (byte) 0b1011_0011, (byte) 0b0000_1101, (byte) 0b1111_1111};
		
//		System.out.println(Integer.toBinaryString(test[0] & (255)));
//		System.out.println(Integer.toBinaryString(test[1] & (255)));
//		System.out.println(Integer.toBinaryString(test[2] & (255)));
//		System.out.println(Integer.toBinaryString(test[3] & (255)));
//		System.out.println(Integer.toBinaryString(test[4] & (255)));

		System.out.println(String.format("%8s", Integer.toBinaryString(test[0] & 0xFF)).replace(' ', '0'));
		System.out.println(String.format("%8s", Integer.toBinaryString(test[1] & 0xFF)).replace(' ', '0'));
		System.out.println(String.format("%8s", Integer.toBinaryString(test[2] & 0xFF)).replace(' ', '0'));
		System.out.println(String.format("%8s", Integer.toBinaryString(test[3] & 0xFF)).replace(' ', '0'));
		System.out.println(String.format("%8s", Integer.toBinaryString(test[4] & 0xFF)).replace(' ', '0'));

		int bits = 25;
		long factor = (1L << bits) - 1;

//		System.out.println("Prob divide factor: " + factor);
//
//		readProb(test, 1, 3, bits, factor);

		bits = 1;
		factor = (1L << bits) - 1;

//		System.out.println("Prob divide factor: " + factor);

//		readProb(test, 0, 4, bits, factor);
		
//		System.out.println("---------");
		
		bits = 23;
		factor = (1L << bits) - 1;

		System.out.println("Prob divide factor: " + factor);

		readProb(test, 8, bits);
		double intOfNBits = getIntOfNBits(test, 8, bits);
		System.out.println("intOfNBits = " + intOfNBits / factor);

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
	private static double readProb(byte[] bytes, int bitOffset, int totalBitsToRead) {

		// Get the byte to start reading from.
		int firstByteIndex = bitOffset / 8;
		// Get the bit within the byte to start reading from.
		int remainingBitOffset = bitOffset % 8;

		int totalBytesMin1 = (remainingBitOffset + totalBitsToRead - 1) / 8;


		//int bitsFromLastByte2 = totalBits - ((totalBytesMin1 - 1) * 8) - (8 - remainingBitOffset);
		//Below is simplification of real formula above
		int nBitsFromLastByte = totalBitsToRead - (totalBytesMin1 * 8) + remainingBitOffset;

		//Below is old way to calculate nBitsFromLastByte but I think above is faster
//		int nBitsFromLastByte = (totalBits + remainingBitOffset) % 8;
//		if (nBitsFromLastByte == 0) {
//			nBitsFromLastByte = 8;
//		}



		int bitShiftAfterFirstByte = 8 - remainingBitOffset;

		System.out.println("Total bytes min 1: " + totalBytesMin1);
		System.out.println("First byte: " + firstByteIndex);
		System.out.println("Index of first bit first byte: " + remainingBitOffset);
		System.out.println("Total bits: " + totalBitsToRead);
		System.out.println("Bits from last byte: " + nBitsFromLastByte);
		System.out.println("First byte mask: " + Integer.toBinaryString(LAST_BYTE_MASK[remainingBitOffset]));
		System.out.println("Last byte mask: " + Integer.toBinaryString(FIRST_BYTE_MASK[nBitsFromLastByte]));
		System.out.println("Bit shift after first byte: " + bitShiftAfterFirstByte);

		long value; // long because values are stored as unsigned int

		//Switch because this is very fast compared to loops or if statements
		switch (totalBytesMin1) {
			case 0:
				value
						= bytes[firstByteIndex] >> (remainingBitOffset) & LAST_BYTE_MASK[8 - totalBitsToRead];
				break;
			case 1:
				//Last byte parsing is different then longer encoding
				value
						= (bytes[firstByteIndex + 1] & LAST_BYTE_MASK[nBitsFromLastByte]) << bitShiftAfterFirstByte
						| (bytes[firstByteIndex] & FIRST_BYTE_MASK[bitShiftAfterFirstByte]) >> remainingBitOffset;
				break;
			case 2:
				value
						= (bytes[firstByteIndex + 2] & LAST_BYTE_MASK[nBitsFromLastByte]) << (8 + bitShiftAfterFirstByte)
						| (bytes[firstByteIndex + 1] & 255) << (bitShiftAfterFirstByte)
						| (bytes[firstByteIndex] & FIRST_BYTE_MASK[bitShiftAfterFirstByte]) >> remainingBitOffset;
				break;
			case 3:
				value
						= ((long) bytes[firstByteIndex + 3] & LAST_BYTE_MASK[nBitsFromLastByte]) << (16 + bitShiftAfterFirstByte)
						| (bytes[firstByteIndex + 2] & 255) << (8 + bitShiftAfterFirstByte)
						| (bytes[firstByteIndex + 1] & 255) << (bitShiftAfterFirstByte)
						| (bytes[firstByteIndex] & FIRST_BYTE_MASK[bitShiftAfterFirstByte]) >> remainingBitOffset;
				break;
			case 4:
				value
						= ((long) bytes[firstByteIndex + 4] & LAST_BYTE_MASK[nBitsFromLastByte]) << (24 + bitShiftAfterFirstByte)
						| (bytes[firstByteIndex + 3] & 255) << (16 + bitShiftAfterFirstByte)
						| (bytes[firstByteIndex + 2] & 255) << (8 + bitShiftAfterFirstByte)
						| (bytes[firstByteIndex + 1] & 255) << (bitShiftAfterFirstByte)
						| (bytes[firstByteIndex] & FIRST_BYTE_MASK[bitShiftAfterFirstByte]) >> remainingBitOffset;
				break;

			default:
				throw new GenotypeDataException("Error parsing bgen file. Debug info: totalBits=" + totalBitsToRead + " remainingBitOffset=" + remainingBitOffset + " totalBits=" + totalBitsToRead);
		}

		System.out.println("Result: " + value);
		System.out.println("Result in bin: " + Long.toBinaryString(value));

		return value;
	}


	/**
	 * Converts a specified range of bits from a byte array into a 'long' value.
	 *
	 * @param bytes The byte array to get the result value from.
	 * @param bitOffset The bit in the byte array to start reading from.
	 * @param totalBitsToRead The total number of bits to read from the byte array.
	 * @return The specified bits converted to a 'long' value.
	 */
	private static long getIntOfNBits(byte[] bytes, long bitOffset, int totalBitsToRead) {
		// Get the byte to start reading from.
		int byteOffset = Math.toIntExact(bitOffset / 8);
		// Get the bit within the byte to start reading from.
		int remainingBitOffset = Math.toIntExact(bitOffset % 8);
		// Get the number of bits to read within the initial byte.
		long nBitsToRead = 8 - remainingBitOffset;

		// Define result value and a counter for the number of read bits.
		long nReadBits = 0;
		long value = 0;

		System.out.println("bitOffset = " + bitOffset);

		// First read the number of bits that should be read from the first byte (8 - remainingBitOffset).
		// After that, if the condition still applies, read all bits of following bytes, placing these bits
		// in before the previous read bits (<< readbits)
		while (nReadBits + 8 < totalBitsToRead) {
			// Read all bits in the current byte, moving them the number of already read bits to the left,
			// and masking all bits that are not within the bits that have previously been read or should be read now.
			System.out.println("nBitsToRead = " + nBitsToRead);
			value |= bytes[byteOffset] << nReadBits & ((1L << (nBitsToRead + nReadBits)) - 1);

			// Update bit reading numbers
			nReadBits += nBitsToRead;
			nBitsToRead = 8;

			String s1 = String.format("%8s", Integer.toBinaryString(bytes[byteOffset] & 0xFF)).replace(' ', '0');
			System.out.println("value = " + s1 + " | " + Long.toString(value,2));

			// Update offset values
			byteOffset += 1;
			remainingBitOffset = 0;
		}

		// Calculate the number of bits that still have to be read.
		nBitsToRead = totalBitsToRead - (nReadBits - remainingBitOffset);
		System.out.println("nBitsToRead = " + nBitsToRead);

		// First, move the bits that have to be read in the last byte to the right boundary of the byte
		// (shifting amount defined by: 8 - nBitsToRead and mask the bits that should not be
		// read.
		long readBits = bytes[byteOffset] >> (8 - nBitsToRead) & ((1 << nBitsToRead) - 1);
		// Secondly, move the read bits to the leftmost part of the result value, and mask all other bits.
		value |= readBits << nReadBits & ((1L << totalBitsToRead) - 1);

		String s1 = String.format("%8s", Integer.toBinaryString(bytes[byteOffset] & 0xFF)).replace(' ', '0');
		System.out.println("value = " + s1 + " | " + Long.toString(value,2) + " | " + value);
		return value;
	}

}

