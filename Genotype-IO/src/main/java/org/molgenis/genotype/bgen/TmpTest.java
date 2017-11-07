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
	static final int[] FIRST_BYTE_MASK = {255, 127, 63, 31, 7, 3, 1};
	/**
	 * Index is the number of first bits used from the last byte
	 */
	static final int[] LAST_BYTE_MASK = {128, 192, 224, 240, 248, 252, 254, 255};

	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) {

		byte[] test = {-1, 12, 2, 0, 0};

		System.out.println(Integer.toBinaryString(test[0]));

		readProb(test, 0, 0, 22);

	}

	/**
	 *
	 * @param bytes
	 * @param firstByte
	 * @param firstBitInFirstByte
	 * @param totalBits
	 * @return
	 */
	private static double readProb(byte[] bytes, int firstByte, int firstBitInFirstByte, int totalBits) {

		int totalBytes = (totalBits + firstBitInFirstByte) / 8;
		int bitsFromLastByteMin1;
		if (totalBytes == 0) {
			totalBytes = 1;
			bitsFromLastByteMin1 = totalBits - 1;
		} else {
			bitsFromLastByteMin1 = (totalBits + firstBitInFirstByte) % 8;
			if (bitsFromLastByteMin1 > 0) {
				++totalBytes;
				bitsFromLastByteMin1 -= 1;
			} else {
				bitsFromLastByteMin1 = 7;
			}
		}

		System.out.println("Total bytes: " + totalBytes);
		System.out.println("Fist byte: " + firstByte);
		System.out.println("First bit first byte: " + firstBitInFirstByte);
		System.out.println("Total bits: " + totalBits);
		System.out.println("Bits from last byte min 1: " + bitsFromLastByteMin1);

		long encodedProb;
		switch (totalBytes) {
			case 5:
				//cast to long to make bit shift work on 5 byte value
				encodedProb
						= ((long) bytes[firstByte + 4] & LAST_BYTE_MASK[bitsFromLastByteMin1]) << (24 + bitsFromLastByteMin1 - firstBitInFirstByte)
						| (bytes[firstByte + 3] & 255) << (24 - firstBitInFirstByte)
						| (bytes[firstByte + 2] & 255) << (16 - firstBitInFirstByte)
						| (bytes[firstByte + 1] & 255) << (8 - firstBitInFirstByte)
						| (bytes[firstByte] & FIRST_BYTE_MASK[firstBitInFirstByte]);
				break;
			case 4:
				encodedProb
						= (bytes[firstByte + 3] & LAST_BYTE_MASK[bitsFromLastByteMin1]) << (16 + bitsFromLastByteMin1 - firstBitInFirstByte)
						| (bytes[firstByte + 2] & 255) << (16 - firstBitInFirstByte)
						| (bytes[firstByte + 1] & 255) << (8 - firstBitInFirstByte)
						| (bytes[firstByte] & FIRST_BYTE_MASK[firstBitInFirstByte]);
				break;
			case 3:
				encodedProb
						= (bytes[firstByte + 2] & LAST_BYTE_MASK[bitsFromLastByteMin1]) << (9 + bitsFromLastByteMin1 - firstBitInFirstByte)
						| (bytes[firstByte + 1] & 255) << (8 - firstBitInFirstByte)
						| (bytes[firstByte] & FIRST_BYTE_MASK[firstBitInFirstByte]);
				break;
			case 2:
				if (bitsFromLastByteMin1 < firstBitInFirstByte) {
					encodedProb
							= (bytes[firstByte + 1] & LAST_BYTE_MASK[bitsFromLastByteMin1]) >> (firstBitInFirstByte - bitsFromLastByteMin1)
							| (bytes[firstByte] & FIRST_BYTE_MASK[firstBitInFirstByte]);
				} else {
					encodedProb
							= (bytes[firstByte + 1] & LAST_BYTE_MASK[bitsFromLastByteMin1]) << (bitsFromLastByteMin1 - firstBitInFirstByte +1)
							| (bytes[firstByte] & FIRST_BYTE_MASK[firstBitInFirstByte]);
				}
				break;
			case 1:

				System.out.println(Long.toBinaryString(((long) (bytes[firstByte])
						& LAST_BYTE_MASK[bitsFromLastByteMin1 + firstBitInFirstByte]
						& FIRST_BYTE_MASK[firstBitInFirstByte])));

				encodedProb
						= (bytes[firstByte]
						& LAST_BYTE_MASK[bitsFromLastByteMin1 + firstBitInFirstByte]
						& FIRST_BYTE_MASK[firstBitInFirstByte])
						>> (7 - firstBitInFirstByte - bitsFromLastByteMin1);
				break;
			default:
				throw new GenotypeDataException("Error parsing bgen file");
		}

		System.out.println(encodedProb);
		System.out.println(Long.toBinaryString(encodedProb));

		return 0d;

	}

}
