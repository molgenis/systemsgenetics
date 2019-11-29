/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.bgen;

import java.io.*;
import java.math.BigInteger;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.channels.Channels;
import java.nio.channels.WritableByteChannel;
import java.nio.charset.Charset;
import java.nio.charset.StandardCharsets;
import java.util.*;

import org.apache.log4j.Logger;
import org.molgenis.genotype.Allele;
import org.molgenis.genotype.GenotypeData;
import org.molgenis.genotype.GenotypeWriter;
import org.molgenis.genotype.Sample;
import org.molgenis.genotype.oxford.OxfordSampleFileWriter;
import org.molgenis.genotype.util.Utils;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variant.NotASnpException;

/**
 *
 * @author patri
 */
public class BgenGenotypeWriter implements GenotypeWriter {

	/**
	 * Index is the number of bits used from the last byte. (these are the rightmost bits)
	 */
	private static final int[] LAST_BYTE_MASK = {255, 1, 3, 7, 15, 31, 63, 127, 255};

	/**
	 * Index is the number of bits used from the first byte. (these are the leftmost bits)
	 */
	private static final int[] FIRST_BYTE_MASK = {0, 128, 192, 224, 240, 248, 252, 254, 255};

	private static final Logger LOGGER = Logger.getLogger(BgenGenotypeWriter.class);
	private static final Charset CHARSET = StandardCharsets.UTF_8;
	private final GenotypeData genotypeData;

	public BgenGenotypeWriter(GenotypeData genotypeData) {
		this.genotypeData = genotypeData;
	}

	@Override
	public void write(String path) throws IOException, NotASnpException {
		throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
	}

	public void write(File bgenFile, File sampleFile) throws IOException {
		LOGGER.info("Writing BGEN file " + bgenFile.getAbsolutePath() + " and sample file "
				+ sampleFile.getAbsolutePath());

		Utils.createEmptyFile(bgenFile, "gen");
		Utils.createEmptyFile(sampleFile, "sample");

		HashMap<Sample, Float> sampleMissingness = writeBgenFile(bgenFile);
		OxfordSampleFileWriter.writeSampleFile(sampleFile, genotypeData, sampleMissingness);
	}

	private HashMap<Sample, Float> writeBgenFile(File bgenFile) throws IOException {

		OutputStream bgenOutputStream = new BufferedOutputStream(new FileOutputStream(bgenFile));
		WritableByteChannel bgenOutputByteChannel = Channels.newChannel(bgenOutputStream);
		// Number of samples
		int sampleCount = genotypeData.getSamples().size();
		int variantCount = genotypeData.getVariantAnnotations().size();

		// Calculate the offset, relative to the fifth byte of the file,
		// of the first byte of the first variant data block
		long freeDataAreaLength = 0;
		int minimumHeaderLength = 20;
		long headerBlockLength = minimumHeaderLength + freeDataAreaLength;

		ByteBuffer headerBytesBuffer = ByteBuffer.allocate(20);
		headerBytesBuffer.order(ByteOrder.LITTLE_ENDIAN);

		// Write the header block.
		// 4 bytes with length of the header block.
		headerBytesBuffer.putInt((int) headerBlockLength);
		// 4 bytes with the number of variant data blocks stored in the file.
		headerBytesBuffer.putInt(variantCount);
		// 4 bytes indicating the number of samples represented in the variant data blocks in the file.
		headerBytesBuffer.putInt(sampleCount);
		// Free data area, leave empty or store with description from genotype data
		// Write flags in 4 bytes

		// Always writing BGEN files with layout version 2
		byte byte_temp = (byte) 8; // 2 << 2 & 28 = 00001000 = 8
		// Use compression type 2
		byte_temp = (byte) (byte_temp | 2);

		// Put bytes
		headerBytesBuffer.put(byte_temp);
		headerBytesBuffer.put((byte) 0);
		headerBytesBuffer.put((byte) 0);

		// Set sample identifier flag
		byte byte_temp_two = (byte) 128; // 1 << 7 & 128 = 10000000 = 128
		headerBytesBuffer.put(byte_temp_two);

		// Write sample identifier block (always; this is recommended)
		ByteBuffer sampleIdentifierBlockHeader = ByteBuffer.allocate(8)
				.order(ByteOrder.LITTLE_ENDIAN);

		// 4 bytes with the length of the sample identifier block
		long sampleIdentifierBlockLength = 8;

		// 4 bytes with the number of samples represented in the file.

		ByteArrayOutputStream sampleBlockOutputStream = new ByteArrayOutputStream();
		WritableByteChannel sampleBlockByteChannel = Channels.newChannel(sampleBlockOutputStream);

		// Write the sample identifiers for every sample
		for (Sample sample : genotypeData.getSamples()) {
			// Get the identifier and its length
			String id = sample.getId();
			ByteBuffer sampleBuffer = getFieldLengthFieldBuffer(id, 2, "sample identifier");
			// Write the buffer to the channel
			sampleIdentifierBlockLength += sampleBlockByteChannel.write(sampleBuffer);
		}

		long offset = headerBlockLength + sampleIdentifierBlockLength;
		ByteBuffer firstFourBytesBuffer = ByteBuffer.allocate(4).order(ByteOrder.LITTLE_ENDIAN);
		firstFourBytesBuffer.putInt((int) offset);

		sampleIdentifierBlockHeader.putInt((int) sampleIdentifierBlockLength);
		sampleIdentifierBlockHeader.putInt(sampleCount);

		// Write the first four bytes to the output channel.
		bgenOutputByteChannel.write(firstFourBytesBuffer);

		// Write header to bgenOutputByteChannel
		bgenOutputByteChannel.write(headerBytesBuffer);

		// Write sample block output stream to the bgen output stream
		sampleBlockOutputStream.writeTo(bgenOutputStream);

		// Loop through variants
		for (GeneticVariant variant : genotypeData) {
			// Write variant data block
			// Start with the length of the variant identifier
			List<String> allIds = variant.getAllIds();
			String alternativeId = allIds.size() > 1 ? allIds.get(1) : "";
			// Write variant identifier
			bgenOutputByteChannel.write(getFieldLengthFieldBuffer(
					alternativeId, 2, "variant identifier"));
			// Write the RSID
			bgenOutputByteChannel.write(getFieldLengthFieldBuffer(
					variant.getPrimaryVariantId(), 2, "rs identifier"));
			// Write the chromosome
			bgenOutputByteChannel.write(getFieldLengthFieldBuffer(
					variant.getSequenceName(), 2, "chromosome"));

			// Write the variant position and the number of alleles
			ByteBuffer variantBuffer = ByteBuffer.allocate(6).order(ByteOrder.LITTLE_ENDIAN);
			variantBuffer.putInt(variant.getStartPos());
			int alleleCount = variant.getAlleleCount();
			if (alleleCount > Short.MAX_VALUE) {
				// do something
			}
			variantBuffer.putShort((short) alleleCount);
			bgenOutputByteChannel.write(variantBuffer);
			// Write alleles
			for (Allele allele : variant.getVariantAlleles()) {
				// Write the allele
				bgenOutputByteChannel.write(getFieldLengthFieldBuffer(
						allele.getAlleleAsString(), 4, "allele"));
			}
			// Genotype data block
			int lengthOfVariantProbabilityData = 0;
			int decompressedLengthOfVariantProbabilityData = 0;
			// Gather probability data storage
			if (variant.getSamplePhasing().contains(false)) {
				variant.getSampleGenotypeProbabilitiesBgenPhased();
			} else {
				variant.getSampleGenotypeProbabilitiesBgen();
			}

			ByteBuffer buffer = ByteBuffer.allocate(8);
			// The number of individuals
			buffer.putInt(sampleCount);
			buffer.putShort((short) alleleCount);
			// Put min ploidy
			// Put max ploidy
		}
		return null;
	}

	private ByteBuffer getFieldLengthFieldBuffer(String fieldValue, int maxFieldBytes, String fieldName) {
		int fieldLength = fieldValue.length();
		byte[] idAsBytes = fieldValue.getBytes(CHARSET);

		// 2 bytes with the length of identifier of sample 1
		if (idAsBytes.length > maxFieldBytes) {
			// The length of the sample identifier exceeds two bytes. Give a warning and truncate.
			LOGGER.warn(String.format(
					"Length of %s %.16s... exceeds the maximum length of 2 bytes. (%d vs %d respectively) " +
							"%nField will be truncated.",
					fieldName, fieldValue, fieldLength, Short.MAX_VALUE));
			// ...truncate...
			fieldLength = Short.MAX_VALUE;
			idAsBytes = Arrays.copyOfRange(idAsBytes, 0, maxFieldBytes + 1);
		}
		// Initialize a byte buffer with the correct size.
		ByteBuffer sampleBuffer = ByteBuffer.allocate(maxFieldBytes + idAsBytes.length)
				.order(ByteOrder.LITTLE_ENDIAN);
		// Put both the length of the sample identifier as a short
		// and the identifier itself int the byte array
		sampleBuffer.putShort((short) fieldLength);
		sampleBuffer.put(idAsBytes);
		return sampleBuffer;
	}

	private static double fractionalPart(double v) {
		return v - Math.floor(v);
	}

	private static void writeProbabilities(byte[] bytes, int bitOffset, double[] probabilities, int probabilitiesLengthInBits) {
		double maxValue = Math.pow(2, probabilitiesLengthInBits) - 1;
		List<Integer> indices = new ArrayList<>();
		double totalFractionalPart = 0;

		for (int i = 0; i < probabilities.length; i++) {
			probabilities[i] *= maxValue;
			indices.add(i);
			totalFractionalPart += fractionalPart(probabilities[i]);
		}

		long roundedSumOfFractionalPart = (int) Math.round(totalFractionalPart);
		assert(roundedSumOfFractionalPart < Integer.MAX_VALUE);

		indices.sort(Comparator.comparingDouble(a -> fractionalPart(probabilities[a])));

		for (int i = 0; i < roundedSumOfFractionalPart; i++) {
			probabilities[indices.get(i)] = Math.ceil(probabilities[indices.get(i)]);
		}
		for (int i = (int) roundedSumOfFractionalPart;
			 i < probabilities.length; i++) {
			probabilities[indices.get(i)] = Math.floor(probabilities[indices.get(i)]);
		}

		for (int i = 0; i < (probabilities.length - 1); i++) {
			probabilityValueToByteArray((long) probabilities[i], bytes, bitOffset, probabilitiesLengthInBits);
			bitOffset += probabilitiesLengthInBits;
		}
	}

	/**
	 * Writes n bits (with a maximum of 32) from a long value, to a byte array starting from the offset value.
	 *
	 * @param value	The probability represented by a long to convert to a byte array.
	 * @param bytes				The byte array to write the probability value to.
	 * @param bitOffset       	The bit in the byte array to start writing from.
	 * @param numberOfBitsToWrite 	The total number of bits to write to represent the probability.
	 */
	private static void probabilityValueToByteArray(long value, byte[] bytes,
													int bitOffset, int numberOfBitsToWrite) {
		// Get the byte to start reading from.
		int firstByteIndex = bitOffset / 8;
		// Get the bit within the byte to start reading from.
		int remainingBitOffset = bitOffset % 8;

		int totalBytesMin1 = (remainingBitOffset + numberOfBitsToWrite - 1) / 8;

		int bitShiftAfterFirstByte = 8 - remainingBitOffset;

		switch (totalBytesMin1) {
			case 4:
				bytes[firstByteIndex + 4] |= (value >>> (24 + bitShiftAfterFirstByte));
			case 3:
				bytes[firstByteIndex + 3] |= (value >>> (16 + bitShiftAfterFirstByte));
			case 2:
				bytes[firstByteIndex + 2] |= (value >>> (8 + bitShiftAfterFirstByte));
			case 1:
				bytes[firstByteIndex + 1] |= (value >>> bitShiftAfterFirstByte);
			case 0:
				bytes[firstByteIndex] |= (value << remainingBitOffset);
		}
	}

}
