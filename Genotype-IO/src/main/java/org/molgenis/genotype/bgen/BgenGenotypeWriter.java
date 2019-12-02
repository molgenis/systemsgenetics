/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.bgen;

import java.io.*;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.channels.Channels;
import java.nio.channels.WritableByteChannel;
import java.nio.charset.Charset;
import java.nio.charset.StandardCharsets;
import java.util.*;

import com.github.luben.zstd.Zstd;
import org.apache.log4j.Logger;
import org.molgenis.genotype.*;
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
	private final double maxValue;
	private final int probabilitiesLengthInBits;
    private final double maxValue32Bits = Math.pow(2, 32) - 1;
	private final double maxValue16Bits = Math.pow(2, 16) - 1;
	private Zstd zstd = new Zstd();

	public BgenGenotypeWriter(GenotypeData genotypeData) {
		this.genotypeData = genotypeData;
		this.probabilitiesLengthInBits = 16;
		this.maxValue = Math.pow(2, probabilitiesLengthInBits) - 1;
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

	/**
	 * Method for writing genotype data to a BGEN file.
	 *
	 * @param bgenFile the File to write to.
	 * @return missingness of genotype data per individual
	 * @throws IOException if an I/O error has occurred.
	 */
	private HashMap<Sample, Float> writeBgenFile(File bgenFile) throws IOException {
        // Initialize an output stream for the .bgen file.
		OutputStream bgenOutputStream = new BufferedOutputStream(new FileOutputStream(bgenFile));

		// A ByteBuffer object holds bytes, and has methods to 'put' Integers, Shorts etc.
        // The ByteBuffer can be set to put these values in little endian byte order or big endian byte order

		// Because we would like to write ByteBuffer objects to the outputfile,
        // we create a new channel to which these byte buffers can be written.
		WritableByteChannel bgenOutputByteChannel = Channels.newChannel(bgenOutputStream);

		// Start writing the file.
        // We do not now the offset values, so we first gather data...

		// Get the number of samples and the number of variants.
		int sampleCount = genotypeData.getSamples().size();
		int variantCount = genotypeData.getVariantAnnotations().size();

		// Calculate the offset, relative to the fifth byte of the file,
		// of the first byte of the first variant data block
		int freeDataAreaLength = 0;
		int minimumHeaderLength = 20;
		int headerBlockLength = minimumHeaderLength + freeDataAreaLength;

		// Initialize a first bytebuffer for the entire header block.
		ByteBuffer headerBytesBuffer = ByteBuffer.allocate(headerBlockLength)
                .order(ByteOrder.LITTLE_ENDIAN);

		// Write the header block.
		// 4 bytes with length of the header block.
		headerBytesBuffer.putInt(headerBlockLength);
		// 4 bytes with the number of variant data blocks stored in the file.
		headerBytesBuffer.putInt(variantCount);
		// 4 bytes indicating the number of samples represented in the variant data blocks in the file.
		headerBytesBuffer.putInt(sampleCount);
		// Free data area, leave empty or store with description from genotype data
		// Write flags in 4 bytes

        // Generate the first byte containing compression value (2, in bits 0-1)
        // and the layout (2, in bits 2-5)
		// Always writing BGEN files with layout version 2
		byte byte_temp = (byte) 8; // 2 << 2 & 28 = 00001000 = 8
		// Use compression type 2
		byte_temp = (byte) (byte_temp | 2);

		// Put this temporary byte in the buffer, and add 2 empty ones
		headerBytesBuffer.put(byte_temp);
		headerBytesBuffer.put((byte) 0);
		headerBytesBuffer.put((byte) 0);

        // Write sample identifier block always; this is recommended

        // Set sample identifier flag (1, on bit 31, last one (zero based index))
		byte byte_temp_two = (byte) 128; // 1 << 7 & 128 = 10000000 = 128
		headerBytesBuffer.put(byte_temp_two);

		// Start with the sample identifier block.

        // 4 bytes with the length of the sample identifier block
        // Initialize as a long, as this stores up to 32^2-1
		long sampleIdentifierBlockLength = 8; // The length is at least 8 bytes.

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

		// Get the total offset
        long offset = headerBlockLength + sampleIdentifierBlockLength;

        // Check if the offset from the start of the file is greater than the max that fits in 32 bits.
        // If the offset fits within the 32 bits, the sample identifier block length will also fit within its
        // respective field.
        if (offset > maxValue32Bits) {
            throw new GenotypeDataException("The length of the sample identifier block " +
                    "is too large for the BGEN file format");
        }
        // We check for the 32bits instead of the maximum value that the signed integer
        // supports because casting it will not change the individual bits (still 32 bits are used.).


		// Create a byte buffer for the first 8 bytes within the sample identifier block.
        ByteBuffer sampleIdentifierBlockHeader = ByteBuffer.allocate(8)
                .order(ByteOrder.LITTLE_ENDIAN);
        // Fill it.
		sampleIdentifierBlockHeader.putInt((int) sampleIdentifierBlockLength);
		sampleIdentifierBlockHeader.putInt(sampleCount);

        // Create a byte buffer for the first four bytes of the entire file and fill it
        ByteBuffer firstFourBytesBuffer = ByteBuffer.allocate(4)
				.order(ByteOrder.LITTLE_ENDIAN);
        firstFourBytesBuffer.putInt((int) offset);

        // Write the first sections of the Bgen file up to the variant sections.
		// Write the first four bytes to the output channel.
		bgenOutputByteChannel.write(firstFourBytesBuffer);

		// Write header to bgenOutputByteChannel
		bgenOutputByteChannel.write(headerBytesBuffer);

		// Write sample block output stream to the bgen output stream
		sampleBlockOutputStream.writeTo(bgenOutputStream);

		// Start writing the variants sections.

		// Initialize missingness array. Implementation mimics that of the gen genotype writer.
		float[] sampleMissingCount = new float[sampleCount];

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

			// Write the number of alleles in 16 bits.
			int alleleCount = variant.getAlleleCount();
			if (alleleCount > maxValue16Bits) {
				throw new GenotypeDataException(String.format("Error, number of alleles for variant %s too large", variant));
			}
			// Can now safely cast to short.
			variantBuffer.putShort((short) alleleCount);
			bgenOutputByteChannel.write(variantBuffer);
			// Write alleles
			for (Allele allele : variant.getVariantAlleles()) {
				// Write the allele
				bgenOutputByteChannel.write(getFieldLengthFieldBuffer(
						allele.getAlleleAsString(), 4, "allele"));
			}

			// Get the genotype data
            ByteBuffer genotypeDataBlockByteBuffer = getGenotypeDataBlock(sampleCount,
                    sampleMissingCount, variant, alleleCount);

            // Write the compressed genotype data to the output channel.
			writeCompressedBgenGenotypeDataBlock(bgenOutputByteChannel, genotypeDataBlockByteBuffer);
		}

		// Return the missingness of the samples.
		HashMap<Sample, Float> sampleMissingness = new HashMap<>();
		for (int i = 0; i < sampleMissingCount.length; ++i) {
			sampleMissingness.put(genotypeData.getSamples().get(i), sampleMissingCount[i] / (float) variantCount);
		}
		return sampleMissingness;
	}

    private ByteBuffer getGenotypeDataBlock(int sampleCount, float[] sampleMissingCount, GeneticVariant variant, int alleleCount) throws IOException {
        // First declare an empty ByteBuffer.
        ByteBuffer genotypeDataBlockByteBuffer;

        // Check if phased data is available for all samples
        if (variant.getSamplePhasing().contains(false)) {
            // Get the phased genotype data block byte buffer if phased data is available.
            genotypeDataBlockByteBuffer = getPhasedGenotypeDataBlockByteBuffer(
                    sampleCount, sampleMissingCount, variant, alleleCount);
        } else {
            // Get the unphased genotype data block byte buffer if phased data is available.
            genotypeDataBlockByteBuffer = getUnphasedGenotypeDataBlockByteBuffer(
                    sampleCount, sampleMissingCount, variant, alleleCount);
        }
        return genotypeDataBlockByteBuffer;
    }

    private ByteBuffer getUnphasedGenotypeDataBlockByteBuffer(int sampleCount,
															  float[] sampleMissingCount,
															  GeneticVariant variant,
															  int alleleCount) throws IOException {
		double[][] sampleGenotypeProbabilitiesBgen = variant.getSampleGenotypeProbabilitiesBgen();

		// Init the minimum ploidy with the max possible value, can only get lower.
		int minimumPloidy = 63;
		// Init the maximum ploidy with the min possible value, can only get higher.
		int maximumPloidy = 0;
		int totalNumberOfProbabilities = 0;
		byte[] ploidyMissingnessBytes = new byte[sampleCount];

		for (int i = 0; i < sampleCount; i++) {
			double[] sampleProbabilities = sampleGenotypeProbabilitiesBgen[i];
			int ploidy = sampleProbabilities.length;
			if (0 >= ploidy || ploidy > 63) {
				throw new GenotypeDataException("Ploidy is not between 0 and 63.");
			}
			if (ploidy < minimumPloidy) {
				minimumPloidy = ploidy;
			}
			if (ploidy > maximumPloidy) {
				maximumPloidy = ploidy;
			}
			boolean missingness = Arrays.stream(sampleProbabilities).sum() == 0;
			sampleMissingCount[i]++;

			// Put the ploidy missingness byte within the buffer.
			ploidyMissingnessBytes[i] = getPloidyMissingnessByte(ploidy, missingness);
			// Add the number of possible probabilities.
			totalNumberOfProbabilities += ploidy * (alleleCount - 1);
		}

		int probabilitiesBytesCount = totalNumberOfProbabilities * (probabilitiesLengthInBits / 8);

		ByteBuffer genotypeDataBlockBuffer = initializeGenotypeDataBlockByteBuffer(
				sampleCount, alleleCount,
				minimumPloidy, maximumPloidy,
				ploidyMissingnessBytes, probabilitiesBytesCount, false);

		int bitOffset = 0;
		byte[] probabilitiesByteArray = new byte[probabilitiesBytesCount];
		for (double[] sampleProbabilities : sampleGenotypeProbabilitiesBgen) {
			writeProbabilitiesSummingToOne(probabilitiesByteArray, bitOffset, sampleProbabilities);
			bitOffset += probabilitiesLengthInBits * (alleleCount - 1);
		}

		genotypeDataBlockBuffer.put(probabilitiesByteArray);
		return genotypeDataBlockBuffer;
	}

	private ByteBuffer initializeGenotypeDataBlockByteBuffer(int sampleCount, int alleleCount,
															 int minimumPloidy, int maximumPloidy,
															 byte[] ploidyMissingnessBytes, int probabilitiesBytesCount,
															 boolean isPhased) {

		ByteBuffer genotypeDataBlockBuffer = ByteBuffer
				.allocate(10 + sampleCount + probabilitiesBytesCount)
				.order(ByteOrder.LITTLE_ENDIAN);

		// The number of individuals
		genotypeDataBlockBuffer.putInt(sampleCount);
		genotypeDataBlockBuffer.putShort((short) alleleCount);

		// Put min ploidy
		genotypeDataBlockBuffer.put((byte) minimumPloidy);
		// Put max ploidy
		genotypeDataBlockBuffer.put((byte) maximumPloidy);
		// Put ploidy & missingness byte
		genotypeDataBlockBuffer.put(ploidyMissingnessBytes);

		// Phased is false...
		genotypeDataBlockBuffer.put((byte) (isPhased ? 1 : 0));
		// Store write probabilities lengths in bits.
		genotypeDataBlockBuffer.put((byte) probabilitiesLengthInBits);
		return genotypeDataBlockBuffer;
	}

	private void writeCompressedBgenGenotypeDataBlock(WritableByteChannel bgenOutputByteChannel, ByteBuffer genotypeDataBlockBuffer) throws IOException {

		ByteBuffer compressedGenotypeDataBuffer = Zstd.compress(genotypeDataBlockBuffer, 3);

		int lengthOfVariantProbabilityData =
				compressedGenotypeDataBuffer.capacity() + 4; // Add 4 since the D field takes 4 bytes
		int decompressedLengthOfVariantProbabilityData = genotypeDataBlockBuffer.capacity();

		ByteBuffer genotypeBlockLengths = ByteBuffer.allocate(8)
				.order(ByteOrder.LITTLE_ENDIAN);

		genotypeBlockLengths.putInt(lengthOfVariantProbabilityData);
		genotypeBlockLengths.putInt(decompressedLengthOfVariantProbabilityData);

		bgenOutputByteChannel.write(genotypeBlockLengths);
		bgenOutputByteChannel.write(compressedGenotypeDataBuffer);
	}

	private ByteBuffer getPhasedGenotypeDataBlockByteBuffer(int sampleCount,
															float[] sampleMissingCount,
															GeneticVariant variant,
															int alleleCount) {
		double[][][] sampleGenotypeProbabilitiesBgenPhased = variant.getSampleGenotypeProbabilitiesBgenPhased();

		// Init the minimum ploidy with the max possible value, can only get lower.
		int minimumPloidy = 63;
		// Init the maximum ploidy with the min possible value, can only get higher.
		int maximumPloidy = 0;
		int totalNumberOfProbabilities = 0;
		byte[] ploidyMissingnessBytes = new byte[sampleCount];

		for (int i = 0; i < sampleCount; i++) {
			double[][] sampleProbabilities = sampleGenotypeProbabilitiesBgenPhased[i];
			int ploidy = sampleProbabilities.length;
			if (0 >= ploidy || ploidy > 63) {
				throw new GenotypeDataException("Ploidy is not between 0 and 63.");
			}
			if (ploidy < minimumPloidy) {
				minimumPloidy = ploidy;
			}
			if (ploidy > maximumPloidy) {
				maximumPloidy = ploidy;
			}
			boolean missingness = Arrays.stream(sampleProbabilities)
					.flatMapToDouble(Arrays::stream).sum() == 0;
			sampleMissingCount[i]++;

			// Put the ploidy missingness byte within the buffer.
			ploidyMissingnessBytes[i] = getPloidyMissingnessByte(ploidy, missingness);
			// Add the number of possible probabilities.
			totalNumberOfProbabilities += ploidy * (alleleCount - 1);
		}

		int probabilitiesBytesCount = totalNumberOfProbabilities * (probabilitiesLengthInBits / 8);

		ByteBuffer genotypeDataBlockBuffer = initializeGenotypeDataBlockByteBuffer(
				sampleCount, alleleCount,
				minimumPloidy, maximumPloidy,
				ploidyMissingnessBytes, probabilitiesBytesCount, true); // data is phased, so true here.

		int bitOffset = 0; // Start with an offset of zero.
		byte[] probabilitiesByteArray = new byte[totalNumberOfProbabilities * (probabilitiesLengthInBits / 8)];
		for (double[][] sampleProbabilities : sampleGenotypeProbabilitiesBgenPhased) {
			for (double[] alleleProbabilities : sampleProbabilities){
				writeProbabilitiesSummingToOne(probabilitiesByteArray, bitOffset, alleleProbabilities);
				bitOffset += probabilitiesLengthInBits * (alleleCount - 1);
			}
		}

		genotypeDataBlockBuffer.put(probabilitiesByteArray);
		return genotypeDataBlockBuffer;
	}

	private byte getPloidyMissingnessByte(int ploidy, boolean missingness) {
		return (byte) ((byte) (ploidy & 63) | (missingness ? 128 : 0));

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

	/**
	 * Method that is responsible for writing an array of probabilities that should sum to one.
	 * The last probability is not included since this can be inferred from the other values.
	 *
	 * The implementation is based on the following c++ reference implementation retrieved on 19-12-02.
	 * https://bitbucket.org/gavinband/bgen/src/default/src/bgen.cpp
	 *
	 * @param bytes A byte array to write the probabilities to.
	 * @param bitOffset The offset within the byte array to write the probabilities from.
	 * @param probabilities The array of probabilities to write.
	 */
	private void writeProbabilitiesSummingToOne(byte[] bytes, int bitOffset, double[] probabilities) {

		// Find sum to renormalize to 1
		double sum = Arrays.stream(probabilities).sum();

		// Calculate the probabilities to store.
		List<Integer> indices = new ArrayList<>();
		double totalFractionalPart = 0;

		// Multiply the probabilities by the maximum value,
		// and calculate the total fractional part of the multiplied probabilities
		for (int i = 0; i < probabilities.length; i++) {
			probabilities[i] /= sum;
			probabilities[i] *= maxValue;
			indices.add(i);
			totalFractionalPart += fractionalPart(probabilities[i]);
		}

		// Round the total fractional part
		long roundedSumOfFractionalPart = (int) Math.round(totalFractionalPart);
		// Make sure things do not proceed if the rounded sum of fractional part is greater than the maximum of an int.
		assert(roundedSumOfFractionalPart < Integer.MAX_VALUE);

		//
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
