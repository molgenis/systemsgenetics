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
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import com.github.luben.zstd.Zstd;
import com.google.common.collect.Iterators;
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

	private static final Logger LOGGER = Logger.getLogger(BgenGenotypeWriter.class);
	private static final Charset CHARSET = StandardCharsets.UTF_8;
	private final GenotypeData genotypeData;
	private final double maxValue;
	private int probabilitiesLengthInBits;
    private final double maxValue32Bits = Math.pow(2, 32) - 1;
	private final double maxValue16Bits = Math.pow(2, 16) - 1;
	private Zstd zstd = new Zstd();

	public BgenGenotypeWriter(GenotypeData genotypeData) {
		this.genotypeData = genotypeData;
		this.probabilitiesLengthInBits = 16;
		this.maxValue = Math.pow(2, probabilitiesLengthInBits) - 1;
	}

	@Override
	public void write(String basePath) throws IOException, NotASnpException {

		if (!genotypeData.isOnlyContaingSaveProbabilityGenotypes()) {
			LOGGER.warn("WARNING!!! writing dosage genotype data to .gen posterior probabilities file. Using heuristic method to convert to probabilities, this is not guaranteed to be accurate. See manual for more details.");
		}

		write(new File(basePath + ".bgen"), new File(basePath + ".sample"));
	}

	public void write(File bgenFile, File sampleFile) throws IOException {
		LOGGER.info("Writing BGEN file " + bgenFile.getAbsolutePath() + " and sample file "
				+ sampleFile.getAbsolutePath());

		Utils.createEmptyFile(bgenFile, "bgen");
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
		// Initialize a bgenix writer.
		File bgenixFile = new File(bgenFile.getAbsolutePath() + ".bgi");
		BgenixWriter bgenixWriter = new BgenixWriter(bgenixFile);


		// A ByteBuffer object holds bytes, and has methods to 'put' Integers, Shorts etc.
        // The ByteBuffer can be set to put these values in little endian byte order or big endian byte order

		// Because we would like to write ByteBuffer objects to the outputfile,
        // we create a new channel to which these byte buffers can be written.
		WritableByteChannel bgenOutputByteChannel = Channels.newChannel(bgenOutputStream);

		// Start writing the file.
        // We do not now the offset values, so we first gather data...

		// Get the number of samples and the number of variants.
		int sampleCount = genotypeData.getSamples().size();
		int variantCount = Iterators.size(genotypeData.iterator());

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

		// Put magic bytes
		headerBytesBuffer.put("b".getBytes(CHARSET));
		headerBytesBuffer.put("g".getBytes(CHARSET));
		headerBytesBuffer.put("e".getBytes(CHARSET));
		headerBytesBuffer.put("n".getBytes(CHARSET));
		// Put this temporary byte in the buffer, and add 2 empty ones
		headerBytesBuffer.put(byte_temp);
		headerBytesBuffer.put((byte) 0);
		headerBytesBuffer.put((byte) 0);

        // Write sample identifier block always; this is recommended

        // Set sample identifier flag (1, on bit 31, last one (zero based index))
		byte byte_temp_two = (byte) 128; // 1 << 7 & 128 = 10000000 = 128
		headerBytesBuffer.put(byte_temp_two);

		headerBytesBuffer.flip(); // reset pointer

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
			sampleIdentifierBlockLength += writeFieldWithFieldLength(sampleBlockByteChannel,
					id, 2, "sample identifier");
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

		sampleIdentifierBlockHeader.flip(); // reset pointer

        // Create a byte buffer for the first four bytes of the entire file and fill it
        ByteBuffer firstFourBytesBuffer = ByteBuffer.allocate(4)
				.order(ByteOrder.LITTLE_ENDIAN);
        firstFourBytesBuffer.putInt((int) offset);
        firstFourBytesBuffer.flip(); // reset pointer

        // Write the first sections of the Bgen file up to the variant sections.
		// Write the first four bytes to the output channel.
		bgenOutputByteChannel.write(firstFourBytesBuffer);

		// Write header to bgenOutputByteChannel
		bgenOutputByteChannel.write(headerBytesBuffer);

		// Write the sample identifier block header
		bgenOutputByteChannel.write(sampleIdentifierBlockHeader);

		// Write sample block output stream to the bgen output stream
		sampleBlockOutputStream.writeTo(bgenOutputStream);

		// Start writing the variants sections.

		// Initialize missingness array. Implementation mimics that of the gen genotype writer.
		float[] sampleMissingCount = new float[sampleCount];

		long variantStartPositionInFile = offset + firstFourBytesBuffer.limit();

		// Loop through variants
		for (GeneticVariant variant : genotypeData) {
			// Write variant data block
			// Initialize counter for the size of the variant data.
			long variantDataSizeInBytes = 0;

			// Start with the length of the variant identifier
			List<String> allIds = variant.getAllIds();
			String alternativeId = allIds.size() > 1 ? allIds.get(1) : "";
			// Write variant identifier
			variantDataSizeInBytes += writeFieldWithFieldLength(
					bgenOutputByteChannel, alternativeId, 2, "variant identifier");
			// Write the RSID
			variantDataSizeInBytes += writeFieldWithFieldLength(
					bgenOutputByteChannel, variant.getPrimaryVariantId(), 2, "rs identifier");
			// Write the chromosome
			variantDataSizeInBytes += writeFieldWithFieldLength(
					bgenOutputByteChannel, variant.getSequenceName(), 2, "chromosome");

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
			variantBuffer.flip(); // reset pointer
			variantDataSizeInBytes += bgenOutputByteChannel.write(variantBuffer);
			// Write alleles
			for (Allele allele : variant.getVariantAlleles()) {
				// Write the allele
				variantDataSizeInBytes += writeFieldWithFieldLength(bgenOutputByteChannel,
						allele.getAlleleAsString(), 4, "allele");
			}

			// Get the genotype data
            ByteBuffer genotypeDataBlockByteBuffer = getGenotypeDataBlock(sampleCount,
                    sampleMissingCount, variant);

            // Write the compressed genotype data to the output channel.
			variantDataSizeInBytes += writeCompressedBgenGenotypeDataBlock(
					bgenOutputByteChannel, genotypeDataBlockByteBuffer);

			LOGGER.debug(String.format("Written %s, %s at %d, of size %d | seq:pos = %s:%d, %d alleles",
					variant.getPrimaryVariantId(),
					!variant.getAlternativeVariantIds().isEmpty() ? variant.getAlternativeVariantIds().get(0) : "-",
					variantStartPositionInFile, variantDataSizeInBytes,
					variant.getSequenceName(), variant.getStartPos(), alleleCount));

			// Add the read variant to the BGENIX file so that it can quickly be retrieved.
			bgenixWriter.addVariantToIndex(
					variant,
					variantStartPositionInFile,
					variantDataSizeInBytes,
					variant.getPrimaryVariantId());

			variantStartPositionInFile += variantDataSizeInBytes;
		}

		bgenOutputStream.close();

		// Finalize bgenix file
		addMetaData(bgenFile, bgenixWriter);
		bgenixWriter.finalizeIndex();

		// Return the missingness of the samples.
		HashMap<Sample, Float> sampleMissingness = new HashMap<>();
		for (int i = 0; i < sampleMissingCount.length; ++i) {
			sampleMissingness.put(genotypeData.getSamples().get(i), sampleMissingCount[i] / (float) variantCount);
		}
		return sampleMissingness;
	}

	private void addMetaData(File bgenFile, BgenixWriter bgenixWriter) throws IOException {
		// Go to the first byte...
		RandomAccessFile randomAccessBgenFile = new RandomAccessFile(bgenFile, "r");
		randomAccessBgenFile.seek(0);

		// Read the first 1000 bytes of the bgen file.
		byte[] firstBytes = new byte[1000];
		randomAccessBgenFile.read(firstBytes, 0, 1000);

		//Add current time in int.
		System.out.println((System.currentTimeMillis() / 1000L));
		// Create and write new metadata.
		BgenixMetadata m = new BgenixMetadata(
				bgenFile.getName(),
				(int) bgenFile.length(),
				(int) (bgenFile.lastModified() / 1000L),
				firstBytes,
				(System.currentTimeMillis() / 1000L));
		bgenixWriter.writeMetadata(m);
	}

	/**
	 * Method that gets the genotype data for a given variant and writes it to a byte buffer.
	 * The sample missing count is updated with +1 when a sample misses data for this variant.
	 *
	 * @param sampleCount The number of samples in the genotype data.
	 * @param sampleMissingCount An array of floats representing the missingness for every sample.
	 * @param variant The variant to get the probabilities from.
	 * @return A ByteBuffer containing the entire probability data storage for the given variant.
	 * @throws IOException if an I/O error has occurred.
	 */
    private ByteBuffer getGenotypeDataBlock(int sampleCount, float[] sampleMissingCount, GeneticVariant variant) throws IOException {
        // First declare an empty ByteBuffer.
        ByteBuffer genotypeDataBlockByteBuffer;

        // Check if phased data is available for all samples
        if (!variant.getSamplePhasing().contains(false)) {
            // Get the phased genotype data block byte buffer if phased data is available.
            genotypeDataBlockByteBuffer = getPhasedGenotypeDataBlockByteBuffer(
                    sampleCount, sampleMissingCount, variant);
        } else {
            // Get the unphased genotype data block byte buffer if phased data is available.
            genotypeDataBlockByteBuffer = getUnphasedGenotypeDataBlockByteBuffer(
                    sampleCount, sampleMissingCount, variant);
        }
        genotypeDataBlockByteBuffer.flip(); // reset pointer
        return genotypeDataBlockByteBuffer;
    }

	/**
	 * Method that gets the <i>unphased</i> genotype data for a given variant and writes it to a byte buffer.
	 * The sample missing count is updated with +1 when a sample misses data for this variant.
	 *
	 * @param sampleCount The number of samples in the genotype data.
	 * @param sampleMissingCount An array of floats representing the missingness for every sample.
	 * @param variant The variant to get the probabilities from.
	 * @return A ByteBuffer containing the entire probability data storage for the given variant.
	 */
    private ByteBuffer getUnphasedGenotypeDataBlockByteBuffer(int sampleCount,
															  float[] sampleMissingCount,
															  GeneticVariant variant) {

    	// Get the unphased bgen probabilities (this can represent polyploidity and multiallelic variants)
		double[][] sampleGenotypeProbabilitiesBgen = variant.getSampleGenotypeProbabilitiesComplex();

		// Get the allele count for the variant.
		int alleleCount = variant.getAlleleCount();

		// Init the minimum ploidy with the max possible value, can only get lower.
		int minimumPloidy = 63;
		// Init the maximum ploidy with the min possible value, can only get higher.
		int maximumPloidy = 0;
		int totalNumberOfProbabilities = 0;

		// Build the combined bytes for ploidy and missingness
		byte[] ploidyMissingnessBytes = new byte[sampleCount];

		for (int i = 0; i < sampleCount; i++) {
			double[] sampleProbabilities = sampleGenotypeProbabilitiesBgen[i];
			int numberOfProbabilities = sampleProbabilities.length;

			int ploidy = getPloidy(numberOfProbabilities, alleleCount);
			// Ploidy is stored in the least significant 6 bits and thus hase a max of 63,
			// Check this.
			if (0 >= ploidy || ploidy > 63) {
				throw new GenotypeDataException("Ploidy is not between 0 and 63.");
			}
			// Update minimum and maximum if conditions apply
			if (ploidy < minimumPloidy) {
				minimumPloidy = ploidy;
			}
			if (ploidy > maximumPloidy) {
				maximumPloidy = ploidy;
			}
			// Get the missingness and update the missingness count.
			boolean missingness = Arrays.stream(sampleProbabilities).sum() == 0;
			sampleMissingCount[i]++;

			// Put the ploidy missingness byte within the buffer.
			ploidyMissingnessBytes[i] = getPloidyMissingnessByte(ploidy, missingness);
			// Add the number of expected probabilities.
			totalNumberOfProbabilities += numberOfProbabilities - 1;
		}

		// Convert the total number of probabilities to the number of bytes
		int probabilitiesBytesCount = (totalNumberOfProbabilities * probabilitiesLengthInBits + 7) / 8;

		// Generate a bytebuffer with samplecount, allelecount etc allready entered, with
		// room for the sample probabilities
		ByteBuffer genotypeDataBlockBuffer = initializeGenotypeDataBlockByteBuffer(
				sampleCount, alleleCount,
				minimumPloidy, maximumPloidy,
				ploidyMissingnessBytes, probabilitiesBytesCount, false);

		// Now start writing the sample probabilities to a byte array.
		int bitOffset = 0; // Initialize at 0
		byte[] probabilitiesByteArray = new byte[probabilitiesBytesCount]; // Initialize the byte array
		for (double[] sampleProbabilities : sampleGenotypeProbabilitiesBgen) {
			// These probabilities should sum to one, add this to the byte buffer.
			// Also add the number of bits appended to the bit offset.
			bitOffset = addProbabilitiesSummingToOne(
					probabilitiesByteArray, bitOffset, sampleProbabilities);
		}

		genotypeDataBlockBuffer.put(probabilitiesByteArray);
		return genotypeDataBlockBuffer;
	}

	/**
	 * Method that initializes a genotype data block with the given values appended
	 * in the order and manner that corresponds to the BGEN specs
	 *
	 * @param sampleCount The number of samples that the BGEN file represents.
	 * @param alleleCount The number of alleles that the variant has.
	 * @param minimumPloidy The minimum ploidy of all samples for this variant.
	 * @param maximumPloidy The maximum ploidy of all samples for this variant.
	 * @param ploidyMissingnessBytes A byte Array of bytes holding missingness
	 *                                  and ploidy information for every sample.
	 * @param probabilitiesBytesCount The number of probabilities that is to be written for this variant
	 * @param isPhased True if the probabilities represent phased data.
	 * @return A byte buffer with the given values and with room for the probabilities.
	 */
	private ByteBuffer initializeGenotypeDataBlockByteBuffer(int sampleCount, int alleleCount,
															 int minimumPloidy, int maximumPloidy,
															 byte[] ploidyMissingnessBytes, int probabilitiesBytesCount,
															 boolean isPhased) {

    	// Initialize a genotype data block buffer so that all probability data for a variant fits.
		ByteBuffer genotypeDataBlockBuffer = ByteBuffer
				.allocateDirect(10 // Fixed value for every genotype block.
						+ sampleCount // Every sample has a single byte containing missingness and ploidy.
						+ probabilitiesBytesCount) // The number of probabilities in bytes.
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

	/**
	 * Method that gets the <i>phased</i> genotype data for a given variant and writes it to a byte buffer.
	 * The sample missing count is updated with +1 when a sample misses data for this variant.
	 *
	 * @param sampleCount The number of samples in the genotype data.
	 * @param sampleMissingCount An array of floats representing the missingness for every sample.
	 * @param variant The variant to get the probabilities from.
	 * @return A ByteBuffer containing the entire probability data storage for the given variant.
	 */
	private ByteBuffer getPhasedGenotypeDataBlockByteBuffer(int sampleCount,
															float[] sampleMissingCount,
															GeneticVariant variant) {
		double[][][] sampleGenotypeProbabilitiesBgenPhased = variant.getSampleGenotypeProbabilitiesPhased();

		// Get the allele count
		int alleleCount = variant.getAlleleCount();

		// Init the minimum ploidy with the max possible value, can only get lower.
		int minimumPloidy = 63;
		// Init the maximum ploidy with the min possible value, can only get higher.
		int maximumPloidy = 0;
		int totalNumberOfProbabilities = 0;
		byte[] ploidyMissingnessBytes = new byte[sampleCount];

		for (int i = 0; i < sampleCount; i++) {
			double[][] sampleProbabilities = sampleGenotypeProbabilitiesBgenPhased[i];
			int ploidy = sampleProbabilities.length;
			// Ploidy is stored in the least significant 6 bits and thus hase a max of 63,
			// Check this.
			if (0 >= ploidy || ploidy > 63) {
				throw new GenotypeDataException("Ploidy is not between 0 and 63.");
			}
			// Update minimum and maximum if conditions apply
			if (ploidy < minimumPloidy) {
				minimumPloidy = ploidy;
			}
			if (ploidy > maximumPloidy) {
				maximumPloidy = ploidy;
			}
			// Get the missingness and update the missingness count.
			boolean missingness = Arrays.stream(sampleProbabilities)
					.flatMapToDouble(Arrays::stream).sum() == 0;
			sampleMissingCount[i]++;

			// Put the ploidy missingness byte within the buffer.
			ploidyMissingnessBytes[i] = getPloidyMissingnessByte(ploidy, missingness);
			// Add the number of possible probabilities.
			totalNumberOfProbabilities += ploidy * (alleleCount - 1);
		}

		// Convert the total number of probabilities to the number of bytes
		int probabilitiesBytesCount = totalNumberOfProbabilities * probabilitiesLengthInBits / 8 + 1;

		// Generate a bytebuffer with samplecount, allelecount etc allready entered, with
		// room for the sample probabilities
		ByteBuffer genotypeDataBlockBuffer = initializeGenotypeDataBlockByteBuffer(
				sampleCount, alleleCount,
				minimumPloidy, maximumPloidy,
				ploidyMissingnessBytes, probabilitiesBytesCount, true); // data is phased, so true here.

		// Now start writing the sample probabilities to a byte array.
		int bitOffset = 0; // Start with an offset of zero.
		byte[] probabilitiesByteArray = new byte[totalNumberOfProbabilities * (probabilitiesLengthInBits / 8)];
		for (double[][] sampleProbabilities : sampleGenotypeProbabilitiesBgenPhased) {
			for (double[] alleleProbabilities : sampleProbabilities){
				// These probabilities should sum to one, add this to the byte buffer.
				// Also add the number of bits appended to the bit offset.
				bitOffset = addProbabilitiesSummingToOne(
						probabilitiesByteArray, bitOffset, alleleProbabilities);
			}
		}

		// Put the probabilities byte array into the genotype block byte buffer and return this.
		genotypeDataBlockBuffer.put(probabilitiesByteArray);
		return genotypeDataBlockBuffer;
	}

	/**
	 * Method that generates a single byte for a sample containing its missigness for a variant
	 * in the most significant bit, and the its ploidy in the 6 least significant bits.
	 *
	 * The ploidy should be a maximum of 63.
	 *
	 * @param ploidy The ploidy of the sample.
	 * @param missingness If genotype data is missing for the sample.
	 * @return The byte with missingness and ploidy.
	 */
	private byte getPloidyMissingnessByte(int ploidy, boolean missingness) {
		return (byte) ((byte) (ploidy & 63) | (missingness ? 128 : 0));
	}

	/**
	 * Write the genotype probability data from the byte buffer in compressed form to the output channel.
	 *
	 * @param bgenOutputByteChannel The channel that is able to pass byte buffers to the output stream.
	 * @param genotypeDataBlockBuffer The byte buffer that holds the genotype probability data.
	 * @throws IOException If an I/O exception has occurred.
	 * @return The number of bytes written.
	 */
	private long writeCompressedBgenGenotypeDataBlock(WritableByteChannel bgenOutputByteChannel,
													  ByteBuffer genotypeDataBlockBuffer) throws IOException {

		// Compress the genotype data using the Zstandard compression method.
		// Compression level 3 is used since this was the defualt when writing this.
		// (The byte buffer compression method does not have have an overloaded method with default compression level)
		ByteBuffer compressedGenotypeDataBuffer = Zstd.compress(genotypeDataBlockBuffer, 3);
		// If this was successful, the limit is equal to the position of the original buffer.
		// Check if the alternative condition applies.
		if (genotypeDataBlockBuffer.position() != genotypeDataBlockBuffer.limit()) {
			throw new GenotypeDataException("Zstd compression of genotype data failed, exiting");
		}

		// Get the length of the compressed variant probability data.
		int lengthOfVariantProbabilityData =
				compressedGenotypeDataBuffer.limit() + 4; // Add 4 since the D field takes 4 bytes
		// Get the decompressed length of variant probability data.
		int decompressedLengthOfVariantProbabilityData = genotypeDataBlockBuffer.capacity();

		// Write the (decompressed) length to a bytebuffer.
		ByteBuffer genotypeBlockLengths = ByteBuffer.allocate(8)
				.order(ByteOrder.LITTLE_ENDIAN);

		genotypeBlockLengths.putInt(lengthOfVariantProbabilityData);
		genotypeBlockLengths.putInt(decompressedLengthOfVariantProbabilityData);

		genotypeBlockLengths.flip(); // reset pointer

		// First write the lengths.
		int numberOfBytesWritten = bgenOutputByteChannel.write(genotypeBlockLengths);
		// Now write the actual data.
		return numberOfBytesWritten + bgenOutputByteChannel.write(compressedGenotypeDataBuffer);
	}

	/**
	 * Method that writes a field with the length of the next field and this next field with the given value.
	 * The size of the this value
	 * is at maximum the largest possible value that fits within the given maximum size of the length field.
	 *
	 * The max length of the length field should be 2 or 4 bytes.
	 *
	 * @param outputByteChannel The output channel to write the fields to
	 * @param fieldValue The value of the field to write.
	 * @param maxLengthFieldSizeInBytes The maximum number of bytes of the length field.
	 * @param fieldName The name of the field, for message formatting.
	 * @return The number of bytes written to the output byte channel
	 */
	private int writeFieldWithFieldLength(WritableByteChannel outputByteChannel, String fieldValue, int maxLengthFieldSizeInBytes, String fieldName) throws IOException {
		assert(maxLengthFieldSizeInBytes == 2 || maxLengthFieldSizeInBytes == 4);
		byte[] fieldValueAsBytes = fieldValue.getBytes(CHARSET);
		long lengthOfFieldInBytes = fieldValueAsBytes.length;

		// Get the maximum number of bytes that the field can have given the maximum the length field can hold.
		int maxLengthOfFieldInBytes = (int) Math.min(Math.pow(2, maxLengthFieldSizeInBytes * 8) - 1, Integer.MAX_VALUE);
		if (lengthOfFieldInBytes > maxLengthOfFieldInBytes) {
			// The length of the field value exceeds the max number of bytes. Give a warning and truncate.
			LOGGER.warn(String.format(
					"Length of %s %.16s... exceeds the maximum number of bytes. (%d vs %d bytes respectively) " +
							"%nField will be truncated.",
					fieldName, fieldValue, fieldValueAsBytes.length, maxLengthOfFieldInBytes));
			// ...truncate...
			lengthOfFieldInBytes = maxLengthOfFieldInBytes;
			fieldValueAsBytes = Arrays.copyOfRange(fieldValueAsBytes, 0,
					maxLengthOfFieldInBytes + 1);
		}

		// Initialize a byte buffer with the correct size.
		ByteBuffer buffer = ByteBuffer.allocate((int) lengthOfFieldInBytes)
				.order(ByteOrder.LITTLE_ENDIAN);
		buffer.put(fieldValueAsBytes);
		ByteBuffer fieldLengthBuffer = ByteBuffer.allocate(maxLengthFieldSizeInBytes)
				.order(ByteOrder.LITTLE_ENDIAN);
		switch (maxLengthFieldSizeInBytes) {
			case 2:
				fieldLengthBuffer.putShort((short) lengthOfFieldInBytes);
				break;
			case 4:
				fieldLengthBuffer.putInt((int) lengthOfFieldInBytes);
				break;
		}
		buffer.flip(); // reset pointer
		fieldLengthBuffer.flip(); // reset pointer

		// Write the buffer to the channel
		int writtenBytes = outputByteChannel.write(fieldLengthBuffer);
		return writtenBytes + outputByteChannel.write(buffer);
	}

	/**
	 * Method returning the fractional part of a given value.
	 *
	 * @param value The value to use.
	 * @return The fracional part.
	 */
	private static double fractionalPart(double value) {
		return value - Math.floor(value);
	}

	/**
	 * Method that should return the a number of times an object is chosen (ploidy),
	 * given a number of ordered combination generated
	 * using a given number of objects (alleles) that are repetitively chosen.
	 *
	 * @param numberOfProbabilities the number of combinations.
	 * @param alleleCount The number of alleles that are counted.
	 * @return the ploidy.
	 */
	public static int getPloidy(int numberOfProbabilities, int alleleCount) {
		List<Integer> pascalRange = IntStream.rangeClosed(1, alleleCount).boxed()
				.collect(Collectors.toList());
		int currentR = 1;
		while (numberOfProbabilities > pascalRange.get(alleleCount - 1)) {
			for (int i = 1; i < pascalRange.size(); i++) {
				pascalRange.set(i, pascalRange.get(i - 1) + pascalRange.get(i));
			}
			currentR++;
		}
		return currentR;
	}

	/**
	 * Method that is responsible for writing an array of probabilities that should sum to one.
	 * The last probability is not included since this can be inferred from the other values.
	 *
	 * The implementation is based on the following c++ reference implementation retrieved on 19-12-02.
	 * https://bitbucket.org/gavinband/bgen/src/default/src/bgen.cpp
	 *  @param bytes A byte array to write the probabilities to.
	 * @param bitOffset The offset within the byte array to write the probabilities from.
	 * @param probabilities The array of probabilities to write.
	 * @return the given bit offset + the number of written bits.
	 */
	private int addProbabilitiesSummingToOne(byte[] bytes, int bitOffset, double[] probabilities) {

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

		// Sort the indices according to the fractional part of their corresponding value.
		indices.sort(Comparator.comparingDouble(a -> fractionalPart(probabilities[a])));

		// Loop through the integers and use ceil and floor to make the sum of these integers 1.
		for (int i = 0; i < roundedSumOfFractionalPart; i++) {
			probabilities[indices.get(i)] = Math.ceil(probabilities[indices.get(i)]);
		}
		for (int i = (int) roundedSumOfFractionalPart;
			 i < probabilities.length; i++) {
			probabilities[indices.get(i)] = Math.floor(probabilities[indices.get(i)]);
		}

		// Write the probabilities in the byte array.
		for (int i = 0; i < (probabilities.length - 1); i++) {
			probabilityValueToByteArray((long) probabilities[i], bytes, bitOffset, probabilitiesLengthInBits);
			bitOffset += probabilitiesLengthInBits;
		}
		return bitOffset;
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
		// Get the total bytes min
		int totalBytesMin1 = (remainingBitOffset + numberOfBitsToWrite - 1) / 8;
		// Find how much bits to shift after the first byte (most significant)
		int bitShiftAfterFirstByte = 8 - remainingBitOffset;
		// Jump into the switch at the correct point and fall-through to fill the remaining bytes.
		switch (totalBytesMin1) {
			case 4: // Bitshift the value at least 3 bytes
				bytes[firstByteIndex + 4] |= (value >>> (24 + bitShiftAfterFirstByte));
			case 3: // Bitshift the value at least 2 bytes
				bytes[firstByteIndex + 3] |= (value >>> (16 + bitShiftAfterFirstByte));
			case 2: // Bitshift the value at least 1 byte
				bytes[firstByteIndex + 2] |= (value >>> (8 + bitShiftAfterFirstByte));
			case 1:
				bytes[firstByteIndex + 1] |= (value >>> bitShiftAfterFirstByte);
			case 0:
				bytes[firstByteIndex] |= (value << remainingBitOffset);
		}
		// >>> and << is used since these both only pad with 0s at the right / left side.
	}

	/**
	 * Setter for the precision of probabilities in number of bits.
	 * @param bgenBitRepresentation The number of bits.
	 */
	public void setProbabilityPrecisionInBits(int bgenBitRepresentation) {
		if (bgenBitRepresentation < 1 || bgenBitRepresentation > 32) {
			throw new GenotypeDataException(String.format("The given probability precision '%d' was " +
					"not in range 1 - 32 bits (both inclusive)", bgenBitRepresentation));
		}
		probabilitiesLengthInBits = bgenBitRepresentation;
	}
}
