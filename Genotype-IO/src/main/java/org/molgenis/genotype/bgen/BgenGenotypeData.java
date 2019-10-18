/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.bgen;

import com.facebook.presto.orc.zstd.ZstdDecompressor;
import org.apache.log4j.Logger;
import org.molgenis.genotype.*;
import org.molgenis.genotype.annotation.Annotation;
import org.molgenis.genotype.annotation.SampleAnnotation;
import org.molgenis.genotype.util.FixedSizeIterable;
import org.molgenis.genotype.variant.*;
import org.molgenis.genotype.variant.range.GeneticVariantRange;
import org.molgenis.genotype.variant.sampleProvider.CachedSampleVariantProvider;
import org.molgenis.genotype.variant.sampleProvider.SampleVariantUniqueIdProvider;
import org.molgenis.genotype.variant.sampleProvider.SampleVariantsProvider;

import java.io.EOFException;
import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.charset.Charset;
import java.util.*;
import java.util.zip.DataFormatException;
import java.util.zip.Inflater;

/**
 *
 * @author Patrick Deelen
 */
public class BgenGenotypeData extends AbstractRandomAccessGenotypeData implements SampleVariantsProvider {

	public enum blockRepresentation {
		compression_0, compression_1, compression_2
	}

	public enum Layout {
		layOut_1, layOut_2
	}

	/**
	 * Index is the number of last bits used from the first byte
	 */
	private static final int[] LAST_BYTE_MASK = {255, 1, 3, 7, 15, 31, 63, 127, 255};
	/**
	 * Index is the number of first bits used from the last byte
	 */
	private static final int[] FIRST_BYTE_MASK = {0, 128, 192, 224, 240, 248, 252, 254, 255};

	private static final double DEFAULT_MINIMUM_POSTERIOR_PROBABILITY_TO_CALL = 0.4f;
	private final RandomAccessFile bgenFile;
	private final byte[] byteArray4 = new byte[4]; //resuable 4 byte array
	private final byte[] byteArray2 = new byte[2]; //resuable 2 byte array
	private static final Logger LOGGER = Logger.getLogger(BgenGenotypeData.class);
	private final boolean sampleIdentifiersPresent;
	private final Inflater gzipInflater = new Inflater();
	private final ZstdDecompressor zstdInflater = new ZstdDecompressor();
	private static final Charset CHARSET = Charset.forName("UTF-8");
	private final blockRepresentation snpBlockRepresentation;
	private final Layout fileLayout;
	private final List<Sample> samples = new ArrayList<>();
	private final LinkedHashSet<String> sequenceNames = new LinkedHashSet<String>();;
	private static final GeneticVariantMeta GP_VARIANT = GeneticVariantMetaMap.getGeneticVariantMetaGp();
	private final SampleVariantsProvider sampleVariantProvider;
	private final int sampleVariantProviderUniqueId;
	private final BgenixReader bgenixReader; // Was previously a final field
	private final long sampleCount;

	public BgenGenotypeData(File bgenFile, File sampleFile) throws IOException {
		this(bgenFile, sampleFile, 1000);
	}

	public BgenGenotypeData(File bgenFile, File sampleFile, int cacheSize) throws IOException {
		this(bgenFile, sampleFile, cacheSize, DEFAULT_MINIMUM_POSTERIOR_PROBABILITY_TO_CALL);
	}

	public BgenGenotypeData(File bgenFile, File sampleFile, int cacheSize, double minimumPosteriorProbabilityToCall) throws IOException {

		this.bgenFile = new RandomAccessFile(bgenFile, "r");

		// Get offset of variants in the first four bytes.
        long snpOffset = readFourBytesAsUInt32(
                "SNP offset",
                "Error reading bgen file header. File is corrupt");

        // Get offset the size of the header in the following four bytes.
        long headerSize = readFourBytesAsUInt32(
                "Header size (Lh)",
                "Error reading bgen file header. File is corrupt");

        // Throw an exception if the size of the header is smaller than.
        if (headerSize > snpOffset) {
			throw new GenotypeDataException(
			        "Error reading bgen file header. Header information is bigger than expected offset.");
		}

        // Get the number of variants in the file.
        long variantCount = readFourBytesAsUInt32(
                "Number of SNPs",
                "Error reading bgen file header. File is corrupt");

//               Not needed when using the bgenix indexing file.  
//		if (snpCount > Integer.MAX_VALUE) {
//			throw new GenotypeDataException("Found more than (2^31)-1 SNPs in bgen file. This is not supported");
//		}

        // Get the number of samples in the file.
        sampleCount = readFourBytesAsUInt32(
                "Number of samples",
                "Error reading bgen file header. File is corrupt");

        // Throw an exception if the number of samples exceeds the maximum value of an integer.
		if (sampleCount > Integer.MAX_VALUE) {
			throw new GenotypeDataException("Found more than (2^31)-1 samples in bgen file. This is not supported");
		}

		// Magic number is in the next four bytes but skipped over.

		// Skip over reserved and free area,
		// this seek works because there are 4 (interesting) bytes
		// in the beginning of the file and in the end of the header block.
		this.bgenFile.seek(headerSize);

		//Read flags
		if (this.bgenFile.read(byteArray4, 0, 4) != 4) {
			throw new GenotypeDataException("Error reading bgen file header. File is corrupt");
		}

		// Changed storing of version and data.
        // Read the SNP block representation (compression state).
        snpBlockRepresentation = readSnpBlockRepresentation();

        // Read the file layout.
		fileLayout = readFileLayout(byteArray4);

        // Throw an exception if the SNP block representation is 2 while layout 1 is used.
        if (snpBlockRepresentation.equals(blockRepresentation.compression_2) & fileLayout.equals(Layout.layOut_1)) {
			throw new GenotypeDataException("Invalid compression method for layout one observed. Trying to use ZSTD compression on layout one file, which is not supported.");
		}

        // Read the sample identifiers presence; set sampleIdentifiersPresent to true if present, false if absent.
        sampleIdentifiersPresent = readSampleIdentifiersPresence(byteArray4[3]);

        if (sampleIdentifiersPresent) {
			// Process sample identifier block.
            processSampleIdentifierBlock(snpOffset, headerSize);
        }

//		System.out.println("offset:"+ (snpOffset + 4));
		File bgenixFile = new File(bgenFile.getAbsolutePath() + ".bgi");
		System.out.println(bgenFile.getAbsolutePath() + ".bgi");

		// Get the start of the variant data block
		long pointerFirstSnp = snpOffset + 4;

        // Check if BGENIX file is present.
        // If it is non-existent, make a new one.
		if (!bgenixFile.exists()) {
			LOGGER.info("Creating bgenix file at: " + bgenixFile.getAbsolutePath());
			BgenixWriter bgenixWriter = new BgenixWriter(bgenixFile);
			createBgenixFile(
			        bgenFile,
                    bgenixWriter,
                    pointerFirstSnp);

			bgenixWriter.finalizeIndex();
		}

		// Not sure what bgenix reader should do in this object other than read an existing BGENIX file.
		bgenixReader = new BgenixReader(bgenixFile);
        readExistingBgenixFile(bgenFile);

		sampleVariantProviderUniqueId = SampleVariantUniqueIdProvider.getNextUniqueId();
		// If the specified cache size is greater than 0, construct a new SampleVariantProvider
		if (cacheSize > 0) {
			sampleVariantProvider = new CachedSampleVariantProvider(this, cacheSize);
		} else {
			sampleVariantProvider = this;
		}
		//Read the first snp to get into genotype-io.
//		readCompleteGeneticVariant(bgenFile, pointerFirstSnp, (int) sampleCount, this.fileLayout, this.snpBlockRepresentation);
	}

	/**
	 * Method responsible for reading a BGENIX file. This BGENIX file stores indexes of variants
	 * for quick random access to genotype data.
	 *
	 * @param bgenFile The BGENIX file to read.
	 * @throws IOException if an I/O error occurs
	 */
    private void readExistingBgenixFile(File bgenFile) throws IOException {
        BgenixMetadata metadata = bgenixReader.getMetadata();
        if (metadata != null) {
            // Check if the metadata of the read BGENIX file is equal to
            // the metadata of the BGEN file.
            if (!metadata.getFileName().equals(bgenFile.getName())) {
                throw new GenotypeDataException("Sample name between bgenix and bgen is not equal. Invalid Bgen and Bgenix combination.");
            }

            // Check if the BGEN file length corresponds to the expected file length in the BGENIX file.
            if (metadata.getFileSize() != bgenFile.length()) {
                throw new GenotypeDataException("Number of expected bytes is different to the numer of observed bytes. Invalid Bgen and Bgenix combination.");
            }

            // Read the first 1000 bytes to check if this is equal to that in the metadata in the
			// BGENIX file.
            this.bgenFile.seek(0);
            byte[] firstBytes = new byte[1000];
            this.bgenFile.read(firstBytes, 0, 1000);
            if (!Arrays.equals(metadata.getFirst1000bytes(), firstBytes)) {
                throw new GenotypeDataException("First 1000 bytes of meta data and actual data are not equal. Invalid Bgen and Bgenix combination.");
            }
        } else {
            // Also throw an exception whenever the BGENIX files metadata is null as returned by the getMetadata method
            // in the bgenixReader object.
            // Not sure if we're ok with null?
            throw new GenotypeDataException("Metadata could not be obtained from the BGENIX file.");
}
    }

	/**
	 * Method that processes the sample identifier block of the BGEN file.
	 *
	 * @param snpOffset The number of bytes the variant data is offset from the start of the file header.
	 * @param headerSize The size of the header within the BGEN file.
	 * @throws IOException if an I/O error has occurred.
	 */
    private void processSampleIdentifierBlock(long snpOffset, long headerSize) throws IOException {
	    // Get the length of the sample identifier data block.
        long byteSizeSampleIds = readFourBytesAsUInt32(
        		"Sample identifier block size",
				"Error in sample identifier block. File is corrupt.");

        // Throw an exception if the sample identifier block size plus the header size is more than the snp offset
        if ((byteSizeSampleIds + headerSize) > snpOffset) {
            throw new GenotypeDataException(
            		"Error reading bgen file header. Combination of header & sample id information is bigger than expected offset.");
        }

        // Read the number of samples represented in this file
        long sampleCountFromSampleIdBlock = readFourBytesAsUInt32(
        		"Number of samples with sample-ids",
				"Error in sample identifier block. File is corrupt.");

        // Throw an exception if the number of samples given by the metadata in the sample identifier block and the
        // header are not equal
        if (sampleCountFromSampleIdBlock != sampleCount) {
            throw new GenotypeDataException("Number of samples in metadata and sample id data is not equal. File is corrupt.");
        }

        // Initialize an array of Strings representing the sample identifiers
        String[] sampleIds = new String[(int) sampleCountFromSampleIdBlock];

        // Read the sample identifiers sample by sample.
        for (int i = 0; i < sampleCountFromSampleIdBlock; i++) {
            // Read the sample id length within the next two bytes
            if (this.bgenFile.read(byteArray2, 0, 2) != 2) {
                throw new GenotypeDataException("Error in sample Id. File is corrupt.");
            }
            int sampleIdLength = getUInt16(byteArray2, 0);
            // Initialize the sample identifier.
            byte[] sampleName = new byte[sampleIdLength];
            // Read the sample identifier.
            this.bgenFile.read(sampleName, 0, sampleIdLength);
            // Append the sample identifier to the array of sample ids.
            sampleIds[i] = new String(sampleName, CHARSET);
//				System.out.println("\t"+sampleIds[i]);
			samples.add(new Sample(sampleIds[i], null, null));
		}
    }

	/**
	 * Read the field that represents the presence of sample identifiers within the BGEN file.
	 *
	 * @param sampleIdentifiersField The field representing the presence of sample identifiers
	 *                               within the BGEN file.
	 * @return true if the smaple identifiers are present, false if not.
	 */
    private boolean readSampleIdentifiersPresence(byte sampleIdentifiersField) {
        if ((sampleIdentifiersField & 128) == 128) {
            LOGGER.debug("SampleIdentifiers present in bgen file");
            return true;
        } else {
            LOGGER.debug("SampleIdentifiers not present in bgen file");
            return false;
        }
    }

	/**
	 * Read the field that indicates which layout is used for this specific BGEN file.
	 *
	 * @return the layout type.
	 * @param flagArray Byte array with flags of length 4 according to BGEN specifications.
	 */
    private Layout readFileLayout(byte[] flagArray) {
    	// Layout flag is located within the first byte.
		byte byte_tmp = flagArray[0];

		// Shift the layout flag, located on bit 2-5 relative to the first bit.
		byte_tmp = (byte) (byte_tmp >> 2); // relocate to bit 0-3.

		// Perform bitwise operation to mask everything outside bit 0-3.
		int layoutFlag = byte_tmp & 7; // mask using 7 (00000111) and convert to an integer.

		// Find the layout that corresponds to the layout flag value.
        switch (layoutFlag) {
            case 1:
                LOGGER.debug("Genotype data is in layout 1 (BGEN 1.1)");
                return Layout.layOut_1;
            case 2:
                LOGGER.debug("Genotype data is in layout 2 (BGEN 1.2 & 1.3)");
                return Layout.layOut_2;
            default:
                throw new GenotypeDataException("Invalid layout, observed: " + (byteArray4[0] & 3));
        }
    }

	/**
	 * Read the flag that indicates how the SNP probability block is represented.
	 * Data can be zlib compressed, zstd compressed (for layout 2), or uncompressed.
	 *
	 * @return the representation of the SNP probabilities block.
	 */
    private blockRepresentation readSnpBlockRepresentation() {

    	// Perform bitwise operation to mask everything outside the 0, 1 bits.
		int blockRepresentationFlag = byteArray4[0] & 3; // mask with 7 (00000011) and convert to an integer.

		// Find the representation that corresponds to the flag.
		switch (blockRepresentationFlag) {
            case 0:
                LOGGER.debug("Genotype data is not compressed.");
                return blockRepresentation.compression_0;
            case 1:
                LOGGER.debug("Genotype data is zlib compressed");
                return blockRepresentation.compression_1;
            case 2:
                LOGGER.debug("Genotype data is zstd compressed");
                return blockRepresentation.compression_2;
            default:
                throw new GenotypeDataException("Invalid compression method, observed: " + (blockRepresentationFlag));
        }
    }

	/**
	 * @param fieldName The name of the field to read.
	 * @param exceptionMessage What to report if the length of bytes read is not equal to four.
	 * @return the bytes as a long.
	 * @throws IOException if an I/O error has occurred.
	 */
    private long readFourBytesAsUInt32(String fieldName, String exceptionMessage) throws IOException {
        // Throw an exception if the read number of bytes is not equal to 4.
    	if (this.bgenFile.read(byteArray4, 0, 4) != 4) {
            throw new GenotypeDataException(exceptionMessage);
        }
    	// Convert the read bytes to a long
        long value = getUInt32(byteArray4, 0);

        LOGGER.debug(fieldName + ": " + value);
        return value;
    }

	/**
	 * Method responsible for creating a BGENIX file using the BGEN file as guide.
	 *
	 * @param bgen The bgen file that is being read.
	 * @param bgenixWriter The bgenixWriter to use for creating this BGENIX file.
	 * @param pointerFirstSnp The byte index of the file to start reading variants.
	 * @throws IOException if an I/O error has occurred.
	 */
    private void createBgenixFile(File bgen,
								  BgenixWriter bgenixWriter,
								  long pointerFirstSnp) throws IOException {
		// Go to the first byte...
	    this.bgenFile.seek(0);

	    // Read the first 1000 bytes of the bgen file.
		byte[] firstBytes = new byte[1000];
		this.bgenFile.read(firstBytes, 0, 1000);

		//Add current time in int.
		System.out.println((System.currentTimeMillis() / 1000L));
		// Create and write new metadata.
		BgenixMetadata m = new BgenixMetadata(
				bgen.getName(),
				(int) this.bgenFile.length(),
				(int) (bgen.lastModified() / 1000L),
				firstBytes,
				(System.currentTimeMillis() / 1000L));
		bgenixWriter.writeMetadata(m);

		//Loop through the start of the file
		long variantReadingPosition = pointerFirstSnp;
		while ((variantReadingPosition) < bgenFile.length()) {
			//Loop through variants.

			// Read the variant data creating a variant.
			GeneticVariant variant = processVariantIdentifyingData(variantReadingPosition);

			// Get the current position of the file pointer.
			long currentPointer = this.bgenFile.getFilePointer();
			// Extract the variant genotype data block info starting from the current position of the
			// file pointer, which should be right after all alleles for a specific variant.
			VariantGenotypeDataBlockInfo variantGenotypeDataBlockInfo =
					extractVariantGenotypeDataBlockInfo(currentPointer);

			// Add the read variant to the BGENIX file so that it can quickly be retrieved.
			bgenixWriter.addVariantToIndex(
					variant,
					pointerFirstSnp,
					(int) variantGenotypeDataBlockInfo.getBlockLength(),
					variant.getVariantId().getPrimairyId());

			readGenotypesFromVariant(currentPointer);
			break;

			// Get the position in the file to start reading the next variant.
//			variantReadingPosition = currentPointer + variantGenotypeDataBlockInfo.getBlockLengthHeaderInclusive();
		}
	}

	/**
	 * Processes variant identifying data from the BGENIX file, using a the variant start
	 * position to read from the right location.
	 *
	 * @param variantStartPosition The position to start reading the variant from.
	 * @return a genetic variant.
	 * @throws IOException if an I/O error has occurred.
	 */
	private GeneticVariant processVariantIdentifyingData(long variantStartPosition) throws IOException {
		LOGGER.debug("variantStartPosition: " + variantStartPosition);

		// If layout is equal to 1 then the variant identifying data starts with 4 bytes describing the
		// number of individuals within the row.
		if (fileLayout == Layout.layOut_1) {
			// We chose to ignore this apparently ...
			variantStartPosition += 4;
		}

		// Go to the start of the variant to begin reading there.
		this.bgenFile.seek(variantStartPosition);

		// //Not sure if we want to do the buffer search here. Or we might be able to take a smaller set.
		byte[] snpInfoBuffer = new byte[8096];
		this.bgenFile.read(snpInfoBuffer, 0, snpInfoBuffer.length);
		int snpInfoBufferPos = 0;

//		if (snpInfoBufferSize < 20) {
//			throw new GenotypeDataException("Error reading bgen snp data. File is corrupt");
//		}
		// Need to check that it is correct with the block in front of the snp id.

		// Get the ids of the variant
		ArrayList<String> VariantIds = new ArrayList<>();
		LOGGER.debug("SNPinfoBufferPos: " + snpInfoBufferPos);
		int fieldLength = getUInt16(snpInfoBuffer, snpInfoBufferPos);
		LOGGER.debug("Snp ID length " + fieldLength);
		snpInfoBufferPos += 2;
		String snpId = new String(snpInfoBuffer, snpInfoBufferPos, fieldLength, CHARSET);
		LOGGER.debug("SNP ID: " + snpId);
//		System.out.println(snpId);
		snpInfoBufferPos += fieldLength; // skip id length and snp id

		fieldLength = getUInt16(snpInfoBuffer, snpInfoBufferPos);
		LOGGER.debug("Snp RS length " + fieldLength);

		snpInfoBufferPos += 2;
		String snpRsId = new String(snpInfoBuffer, snpInfoBufferPos, fieldLength, CHARSET);
		LOGGER.debug("SNP RSID: " + snpRsId);

//		System.out.println(snpRsId);
		snpInfoBufferPos += fieldLength;
		VariantIds.add(snpRsId);
		VariantIds.add(snpId);

		fieldLength = getUInt16(snpInfoBuffer, snpInfoBufferPos);
		snpInfoBufferPos += 2;
		String seqName = new String(snpInfoBuffer, snpInfoBufferPos, fieldLength, CHARSET);
//		System.out.println(seqName);
		snpInfoBufferPos += fieldLength;

		// Get the position of the variant.
		int variantPosition = getVariantPosition(snpInfoBuffer, snpInfoBufferPos);
		snpInfoBufferPos += 4;
//		System.out.println("SNP pos " + variantPosition);

		// Get the alleles for this variant.
		int numberOfAlleles = 2;
		if (fileLayout.equals(Layout.layOut_2)) {
			numberOfAlleles = getUInt16(snpInfoBuffer, snpInfoBufferPos);
			snpInfoBufferPos += 2;
//			System.out.println("SNP Alleles " + numberOfAlleles);
		}

		ArrayList<String> alleles = new ArrayList<>();
		for (int i = 0; i < numberOfAlleles; i++) {
			snpInfoBufferPos = readAllele(snpInfoBuffer, snpInfoBufferPos, alleles);
		}

		// Seek to the end of the read data.
		this.bgenFile.seek(snpInfoBufferPos + variantStartPosition);

		return ReadOnlyGeneticVariant.createVariant(
				GP_VARIANT, VariantIds, variantPosition, seqName, sampleVariantProvider, alleles, alleles.get(0));
	}

	/**
	 * Method for extracting the position of the current variant.
	 *
	 * @param snpInfoBuffer The byte array buffer starting from the start of the variant block.
	 * @param snpInfoBufferPos The position of the variant position field in the snp info buffer.
	 * @return the position of the variant.
	 */
	private int getVariantPosition(byte[] snpInfoBuffer, int snpInfoBufferPos) {
		long variantPosLong = getUInt32(snpInfoBuffer, snpInfoBufferPos);
		// Throw an exception if the variant position exceeds the maximum value of an integer.
		if (variantPosLong > Integer.MAX_VALUE) {
			throw new GenotypeDataException("SNP pos larger than (2^31)-1 not supported");
		}
		return (int) variantPosLong;
	}

	/**
	 * Read an allele and add this to the list of alleles
	 *
	 * @param snpInfoBuffer A byte array buffer starting from the start of the variant block.
	 * @param snpInfoBufferPos The position of an allele block in the snp info buffer.
	 * @param alleles A list of alleles.
	 * @return the end position of the allele block.
	 */
	private int readAllele(byte[] snpInfoBuffer, int snpInfoBufferPos, List<String> alleles) {

		// Length of the allele
		long fieldLengthLong = getUInt32(snpInfoBuffer, snpInfoBufferPos);
		snpInfoBufferPos += 4;

		// Throw an exception if the allele length is longer
		if (fieldLengthLong > Integer.MAX_VALUE) {
			throw new GenotypeDataException("SNP with allele longer than (2^31)-1 characters not supported");
		}

		// Get the allele from the buffer.
		String allele = new String(snpInfoBuffer, snpInfoBufferPos, (int) fieldLengthLong, CHARSET);
		snpInfoBufferPos += ((int) fieldLengthLong);
		alleles.add(allele);
		return snpInfoBufferPos;
	}

	private float[][] readGenotypesFromVariant(long filePointer) throws IOException {
		this.bgenFile.seek(filePointer);

		//Not sure if we want to do the buffer search here. Or we might be able to take a smaller set.
		byte[] snpInfoBuffer = new byte[8096];
//		int snpInfoBufferSize = 
		this.bgenFile.read(snpInfoBuffer, 0, snpInfoBuffer.length);
		int snpInfoBufferPos = 0;

		float[][] probabilities = new float[(int) sampleCount][3];

		// Read
		if (fileLayout == Layout.layOut_1) {
			VariantGenotypeDataBlockInfo blockInfo = getVariantGenotypeDataBlockInfoForLayoutOne(
					snpInfoBuffer,
					snpInfoBufferPos);

            byte[] variantBlockData = getDecompressedBlockData(filePointer, blockInfo);
			for (int sampleIndex = 0; sampleIndex < sampleCount; sampleIndex++) {
				float[] sampleProbabilities = new float[3];

//				System.out.println(getUInt16(variantBlockData, 0) / 32768f + " " + getUInt16(variantBlockData, 2) / 32768f + " " + getUInt16(variantBlockData, 4) / 32768f);
				sampleProbabilities[0] = getUInt16(variantBlockData, 0) / 32768f;
				sampleProbabilities[1] = getUInt16(variantBlockData, 2) / 32768f;
				sampleProbabilities[2] = getUInt16(variantBlockData, 4) / 32768f;
				probabilities[sampleIndex] = sampleProbabilities;
			}
			return probabilities;

        } else if (fileLayout.equals(Layout.layOut_2)) {
			VariantGenotypeDataBlockInfo blockInfo = determineVariantGenotypeDataBlockSizeForLayoutTwo(
					snpInfoBuffer,
					snpInfoBufferPos);

            System.out.println("Snp block size: " + blockInfo.getBlockLength());
			System.out.println("Snp block size decompressed: " + blockInfo.getDecompressedBlockLength());

			// Get the decompressed variant data of length D
            byte[] variantBlockData = getDecompressedBlockData(filePointer, blockInfo);
			System.out.println("variantBlockData.length = " + variantBlockData.length);

            int blockBufferOffset = 0;
			//must equal data before.
			if ((int) getUInt32(variantBlockData, blockBufferOffset) != sampleCount) {
                throw new GenotypeDataException(
                        "BGEN file format error. Sample counts do not match.");
            }
			System.out.println("Number of individuals: " + sampleCount);
			blockBufferOffset += 4;

			//must equal data before.
			int numberOfAlleles = getUInt16(variantBlockData, blockBufferOffset);
			System.out.println("Number of Alleles: " + numberOfAlleles);
			blockBufferOffset += 2;

			// Read the minimum ploidy
			int minPloidy = getUInt8(variantBlockData, blockBufferOffset);
			System.out.println("Min ploidy: " + minPloidy);
			blockBufferOffset += 1;

			// Read the maximum ploidy
			int maxPloidy = getUInt8(variantBlockData, blockBufferOffset);
			System.out.println("Max ploidy: " + maxPloidy);
			blockBufferOffset += 1;

			// Loop through every individual
			List<Boolean> isMissing = new ArrayList<>();
			List<Integer> ploidies = new ArrayList<>();
			ReadPloidiesAndMissingnessByte(variantBlockData, blockBufferOffset, isMissing, ploidies);
			// Add the number of individuals to the buffer.
			blockBufferOffset += sampleCount;

			boolean phased = isPhased(variantBlockData, blockBufferOffset);
			blockBufferOffset += 1;

			int probabilitiesLengthInBits = getUInt8(variantBlockData, blockBufferOffset);
			System.out.println("Bit representation of probability: " + probabilitiesLengthInBits);
			blockBufferOffset += 1;

			byte[] probabilitiesArray = Arrays.copyOfRange(
					variantBlockData, blockBufferOffset,
					variantBlockData.length);

			if (phased) {
				readHaplotypeProbabilities(
						probabilitiesArray,
						probabilitiesLengthInBits,
						numberOfAlleles,
						isMissing, ploidies);
			} else {
				double[][] completeGenotypeProbabilities = readGenotypeProbabilities(
						probabilitiesArray,
						probabilitiesLengthInBits,
						numberOfAlleles,
						isMissing, ploidies);
			}
		}
		return probabilities;
	}

	private double calculateProbability(byte[] probabilitiesArray, int bitOffset, int probabilitiesLengthInBits) {
		// Each probability is stored in B bits.
		// Values are interpreted by linear interpolation between 0 and 1;
		// value b corresponds to probability b / ((2^B)-1).

		// To interpret a stored value x as a probability:
		// Convert x to an integer in floating point representation
		double probability = readProb(probabilitiesArray, bitOffset, probabilitiesLengthInBits);

		// Divide by (2^B)-1
		double maxValue = Math.pow(2, probabilitiesLengthInBits) - 1;
		probability = probability / maxValue;
		return probability;
	}

	private double[] computeApproximateProbabilities(byte[] probabilitiesArray, int probabilitiesLengthInBits, int bitOffset, int numberOfProbabilities) {
		// Keep track of the probabilities
		double[] probabilities = new double[numberOfProbabilities];
		// Keep track of the sum of all probabilities as the last value is calculated by
		// subtracting this sum from 1.
		double sumOfProbabilities = 0;
		for (int j = 0; j < numberOfProbabilities - 1; j++) {
			// Get the probability
			double probability = calculateProbability(probabilitiesArray, bitOffset, probabilitiesLengthInBits);
			probabilities[j] = probability;
			sumOfProbabilities += probability;
			// Update the offset of the current bit
			bitOffset += probabilitiesLengthInBits;
		}
		// Calculate probability of Kth allele (P_iK)
		probabilities[numberOfProbabilities - 1] = 1 - sumOfProbabilities;
		return probabilities;
	}

	private double[][] readGenotypeProbabilities(
			byte[] probabilitiesArray,
			int probabilitiesLengthInBits,
			int numberOfAlleles, List<Boolean> isMissing,
			List<Integer> ploidies) {

		// Get bit offset
		int bitOffset = 0;

		// Initialize an array of probabilities.
		double[][] probabilities = new double[getSamples().size()][];

		for (int sampleIndex = 0; sampleIndex < sampleCount; sampleIndex++) {
			if (ploidies.get(sampleIndex) != 2) {
				throw new UnsupportedOperationException(String.format(
						"The genotype of a sample with a ploidy of %d was requested, " +
								"but only a ploidy of 2 is supported", ploidies.get(sampleIndex)));
			}

			// Get the list of combinations that should have probabilities in this probability-block.
			List<List<Integer>> combinations = combinations(
					numberOfAlleles,
					ploidies.get(sampleIndex));

			// If the probabilities are missing for this sample, read zero and continue with the
			// next sample.
			if (isMissing.get(sampleIndex)) {
				// If this is missing, the probability is zero.
				bitOffset += probabilitiesLengthInBits * (combinations.size() - 1);
				continue;
			}

			// Get all probabilities for this sample.
			double[] genotypeProbabilities = computeApproximateProbabilities(
					probabilitiesArray, probabilitiesLengthInBits, bitOffset, combinations.size());
			probabilities[sampleIndex] = genotypeProbabilities;

			// Update the bit to read the next probabilities from.
			bitOffset += probabilitiesLengthInBits * (combinations.size() - 1);

			System.out.println("genotypeProbabilities = " + Arrays.toString(genotypeProbabilities));
		}
		return probabilities;
	}

	private void readHaplotypeProbabilities(
			byte[] probabilitiesArray,
			int probabilitiesLengthInBits,
			int numberOfAlleles, List<Boolean> isMissing,
			List<Integer> ploidies) {

    	// Define an array consisting of an array of posterior probabilities for each genotype
		float[][] probabilities = new float[getSamples().size()][3];

		// Get bit offset
		int bitOffset = 0;

		// Each probability is stored in B bits.
		// Values are interpreted by linear interpolation between 0 and 1;
		// value b corresponds to probability b / ((2^B)-1).

		for (int sampleIndex = 0; sampleIndex < sampleCount; sampleIndex++) {
			// If the probabilities are missing for this sample, read zero and continue with the
			// next sample.
			if (isMissing.get(sampleIndex)) {
				// If this is missing, the probability is zero.
				bitOffset += probabilitiesLengthInBits;
				continue;
			}

			float[] sampleProbabilities = new float[3];
			// Calculate the probability for homozygous genotype 'AA',
			// the probability for 'AB', and the probability for 'BB'

			Integer ploidy = ploidies.get(sampleIndex);
			for (int i = 0; i < ploidy; i++) {
				// Get the probabilities for every allele in this haplotype.
				double[] alleleProbabilities = computeApproximateProbabilities(
						probabilitiesArray, probabilitiesLengthInBits, bitOffset, numberOfAlleles);

				// What now?

				// We have per sample probabilities for A1, B1, A2, B2
				// Which can be converted to AA, AB, BA, BB
			}

		}
	}

	/**
	 * Read the ploidity and missingness.
	 *  @param snpBlockData The byte buffer array to read the genotype data.
	 * @param byteOffset the byte number to start reading ploidies from within the input byte array.
	 * @param isMissing A list of booleans representing if the corresponding sample has missing probabilities
	 * @param ploidity A list representing the ploidy of corresponding samples.
	 */
	private void ReadPloidiesAndMissingnessByte(
			byte[] snpBlockData, int byteOffset,
			List<Boolean> isMissing, List<Integer> ploidity) {

		for (int i = byteOffset; i < sampleCount + byteOffset; i++) {
			// Here we need to handle ploidity and missigness of prababilities.
			// Missingness is encoded by the most significant bit (MSB);
			// thus a value of 1 for the most significant bit
			// indicates that no probability data is stored for this sample.
			// Ploidity is encoded by the 6 least significant bits.

			// Get the MSB and check if the result equals 128.
			isMissing.add((snpBlockData[i] & 0x80) == 128);
			// Get the 6 least significant bits.
			ploidity.add(snpBlockData[i] & 0x3F);

//			String s1 = String.format("%8s", Integer.toBinaryString(snpBlockData[i] & 0xFF)).replace(' ', '0');
//			System.out.println("value = " + s1 + " | " + (snpBlockData[i] & 0x3F) + " | " + (snpBlockData[i] & 0x80));
		}
	}

	/**
	 * Read if the variant data is phased or not.
	 *
	 * @param snpBlockData The byte buffer array to read the genotype data.
	 * @param isPhasedFieldPosition The position in the byte buffer array of the field representing the phasing
	 * @return true if the data is phased, false if not.
	 */
	private boolean isPhased(byte[] snpBlockData, int isPhasedFieldPosition) {
		boolean phased;
		// Get the field with the phasing flag.
		switch (getUInt8(snpBlockData, isPhasedFieldPosition)) {
			case 0:
				phased = false;
				break;
			case 1:
				phased = true;
				break;
			default:
				throw new GenotypeDataException(
						"Bgen file format error. Unsupported value for the phased flag observed."
				);
		}
		System.out.println("phased: " + phased);
		return phased;
	}

	/**
	 * Method for obtaining decompressed data for layout 2.
	 *
	 * @param variantDataBlockOffset The position of the variant data block in bytes.
	 * @param blockInfo The info of the variant genotype data block.
	 */
    private byte[] getDecompressedBlockData(long variantDataBlockOffset, VariantGenotypeDataBlockInfo blockInfo) throws IOException {
		int variantBlockLength = Math.toIntExact(blockInfo.getBlockLength());
		int decompressedVariantBlockLength = Math.toIntExact(blockInfo.getDecompressedBlockLength());

		// Initialize byte arrays.
		byte[] compressedBlockData = new byte[variantBlockLength];
		byte[] decompressedBlockData = new byte[decompressedVariantBlockLength];

		// Read the compressed / uncompressed data starting from the correct location.
		this.bgenFile.seek(variantDataBlockOffset + blockInfo.getBlockOffset());
		this.bgenFile.read(compressedBlockData, 0, variantBlockLength);

        switch (snpBlockRepresentation) {

            case compression_1:
                decompressVariantBlockGzip(
                		compressedBlockData,
						decompressedBlockData);
                break;

            case compression_2:
                zstdInflater.decompress(
                		compressedBlockData,
						0,
						variantBlockLength,
						decompressedBlockData,
						0,
						decompressedVariantBlockLength);
                break;

            default:
				decompressedBlockData = compressedBlockData;
                break;
        }
        return decompressedBlockData;
    }

	/**
	 * Method that decompresses data using the gzipInflater.
	 *  @param compressedVariantDataBlock The input byte array to decompress.
	 * @param outputVariantDataBlock The decompressed output byte array.
	 */
    private void decompressVariantBlockGzip(
			byte[] compressedVariantDataBlock, byte[] outputVariantDataBlock) {

    	// Set the input for the gzip inflater.
		gzipInflater.setInput(compressedVariantDataBlock);

		// Try to decompress the data.
		try {
			gzipInflater.inflate(outputVariantDataBlock);
		} catch (DataFormatException e) {
			throw new GenotypeDataException("Error decompressing bgen data", e);
		}
		gzipInflater.reset();
	}

	/**
	 * Method that extracts info about the variant genotype data block.
	 *
	 * @param filePointer The position in the BGEN file where the variant genotype data block starts.
	 * @return an object with info of this block.
	 * @throws IOException if an I/O error has occurred.
	 */
	private VariantGenotypeDataBlockInfo extractVariantGenotypeDataBlockInfo(
			long filePointer) throws IOException {

		// Seek towards the file pointer (bgenFile.getFilePointer();)
		this.bgenFile.seek(filePointer);

		// Not sure if we want to do the buffer search here. Or we might be able to take a smaller set.
		byte[] snpInfoBuffer = new byte[8096];
		this.bgenFile.read(snpInfoBuffer, 0, snpInfoBuffer.length);
		int snpInfoBufferPos = 0;

		VariantGenotypeDataBlockInfo variantGenotypeDataBlockInfo = null;

		// Separate methods for layout 1 and 2.
		// Separation like this indicates that subclasses of this reader class should probably be made...
		if (fileLayout == Layout.layOut_1) {
			variantGenotypeDataBlockInfo = getVariantGenotypeDataBlockInfoForLayoutOne(
					snpInfoBuffer, snpInfoBufferPos);
		} else if (fileLayout == Layout.layOut_2) {
			variantGenotypeDataBlockInfo = determineVariantGenotypeDataBlockSizeForLayoutTwo(
					snpInfoBuffer, snpInfoBufferPos);
		}
		return variantGenotypeDataBlockInfo;
	}

	/**
	 * Method for obtaining info of a genotype data block for layout two.
	 *
	 * @param snpInfoBuffer A byte array buffer starting from the start of the variant block.
	 * @param snpInfoBufferPos The position of an genotype data block in the snp info buffer.
	 * @return an object with info of this block.
	 */
	private VariantGenotypeDataBlockInfo determineVariantGenotypeDataBlockSizeForLayoutTwo(
			byte[] snpInfoBuffer, int snpInfoBufferPos) {
		long variantBlockSize;
		long snpBlockSizeDecompressed;
		// If the file layout is 2,
		// The genotype data for this variant is equal to the next 4 bytes
		variantBlockSize = getUInt32(snpInfoBuffer, snpInfoBufferPos);
		// Shift the buffer position 4 bytes
		snpInfoBufferPos += 4;

		// If the probabilities block is compressed read the size of the decompressed data within
		// the next four bytes.
		if (!snpBlockRepresentation.equals(blockRepresentation.compression_0)) {
			snpBlockSizeDecompressed = getUInt32(snpInfoBuffer, snpInfoBufferPos);
			return new VariantGenotypeDataBlockInfo(
					fileLayout,
					variantBlockSize,
					snpBlockSizeDecompressed, true); // The data is compressed.
		}
		return new VariantGenotypeDataBlockInfo(
				fileLayout,
				variantBlockSize,
				false); // The data is not compressed.
	}

	/**
	 * Method for obtaining info of a genotype data block for layout one.
	 *
	 * @param snpInfoBuffer A byte array buffer starting from the start of the variant block.
	 * @param snpInfoBufferPos The position of an genotype data block in the snp info buffer.
	 * @return an object with info of this block.
	 */
	private VariantGenotypeDataBlockInfo getVariantGenotypeDataBlockInfoForLayoutOne(
			byte[] snpInfoBuffer, int snpInfoBufferPos) {
		long variantBlockSize;
		if (snpBlockRepresentation.equals(blockRepresentation.compression_1)) {
			// If the file layout is 1 and the variants are zlib compressed
			// The snp block starts in the 4 bytes after the block size field
			variantBlockSize = getUInt32(snpInfoBuffer, snpInfoBufferPos);
			return new VariantGenotypeDataBlockInfo(fileLayout, variantBlockSize, true);
		} else {
			// If the file layout is 1 and the variants are not zlib compressed,
			// the genotype data for this variant is 6 times the number of samples
			variantBlockSize = 6 * sampleCount;
			return new VariantGenotypeDataBlockInfo(fileLayout, variantBlockSize, false);
		}
	}

	@Override
	public List<Sample> getSamples() {
		return Collections.unmodifiableList(samples);
	}

	@Override
	public Map<String, Annotation> getVariantAnnotationsMap() {
		throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
	}

	@Override
	public Map<String, SampleAnnotation> getSampleAnnotationsMap() {
		throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
	}

	@Override
	public boolean isOnlyContaingSaveProbabilityGenotypes() {
		return true;
	}

	@Override
	public void close() throws IOException {
		bgenFile.close();
	}

	@Override
	public List<String> getSeqNames() {
		return new ArrayList<>(sequenceNames);
	}

	@Override
	public Iterable<Sequence> getSequences() {
		List<Sequence> sequences = new ArrayList<>();
		for (String seqName : getSeqNames()) {
			sequences.add(new SimpleSequence(seqName, null, this));
		}
		return sequences;
	}

	@Override
	public Iterator<GeneticVariant> iterator() {
		return getGeneticVariants(bgenixReader.getVariants()).iterator();
	}

	@Override
	public Iterable<GeneticVariant> getVariantsByPos(String seqName, int startPos) {
		return getGeneticVariants(bgenixReader.getVariantsPostion(seqName, startPos));
	}

	@Override
	public Iterable<GeneticVariant> getSequenceGeneticVariants(String seqName) {
		return getGeneticVariants(bgenixReader.getVariantsChromosome(seqName));
	}

	@Override
	public Iterable<GeneticVariant> getVariantsByRange(String seqName, int rangeStart, int rangeEnd) {
		return getGeneticVariants(bgenixReader.getVariantsRange(seqName, rangeStart, rangeEnd));
	}

	@Override
	public List<Alleles> getSampleVariants(GeneticVariant variant) {
		throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
	}

	@Override
	public FixedSizeIterable<GenotypeRecord> getSampleGenotypeRecords(GeneticVariant variant) {
		throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
	}

	@Override
	public List<Boolean> getSamplePhasing(GeneticVariant variant) {
		throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
	}

	@Override
	public int cacheSize() {
		return 0;
	}

	@Override
	public int getSampleVariantProviderUniqueId() {
		return sampleVariantProviderUniqueId;
	}

	@Override
	public byte[] getSampleCalledDosage(GeneticVariant variant) {
		throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
	}

	@Override
	public float[] getSampleDosage(GeneticVariant variant) {
		throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
	}

	@Override
	public float[][] getSampleProbilities(GeneticVariant variant) {
//		throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
		long variantDataPosition = getVariantDataPosition(variant);
//		if (fileLayout == Layout.layOut_1) {
//			// Just read the layout 1 for probabilities
//		} else if (fileLayout == Layout.layOut_2) {
//			//
//		}
		try {
			float[][] probabilities = readGenotypesFromVariant(variantDataPosition);
		} catch (IOException e) {
			e.printStackTrace();
		}
		return null;
	}

	private long getVariantDataPosition(GeneticVariant variant) {
		BgenixVariantQueryResult variantQueryResult = bgenixReader.getVariantsPostion(
				variant.getSequenceName(), variant.getStartPos());
		for (BgenixVariantQueryResult it = variantQueryResult; it.hasNext(); ) {
			BgenixVariantData variantData = it.next();
			if (variantData.getRsid().equals(variant.getPrimaryVariantId())) {
				return variantData.getFile_start_position();
			}
		}
		throw new GenotypeDataException(String.format("Variant with primary ID '%s' does not match any indexed variant",
				variant.getPrimaryVariantId()));
	}

	private Iterable<GeneticVariant> getGeneticVariants(BgenixVariantQueryResult variantQueryResult) {
		GeneticVariantRange.GeneticVariantRangeCreate variantRangeFactory = GeneticVariantRange.createRangeFactory();

		for (BgenixVariantQueryResult it = variantQueryResult; it.hasNext(); ) {
			BgenixVariantData variantData = it.next();
			long variantDataPosition = variantData.getFile_start_position();
			try {
				variantRangeFactory.addVariant(processVariantIdentifyingData(variantDataPosition));
			} catch (IOException e) {
				throw new GenotypeDataException(String.format(
						"Could not read variant data %s at position %d%n",
						variantData.getRsid(), variantDataPosition));
			}
		}
		return variantRangeFactory.createRange();
	}

	/**
	 * Method that returns the list of all possible allele combinations, for the given number of ploidity.
	 * These allele combinations, or genotypes represent the way in which genotypes are stored in a layout 2 BGEN file
	 * for a specific sample.
	 *
	 * For example, the method returns 4 combinations with 2 alleles and a ploidity of 2:
	 * [1, 1]
	 * [2, 1]
	 * [1, 2]
	 * [2, 2]
	 *
	 * @param numberOfAlleles The number of alleles for the variant to create combinations for.
	 * @param ploidy The ploidity of the sample to create combinations for.
	 * @return the combinations of alleles, genotypes, that should be stored in a BGEN file
	 * for an unphased variant and sample, in this specified order.
	 */
	private static List<List<Integer>> combinations(int numberOfAlleles, int ploidy) {
		// Construct nested lists
		List<List<Integer>> combinations = new ArrayList<>();
		// Get the combinations
		getCombinationsRecursively(combinations, new ArrayList<>(), numberOfAlleles, ploidy);
		return combinations;
	}

	/**
	 * Method that recursively fills a list of combinations.
	 *
	 * @param combinations The list of combinations to fill.
	 * @param combination The current combination that is being constructed.
	 * @param numberOfAlleles The number of different values to fit into a combinations.
	 * @param ploidy The size of a combination.
	 */
	private static void getCombinationsRecursively(List<List<Integer>> combinations, List<Integer> combination, int numberOfAlleles, int ploidy) {
		// If the combination is complete, the size of the combination equals the required size.
		// Add the combination and return
		if (combination.size() == ploidy) {
			combinations.add(0, combination);
			return;
		}
		// Get the first value from the current combination (which is the last one that is inserted),
		// if it exists, otherwise get the max value.
		int i = (combination.size() > 0) ? combination.get(0) : numberOfAlleles;
		// This prevents higher values than the previous value being inserted into the combination, which will cause
		// duplicate combinations (although in reverse order).

		// Loop through the possible values from high to low.
		for (; i > 0; i--) {
			// Copy the preliminary combination
			List<Integer> newCombination = new ArrayList<>(combination);
			// Add a new value to the combination
			newCombination.add(0, i);
			getCombinationsRecursively(combinations, newCombination, numberOfAlleles, ploidy);
		}
	}

//	/**
//	 * Converts a specified range of bits from a byte array into a 'long' value.
//	 *
//	 * @param bytes The byte array to get the result value from.
//	 * @param bitOffset The bit in the byte array to start reading from.
//	 * @param totalBitsToRead The total number of bits to read from the byte array.
//	 * @return The specified bits converted to a 'long' value.
//	 */
//	private static long getIntOfNBits(byte[] bytes, long bitOffset, int totalBitsToRead) {
//		// Get the byte to start reading from.
//		int byteOffset = Math.toIntExact(bitOffset / 8);
//		// Get the bit within the byte to start reading from.
//		int remainingBitOffset = Math.toIntExact(bitOffset % 8);
//		// Get the number of bits to read within the initial byte.
//		long nBitsToRead = 8 - remainingBitOffset;
//
//		// Define result value and a counter for the number of read bits.
//		long nReadBits = 0;
//		long value = 0;
//
//		System.out.println("bitOffset = " + bitOffset);
//
//		// First read the number of bits that should be read from the first byte (8 - remainingBitOffset).
//		// After that, if the condition still applies, read all bits of following bytes, placing these bits
//		// in before the previous read bits (<< readbits)
//		while (nReadBits + 8 < totalBitsToRead) {
//			// Read all bits in the current byte, moving them the number of already read bits to the left,
//			// and masking all bits that are not within the bits that have previously been read or should be read now.
//			System.out.println("nBitsToRead = " + nBitsToRead);
//			value |= bytes[byteOffset] >> (8 - nBitsToRead) & ((1L << (nBitsToRead + nReadBits)) - 1);
//
//			// Update bit reading numbers
//			nReadBits += nBitsToRead;
//			nBitsToRead = 8;
//
//			String s1 = String.format("%8s", Integer.toBinaryString(bytes[byteOffset] & 0xFF)).replace(' ', '0');
//			System.out.println("value = " + s1 + " | " + Long.toString(value,2));
//
//			// Update offset values
//			byteOffset += 1;
//			remainingBitOffset = 0;
//		}
//
//		// Calculate the number of bits that still have to be read.
//		nBitsToRead = totalBitsToRead - (nReadBits - remainingBitOffset);
//		System.out.println("nBitsToRead = " + nBitsToRead);
//
//		// First, move the bits that have to be read in the last byte to the right boundary of the byte
//		// (shifting amount defined by: 8 - nBitsToRead and mask the bits that should not be
//		// read.
//		long readBits = bytes[byteOffset] >> (8 - nBitsToRead) & ((1 << nBitsToRead) - 1);
//		// Secondly, move the read bits to the leftmost part of the result value, and mask all other bits.
//		value |= readBits << nReadBits & ((1L << totalBitsToRead) - 1);
//
//		String s1 = String.format("%8s", Integer.toBinaryString(bytes[byteOffset] & 0xFF)).replace(' ', '0');
//		System.out.println("value = " + s1 + " | " + Long.toString(value,2) + " | " + value);
//		return value;
//	}

	/**
	 * Convert 4 bytes to unsigned 32 bit int from index. Returns long since
	 * java does not have unsigned int
	 *
	 * https://stackoverflow.com/questions/13203426/convert-4-bytes-to-an-unsigned-32-bit-integer-and-storing-it-in-a-long
	 *
	 * @return
	 * @throws EOFException
	 * @throws IOException
	 */
	private long getUInt32(byte[] bytes, int startIndex) {
		long value = bytes[0 + startIndex] & 0xFF;
		value |= (bytes[1 + startIndex] << 8) & 0xFFFF;
		value |= (bytes[2 + startIndex] << 16) & 0xFFFFFF;
		value |= (bytes[3 + startIndex] << 24) & 0xFFFFFFFF;
		return value;
	}

	/**
	 * Convert 2 bytes to unsigned 16 bit int from start index.
	 *
	 * https://stackoverflow.com/questions/13203426/convert-4-bytes-to-an-unsigned-32-bit-integer-and-storing-it-in-a-long
	 *
	 * @return
	 * @throws EOFException
	 * @throws IOException
	 */
	private int getUInt16(byte[] bytes, int startIndex) {
		int value = bytes[0 + startIndex] & 0xFF;
		value |= (bytes[1 + startIndex] << 8) & 0xFFFF;
		return value;
	}

	/**
	 * Convert 1 byte to unsigned 8 bit int from start index.
	 *
	 * https://stackoverflow.com/questions/13203426/convert-4-bytes-to-an-unsigned-32-bit-integer-and-storing-it-in-a-long
	 *
	 * @return
	 * @throws EOFException
	 * @throws IOException
	 */
	private int getUInt8(byte[] bytes, int startIndex) {
		int value = bytes[0 + startIndex] & 0xFF;
		return value;
	}

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

		System.out.println("bitOffset = " + bitOffset);
		System.out.println("Total bytes min 1: " + totalBytesMin1);
		System.out.println("First byte: " + firstByteIndex);
		System.out.println("Index of first bit first byte: " + remainingBitOffset);
		System.out.println("Total bits: " + totalBitsToRead);
		System.out.println("Bits from last byte: " + nBitsFromLastByte);
		System.out.println("Last byte mask: " + Integer.toBinaryString(LAST_BYTE_MASK[nBitsFromLastByte]));
		System.out.println("First byte mask: " + Integer.toBinaryString(FIRST_BYTE_MASK[bitShiftAfterFirstByte]));
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
				System.out.println(String.format("%8s", Integer.toBinaryString(bytes[firstByteIndex] & 0xFF)).replace(' ', '0'));
				System.out.println(String.format("%8s", Integer.toBinaryString(bytes[firstByteIndex + 1] & 0xFF)).replace(' ', '0'));
				System.out.println(String.format("%8s", Integer.toBinaryString(bytes[firstByteIndex + 2] & 0xFF)).replace(' ', '0'));
				System.out.println(String.format("%8s", Integer.toBinaryString(bytes[firstByteIndex + 3] & 0xFF)).replace(' ', '0'));
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
}
