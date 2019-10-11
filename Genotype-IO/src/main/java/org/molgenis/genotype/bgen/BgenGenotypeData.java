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
import org.molgenis.genotype.variant.sampleProvider.CachedSampleVariantProvider;
import org.molgenis.genotype.variant.sampleProvider.SampleVariantsProvider;

import java.io.EOFException;
import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
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
	private static final GeneticVariantMeta GP_VARIANT = GeneticVariantMetaMap.getGeneticVariantMetaGp();
	private final SampleVariantsProvider sampleVariantProvider;
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

//			readGenotypesFromVariant(currentPointer);

			// Get the position in the file to start reading the next variant.
			variantReadingPosition = currentPointer + variantGenotypeDataBlockInfo.getBlockLengthHeaderInclusive();
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

	private void readGenotypesFromVariant(long filePointer) throws IOException {
		this.bgenFile.seek(filePointer);

		//Not sure if we want to do the buffer search here. Or we might be able to take a smaller set.
		byte[] snpInfoBuffer = new byte[8096];
//		int snpInfoBufferSize = 
		this.bgenFile.read(snpInfoBuffer, 0, snpInfoBuffer.length);
		int snpInfoBufferPos = 0;

		float[][] probabilities = new float[getSamples().size()][3];

		// Read
		if (fileLayout == Layout.layOut_1) {
			VariantGenotypeDataBlockInfo blockInfo = getVariantGenotypeDataBlockInfoForLayoutOne(
					snpInfoBuffer,
					snpInfoBufferPos);

            byte[] snpBlockData = getDecompressedBlockData(snpInfoBuffer, blockInfo);
			for (int sampleIndex = 0; sampleIndex < sampleCount; sampleIndex++) {
				float[] sampleProbabilities = new float[3];

//				System.out.println(getUInt16(snpBlockData, 0) / 32768f + " " + getUInt16(snpBlockData, 2) / 32768f + " " + getUInt16(snpBlockData, 4) / 32768f);
				sampleProbabilities[0] = getUInt16(snpBlockData, 0) / 32768f;
				sampleProbabilities[1] = getUInt16(snpBlockData, 2) / 32768f;
				sampleProbabilities[2] = getUInt16(snpBlockData, 4) / 32768f;
				probabilities[sampleIndex] = sampleProbabilities;
			}

        } else if (fileLayout.equals(Layout.layOut_2)) {
			VariantGenotypeDataBlockInfo blockInfo = determineVariantGenotypeDataBlockSizeForLayoutTwo(
					snpInfoBuffer,
					snpInfoBufferPos);

            System.out.println("Snp block size: " + blockInfo.getBlockLength());
			System.out.println("Snp block size decompressed: " + blockInfo.getDecompressedBlockLength());

			// Get the decompressed variant data of length D
            byte[] variantBlockData = getDecompressedBlockData(
                    snpInfoBuffer,
					blockInfo.getBlockOffset(),
                    blockInfo.getBlockLength(),
                    (int) blockInfo.getDecompressedBlockLength());

            int blockBufferOffset = 0;
			//must equal data before.
			if ((int) getUInt32(variantBlockData, blockBufferOffset) != sampleCount) {
                throw new GenotypeDataException(
                        "BGEN file format error. Sample counts do not match.");
            }
			System.out.println("Number of individuals: " + sampleCount);
			blockBufferOffset += 4;

			//must equal data before.
			int numberOfAlleles = (int) getUInt16(variantBlockData, blockBufferOffset);
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
			List<Integer> ploidity = new ArrayList<>();
			ReadPloidity(variantBlockData, isMissing, ploidity);
			// Add the number of individuals to the buffer.
			blockBufferOffset += sampleCount;

			boolean phased = isPhased(variantBlockData, blockBufferOffset);
			blockBufferOffset += 1;

			int probabilitiesLengthInBits = getUInt8(variantBlockData, blockBufferOffset);
			System.out.println("Bit representation of probability: " + probabilitiesLengthInBits);
			blockBufferOffset += 1;

			byte[] probabilitiesArray = Arrays.copyOfRange(variantBlockData, blockBufferOffset,
					variantBlockData.length - blockBufferOffset);

			if (phased) {
				readHaplotypeProbabilities(
						probabilitiesArray,
						probabilitiesLengthInBits,
						numberOfAlleles,
						isMissing, ploidity);
			} else {
				readGenotypeProbabilities(
						probabilitiesArray,
						probabilitiesLengthInBits,
						isMissing, ploidity);
			}

		}
	}

	private void readGenotypeProbabilities(
			byte[] probabilities,
			int probabilitiesLengthInBits,
			List<Boolean> isMissing,
			List<Integer> ploidity) {

	}

	private void readHaplotypeProbabilities(
			byte[] probabilitiesArray,
			int probabilitiesLengthInBits,
			int numberOfAlleles, List<Boolean> isMissing,
			List<Integer> ploidity) {

    	// Define an array consisting of an array of posterior probabilities for each genotype
		float[][] probabilities = new float[getSamples().size()][3];


		// Each probability is stored in B bits.
		// Values are interpreted by linear interpolation between 0 and 1;
		// value b corresponds to probability b / ((2^B)-1).

		for (int sampleIndex = 0; sampleIndex < sampleCount; sampleIndex++) {

			float[] sampleProbabilities = new float[3];
			// Calculate the probability for homozygous genotype 'AA',
			// the probability for 'AB', and the probability for 'BB'

			// Read Z * (K-1) probabilities of length B.
			Integer ploidy = ploidity.get(sampleIndex);
			for (int i = 0; i < ploidy; i++) {
				for (int j = 0; j < numberOfAlleles - 1; j++) {
					// Get the P_ij probability
					int probability = 0;

				}
				// Calculate probability of Kth allele (P_iK)
			}

		}
	}

	/**
	 * Read the ploidity and missingness.
	 *
	 * @param snpBlockData The byte buffer array to read the genotype data.
	 * @param isMissing A list of booleans representing if the corresponding sample has missing probabilities
	 * @param ploidity A list representing the ploidy of corresponding samples.
	 */
	private void ReadPloidity(byte[] snpBlockData, List<Boolean> isMissing, List<Integer> ploidity) {
		for (int i = 0; i < sampleCount; i++) {
			// Here we need to handle ploidity and missigness of prababilities.
			// Missingness is encoded by the most significant bit (MSB);
			// thus a value of 1 for the most significant bit
			// indicates that no probability data is stored for this sample.
			// Ploidity is encoded by the 6 least significant bits.

			// Get the MSB and check if the result equals 128.
			isMissing.add((snpBlockData[i] & 0x80) == 128);
			// Get the 6 least significant bits.
			ploidity.add(snpBlockData[i] & 0x3F);

//				String s1 = String.format("%8s", Integer.toBinaryString(snpBlockData[i] & 0xFF)).replace(' ', '0');
//                System.out.println("value = " + s1 + " | " + value + " | " + isMissing);
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
	 * Method for obtaining decompressed data for layout 1.
	 *
	 * @param snpInfoBuffer The input byte array to decompress.
	 * @param blockInfo The info of the variant genotype data block.
	 */
	private byte[] getDecompressedBlockData(byte[] snpInfoBuffer, VariantGenotypeDataBlockInfo blockInfo) {

		// Initialize byte array of length 6 times the number of samples.
        byte[] variantDataBlock = new byte[6 * (int) sampleCount];
		long variantBlockSize = blockInfo.getBlockLength();
		int variantBlockOffset = blockInfo.getBlockOffset();
		if (snpBlockRepresentation.equals(blockRepresentation.compression_1)) {
			// GZIP decompression
            decompressVariantBlockGzip(snpInfoBuffer,
					variantBlockOffset,
                    (int) variantBlockSize,
                    variantDataBlock);
		} else {
			// No decompression, copy variant data block.
            variantDataBlock = Arrays.copyOfRange(
                    snpInfoBuffer,
					variantBlockOffset,
                    (int) (variantBlockOffset + variantBlockSize));
        }
        return variantDataBlock;
    }

	/**
	 * Method for obtaining decompressed data for layout 2.
	 *
	 * @param snpInfoBuffer The input byte array to decompress.
	 * @param blockOffset The offset of the block size from the start of the input byte array.
	 * @param variantBlockSize The uncompressed block size.
	 * @param snpBlockSizeDecompressed The block size after
	 */
    private byte[] getDecompressedBlockData(byte[] snpInfoBuffer, int blockOffset, long variantBlockSize, int snpBlockSizeDecompressed) {
        // Initialize byte array.
        byte[] snpBlockData = new byte[snpBlockSizeDecompressed];
        switch (snpBlockRepresentation) {

            case compression_1:
                decompressVariantBlockGzip(snpInfoBuffer,
						blockOffset,
                        (int) variantBlockSize,
                        snpBlockData);
                break;

            //At genotype / haplotype data
            case compression_2:
                zstdInflater.decompress(snpInfoBuffer, blockOffset, (int) variantBlockSize, snpBlockData, 0, snpBlockSizeDecompressed);
                //Is this enough?
                break;

            default:
                snpBlockData = Arrays.copyOfRange(
                        snpInfoBuffer,
						blockOffset,
                        (int) (blockOffset + variantBlockSize));
                break;
//					snpInfoBufferPos = ;
            //Not compressed.
        }
        return snpBlockData;
    }

	/**
	 * Method that decompresses data using the gzipInflater.
	 *
	 * @param snpInfoBuffer The input byte array to decompress.
	 * @param blockOffset The offset of the block size from the start of the input byte array.
	 * @param snpBlockSize The uncompressed block size.
	 * @param snpBlockData The decompressed output byte array.
	 */
    private void decompressVariantBlockGzip(
    		byte[] snpInfoBuffer, int blockOffset, int snpBlockSize, byte[] snpBlockData) {

    	// Set the input for the gzip inflater.
		gzipInflater.setInput(snpInfoBuffer, blockOffset, snpBlockSize);

		// Try to decompress the data.
		try {
			gzipInflater.inflate(snpBlockData);
		} catch (DataFormatException e) {
			throw new GenotypeDataException("Error decompressing bgen data", e);
		}
		gzipInflater.reset();
	}

//	private void readCompleteGeneticVariant(File bgen, long lastSnpStart, int sampleCount, layout currentFileLayout, blockRepresentation currentBlockRepresentation) throws IOException {
//		//Binary index writer.
////		for (int snpI = 0; snpI < variantCount; ++snpI) {
//
//		long currentPointer= this.bgenFile.getFilePointer();
//		long stepSize=0L;
//		
//		while((currentPointer+stepSize)<bgen.length()){
//			GeneticVariant var = readSnpInfo((currentPointer+stepSize));
//			System.out.println(var.getAllIds().get(0));
//			currentPointer= this.bgenFile.getFilePointer();
//			
//			stepSize = DetermineStepSize(currentPointer, currentFileLayout, currentBlockRepresentation, sampleCount);
//			if(currentFileLayout.equals(layout.layOut_2)){
//				stepSize+=4;
//			}
//			this.bgenFile.seek(currentPointer);
//			//Here we can process actual genotype info.
//			
//			
////			System.out.println(bgen.length() - (currentPointer+stepSize));
//			System.out.println("Step to skip, to next variant: " + (int) stepSize);
//		}
//
//	}

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
		throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
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
		throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
	}

	@Override
	public List<String> getSeqNames() {
	    // Corresponds to
		throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
	}

	@Override
	public Iterable<Sequence> getSequences() {
		throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
	}

	@Override
	public Iterable<GeneticVariant> getVariantsByPos(String seqName, int startPos) {
		throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
	}

	@Override
	public Iterable<GeneticVariant> getSequenceGeneticVariants(String seqName) {
		BgenixVariantQueryResult bgenixVariants = bgenixReader.getVariantsChromosome(seqName);

		return null;
	}

	@Override
	public Iterable<GeneticVariant> getVariantsByRange(String seqName, int rangeStart, int rangeEnd) {
		throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
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
		throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
	}

	@Override
	public int getSampleVariantProviderUniqueId() {
		throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
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
		throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
	}

	private long getUIntUpTo32Bits(byte[] bytes, int bitOffset, int totalBitsToRead) {
		int byteOffset = bitOffset / 8;
		int remainingBitOffset = bitOffset % 8;

		int numberOfBytes = totalBitsToRead / 8;
		int remainingBits = totalBitsToRead % 8;

		int bitsToRead = 8 - remainingBitOffset;
		int readBits = 0;

		int value = 0;

		if (bitsToRead <= totalBitsToRead) {
			value |= bytes[byteOffset] & ((1 << bitsToRead) - 1);
			byteOffset += 1;
		}

//		if (11)

		value |= bytes[byteOffset] >> (8 - bitsToRead) & ((1 << bitsToRead)-1);
		return value;
	}

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
}
