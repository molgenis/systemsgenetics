/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.bgen;

import java.io.*;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.channels.Channel;
import java.nio.channels.Channels;
import java.nio.channels.WritableByteChannel;
import java.nio.charset.Charset;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.HashMap;

import org.apache.log4j.Logger;
import org.molgenis.genotype.GenotypeData;
import org.molgenis.genotype.GenotypeWriter;
import org.molgenis.genotype.Sample;
import org.molgenis.genotype.oxford.OxfordSampleFileWriter;
import org.molgenis.genotype.util.Utils;
import org.molgenis.genotype.variant.NotASnpException;

/**
 *
 * @author patri
 */
public class BgenGenotypeWriter implements GenotypeWriter {

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
			int sampleIdentifierLength = id.length();
			byte[] idAsBytes = id.getBytes(CHARSET);

			// 2 bytes with the length of identifier of sample 1
			if (idAsBytes.length > 2) {
				// The length of the sample identifier exceeds two bytes. Give a warning and truncate.
				LOGGER.warn(String.format(
						"Length of variant identifier %.16s... exceeds the maximum length of 2 bytes. " +
								"(%d vs %d respectively) %nSample id gets truncated.",
						id, sampleIdentifierLength, Short.MAX_VALUE));
				// ...truncate...
				sampleIdentifierLength = Short.MAX_VALUE;
				idAsBytes = Arrays.copyOfRange(idAsBytes, 0, 3);
			}
			// Initialize a byte buffer with the correct size.
			ByteBuffer sampleBuffer = ByteBuffer.allocate(2 + idAsBytes.length)
					.order(ByteOrder.LITTLE_ENDIAN);
			// Put both the length of the sample identifier as a short
			// and the identifier itself int the byte array
			sampleBuffer.putShort((short) sampleIdentifierLength);
			sampleBuffer.put(idAsBytes);
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

		return null;
	}

}
