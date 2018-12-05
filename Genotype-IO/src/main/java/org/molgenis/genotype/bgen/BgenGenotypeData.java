/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.bgen;

import java.io.EOFException;
import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.zip.DataFormatException;
import java.util.zip.Inflater;
import com.facebook.presto.orc.zstd.ZstdDecompressor;
import java.util.List;
import java.util.Map;
import org.apache.log4j.Logger;
import org.molgenis.genotype.AbstractRandomAccessGenotypeData;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.GenotypeDataException;
import org.molgenis.genotype.Sample;
import org.molgenis.genotype.Sequence;
import org.molgenis.genotype.annotation.Annotation;
import org.molgenis.genotype.annotation.SampleAnnotation;
import org.molgenis.genotype.util.FixedSizeIterable;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variant.GeneticVariantMeta;
import org.molgenis.genotype.variant.GeneticVariantMetaMap;
import org.molgenis.genotype.variant.GenotypeRecord;
import org.molgenis.genotype.variant.ReadOnlyGeneticVariant;
import org.molgenis.genotype.variant.sampleProvider.CachedSampleVariantProvider;
import org.molgenis.genotype.variant.sampleProvider.SampleVariantsProvider;

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
	private final BgenixReader bgenixReader;
	private final long sampleCount;

	public BgenGenotypeData(File bgenFile, File sampleFile) throws IOException {
		this(bgenFile, sampleFile, 1000);
	}

	public BgenGenotypeData(File bgenFile, File sampleFile, int cacheSize) throws IOException {
		this(bgenFile, sampleFile, cacheSize, DEFAULT_MINIMUM_POSTERIOR_PROBABILITY_TO_CALL);
	}

	public BgenGenotypeData(File bgenFile, File sampleFile, int cacheSize, double minimumPosteriorProbabilityToCall) throws IOException {

		this.bgenFile = new RandomAccessFile(bgenFile, "r");
		
		
		
		if (this.bgenFile.read(byteArray4, 0, 4) != 4) {
			throw new GenotypeDataException("Error reading bgen file header. File is corrupt");
		}
		long snpOffset = getUInt32(byteArray4, 0);
		LOGGER.debug("SNP offset: " + snpOffset);

		if (this.bgenFile.read(byteArray4, 0, 4) != 4) {
			throw new GenotypeDataException("Error reading bgen file header. File is corrupt");
		}
		long headerSize = getUInt32(byteArray4, 0);
		LOGGER.debug("Header size (Lh): " + headerSize);
		if (headerSize > snpOffset) {
			throw new GenotypeDataException("Error reading bgen file header. Header information is bigger than expected offset.");
		}
		//add checks for header size and header size must be smaller than offset.
		if (this.bgenFile.read(byteArray4, 0, 4) != 4) {
			throw new GenotypeDataException("Error reading bgen file header. File is corrupt");
		}
		long variantCount = getUInt32(byteArray4, 0);
		LOGGER.debug("Number of SNPs: " + variantCount);

//               Not needed when using the bgenix indexing file.  
//		if (snpCount > Integer.MAX_VALUE) {
//			throw new GenotypeDataException("Found more than (2^31)-1 SNPs in bgen file. This is not supported");
//		}
		if (this.bgenFile.read(byteArray4, 0, 4) != 4) {
			throw new GenotypeDataException("Error reading bgen file header. File is corrupt");
		}
		sampleCount = getUInt32(byteArray4, 0);
		LOGGER.debug("Number of samples: " + sampleCount);

		if (sampleCount > Integer.MAX_VALUE) {
			throw new GenotypeDataException("Found more than (2^31)-1 samples in bgen file. This is not supported");
		}
		//Magic number is in the next four bytes but skipped over.

		//Skip over reserved and free area, 
		//this seek works because there are 4 (interesting) bytes 
		//in the beginning of the file and in the end of the header block.
		this.bgenFile.seek(headerSize);

		//Read flags
		if (this.bgenFile.read(byteArray4, 0, 4) != 4) {
			throw new GenotypeDataException("Error reading bgen file header. File is corrupt");
		}

		//Changed storing of version and data.
		switch (byteArray4[0] & 3) {
			case 0:
				LOGGER.debug("Genotype data is not compressed.");
				snpBlockRepresentation = blockRepresentation.compression_0;
				break;
			case 1:
				LOGGER.debug("Genotype data is zlib compressed");
				snpBlockRepresentation = blockRepresentation.compression_1;
				break;
			case 2:
				LOGGER.debug("Genotype data is zstd compressed");
				snpBlockRepresentation = blockRepresentation.compression_2;
				break;
			default:
				throw new GenotypeDataException("Invalid compression method, observed: " + (byteArray4[0] & 3));
		}
		byte byte_tmp = byteArray4[0];
		byte_tmp = (byte) (byte_tmp >> 2);
		switch (byte_tmp & 7) {
			case 1:
				LOGGER.debug("Genotype data is in layout 1 (BGEN 1.1)");
				fileLayout = Layout.layOut_1;
				break;
			case 2:
				LOGGER.debug("Genotype data is in layout 2 (BGEN 1.2 & 1.3)");
				fileLayout = Layout.layOut_2;
				break;
			default:
				throw new GenotypeDataException("Invalid layout, observed: " + (byteArray4[0] & 3));
		}

		if (snpBlockRepresentation.equals(blockRepresentation.compression_2) & fileLayout.equals(Layout.layOut_1)) {
			throw new GenotypeDataException("Invalid compression method for layout one observed. Trying to use ZSTD compression on layout one file, which is not supported.");
		}

		if ((byteArray4[3] & 128) == 128) {
			LOGGER.debug("SampleIdentifiers present in bgen file");
			sampleIdentifiersPresent = true;
		} else {
			LOGGER.debug("SampleIdentifiers not present in bgen file");
			sampleIdentifiersPresent = false;
		}

		if (sampleIdentifiersPresent) {
			// Proccess sample indentifier block
			if (this.bgenFile.read(byteArray4, 0, 4) != 4) {
				throw new GenotypeDataException("Error in sample identifier block. File is corrupt.");
			}
			long byteSizeSampleIds = getUInt32(byteArray4, 0);
			LOGGER.debug("LSI size: " + byteSizeSampleIds);
			if ((byteSizeSampleIds + headerSize) > snpOffset) {
				throw new GenotypeDataException("Error reading bgen file header. Combination of header & sample id information is bigger than expected offset.");
			}

			// Proccess sample indentifier block
			if (this.bgenFile.read(byteArray4, 0, 4) != 4) {
				throw new GenotypeDataException("Error in sample identifier block. File is corrupt.");
			}
			long sampleCount2 = getUInt32(byteArray4, 0);
			LOGGER.debug("Number of samples with sample-ids: " + sampleCount2);
			if (sampleCount2 != sampleCount) {
				throw new GenotypeDataException("Numer of samples in metadata and sample id data is not equal. File is corrupt.");
			}
			String[] sampleIds = new String[(int) sampleCount2];
			for (int i = 0; i < sampleCount2; i++) {
				if (this.bgenFile.read(byteArray2, 0, 2) != 2) {
					throw new GenotypeDataException("Error in sample Id. File is corrupt.");
				}
				int sampleIdLength = getUInt16(byteArray2, 0);
				byte[] sampleName = new byte[sampleIdLength];
				this.bgenFile.read(sampleName, 0, sampleIdLength);
				sampleIds[i] = new String(sampleName, CHARSET);
//				System.out.println("\t"+sampleIds[i]);
			}
		}
//		System.out.println("offset:"+ (snpOffset + 4));
		//Check if bgenix file is present.
		File bgenixFile = new File(bgenFile.getAbsolutePath() + ".bgi");
		System.out.println(bgenFile.getAbsolutePath() + ".bgi");

		long lastSnpStart = snpOffset + 4;

		if (!bgenixFile.exists()) {
			LOGGER.info("Creating bgenix file at: " + bgenixFile.getAbsolutePath());
			BgenixWriter b = new BgenixWriter(bgenixFile);
			createBgenixFile(bgenFile, b, lastSnpStart, (int) sampleCount, this.fileLayout, this.snpBlockRepresentation);
			b.finalizeIndex();
		}

		bgenixReader = new BgenixReader(bgenixFile);
		BgenixMetadata metadata = bgenixReader.getMetadata();
		if (metadata != null) {
			if (!metadata.getFileName().equals(bgenFile.getName())) {
				throw new GenotypeDataException("Sample name between bgenix and bgen is not equal. Invalid Bgen and Bgenix combination.");
			}

			if (metadata.getFileSize() != bgenFile.length()) {
				throw new GenotypeDataException("Number of expected bytes is different to the numer of observed bytes. Invalid Bgen and Bgenix combination.");
			}

			this.bgenFile.seek(0);
			byte[] firstBytes = new byte[1000];
			this.bgenFile.read(firstBytes, 0, 1000);
			if (!Arrays.equals(metadata.getFirst1000bytes(), firstBytes)) {
				throw new GenotypeDataException("First 1000 bytes of meta data and actual data are not equal. Invalid Bgen and Bgenix combination.");
			}
		}

		if (cacheSize > 0) {
			sampleVariantProvider = new CachedSampleVariantProvider(this, cacheSize);
		} else {
			sampleVariantProvider = this;
		}
		//Read the first snp to get into genotype-io.
//		readCompleteGeneticVariant(bgenFile, lastSnpStart, (int) sampleCount, this.fileLayout, this.snpBlockRepresentation);
	}

	private void createBgenixFile(File bgen, BgenixWriter b, long pointerFirstSnp, int nSamples, Layout fileLayout, blockRepresentation fileBlockRepresentation) throws IOException {
		this.bgenFile.seek(0);

		byte[] firstBytes = new byte[1000];
		this.bgenFile.read(firstBytes, 0, 1000);

		//Add current time in int.
		System.out.println((System.currentTimeMillis() / 1000L));
		BgenixMetadata m = new BgenixMetadata(bgen.getName(), (int) this.bgenFile.length(), (int) (bgen.lastModified() / 1000L), firstBytes, (System.currentTimeMillis() / 1000L));
		b.writeMetadata(m);

		//Loop through the start of the file
//		int stepToNextVariant = 0;
//		if (fileLayout.equals(layout.layOut_1)) {
//			stepToNextVariant += 4;
//		}
		long startSize = pointerFirstSnp;
		while ((startSize) < bgenFile.length()) {
			//Loop through variants.
//			long currentStart = pointerFirstSnp + stepToNextVariant;

			GeneticVariant var = readSnpInfo(startSize);
			long currentPointer = this.bgenFile.getFilePointer();
			long stepSize = DetermineStepSize(currentPointer, fileLayout, fileBlockRepresentation, nSamples);

			if (fileLayout.equals(fileLayout.layOut_2)) {
				stepSize += 4;
			}
//			this.bgenFile.seek(currentPointer);
			startSize = currentPointer + stepSize;

			b.addVariantToIndex(var, pointerFirstSnp, (int) stepSize, var.getVariantId().getPrimairyId());
		}
	}

	private GeneticVariant readSnpInfo(long filePointer) throws IOException {
		long lastSnpStart = filePointer;
		if (fileLayout == Layout.layOut_1) {
			lastSnpStart += 4;
		}
//		System.out.println("Seeking to:" +lastSnpStart);
		this.bgenFile.seek(lastSnpStart);
		//Not sure if we want to do the buffer search here. Or we might be able to take a smaller set.
		byte[] snpInfoBuffer = new byte[8096];
//		int snpInfoBufferSize = 
		this.bgenFile.read(snpInfoBuffer, 0, snpInfoBuffer.length);
		int snpInfoBufferPos = 0;

//		if (snpInfoBufferSize < 20) {
//			throw new GenotypeDataException("Error reading bgen snp data. File is corrupt");
//		}
		// Need to check that it is correct with the block infront of the snp id.
		ArrayList<String> snpIds = new ArrayList<>();
		int fieldLength = getUInt16(snpInfoBuffer, snpInfoBufferPos);
		LOGGER.debug("Snp ID length " + fieldLength);
		snpInfoBufferPos += 2;
		String snpId = new String(snpInfoBuffer, snpInfoBufferPos, fieldLength, CHARSET);
//		System.out.println(snpId);
		snpInfoBufferPos += fieldLength; // skip id length and snp id

		fieldLength = getUInt16(snpInfoBuffer, snpInfoBufferPos);
		LOGGER.debug("Snp RS length " + fieldLength);

		snpInfoBufferPos += 2;
		String snpRs = new String(snpInfoBuffer, snpInfoBufferPos, fieldLength, CHARSET);
//		System.out.println(snpRs);
		snpInfoBufferPos += fieldLength;
		snpIds.add(snpRs);
		snpIds.add(snpId);

		fieldLength = getUInt16(snpInfoBuffer, snpInfoBufferPos);
		snpInfoBufferPos += 2;
		String seqName = new String(snpInfoBuffer, snpInfoBufferPos, fieldLength, CHARSET);
//		System.out.println(seqName);
		snpInfoBufferPos += fieldLength;

		long snpPosLong = getUInt32(snpInfoBuffer, snpInfoBufferPos);
		snpInfoBufferPos += 4;
		if (snpPosLong > Integer.MAX_VALUE) {
			throw new GenotypeDataException("SNP pos larger than (2^31)-1 not supported");
		}

		int snpPos = (int) snpPosLong;

//		System.out.println("SNP pos " + snpPos);
		int numberOfAlleles = 2;
		if (fileLayout.equals(Layout.layOut_2)) {
			numberOfAlleles = getUInt16(snpInfoBuffer, snpInfoBufferPos);
			snpInfoBufferPos += 2;
//			System.out.println("SNP Alleles " + numberOfAlleles);
		}

		ArrayList<String> alleles = new ArrayList<>();
		for (int i = 0; i < numberOfAlleles; i++) {
			long fieldLengthLong = getUInt32(snpInfoBuffer, snpInfoBufferPos);
			snpInfoBufferPos += 4;
			if (fieldLengthLong > Integer.MAX_VALUE) {
				throw new GenotypeDataException("SNP with allele longer than (2^31)-1 characters not supported");
			}
			String a = new String(snpInfoBuffer, snpInfoBufferPos, (int) fieldLengthLong, CHARSET);
			snpInfoBufferPos += ((int) fieldLengthLong);
			alleles.add(a);
		}

//		System.out.println(alleles.toString());
//			System.out.println("Location where genotypes start: "+this.bgenFile.getFilePointer());
		this.bgenFile.seek(snpInfoBufferPos + lastSnpStart);

		GeneticVariant var = ReadOnlyGeneticVariant.createVariant(GP_VARIANT, snpIds, snpPos, seqName, sampleVariantProvider, alleles, alleles.get(0));

		return var;
	}

	private void readGenotypesFromVariant(long filePointer) throws IOException {
		this.bgenFile.seek(filePointer);

		//Not sure if we want to do the buffer search here. Or we might be able to take a smaller set.
		byte[] snpInfoBuffer = new byte[8096];
//		int snpInfoBufferSize = 
		this.bgenFile.read(snpInfoBuffer, 0, snpInfoBuffer.length);
		int snpInfoBufferPos = 0;

		if (fileLayout == fileLayout.layOut_1) {
			long snpBlockSize;
			if (snpBlockRepresentation.equals(blockRepresentation.compression_1)) {
				snpBlockSize = getUInt32(snpInfoBuffer, snpInfoBufferPos);
				snpBlockSize += 4;
			} else {
				snpBlockSize = 6 * sampleCount;
			}
			System.out.println("Snp block size: " + snpBlockSize);
			byte[] snpBlockData = new byte[6 * (int) sampleCount];
			if (snpBlockRepresentation.equals(blockRepresentation.compression_1)) {
				gzipInflater.setInput(snpInfoBuffer, snpInfoBufferPos + 4, (int) snpBlockSize - 4);
				try {
					gzipInflater.inflate(snpBlockData);
				} catch (DataFormatException ex) {
					throw new GenotypeDataException("Error decompressing bgen data", ex);
				}
				gzipInflater.reset();

				//Here we can parse from the inflater block for all the samples.
				System.out.println(getUInt16(snpBlockData, 0) / 32768f + " " + getUInt16(snpBlockData, 2) / 32768f + " " + getUInt16(snpBlockData, 4) / 32768f);
			} else {
				//Here we can directly parse from the snpInfoBuffer.
			}

		} else if (fileLayout.equals(fileLayout.layOut_2)) {
			long snpBlockSize;
			long snpBlockSizeDecompressed;
			if (!snpBlockRepresentation.equals(blockRepresentation.compression_0)) {
				snpBlockSize = getUInt32(snpInfoBuffer, snpInfoBufferPos);
				snpInfoBufferPos += 4;
				snpBlockSizeDecompressed = getUInt32(snpInfoBuffer, snpInfoBufferPos);
				snpInfoBufferPos += 4;
			} else {
				snpBlockSize = getUInt32(snpInfoBuffer, snpInfoBufferPos);
				snpInfoBufferPos += 4;
				snpBlockSizeDecompressed = snpBlockSize;
			}
			System.out.println("Snp block size: " + snpBlockSize);
			System.out.println("Snp block size decompressed: " + snpBlockSizeDecompressed);

			byte[] snpBlockData = new byte[(int) snpBlockSizeDecompressed];
			switch (snpBlockRepresentation) {

				case compression_1:
					gzipInflater.setInput(snpInfoBuffer, snpInfoBufferPos, (int) snpBlockSize - 4);
					try {
						gzipInflater.inflate(snpBlockData);
					} catch (DataFormatException ex) {
						throw new GenotypeDataException("Error decompressing bgen data", ex);
					}
					gzipInflater.reset();

				//At genotype / haplotype data
				case compression_2:
//					zstdInflater.decompress(snpInfoBuffer, snpInfoBufferPos,(int) snpBlockSize - 4,snpBlockData, 0, (int) snpBlockSizeDecompressed);
				//Is this enough?

				default:
					break;
//					snpInfoBufferPos = ;
				//Not compressed.
			}
			int blockBuffer = 0;
			//must equal data before.
			int numberOfIndividuals = (int) getUInt32(snpBlockData, blockBuffer);
			System.out.println("Number of individuals: " + numberOfIndividuals);
			blockBuffer += 4;
			//must equal data before.
			int numberOfAlleles = (int) getUInt16(snpBlockData, blockBuffer);
			System.out.println("Number of Alleles: " + numberOfAlleles);
			blockBuffer += 2;

			int minPloidy = getUInt8(snpBlockData, blockBuffer);
			System.out.println("Min ploidy: " + minPloidy);
			blockBuffer += 1;
			int maxPloidy = getUInt8(snpBlockData, blockBuffer);
			System.out.println("Max ploidy: " + maxPloidy);
			blockBuffer += 1;
			for (int i = 0; i < numberOfIndividuals; i++) {
				//Here we need to handle missing ploidity.
				//Missingness is encoded by the most significant bit; thus a value of 1 for the most significant bit indicates that no probability data is stored for this sample.
//				System.out.println("ploidity: " + getUInt8(snpBlockData, blockBuffer));
				blockBuffer += 1;
			}
			int phased = getUInt8(snpBlockData, blockBuffer);
			System.out.println("phased: " + phased);
			blockBuffer += 1;
			if (phased > 1) {
				throw new GenotypeDataException("Bgen file format error. Unsupported value for the phased flag observed.");
			}
			int bitProbabilityRepresentation = getUInt8(snpBlockData, blockBuffer);
			System.out.println("Bit representation of probability: " + bitProbabilityRepresentation);
			blockBuffer += 1;

			if (phased == 1) {

			} else {

			}

		}
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
	private Long DetermineStepSize(long filePointer, Layout currentFileLayout, blockRepresentation currentBlockRepresentation, int sampleCount) throws IOException {
		this.bgenFile.seek(filePointer);

		//Not sure if we want to do the buffer search here. Or we might be able to take a smaller set.
		byte[] snpInfoBuffer = new byte[8096];
//		int snpInfoBufferSize = 
		this.bgenFile.read(snpInfoBuffer, 0, snpInfoBuffer.length);
		int snpInfoBufferPos = 0;

		Long snpBlockSize = null;
		long snpBlockSizeDecompressed;

		if (currentFileLayout == fileLayout.layOut_1) {
			if (currentBlockRepresentation.equals(blockRepresentation.compression_1)) {
				snpBlockSize = getUInt32(snpInfoBuffer, snpInfoBufferPos);
				snpBlockSize += 4;
			} else {
				snpBlockSize = new Long((long) (6 * sampleCount));
			}
			System.out.println("Snp block size: " + snpBlockSize);

		} else if (currentFileLayout == fileLayout.layOut_2) {

			if (!currentBlockRepresentation.equals(blockRepresentation.compression_0)) {
				snpBlockSize = getUInt32(snpInfoBuffer, snpInfoBufferPos);
				snpInfoBufferPos += 4;
				snpBlockSizeDecompressed = getUInt32(snpInfoBuffer, snpInfoBufferPos);
				snpInfoBufferPos += 4;
			} else {
				snpBlockSize = getUInt32(snpInfoBuffer, snpInfoBufferPos);
				snpInfoBufferPos += 4;
				snpBlockSizeDecompressed = snpBlockSize;
			}
			System.out.println("Snp block size: " + snpBlockSize);
			System.out.println("Snp block size decompressed: " + snpBlockSizeDecompressed);
		}
		return snpBlockSize;
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
