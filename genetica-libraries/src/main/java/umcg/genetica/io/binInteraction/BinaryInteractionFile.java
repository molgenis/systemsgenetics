package umcg.genetica.io.binInteraction;

import com.google.common.io.CountingInputStream;
import java.io.BufferedInputStream;
import java.io.DataInputStream;
import java.io.EOFException;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.HashSet;
import org.apache.commons.lang.NotImplementedException;

/**
 *
 * @author Patrick Deelen
 */
public class BinaryInteractionFile {

	private static final int GENERAL_HEADER_SIZE = 10;//bytes
	private static final byte MAGIC_1 = 81;
	private static final byte MAGIC_2 = 73;
	private static final byte MAJOR_VERSION = 1;
	private static final byte MINOR_VERSION = 0;
	private final File interactionFile;
	private final boolean readOnly;
	private final BinaryInteractionCohort[] cohorts;
	private final BinaryInteractionGene[] genes;
	private final BinaryInteractionVariant[] variants;
	private final String[] covariats;
	private final String[] chrDictionary;
	private final String[] alleleDictionary;
	private final long timeStamp;
	private final boolean allCovariants;
	private final boolean metaAnalysis;
	private final boolean normalQtlStored;
	private final String fileDescription;
	private final int interactions;
	private final int startQtlBlock;
	private final int startInteractionBlock;

	protected BinaryInteractionFile(File interactionFile, boolean readOnly, BinaryInteractionCohort[] cohorts, BinaryInteractionGene[] genes, BinaryInteractionVariant[] variants, String[] covariats, String[] chrDictionary, String[] alleleDictionary, long timeStamp, boolean allCovariants, boolean metaAnalysis, boolean normalQtlStored, String fileDescription, int interactions, int startQtlBlock, int startInteractionBlock) {
		this.interactionFile = interactionFile;
		this.readOnly = readOnly;
		this.cohorts = cohorts;
		this.genes = genes;
		this.variants = variants;
		this.covariats = covariats;
		this.chrDictionary = chrDictionary;
		this.alleleDictionary = alleleDictionary;
		this.timeStamp = timeStamp;
		this.allCovariants = allCovariants;
		this.metaAnalysis = metaAnalysis;
		this.normalQtlStored = normalQtlStored;
		this.fileDescription = fileDescription;
		this.interactions = interactions;
		this.startQtlBlock = startQtlBlock;
		this.startInteractionBlock = startInteractionBlock;
	}

	public static BinaryInteractionFile createNew() {
		//todo
		throw new NotImplementedException();
	}

	public static BinaryInteractionFile load(File interactionFile) throws FileNotFoundException, IOException, BinaryInteractionFileException {
		return load(interactionFile, true);
	}

	public static BinaryInteractionFile load(File interactionFile, boolean readOnly) throws FileNotFoundException, IOException, BinaryInteractionFileException {

		BinaryInteractionFileConstructorBuilder builder = new BinaryInteractionFileConstructorBuilder();

		builder.setInteractionFile(interactionFile);
		builder.setReadOnly(readOnly);

		CountingInputStream inputStreamCounted = new CountingInputStream(new BufferedInputStream(new FileInputStream(interactionFile)));
		DataInputStream inputStream = new DataInputStream(inputStreamCounted);

		try {

			inputStream.read();
			//First parse heaeder		
			if (inputStream.readByte() != MAGIC_1 || inputStream.readByte() != MAGIC_2) {
				throw new BinaryInteractionFileException("No a valid interaction file");
			}

			if (inputStream.readByte() != MAJOR_VERSION || inputStream.readByte() != MINOR_VERSION) {
				throw new BinaryInteractionFileException("Interaction file version not supported");
			}

			builder.setTimeStamp(inputStream.readLong());
			builder.setAllCovariants(inputStream.readBoolean());
			builder.setMetaAnalysis(inputStream.readBoolean());
			builder.setNormalQtlStored(inputStream.readBoolean());

			inputStream.skip(4);//Skip reserved

			builder.setFileDescription(readString(inputStream));

			final int chortsCount = inputStream.readInt();
			final int chrsCount = inputStream.readInt();
			final int allelesCount = inputStream.readInt();
			final int genesCount = inputStream.readInt();
			final int snpsCount = inputStream.readInt();
			final int covariatsCount = inputStream.readInt();

			builder.setInteractions(inputStream.readInt());
			
			builder.setCohorts(readCohorts(inputStream, chortsCount));
			builder.setChrDictionary(readDictionary(inputStream, chrsCount, "Chromosomes"));
			builder.setChrDictionary(readDictionary(inputStream, allelesCount, "Alleles"));
			

			return builder.createBinaryInteractionFile();

		} catch (EOFException ex) {
			throw new BinaryInteractionFileException("Error parsing header, unexpected EOF", ex);
		}

	}

	private static BinaryInteractionCohort[] readCohorts(DataInputStream inputStream, int cohortsCount) throws EOFException, IOException, BinaryInteractionFileException {

		final BinaryInteractionCohort[] cohorts = new BinaryInteractionCohort[cohortsCount];
		final HashSet<String> cohortNames = new HashSet<String>(cohortsCount);

		for (int i = 0; i < cohortsCount; ++i) {

			String name = readString(inputStream);

			if (!cohortNames.add(name)) {
				throw new BinaryInteractionFileException("Cohort names must be unique. " + name + " has been found multiple times");
			}

			cohorts[i] = new BinaryInteractionCohort(name, inputStream.readInt());

		}
		return cohorts;
	}

	private static String[] readDictionary(DataInputStream inputStream, int size, String name) throws BinaryInteractionFileException, EOFException, IOException {
		
		final String[] dictrionary = new String[size];
		final HashSet<String> cohortNames = new HashSet<String>(size);

		for (int i = 0; i < size; ++i) {

			dictrionary[i] = readString(inputStream);
			if (!cohortNames.add(dictrionary[i])) {
				throw new BinaryInteractionFileException(name + " entries must be unique. " + dictrionary[i] + " has been found multiple times");
			}

		}
		return dictrionary;
		
	}

	/**
	 * Will first read int to determine length. Will then read length chars and
	 * convert to string.
	 *
	 * @param inputStream
	 * @return
	 */
	private static String readString(DataInputStream inputStream) throws EOFException, IOException {

		int lenght = inputStream.readInt();
		char[] chars = new char[lenght];

		for (int i = 0; i < lenght; ++i) {
			chars[i] = inputStream.readChar();
		}

		return new String(chars);

	}
}
