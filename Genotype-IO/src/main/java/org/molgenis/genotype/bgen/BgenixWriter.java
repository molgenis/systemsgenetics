/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.bgen;

import java.io.File;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.SQLException;
import java.util.List;
import org.molgenis.genotype.GenotypeDataException;
import org.molgenis.genotype.variant.GeneticVariant;

/**
 *
 * @author Patrick Deelen
 */
public class BgenixWriter {

	private final Connection dbConnection;
	private final PreparedStatement addVariantStatement;
	private int addVariantBufferCounter = 0;
	private boolean finalized;
	private static final int ADD_BUFFER_MAX = 5000;

	public BgenixWriter(File bgenixFile) {

		if (bgenixFile.exists()) {
			throw new GenotypeDataException("Bgen index file already exisist: " + bgenixFile.getAbsolutePath());
		}

		try {
			dbConnection = DriverManager.getConnection("jdbc:sqlite: " + bgenixFile.getPath());
			dbConnection.setAutoCommit(false);
		} catch (SQLException ex) {
			throw new GenotypeDataException("Unable to create bgenix file. Error: " + ex.getMessage(), ex);
		}

		try {

			dbConnection.createStatement().execute("CREATE TABLE Variant ( chromosome TEXT NOT NULL, position INT NOT NULL, rsid TEXT NOT NULL, number_of_alleles INT NOT NULL, allele1 TEXT NOT NULL, allele2 TEXT NULL, file_start_position INT NOT NULL, size_in_bytes INT NOT NULL, PRIMARY KEY (chromosome, position, rsid, allele1, allele2, file_start_position )) WITHOUT ROWID;");

			addVariantStatement = dbConnection.prepareStatement("INSERT INTO Metadata(chromosome, position, rsid, number_of_alleles, allele1, allele2, file_start_position, size_in_bytes) VALUES(?,?,?,?,?,?,?,?)");

		} catch (SQLException ex) {
			throw new GenotypeDataException("Unable to initialize bgenix file. Error: " + ex.getMessage(), ex);
		}

	}

	/**
	 *
	 *
	 * @param variant
	 * @param startPos
	 * @param sizeInBytesInBgenFile
	 * @param variantId variant ID as is written to bgen file. Variant ID is
	 * requered by bgen so might be different from the null ID that can be
	 * sorted in Genetic variant object
	 */
	public synchronized void addVariantToIndex(GeneticVariant variant, long startPos, int sizeInBytesInBgenFile, String variantId) {

		if (finalized) {
			throw new GenotypeDataException("Cannot add variants to bgenix after finalizing");
		}

		++addVariantBufferCounter;

		try {

			addVariantStatement.setString(1, variant.getSequenceName());
			addVariantStatement.setInt(2, variant.getStartPos());
			addVariantStatement.setString(3, variantId);

			List<String> alleles = variant.getVariantAlleles().getAllelesAsString();
			
			addVariantStatement.setInt(4,alleles.size());
			addVariantStatement.setString(5, alleles.get(0));
			addVariantStatement.setString(6, alleles.size() > 1 ? alleles.get(1) : null);
			
			addVariantStatement.setLong(7, startPos);
			addVariantStatement.setInt(8, sizeInBytesInBgenFile);

			addVariantStatement.addBatch();
			
		} catch (SQLException ex) {
			throw new GenotypeDataException("Unable to add variant to bgenix write buffer. Error: " + ex.getMessage(), ex);
		}

		if (addVariantBufferCounter > ADD_BUFFER_MAX) {

			try {
				addVariantStatement.executeBatch();
				dbConnection.commit();
			} catch (SQLException ex) {
				throw new GenotypeDataException("Unable to write variants from buffer to bgenix file. Error: " + ex.getMessage(), ex);
			}
			addVariantBufferCounter = 0;
		}

	}

	/**
	 * Finalize writing variants and meta data to bgenix file.
	 *
	 * @return db connection for use as index
	 */
	public BgenixReader finalizeIndex() {

		try {

			if (addVariantBufferCounter > 0) {
				addVariantStatement.executeBatch();
				dbConnection.commit();
			}

			dbConnection.setAutoCommit(true);//this is the default. Not sure if needed to set to true but just to be sure
			return new BgenixReader(dbConnection);

		} catch (SQLException ex) {
			throw new GenotypeDataException("Unable to finalize bgenix file. Error: " + ex.getMessage(), ex);
		}

	}

	public void writeMetadata(BgenixMetadata metadata) {

		try {

			dbConnection.createStatement().execute("CREATE TABLE Metadata (filename TEXT NOT NULL, file_size INT NOT NULL, last_write_time INT NOT NULL, first_1000_bytes BLOB NOT NULL, index_creation_time INT NOT NULL);");

			PreparedStatement metadataStatement = dbConnection.prepareStatement("INSERT INTO Metadata(filename, file_size, last_write_time, first_1000_bytes, index_creation_time) VALUES(?,?,?,?,?)");

			metadataStatement.setString(1, metadata.getFileName());
			metadataStatement.setInt(2, metadata.getFileSize());
			metadataStatement.setInt(3, metadata.getLastWriteTime());
			metadataStatement.setBytes(4, metadata.getFirst1000bytes());
			metadataStatement.setInt(3, metadata.getIndexCreationTime());

			metadataStatement.executeUpdate();

			dbConnection.commit();

		} catch (SQLException ex) {
			throw new GenotypeDataException("Unable to save metadata in bgenix file. Error: " + ex.getMessage(), ex);
		}

	}

}
