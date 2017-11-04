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
import java.util.logging.Level;
import java.util.logging.Logger;
import org.molgenis.genotype.GenotypeDataException;
import org.molgenis.genotype.variant.GeneticVariant;

/**
 *
 * @author patri
 */
public class BgenixWriter {
	
	private final Connection dbConnection;
	private final PreparedStatement addVariantStatement;

	public BgenixWriter(File bgenixFile) {
		
		if(bgenixFile.exists()){
			throw new GenotypeDataException("Bgen index file already exisist: " + bgenixFile.getAbsolutePath());
		}
		
		try {
			dbConnection = DriverManager.getConnection("jdbc:sqlite: " + bgenixFile.getPath());
		} catch (SQLException ex) {
			throw new GenotypeDataException("Unable to create bgenix file. Error: " + ex.getMessage());
		}
		
		try {
			
			dbConnection.createStatement().execute("CREATE TABLE Variant ( chromosome TEXT NOT NULL, position INT NOT NULL, rsid TEXT NOT NULL, number_of_alleles INT NOT NULL, allele1 TEXT NOT NULL, allele2 TEXT NULL, file_start_position INT NOT NULL, size_in_bytes INT NOT NULL, PRIMARY KEY (chromosome, position, rsid, allele1, allele2, file_start_position )) WITHOUT ROWID;");
			
			addVariantStatement = dbConnection.prepareStatement("INSERT INTO Metadata(chromosome, position, rsid, number_of_alleles, allele1, allele2, file_start_position, size_in_bytes) VALUES(?,?,?,?,?,?,?,?)");
			
		} catch (SQLException ex) {
			throw new GenotypeDataException("Unable to initialize bgenix file. Error: " + ex.getMessage());
		}
		
	}
	
	public void addVariantToIndex(GeneticVariant variant, long startPos, int sizeInBytesInBgenFile){
		
	}
	
	public void writeMetadata(BgenixMetadata metadata){
		
		try {
			
			dbConnection.createStatement().execute("CREATE TABLE Metadata (filename TEXT NOT NULL, file_size INT NOT NULL, last_write_time INT NOT NULL, first_1000_bytes BLOB NOT NULL, index_creation_time INT NOT NULL);");
			
			PreparedStatement metadataStatement = dbConnection.prepareStatement("INSERT INTO Metadata(filename, file_size, last_write_time, first_1000_bytes, index_creation_time) VALUES(?,?,?,?,?)");
			
			metadataStatement.setString(1, metadata.getFileName());
			metadataStatement.setInt(2, metadata.getFileSize());
			metadataStatement.setInt(3, metadata.getLastWriteTime());
			metadataStatement.setBytes(4, metadata.getFirst1000bytes());
			metadataStatement.setInt(3, metadata.getIndexCreationTime());
			
			metadataStatement.executeUpdate();
			
		} catch (SQLException ex) {
			throw new GenotypeDataException("Unable to save metadata in bgenix file. Error: " + ex.getMessage());
		}
		
		
		
	}
	
	
}
