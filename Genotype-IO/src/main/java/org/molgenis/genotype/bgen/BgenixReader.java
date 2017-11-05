package org.molgenis.genotype.bgen;

import java.io.File;
import java.io.IOException;
import java.sql.Connection;
import java.sql.DatabaseMetaData;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import org.molgenis.genotype.GenotypeDataException;

/**
 *
 * @author Patrick Deelen
 */
public class BgenixReader {

	private final Connection dbConnection;
	private final PreparedStatement queryByChromosome;
	private final PreparedStatement queryByPosition;
	private final PreparedStatement queryByRange;
	
	public BgenixReader(File bgenixFile) {
		this(createNewConnection(bgenixFile));
	}

	public BgenixReader(Connection dbConnection) {
		this.dbConnection = dbConnection;
		
		try {
			queryByChromosome = dbConnection.prepareStatement("SELECT * FROM Variant WHERE chromosome = ?");
			queryByPosition = dbConnection.prepareStatement("SELECT * FROM Variant WHERE chromosome = ? AND position = ?");
			queryByRange = dbConnection.prepareStatement("SELECT * FROM Variant WHERE chromosome = ? AND (position BETWEEN ? AND ?)");
		} catch (SQLException ex) {
			throw new GenotypeDataException("Unable to load bgenix file. Error: " + ex.getMessage(), ex);
		}
		
	}
	
	public static Connection createNewConnection(File bgenixFile){
		Connection newDbConnection;
		try {
			newDbConnection = DriverManager.getConnection("jdbc:sqlite:"+bgenixFile.getPath());
		} catch (SQLException ex) {
			throw new GenotypeDataException("Unable to load bgenix file. Error: " + ex.getMessage(), ex);
		}
		return newDbConnection;
	}
	
	public synchronized BgenixVariantQueryResult getVariantsChromosome(String chr){
		try {
			queryByChromosome.setString(1, chr);
			return new BgenixVariantQueryResult(queryByChromosome.executeQuery());
		} catch (SQLException ex) {
			throw new GenotypeDataException("Unable to query bgenix file. Error: " + ex.getMessage(), ex);
		}
	}
	
	public synchronized BgenixVariantQueryResult getVariantsPostion(String chr, int position){
		try {
			queryByPosition.setString(1, chr);
			queryByPosition.setInt(2, position);
			return new BgenixVariantQueryResult(queryByPosition.executeQuery());
		} catch (SQLException ex) {
			throw new GenotypeDataException("Unable to query bgenix file. Error: " + ex.getMessage(), ex);
		}
	}
	
	public synchronized BgenixVariantQueryResult getVariantsRange(String chr, int from, int to){
		try {
			queryByRange.setString(1, chr);
			queryByRange.setInt(2, from);
			queryByRange.setInt(3, to);
			return new BgenixVariantQueryResult(queryByRange.executeQuery());
		} catch (SQLException ex) {
			throw new GenotypeDataException("Unable to query bgenix file. Error: " + ex.getMessage(), ex);
		}
	}
	
	/**
	 * Returns null if no meta data is found in bgenix file
	 * 
	 * @return 
	 */
	public BgenixMetadata getMetadata() {
		
		try {
			DatabaseMetaData md = dbConnection.getMetaData();
			
			//Check if metadata present
			ResultSet rs = md.getTables(null, null, "Metadata", null);
			if(!rs.next()){
				return null;
			}
			
			ResultSet metaDataResultSet = dbConnection.createStatement().executeQuery("SELECT * FROM Metadata");
			
			if(!metaDataResultSet.next()){
				return null;
			}
			
			String fileName = metaDataResultSet.getString("filename");
			int fileSize = metaDataResultSet.getInt("file_size");
			int lastWriteTime = metaDataResultSet.getInt("last_write_time");
			byte[] first1000bytes = new byte[1000];
			metaDataResultSet.getBinaryStream("first_1000_bytes").read(first1000bytes);
			int indexCreationTime = metaDataResultSet.getInt("index_creation_time");
			
			return new BgenixMetadata(fileName, fileSize, lastWriteTime, first1000bytes, indexCreationTime);
			
		} catch (SQLException | IOException ex) {
			throw new GenotypeDataException("Unable to read bgenix metadata. Error: " + ex.getMessage(), ex);
		}
		
	}

}
