package org.molgenis.genotype.bgen;

import java.io.File;
import java.sql.Connection;
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
			throw new GenotypeDataException("Unable to load bgenix file. Error: " + ex.getMessage());
		}
		
	}
	
	public static Connection createNewConnection(File bgenixFile){
		Connection newDbConnection;
		try {
			newDbConnection = DriverManager.getConnection("jdbc:sqlite:"+bgenixFile.getPath());
		} catch (SQLException ex) {
			throw new GenotypeDataException("Unable to load bgenix file. Error: " + ex.getMessage());
		}
		return newDbConnection;
	}
	
	public synchronized ResultSet getVariantsChromosome(String chr){
		try {
			queryByChromosome.setString(1, chr);
			return queryByChromosome.executeQuery();
		} catch (SQLException ex) {
			throw new GenotypeDataException("Unable to load bgenix file. Error: " + ex.getMessage());
		}
	}
	
	public synchronized ResultSet getVariantsPostion(String chr, int position){
		try {
			queryByPosition.setString(1, chr);
			queryByPosition.setInt(1, position);
			return queryByPosition.executeQuery();
		} catch (SQLException ex) {
			throw new GenotypeDataException("Unable to load bgenix file. Error: " + ex.getMessage());
		}
	}

}
