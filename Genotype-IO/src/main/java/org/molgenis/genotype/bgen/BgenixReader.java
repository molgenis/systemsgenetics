package org.molgenis.genotype.bgen;

import java.io.File;
import java.sql.Connection;
import java.sql.DriverManager;
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
	
	public BgenixReader(File bgenixFile) {
		this(createNewConnection(bgenixFile));
	}

	public BgenixReader(Connection dbConnection) {
		this.dbConnection = dbConnection;
		
		
	}
	
	public static Connection createNewConnection(File bgenixFile){
		Connection newDbConnection;
		try {
			newDbConnection = DriverManager.getConnection(bgenixFile.getPath());
		} catch (SQLException e) {
			throw new GenotypeDataException("Unable to load bgenix file. Error: " + e.getMessage());
		}
		return newDbConnection;
	}
	
	public ResultSet getVariantsChromosome(String chr){
		return null;
	}

}
