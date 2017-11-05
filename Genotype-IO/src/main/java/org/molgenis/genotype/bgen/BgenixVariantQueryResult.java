/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.bgen;

import java.io.Closeable;
import java.io.IOException;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.Iterator;
import org.molgenis.genotype.GenotypeDataException;

/**
 *
 * @author Patrick Deelen
 */
public class BgenixVariantQueryResult implements Iterator<BgenixVariantData>, Closeable {

	private final ResultSet result;
	private boolean resultClosed;
	private boolean hasNext;
	private boolean didNext;

	public BgenixVariantQueryResult(ResultSet resultSet) {
		this.result = resultSet;
	}

	@Override
	public boolean hasNext() {
		
		if(resultClosed){
			return false;
		}
		
		try {
			if (!didNext) {

				hasNext = result.next();
				didNext = true;
			}

			if (!hasNext & !resultClosed) {
				result.close();
				resultClosed = true;
			}
		} catch (SQLException ex) {
			throw new GenotypeDataException("Unable to query bgenix file. Error: " + ex.getMessage(), ex);
		}

		return hasNext;
	}

	@Override
	public BgenixVariantData next() {
		
		if(resultClosed){
			throw new GenotypeDataException("Connection to bgenix already closed");
		}
		
		try {
			if (!didNext) {

				hasNext = result.next();

			}
			
			if(!hasNext){
				throw new GenotypeDataException("No more variants in query");
			}
			
			didNext = false;

			return new BgenixVariantData(
					result.getString("chromosome"),
					result.getInt("position"),
					result.getString("rsid"),
					result.getInt("number_of_alleles"),
					result.getString("allele1"),
					result.getString("allele2"),
					result.getLong("file_start_position"),
					result.getInt("size_in_bytes")
			);

		} catch (SQLException ex) {
			throw new GenotypeDataException("Unable to query bgenix file. Error: " + ex.getMessage(), ex);
		}

	}

	@Override
	public void close() throws IOException {
		try {
			if (!resultClosed) {
				result.close();
				resultClosed = true;
			}
		} catch (SQLException ex) {
			throw new GenotypeDataException("Unable to to close connection with bgenix file. Error: " + ex.getMessage(), ex);
		}
	}

}
