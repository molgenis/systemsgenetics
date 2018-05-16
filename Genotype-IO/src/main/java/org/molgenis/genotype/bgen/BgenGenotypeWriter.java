/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.bgen;

import java.io.IOException;
import org.apache.log4j.Logger;
import org.molgenis.genotype.GenotypeData;
import org.molgenis.genotype.GenotypeWriter;
import org.molgenis.genotype.variant.NotASnpException;

/**
 *
 * @author patri
 */
public class BgenGenotypeWriter implements GenotypeWriter {

	private static final Logger LOGGER = Logger.getLogger(BgenGenotypeWriter.class);
	private final GenotypeData genotypeData;

	public BgenGenotypeWriter(GenotypeData genotypeData) {
		this.genotypeData = genotypeData;
	}

	@Override
	public void write(String path) throws IOException, NotASnpException {
		throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
	}

}
