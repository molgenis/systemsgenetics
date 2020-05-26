package org.molgenis.genotype;

import java.io.IOException;

import org.molgenis.genotype.variant.NotASnpException;

public interface GenotypeWriter
{

	/**
	 * 
	 * @param path
	 * @throws IOException
	 * @throws NotASnpException
	 *             when writer only supports SNPs but is given a non SNP variant
	 */
	public void write(String path) throws IOException, NotASnpException;

}
