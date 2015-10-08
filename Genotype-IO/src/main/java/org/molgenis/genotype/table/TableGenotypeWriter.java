/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.table;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.GenotypeData;
import org.molgenis.genotype.GenotypeDataException;
import org.molgenis.genotype.GenotypeWriter;
import org.molgenis.genotype.variant.GeneticVariant;

/**
 *
 * @author Patrick Deelen
 */
public class TableGenotypeWriter implements GenotypeWriter {

	private final GenotypeData genotypeData;

	public TableGenotypeWriter(GenotypeData genotypeData) {
		this.genotypeData = genotypeData;
	}

	@Override
	public void write(String path) {
		try {
			BufferedWriter dosageWriter = new BufferedWriter(new FileWriter(path + ".dosages.txt"));
			BufferedWriter genotypeWriter = new BufferedWriter(new FileWriter(path + ".genotypes.txt"));

			for (String sample : genotypeData.getSampleNames()) {
				dosageWriter.append('\t');
				dosageWriter.append(sample);
				genotypeWriter.append('\t');
				genotypeWriter.append(sample);
			}
			genotypeWriter.append('\n');
			dosageWriter.append('\n');

			for (GeneticVariant variant : genotypeData) {

				String variantId = variant.getPrimaryVariantId();
				if(variantId == null){
					variantId = variant.getSequenceName() + ":" + variant.getStartPos();
				}
								
				dosageWriter.append(variantId);
				genotypeWriter.append(variantId);

				for (float dosage : variant.getSampleDosages()) {
					dosageWriter.append('\t');
					dosageWriter.append(String.valueOf(dosage));
				}

				for (Alleles alleles : variant.getSampleVariants()) {

					genotypeWriter.append('\t');

					boolean notFirst = false;
					for (String allele : alleles.getAllelesAsString()) {
						if (notFirst) {
							genotypeWriter.append('/');
						}
						genotypeWriter.append(allele);
						notFirst = true;
					}

				}

				genotypeWriter.append('\n');
				dosageWriter.append('\n');

			}

			genotypeWriter.close();
			dosageWriter.close();
		} catch (IOException ex) {
			throw new GenotypeDataException("Error writing genotype tables: " + ex);
		}

	}
}
