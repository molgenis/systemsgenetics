package org.molgenis.genotype.util;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import org.apache.log4j.Logger;
import org.molgenis.genotype.GenotypeData;
import org.molgenis.genotype.GenotypeDataException;

/**
 * General Uilities
 *
 * @author erwin
 *
 */
public class Utils {

	public static final Map<Character, Character> COMPLEMENTAL_NUCLEOTIDES = new HashMap<Character, Character>();
	private static final Logger LOGGER = Logger.getLogger(GenotypeData.class);

	static {
		COMPLEMENTAL_NUCLEOTIDES.put('A', 'T');
		COMPLEMENTAL_NUCLEOTIDES.put('a', 't');
		COMPLEMENTAL_NUCLEOTIDES.put('T', 'A');
		COMPLEMENTAL_NUCLEOTIDES.put('t', 'a');
		COMPLEMENTAL_NUCLEOTIDES.put('C', 'G');
		COMPLEMENTAL_NUCLEOTIDES.put('c', 'g');
		COMPLEMENTAL_NUCLEOTIDES.put('G', 'C');
		COMPLEMENTAL_NUCLEOTIDES.put('g', 'c');
		COMPLEMENTAL_NUCLEOTIDES.put('0', '0');
	}

	/**
	 * Add all objects from the iterator in a ArrayList
	 *
	 * @param iterator
	 * @return
	 */
	public static <T> List<T> iteratorToList(Iterator<T> iterator) {
		List<T> result = new ArrayList<T>();
		while (iterator.hasNext()) {
			result.add(iterator.next());
		}

		return result;
	}

	/**
	 * Variant is a SNP if all alleles are one nucleotide long
	 */
	public static boolean isSnp(List<String> alleles) {
		for (String allele : alleles) {
			if ((allele != null) && allele.length() != 1) {
				return false;
			}
		}

		return true;
	}

	public static char getComplementNucleotide(char allele) {
		Character complementAllele = COMPLEMENTAL_NUCLEOTIDES.get(allele);
		if (complementAllele == null) {
			throw new GenotypeDataException("Failed to get comlpement for allele: " + allele);
		}
		return complementAllele;
	}

	public static void createEmptyFile(File file, String fileName) {

		if (file.exists()) {
			if (file.isDirectory()) {
				throw new GenotypeDataException("Cannot overwrite dir with " + fileName + " file:" + file.getAbsolutePath());
			}
			if (file.isFile()) {
				LOGGER.warn("Overriding " + fileName + " file " + file.getAbsolutePath());
                                
				if (!file.delete()) {
					throw new GenotypeDataException("Failed to overwrite " + fileName + " file: " + file.getAbsolutePath());
				}
			}
		}

		if (file.getParentFile() != null && !file.getParentFile().exists()) {
			if (!file.getParentFile().mkdirs()) {
				throw new GenotypeDataException("Failed to create parent dirs for " + fileName + " file: " + file.getAbsolutePath());
			}
		}

		try {
			if (!file.createNewFile()) {
				throw new GenotypeDataException("Error creating " + fileName + " file: " + file.getAbsolutePath());
			}
		} catch (IOException ex) {
			throw new GenotypeDataException("Error creating " + fileName + " file: " + file.getAbsolutePath(), ex);
		}

		if (!file.canWrite()) {
			throw new GenotypeDataException("Created " + fileName + " file but can not write to file: " + file.getAbsolutePath());
		}

	}
}
