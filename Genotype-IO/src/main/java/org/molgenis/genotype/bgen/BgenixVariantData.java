/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.bgen;

/**
 *
 * @author Patrick Deelen
 */
public class BgenixVariantData {

	private final String chromosome;
	private final int position;
	private final String rsid;
	private final int number_of_alleles;
	private final String allele1;
	private final String allele2;
	private final long file_start_position;
	private final int size_in_bytes;

	public BgenixVariantData(String chromosome, int position, String rsid, int number_of_alleles, String allele1, String allele2, long file_start_position, int size_in_bytes) {
		this.chromosome = chromosome;
		this.position = position;
		this.rsid = rsid;
		this.number_of_alleles = number_of_alleles;
		this.allele1 = allele1;
		this.allele2 = allele2;
		this.file_start_position = file_start_position;
		this.size_in_bytes = size_in_bytes;
	}

	public String getChromosome() {
		return chromosome;
	}

	public int getPosition() {
		return position;
	}

	public String getRsid() {
		return rsid;
	}

	public int getNumber_of_alleles() {
		return number_of_alleles;
	}

	public String getAllele1() {
		return allele1;
	}

	public String getAllele2() {
		return allele2;
	}

	public long getFile_start_position() {
		return file_start_position;
	}

	public int getSize_in_bytes() {
		return size_in_bytes;
	}

}
