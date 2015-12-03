package org.molgenis.genotype.trityper;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashSet;
import org.apache.log4j.Logger;
import org.molgenis.genotype.Allele;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.GenotypeData;
import org.molgenis.genotype.GenotypeDataException;
import org.molgenis.genotype.GenotypeWriter;
import org.molgenis.genotype.Sample;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variant.NotASnpException;
import org.molgenis.genotype.variant.id.GeneticVariantId;

/**
 *
 * @author Patrick Deelen
 */
public class TriTyperGenotypeWriter implements GenotypeWriter {

	private static Logger LOGGER = Logger.getLogger(TriTyperGenotypeWriter.class);
	private final GenotypeData genotypeData;

	public TriTyperGenotypeWriter(GenotypeData genotypeData) {
		this.genotypeData = genotypeData;
	}

	@Override
	public void write(String path) throws IOException, NotASnpException {
		write(new File(path));
	}

	public void write(File folder) throws IOException {

		if (!folder.isDirectory()) {
			if (!folder.mkdirs()) {
				throw new GenotypeDataException("Failed to create trityper dir at: " + folder.getAbsolutePath());
			}
		}

		File genotypeDataFile = new File(folder, "GenotypeMatrix.dat");
		File imputedDosageDataFile = new File(folder, "ImputedDosageMatrix.dat");
		File snpFile = new File(folder, "SNPs.txt");
		File snpMapFile = new File(folder, "SNPMappings.txt");
		File individualFile = new File(folder, "Individuals.txt");
		File phenotypeAnnotationFile = new File(folder, "PhenotypeInformation.txt");
        File allelRecodingFile = new File(folder, "AlleleRecodingInformation.txt");

		writeSnps(snpFile, snpMapFile);
		writeSamples(individualFile, phenotypeAnnotationFile);
		writeGenotypes(genotypeDataFile, imputedDosageDataFile, allelRecodingFile);


	}

	private void writeSnps(File snpFile, File snpMapFile) throws IOException {

		BufferedWriter snpFileWriter = new BufferedWriter(new FileWriter(snpFile));
		BufferedWriter snpMapFileWriter = new BufferedWriter(new FileWriter(snpMapFile));

		for (GeneticVariant variant : genotypeData) {

//			if (!variant.isSnp()) {
//				LOGGER.warn("Skipping variant: " + variant.getPrimaryVariantId() + ", it is not a SNP.\n");
//                LOGGER.warn("Will create a remapping file to use the data in TriTyper file format.");
//				continue;
//			}

			final String snpName = createTriTyperVariantId(variant);
			
			snpFileWriter.append(snpName);
			snpFileWriter.append('\n');

			snpMapFileWriter.append(variant.getSequenceName());
			snpMapFileWriter.append('\t');
			snpMapFileWriter.append(String.valueOf(variant.getStartPos()));
			snpMapFileWriter.append('\t');
			snpMapFileWriter.append(snpName);
			snpMapFileWriter.append('\n');
		}

		snpFileWriter.close();
		snpMapFileWriter.close();

	}

	private void writeSamples(File individualFile, File phenotypeAnnotationFile) throws IOException {

		BufferedWriter individualFileWriter = new BufferedWriter(new FileWriter(individualFile));
		BufferedWriter phenotypeAnnotationFileWriter = new BufferedWriter(new FileWriter(phenotypeAnnotationFile));

		for (Sample sample : genotypeData.getSamples()) {

			individualFileWriter.append(sample.getId());
			individualFileWriter.append('\n');

			phenotypeAnnotationFileWriter.append(sample.getId());
			phenotypeAnnotationFileWriter.append('\t');
			phenotypeAnnotationFileWriter.append(sample.getCaseControlAnnotation().getTriTyperName());
			phenotypeAnnotationFileWriter.append('\t');
			phenotypeAnnotationFileWriter.append(sample.isIncluded() ? "include" : "exclude");
			phenotypeAnnotationFileWriter.append('\t');
			phenotypeAnnotationFileWriter.append(sample.getSex().getGender().toLowerCase());
			phenotypeAnnotationFileWriter.append('\n');

		}

		individualFileWriter.close();
		phenotypeAnnotationFileWriter.close();

	}

	private void writeGenotypes(File genotypeDataFile, File imputedDosageDataFile, File allelRecodingFile) throws IOException {

		// no need for buffered stream writer. data we write per SNP.
		FileOutputStream genotypeDataFileWriter = new FileOutputStream(genotypeDataFile);
		FileOutputStream genotypeDosageDataFileWriter = new FileOutputStream(imputedDosageDataFile);
        
        HashSet<String> snpRecodingInfo = new HashSet<String>();
        
        //Should we skip writing the genotypes?
        
		String[] samples = genotypeData.getSampleNames();
		int sampleCount = samples.length;

		byte[] snpBuffer = new byte[sampleCount * 2];
		byte[] dosageBuffer = new byte[sampleCount];

		for (GeneticVariant variant : genotypeData) {

			float[] dosageValues = variant.getSampleDosages();
			int i = 0;
			for (Alleles sampleAlleles : variant.getSampleVariants()) {
                
				if (sampleAlleles.getAlleleCount() != 2) {
					LOGGER.debug("variant at: " + variant.getSequenceName() + ":" + variant.getStartPos() + " set to missing for " + samples[i]);
					sampleAlleles = Alleles.BI_ALLELIC_MISSING;
                    //ToDo should we continue here?
				}

				try {
                    byte a;
                    byte b;
                    if (variant.isSnp()){
                        a = sampleAlleles.get(0).isSnpAllele() && sampleAlleles.get(0) != Allele.ZERO ? (byte) sampleAlleles.get(0).getAlleleAsSnp() : 0;
                        b = sampleAlleles.get(1).isSnpAllele() && sampleAlleles.get(1) != Allele.ZERO ? (byte) sampleAlleles.get(1).getAlleleAsSnp() : 0;
                    } else {
                        snpRecodingInfo.add(createTriTyperVariantId(variant)+"\t"+variant.getSequenceName()+"\t"+variant.getStartPos()+"\t"+variant.getVariantAlleles().get(0)+"\t"+variant.getVariantAlleles().get(1));
                        
                        if(sampleAlleles.get(0).equals(variant.getVariantAlleles().get(0))){
                            a = (byte) 'A';
                        } else if(sampleAlleles.get(0).equals(variant.getVariantAlleles().get(1))){
                            a = (byte) 'C';
                        } else {
                            a = 0;
                        }
                        
                        if(sampleAlleles.get(1).equals(variant.getVariantAlleles().get(0))){
                            b = (byte) 'A';
                        } else if(sampleAlleles.get(1).equals(variant.getVariantAlleles().get(1))){
                            b = (byte) 'C';
                        } else {
                            b = 0;
                        }

                    }
                    
                    snpBuffer[i] = a;
                    snpBuffer[i + sampleCount] = b;

				} catch (Exception e) {
					throw new GenotypeDataException("Error writing TriTyper data: " + e.getMessage(), e);
				}

				float dosage = dosageValues[i];
				if (dosage == -1) {
					dosageBuffer[i] = (byte) 127;
				} else {
					int dosageInt = (int) Math.round(dosage * 100d);
					dosageBuffer[i] = (byte) (Byte.MIN_VALUE + dosageInt);
				}

				++i;
			}

			genotypeDataFileWriter.write(snpBuffer);
			genotypeDosageDataFileWriter.write(dosageBuffer);

		}
        
		genotypeDataFileWriter.close();
		genotypeDosageDataFileWriter.close();
        
        if(!snpRecodingInfo.isEmpty()){
            BufferedWriter allelRecodingFileWriter = new BufferedWriter(new FileWriter(allelRecodingFile));
            
            allelRecodingFileWriter.write("Variant_ID\tChr\tPos\tAllele1\tAllele2\n");
            for(String s : snpRecodingInfo){
                allelRecodingFileWriter.write(s+"\n");
            }
            
            allelRecodingFileWriter.close();
        }
	}

	private String createTriTyperVariantId(GeneticVariant variant) {
		final GeneticVariantId snpId = variant.getVariantId();
		String snpName;
		if(snpId.containsId()){
			snpName = snpId.getPrimairyId();
		} else {
			snpName = variant.getSequenceName() + ':' + String.valueOf(variant.getStartPos());
			if(!variant.isSnp()){
				for(Allele allele : variant.getVariantAlleles()){
					snpName = snpName + "_" + allele.getAlleleAsString();
				}
				
			}
		}
		return snpName;
	}
}
