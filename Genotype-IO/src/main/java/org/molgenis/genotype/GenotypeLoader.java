/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype;

import java.io.File;
import java.util.Arrays;
import java.util.EnumMap;
import java.util.EnumSet;
import java.util.List;
import org.apache.log4j.Logger;
import static org.molgenis.genotype.GenotypeFileType.*;

/**
 *
 * @author Patrick Deelen
 */
public class GenotypeLoader {

	private static final Logger LOGGER = Logger.getLogger(GenotypeLoader.class);

	public static RandomAccessGenotypeData load(String... paths) {
		return load(Arrays.asList(paths));
	}

	private static EnumMap<GenotypeFileType, File> deterimeFileList(List<String> paths){
		
		EnumMap<GenotypeFileType, File> files = new EnumMap<GenotypeFileType, File>(GenotypeFileType.class);
		EnumSet<GenotypeFileType> overridableFiles = EnumSet.copyOf(GenotypeFileType.getTrityperFileTypes());
		
		if (paths.isEmpty()) {
			LOGGER.debug("No paths to genotype data specified");
			throw new GenotypeDataException("No path to genotype data specified");
		} else if (paths.size() == 1) {
			//potenially were dealing with a prefix and not individual files
			LOGGER.debug("One path for genotype data specified, checking if it is a prefix");

			String path = paths.get(0);
			File file = new File(path);

			if (!file.isDirectory()) {

				//We deal with dirs later

				for (GenotypeFileType type : GenotypeFileType.values()) {

					File prefixedFile = new File(path + type.getSuffix());

					if (prefixedFile.exists()) {
						LOGGER.debug("Detected " + type.getFriendlyName() + " file at: " + prefixedFile.getAbsolutePath());
						files.put(type, prefixedFile);
					}


				}


			}

		}

		for (String path : paths) {


			File file = new File(path);

			if (file.isDirectory()) {
				LOGGER.debug("Genotype data path is directory");

				folderFiles:
				for (File folderFile : file.listFiles()) {

					if (folderFile.getName().endsWith(".vcf.gz")) {
						LOGGER.info("Genotype loader detected vcf.gz files in specified folder. Assuming " + RandomAccessGenotypeDataReaderFormats.VCF_FOLDER.getName() + ", " + RandomAccessGenotypeDataReaderFormats.VCF_FOLDER.getDescription());
						if (files.containsKey(VCF_FOLDER)) {
							throw new GenotypeDataException("Detected second VCF folder. That is currently not supported");
						}
						files.put(VCF_FOLDER, file);
						continue folderFiles;
					}

					for (GenotypeFileType trityerType : GenotypeFileType.getTrityperFileTypes()) {

						if(files.containsKey(trityerType)){
								//ignore if already set to allow overwrite of defaults
								continue;
							}
						
						if (folderFile.getName().equals(trityerType.getSuffix())) {	
							files.put(trityerType, folderFile);
							continue folderFiles;
						}

						if (trityerType != TRITYPER_GENOTYPE && trityerType != TRITYPER_DOSAGE) {
							if (folderFile.getName().equals(trityerType.getSuffix() + ".gz")) {
								files.put(trityerType, folderFile);
								continue folderFiles;
							}

						}


					}


				}


			} else if (file.exists()) {
				
				for (GenotypeFileType type : GenotypeFileType.values()) {
					
					if(type.matches(file)){
						
						LOGGER.debug("Found " + type.getFriendlyName() + " file at: " + file.getAbsolutePath());
						
						if(files.containsKey(type) && !overridableFiles.contains(type)){
							
							throw new GenotypeDataException("Found " + type.getFriendlyName() + " file twice. First at: " + files.get(type).getAbsolutePath() + " and then again at: " + file.getAbsolutePath());
							
						} else {
							
							files.put(type, file);
							overridableFiles.remove(type);//only overwrite once.
							
						}
					
					}
					
				}
				
			} else {
				//file does not exist. That is fine if we only specified one path as a prefix.

				if (paths.size() == 1) {
					continue;
				} else {
					throw new GenotypeDataException("No genotype file found at: " + file.getAbsolutePath());
				}

			}


		}
		
		if(files.containsKey(SAMPLE) && !files.containsKey(HAPS) && !files.containsKey(GEN) && files.containsKey(UNKNOWN)){
			LOGGER.info("Found a file without an extention and and an " + SAMPLE.getFriendlyName() + " file. Assuming that this file: " + files.get(UNKNOWN) + " is a " + GEN.getFriendlyName() + " file");
			files.put(GEN, files.remove(UNKNOWN));
		}
		
		return files;
		
	}
	
	public static RandomAccessGenotypeData load(List<String> paths) {
		
		EnumMap<GenotypeFileType, File> files = deterimeFileList(paths);
		//By this point we should have an invatory of all matched files
	
		RandomAccessGenotypeDataReaderFormats identifiedFormat = null;
		
		for(RandomAccessGenotypeDataReaderFormats format : RandomAccessGenotypeDataReaderFormats.values()){
			if(files.keySet().containsAll(format.getRequiredFiles())){
				
				
				
				if(identifiedFormat != null){
					throw new GenotypeDataException("Unable to unambiguously identify the genotype file format. Could be: " + identifiedFormat.getName() + " or " + format.getName() + " based on the identified files");
				} else if(false) {
					//TODO
					identifiedFormat = format;
				}
			}
		}

		if(identifiedFormat != null){
			
		}
		
		return null;

	}

	public static RandomAccessGenotypeData loadFormat(String format, String... paths) {
		return null;
	}

	public static RandomAccessGenotypeData loadFormat(RandomAccessGenotypeDataReaderFormats format, String... paths) {
		return null;
	}
}
