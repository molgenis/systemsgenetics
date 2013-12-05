package umcg.genetica.io.trityper.util;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.bin.BinaryFile;
import umcg.genetica.io.text.TextFile;

/**
 * This class can be used to concatenate several trityper datasets. All datasets
 * must contain the exact same number of individuals in the same order. If first
 * dataset to concat contains dosage it is assumed that other datasets also have
 * dosage. Individuals.txt of first dataset to concat is used.
 *
 * @author Patrick Deelen
 */
public class TriTyperConcatDatasets {

    private final ArrayList<String> datasetsToConcatFolders;
    private final String outputFolder;
    private final HashSet<String> totalSnps;

    public TriTyperConcatDatasets(String outputFolderPath) throws IOException {
        this.outputFolder = outputFolderPath;
        this.datasetsToConcatFolders = new ArrayList<String>();
        this.totalSnps = new HashSet<String>();

        Gpio.createDir(outputFolderPath);
    }

    public void addDatasetToConcat(String datasetFolderPath) throws IOException {

        if (!Gpio.exists(datasetFolderPath)) {
            throw new FileNotFoundException("Folder with dataset to concat does not exists: " + datasetFolderPath);
        }

        if (!Gpio.isDir(datasetFolderPath)) {
            throw new IOException("Path to folder with dataset to concat is not a folder: " + datasetFolderPath);
        }

        if (!Gpio.canRead(datasetFolderPath)) {
            throw new IOException("Can not read folder of dataset to concat: " + datasetFolderPath);
        }

        if (!datasetFolderPath.endsWith(Gpio.getFileSeparator())) {
            datasetFolderPath += Gpio.getFileSeparator();
        }
        this.datasetsToConcatFolders.add(datasetFolderPath);

    }

    public void writeConcatedDataset() throws FileNotFoundException, IOException, Exception {

        if (this.datasetsToConcatFolders.isEmpty()) {
            System.out.println("No data to concatenate");
            return;
        }

        System.out.println("Concatenating data to: " + this.outputFolder);

        TextFile outputSnps = new TextFile(this.outputFolder + "SNPs.txt", TextFile.W);
        BinaryFile outputGenotypeMatrix = new BinaryFile(this.outputFolder + "GenotypeMatrix.dat", BinaryFile.W, (8 * 1024));
        BinaryFile outputImputedDosageMatrix = null;
		TextFile outputMappings = null;
		
        final boolean concatImputedDosageMatrix = hasDatasetFolderImputedDosageMatrix(this.datasetsToConcatFolders.get(0));
		final boolean concatSnpMappings = hasDatasetFolderSnpMapping(this.datasetsToConcatFolders.get(0));

        //This will be filled be first dataset.
        List<String> individuals = null;

        if (concatImputedDosageMatrix) {
            outputImputedDosageMatrix = new BinaryFile(this.outputFolder + "ImputedDosageMatrix.dat", BinaryFile.W, (8 * 1024));
        }
		
		if(concatSnpMappings) {
			System.out.println("SNP mappigns detected and will be concatted");
			outputMappings = new TextFile(this.outputFolder + "SNPMappings.txt", TextFile.W);
		}

        for (String datasetToConcatFolder : this.datasetsToConcatFolders) {

            System.out.println("Adding data from: " + datasetToConcatFolder);

            if (individuals == null) {
                individuals = readIndividuals(datasetToConcatFolder);
            } else {
                if (!individuals.equals(readIndividuals(datasetToConcatFolder))) {
                    throw new Exception("Individuals.txt does not match for dataset (check whether the order, number and identifiers of samples are identical): " + datasetToConcatFolder);
                }
            }

            TextFile intputSnps = new TextFile(datasetToConcatFolder + "SNPs.txt", TextFile.R);
            addSnpsToTotalList(intputSnps, outputSnps);
            intputSnps.close();

            BinaryFile intputGenotypeMatrix = new BinaryFile(datasetToConcatFolder + "GenotypeMatrix.dat", BinaryFile.R, (8 * 1024));
            addInputStreamToOutputStream(intputGenotypeMatrix, outputGenotypeMatrix);
            intputGenotypeMatrix.close();

            if (concatImputedDosageMatrix) {
                BinaryFile intputImputedDosageMatrix = new BinaryFile(datasetToConcatFolder + "ImputedDosageMatrix.dat", BinaryFile.R, (8 * 1024));
                addInputStreamToOutputStream(intputImputedDosageMatrix, outputImputedDosageMatrix);
                intputImputedDosageMatrix.close();
            }
			
			if(concatSnpMappings) {
				TextFile inputSnpMappings = new TextFile(datasetToConcatFolder + "SNPMappings.txt", TextFile.R);
				copyContent(inputSnpMappings, outputMappings);
				inputSnpMappings.close();
			}

        }	

        outputSnps.close();
        outputGenotypeMatrix.close();
        if (concatImputedDosageMatrix) {
            outputImputedDosageMatrix.close();
        }
		
		if(concatSnpMappings) {
			outputMappings.close();
		}

        writeIndividuals(this.outputFolder, individuals);

		if(hasDatasetFolderPhenoInfo(this.datasetsToConcatFolders.get(0))){
			System.out.println("Detected pheno information this is copied");
			
			TextFile outputPheno = new TextFile(this.outputFolder + "PhenotypeInformation.txt", TextFile.W);
			TextFile inputPheno = new TextFile(this.datasetsToConcatFolders.get(0) + "PhenotypeInformation.txt", TextFile.R);
			
			copyContent(inputPheno, outputPheno);
			outputPheno.close();
			inputPheno.close();
		}
		
        System.out.println("Concatenated data written");

    }

    private void addInputStreamToOutputStream(BinaryFile input, BinaryFile output) throws IOException {
        int i;
        while ((i = input.read()) != -1) {
            output.write(i);
        }
    }

    private void addSnpsToTotalList(TextFile snpFileReader, TextFile snpFileWriter) throws IOException, Exception {

        String line;
        while ((line = snpFileReader.readLine()) != null) {

            if (this.totalSnps.contains(line)) {
               // throw new Exception("Snp: " + line + " is present twice");
            }
            totalSnps.add(line);
            snpFileWriter.writeln(line);
        }
    }
	
		private void copyContent(TextFile source, TextFile target) throws IOException, Exception {

        String line;
        while ((line = source.readLine()) != null) {
            target.writeln(line);
        }
		
    }

    private boolean hasDatasetFolderImputedDosageMatrix(String trityperDatasetFolder) {
        return Gpio.exists(trityperDatasetFolder + "ImputedDosageMatrix.dat");
    }
	
	private boolean hasDatasetFolderSnpMapping(String trityperDatasetFolder) {
        return Gpio.exists(trityperDatasetFolder + "SNPMappings.txt");
    }
	
	private boolean hasDatasetFolderPhenoInfo(String trityperDatasetFolder) {
        return Gpio.exists(trityperDatasetFolder + "PhenotypeInformation.txt");
    }


    private List<String> readIndividuals(String trityperDatasetFolder) throws FileNotFoundException, IOException {

        ArrayList<String> individuals = null;

        TextFile sampleFile = new TextFile(trityperDatasetFolder + "Individuals.txt", TextFile.R);
        individuals = sampleFile.readAsArrayList();
        sampleFile.close();
        return Collections.unmodifiableList(individuals);

    }

    private void writeIndividuals(String targetFolder, List<String> individuals) throws IOException {

        TextFile individualsWriter = new TextFile(targetFolder + "Individuals.txt", TextFile.W);
        individualsWriter.writeList(individuals);
        individualsWriter.close();

    }
}
