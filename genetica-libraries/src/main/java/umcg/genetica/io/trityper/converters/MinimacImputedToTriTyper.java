/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.io.trityper.converters;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.WGAFileMatrixGenotype;
import umcg.genetica.io.trityper.WGAFileMatrixImputedDosage;

/**
 *
 * @author lude
 */
public class MinimacImputedToTriTyper {

	private static class MinimacSnp {

		private final char a1;
		private final char a2;
		private final String snpId;

		public MinimacSnp(char a1, char a2, String snpId) {
			this.a1 = a1;
			this.a2 = a2;
			this.snpId = snpId;
		}

		public char getA2() {
			return a2;
		}

		public char getA1() {
			return a1;
		}

		public String getSnpId() {
			return snpId;
		}

		public char[] getAllelesForDossage(float dosage) {

			char[] alleles = new char[2];
			
			switch (Math.round(dosage)) {
				case 0:
					alleles[0] = a2;
					alleles[1] = a2;
					break;
				case 1:
					alleles[0] = a1;
					alleles[1] = a2;
					break;
				case 2:
					alleles[0] = a1;
					alleles[1] = a1;
					break;

			}
			//System.out.println(dosage + " -> " + alleles[0] + "/" + alleles[1]);
			return alleles;

		}
	}

	public static void convertMinimacImputedToTriTyper(String dir, String outputDir) throws IOException {
		Gpio.createDir(outputDir);
		HashMap<String, Integer> hashSNPs = new HashMap<String, Integer>();
		ArrayList<String> vecSNPs = new ArrayList<String>();
		HashMap<Integer, MinimacSnp> minimacSnpInfo = new HashMap<Integer, MinimacSnp>();


		//Load all imputed SNPs:
		System.out.println("Determining number of unique SNPs:");
		TextFile outSNPs = new TextFile(outputDir + "/SNPs.txt", TextFile.W);

		String[] fileList = Gpio.getListOfFiles(dir, "info");
		TextFile in = null;
		String[] str = null;
		for (String filename : fileList) {
			System.out.println("Processing file:\t" + filename);

			in = new TextFile(filename, TextFile.R);
			str = in.readLineElemsReturnReference(TextFile.tab);

			while (str != null) {
				String snp = new String(str[0].getBytes());
				
				if (!hashSNPs.containsKey(snp)) {
					hashSNPs.put(snp, vecSNPs.size());
					minimacSnpInfo.put(vecSNPs.size(), new MinimacSnp(str[1].charAt(0), str[2].charAt(0), snp));
					vecSNPs.add(snp);
				}
				
				
				
				outSNPs.write(snp + "\n");
				str = in.readLineElemsReturnReference(TextFile.tab);
			}
			in.close();
			in = null;
		}
		outSNPs.close();

		//Load all individuals:
		System.out.println("Determining unique number of individuals:");
		HashMap<String, Integer> hashInds = new HashMap<String, Integer>();
		ArrayList<String> vecInds = new ArrayList<String>();
		fileList = Gpio.getListOfFiles(dir, "dose");
		for (String filename : fileList) {
			in = new TextFile(filename, TextFile.R);
			str = in.readLineElemsReturnReference(TextFile.tab);
			while (str != null) {
				if (!hashInds.containsKey(str[0])) {
					hashInds.put(str[0], vecInds.size());
					vecInds.add(str[0]);
				}
				str = in.readLineElemsReturnReference(TextFile.tab);
			}
		}


		TextFile outInds = new TextFile(outputDir + "/Individuals.txt", TextFile.W);
		for (String ind : vecInds) {
			String individual = ind.split(">")[1];
			outInds.write(individual + "\n");
		}
		outInds.close();

		System.out.println("Number of unique SNPs:\t" + vecSNPs.size());
		System.out.println("Number of unique individuals:\t" + vecInds.size());

		//Process genotypes:
		System.out.println("Importing genotypes:");
		File fileGenotypeMatrix = new File(outputDir + "/GenotypeMatrix.dat");
		WGAFileMatrixGenotype matrixGenotype = new WGAFileMatrixGenotype(vecSNPs.size(), vecInds.size(), fileGenotypeMatrix, false);

		File fileImputedDosageMatrix = new File(outputDir + "/ImputedDosageMatrix.dat");
		WGAFileMatrixImputedDosage matrixImputedDosage = new WGAFileMatrixImputedDosage(vecSNPs.size(), vecInds.size(), fileImputedDosageMatrix, false);

		fileList = Gpio.getListOfFiles(dir, "info");

		for (String filename : fileList) {
			System.out.println("Processing " + filename);
			int[] snpIDArray = new int[vecSNPs.size()];
			in = new TextFile(filename, TextFile.R);
			str = in.readLineElemsReturnReference(TextFile.tab);
			int counter = 0;
			while (str != null) {
				String snp = new String(str[0].getBytes());
				snpIDArray[counter] = hashSNPs.get(snp);
				counter++;
				str = in.readLineElemsReturnReference(TextFile.tab);
			}
			in.close();

			String genofilename = filename.substring(0, filename.length() - 5);
			genofilename += ".dose";
			in = new TextFile(genofilename, TextFile.R);
	
			while ( (str = in.readLineElemsReturnReference(TextFile.tab))!= null){

				int indID = hashInds.get(str[0]);
				for (int c = 2; c < str.length; c++) {
					int snpID = snpIDArray[c - 2];
					float dosageFloat = Float.parseFloat(str[c]);
					int dosageInt = Math.round(dosageFloat * 100f);
					if (dosageInt < 0 || dosageInt > 200) {
						System.out.println("Warning, incorrect dosage!:\t" + dosageInt + "\t" + snpID + "\t" + str[c]);
					}
					byte dosageByte = (byte) (Byte.MIN_VALUE + dosageInt);
					byte[] dosage = new byte[1];
					dosage[0] = (byte) dosageByte;
					matrixImputedDosage.setDosage(snpID, indID, dosage);

					char[] alleles = minimacSnpInfo.get(snpID).getAllelesForDossage(dosageFloat);

					byte[] allele1 = new byte[1];
					allele1[0] = (byte) alleles[0];
					byte[] allele2 = new byte[1];
					allele2[0] = (byte) alleles[1];
					
					

					matrixGenotype.setAllele1(snpID, indID, allele1);
					matrixGenotype.setAllele2(snpID, indID, allele2);


				}
				System.out.println("Sample\t" + str[0] + "\thas been processed.");
			}
			in.close();
		}

		matrixGenotype.close();
		matrixImputedDosage.close();


		System.out.println("\n\n\n");
		System.out.println("Minimac imputed data has been imported.");
		System.out.println("\n\n\n");

	}
	
	
	public static void convertSingleFileMinimacImputedToTriTyper(String baseFilePath, String outputDir) throws IOException {
		//quick hack to do this one file based on: convertMinimacImputedToTriTyper
		
		Gpio.createDir(outputDir);
		HashMap<String, Integer> hashSNPs = new HashMap<String, Integer>();
		ArrayList<String> vecSNPs = new ArrayList<String>();
		HashMap<Integer, MinimacSnp> minimacSnpInfo = new HashMap<Integer, MinimacSnp>();


		//Load all imputed SNPs:
		System.out.println("Determining number of unique SNPs:");
		TextFile outSNPs = new TextFile(outputDir + "/SNPs.txt", TextFile.W);

		
		TextFile in = null;
		String[] str = null;
		
		String filename = baseFilePath + ".info";
			System.out.println("Processing file:\t" + filename);

			in = new TextFile(filename, TextFile.R);
			//skip first line in the file
			in.readLineElemsReturnReference(TextFile.tab);
			
			
			str = in.readLineElemsReturnReference(TextFile.tab);

			while (str != null) {
				String snp = new String(str[0].getBytes());
				
				if (!hashSNPs.containsKey(snp)) {
					hashSNPs.put(snp, vecSNPs.size());
					minimacSnpInfo.put(vecSNPs.size(), new MinimacSnp(str[1].charAt(0), str[2].charAt(0), snp));
					vecSNPs.add(snp);
				}
				
				outSNPs.write(snp + "\n");
				str = in.readLineElemsReturnReference(TextFile.tab);
			}
			in.close();
			in = null;
		
		outSNPs.close();

		//Load all individuals:
		System.out.println("Determining unique number of individuals:");
		HashMap<String, Integer> hashInds = new HashMap<String, Integer>();
		ArrayList<String> vecInds = new ArrayList<String>();
	
		
			
			filename = baseFilePath + ".dose";
			in = new TextFile(filename, TextFile.R);
			str = in.readLineElemsReturnReference(TextFile.tab);
			while (str != null) {
				if (!hashInds.containsKey(str[0])) {
					hashInds.put(str[0], vecInds.size());
					vecInds.add(str[0]);
				}
				str = in.readLineElemsReturnReference(TextFile.tab);
			}
		


		TextFile outInds = new TextFile(outputDir + "/Individuals.txt", TextFile.W);
		for (String ind : vecInds) {
			String individual = ind.split(">")[1];
			outInds.write(individual + "\n");
		}
		outInds.close();

		System.out.println("Number of unique SNPs:\t" + vecSNPs.size());
		System.out.println("Number of unique individuals:\t" + vecInds.size());

		//Process genotypes:
		System.out.println("Importing genotypes:");
		File fileGenotypeMatrix = new File(outputDir + "/GenotypeMatrix.dat");
		WGAFileMatrixGenotype matrixGenotype = new WGAFileMatrixGenotype(vecSNPs.size(), vecInds.size(), fileGenotypeMatrix, false);

		File fileImputedDosageMatrix = new File(outputDir + "/ImputedDosageMatrix.dat");
		WGAFileMatrixImputedDosage matrixImputedDosage = new WGAFileMatrixImputedDosage(vecSNPs.size(), vecInds.size(), fileImputedDosageMatrix, false);

		
			filename = baseFilePath + ".info";
			System.out.println("Processing " + filename);
			int[] snpIDArray = new int[vecSNPs.size()];
			in = new TextFile(filename, TextFile.R);
			
			//skip first line
			in.readLineElemsReturnReference(TextFile.tab);
			
			str = in.readLineElemsReturnReference(TextFile.tab);
			int counter = 0;
			while (str != null) {
				String snp = new String(str[0].getBytes());
				snpIDArray[counter] = hashSNPs.get(snp);
				counter++;
				str = in.readLineElemsReturnReference(TextFile.tab);
			}
			in.close();

			String genofilename = filename.substring(0, filename.length() - 5);
			genofilename += ".dose";
			in = new TextFile(genofilename, TextFile.R);
	
			while ( (str = in.readLineElemsReturnReference(TextFile.tab))!= null){

				int indID = hashInds.get(str[0]);
//				System.out.println("Sample id: " + indID + " sample name " + str[0]);
				for (int c = 2; c < str.length; c++) {
					int snpID = snpIDArray[c - 2];
					float dosageFloat = Float.parseFloat(str[c]);
					int dosageInt = Math.round(dosageFloat * 100f);
					if (dosageInt < 0 || dosageInt > 200) {
						System.out.println("Warning, incorrect dosage!:\t" + dosageInt + "\t" + snpID + "\t" + str[c]);
					}
					byte dosageByte = (byte) (Byte.MIN_VALUE + dosageInt);
					byte[] dosage = new byte[1];
					dosage[0] = (byte) dosageByte;
					matrixImputedDosage.setDosage(snpID, indID, dosage);

					char[] alleles = minimacSnpInfo.get(snpID).getAllelesForDossage(dosageFloat);

					byte[] allele1 = new byte[1];
					allele1[0] = (byte) alleles[0];
					byte[] allele2 = new byte[1];
					allele2[0] = (byte) alleles[1];
//					
//					if(minimacSnpInfo.get(snpID).getSnpId().equals("12:30826335")){
//						
//						System.out.println(minimacSnpInfo.get(snpID).getSnpId() + "(" + snpID + "): " + dosageFloat + " (" + dosageInt + "  " + dosageByte + ") " + (char) allele1[0] + "/" + (char) allele2[0]);
//						
//					}
					
					

					matrixGenotype.setAllele1(snpID, indID, allele1);
					matrixGenotype.setAllele2(snpID, indID, allele2);


				}
				System.out.println("Sample\t" + str[0] + "\thas been processed.");
			}
			in.close();
		

		matrixGenotype.close();
		matrixImputedDosage.close();


		System.out.println("\n\n\n");
		System.out.println("Minimac imputed data has been imported.");
		System.out.println("\n\n\n");

	}
	
}
