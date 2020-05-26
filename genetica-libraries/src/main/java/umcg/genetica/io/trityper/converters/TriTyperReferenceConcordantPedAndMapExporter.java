/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.io.trityper.converters;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.containers.SortableSNP;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.SNPLoader;
import umcg.genetica.io.trityper.TriTyperGenotypeData;
import umcg.genetica.io.trityper.util.DetermineLD;

/**
 *
 * @author harmjan
 */
public class TriTyperReferenceConcordantPedAndMapExporter {

    private Double CRTHRESHOLD = 0.95;
    private Double HWEPTHRESHOLD = 0.0001;
    private Double MAFTHRESHOLD = 0.01;
    int numDatasetIndividuals = 0;
    int numReferenceIndividuals = 0;
    private int[] alleleIndex;
    String[] alleleNames = {"A", "C", "G", "T"};
    private TriTyperGenotypeData referenceGenotypeDataset;
    private TriTyperGenotypeData dataGenotypeDataset;
    private TextFile log;
    private DetermineLD ldcalcdataset, ldcalchapmap;
    private HashSet<String> uniquebatches = null;
    private HashMap<String, String> sampleInBatch = null;
    private Object familyData;
    private HashMap<String, String> referenceAlleles = new HashMap<String, String>();
    private HashSet<String> snpsExcluded = new HashSet<String>();
    private static int selectedchr = -1;

    public TriTyperReferenceConcordantPedAndMapExporter() {


		alleleIndex = new int[256];
		for (int a = 0; a < 256; a++) {
			alleleIndex[a] = 4;
		}
		alleleIndex[65] = 0;
		alleleIndex[67] = 1;
		alleleIndex[71] = 2;
		alleleIndex[84] = 3;
		ldcalcdataset = new DetermineLD();
		ldcalchapmap = new DetermineLD();
    }

    public void export(String referenceLoc, String datasetLoc, String batches, String excludesnps, String output) throws IOException {

		if (excludesnps != null) {
			loadExcludeSNPs(excludesnps);
		}
		if (batches != null) {
			loadBatches(batches);
		}


		referenceGenotypeDataset = new TriTyperGenotypeData();
		referenceGenotypeDataset.load(referenceLoc);
		numReferenceIndividuals = referenceGenotypeDataset.getIndividuals().length;

		dataGenotypeDataset = new TriTyperGenotypeData();
		dataGenotypeDataset.load(datasetLoc);
		numDatasetIndividuals = dataGenotypeDataset.getIndividuals().length;

		String[] referenceSNPs = referenceGenotypeDataset.getSNPs();
		String[] datasetSNPs = dataGenotypeDataset.getSNPs();

		ArrayList<SortableSNP> sortedOverlapSNPList = new ArrayList<SortableSNP>();

		int nrnotinref = 0;
		int nrrefwrongmapping = 0;
		int nrwrongmapping = 0;
		int nrdiffchrpos = 0;

		if (!output.endsWith("/")) {
			output += "/";
		}
		log = new TextFile(output + "exportlog.txt", TextFile.W);


		for (int i = 0; i < datasetSNPs.length; i++) {
			String s = datasetSNPs[i];
			if (snpsExcluded.contains(s)) {
			log.writeln(s + "\t excluded by user");
			} else {
			int refSNPId = referenceGenotypeDataset.getSnpToSNPId().get(s);
			if (refSNPId == -9) {
				// SNP IS NOT IN REFERENCE. EXCLUDE!
				log.writeln(s + "\t is not in reference");
				nrnotinref++;
			} else {
				int refsnpchr = referenceGenotypeDataset.getChr(refSNPId);
				int refsnpchrpos = referenceGenotypeDataset.getChrPos(refSNPId);
				if (refsnpchr < 1 || refsnpchrpos < 1 || refsnpchr > 22) {
				// REFERENCE SNP IS NOT PROPERLY MAPPED
				log.writeln(s + "\t is wrongly mapped or on sex chr in reference (chr: " + refsnpchr + ", chrpos: " + refsnpchrpos + ")");
				nrrefwrongmapping++;
				} else {
	// SNP IS IN REFERENCE.. WE NEED TO CHECK LD, CR, ALLELES, ETC
				byte chr = dataGenotypeDataset.getChr(i);
				int chrpos = dataGenotypeDataset.getChrPos(i);
				if (chr < 1 || chrpos < 1) {
					// SNP has no proper mapping
					log.writeln(s + "\t is wrongly mapped in dataset (chr: " + chr + ", chrpos: " + chrpos + ")");
					nrwrongmapping++;
				} else {

					if (refsnpchr == chr && refsnpchrpos == chrpos && (selectedchr < 1 || selectedchr == chr)) {
					sortedOverlapSNPList.add(new SortableSNP(s, i, chr, chrpos, SortableSNP.SORTBY.ID));
					} else {
					// REFERENCE AND DATASET SNP MAP TO DIFFERENT POSITIONS
					log.writeln(s + "\t maps to a different position in reference (dataset chr: " + chr + ", chrpos: " + chrpos + ", reference chr: " + refsnpchr + ", chrpos: " + refsnpchrpos + ")");
					nrdiffchrpos++;
					}
				}
				}
			}
			}

		}

		System.out.println(sortedOverlapSNPList.size() + "\t SNPs overlap between reference and dataset");
		System.out.println(nrnotinref + "\tnot in reference");
		System.out.println(nrrefwrongmapping + "\tref wrongly mapped (chr < 1 || chrpos < 1)");
		System.out.println(nrwrongmapping + "\tdataset wrongly mapped");
		System.out.println(nrdiffchrpos + "\tref and data have different mapping positions");

		Collections.sort(sortedOverlapSNPList);

		SNPLoader referenceLoader = referenceGenotypeDataset.createSNPLoader();
		SNPLoader datasetLoader = dataGenotypeDataset.createSNPLoader();

		// we now have a set of overlapping SNPs, sorted by ID (of the dataset, not the reference)
		int buffersize = 10000;

		int snpcounter = 0;

		int nrExcludedByAllelicDirectionFilter = 0;
		int nrWithSwappedDirection = 0;
		ArrayList<SortableSNP> includedSortedDatasetSNPList = new ArrayList<SortableSNP>();
		HashMap<String, Boolean> alleleswapped = new HashMap<String, Boolean>();

		System.out.println("Checking allelic directions and SNP quality checks");
		ProgressBar pb = new ProgressBar(sortedOverlapSNPList.size());

		int chr1snpspassingqc = 0;
		int chr1snpsnotpassingqc = 0;
		int chr1snpsnotalleleqc = 0;
		int chr1snpstotal = 0;
		int chr1snpspassingalleleqc = 0;

		Boolean[] snpIncluded = new Boolean[sortedOverlapSNPList.size()];

		while (snpcounter < sortedOverlapSNPList.size()) {

			if (snpcounter + buffersize > sortedOverlapSNPList.size()) {
				buffersize = (sortedOverlapSNPList.size() - snpcounter);
			}

			SNP[] datasetSNPBuffer = new SNP[buffersize];
			SNP[] referenceSNPBuffer = new SNP[buffersize];
			boolean[] snpPassesQC = new boolean[buffersize];

			ArrayList<SortableSNP> refbufferSortedSNPs = new ArrayList<SortableSNP>();

			HashMap<String, SortableSNP> datasetToReference = new HashMap<String, SortableSNP>();
			HashMap<String, Integer> referenceToDataset = new HashMap<String, Integer>();


			// load the dataset SNPs
			for (int i = 0; i < buffersize; i++) {
				snpIncluded[i + snpcounter] = false;
				SortableSNP snp = sortedOverlapSNPList.get(i + snpcounter);
				datasetSNPBuffer[i] = dataGenotypeDataset.getSNPObject(snp.id);
				datasetLoader.loadGenotypes(datasetSNPBuffer[i]);

				SortableSNP refSNP = new SortableSNP(snp.name, referenceGenotypeDataset.getSnpToSNPId().get(snp.name), snp.chr, snp.chrpos, SortableSNP.SORTBY.ID);
				datasetToReference.put(snp.name, refSNP);
				referenceToDataset.put(snp.name, i);
				refbufferSortedSNPs.add(refSNP);
			}



			for (int i = 0; i < buffersize; i++) {
				SNP datasetSNP = datasetSNPBuffer[i];
				if (datasetSNP.getChr() == 1) {
					chr1snpstotal++;
				}
				if (datasetSNP.getCR() > CRTHRESHOLD && datasetSNP.getHWEP() > HWEPTHRESHOLD && datasetSNP.getMAF() > MAFTHRESHOLD) {
					snpPassesQC[i] = true;
					if (datasetSNP.getChr() == 1) {
					chr1snpspassingqc++;
					}
				} else {
					if (datasetSNP.getChr() == 1) {
					chr1snpsnotpassingqc++;
					}

					datasetSNP.clearGenotypes();
					datasetSNPBuffer[i] = null;
					// EXCLUDE SNP!
					snpPassesQC[i] = false;
					log.writeln(datasetSNP.getName() + "\t does not pass QC. CR: " + datasetSNP.getCR() + ", HWEP: " + datasetSNP.getHWEP() + ", MAF: " + datasetSNP.getMAF());
				}
			}

			Collections.sort(refbufferSortedSNPs);

			for (int i = 0; i < buffersize; i++) {
				Integer id = refbufferSortedSNPs.get(i).id;
				SNP refSNP = referenceGenotypeDataset.getSNPObject(id);
				referenceLoader.loadGenotypes(refSNP);
				Integer position = referenceToDataset.get(refSNP.getName());
				referenceSNPBuffer[position] = refSNP;
			}

			// compare allelic directions
			for (int i = 0; i < buffersize; i++) {
				SNP datasetSNP = datasetSNPBuffer[i];
				SNP referenceSNP = referenceSNPBuffer[i];

				if (datasetSNP != null && referenceSNP != null) {
					if (datasetSNP.getName().equals(referenceSNP.getName())) {
					if (snpPassesQC[i]) {
						Boolean swapalleles = compareAllelicDirection(referenceSNP, datasetSNP);

						if (swapalleles == null) {
							if (datasetSNP.getChr() == 1) {
								chr1snpsnotalleleqc++;
							}

							nrExcludedByAllelicDirectionFilter++;
							chr1snpsnotpassingqc++;
						} else {

							if (datasetSNP.getChr() == 1) {
								chr1snpspassingalleleqc++;
							}

							// THIS SNP PASSES QC: NOW ON TO LD COMPARISONS
							SortableSNP selected = sortedOverlapSNPList.get(i + snpcounter);
							selected.s = SortableSNP.SORTBY.CHR;
							includedSortedDatasetSNPList.add(sortedOverlapSNPList.get(i + snpcounter));
							alleleswapped.put(datasetSNP.getName(), swapalleles);
							if (swapalleles) {
								nrWithSwappedDirection++;
							}

						}

					}
					} else {
						System.out.println("ERROR names of SNPs not equal:" + datasetSNP.getName() + "\t" + referenceSNP.getName());
						System.exit(0);
					}
				} else if (datasetSNP != null && referenceSNP == null) {
					System.out.println("ERROR: REFERENCE SNP IS NULL WHILE DATASET IS NOT");
				}
			}

			// clear the SNP Buffers
			for (int i = 0; i < buffersize; i++) {

				if (datasetSNPBuffer[i] != null) {
					datasetSNPBuffer[i].clearGenotypes();
				}

				if (referenceSNPBuffer[i] != null) {
					referenceSNPBuffer[i].clearGenotypes();
				}
				datasetSNPBuffer[i] = null;
				referenceSNPBuffer[i] = null;
			}

			snpcounter += buffersize;

			pb.set(snpcounter);
		}

		pb.close();


		System.out.println("pqc " + chr1snpspassingqc + ", npqc " + chr1snpsnotpassingqc + ", naqc " + chr1snpsnotalleleqc + ", tot " + chr1snpstotal + ", totpaqc " + chr1snpspassingalleleqc + "");



		System.out.println(includedSortedDatasetSNPList.size() + "/" + sortedOverlapSNPList.size() + "\t pass allelic direction check.\t" + nrExcludedByAllelicDirectionFilter + " excluded");
		System.out.println(nrWithSwappedDirection + "\t snps with swapped direction alleles");
		Collections.sort(includedSortedDatasetSNPList);

		// compare LD structure

		for (int chr = 1; chr < 23; chr++) {
			ArrayList<SortableSNP> snpsForChr = new ArrayList<SortableSNP>();
			for (int snp = 0; snp < includedSortedDatasetSNPList.size(); snp++) {
			if (includedSortedDatasetSNPList.get(snp).chr == chr) {
				SortableSNP selected = includedSortedDatasetSNPList.get(snp);
				selected.s = SortableSNP.SORTBY.CHRPOS;
				snpsForChr.add(selected);
			}
			}

			Collections.sort(snpsForChr);

			System.out.println(snpsForChr.size() + "\t SNPs for chr " + chr);
			for (int snp = 0; snp < snpsForChr.size() - 1; snp++) {
			if (snpsForChr.get(snp).chrpos > snpsForChr.get(snp + 1).chrpos) {
				System.out.println("SNP LIST NOT PROPERLY SORTED");
				System.exit(0);
			}
			}

			SNP[] datasetSNPObjs = new SNP[snpsForChr.size()];
			boolean[] takeComplement = new boolean[datasetSNPObjs.length];
			SNP[] referenceSNPObjs = new SNP[snpsForChr.size()];

			for (int i = 0; i < datasetSNPObjs.length; i++) {
			SortableSNP datasetSNP = snpsForChr.get(i);
			Integer refSNPId = referenceGenotypeDataset.getSnpToSNPId().get(datasetSNP.name);
			referenceSNPObjs[i] = referenceGenotypeDataset.getSNPObject(refSNPId);
			datasetSNPObjs[i] = dataGenotypeDataset.getSNPObject(datasetSNP.id);
			datasetLoader.loadGenotypes(datasetSNPObjs[i]);

			takeComplement[i] = alleleswapped.get(datasetSNP.name);
			}

			// load the reference snps
			for (int i = 0; i < datasetSNPObjs.length; i++) {
			referenceLoader.loadGenotypes(referenceSNPObjs[i]);
			}

			// test the LD for this Chr, return list of SNPs passing check
			System.out.println("Now testing LD structure for chromosome: " + chr);
			boolean[] snpsPassingQC = testLDStructure(referenceSNPObjs, datasetSNPObjs, takeComplement);

			// clear the reference snps
			for (int i = 0; i < datasetSNPObjs.length; i++) {
			referenceSNPObjs[i].clearGenotypes();
			referenceSNPObjs[i] = null;
			}

			// print the results for this chr
			writePedMapDat(output, chr, datasetSNPObjs, takeComplement, snpsPassingQC);

			// clear the buffer..
			for (int i = 0; i < datasetSNPObjs.length; i++) {
			datasetSNPObjs[i].clearGenotypes();
			datasetSNPObjs[i] = null;
			}

		}
		log.close();

    }

    private Boolean compareAllelicDirection(SNP snpDataObjectHapMap, SNP snpDataObject) throws IOException {

		//SNP has been typed both in dataset and HapMap, do we need to take complementary alleles?

		byte[] allelesbytesreference = snpDataObjectHapMap.getAlleles();
		String alleles = new String(snpDataObject.getAlleles());
		String allelesreference = null;
		Boolean takeComplement = false;
		try {
			allelesreference = new String(snpDataObjectHapMap.getAlleles(), "UTF-8");
		} catch (Exception e) {
		}


		referenceAlleles.put(snpDataObject.getName(), allelesreference);
		boolean allelesOk = true;

		if (allelesbytesreference[0] == 00 || allelesbytesreference[1] == 00) {
	//	    exclude[s] = true;
	//	    excludeReason[s] += "\tHapMap has a null allele HapMap=" + allelesHapMap;
			log.writeln(snpDataObject.getName() + "\t has a null allele in reference: " + allelesreference);
			return null;
		}


		if (allelesreference == null) {
	//	    exclude[s] = true;
	//	    excludeReason[s] += "\tHapMap alleles are NULL ";
			System.out.println("Alleles in HAPMAP are NULL");
			System.exit(1);
		}

		boolean strandForward = true;
		int[] alleleIndex1 = new int[5];

		for (int ind = 0; ind < numDatasetIndividuals; ind++) {
			if (dataGenotypeDataset.getIsIncluded()[ind]) {
				byte[] snpallele1 = snpDataObject.getAllele1();
				byte[] snpallele2 = snpDataObject.getAllele2();

				byte allele1Byte = snpallele1[ind];
				alleleIndex1[alleleIndex[allele1Byte]]++;
				byte allele2Byte = snpallele2[ind];
				alleleIndex1[alleleIndex[allele2Byte]]++;
			}
		}
		int[] alleleIndex2 = new int[5];
		for (int ind = 0; ind < numReferenceIndividuals; ind++) {
			if (referenceGenotypeDataset.getIsIncluded()[ind]) {
			byte[] snpallele1 = snpDataObjectHapMap.getAllele1();
			byte[] snpallele2 = snpDataObjectHapMap.getAllele2();

			byte allele1Byte = snpallele1[ind];
			alleleIndex2[alleleIndex[allele1Byte]]++;
			byte allele2Byte = snpallele2[ind];
			alleleIndex2[alleleIndex[allele2Byte]]++;
			}
		}

		double[] alleleIndexFreq1 = new double[4];
		double[] alleleIndexFreq2 = new double[4];
		int itr = 0;

		boolean issueResolved = false;

		while (!issueResolved) {
			//Take complement alleles when necessary:
			if (!strandForward) {
				int[] alleleIndex1Copy = new int[4];
				System.arraycopy(alleleIndex1, 0, alleleIndex1Copy, 0, 4);
				alleleIndex1[0] = alleleIndex1Copy[3];
				alleleIndex1[1] = alleleIndex1Copy[2];
				alleleIndex1[2] = alleleIndex1Copy[1];
				alleleIndex1[3] = alleleIndex1Copy[0];
			}
			//Determine total number of called alleles:
			int totalCalled1 = 0;
			int totalCalled2 = 0;
			for (int a = 0; a < 4; a++) {
				totalCalled1 += alleleIndex1[a];
				totalCalled2 += alleleIndex2[a];
			}
			//Calculate allele freq:
			for (int a = 0; a < 4; a++) {
				alleleIndexFreq1[a] = (double) alleleIndex1[a] / (double) totalCalled1;
				alleleIndexFreq2[a] = (double) alleleIndex2[a] / (double) totalCalled2;
			}
			//Check whether alleles are identical:
			int nrDifferentAllelesPresent = 0;
			for (int a = 0; a < 4; a++) {
				if (alleleIndexFreq1[a] > 0 || alleleIndexFreq2[a] > 0) {
					nrDifferentAllelesPresent++;
				}
			}
			if (nrDifferentAllelesPresent > 2) {
				strandForward = !strandForward;
			} else {
			//break;
				issueResolved = true;
			}
			itr++;

			if (itr >= 2) {
				if (!issueResolved) {
					//Taking complementary alleles does not resolve anything:
		//		    exclude[s] = true;
		//		    excludeReason[s] += "\tIncompatibleAlleles:Dataset=" + alleles + ",HapMap=" + allelesHapMap;
					log.writeln(snpDataObject.getName() + "\tIncompatibleAlleles:Dataset=" + alleles + ",HapMap=" + allelesreference);

					issueResolved = true;
					return null;
				}
			}
		}

		takeComplement = !strandForward;

		//Check whether allele freq is comparable:
		boolean concordant = true;
		for (int a = 0; a < 4; a++) {
			if (alleleIndexFreq1[a] > 0 && alleleIndexFreq2[a] > 0) {
				if (alleleIndexFreq1[a] > 0.5 && alleleIndexFreq2[a] < 0.5) {
					concordant = false;
				}
				if (alleleIndexFreq1[a] < 0.5 && alleleIndexFreq2[a] > 0.5) {
					concordant = false;
				}
			}
		}

		//If SNP is AT or CG SNP, it can be we had to take complementary allele:
		byte[] snpAlleles = snpDataObject.getAlleles();
		if (snpAlleles[0] + snpAlleles[1] == 65 + 84 || snpAlleles[0] + snpAlleles[1] == 67 + 71) {
			if (!concordant) {
			takeComplement = !takeComplement;
			concordant = true;
			}
		}




		//Output SNPs that have strongly different allele frequencies:
		if (!concordant) {

			String output = "";
			double maxAlleleFrequencyDifference = 0;
			for (int a = 0; a < 4; a++) {
				if (alleleIndexFreq1[a] > 0 || alleleIndexFreq2[a] > 0) {
					int freqData = (int) Math.round(alleleIndexFreq1[a] * 100.0d);
					int freqHapMap = (int) Math.round(alleleIndexFreq2[a] * 100.0d);
					double absDiffFreq = Math.abs(freqData - freqHapMap);
					if (absDiffFreq > maxAlleleFrequencyDifference) {
					maxAlleleFrequencyDifference = absDiffFreq;
					}
					if (output.length() > 0) {
					output += " / ";
					}
					output += alleleNames[a] + ": " + freqData + "%, " + freqHapMap + "%";
				}
			}
			if (maxAlleleFrequencyDifference > 25) {

			log.writeln(snpDataObject.getName() + "\t allele frequencies differ too much with reference: " + output);
			return null;
	//		exclude[s] = true;
	//		excludeReason[s] += "\tAlleleFrequenciesDifferTooMuch:" + output;
			}
	//	    warning[s] += "AlleleFrequenciesDiffer:" + output;
		}

		return takeComplement;
    }

    private boolean[] testLDStructure(SNP[] referenceSNPObjs, SNP[] datasetSNPObjs, boolean[] takeComplement) throws IOException {
		//Complement allele information:
		ProgressBar pb = new ProgressBar(datasetSNPObjs.length);
		boolean[] output = new boolean[datasetSNPObjs.length];
		byte[] complementAllele = new byte[256];
		complementAllele[0] = 0; // Nothing
		complementAllele[84] = 65; //T > A
		complementAllele[65] = 84; //A > T
		complementAllele[48] = 48; //- > -
		complementAllele[67] = 71; //C > G
		complementAllele[71] = 67; //G > C

		int numSNPs = datasetSNPObjs.length;
		int passingLDQC = 0;
		for (int s = 0; s < numSNPs; s++) {

			String rsName = datasetSNPObjs[s].getName();

			SNP datasetSNP1 = datasetSNPObjs[s];
			SNP referenceSNP1 = referenceSNPObjs[s];

			int nrSNPsLDAssessed = 0;
			int nrSNPsLDAssessedPos = 0;
			int nrSNPsLDAssessedNeg = 0;
			String nrSNPsLDAssessedNegLog = "";

			for (int s2 = s - 50; s2 < s + 50; s2++) {
			if (s2 >= 0 && s2 < numSNPs && s2 != s) {

				SNP datasetSNP2 = datasetSNPObjs[s2];
				SNP referenceSNP2 = referenceSNPObjs[s2];

				String rsName2 = datasetSNP2.getName();
				//Do not check flanking SNPS that have A/T or C/G alleles, as this might pose problems:
				byte[] snp2Alleles = referenceSNP2.getAlleles();
				if (snp2Alleles[0] + snp2Alleles[1] != 65 + 84 && snp2Alleles[0] + snp2Alleles[1] != 67 + 71) {

				double rSquared = ldcalcdataset.getRSquared(datasetSNP1, datasetSNP2, dataGenotypeDataset, ldcalcdataset.RETURN_R_SQUARED, ldcalcdataset.INCLUDE_CASES_AND_CONTROLS, false);
				double rSquaredHapMap = ldcalchapmap.getRSquared(referenceSNP1, referenceSNP2, referenceGenotypeDataset, ldcalcdataset.RETURN_R_SQUARED, ldcalcdataset.INCLUDE_CASES_AND_CONTROLS, false);

				if (rSquared > 0.1 && rSquaredHapMap > 0.1) {
					double[] h = new double[4];
					h[0] = ldcalcdataset.h11;
					h[1] = ldcalcdataset.h21;
					h[2] = ldcalcdataset.h12;
					h[3] = ldcalcdataset.h22;
					double[] hHapMap = new double[4];
					hHapMap[0] = ldcalchapmap.h11;
					hHapMap[1] = ldcalchapmap.h21;
					hHapMap[2] = ldcalchapmap.h12;
					hHapMap[3] = ldcalchapmap.h22;
					byte[][] hAlleles = new byte[4][2];

					byte[] snpDataAlleles1 = datasetSNP1.getAlleles();
					byte[] snpDataAlleles2 = datasetSNP2.getAlleles();

					hAlleles[0][0] = snpDataAlleles1[0];
					hAlleles[1][0] = snpDataAlleles1[0];
					hAlleles[2][0] = snpDataAlleles1[1];
					hAlleles[3][0] = snpDataAlleles1[1];
					hAlleles[0][1] = snpDataAlleles2[0];
					hAlleles[1][1] = snpDataAlleles2[1];
					hAlleles[2][1] = snpDataAlleles2[0];
					hAlleles[3][1] = snpDataAlleles2[1];
					if (takeComplement[s]) {
					for (int a = 0; a < 4; a++) {
						hAlleles[a][0] = complementAllele[hAlleles[a][0]];
					}
					}
					if (takeComplement[s2]) {
					for (int a = 0; a < 4; a++) {
						hAlleles[a][1] = complementAllele[hAlleles[a][1]];
					}
					}
					byte[][] hAllelesHapMap = new byte[4][2];

					byte[] snpHapMapAlleles1 = referenceSNP1.getAlleles();
					byte[] snpHapMapAlleles2 = referenceSNP2.getAlleles();

					hAllelesHapMap[0][0] = snpHapMapAlleles1[0];
					hAllelesHapMap[1][0] = snpHapMapAlleles1[0];
					hAllelesHapMap[2][0] = snpHapMapAlleles1[1];
					hAllelesHapMap[3][0] = snpHapMapAlleles1[1];
					hAllelesHapMap[0][1] = snpHapMapAlleles2[0];
					hAllelesHapMap[1][1] = snpHapMapAlleles2[1];
					hAllelesHapMap[2][1] = snpHapMapAlleles2[0];
					hAllelesHapMap[3][1] = snpHapMapAlleles2[1];
					double[] xVals = new double[4];
					double[] yVals = new double[4];
					int itr = 0;
					for (int a = 0; a < 4; a++) {
					for (int b = 0; b < 4; b++) {
						if (hAlleles[a][0] == hAllelesHapMap[b][0] && hAlleles[a][1] == hAllelesHapMap[b][1]) {
						xVals[itr] = h[a];
						yVals[itr] = hHapMap[b];
						itr++;
						}
					}
					}

					double correlation = JSci.maths.ArrayMath.correlation(xVals, yVals);
					nrSNPsLDAssessed++;
					if (correlation < 0) {
					nrSNPsLDAssessedNeg++;
					nrSNPsLDAssessedNegLog += "\t" + rsName2;
					}
					if (correlation > 0) {
					nrSNPsLDAssessedPos++;
					}
				}
				}
			}
			}

			if (nrSNPsLDAssessedNeg > nrSNPsLDAssessedPos) {
			//There are SNPs with whom this SNP is in LD, but the haplotype frequencies are completeley different:
			output[s] = false;
			log.writeln(rsName + "\tFlankingSNPsLDPatternsDifferent.NrAssessed:" + nrSNPsLDAssessed + ",Pos:" + nrSNPsLDAssessedPos + ",Neg:" + nrSNPsLDAssessedNeg + ",NegSNPs:" + nrSNPsLDAssessedNegLog);
			} else {
			passingLDQC++;
			output[s] = true;
			}
			pb.iterate();
		}

		pb.close();

		System.out.println(passingLDQC+" out of "+datasetSNPObjs.length+" snps pass LD QC");
		return output;
	}

    private void writePedMapDat(String output, int chr, SNP[] datasetSNPObjs, boolean[] takeComplement, boolean[] snpsPassingQC) throws IOException {

	//if there are no snps do not write output. Allows per chromosome conversion
	if(datasetSNPObjs.length == 0){
		System.out.println("No SNPs for chr " + chr);
		return;
	}
		
	byte[] complementAllele = new byte[256];
	complementAllele[0] = 0; // Nothing
	complementAllele[84] = 65; //T > A
	complementAllele[65] = 84; //A > T
	complementAllele[48] = 48; //- > -
	complementAllele[67] = 71; //C > G
	complementAllele[71] = 67; //G > C

	int numbatches = 0;
	HashMap<String, Integer> batchtofile = new HashMap<String, Integer>();

	String[] batchnames = null;

	if (uniquebatches != null) {
	    numbatches = uniquebatches.size();
	    batchnames = new String[numbatches];
	    uniquebatches.toArray(batchnames);
	    for (int i = 0; i < batchnames.length; i++) {
			batchtofile.put(batchnames[i], i);
			System.out.println(i + "\t" + batchnames[i]);
	    }
	} else {
		//only one single batch per chromosome
		numbatches = 1;
	}

//	System.exit(0);

	TextFile[] batchpedfiles = new TextFile[numbatches];

	if(uniquebatches != null){
		for (int i = 0; i < batchpedfiles.length; i++) {
			String batchname = batchnames[i];
			batchpedfiles[i] = new TextFile(output + "/chr" + chr + "-" + batchname + ".ped", TextFile.W);
		}
	} else {
		//no batches. Write output per chromosome
		batchpedfiles[0] = new TextFile(output + "/chr" + chr + ".ped", TextFile.W);
	}



	// write PED files
	String[] individuals = dataGenotypeDataset.getIndividuals();
	
	individuals:
	for (int i = 0; i < individuals.length; i++) {
		
		String individual = individuals[i];
		Integer filehandle;
		
		if(uniquebatches != null){
			
			String batch = sampleInBatch.get(individual);
			if (batch == null) {
				System.out.println(individual + "\t not included in any batch");
				continue individuals;
			} else {
				filehandle = batchtofile.get(batch);
				if (filehandle == null) {
					System.out.println(batch + "\t has no handle?");
					continue individuals;
				}
			}
			
		} else {
			//no batches so use first that contains handle for whole chromosome
			filehandle = 0;
		}
		
		StringBuilder individualOutputString = new StringBuilder();

		String affectionStatus = "-9";
		if (dataGenotypeDataset.getIsCase()[i] == null) {
		affectionStatus = "-9";
		} else if (dataGenotypeDataset.getIsCase()[i]) {
		affectionStatus = "2";
		} else if (!dataGenotypeDataset.getIsCase()[i]) {
		affectionStatus = "1";
		}

		String sex = "-9";

		if (dataGenotypeDataset.getIsFemale()[i] == null) {
		sex = "-9";
		} else if (dataGenotypeDataset.getIsFemale()[i]) {
		sex = "2";
		} else if (!dataGenotypeDataset.getIsFemale()[i]) {
		sex = "1";
		}

		if (individual.contains(" ")) {
		individual = individual.replaceAll(" ", "-");
		}

		if (familyData != null) {
//		if (familyData.get(individual) != null) {
//		    outputInd[indCounter] += familyData.get(individual);
//		} else {
//		    outputInd[indCounter] += "1 " + individual + " 0 0 " + sex + " " + affectionStatus;
//		}
		} else {
		individualOutputString.append("1 ").append(individual).append(" 0 0 ").append(sex).append(" ").append(affectionStatus);
		}
		// build string


		for (int snp = 0; snp < datasetSNPObjs.length; snp++) {
			
			if (snpsPassingQC[snp]) {
				byte[] content = new byte[4];
				SNP snpDataObject = datasetSNPObjs[snp];

				content[0] = 32;
				byte[] snpAlleles1 = snpDataObject.getAllele1();
				byte[] snpAlleles2 = snpDataObject.getAllele2();

				byte value1 = snpAlleles1[i];
				byte value2 = snpAlleles2[i];
				if (value1 == 0) {
				value1 = 48;
				}
				if (value2 == 0) {
				value2 = 48;
				}

	//		byte[] rValues = snpDataObject.getRValues();
	//		byte[] thetaValues = snpDataObject.getThetaValues();
	//		if (genotypeDataset.rawDataAvailable() && rValues[ind] == 0 && thetaValues[ind] == 0) {
	//		    value1 = 48;
	//		    value2 = 48;
	//		}

				if (takeComplement[snp]) {
				value1 = complementAllele[value1];
				value2 = complementAllele[value2];
				}

				content[1] = value1;
				content[2] = 32;
				content[3] = value2;

				individualOutputString.append(new String(content));
			}

		}



		batchpedfiles[filehandle].writeln(individualOutputString.toString());

	}

	for (int i = 0; i < batchpedfiles.length; i++) {
	    batchpedfiles[i].close();
	}

	// write MAP files

	// write FAM files

	TextFile mapfile = new TextFile(output + "/chr" + chr + ".map", TextFile.W);
	TextFile datfile = new TextFile(output + "/chr" + chr + ".dat", TextFile.W);


	datfile.writeln("A Status");

	for (int s = 0; s < datasetSNPObjs.length; s++) {
	    if (snpsPassingQC[s]) {
		SNP datasetSNP = datasetSNPObjs[s];
		String rsName = datasetSNP.getName();
		Integer snpID = datasetSNP.getId();

		byte snpchr = datasetSNP.getChr();
		int snpchrpos = datasetSNP.getChrPos();
		mapfile.writeln(chr + "\t" + rsName + "\t0\t" + snpchrpos);
		datfile.writeln("M " + rsName);
		//out.write(chr + " " + rsName + " " + pos + " " + strand + "\n");
	    }
	}
	datfile.close();
	mapfile.close();



	ArrayList<SortableSNP> sortedSNPs = new ArrayList<SortableSNP>();
	for (int s = 0; s < datasetSNPObjs.length; s++) {
	    if (snpsPassingQC[s]) {
		SNP datasetSNP = datasetSNPObjs[s];
		sortedSNPs.add(new SortableSNP(datasetSNP.getName(), s, datasetSNP.getChr(), datasetSNP.getChrPos(), SortableSNP.SORTBY.CHRPOS));
	    }
	}

//	Collections.sort(sortedSNPs);

	TextFile bglfile = new TextFile(output + "/chr" + chr + ".markersbeagleformat", TextFile.W);
	for (int s = 0; s < sortedSNPs.size(); s++) {
	    SortableSNP sortedSNP = sortedSNPs.get(s);

	    String rsName = sortedSNP.name;
	    int pos = sortedSNP.chrpos;
	    String allelesHapMap = referenceAlleles.get(rsName);
	    if (allelesHapMap != null) {
		bglfile.writeln(rsName + " " + pos + " " + allelesHapMap.substring(0, 1) + " " + allelesHapMap.substring(1, 2));
	    }
	}
	bglfile.close();

    }

    private void loadBatches(String batches) throws IOException {
	TextFile tf = new TextFile(batches, TextFile.R);
	String[] elems = tf.readLineElemsReturnReference(TextFile.tab);
	uniquebatches = new HashSet<String>();
	sampleInBatch = new HashMap<String, String>();
	while (elems != null) {
	    uniquebatches.add(elems[0]);
	    sampleInBatch.put(elems[1], elems[0]);
	    elems = tf.readLineElemsReturnReference(TextFile.tab);
	}
	tf.close();
	System.out.println("Loaded " + uniquebatches.size() + " batches," + sampleInBatch.size() + " samples");
    }

    private void loadExcludeSNPs(String excludesnps) throws IOException {
	TextFile tf = new TextFile(excludesnps, TextFile.R);
	ArrayList<String> snps = tf.readAsArrayList();
	snpsExcluded.addAll(snps);
	tf.close();
    }
}
