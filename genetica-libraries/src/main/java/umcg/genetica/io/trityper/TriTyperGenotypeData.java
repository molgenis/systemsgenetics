/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.io.trityper;

import gnu.trove.map.hash.TObjectIntHashMap;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.util.ChrAnnotation;
import umcg.genetica.text.Strings;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;

/**
 * @author harmjan
 */
public class TriTyperGenotypeData {
	
	private String[] SNPs;
	private Boolean[] isFemale;
	private Boolean[] isCase;
	private Boolean[] isIncluded;
	private String[] individuals;
	private TObjectIntHashMap<String> individualToId;
	private TObjectIntHashMap<String> snpToSNPId;
	private LinkedList<String[]> alleles;
	private String genotypeFileName;
	private String dosageFileName;
	
	private byte[] chr;
	private int[] chrPos;
	public boolean displayWarnings = true;
	
	public TriTyperGenotypeData(String loc) throws IOException {
		this.load(loc);
	}
	
	public TriTyperGenotypeData() {
	}
	
	public final void load(String loc) throws IOException {
		
		loc = Gpio.formatAsDirectory(loc);
		
		checkFiles(loc);
		
		TextFile t = new TextFile(loc + "Individuals.txt", TextFile.R);
		String[] lineElems = t.readLineElemsReturnReference(TextFile.tab);
		ArrayList<String> alInd = new ArrayList<String>();
		int i = 0;
		individualToId = new TObjectIntHashMap<String>(lineElems.length, 1f, -9);
		while (lineElems != null) {
			String individual = new String(lineElems[0].getBytes("UTF-8"));
			individualToId.put(individual, i);
			alInd.add(individual);
			i++;
			lineElems = t.readLineElemsReturnReference(TextFile.tab);
		}
		t.close();
		
		if (alInd.isEmpty()) {
			System.err.println("ERROR: your dataset does not contain any individuals.");
			System.exit(-1);
		}
		
		
		int numInds = individualToId.size();
		
		setIsFemale(new Boolean[numInds]);
		setIsCase(new Boolean[numInds]);
		setIsIncluded(new Boolean[numInds]);
		
		individuals = new String[numInds];
		
		for (i = 0; i < numInds; i++) {
			getIsFemale()[i] = null;
			getIsCase()[i] = null;
			getIsIncluded()[i] = null;
			individuals[i] = alInd.get(i);
		}
		
		t = new TextFile(loc + "PhenotypeInformation.txt", TextFile.R);
		
		int numMales = 0;
		int numFemales = 0;
		int numCases = 0;
		int numControls = 0;
		int numIncluded = 0;
		
		lineElems = t.readLineElemsReturnReference(TextFile.tab);
		int nrexcludestatusunparseable = 0;
		int nrgenderstatusunparseable = 0;
		int nrcasecontrolstatusunparseable = 0;
		while (lineElems != null) {
			String ind = lineElems[0];
			Integer indId = individualToId.get(ind);
			if (indId != -9) {
				if (lineElems[1].equals("control")) {
					isCase[indId] = false;
					numControls++;
				} else if (lineElems[1].equals("case")) {
					isCase[indId] = true;
					numCases++;
				} else {
					nrcasecontrolstatusunparseable++;
//					if (displayWarnings) {
//						System.err.println("Warning: case/control status unparseable for\t" + lineElems[1] + "\tfor\t" + indId);
//					}
				}
				
				if (lineElems[2].equals("exclude")) {
					isIncluded[indId] = false;
					
				} else if (lineElems[2].equals("include")) {
					isIncluded[indId] = true;
					numIncluded++;
				} else {
					nrexcludestatusunparseable++;
//					if (displayWarnings) {
//						System.err.println("Warning: include/exclude status unparseable\t" + lineElems[2] + "\tfor\t" + indId);
//					}
				}
				
				if (lineElems[3].equals("male")) {
					isFemale[indId] = false;
					numMales++;
				} else if (lineElems[3].equals("female")) {
					isFemale[indId] = true;
					numFemales++;
				} else {
					nrgenderstatusunparseable++;
//					if (displayWarnings) {
//						System.err.println("Warning: gender status unparseable\t" + lineElems[3] + "\tfor\t" + indId);
//					}
				}
			}
			lineElems = t.readLineElemsReturnReference(TextFile.tab);
		}
		t.close();
		
		
		if (nrcasecontrolstatusunparseable > 0) {
			System.out.println("Warning: case/control status unparseable for " + nrcasecontrolstatusunparseable + " individuals.");
		}
		if (nrexcludestatusunparseable > 0) {
			System.out.println("Warning: include/exclude status unparseable for " + nrexcludestatusunparseable + " individuals.");
		}
		if (nrgenderstatusunparseable > 0) {
			System.out.println("Warning: male/female status unparseable for " + nrgenderstatusunparseable + " individuals.");
		}
		
		System.out.println(numInds + " individuals detected, " + numMales + " males, " + numFemales + " females, " + numCases + " cases, " + numControls + " controls, " + numIncluded + " included");
		
		if (numIncluded == 0) {
			System.err.println("ERROR: None of the samples in your dataset will be included. Please check your PhenotypeInformation.txt");
			System.exit(-1);
		}
		
		if (Gpio.exists(loc + "SNPs.txt")) {
			t = new TextFile(loc + "SNPs.txt", TextFile.R, 32 * 1024);
		} else if (Gpio.exists(loc + "SNPs.txt.gz")) {
			t = new TextFile(loc + "SNPs.txt.gz", TextFile.R, 32 * 1024);
		} else {
			throw new FileNotFoundException("SNPs file not found");
		}
		
		
		LinkedList<String> tmpSNP = new LinkedList<String>();
		LinkedList<String[]> tmpSNPAlleles = null; //
		System.out.println("Reading: " + t.getFileName());
		int ctr = 0;
		String[] elems = t.readLineElems(TextFile.tab);
		HashMap<String, String> alleleCache = new HashMap<String, String>();
		
		while (elems != null) {
			if (elems.length > 1) {
				String snp = elems[0];
				String alleleStr = elems[1];
				String[] alleleElems = Strings.comma.split(alleleStr);
				
				for (int a = 0; a < alleleElems.length; a++) {
					String tmpa = alleleElems[a];
					String allele = alleleCache.get(tmpa);
					if (allele == null) {
						String str = new String(tmpa);
						alleleCache.put(alleleElems[a], str);
					}
					alleleElems[a] = allele;
				}
				if (tmpSNPAlleles == null) {
					tmpSNPAlleles = new LinkedList<String[]>();
				}
				tmpSNPAlleles.add(alleleElems);
				tmpSNP.add(snp);
			} else {
				String line = elems[0];
				if (line.trim().length() > 0) {
					tmpSNP.add(line.intern());
				}
			}
			
			elems = t.readLineElems(TextFile.tab);
			ctr++;
			if (ctr % 100000 == 0) {
				System.out.print(ctr + " snps read so far.\r");
			}
		}
		System.out.println();
		t.close();
		
		
		//value if absent now will be -9
		snpToSNPId = new TObjectIntHashMap<String>(tmpSNP.size(), 1f, -9);
		int snpId = 0;
		for (String s : tmpSNP) {
			if (snpToSNPId.containsKey(s)) {
				System.err.println("Warning: duplicate snp " + s + " in " + t.getFileName());
			}
			snpToSNPId.put(s, snpId);
			snpId++;
		}
		
		SNPs = tmpSNP.toArray(new String[0]);
		alleles = tmpSNPAlleles;
		
		System.out.println(SNPs.length + " snps loaded");
		
		if (SNPs.length == 0) {
			System.err.println("ERROR: no SNPs have been detected. Please check your dataset.");
			System.exit(-1);
		}
		
		if (Gpio.exists(loc + "SNPMappings.txt")) {
			t = new TextFile(loc + "SNPMappings.txt", TextFile.R, 100 * 1024);
		} else if (Gpio.exists(loc + "SNPMappings.txt.gz")) {
			t = new TextFile(loc + "SNPMappings.txt.gz", TextFile.R, 100 * 1024);
		} else {
			throw new FileNotFoundException("SNPMappings not found");
		}
		
		System.out.println("Reading: " + t.getFileName());
		lineElems = t.readLineElemsReturnReference(TextFile.tab);
		
		chr = new byte[SNPs.length];
		chrPos = new int[SNPs.length];
		int linenr = 0;
		int nrWoAnnotation = 0;
		while (lineElems != null) {
			if (lineElems.length > 2) {
				int snpNum = snpToSNPId.get(lineElems[2]);
				boolean unknownchr = false;
				boolean unknownchrpos = false;
				if (snpNum != -9 && !"null".equals(lineElems[0]) && !"null".equals(lineElems[1])) {
					
					if (lineElems[0].length() > 0 && !lineElems[0].equals("-1")) {
						chr[snpNum] = ChrAnnotation.parseChr(lineElems[0]);
					} else {
						chr[snpNum] = -1;
						unknownchr = true;
					}
					
					if (lineElems[1].length() > 0 && !lineElems[1].equals("-1")) {
						chrPos[snpNum] = Integer.parseInt(lineElems[1]);
					} else {
						chrPos[snpNum] = -1;
						unknownchrpos = true;
					}
					
					if (unknownchr || unknownchrpos) {
						nrWoAnnotation++;
					}
				}
				
			}
			linenr++;
			if (linenr % 100000 == 0) {
				System.out.print(linenr + " snp mappings read so far.\r");
			}
			lineElems = t.readLineElemsReturnReference(TextFile.tab);
			
		}
		t.close();
		System.out.println();
		System.out.println("Nr of SNPs without annotation: " + nrWoAnnotation + " / " + snpToSNPId.size());
		// open random access file
		setGenotypeFileName(loc + "GenotypeMatrix.dat");
		setDosageFileName(loc + "ImputedDosageMatrix.dat");
		
		checkFileSizes();
		
	}
	
	/**
	 * @return the isFemale
	 */
	public Boolean[] getIsFemale() {
		return isFemale;
	}
	
	/**
	 * @param isFemale the isFemale to set
	 */
	public void setIsFemale(Boolean[] isFemale) {
		this.isFemale = isFemale;
	}
	
	/**
	 * @return the isCase
	 */
	public Boolean[] getIsCase() {
		return isCase;
	}
	
	/**
	 * @param isCase the isCase to set
	 */
	public void setIsCase(Boolean[] isCase) {
		this.isCase = isCase;
	}
	
	/**
	 * @return the isIncluded
	 */
	public Boolean[] getIsIncluded() {
		return isIncluded;
	}
	
	/**
	 * @param isIncluded the isIncluded to set
	 */
	public void setIsIncluded(Boolean[] isIncluded) {
		this.isIncluded = isIncluded;
	}
	
	/**
	 * @return the individualToId
	 */
	public TObjectIntHashMap<String> getIndividualToId() {
		return individualToId;
	}
	
	/**
	 * @param individualToId the individualToId to set
	 */
	public void setIndividualToId(TObjectIntHashMap<String> individualToId) {
		this.individualToId = individualToId;
	}
	
	/**
	 * @return the snpToSNPObject
	 */
	public TObjectIntHashMap<String> getSnpToSNPId() {
		return snpToSNPId;
	}
	
	/**
	 * @param snpToSNPObject the snpToSNPObject to set
	 */
	public void setSnpToSNPObject(TObjectIntHashMap<String> snpToSNPObject) {
		this.snpToSNPId = snpToSNPObject;
	}
	
	/**
	 * @return the genotypeFileName
	 */
	public String getGenotypeFileName() {
		return genotypeFileName;
	}
	
	/**
	 * @param genotypeFileName the genotypeFileName to set
	 */
	public void setGenotypeFileName(String genotypeFileName) {
		this.genotypeFileName = genotypeFileName;
	}
	
	/**
	 * @return the dosageFileName
	 */
	public String getDosageFileName() {
		return dosageFileName;
	}
	
	/**
	 * @param dosageFileName the dosageFileName to set
	 */
	public void setDosageFileName(String dosageFileName) {
		this.dosageFileName = dosageFileName;
	}
	
	/**
	 * @return the SNPs
	 */
	public String[] getSNPs() {
		return SNPs;
	}
	
	/**
	 * @param SNPs the SNPs to set
	 */
	public void setSNPs(String[] SNPs) {
		this.SNPs = SNPs;
	}
	
	public String[] getIndividuals() {
		return individuals;
	}
	
	private void checkFiles(String loc) throws IOException {
		if (!Gpio.exists(loc)) {
			throw new IOException("Error: Directory " + loc + " does not exist.");
		} else if (!Gpio.exists(loc + "PhenotypeInformation.txt")) {
			throw new IOException("Error: Required file " + loc + "PhenotypeInformation.txt does not exist.");
		} else if (!Gpio.exists(loc + "Individuals.txt")) {
			throw new IOException("Error: Required file " + loc + "Individuals.txt does not exist.");
		} else if (!Gpio.exists(loc + "SNPMappings.txt") && !Gpio.exists(loc + "SNPMappings.txt.gz")) {
			throw new IOException("Error: Required file " + loc + "SNPMappings.txt (or SNPMappings.txt.gz) does not exist.");
		} else if (!Gpio.exists(loc + "SNPs.txt") && !Gpio.exists(loc + "SNPs.txt.gz")) {
			throw new IOException("Error: Required file " + loc + "SNPs.txt (or SNPs.txt.gz) does not exist.");
		}
		// GenotypeMatrix.dat
	}
	
	public Integer getIndividualId(String key) {
		return individualToId.get(key);
	}
	
	public Byte getChr(int s) {
		return chr[s];
	}
	
	public int getChrPos(int s) {
		return chrPos[s];
	}
	
	public SNP getSNPObject(int d) {
		SNP out = new SNP();
		out.setId(d);
		out.setChr(chr[d]);
		out.setChrPos(chrPos[d]);
		out.setName(SNPs[d]);
		out.setAlleleEncoding(alleles.get(d));
		return out;
	}
	
	public SNPLoader createSNPLoader() throws IOException {
		return createSNPLoader(1);
	}
	
	public SNPLoader createSNPLoader(int numberOfVariantsToBuffer) throws IOException {
		RandomAccessFile dosageHandle = null;
		RandomAccessFile genotypeHandle = new RandomAccessFile(genotypeFileName, "r");
		if (Gpio.exists(dosageFileName)) {
			dosageHandle = new RandomAccessFile(dosageFileName, "r");
		}
		
		SNPLoader s = new SNPLoader(genotypeHandle, dosageHandle, isIncluded, isFemale, numberOfVariantsToBuffer);
		s.setNumIndividuals(individuals.length);
		return s;
	}
	
	private void checkFileSizes() throws IOException {
		if (Gpio.exists(genotypeFileName)) {
			long expectedfilesize = (long) (SNPs.length * 2) * (long) individuals.length;
			long detectedsize = Gpio.getFileSize(genotypeFileName);
			if (expectedfilesize != detectedsize) {
				throw new IOException("Size of GenotypeMatrix.dat does not match size defined by Individuals.txt and SNPs.txt. Expected size: " + expectedfilesize + " (" + Gpio.humanizeFileSize(expectedfilesize) + ")\tDetected size: " + detectedsize + " (" + Gpio.humanizeFileSize(detectedsize) + ")\tDiff: " + Math.abs(expectedfilesize - detectedsize));
			}
		} else {
			throw new IOException("GenotypeMatrix.dat not detected at location: " + genotypeFileName);
		}
		
		if (Gpio.exists(dosageFileName)) {
			long expectedfilesize = (long) (SNPs.length) * (long) individuals.length;
			long detectedsize = Gpio.getFileSize(dosageFileName);
			
			if (expectedfilesize != detectedsize) {
				throw new IOException("Size of ImputedDosageMatrix.dat does not match size defined by Indivuals.txt and SNPs.txt. Expected size: " + expectedfilesize + " (" + Gpio.humanizeFileSize(expectedfilesize) + ")\tDetected size: " + detectedsize + " (" + Gpio.humanizeFileSize(detectedsize) + ")\tDiff: " + Math.abs(expectedfilesize - detectedsize));
			}
		}
	}
	
	
}
