/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.io.trityper.converters;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.WGAFileMatrixGenotype;
import umcg.genetica.io.trityper.util.BaseAnnot;
import umcg.genetica.io.trityper.util.ChrAnnotation;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class VCFToTriTyper {

	
	Integer[] colToIndId = null;
	private HashMap<String, Integer> individualMap;
	private ArrayList<String> individuals;
	private HashMap<String, Integer> snpMap;
	private HashMap<String, Byte> snpChrMap;
	private HashMap<String, Integer> snpChrPosMap;
	private ArrayList<String> snpList;
	private int multiallelicSNPsExcluded;
	private int snpctr;
	private int finalNrInds;
	private Pattern snppattern; 
	
	/**
	 * Parse and assume all variants are SNPs
	 * 
	 * @param dir
	 * @param oudir 
	 */
	public void parse(String dir, String oudir) throws IOException{
		parse(dir, oudir, null);
	}

	/**
	 * Parser that can filter variants to only include SNPs
	 * 
	 * @param dir
	 * @param outdir
	 * @param snpPatternString pattern INFO column of VCF must contain. Use null for no check
	 * @throws IOException 
	 */
	public void parse(String dir, String outdir, String snpPatternString) throws IOException {

		if(snpPatternString == null){
			this.snppattern = null;
			System.out.println("All variants are assumed to be SNPs");
		} else {
			this.snppattern = Pattern.compile(snpPatternString);
			System.out.println("Variants are filtered for SNPs using this pattern on the INFO column: " + snpPatternString);
		}
		
		
		if (!Gpio.exists(dir)) {
			throw new IOException("Error: could not find dir: " + dir);
		}


		if (!outdir.endsWith("/")) {
			outdir += "/";
		}

		Gpio.createDir(outdir);


		String[] files = Gpio.getListOfFiles(dir);
		ArrayList<String> finalFiles = new ArrayList<String>();
		for (String file : files) {
			if (file.endsWith(".vcf")) {
				finalFiles.add(file);
			}
		}

		files = finalFiles.toArray(new String[0]);



		System.out.println("Found " + files.length + " vcf files");
		if (files.length == 0) {
			System.exit(0);
		}


		for (String file : files) {
			parseFile(dir, file, outdir);
		}



	}

	private void parseFile(String dir, String file, String outdir) throws IOException {
		String[] fileNameElems = Strings.dot.split(file);

		String chrNum = "";
		for (String f : fileNameElems) {
			if (f.startsWith("chr")) {
				chrNum = "Chr" + f.substring(3);
			}

		}

		file = dir + "/" + file;


		outdir += chrNum + "/";
		Gpio.createDir(outdir);

		System.out.println("Will write to " + outdir);

		snpMap = new HashMap<String, Integer>();
		snpChrMap = new HashMap<String, Byte>();
		snpChrPosMap = new HashMap<String, Integer>();

		snpctr = 0;

		int filecounter = 0;
		HashMap<String, Integer> individualsInFilesCounter = new HashMap<String, Integer>();

		HashSet<String> uniqueIndividuals = new HashSet<String>();


		System.out.println("Parsing file: " + file);
		// first get an overview of the number of individuals in each of the files...
		TextFile tf = new TextFile(file, TextFile.R);
		int indcounter = 0;
		String[] elems = tf.readLineElems(TextFile.tab);

		while (elems != null) {
			if (elems[0].startsWith("##")) {
				// comment line
			} else if (elems[0].startsWith("#CHROM")) {
				// this is the header
				if (elems.length > 9) {
					for (int i = 9; i < elems.length; i++) {
						String sample = elems[i];
						Integer numFilesPresent = individualsInFilesCounter.get(sample);
						if (numFilesPresent == null) {
							numFilesPresent = 0;
						}

						numFilesPresent++;
						individualsInFilesCounter.put(sample, numFilesPresent);

						uniqueIndividuals.add(sample);
					}
				}
				break;
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();


		String[] uniqueIndividualsArray = uniqueIndividuals.toArray(new String[0]);

		individualMap = new HashMap<String, Integer>();
		individuals = new ArrayList<String>();
		int indIdCtr = 0;

		for (String individual : uniqueIndividualsArray) {
			Integer ctr = individualsInFilesCounter.get(individual);
//	    if (ctr != null && ctr == files.length) {
			// include individual
			individualMap.put(individual, indIdCtr);
			individuals.add(individual);
			indIdCtr++;
//	    } else {
//		System.out.println("Individual " + individual + " not present in all files! This individual will not be included in the final GenotypeMatrix.dat");
//	    }
		}

		System.out.println("Total number of detected individuals: " + individuals.size());

		System.out.println("Now writing individuals to output directory");
		TextFile indOut = new TextFile(outdir + "Individuals.txt", TextFile.W);
		TextFile indPhenoOut = new TextFile(outdir + "PhenotypeInformation.txt", TextFile.W);
		for (String ind : individuals) {
			indOut.writeln(ind);
			indPhenoOut.writeln(ind + "\tunknown\tinclude\tunknown");
		}
		indOut.close();
		indPhenoOut.close();



		multiallelicSNPsExcluded = 0;

		snpList = new ArrayList<String>();

		System.out.println("Parsing file: " + file);
		int lnctr = 0;
		tf = new TextFile(file, TextFile.R);

		colToIndId = null;
		// #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  HG00096
		elems = tf.readLineElemsReturnObjects(TextFile.tab);
		int nrSNPs = 0;
		while (elems != null) {
			if (elems[0].startsWith("##")) {
				// comment line
				// parse comment
			} else if (elems[0].startsWith("#CHROM")) {
				// parse chromline
				parseHeaderLine(elems);
			} else {
				// parse SNP line
				parseVCFSNPLine(elems, true);
			}
			elems = tf.readLineElemsReturnObjects(TextFile.tab);

			lnctr++;

			if (lnctr % 500000 == 0) {
				System.out.println("Parsed\t" + snpList.size() + "\tsnps.");
			}
		}

		tf.close();

		filecounter++;
		System.out.println(snpMap.size() + "\tsnps detected");
		System.out.println(multiallelicSNPsExcluded + "\tmulti allelic SNPs excluded.");



		System.out.println("Final totals: ");
		System.out.println(snpMap.size() + "\tsnps detected");
		System.out.println(multiallelicSNPsExcluded + "\tmulti allelic SNPs excluded.");

		System.out.println("Now writing snps to output directory!");
		String[] availableSNPs = snpList.toArray(new String[snpList.size()]);
		TextFile snpout = new TextFile(outdir + "SNPs.txt", TextFile.W);
		TextFile snpmapout = new TextFile(outdir + "SNPMappings.txt", TextFile.W);
		for (String snp : availableSNPs) {
			snpout.writeln(snp);
			snpmapout.writeln(snpChrMap.get(snp) + "\t" + snpChrPosMap.get(snp) + "\t" + snp);
		}
		snpout.close();
		snpmapout.close();

		finalNrInds = individuals.size();
		// create genotypeMatrix file
		String outfilename = outdir + "GenotypeMatrix.dat";
		WGAFileMatrixGenotype genotypefile = new WGAFileMatrixGenotype(snpctr, individuals.size(), new File(outfilename), false);

		filecounter = 0;

		ProgressBar pb = new ProgressBar(lnctr, "writing genotypes from file: " + file);
		tf = new TextFile(file, TextFile.R);
		elems = tf.readLineElems(TextFile.tab);

		// create a map, mapping individual to certain column
		colToIndId = null;

		lnctr = 0;
		// #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  HG00096
		snpctr = 0;
		while (elems != null) {

			if (elems[0].startsWith("##")) {
				// comment line
			} else if (elems[0].startsWith("#CHROM")) {
				// this is the header
				parseHeaderLine(elems);
			} else {
				// snp line
				String snp = elems[2];
				if (snp.equals(".")) {
					snp = elems[0] + "_" + elems[1];
				}
				Integer snpid = snpMap.get(snp);
				if (snpid != null) {
					Pair<byte[], byte[]> p = parseVCFSNPLine(elems, false);
					if (p != null) {
						genotypefile.setAlleles(snpid, p.getLeft(), p.getRight());
						snpctr++;
					}

				}

			}
			elems = tf.readLineElems(TextFile.tab);
			pb.set(lnctr);
			lnctr++;
		}

		pb.close();
		filecounter++;
		tf.close();
		System.out.println("");


		genotypefile.close();
	}
	
	private static final Pattern zero = Pattern.compile("0");
	private static final Pattern one = Pattern.compile("1");

	private Pair<byte[], byte[]> parseVCFSNPLine(String[] elems, boolean inventory) {
		// snp line
		String info = elems[7];
		String[] infoelems = Strings.semicolon.split(info);

		boolean isSNP;
		if(snppattern == null ){
			isSNP = true;
		} else {
			isSNP = false;
			for (String infoelem : infoelems) {
				if (snppattern.matcher(infoelem).matches()) {
					isSNP = true;
				}
			}
		}

		if (isSNP) {
			String snp = elems[2];

			if (Strings.dot.matcher(snp).matches()) {
				snp = elems[0] + "_" + elems[1];
			}
			String ref = elems[3];
			String alt = elems[4];
			if (Strings.comma.split(alt).length == 1 && Strings.comma.split(ref).length == 1) {

				byte refb = BaseAnnot.toByte(elems[3]);
				byte altb = BaseAnnot.toByte(elems[4]);

				if (refb == 0 || altb == 0) {
					System.err.println("WARNING: could not properly parse reference or alternative allele for snp\t" + snp + "\t" + elems[3] + "-" + elems[4]);
				} else if (!inventory && snpMap.containsKey(snp)) {
					byte[] allele1 = new byte[finalNrInds];
					byte[] allele2 = new byte[finalNrInds];
					for (int i = 9; i < elems.length; i++) {
						Integer indId = colToIndId[i];
						if (indId != null) {
							// write alleles
							String[] gtElems = Strings.colon.split(elems[i]);
							String[] genotypes = Strings.pipe.split(gtElems[0]);

							if (genotypes.length == 1) {
								genotypes = Strings.forwardslash.split(gtElems[0]);
							}
							if (genotypes.length == 1 || genotypes.length > 2) {
								System.err.println("WARNING: genotype could not be parsed for sample " + individuals.get(indId) + "\t" + gtElems[0]);
							} else {
								if (zero.matcher(genotypes[0]).matches()) {
									allele1[indId] = refb;
								} else if (one.matcher(genotypes[0]).matches()) {
									allele1[indId] = altb;
								} else if (Strings.dot.matcher(genotypes[0]).matches())  {
									allele1[indId] = 0;						
								}else {
									System.err.println("Could not parse allele1 of genotype for sample " + individuals.get(indId) + "\t" + gtElems[0]);
								}

								if (zero.matcher(genotypes[1]).matches()) {
									allele2[indId] = refb;
								} else if (one.matcher(genotypes[1]).matches()) {
									allele2[indId] = altb;
								} else if (Strings.dot.matcher(genotypes[1]).matches())  {
									allele2[indId] = 0;		
								} else {
									System.err.println("Could not parse allele2 of genotype for sample " + individuals.get(indId) + "\t" + gtElems[0]);
								}


							}
						}
					}
					return new Pair<byte[], byte[]>(allele1, allele2);
				} else if (inventory) {
					if (snpMap.containsKey(snp)) {
						System.err.println("WARNING: " + snp + " already parsed?\n" + Strings.concat(elems, Strings.tab));
					} else {
						Byte chr = ChrAnnotation.parseChr(elems[0]);
						Integer pos = Integer.parseInt(elems[1]);

						snpList.add(snp);
						snpMap.put(snp, snpctr);
						snpChrMap.put(snp, chr);
						snpChrPosMap.put(snp, pos);
						snpctr++;
					}
				}
			} else {
				multiallelicSNPsExcluded++;
				System.out.println("SNP " + snp + " is multi-allelic, therefore exlcuding it! " + ref + "\t" + alt);
			}
		}


		return null;
	}

	private void parseHeaderLine(String[] elems) {

		colToIndId = new Integer[elems.length];
		for (int i = 9; i < elems.length; i++) {
			String sample = elems[i];
			Integer indId = individualMap.get(sample);
			colToIndId[i] = indId;
		}

	}
}
