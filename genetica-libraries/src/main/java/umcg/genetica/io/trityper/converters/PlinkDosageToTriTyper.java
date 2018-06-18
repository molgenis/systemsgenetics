package umcg.genetica.io.trityper.converters;

import umcg.genetica.console.ProgressBar;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.WGAFileMatrixGenotype;
import umcg.genetica.io.trityper.WGAFileMatrixImputedDosage;
import umcg.genetica.io.trityper.util.BaseAnnot;
import umcg.genetica.text.Strings;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

public class PlinkDosageToTriTyper {
	
	
	private ArrayList<String> vArrayListInd;
	private HashMap<String, Integer> vhashInd;
	private HashMap<String, String> familyData;
	private boolean familyDataLoaded;
	
	public void convert(String famfile, String dosefile, String output) throws IOException {
		
		loadFamFile(famfile);
		System.out.println(vArrayListInd.size() + " individuals loaded");
		
		// make list of variants
		System.out.println("Parsing file: " + dosefile);
		ArrayList<String> variants = new ArrayList<>();
		ArrayList<Pair<String, String>> alleles = new ArrayList<Pair<String, String>>();
		TextFile tf = new TextFile(dosefile, TextFile.R);
		String ln = tf.readLine();
		while (ln != null) {
			// chr22:16050036 C A 1.53739 1.53962
			String[] elems = Strings.subsplit(ln, Strings.whitespace, 0, 4);
			String variant = elems[0];
			String allele0 = elems[1];
			String allele1 = elems[2];
			
			
			variants.add(variant);
			if (variants.size() % 10000 == 0) {
				System.out.print("\r" + variants.size() + " variants loaded so far");
			}
			ln = tf.readLine();
		}
		tf.close();
		
		System.out.println(variants.size() + " total");
		
		// write individuals
		TextFile outPhe = new TextFile(output + "PhenotypeInformation.txt", TextFile.W);
		TextFile outinds = new TextFile(output + "Individuals.txt", TextFile.W);
		for (String ind : vArrayListInd) {
			outinds.writeln(ind);
			if (familyData.get(ind) != null) {
				outPhe.write(familyData.get(ind));
			} else {
				outPhe.write(ind + "\tcontrol\tinclude\tunknown" + "\n");
			}
		}
		outinds.close();
		outPhe.close();
		
		TextFile outsnps = new TextFile(output + "SNPs.txt.gz", TextFile.W);
		for (String var : variants) {
			outsnps.writeln(var);
		}
		outsnps.close();
		
		//Process genotypes:
		System.out.println("");
		File fileGenotypeMatrix = new File(output + "/GenotypeMatrix.dat");
		WGAFileMatrixGenotype matrixGenotype = new WGAFileMatrixGenotype(variants.size(), vArrayListInd.size(), fileGenotypeMatrix, false);
		File fileImputedDosageMatrix = new File(output + "/ImputedDosageMatrix.dat");
		WGAFileMatrixImputedDosage matrixImputedDosage = new WGAFileMatrixImputedDosage(variants.size(), vArrayListInd.size(), fileImputedDosageMatrix, false);
		tf = new TextFile(dosefile, TextFile.R);
		String[] elems = tf.readLineElems(Strings.whitespace);
		int snpctr = 0;
		ProgressBar pb = new ProgressBar(variants.size(), "Importing genotypes:");
		while (elems != null) {
			// chr22:16050036 C A 1.53739 1.53962
			int nrinds = elems.length - 3;
			if (nrinds != vArrayListInd.size()) {
				System.out.println("Dosage file has " + nrinds + " while fam has " + vArrayListInd.size());
				System.exit(-1);
			}
			String variant = elems[0];
			String allele0 = elems[1];
			String allele1 = elems[2];
			byte allele0b = BaseAnnot.toByte(allele0);
			byte allele1b = BaseAnnot.toByte(allele1);
			byte[] alleles0 = new byte[vArrayListInd.size()];
			byte[] alleles1 = new byte[vArrayListInd.size()];
			byte[] dosagevals = new byte[vArrayListInd.size()];
			for (int e = 3; e < elems.length; e++) {
				Double dosageval = Double.parseDouble(elems[e]);
				byte indAllele0 = -1;
				byte indAllele1 = -1;
				if (dosageval < 0.5) {
					indAllele0 = allele0b;
					indAllele1 = allele0b;
				} else if (dosageval > 1.5) {
					indAllele0 = allele1b;
					indAllele1 = allele1b;
				} else {
					indAllele0 = allele0b;
					indAllele1 = allele1b;
				}
				
				
				int dosageInt = (int) Math.round(dosageval * 100d);
				if (dosageInt < 0 || dosageInt > 200) {
					System.out.println("Warning, incorrect dosage!:\t" + dosageInt + "\t" + snpctr + "\t" + elems[e]);
					System.exit(-1);
				}
				byte dosageByte = (byte) (Byte.MIN_VALUE + dosageInt);
				
				
				dosagevals[e - 3] = dosageByte;
				alleles0[e - 3] = indAllele0;
				alleles1[e - 3] = indAllele1;
			}
			
			matrixGenotype.setAlleles(snpctr, alleles0, alleles1);
			matrixImputedDosage.setDosages(snpctr, dosagevals);
			snpctr++;
			pb.iterate();
			elems = tf.readLineElems(Strings.whitespace);
		}
		tf.close();
		pb.close();
		matrixGenotype.close();
		matrixImputedDosage.close();
		
		System.out.println("done");
	}
	
	public void loadFamFile(String file) throws IOException {
		vArrayListInd = new ArrayList();
		vhashInd = new HashMap();
		System.out.println("Loading FAM file:\t" + file);
		
		TextFile in = new TextFile(file, TextFile.R);
		
		String line = "";
		familyData = new HashMap<String, String>();
		
		while ((line = in.readLine()) != null) {
			String[] elems = line.split(" ");
			
			if (elems.length >= 6) {
				String sample = elems[1];
				
				String fid = elems[0];
				String pid = elems[2];
				String mid = elems[3];
				String sex = elems[4];
				String phe = elems[5];
				
				if (sex.equals("1")) {
					sex = "male";
				} else if (sex.equals("2")) {
					sex = "female";
				} else {
					sex = "unknown";
				}
				
				if (phe.equals("1")) {
					phe = "control";
				} else if (phe.equals("2")) {
					phe = "case";
				} else {
					phe = "unknown";
				}
				
				familyData.put(sample, sample + "\t" + phe + "\tinclude\t" + sex + "\t" + fid + "\t" + pid + "\t" + mid + "\n");
				
				if (!vhashInd.containsKey(sample)) {
					vhashInd.put(sample, vArrayListInd.size());
					vArrayListInd.add(sample);
				}
				
			}
			
		}
		
		in.close();
		
		familyDataLoaded = true;
		
		
	}
	
}
