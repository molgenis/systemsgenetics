/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.io.trityper.converters;


import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Pattern;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.SNPLoader;
import umcg.genetica.io.trityper.TriTyperGenotypeData;

/**
 *
 * @author harmjan
 */
public class TriTyperToDosageConverter {

    HashMap<String, String> sampleToFamilyId = null;

    public void exportAllSNPs(String datadir, String outputdir, String famfile) throws IOException, Exception {
	TriTyperGenotypeData ds = new TriTyperGenotypeData();
	ds.load(datadir);

	if (!outputdir.endsWith("/")) {
	    outputdir += "/";
	}
	if (!Gpio.isDir(outputdir)) {
	    Gpio.createDir(outputdir);
	}

	if (famfile != null) {
	    loadFamFile(famfile);
	}


	String[] snps = ds.getSNPs();
	for (int chr = 1; chr < 23; chr++) {
	    ArrayList<Integer> snpsToExport = new ArrayList<Integer>();

	    for (int i = 0; i < snps.length; i++) {
		SNP s = ds.getSNPObject(i);
		if (s.getChr() == chr) {
		    snpsToExport.add(i);
		}
	    }

	    TextFile out = new TextFile(outputdir + "output." + chr + ".dose", TextFile.W);
	    exportSNPs(snpsToExport, ds, out);
	    out.close();
	}
    }

    public void exportSubsetOfSNPs(String datadir, String outputdir, String snpSubsetFile, String famfile) throws IOException, Exception {

	TriTyperGenotypeData ds = new TriTyperGenotypeData();
	ds.load(datadir);

	if (!outputdir.endsWith("/")) {
	    outputdir += "/";
	}
	if (!Gpio.isDir(outputdir)) {
	    Gpio.createDir(outputdir);
	}

	if (famfile != null) {
	    loadFamFile(famfile);
	}

	TextFile snpfile = new TextFile(snpSubsetFile, TextFile.R);
	String[] requestedSNPs = snpfile.readAsArray();
	snpfile.close();

	String[] snps = ds.getSNPs();

	HashSet<String> availableSNPs = new HashSet<String>();
	availableSNPs.addAll(Arrays.asList(snps));

	ArrayList<Integer> snpsToExport = new ArrayList<Integer>();

	for (String s : requestedSNPs) {
	    if (availableSNPs.contains(s)) {
            snpsToExport.add(ds.getSnpToSNPId().get(s));
	    }
	}

	TextFile out = new TextFile(outputdir + "output.dose", TextFile.W);
	exportSNPs(snpsToExport, ds, out);
	out.close();

	TriTyperToPedAndMapConverter tgd = new TriTyperToPedAndMapConverter();
	tgd.exportFamFile(ds, outputdir);
    }

    private void exportSNPs(ArrayList<Integer> snpsToExport, TriTyperGenotypeData ds, TextFile out) throws IOException {
	System.out.println("Exporting "+snpsToExport.size()+" snps");
	String header = "SNP\tA1\tA2";
	String[] individuals = ds.getIndividuals();
	int nrSamples = individuals.length;

	int nrIncluded = 0;
	for (int i = 0; i < nrSamples; i++) {
	    if (ds.getIsIncluded()[i]) {
		nrIncluded++;
		String familyId = "1";
		String sample = individuals[i];
		if (sampleToFamilyId != null) {
		    familyId = sampleToFamilyId.get(sample);
		    if (familyId == null) {
			System.out.println("WARNING: sample\t" + sample + "\tis not present in famfile!. Setting family ID to -1!");
			familyId = "-1";
		    }

		}
		header += "\t" + familyId + "\t" + sample;
	    }
	}

	java.text.DecimalFormat df = new java.text.DecimalFormat("0.00", new java.text.DecimalFormatSymbols(java.util.Locale.US));

	SNPLoader loader = ds.createSNPLoader();

	Boolean[] indIncluded = ds.getIsIncluded();
	int numsnps = snpsToExport.size();
	ProgressBar pb = new ProgressBar(numsnps);
	for (int s = 0; s < numsnps; s++) {
	    int snpId = snpsToExport.get(s);
	    SNP snpObj = ds.getSNPObject(snpId);
	    loader.loadGenotypes(snpObj);
	    if (snpObj.getMAF() > 0) {
		loader.loadDosage(snpObj);

		boolean takeComplement = false;
		for (int i = 0; i < individuals.length; i++) {
		    if (indIncluded[i]) {
			double dosagevalue = snpObj.getDosageValues()[i];
			short genotype = snpObj.getGenotypes()[i];
			if (genotype == 0 && dosagevalue > 1) {
			    takeComplement = true;
			    break;
			}
			if (genotype == 2 && dosagevalue < 1) {
			    takeComplement = true;
			    break;
			}
		    }
		}

		int nrChars = snpObj.getName().length() + 2 + 2 + nrIncluded + (nrIncluded*4);

		StringBuilder sb = new StringBuilder(nrChars);
		sb.append(snpObj.getName());
		sb.append("\t");
		byte[] alleles = snpObj.getAlleles();
		if (takeComplement) {
		    sb.append((char) alleles[1]);
		    sb.append("\t");
		    sb.append((char) alleles[0]);
		} else {
		    sb.append((char) alleles[0]);
		    sb.append("\t");
		    sb.append((char) alleles[1]);
		}

		for (int i = 0; i < individuals.length; i++) {
		    if (indIncluded[i]) {
			double dosage = 2 - snpObj.getDosageValues()[i];
			String dosageString = df.format(dosage);
			sb.append("\t");
			sb.append(dosageString);
		    }
		}

		out.writeln(sb.toString());
	    }
	    snpObj.clearGenotypes();
	    pb.iterate();
	}
	loader.close();
    }

    private void loadFamFile(String famfile) throws IOException {

	TextFile fam = new TextFile(famfile, TextFile.R);

	Pattern usedPattern = TextFile.space;
	String[] elems = fam.readLineElems(TextFile.space);
	if (elems.length == 1) {
	    fam.open();
	    elems = fam.readLineElems(TextFile.tab);
	    usedPattern = TextFile.tab;
	    if(elems.length == 1){
		System.out.println("ERROR! parsing FAM file:\t"+famfile+"\t! It should be either space or tab delimited, but both don't seem to fit!");
		throw new IOException();
	    }
	}

	if(elems.length == 0){
	    System.out.println("ERROR! parsing fam file:\t"+famfile+"\t! Are you sure it contains any data?");
	    throw new IOException();
	}

	sampleToFamilyId = new HashMap<String, String>();
	while (elems != null) {
	    sampleToFamilyId.put(elems[1], elems[0]);
	    elems = fam.readLineElems(usedPattern);
	}

	fam.close();

    }
}
