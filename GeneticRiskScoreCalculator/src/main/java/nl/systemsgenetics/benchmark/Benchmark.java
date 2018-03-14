package nl.systemsgenetics.benchmark;

import nl.systemsgenetics.simplegeneticriskscorecalculator.Main;
import org.apache.commons.cli.*;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.SNPLoader;
import umcg.genetica.io.trityper.TriTyperGenotypeData;
import umcg.genetica.text.Strings;
import umcg.genetica.util.Primitives;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.logging.Level;
import java.util.logging.Logger;

public class Benchmark {
	
	enum MODE {
		CORRELATE, PARSELOG, COMBINE, EXTRACT, CONCATLOG, FILTER
	}
	
	public static void main(String[] args) {
		Benchmark b = new Benchmark();
		
		CommandLineParser parser = new GnuParser();
		Options options = new Options();
		
		Option fileOut = OptionBuilder.withArgName("path").hasArg().withDescription("Location (folder) for the output.").withLongOpt("output").create("o");
		Option mode = OptionBuilder.withArgName("path").hasArg().withDescription("Set mode (filter | correlate | parselog | concatlog | combine | extract)").withLongOpt("mode").create("m");
		Option fileIn = OptionBuilder.withArgName("path").hasArg().withDescription("Location (folder) for the input.").withLongOpt("input").create("i");
		Option fileIn2 = OptionBuilder.withArgName("path").hasArg().withDescription("Location (folder) for the input 2.").withLongOpt("input2").create("i2");
		Option genotypeIn = OptionBuilder.withArgName("path").hasArg().withDescription("Location for the reference data.").withLongOpt("genotypes").create("g");
		Option traitIn = OptionBuilder.withArgName("path").hasArg().withDescription("Location for the PRS data.").withLongOpt("traits").create("t");
		Option traitToSNPCpl = OptionBuilder.withArgName("path").hasArg().withDescription("Location for the file that links traits to snps.").withLongOpt("traittosnpcoupling").create("ttsc");
		
		options.addOption(fileOut)
				.addOption(mode)
				.addOption(fileIn)
				.addOption(fileIn2)
				.addOption(genotypeIn)
				.addOption(traitIn)
				.addOption(traitToSNPCpl);
		
		CommandLine cmd;
		try {
			cmd = parser.parse(options, args);
			HelpFormatter formatter = new HelpFormatter();
			
			String in = null;
			String out = null;
			String gt = null;
			String ttsc = null;
			String tr = null;
			
			MODE m = null;
			
			if (cmd.hasOption("m")) {
				// initialise the member variable
				String modestr = cmd.getOptionValue("m");
				if (modestr.equals("filter")) {
					m = MODE.FILTER;
				} else if (modestr.equals("correlate")) {
					m = MODE.CORRELATE;
				} else if (modestr.equals("parselog")) {
					m = MODE.PARSELOG;
				} else if (modestr.equals("combine")) {
					m = MODE.COMBINE;
				} else if (modestr.equals("extract")) {
					m = MODE.EXTRACT;
				} else if (modestr.equals("concatlog")) {
					m = MODE.CONCATLOG;
				}
			}
			if (cmd.hasOption("i")) {
				in = cmd.getOptionValue("i");
			}
			String in2 = null;
			if (cmd.hasOption("i2")) {
				in2 = cmd.getOptionValue("i2");
			}
			if (cmd.hasOption("o")) {
				out = cmd.getOptionValue("o");
			}
			if (cmd.hasOption("g")) {
				gt = cmd.getOptionValue("g");
			}
			if (cmd.hasOption("t")) {
				tr = cmd.getOptionValue("t");
			}
			if (cmd.hasOption("ttsc")) {
				ttsc = cmd.getOptionValue("ttsc");
			}
			
			if (m == null) {
				formatter.printHelp("benchmark", options);
			} else if (m.equals(MODE.FILTER)) {
				b.traitfilter(gt, in, out);
			} else if (m.equals(MODE.CORRELATE)) {
				b.correlate(tr, gt, ttsc, out);
			} else if (m.equals(MODE.PARSELOG)) {
				b.logparser(ttsc, in, out);
			} else if (m.equals(MODE.COMBINE)) {
				b.combinefiles(in, in2, out);
			} else if (m.equals(MODE.EXTRACT)) {
				b.extractSNP(gt, in, out);
			} else if (m.equals(MODE.CONCATLOG)) {
				b.logconcat(in, in2, out);
			}
			
		} catch (ParseException ex) {
			Logger.getLogger(Main.class.getName()).log(Level.SEVERE, null, ex);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public void logconcat(String logdir, String traitlist, String out) throws IOException {
		TextFile tf = new TextFile(traitlist, TextFile.R);
		ArrayList<String> listoftraits = tf.readAsArrayList();
		tf.close();
		
		TextFile outf = new TextFile(out, TextFile.W);
		int tctr = 0;
		for (String t : listoftraits) {
			ArrayList<String> snpids = new ArrayList<>();
			
			for (int c = 1; c < 23; c++) {
				String logfile = logdir + t + "Chr" + c + ".log";
				if (Gpio.exists(logfile)) {
					System.out.print(tctr + "/" + listoftraits.size() + "\t" + logfile + "\r");
					TextFile tfl = new TextFile(logfile, TextFile.R);
					String ln = tfl.readLine();
					while (ln != null) {
						if (ln.startsWith("rs")) {
							String[] logelems = ln.split("\t");
							String rs = logelems[0];
							snpids.add(rs);
						}
						ln = tfl.readLine();
					}
					tfl.close();
				}
				
			}
			String snpstr = "NA";
			if (!snpids.isEmpty()) {
				snpstr = Strings.concat(snpids, Strings.semicolon);
			}
			
			outf.writeln(t + "\t" + snpids.size() + "\t" + snpstr);
			
			tctr++;
		}
		outf.close();
	}
	
	public void logparser(String traitToSNPCoupling, String logdir, String out) throws IOException {
		TextFile tf = new TextFile(traitToSNPCoupling, TextFile.R);
		String[] elems = tf.readLineElems(TextFile.tab);
		HashMap<String, String> traitToSNP = new HashMap<>();
		ArrayList<String> traits = new ArrayList<>();
		while (elems != null) {
			String trait = elems[0];
			String snp = elems[1];
			traits.add(trait);
			if (!traitToSNP.containsKey(trait)) {
				traitToSNP.put(trait, snp);
			} else {
				System.out.println("Trait " + trait + " is a duplicate?");
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		
		TextFile outf = new TextFile(out, TextFile.W);
		outf.writeln("Trait\tInputSNP\tOutputSNP");
		for (String t : traits) {
			ArrayList<String> snpids = new ArrayList<>();
			for (int c = 1; c < 23; c++) {
				String logfile = logdir + t + "Chr" + c + ".log";
				if (Gpio.exists(logfile)) {
					TextFile tfl = new TextFile(logfile, TextFile.R);
					String ln = tfl.readLine();
					while (ln != null) {
						if (ln.startsWith("rs")) {
							String[] logelems = ln.split("\t");
							String rs = logelems[0];
							snpids.add(rs);
						}
						ln = tfl.readLine();
					}
					tfl.close();
				}
			}
			
			
			String snpstr = "NA";
			if (!snpids.isEmpty()) {
				snpstr = Strings.concat(snpids, Strings.semicolon);
			}
			
			String ln = t + "\t" + traitToSNP.get(t) + "\t" + snpstr;
			outf.writeln(ln);
			
			
		}
		outf.close();
		
	}
	
	public void combinefiles(String correlfile, String log, String out) throws IOException {
		TextFile tf = new TextFile(correlfile, TextFile.R);
		String[] elems = tf.readLineElems(TextFile.tab);
		HashMap<String, String> traitToStr = new HashMap<>();
		while (elems != null) {
			String trait = elems[0];
			traitToStr.put(trait, Strings.concat(elems, Strings.tab));
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		
		TextFile tf2 = new TextFile(log, TextFile.R);
		TextFile outf = new TextFile(out, TextFile.W);
		elems = tf2.readLineElems(TextFile.tab);
		while (elems != null) {
			String trait = elems[0];
			String input = elems[1];
			boolean containsinput = false;
			if (elems[2].contains(input)) {
				containsinput = true;
			}
			String[] split = elems[2].split(";");
			
			outf.writeln(Strings.concat(elems, Strings.tab) + "\t" + traitToStr.get(trait) + "\t" + containsinput + "\t" + split.length);
			elems = tf2.readLineElems(TextFile.tab);
		}
		tf2.close();
		outf.close();
	}
	
	public void extractSNP(String ttdir, String snpselect, String out) throws IOException {
		
		TriTyperGenotypeData ds = new TriTyperGenotypeData();
		ds.displayWarnings = false;
		ds.load(ttdir);
		
		TextFile tf = new TextFile(snpselect, TextFile.R);
		ArrayList<String> snps = tf.readAsArrayList();
		tf.close();
		
		TextFile tfout = new TextFile(out + "mat.txt", TextFile.W);
		String header = "Variant\t" + Strings.concat(ds.getIndividuals(), Strings.tab);
		tfout.writeln(header);
		SNPLoader loader = ds.createSNPLoader();
		for (String s : snps) {
			Integer id = ds.getSnpToSNPId().get(s);
			if (id != null && id > -1) {
				SNP obj = ds.getSNPObject(id);
				loader.loadGenotypes(obj);
				double[] arr = Primitives.toDoubleArr(obj.getGenotypes());
				String outln = s + "-gt\t" + Strings.concat(arr, Strings.tab);
				tfout.writeln(outln);
				if (loader.hasDosageInformation()) {
					loader.loadDosage(obj);
					outln = s + "-ds\t" + Strings.concat(obj.getDosageValues(), Strings.tab);
					tfout.writeln(outln);
					
					arr = Primitives.toDoubleArr(obj.getGenotypes());
					outln = s + "-gtds\t" + Strings.concat(arr, Strings.tab);
					tfout.writeln(outln);
				}
				
			}
		}
		loader.close();
		tfout.close();
		
	}
	
	public void traitfilter(String ttdataset, String traitlist, String traitlistout) throws IOException {
		TextFile tf = new TextFile(ttdataset + "SNPs.txt", TextFile.R);
		ArrayList<String> snplist = tf.readAsArrayList();
		HashSet<String> snpset = new HashSet<String>();
		snpset.addAll(snplist);
		tf.close();
		
		System.out.println(snplist.size() + " snps in file");
		
		TextFile tf2 = new TextFile(traitlist, TextFile.R);
		TextFile out = new TextFile(traitlistout, TextFile.W);
		String[] elems = tf2.readLineElems(TextFile.tab);
		while (elems != null) {
			
			String snp = elems[1];
			if (snpset.contains(snp)) {
				out.writeln(Strings.concat(elems, Strings.tab));
			}
			elems = tf2.readLineElems(TextFile.tab);
		}
		tf2.close();
		out.close();
	}
	
	
	public void correlate(String traitdsloc, String gtdsloc, String traitToSNPCoupling, String out) throws IOException {
		
		
		HashMap<String, String> traitToSNP = new HashMap<>();
		TextFile tf = new TextFile(traitToSNPCoupling, TextFile.R);
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			String trait = elems[0];
			String snp = elems[1];
			traitToSNP.put(trait, snp);
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		
		System.out.println(traitToSNP.size() + "  combos loaded from: " + traitToSNPCoupling);
		
		TriTyperGenotypeData traitDs = new TriTyperGenotypeData();
		traitDs.displayWarnings = false;
		traitDs.load(traitdsloc);
		
		String[] traits = traitDs.getSNPs();
		int nroverlap = 0;
		for (String trait : traits) {
			if (traitToSNP.containsKey(trait)) {
				nroverlap++;
			}
		}
		if (nroverlap == 0) {
			System.out.println("No overlap found between " + traitdsloc + " and " + traitToSNPCoupling);
			System.exit(-1);
		}
		
		TriTyperGenotypeData gtDs = new TriTyperGenotypeData();
		gtDs.displayWarnings = false;
		gtDs.load(gtdsloc);
		String[] gtInds = gtDs.getIndividuals();
		
		// determine if the individuals overlap
		if (gtInds.length != traitDs.getIndividuals().length) {
			System.out.println("Numver of individuals is not equal");
			System.exit(-1);
		}
		
		int overlap = 0;
		for (String gtInd : gtInds) {
			Integer id = traitDs.getIndividualId(gtInd);
			if (id != null && id > -1) {
				overlap++;
			} else {
				System.out.println("Individual: " + gtInd + "  not found in " + gtdsloc);
			}
		}
		
		if (overlap < gtInds.length) {
			System.exit(-1);
		}
		
		int[] traitIndToGtInd = new int[gtInds.length];
		String[] traitInds = traitDs.getIndividuals();
		for (int i = 0; i < traitInds.length; i++) {
			Integer id = gtDs.getIndividualId(traitInds[i]);
			traitIndToGtInd[i] = id;
		}
		
		
		SNPLoader gtDsLoader = gtDs.createSNPLoader();
		SNPLoader traitDsLoader = traitDs.createSNPLoader();
		TextFile outf = new TextFile(out, TextFile.W);
		String header = "trait\ttraitmaf\ttraithwep\tsnp\tsnpmaf\tsnphwep\tgtCorrel\tdsCorrel";
		outf.writeln(header);
		for (int t = 0; t < traits.length; t++) {
			String trait = traits[t];
			String snp = traitToSNP.get(trait);
			int snpId = gtDs.getSnpToSNPId().get(snp);
			if (snpId > -1) {
				System.out.println("Matched " + trait + " to snp " + snp);
				SNP snpobj = gtDs.getSNPObject(snpId);
				SNP traitobj = traitDs.getSNPObject(t);
				
				gtDsLoader.loadGenotypes(snpobj);
				if (gtDsLoader.hasDosageInformation()) {
					gtDsLoader.loadDosage(snpobj);
				}
				
				traitDsLoader.loadGenotypes(traitobj);
				if (traitDsLoader.hasDosageInformation()) {
					traitDsLoader.loadDosage(traitobj);
				}
				
				double[] xgt = new double[gtInds.length];
				double[] ygt = new double[gtInds.length];
				double[] xdos = new double[gtInds.length];
				double[] ydos = new double[gtInds.length];
				
				PearsonsCorrelation pc = new PearsonsCorrelation();
				
				
				double[] gtdos = snpobj.getDosageValues();
				byte[] gtgt = snpobj.getGenotypes();
				double[] trdos = traitobj.getDosageValues();
				byte[] trgt = traitobj.getGenotypes();
				
				for (int i = 0; i < traitInds.length; i++) {
					int gtId = traitIndToGtInd[i];
					
					xdos[i] = gtdos[gtId];
					ydos[i] = trdos[i];
					xgt[i] = gtgt[gtId];
					ygt[i] = trgt[i];
				}
				
				double r1 = pc.correlation(xdos, ydos);
				double r2 = pc.correlation(xgt, ygt);
				
				String outln = trait + "\t" + traitobj.getMAF() + "\t" + traitobj.getHWEP() + "\t" + snp + "\t" + snpobj.getMAF() + "\t" + snpobj.getHWEP() + "\t" + r2 + "\t" + r1;
				outf.writeln(outln);
				System.out.println(outln);
				snpobj.clearGenotypes();
				traitobj.clearGenotypes();
				
			}
			
			
		}
		gtDsLoader.close();
		traitDsLoader.close();
		outf.close();
	}
	
}
