package eqtlmappingpipeline.util;

import com.mastfrog.util.streams.HashingOutputStream;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.bin.BinaryFile;
import umcg.genetica.io.text.TextFile;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;

public class TriTyperSplitter {
	
	
	public void split(String in, String dirOut) throws IOException {
		
		System.out.println("Splitting: " + in);
		TextFile tf = new TextFile(in + "SNPs.txt", TextFile.R);
		ArrayList<String> snps = tf.readAsArrayList();
		new ArrayList<String>();
		HashMap<String, Integer> snpToInt = new HashMap<String, Integer>();
		for (int i = 0; i < snps.size(); i++) {
			snpToInt.put(snps.get(i), i);
		}
		tf.close();
		System.out.println(snps.size() + " snps loaded. ");
		
		tf = new TextFile(in + "Individuals.txt", TextFile.R);
		ArrayList<String> individuals = tf.readAsArrayList();
		tf.close();
		System.out.println(individuals.size() + " individuals.");
		
		int[] chrs = new int[snps.size()];
		int[] pos = new int[snps.size()];
		TextFile tf2 = new TextFile(in + "SNPMappings.txt", TextFile.R);
		String[] elems = tf2.readLineElems(TextFile.tab);
		while (elems != null) {
			Integer id = snpToInt.get(elems[2]);
			if (id != null) {
				chrs[id] = Integer.parseInt(elems[0]);
				pos[id] = Integer.parseInt(elems[1]);
				if (chrs[id] < 0 || chrs[id] > 22) {
					chrs[id] = -1;
				}
			}
			elems = tf2.readLineElems(TextFile.tab);
		}
		tf2.close();
		
		
		int nrBytesPerGt = 2;
		int nrBytesPerGtDosage = 1;
		int nrBytes = individuals.size() * nrBytesPerGt;
		int nrBytesDosage = individuals.size() * nrBytesPerGtDosage;
		{
			BufferedInputStream bf = new BufferedInputStream(new FileInputStream(in + "GenotypeMatrix.dat"));
			
			BufferedOutputStream[] bfout = new BufferedOutputStream[22];
			TextFile[] outSNPs = new TextFile[22];
			TextFile[] outSNPAnnot = new TextFile[22];
			for (int i = 0; i < bfout.length; i++) {
				String dir = dirOut + "chr" + (i + 1) + "/";
				Gpio.createDir(dir);
				bfout[i] = new BufferedOutputStream(new FileOutputStream(dir + "GenotypeMatrix.dat"));
				
				outSNPs[i] = new TextFile(dir + "SNPs.txt.gz", TextFile.W);
				outSNPAnnot[i] = new TextFile(dir + "SNPMappings.txt.gz", TextFile.W);
				Gpio.copyFile(in + "Individuals.txt", dir + "Individuals.txt");
				Gpio.copyFile(in + "PhenotypeInformation.txt", dir + "PhenotypeInformation.txt");
			}
			
			byte[] b = new byte[nrBytes];
			ProgressBar pb = new ProgressBar(snps.size(), "Processing GenotypeMatrix.dat");
			for (int i = 0; i < snps.size(); i++) {
				bf.read(b);
				
				int chr = chrs[i];
				if (chr > 0) {
					int chrindex = chr - 1;
					String snp = snps.get(i);
					bfout[chrindex].write(b);
					outSNPs[chrindex].writeln(snp);
					outSNPAnnot[chrindex].writeln(chr + "\t" + pos[i] + "\t" + snp);
				}
				pb.iterate();
			}
			
			for (int i = 0; i < bfout.length; i++) {
				bfout[i].close();
				outSNPs[i].close();
				outSNPAnnot[i].close();
			}
			pb.close();
		}
		
		// process dosages, if any
		if (Gpio.exists(in + "ImputedDosageMatrix.dat")) {
			
			BufferedInputStream bf = new BufferedInputStream(new FileInputStream(in + "ImputedDosageMatrix.dat"));
			
			BufferedOutputStream[] bfout = new BufferedOutputStream[22];
			
			for (int i = 0; i < bfout.length; i++) {
				String dir = dirOut + "chr" + (i + 1) + "/";
				Gpio.createDir(dir);
				bfout[i] = new BufferedOutputStream(new FileOutputStream(dir + "ImputedDosageMatrix.dat"));
				
				
			}
			
			byte[] b = new byte[nrBytesDosage];
			ProgressBar pb = new ProgressBar(snps.size(), "Processing ImputedDosageMatrix.dat");
			for (int i = 0; i < snps.size(); i++) {
				bf.read(b);
				
				int chr = chrs[i];
				if (chr > 0) {
					int chrindex = chr - 1;
					bfout[chrindex].write(b);
				}
				pb.iterate();
			}
			
			for (int i = 0; i < bfout.length; i++) {
				bfout[i].close();
			}
			pb.close();
		}
		
	}
	
}


