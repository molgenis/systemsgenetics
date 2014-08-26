/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.binarymeta.util;

import com.lowagie.text.DocumentException;
import eqtlmappingpipeline.metaqtl3.FDR;
import eqtlmappingpipeline.metaqtl3.FDR.FDRMethod;
import eqtlmappingpipeline.metaqtl3.graphics.EQTLDotPlot;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.logging.Level;
import java.util.logging.Logger;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author harmjan
 */
public class Filter {

    public void run(String indir, String probetranslationfile, String requestedAnnotation) throws IOException {

	if (!indir.endsWith("/")) {
	    indir += "/";
	}
	if (!Gpio.isDir(indir)) {
	    throw new IOException("Directory " + indir + " does not exist.");
	}
	if (!Gpio.exists(indir + "eQTLs.txt")) {
	    throw new IOException("Expected to find eQTLs.txt in directory " + indir);
	}


	String outdir = indir + "filtered/";
	Gpio.createDir(outdir);

	System.out.println("In:\t" + indir);
	System.out.println("Out:\t" + outdir);

	String[] filelist = Gpio.getListOfFiles(indir, "gz");

	if (filelist == null) {
	    throw new IOException("No files found ending with .gz in folder " + indir);
	}

	ArrayList<String> permutedFileList = new ArrayList<String>();


	for (String f : filelist) {
	    if (f.contains("PermutedEQTLsPermutationRound") && !f.contains("OppositeEffects")) {
		permutedFileList.add(f);
	    }
	}

	System.out.println("Detected " + permutedFileList.size() + " permutations in directory");

	HashSet<String> allowedProbes = readListOfAllowedProbes(probetranslationfile, requestedAnnotation);




	int minimumNrOfEQTLs = filterFile(indir + "eQTLs.txt", outdir + "eQTLs.txt", allowedProbes);



	for (int i = 0; i < permutedFileList.size(); i++) {
	    String file = permutedFileList.get(i);
	    String filename = file.substring(indir.length());

	    int nrlines = filterFile(file, outdir + filename, allowedProbes);
	    if(nrlines < minimumNrOfEQTLs){
		minimumNrOfEQTLs = nrlines;
	    }
	}


	FDR f = new FDR();
	f.calculateFDR(outdir, permutedFileList.size(), minimumNrOfEQTLs, 0.05, true, null, null, FDRMethod.ALL, false);
	EQTLDotPlot edp = new EQTLDotPlot();
        try {
            edp.draw(outdir + "/eQTLsFDR" + 0.05 + ".txt", outdir + "/DotPlot-FDR" + 0.05 + ".png", EQTLDotPlot.Output.PDF); // "/eQTLsFDR" + fdrCutOff + ".txt", outputReportsDir + "/eQTLsFDR" + fdrCutOff + "DotPlot.png"
        } catch (DocumentException ex) {
            Logger.getLogger(Filter.class.getName()).log(Level.SEVERE, null, ex);
        }
	edp = null;

    }

    private int filterFile(String file, String outfile, HashSet<String> allowedProbes) throws IOException {

	System.out.println("Filtering " + file);
	TextFile in = new TextFile(file, TextFile.R);
	TextFile out = new TextFile(outfile, TextFile.W);

	String line = in.readLine();
	out.writeln(line);

	int nrLinesIncluded = 0;
	while (line != null) {

	    String[] elems = line.split("\t");
	    String probe = elems[4];
	    if (allowedProbes.contains(probe)) {
		out.writeln(line);
		nrLinesIncluded++;
	    }
	    line = in.readLine();
	}

	in.close();
	out.close();
	System.out.println(nrLinesIncluded+" remaining");
	return nrLinesIncluded;
    }

    private HashSet<String> readListOfAllowedProbes(String probetranslationfile, String requestedAnnotation) throws IOException {
	TextFile tf = new TextFile(probetranslationfile, TextFile.R);
	String[] elems = tf.readLineElems(TextFile.tab);

	int colWithProbeAnnotation = -1;
	int i = 0;
	for (String s : elems) {
	    if (s.equals(requestedAnnotation)) {
		colWithProbeAnnotation = i;
	    }
	    i++;
	}
	if (colWithProbeAnnotation == -1) {
	    System.out.println("Requested probe annotation " + requestedAnnotation + " not found in " + probetranslationfile);
	    return null;
	}

	HashSet<String> allowedProbes = new HashSet<String>();
	while (elems != null) {

	    String probe = elems[0];
	    String arrayAddress = elems[colWithProbeAnnotation];

	    if (!arrayAddress.equals("-")) {
		allowedProbes.add(probe);
	    }


	    elems = tf.readLineElems(TextFile.tab);
	}

	tf.close();

	System.out.println(allowedProbes.size() + " probes allowed");

	return allowedProbes;
    }
}
